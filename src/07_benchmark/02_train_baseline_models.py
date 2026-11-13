"""Train the baseline models of every benchmark task under the naive splits and the patient-safe splits, and score them within every test unit."""
import hashlib
import importlib
import warnings

C = importlib.import_module("01_audit_split_leakage")

REASONS = {"value_scale_not_comparable": "Its gene matrix is not on a log scale, so the correction stage gave it no aligned values.",
           "study_has_too_few_samples": "Its study unit has too few estimation samples for the correction stage, which did not align it.",
           "sample_type_not_corrected": "Its sample type is in no group of the correction stage, so it has no aligned values."}
UNIT_REASONS = {"fewer_than_minimum_of_a_class": "The unit holds fewer than splits.minimum_per_class records with values of one of the two classes.",
                "fewer_than_minimum_with_a_reference_score": "Fewer than splits.minimum_per_class records of a class of the unit have a Sepsis MetaScore "
                                                             "from the validation stage, which scores only strata that hold both of its groups."}
REFERENCE = "sepsis_metascore"
REFERENCE_TASK = "sepsis_vs_healthy_control"
RECORD_COLUMNS = ["task", *C.KEYS, "sample_id", "patient_key", "patient_identified", "patient_safe_group", "label", "in_gts", "has_values", "reason"]
UNIT_COLUMNS = ["task", "model", "population", "scheme", "fold", "unit", "train_records", "train_positives", "train_negatives", "genes", "records",
                "predicted", "positives", "negatives", "status", "reason", "auroc", "auroc_low", "auroc_high", "auprc", "auprc_low", "auprc_high", "prevalence"]
PREDICTION_COLUMNS = ["task", "model", "scheme", "fold", "unit", *C.KEYS, "label", "score"]
GENE_COLUMNS = ["task", "scheme", "fold", "rank", "hgnc_id", "statistic", "present_share"]
SUMMARY_COLUMNS = ["task", "model", "population", "scheme", "units", "units_scored", "records_scored", "auroc_mean", "auroc_median", "auprc_mean",
                   "auprc_median", "prevalence_mean", "inflation_auroc", "inflation_auprc"]
METRICS = ("auroc", "auroc_low", "auroc_high", "auprc", "auprc_low", "auprc_high", "prevalence")
SUBSETS = ("all", "without_twins")
FLAGS = ("in_additive_set", "in_twins_set", "in_seen_study_set", "in_twins_across_set")
STUDY_COLUMNS = ["task", "model", "scheme", "fold", "unit", "study_accession", "subset", "records", "positives", "negatives", "status", "reason", "auroc", *FLAGS]
EFFECTS = {"total_inflation": ("pooled AUROC of the naive random split minus the patient-safe leave-one-group-out AUROC, as in model_summary", False),
           "study_shortcut": ("A - B: pooled minus within-study AUROC of the naive random split, over every study of its test sets", False),
           "twins": ("B - C: within-study AUROC of the naive random split with minus without twins, over the study repeats scorable both ways", True),
           "seen_study": ("C - D: within-study AUROC without twins minus that of leave-one-group-out, over the studies scorable in both", True),
           "twins_across_studies": ("E - F: leave-one-study-out within the held-out study with minus without twins, over the studies scorable both ways", True)}
EFFECT_COLUMNS = ["task", "model", "effect", "comparison", "matched", "status", "reason", "first", "second", "estimate", "low", "high", "studies",
                  "first_records", "second_records"]
VALUES = ["a_pooled_random_split", "b_within_study_random_split", "c_within_study_random_split_without_twins", "d_within_study_leave_one_group_out",
          "study_confounding", "twins", "own_study_in_training", "total", "e_within_study_leave_one_study_out", "f_within_study_leave_one_study_out_without_twins",
          "twins_across_studies"]
MATCHED_COUNTS = ["matched_studies", "a_records", "b_records", "c_records", "d_records", "loso_matched_studies", "e_records", "f_records"]
ADDITIVE_COLUMNS = ["task", "model", "status", "reason", "loso_status", "loso_reason", *VALUES, *MATCHED_COUNTS]
EFFECT_REASONS = {"fewer_than_minimum_studies": "Fewer than metrics.inflation_effects.minimum_studies studies are in the set of the comparison."}


def task_records(populations, correction):
    """Every record of the naive population of every task, whether it is also a sample of the task and whether it has aligned values, or why not."""
    import numpy as np

    group = C.CONFIG["features"]["group"]
    naive = populations[populations["population"] == "naive"]
    safe = populations[populations["population"] == "gts"]
    marked = naive.merge(safe[["task", *C.KEYS]].assign(in_gts=True), on=["task", *C.KEYS], how="left", validate="one_to_one")
    C.require(len(marked) == len(naive) and int(marked["in_gts"].notna().sum()) == len(safe), "Every sample of a task must be in its naive population")
    marked["in_gts"] = marked["in_gts"].eq(True)
    held = correction[[*C.KEYS, "group", "status", "reason"]].rename(columns={"reason": "correction_reason"})
    marked = marked.merge(held, on=C.KEYS, how="left", validate="many_to_one")
    C.require(marked["status"].notna().all(), "Every record of a task must be in the records of the correction stage")
    marked["has_values"] = (marked["group"] == group) & (marked["status"] == "in_unit")
    marked["reason"] = np.where(marked["has_values"], "", np.where(marked["correction_reason"] != "", marked["correction_reason"], "sample_type_not_corrected"))
    C.require(set(marked["reason"]) <= set(REASONS) | {""}, f"Unknown reason for a record without values; name it in REASONS: {sorted(set(marked['reason']))}")
    return marked[RECORD_COLUMNS].reset_index(drop=True)


def feature_values(frame):
    """The aligned values of the records of a frame, records in rows and the genes of the feature group in columns, as float32, and the genes."""
    import numpy as np
    import pandas as pd

    group = C.CONFIG["features"]["group"]
    units = frame.groupby(["study_accession", "platform"]).indices
    values, genes = None, None
    for study, platform in sorted(units):
        positions = units[(study, platform)]
        columns = list(frame["sample_accession"].to_numpy()[positions])
        part = pd.read_parquet(C.CORRECTION / "aligned" / group / f"{study}__{platform}.parquet", columns=["hgnc_id", *columns])
        if genes is None:
            genes = list(part["hgnc_id"])
            values = np.full((len(frame), len(genes)), np.nan, dtype="float32")
        C.require(list(part["hgnc_id"]) == genes, f"The aligned table of {study} {platform} holds other genes than the other units of {group}")
        values[positions] = part[columns].to_numpy(dtype="float32").T
    return values, genes


def select_genes(values, rows, labels):
    """The genes of one training set: of the genes with values in at least minimum_present_share of its records, the first `genes` by the absolute Welch
    t statistic between its classes, ties in the order of the genes; returns their columns, statistics and present shares."""
    import numpy as np

    policy = C.CONFIG["features"]
    statistic, present = np.full(values.shape[1], np.nan), np.zeros(values.shape[1])
    cases = np.asarray(labels) == 1
    for start in range(0, values.shape[1], policy["chunk_genes"]):
        part = values[rows, start:start + policy["chunk_genes"]].astype("float64")
        held = ~np.isnan(part)
        present[start:start + part.shape[1]] = held.mean(axis=0)
        moments = []
        for side in (cases, ~cases):
            own, mask = part[side], held[side]
            count = mask.sum(axis=0)
            mean = np.where(mask, own, 0.0).sum(axis=0) / np.maximum(count, 1)
            spread = np.where(mask, (own - mean) ** 2, 0.0).sum(axis=0) / np.maximum(count - 1, 1)
            moments.append((count, mean, spread))
        (first, centre, spread), (second, other, scatter) = moments
        error = np.sqrt(spread / np.maximum(first, 1) + scatter / np.maximum(second, 1))
        valid = (first >= 2) & (second >= 2) & (error > 0)
        statistic[start:start + part.shape[1]] = np.where(valid, (centre - other) / np.where(valid, error, 1.0), np.nan)
    candidate = (present >= policy["minimum_present_share"]) & np.isfinite(statistic)
    C.require(int(candidate.sum()) >= policy["genes"], "Too few genes with values in the training records to select features.genes from")
    order = np.argsort(-np.where(candidate, np.abs(statistic), -1.0), kind="stable")[:policy["genes"]]
    return order, statistic[order], present[order]


def build_model(name):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    policy = C.CONFIG["models"][name]
    if name == "logistic_regression":
        return make_pipeline(SimpleImputer(strategy="mean"), StandardScaler(), LogisticRegression(C=policy["C"], max_iter=policy["max_iter"], tol=policy["tol"]))
    return make_pipeline(SimpleImputer(strategy="mean"), HistGradientBoostingClassifier(
        learning_rate=policy["learning_rate"], max_iter=policy["max_iter"], max_leaf_nodes=policy["max_leaf_nodes"], min_samples_leaf=policy["min_samples_leaf"],
        l2_regularization=policy["l2_regularization"], early_stopping=False, random_state=policy["random_state"]))


def fit_predict(name, train_x, train_y, test_x):
    """Fit one model on a training set and return the probability of the positive class for the test records; a model that does not converge stops the run."""
    from sklearn.exceptions import ConvergenceWarning

    model = build_model(name)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        try:
            model.fit(train_x, train_y)
        except ConvergenceWarning as error:
            raise ValueError(f"{name} did not converge ({error}); raise models.{name}.max_iter in benchmark_config.json") from error
    return model.predict_proba(test_x)[:, 1].astype("float64")


def task_folds(rows, splits, task):
    """Every fold of one task as (population, scheme, fold, training mask, [(unit, test mask)]) over its naive records; masks keep only records with values."""
    import numpy as np

    held = rows["has_values"].to_numpy()
    studies, groups, safe = rows["study_accession"].to_numpy(), rows["patient_safe_group"].to_numpy(), rows["in_gts"].to_numpy()
    for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1):
        test = C.record_split(rows, repeat)
        yield "naive", "naive_record_split", C.repeat_name(repeat), ~test & held, [("pooled", test & held)]
    for study in sorted(set(studies)):
        yield "naive", "naive_leave_one_study_out", study, (studies != study) & held, [(study, (studies == study) & held)]
    for scheme, fold, train_groups, test_groups, _ in C.gts_folds(splits, task):
        yield "gts", scheme, fold, np.isin(groups, train_groups) & safe & held, [(group, (groups == group) & safe & held) for group in test_groups]


def study_scores(records, predictions, splits):
    """The AUROC within every study of the test records of every unit of the naive random splits, the naive leave-one-study-out split and the patient-safe
    scheme of metrics.inflation_reference: over all of them, and for a naive scheme over those without a linked record among the training records of
    their fold."""
    import numpy as np
    import pandas as pd
    from sklearn.metrics import roc_auc_score

    grouped, rows = predictions.groupby(["task", "scheme", "fold"], sort=False), []
    schemes = (*C.NAIVE_SCHEMES, C.CONFIG["metrics"]["inflation_reference"])
    for task, part in records.groupby("task", sort=False):
        part = part.reset_index(drop=True)
        pairs, where = C.twin_pairs(part), pd.MultiIndex.from_frame(part[["study_accession", "sample_accession"]])
        for _, scheme, fold, train, _ in task_folds(part, splits, task):
            if scheme not in schemes or (task, scheme, fold) not in grouped.groups:
                continue
            found = grouped.get_group((task, scheme, fold))
            place = where.get_indexer(pd.MultiIndex.from_frame(found[["study_accession", "sample_accession"]]))
            C.require((place >= 0).all(), f"Every prediction of {task} must be a record of its task")
            test = np.zeros(len(part), dtype=bool)
            test[place] = True
            # Exposure is that of the leakage audit: a deposit of the same array, or a sample of the same identified patient, among the fold's training records.
            twin = C.exposed(pairs, test, train)["any"][place]
            subsets = (("all", np.ones(len(found), dtype=bool)), ("without_twins", ~twin)) if scheme in C.NAIVE_SCHEMES else (("all", np.ones(len(found), dtype=bool)),)
            for (model, unit), own in found.assign(row=np.arange(len(found))).groupby(["model", "unit"], sort=False):
                for study in sorted(set(own["study_accession"])):
                    for subset, keep in subsets:
                        chosen = own[(own["study_accession"] == study) & keep[own["row"].to_numpy()]]
                        labels = chosen["label"].to_numpy(dtype="int64")
                        ready = C.scorable(labels)
                        rows.append(dict(task=task, model=model, scheme=scheme, fold=fold, unit=unit, study_accession=study, subset=subset, records=len(chosen),
                                         positives=int(labels.sum()), negatives=int((labels == 0).sum()), status="scored" if ready else "not_scored",
                                         reason="" if ready else "fewer_than_minimum_of_a_class",
                                         auroc=float(roc_auc_score(labels, chosen["score"])) if ready else np.nan, **dict.fromkeys(FLAGS, False)))
    return pd.DataFrame(rows, columns=STUDY_COLUMNS)


def effect_seed(task, model, effect):
    """The seed of the bootstrap of one effect, from the configured seed and the names of the task, the model and the effect."""
    text = f"{C.CONFIG['metrics']['inflation_effects']['bootstrap']['seed']}:{task}:{model}:{effect}"
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:16], 16)


def weighted_auroc(labels, groups, negatives, weights):
    """AUROC of records that count weights times each, ties as half; groups are the indices of the distinct scores in increasing order, and negatives the
    zero weights of the positives."""
    import numpy as np

    below = np.bincount(groups, weights=weights * negatives, minlength=groups.max() + 1)
    under = np.cumsum(below) - below
    positive = weights * (labels == 1)
    return float((positive * (under[groups] + 0.5 * below[groups])).sum() / (positive.sum() * (weights * negatives).sum()))


def interval(draws):
    """A percentile interval of the bootstrap replicates of an effect; a replicate whose studies leave a side without a scorable value is left out."""
    import numpy as np

    tail = (1 - C.CONFIG["metrics"]["inflation_effects"]["bootstrap"]["level"]) / 2
    draws = np.asarray(draws, dtype="float64")
    draws = draws[np.isfinite(draws)]
    return float(np.quantile(draws, tail)), float(np.quantile(draws, 1 - tail))


def shortcut_effect(scored, pooled, seed):
    """The study shortcut, A - B: the mean over the repeats of the pooled AUROC over all test records, minus the mean over the repeats of the mean
    within-study AUROC over the scored studies; the bootstrap resamples the studies of the test sets and recomputes both."""
    import numpy as np

    policy, minimum = C.CONFIG["metrics"]["inflation_effects"]["bootstrap"], C.CONFIG["splits"]["minimum_per_class"]
    studies = sorted(set(pooled["study_accession"]))
    index = {study: position for position, study in enumerate(studies)}
    repeats = sorted(set(pooled["fold"]))
    within = scored[(scored["scheme"] == "naive_record_split") & (scored["subset"] == "all")]
    values = np.full((len(studies), len(repeats)), np.nan)
    for row in within.itertuples(index=False):
        values[index[row.study_accession], repeats.index(row.fold)] = row.auroc
    held = np.isfinite(values)
    parts = []
    for repeat in repeats:
        own = pooled[pooled["fold"] == repeat]
        labels = own["label"].to_numpy(dtype="int64")
        parts.append((labels, np.unique(own["score"].to_numpy(dtype="float64"), return_inverse=True)[1], (labels == 0).astype("float64"),
                      own["study_accession"].map(index).to_numpy()))

    def effect(counts):
        found = []
        for labels, groups, negatives, owner in parts:
            weights = counts[owner].astype("float64")
            if (weights * (labels == 1)).sum() >= minimum and (weights * negatives).sum() >= minimum:
                found.append(weighted_auroc(labels, groups, negatives, weights))
        weight = (counts[:, None] * held).sum(axis=0)
        means = (counts[:, None] * np.where(held, values, 0.0)).sum(axis=0)[weight > 0] / weight[weight > 0]
        return float(np.mean(found)) if found else np.nan, float(np.mean(means)) if len(means) else np.nan

    first, second = effect(np.ones(len(studies), dtype="int64"))
    generator = np.random.default_rng(seed)
    draws = generator.integers(0, len(studies), (policy["replicates"], len(studies)))
    replicates = [np.subtract(*effect(np.bincount(draw, minlength=len(studies)))) for draw in draws]
    return dict(first=first, second=second, estimate=first - second, studies=len(studies), first_records=len(pooled),
                second_records=int(within["records"].sum())), interval(replicates)


def paired_effect(first, second, seed):
    """An effect measured on one set of studies: the mean over its studies of the first minus the second value, with a bootstrap that resamples them."""
    import numpy as np

    policy = C.CONFIG["metrics"]["inflation_effects"]["bootstrap"]
    gaps = np.asarray(first, dtype="float64") - np.asarray(second, dtype="float64")
    generator = np.random.default_rng(seed)
    draws = generator.integers(0, len(gaps), (policy["replicates"], len(gaps)))
    return float(np.mean(gaps)), interval(gaps[draws].mean(axis=1))


def effect_sets(scored):
    """The set of every matched effect of one task and model, with its per-study values and records: twins over the study repeats scorable with and
    without twins, the seen study over the studies scorable without twins in a repeat and in leave-one-group-out, and twins across studies over the
    held-out studies scorable both ways; and the per-study rows each set uses."""
    reference = C.CONFIG["metrics"]["inflation_reference"]

    def pick(scheme, subset):
        return scored[(scored["scheme"] == scheme) & (scored["subset"] == subset)].set_index(["fold", "study_accession"])

    first, second = pick("naive_record_split", "all"), pick("naive_record_split", "without_twins")
    pairs = first.index.intersection(second.index)
    twins = (first.loc[pairs, "auroc"].groupby(level="study_accession").mean(), second.loc[pairs, "auroc"].groupby(level="study_accession").mean(),
             int(first.loc[pairs, "records"].sum()), int(second.loc[pairs, "records"].sum()))
    safe = scored[(scored["scheme"] == reference) & (scored["subset"] == "all")]
    # A study whose platforms fall in two patient-safe groups is held out in two folds; its value is the mean over the folds that score it.
    held = safe.groupby("study_accession").agg(auroc=("auroc", "mean"), records=("records", "sum"))
    without = second.reset_index()
    seen = sorted(set(without["study_accession"]) & set(held.index))
    kept = without[without["study_accession"].isin(seen)]
    seen_study = (kept.groupby("study_accession")["auroc"].mean(), held.loc[seen, "auroc"], int(kept["records"].sum()), int(held.loc[seen, "records"].sum()))
    one, other = pick("naive_leave_one_study_out", "all"), pick("naive_leave_one_study_out", "without_twins")
    shared = one.index.intersection(other.index)
    across = (one.loc[shared, "auroc"].droplevel("fold"), other.loc[shared, "auroc"].droplevel("fold"), int(one.loc[shared, "records"].sum()),
              int(other.loc[shared, "records"].sum()))
    used = {"in_twins_set": {("naive_record_split", fold, study) for fold, study in pairs},
            "in_seen_study_set": {("naive_record_split", fold, study) for fold, study in zip(kept["fold"], kept["study_accession"])}
                                 | {(reference, fold, study) for fold, study in zip(safe["fold"], safe["study_accession"]) if study in seen},
            "in_twins_across_set": {("naive_leave_one_study_out", fold, study) for fold, study in shared}}
    return {"twins": twins, "seen_study": seen_study, "twins_across_studies": across}, used


def additive_matched(scored, pooled):
    """The sensitivity analysis: an additive split over one set of studies S, those scorable within study in their leave-one-group-out fold and, in at least
    one repeat of the random split, both with and without twins. B_s and C_s are the means over the repeats in which the study is scorable in both, D_s
    its value in its leave-one-group-out folds, and A each repeat's pooled AUROC over the test records of S; leave-one-study-out is matched on the held-out
    studies scorable with and without twins. Returns the values, their records and the per-study rows used."""
    import numpy as np
    from sklearn.metrics import roc_auc_score

    reference, minimum = C.CONFIG["metrics"]["inflation_reference"], C.CONFIG["metrics"]["inflation_effects"]["minimum_studies"]

    def pick(scheme, subset):
        return scored[(scored["scheme"] == scheme) & (scored["subset"] == subset)].set_index(["fold", "study_accession"])

    first, second = pick("naive_record_split", "all"), pick("naive_record_split", "without_twins")
    both = first.index.intersection(second.index)
    safe = scored[(scored["scheme"] == reference) & (scored["subset"] == "all")]
    held = safe.groupby("study_accession").agg(auroc=("auroc", "mean"), records=("records", "sum"))
    chosen = sorted(set(both.get_level_values("study_accession")) & set(held.index))
    pairs = both[both.get_level_values("study_accession").isin(chosen)]
    one, other = pick("naive_leave_one_study_out", "all"), pick("naive_leave_one_study_out", "without_twins")
    shared = one.index.intersection(other.index)
    values = dict(dict(zip(VALUES, [np.nan] * len(VALUES))), matched_studies=len(chosen), loso_matched_studies=len(shared), a_records=0, b_records=0,
                  c_records=0, d_records=0, e_records=0, f_records=0, status="estimable", reason="", loso_status="estimable", loso_reason="")
    used = set()
    if len(chosen) >= minimum:
        b = first.loc[pairs, "auroc"].groupby(level="study_accession").mean()
        c = second.loc[pairs, "auroc"].groupby(level="study_accession").mean()
        d = held.loc[chosen, "auroc"]
        repeats, records = [], 0
        for _, part in pooled[pooled["study_accession"].isin(chosen)].groupby("fold", sort=True):
            if C.scorable(part["label"].to_numpy()):
                repeats.append(float(roc_auc_score(part["label"], part["score"])))
                records += len(part)
        a, b, c, d = float(np.mean(repeats)) if repeats else np.nan, float(b.mean()), float(c.mean()), float(d.mean())
        values.update(a_pooled_random_split=a, b_within_study_random_split=b, c_within_study_random_split_without_twins=c, d_within_study_leave_one_group_out=d,
                      study_confounding=a - b, twins=b - c, own_study_in_training=c - d, total=a - d, a_records=records,
                      b_records=int(first.loc[pairs, "records"].sum()), c_records=int(second.loc[pairs, "records"].sum()),
                      d_records=int(held.loc[chosen, "records"].sum()))
        used |= {("naive_record_split", fold, study) for fold, study in pairs}
        used |= {(reference, fold, study) for fold, study in zip(safe["fold"], safe["study_accession"]) if study in chosen}
    else:
        values.update(status="not_estimable", reason="fewer_than_minimum_studies")
    if len(shared) >= minimum:
        e, f = float(one.loc[shared, "auroc"].mean()), float(other.loc[shared, "auroc"].mean())
        values.update(e_within_study_leave_one_study_out=e, f_within_study_leave_one_study_out_without_twins=f, twins_across_studies=e - f,
                      e_records=int(one.loc[shared, "records"].sum()), f_records=int(other.loc[shared, "records"].sum()))
        used |= {("naive_leave_one_study_out", fold, study) for fold, study in shared}
    else:
        values.update(loso_status="not_estimable", loso_reason="fewer_than_minimum_studies")
    return values, used


def inflation_effects(summary, studies, predictions):
    """The total inflation of every trained model, as in the summary, and three effects, each measured on the comparison that isolates it and not added up:
    the study shortcut over every study, twins and the seen study over their matched sets, and twins across studies over the held-out studies scorable
    both ways, each with a bootstrap interval over its studies; and the additive sensitivity analysis. Returns the effect rows, the additive rows and, for
    every flag of the per-study table, the rows it marks."""
    import numpy as np

    minimum = C.CONFIG["metrics"]["inflation_effects"]["minimum_studies"]
    effects, additive, flags = [], [], {name: set() for name in FLAGS}
    keyed = summary.set_index(["task", "model", "scheme"])
    for (task, model), part in studies[studies["model"].isin(C.MODELS)].groupby(["task", "model"], sort=False):
        scored = part[part["status"] == "scored"]
        pooled = predictions[(predictions["task"] == task) & (predictions["model"] == model) & (predictions["scheme"] == "naive_record_split")]

        def row(effect, **values):
            comparison, matched = EFFECTS[effect]
            base = dict(task=task, model=model, effect=effect, comparison=comparison, matched=matched, status="estimable", reason="", first=np.nan,
                        second=np.nan, estimate=np.nan, low=np.nan, high=np.nan, studies=0, first_records=0, second_records=0)
            return dict(base, **values)

        total = keyed.loc[(task, model, "naive_record_split")]
        safe = keyed.loc[(task, model, C.CONFIG["metrics"]["inflation_reference"])]
        effects.append(row("total_inflation", status="reported", first=total["auroc_mean"], second=safe["auroc_mean"], estimate=total["inflation_auroc"],
                           studies=int(safe["units_scored"]), first_records=int(total["records_scored"]), second_records=int(safe["records_scored"])))
        found, (low, high) = shortcut_effect(scored, pooled, effect_seed(task, model, "study_shortcut"))
        effects.append(row("study_shortcut", low=low, high=high, **found) if found["studies"] >= minimum
                       else row("study_shortcut", status="not_estimable", reason="fewer_than_minimum_studies", studies=found["studies"]))
        sets, used = effect_sets(scored)
        for effect, (first, second, first_records, second_records) in sets.items():
            if len(first) < minimum:
                effects.append(row(effect, status="not_estimable", reason="fewer_than_minimum_studies", studies=len(first)))
                continue
            # The studies are drawn in the order of their accessions, so that the bootstrap does not depend on the order of the rows.
            first = first.sort_index()
            second = second.loc[first.index]
            estimate, (low, high) = paired_effect(first.to_numpy(), second.to_numpy(), effect_seed(task, model, effect))
            effects.append(row(effect, first=float(first.mean()), second=float(second.mean()), estimate=estimate, low=low, high=high, studies=len(first),
                               first_records=first_records, second_records=second_records))
        values, own = additive_matched(scored, pooled)
        additive.append(dict(task=task, model=model, **values))
        for name, keys in {**used, "in_additive_set": own}.items():
            flags[name] |= {(task, model, *key) for key in keys}
    return effects, additive, flags


def unit_row(task, model, population, scheme, fold, unit, train, labels, genes, test, predicted, status, reason, metrics=None):
    import numpy as np

    taken = labels[predicted] if predicted is not None else labels[:0]
    return dict(task=task, model=model, population=population, scheme=scheme, fold=fold, unit=unit, train_records=int(train.sum()) if model != REFERENCE else 0,
                train_positives=int((labels[train] == 1).sum()) if model != REFERENCE else 0,
                train_negatives=int((labels[train] == 0).sum()) if model != REFERENCE else 0, genes=genes, records=int(test.sum()),
                predicted=int(len(taken)), positives=int((taken == 1).sum()), negatives=int((taken == 0).sum()), status=status, reason=reason,
                **(metrics or {name: np.nan for name in METRICS}))


def main():
    import numpy as np
    import pandas as pd

    print(f"[02] Training baseline models v{C.VERSION}", flush=True)
    start = C.begin(2)
    policy = C.cohort_policy()
    populations = pd.read_parquet(C.DATA / "split_populations.parquet")
    correction = pd.read_parquet(C.CORRECTION / "correction_records.parquet")
    splits = pd.read_parquet(C.COHORT / "benchmark_splits.parquet")
    reference = pd.read_parquet(C.VALIDATION / "sepsis_score_samples.parquet", columns=["study_accession", "sample_accession", "platform", "score"])
    C.require(not reference.duplicated(["study_accession", "sample_accession", "platform"]).any(), "Every record has at most one Sepsis MetaScore")
    reference = reference.set_index(["study_accession", "sample_accession", "platform"])["score"]
    records = task_records(populations, correction)
    wanted = records[records["has_values"]].drop_duplicates(C.KEYS).sort_values(["study_accession", "platform", "sample_accession"], kind="stable")
    wanted = wanted.reset_index(drop=True)
    values, genes = feature_values(wanted)
    where = pd.MultiIndex.from_frame(wanted[C.KEYS])
    print(f"[02] {len(wanted):,} records with aligned values over {len(genes):,} genes ({values.nbytes / 1e9:.2f} GB)", flush=True)

    units, predictions, chosen = [], [], []
    for task in [entry["task"] for entry in policy["tasks"]]:
        rows = records[records["task"] == task].reset_index(drop=True)
        labels = rows["label"].to_numpy(dtype="int64")
        place = where.get_indexer(pd.MultiIndex.from_frame(rows[C.KEYS]))
        C.require(((place >= 0) == rows["has_values"].to_numpy()).all(), f"Every record of {task} with values must be in the feature table")
        keyed = list(zip(rows["study_accession"], rows["sample_accession"], rows["platform"]))
        scored_by = rows.assign(score=reference.reindex(keyed).to_numpy()) if task == REFERENCE_TASK else None
        for population, scheme, fold, train, tests in task_folds(rows, splits, task):
            ready = [(unit, test) for unit, test in tests if C.scorable(labels[test])]
            for unit, test in tests:
                if not C.scorable(labels[test]):
                    for model in (*C.MODELS, *([REFERENCE] if scored_by is not None else [])):
                        units.append(unit_row(task, model, population, scheme, fold, unit, train, labels, 0, test, None, "not_scored",
                                              "fewer_than_minimum_of_a_class"))
            if not ready:
                continue
            C.require(C.scorable(labels[train]), f"The training records of fold {fold} of {scheme} for {task} need the minimum of each class")
            columns, statistic, present = select_genes(values, place[train], labels[train])
            chosen.append(pd.DataFrame({"task": task, "scheme": scheme, "fold": fold, "rank": np.arange(1, len(columns) + 1),
                                        "hgnc_id": [genes[column] for column in columns], "statistic": statistic, "present_share": present}))
            train_x = values[np.ix_(place[train], columns)]
            test_rows = np.concatenate([np.flatnonzero(test) for _, test in ready])
            test_x = values[np.ix_(place[test_rows], columns)]
            for model in C.MODELS:
                scores = dict(zip(test_rows, fit_predict(model, train_x, labels[train], test_x)))
                for unit, test in ready:
                    own = np.flatnonzero(test)
                    found = np.array([scores[index] for index in own])
                    metrics = C.score_unit(labels[own], found, C.unit_seed(task, model, scheme, fold, unit))
                    units.append(unit_row(task, model, population, scheme, fold, unit, train, labels, len(columns), test, own, "scored", "", metrics))
                    predictions.append(rows.loc[own, C.KEYS].assign(task=task, model=model, scheme=scheme, fold=fold, unit=unit, label=labels[own], score=found))
            if scored_by is None:
                continue
            for unit, test in ready:
                own = np.flatnonzero(test & scored_by["score"].notna().to_numpy())
                if not C.scorable(labels[own]):
                    units.append(unit_row(task, REFERENCE, population, scheme, fold, unit, train, labels, 0, test, own, "not_scored",
                                          "fewer_than_minimum_with_a_reference_score"))
                    continue
                found = scored_by["score"].to_numpy(dtype="float64")[own]
                metrics = C.score_unit(labels[own], found, C.unit_seed(task, REFERENCE, scheme, fold, unit))
                units.append(unit_row(task, REFERENCE, population, scheme, fold, unit, train, labels, 0, test, own, "scored", "", metrics))
                predictions.append(rows.loc[own, C.KEYS].assign(task=task, model=REFERENCE, scheme=scheme, fold=fold, unit=unit, label=labels[own], score=found))
        done = pd.DataFrame([row for row in units if row["task"] == task])
        print(f"[02] {task}: " + "; ".join(f"{model} {scheme} {part['auroc'].mean():.3f} over {len(part)} units" for (model, scheme), part in
                                          done[done["status"] == "scored"].groupby(["model", "scheme"], sort=False)), flush=True)

    units = pd.DataFrame(units, columns=UNIT_COLUMNS)
    predictions = pd.concat(predictions, ignore_index=True)[PREDICTION_COLUMNS]
    chosen = pd.concat(chosen, ignore_index=True)[GENE_COLUMNS]
    reference_scheme = C.CONFIG["metrics"]["inflation_reference"]
    summary = []
    for (task, model, population, scheme), part in units.groupby(["task", "model", "population", "scheme"], sort=False):
        scored = part[part["status"] == "scored"]
        summary.append(dict(task=task, model=model, population=population, scheme=scheme, units=len(part), units_scored=len(scored),
                            records_scored=int(scored["predicted"].sum()), auroc_mean=float(scored["auroc"].mean()), auroc_median=float(scored["auroc"].median()),
                            auprc_mean=float(scored["auprc"].mean()), auprc_median=float(scored["auprc"].median()),
                            prevalence_mean=float(scored["prevalence"].mean())))
    summary = pd.DataFrame(summary)
    base = summary[summary["scheme"] == reference_scheme].set_index(["task", "model"])
    for name in ("auroc", "auprc"):
        against = base[f"{name}_mean"].reindex(list(zip(summary["task"], summary["model"]))).to_numpy()
        summary[f"inflation_{name}"] = np.where(summary["population"] == "naive", summary[f"{name}_mean"].to_numpy() - against, np.nan)
    summary = summary[SUMMARY_COLUMNS]
    studies = study_scores(records, predictions, splits)
    found, matched, flags = inflation_effects(summary, studies, predictions)
    effects = pd.DataFrame(found, columns=EFFECT_COLUMNS)
    additive = pd.DataFrame(matched, columns=ADDITIVE_COLUMNS).sort_values(["task", "model"], kind="stable").reset_index(drop=True)
    keys = list(zip(studies["task"], studies["model"], studies["scheme"], studies["fold"], studies["study_accession"]))
    for name in FLAGS:
        studies[name] = [key in flags[name] for key in keys]
    studies = studies.sort_values(["task", "model", "scheme", "fold", "study_accession", "subset"], kind="stable").reset_index(drop=True)
    effects = effects.sort_values(["task", "model"], kind="stable").reset_index(drop=True)
    order = ["task", "model", "scheme", "fold", "unit"]
    units = units.sort_values(order, kind="stable").reset_index(drop=True)
    predictions = predictions.sort_values([*order, "study_accession", "sample_accession"], kind="stable").reset_index(drop=True)
    chosen = chosen.sort_values(["task", "scheme", "fold", "rank"], kind="stable").reset_index(drop=True)
    summary = summary.sort_values(["task", "model", "population", "scheme"], kind="stable").reset_index(drop=True)
    records = records.sort_values(["task", "study_accession", "sample_accession"], kind="stable").reset_index(drop=True)
    outputs = {name: C.DATA / f"{name}.parquet" for name in ("model_records", "model_units", "model_predictions", "model_genes", "model_summary",
                                                             "model_study_scores", "model_inflation_effects", "model_additive_matched")}
    for name, frame in zip(outputs, (records, units, predictions, chosen, summary, studies, effects, additive)):
        C.save_table(frame, outputs[name])
    C.save_review(summary, C.OUT / "model_summary_review.csv")
    C.save_review(effects, C.OUT / "model_inflation_effects_review.csv")
    sources, sensitivity = {}, {}
    for row in effects.itertuples(index=False):
        sources.setdefault(row.task, {}).setdefault(row.model, {})[row.effect] = dict(
            status=row.status, estimate=C.rounded(row.estimate), low=C.rounded(row.low), high=C.rounded(row.high), studies=int(row.studies),
            first_records=int(row.first_records), second_records=int(row.second_records))
    for row in additive.itertuples(index=False):
        sensitivity.setdefault(row.task, {})[row.model] = dict(status=row.status, loso_status=row.loso_status, matched_studies=int(row.matched_studies),
                                                               loso_matched_studies=int(row.loso_matched_studies),
                                                               **{name: C.rounded(getattr(row, name)) for name in VALUES})
    by_task = {task: dict(records=len(part), with_values=int(part["has_values"].sum()), samples_of_the_task=int(part["in_gts"].sum()),
                          samples_of_the_task_with_values=int((part["in_gts"] & part["has_values"]).sum()),
                          left_out_by_reason=C.counts(part.loc[~part["has_values"], "reason"])) for task, part in records.groupby("task", sort=False)}
    results = {}
    for row in summary.itertuples(index=False):
        results.setdefault(row.task, {}).setdefault(row.model, {})[row.scheme] = dict(
            units_scored=int(row.units_scored), auroc_mean=C.rounded(row.auroc_mean), auprc_mean=C.rounded(row.auprc_mean),
            inflation_auroc=C.rounded(row.inflation_auroc), inflation_auprc=C.rounded(row.inflation_auprc))
    C.finish(2, start, [*outputs.values(), C.OUT / "model_summary_review.csv", C.OUT / "model_inflation_effects_review.csv"], method=C.CONFIG["models"]["method"],
             features=dict(method=C.CONFIG["features"]["method"], group=C.CONFIG["features"]["group"], genes_in_group=len(genes),
                           records_with_values=len(wanted), genes_selected=C.CONFIG["features"]["genes"], by_task=by_task),
             reasons=REASONS, unit_reasons=UNIT_REASONS, folds_fitted=int(chosen.groupby(["task", "scheme", "fold"]).ngroups),
             units=dict(total=len(units), scored=int((units["status"] == "scored").sum()),
                        not_scored_by_reason=C.counts(units.loc[units["status"] != "scored", "reason"])),
             predictions=len(predictions), summary=results, inflation_effects=sources, additive_matched=sensitivity, effect_reasons=EFFECT_REASONS,
             effects={name: comparison for name, (comparison, _) in EFFECTS.items()},
             study_scores=dict(total=len(studies), scored=int((studies["status"] == "scored").sum()), **{name: int(studies[name].sum()) for name in FLAGS}),
             models={name: C.CONFIG["models"][name] for name in C.MODELS}, metrics=C.CONFIG["metrics"], protocol=C.CONFIG["protocol"])
    for (task, model), part in effects.groupby(["task", "model"], sort=False):
        found = part.set_index("effect")
        shown = [f"total inflation {found.at['total_inflation', 'estimate']:+.3f}"]
        for effect in list(EFFECTS)[1:]:
            entry = found.loc[effect]
            shown.append(f"{effect.replace('_', ' ')} {entry['estimate']:+.3f} ({entry['low']:+.3f} to {entry['high']:+.3f}, {int(entry['studies'])} studies)"
                         if entry["status"] == "estimable" else f"{effect.replace('_', ' ')} not estimable ({int(entry['studies'])} studies)")
        print(f"[02] {task} {model}: " + "; ".join(shown), flush=True)
    print("[02] Complete: " + "; ".join(f"{row.task} {row.model} inflation {row.inflation_auroc:+.3f} ({row.scheme})" for row in summary.itertuples(index=False)
                                        if row.population == "naive" and np.isfinite(row.inflation_auroc))
          + f"; review sheet: {(C.OUT / 'model_summary_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
