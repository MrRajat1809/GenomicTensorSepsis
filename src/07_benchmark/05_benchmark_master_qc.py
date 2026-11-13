"""Verify the leakage audit, the baseline models, the study shortcuts and the signature tests of the benchmark stage from their files."""
import hashlib
import importlib
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_audit_split_leakage")
M = importlib.import_module("02_train_baseline_models")
S = importlib.import_module("03_test_study_shortcuts")
G = importlib.import_module("04_retest_published_signatures")

RECORD = ["study_accession", "sample_accession"]
MEANING = ("PASS means that every population of the leakage audit is the set of records its rule gives, with its links to other deposits of an array and "
           "other samples of an identified patient counted again by joining the cohort table with itself; that the naive splits, drawn again from their "
           "seeds, expose the counted test records, and that no fold of the cohort stage's patient-safe schemes holds a linked record on both sides; that "
           "every model predicted exactly the records of its test unit, that no record of a patient-safe test unit shares an array or an identified "
           "patient with a record of its training groups, that every metric and interval follows from the stored predictions and seeds, and that one fold "
           "of each model, fitted again, gives the stored genes and predictions; that every within-study AUROC, over all test records and over those "
           "without a twin among their fold's training records found again by joining the records with themselves, follows from the stored predictions, "
           "that the study shortcut, the twins, the seen study and the twins across studies follow from their rebuilt sets with intervals replayed from "
           "their seeds, and that the additive split over one matched set of studies adds up; that the study-only predictions and the study identification follow from "
           "their populations and folds, with every held-out study scoring one half; and that every signature stratum was tested or skipped by the first "
           "reason that applies, with the scores of a sample of strata computed again by gene symbol, every AUROC counted again from the scores, every "
           "stratum marked independent by the validation stage's rule, and the Sepsis MetaScore and its independent strata equal to those of the "
           "validation stage. It verifies that the benchmark was carried out as described; what its numbers show is in the tables and review sheets.")


class QC:
    def __init__(self):
        self.checks = []
        self.section = ""

    def check(self, name, passed, detail=None):
        self.checks.append(dict(section=self.section, name=name, passed=bool(passed), detail=detail))
        if not passed:
            raise ValueError(f"QC failed: {name}")

    def phase(self, name, action):
        self.section = name
        print(f"[05] {name}", flush=True)
        action()


def close(found, stored, tolerance=None):
    """Whether recomputed numbers equal stored ones within the metric tolerance, missing exactly where the stored ones are missing."""
    found, stored = np.asarray(found, dtype="float64"), np.asarray(stored, dtype="float64")
    tolerance = C.CONFIG["quality_control"]["metric_tolerance"] if tolerance is None else tolerance
    same = np.isnan(found) == np.isnan(stored)
    return bool(same.all() and (np.abs(np.where(np.isnan(found), 0, found) - np.where(np.isnan(stored), 0, stored)) <= tolerance).all())


def pairwise(cases, others):
    """AUROC counted pair by pair, ties as half: a second route to the area under the curve."""
    cases, others = np.asarray(cases, dtype="float64"), np.asarray(others, dtype="float64")
    return float(((cases[:, None] > others[None, :]).sum() + 0.5 * (cases[:, None] == others[None, :]).sum()) / (len(cases) * len(others)))


def link_table(frame):
    """Every ordered pair of linked records of a frame, found by joining it with itself on the array and on the identified patient."""
    keyed = frame[["study_accession", "sample_id", "patient_key", "patient_identified"]].assign(row=np.arange(len(frame)))
    arrays = keyed.merge(keyed, on="sample_id", suffixes=("", "_other"))
    arrays = arrays[arrays["row"] != arrays["row_other"]].assign(kind="same_array")
    known = keyed[keyed["patient_identified"]]
    patients = known.merge(known, on="patient_key", suffixes=("", "_other"))
    patients = patients[patients["sample_id"] != patients["sample_id_other"]].assign(kind="same_patient")
    columns = ["row", "row_other", "study_accession", "study_accession_other", "kind"]
    return pd.concat([arrays[columns], patients[columns]], ignore_index=True)


def exposure(links, test, train):
    hit = links[test[links["row"].to_numpy()] & train[links["row_other"].to_numpy()]]
    return (hit["row"].nunique(), hit.loc[hit["kind"] == "same_array", "row"].nunique(), hit.loc[hit["kind"] == "same_patient", "row"].nunique())


def random_test(frame, repeat):
    """The test records of one repeat of the naive random split, drawn again with Python's sort of the digests."""
    policy = C.CONFIG["splits"]
    digests = sorted((hashlib.sha256(f"{policy['seed']}:{repeat}:{study}:{sample}".encode("utf-8")).hexdigest(), study, sample)
                     for study, sample in zip(frame["study_accession"], frame["sample_accession"]))
    count = int(policy["test_fraction"] * len(frame) + 0.5)
    chosen = {(study, sample) for _, study, sample in digests[:count]}
    return np.array([key in chosen for key in zip(frame["study_accession"], frame["sample_accession"])])


def provenance_checks(qc, state):
    state["start"] = C.begin(5)
    qc.check("all_stage_receipts_and_code_bindings", True)
    state["policy"] = C.cohort_policy()
    state["table"] = pd.read_parquet(C.COHORT / "cohort_samples.parquet")
    state["tasks"] = pd.read_parquet(C.COHORT / "benchmark_tasks.parquet")
    state["splits"] = pd.read_parquet(C.COHORT / "benchmark_splits.parquet")
    state["summary"] = {}


def leakage_checks(qc, state):
    policy, table, tasks, splits = state["policy"], state["table"], state["tasks"], state["splits"]
    populations = pd.read_parquet(C.DATA / "split_populations.parquet")
    leakage = pd.read_parquet(C.DATA / "split_leakage.parquet")
    qc.check("tables_have_their_columns", list(populations.columns) == C.POPULATION_COLUMNS and list(leakage.columns) == C.LEAKAGE_COLUMNS)
    # The populations are rebuilt here with one selection each, apart from the helpers of stage 01.
    expected = {("all_records", ""): (set(zip(table["study_accession"], table["sample_accession"])), None)}
    for task in policy["tasks"]:
        value = table[task["variable"]]
        chosen = table["tissue"].isin(policy["tissues"]) & value.isin([task["positive"], task["negative"]])
        chosen &= (table["sample_class"] == task["within_class"]) if "within_class" in task else True
        expected[("naive", task["task"])] = (set(zip(table.loc[chosen, "study_accession"], table.loc[chosen, "sample_accession"])),
                                             dict(zip(zip(table.loc[chosen, "study_accession"], table.loc[chosen, "sample_accession"]),
                                                      (value[chosen] == task["positive"]).astype(int))))
        safe = tasks[(tasks["task"] == task["task"]) & tasks["eligible"]]
        expected[("gts", task["task"])] = (set(zip(safe["study_accession"], safe["sample_accession"])),
                                           dict(zip(zip(safe["study_accession"], safe["sample_accession"]), safe["label"].astype(int))))
    right, links_right, frames = True, True, {}
    for (population, task), (keys, labels) in expected.items():
        part = populations[(populations["population"] == population) & (populations["task"] == task)].reset_index(drop=True)
        found = list(zip(part["study_accession"], part["sample_accession"]))
        right = right and len(found) == len(keys) and set(found) == keys and (labels is None and part["label"].isna().all()
                                                                               or labels is not None and [labels[key] for key in found] == list(part["label"]))
        links = link_table(part)
        arrays, patients = (np.bincount(links.loc[links["kind"] == kind, "row"].to_numpy(dtype="int64"), minlength=len(part)) for kind in ("same_array", "same_patient"))
        other = links[links["study_accession"] != links["study_accession_other"]]["row"].unique()
        links_right = (links_right and (arrays == part["same_array_links"].to_numpy()).all() and (patients == part["same_patient_links"].to_numpy()).all()
                       and (arrays + patients == part["linked_records"].to_numpy()).all()
                       and (np.isin(np.arange(len(part)), other) == part["linked_in_other_study"].to_numpy()).all())
        frames[(population, task)] = (part, links)
    qc.check("populations_are_the_records_their_rules_give", right, {f"{population}:{task}" if task else population: len(keys)
                                                                   for (population, task), (keys, _) in expected.items()})
    qc.check("links_counted_again_by_joining_the_table_with_itself", links_right)
    # The naive splits are drawn again and their exposures counted again from the joined links.
    exposures, share = True, C.CONFIG["splits"]["test_fraction"]
    for (population, task), (part, links) in frames.items():
        if population == "gts":
            continue
        stored = leakage[(leakage["population"] == population) & (leakage["task"] == task)].set_index(["scheme", "fold"])
        studies = part["study_accession"].to_numpy()
        folds = [("naive_record_split", C.repeat_name(repeat), random_test(part, repeat)) for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1)]
        folds += [("naive_leave_one_study_out", study, studies == study) for study in sorted(set(studies))]
        expected_share = float(np.mean(1 - share ** part["linked_records"].to_numpy()))
        for scheme, fold, test in folds:
            row = stored.loc[(scheme, fold)]
            counts = exposure(links, test, ~test)
            exposures = (exposures and row["test_records"] == test.sum() and row["train_records"] == (~test).sum()
                         and (row["exposed"], row["exposed_same_array"], row["exposed_same_patient"]) == counts
                         and (scheme != "naive_record_split" or abs(row["expected_share"] - expected_share) < 1e-12))
        exposures = exposures and len(stored) == len(folds)
    qc.check("naive_splits_drawn_again_expose_the_counted_records", exposures)
    # Every record of the cohort, not only the samples of a task, is placed by its group, so a patient's other samples count as well.
    whole, safe_right = link_table(table), True
    groups = table["patient_safe_group"].to_numpy()
    for (task, scheme, fold), part in splits.groupby(["task", "scheme", "fold"], sort=False):
        train, held = part.loc[part["role"] == "train", "patient_safe_group"], part.loc[part["role"] != "train", "patient_safe_group"]
        crossing = exposure(whole, np.isin(groups, held), np.isin(groups, train))[0]
        row = leakage[(leakage["population"] == "gts") & (leakage["task"] == task) & (leakage["scheme"] == scheme) & (leakage["fold"] == fold)]
        safe_right = safe_right and crossing == 0 and len(row) == 1 and int(row["exposed"].iloc[0]) == 0 and int(row["linked_cohort_records"].iloc[0]) == 0
    qc.check("no_patient_safe_fold_links_records_across_its_split", safe_right and int((leakage["population"] == "gts").sum()) == splits.groupby(
        ["task", "scheme", "fold"]).ngroups, int((leakage["population"] == "gts").sum()))
    state["populations"] = populations
    every = leakage[(leakage["population"] == "all_records")]
    state["summary"].update(records=len(table), share_exposed_random_split=round(float(every.loc[every["scheme"] == "naive_record_split", "exposed"].sum()
                                                                                    / every.loc[every["scheme"] == "naive_record_split", "test_records"].sum()), 6),
                            share_exposed_leave_one_study_out=round(float(every.loc[every["scheme"] == "naive_leave_one_study_out", "exposed"].sum() / len(table)), 6),
                            patient_safe_folds=int((leakage["population"] == "gts").sum()))


def model_checks(qc, state):
    populations, splits, table = state["populations"], state["splits"], state["table"]
    records = pd.read_parquet(C.DATA / "model_records.parquet")
    units = pd.read_parquet(C.DATA / "model_units.parquet")
    predictions = pd.read_parquet(C.DATA / "model_predictions.parquet")
    chosen = pd.read_parquet(C.DATA / "model_genes.parquet")
    summary = pd.read_parquet(C.DATA / "model_summary.parquet")
    qc.check("model_tables_have_their_columns", list(records.columns) == M.RECORD_COLUMNS and list(units.columns) == M.UNIT_COLUMNS
             and list(predictions.columns) == M.PREDICTION_COLUMNS and list(chosen.columns) == M.GENE_COLUMNS and list(summary.columns) == M.SUMMARY_COLUMNS)
    correction = pd.read_parquet(C.CORRECTION / "correction_records.parquet").set_index(RECORD)
    group = C.CONFIG["features"]["group"]
    keyed = list(zip(records["study_accession"], records["sample_accession"]))
    has = ((correction.loc[keyed, "group"] == group) & (correction.loc[keyed, "status"] == "in_unit")).to_numpy()
    naive = populations[populations["population"] == "naive"]
    safe = set(zip(*(populations.loc[populations["population"] == "gts", name] for name in ("task", "study_accession", "sample_accession"))))
    qc.check("model_records_are_the_naive_populations_with_their_values", len(records) == len(naive)
             and set(zip(records["task"], records["study_accession"], records["sample_accession"])) == set(zip(naive["task"], naive["study_accession"], naive["sample_accession"]))
             and (records["has_values"].to_numpy() == has).all()
             and (records["in_gts"].to_numpy() == np.array([key in safe for key in zip(records["task"], records["study_accession"], records["sample_accession"])])).all()
             and (records.loc[~records["has_values"], "reason"].to_numpy() == correction.loc[[key for key, held in zip(keyed, has) if not held], "reason"].replace(
                 "", "sample_type_not_corrected").to_numpy()).all() and (records.loc[records["has_values"], "reason"] == "").all(),
             C.counts(records.loc[~records["has_values"], "reason"]))
    # The test records of every unit are rebuilt from the rules of its scheme, apart from the fold generator of stage 02.
    indexed = table.set_index(RECORD)
    with_score = set(map(tuple, pd.read_parquet(C.VALIDATION / "sepsis_score_samples.parquet", columns=["study_accession", "sample_accession", "platform"]).to_numpy()))
    integrity, crossing, same_sets = True, [], True
    for task, rows in records.groupby("task", sort=False):
        rows = rows.reset_index(drop=True)
        held, inside = rows["has_values"].to_numpy(), rows["in_gts"].to_numpy()
        has_score = np.array([key in with_score for key in zip(rows["study_accession"], rows["sample_accession"], rows["platform"])])
        expected = {("naive_record_split", C.repeat_name(repeat), "pooled"): random_test(rows, repeat) & held for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1)}
        expected.update({("naive_leave_one_study_out", study, study): (rows["study_accession"] == study).to_numpy() & held
                          for study in sorted(set(rows["study_accession"]))})
        trains = {}
        for (scheme, fold), part in splits[splits["task"] == task].groupby(["scheme", "fold"], sort=False):
            trains[(scheme, fold)] = set(part.loc[part["role"] == "train", "patient_safe_group"])
            for name in part.loc[part["role"] == "test", "patient_safe_group"]:
                expected[(scheme, fold, name)] = (rows["patient_safe_group"] == name).to_numpy() & inside & held
        own_units = units[units["task"] == task]
        listed = set(zip(own_units["scheme"], own_units["fold"], own_units["unit"]))
        integrity = integrity and listed == set(expected)
        labels = rows["label"].to_numpy()
        record_of = {key: index for index, key in enumerate(zip(rows["study_accession"], rows["sample_accession"]))}
        for (scheme, fold, unit), mask in expected.items():
            states = own_units[(own_units["scheme"] == scheme) & (own_units["fold"] == fold) & (own_units["unit"] == unit)]
            models, ready = states[states["model"] != M.REFERENCE], C.scorable(labels[mask])
            integrity = (integrity and len(models) == len(C.MODELS) and (models["status"] == "scored").all() == ready and (states["records"] == mask.sum()).all()
                         and len(states) == len(C.MODELS) + (task == M.REFERENCE_TASK))
            if task == M.REFERENCE_TASK:
                own = states[states["model"] == M.REFERENCE].iloc[0]
                integrity = integrity and (own["status"] == "scored") == (ready and C.scorable(labels[mask & has_score])) and (
                    own["status"] == "scored" or own["reason"] == ("fewer_than_minimum_of_a_class" if not ready else "fewer_than_minimum_with_a_reference_score"))
            if not ready:
                integrity = integrity and (states["reason"] == "fewer_than_minimum_of_a_class").all()
                continue
            wanted = set(zip(rows.loc[mask, "study_accession"], rows.loc[mask, "sample_accession"]))
            found = predictions[(predictions["task"] == task) & (predictions["scheme"] == scheme) & (predictions["fold"] == fold) & (predictions["unit"] == unit)]
            for model in C.MODELS:
                got = found[found["model"] == model]
                same_sets = same_sets and set(zip(got["study_accession"], got["sample_accession"])) == wanted and len(got) == len(wanted)
            if scheme in C.GTS_SCHEMES:
                test_keys = list(wanted)
                train_records = table[table["patient_safe_group"].isin(trains[(scheme, fold)])]
                test_records = indexed.loc[test_keys]
                known = test_records["patient_key"][test_records["patient_key"] != ""]
                if (set(train_records["sample_id"]) & set(test_records["sample_id"]) or set(train_records.loc[train_records["patient_key"] != "", "patient_key"]) & set(known)
                        or set(test_records["patient_safe_group"]) & trains[(scheme, fold)]):
                    crossing.append(f"{task} {scheme} {fold} {unit}")
            integrity = integrity and all(labels[record_of[key]] == label for key, label in zip(zip(found["study_accession"], found["sample_accession"]), found["label"]))
    qc.check("every_unit_holds_the_records_of_its_scheme", integrity, len(units))
    qc.check("both_models_predict_exactly_the_records_of_every_scored_unit", same_sets)
    qc.check("no_patient_safe_test_record_shares_an_array_or_patient_with_training", not crossing, crossing[:20])
    # Every metric is computed again from the stored predictions with scikit-learn, and every interval drawn again from its seed.
    from sklearn.metrics import average_precision_score, roc_auc_score

    # Summed in floating point, the average precision of a perfect unit can exceed one by an ulp.
    scored, tolerance = units[units["status"] == "scored"].reset_index(drop=True), C.CONFIG["quality_control"]["metric_tolerance"]
    grouped = predictions.groupby(["task", "model", "scheme", "fold", "unit"], sort=False)
    again, counted = [], True
    for row in scored.itertuples(index=False):
        part = grouped.get_group((row.task, row.model, row.scheme, row.fold, row.unit))
        labels, scores = part["label"].to_numpy(dtype="int64"), part["score"].to_numpy(dtype="float64")
        replay = C.score_unit(labels, scores, C.unit_seed(row.task, row.model, row.scheme, row.fold, row.unit))
        again.append([roc_auc_score(labels, scores), average_precision_score(labels, scores), *(replay[name] for name in ("auroc_low", "auroc_high", "auprc_low",
                                                                                                                          "auprc_high", "prevalence"))])
        counted = counted and row.predicted == len(part) and row.positives == int(labels.sum()) and row.negatives == int((labels == 0).sum())
    again = np.array(again)
    qc.check("metrics_and_intervals_follow_from_the_predictions_and_seeds", counted and len(again) == len(scored) and grouped.ngroups == len(scored)
             and all(close(again[:, index], scored[name]) for index, name in enumerate(("auroc", "auprc", "auroc_low", "auroc_high", "auprc_low", "auprc_high",
                                                                                         "prevalence")))
             and (scored["auroc_low"] <= scored["auroc_high"]).all() and (scored["auprc_low"] <= scored["auprc_high"]).all()
             and scored[list(M.METRICS)].stack().between(-tolerance, 1 + tolerance).all(), len(scored))
    reference = pd.read_parquet(C.VALIDATION / "sepsis_score_samples.parquet").set_index(["study_accession", "sample_accession", "platform"])["score"]
    own = predictions[predictions["model"] == M.REFERENCE]
    qc.check("reference_scores_are_those_of_the_validation_stage", (own["task"] == M.REFERENCE_TASK).all()
             and close(reference.loc[list(zip(own["study_accession"], own["sample_accession"], own["platform"]))].to_numpy(), own["score"]), len(own))
    rebuilt, base = [], C.CONFIG["metrics"]["inflation_reference"]
    for (task, model, population, scheme), part in units.groupby(["task", "model", "population", "scheme"], sort=False):
        good = part[part["status"] == "scored"]
        rebuilt.append((task, model, scheme, len(part), len(good), good["auroc"].mean(), good["auprc"].mean()))
    rebuilt = pd.DataFrame(rebuilt, columns=["task", "model", "scheme", "units", "units_scored", "auroc_mean", "auprc_mean"]).set_index(["task", "model", "scheme"])
    keyed = summary.set_index(["task", "model", "scheme"])
    inflation = keyed["auroc_mean"] - [rebuilt.loc[(task, model, base), "auroc_mean"] for task, model, _ in keyed.index]
    qc.check("summary_and_inflation_follow_from_the_units", len(keyed) == len(rebuilt) and close(rebuilt.loc[keyed.index, "auroc_mean"], keyed["auroc_mean"])
             and close(rebuilt.loc[keyed.index, "auprc_mean"], keyed["auprc_mean"]) and (rebuilt.loc[keyed.index, "units_scored"].to_numpy() == keyed["units_scored"].to_numpy()).all()
             and close(np.where(keyed["population"] == "naive", inflation, np.nan), keyed["inflation_auroc"]))
    fitted = chosen.groupby(["task", "scheme", "fold"])
    qc.check("every_fitted_fold_keeps_its_genes_from_its_training_records", (fitted.size() == C.CONFIG["features"]["genes"]).all()
             and (chosen["present_share"] >= C.CONFIG["features"]["minimum_present_share"]).all() and not chosen.duplicated(["task", "scheme", "fold", "hgnc_id"]).any()
             and set(fitted.groups) == set(zip(scored.loc[scored["model"] != M.REFERENCE, "task"], scored.loc[scored["model"] != M.REFERENCE, "scheme"],
                                                scored.loc[scored["model"] != M.REFERENCE, "fold"])), fitted.ngroups)
    # One fold of each model is fitted again: its genes from its training records alone, and its predictions.
    refit = C.CONFIG["quality_control"]["refit"]
    task_rows = records[records["task"] == refit["task"]].reset_index(drop=True)
    fold = sorted(scored.loc[(scored["task"] == refit["task"]) & (scored["scheme"] == refit["scheme"]), "fold"])[0]
    _, scheme, name, train, tests = next(entry for entry in M.task_folds(task_rows, splits, refit["task"]) if entry[1] == refit["scheme"] and entry[2] == fold)
    trains = set(splits.loc[(splits["task"] == refit["task"]) & (splits["scheme"] == scheme) & (splits["fold"] == fold) & (splits["role"] == "train"), "patient_safe_group"])
    qc.check("refit_fold_trains_on_the_samples_of_its_training_groups", (train == (task_rows["patient_safe_group"].isin(trains) & task_rows["in_gts"]
                                                                                  & task_rows["has_values"]).to_numpy()).all())
    wanted = task_rows[task_rows["has_values"]].reset_index(drop=True)
    values, genes = M.feature_values(wanted)
    place = pd.MultiIndex.from_frame(wanted[C.KEYS]).get_indexer(pd.MultiIndex.from_frame(task_rows[C.KEYS]))
    labels = task_rows["label"].to_numpy(dtype="int64")
    columns, statistic, _ = M.select_genes(values, place[train], labels[train])
    stored = chosen[(chosen["task"] == refit["task"]) & (chosen["scheme"] == scheme) & (chosen["fold"] == fold)].sort_values("rank")
    ready = [(unit, test) for unit, test in tests if C.scorable(labels[test])]
    test_rows = np.concatenate([np.flatnonzero(test) for _, test in ready])
    worst = 0.0
    for model in C.MODELS:
        found = M.fit_predict(model, values[np.ix_(place[train], columns)], labels[train], values[np.ix_(place[test_rows], columns)])
        kept = predictions[(predictions["task"] == refit["task"]) & (predictions["model"] == model) & (predictions["scheme"] == scheme) & (predictions["fold"] == fold)]
        kept = kept.set_index(RECORD)["score"].loc[list(zip(task_rows.loc[test_rows, "study_accession"], task_rows.loc[test_rows, "sample_accession"]))]
        worst = max(worst, float(np.abs(found - kept.to_numpy()).max()))
    qc.check("one_fold_of_each_model_fitted_again_gives_the_stored_genes_and_predictions", [genes[column] for column in columns] == list(stored["hgnc_id"])
             and close(statistic, stored["statistic"], 1e-9) and worst <= C.CONFIG["quality_control"]["refit_tolerance"], dict(fold=fold, largest_difference=worst))
    state["summary"].update(model_units_scored=len(scored), inflation={f"{row.task}:{row.model}:{row.scheme}": C.rounded(row.inflation_auroc)
                                                                         for row in summary.itertuples(index=False) if row.population == "naive"})


def inflation_checks(qc, state):
    splits = state["splits"]
    records = pd.read_parquet(C.DATA / "model_records.parquet")
    predictions = pd.read_parquet(C.DATA / "model_predictions.parquet")
    units = pd.read_parquet(C.DATA / "model_units.parquet")
    studies = pd.read_parquet(C.DATA / "model_study_scores.parquet")
    effects = pd.read_parquet(C.DATA / "model_inflation_effects.parquet")
    additive = pd.read_parquet(C.DATA / "model_additive_matched.parquet")
    summary = pd.read_parquet(C.DATA / "model_summary.parquet").set_index(["task", "model", "scheme"])
    qc.check("inflation_tables_have_their_columns", list(studies.columns) == M.STUDY_COLUMNS and list(effects.columns) == M.EFFECT_COLUMNS
             and list(additive.columns) == M.ADDITIVE_COLUMNS)
    # The twin status of every test record is decided again by joining the task's records with themselves, against the training records of its own fold,
    # which are rebuilt from the rules of its scheme; every within-study AUROC is then counted again, pair by pair, from the stored predictions.
    grouped, expected, minimum = predictions.groupby(["task", "scheme", "fold"]), {}, C.CONFIG["splits"]["minimum_per_class"]
    reference = C.CONFIG["metrics"]["inflation_reference"]
    for task, rows in records.groupby("task", sort=False):
        rows = rows.reset_index(drop=True)
        links, held, inside = link_table(rows), rows["has_values"].to_numpy(), rows["in_gts"].to_numpy()
        studies_of, groups = rows["study_accession"].to_numpy(), rows["patient_safe_group"].to_numpy()
        trains = {("naive_record_split", C.repeat_name(repeat)): ~random_test(rows, repeat) & held for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1)}
        trains.update({("naive_leave_one_study_out", study): (studies_of != study) & held for study in sorted(set(studies_of))})
        for fold, part in splits[(splits["task"] == task) & (splits["scheme"] == reference)].groupby("fold"):
            trains[(reference, fold)] = np.isin(groups, part.loc[part["role"] == "train", "patient_safe_group"].to_numpy()) & inside & held
        position = {key: index for index, key in enumerate(zip(rows["study_accession"], rows["sample_accession"]))}
        for (scheme, fold), train in trains.items():
            if (task, scheme, fold) not in grouped.groups:
                continue
            found = grouped.get_group((task, scheme, fold))
            place = np.array([position[key] for key in zip(found["study_accession"], found["sample_accession"])])
            hit = links[train[links["row_other"].to_numpy()]]["row"].unique()
            twin = np.isin(place, hit)
            for (model, unit, study), own in found.assign(twin=twin).groupby(["model", "unit", "study_accession"]):
                for subset in (("all", "without_twins") if scheme in C.NAIVE_SCHEMES else ("all",)):
                    chosen = own if subset == "all" else own[~own["twin"]]
                    cases, others = chosen.loc[chosen["label"] == 1, "score"], chosen.loc[chosen["label"] == 0, "score"]
                    ready = len(cases) >= minimum and len(others) >= minimum
                    expected[(task, model, scheme, fold, unit, study, subset)] = (len(chosen), len(cases), len(others), "scored" if ready else "not_scored",
                                                                                 pairwise(cases, others) if ready else np.nan)
    keyed = studies.set_index(["task", "model", "scheme", "fold", "unit", "study_accession", "subset"])
    follows = len(keyed) == len(expected) and set(keyed.index) == set(expected) and all(
        (keyed.at[key, "records"], keyed.at[key, "positives"], keyed.at[key, "negatives"], keyed.at[key, "status"]) == value[:4]
        and close([keyed.at[key, "auroc"]], [value[4]]) for key, value in expected.items())
    qc.check("within_study_aurocs_and_twins_follow_from_the_predictions_and_each_folds_training_records", follows, len(expected))
    # The sets of the effects are rebuilt from the study rows as tables of study by repeat, and every estimate is computed again.
    scored, reference, policy = studies[studies["status"] == "scored"], C.CONFIG["metrics"]["inflation_reference"], C.CONFIG["metrics"]["inflation_effects"]
    least, flags, estimates, intervals, replays = policy["minimum_studies"], {name: set() for name in M.FLAGS}, True, True, 0
    keyed, stored = summary, effects.set_index(["task", "model", "effect"])
    tail = (1 - policy["bootstrap"]["level"]) / 2

    def bounds(values):
        values = np.asarray(values, dtype="float64")
        values = values[np.isfinite(values)]
        return [float(np.quantile(values, tail)), float(np.quantile(values, 1 - tail))]

    def rank_auroc(labels, scores):
        from scipy.stats import rankdata

        ranks = rankdata(scores)
        cases = int((labels == 1).sum())
        return float((ranks[labels == 1].sum() - cases * (cases + 1) / 2) / (cases * (len(labels) - cases)))

    for (task, model), own in scored[scored["model"].isin(C.MODELS)].groupby(["task", "model"]):
        split = own[own["scheme"] == "naive_record_split"]
        table = split.pivot_table(index=["study_accession", "fold"], columns="subset", values="auroc").reindex(columns=list(M.SUBSETS))
        sizes = split.pivot_table(index=["study_accession", "fold"], columns="subset", values="records", aggfunc="sum").reindex(columns=list(M.SUBSETS))
        paired = table.dropna()
        safe = own[(own["scheme"] == reference) & (own["subset"] == "all")]
        held = safe.groupby("study_accession").agg(auroc=("auroc", "mean"), records=("records", "sum"))
        without = table["without_twins"].dropna()
        seen = sorted(set(without.index.get_level_values("study_accession")) & set(held.index))
        loso = own[own["scheme"] == "naive_leave_one_study_out"]
        across = loso.pivot_table(index="fold", columns="subset", values="auroc").reindex(columns=list(M.SUBSETS)).dropna()
        counts = loso.pivot_table(index="fold", columns="subset", values="records", aggfunc="sum").reindex(columns=list(M.SUBSETS))
        flags["in_twins_set"] |= {(task, model, "naive_record_split", fold, study) for study, fold in paired.index}
        chosen = without[without.index.get_level_values("study_accession").isin(seen)]
        flags["in_seen_study_set"] |= {(task, model, "naive_record_split", fold, study) for study, fold in chosen.index}
        flags["in_seen_study_set"] |= {(task, model, reference, fold, study) for fold, study in zip(safe["fold"], safe["study_accession"]) if study in seen}
        flags["in_twins_across_set"] |= {(task, model, "naive_leave_one_study_out", fold, fold) for fold in across.index}
        first_of = paired["all"].groupby(level="study_accession").mean()
        sets = {"twins": (first_of, paired["without_twins"].groupby(level="study_accession").mean(), int(sizes.loc[paired.index, "all"].sum()),
                          int(sizes.loc[paired.index, "without_twins"].sum())),
                "seen_study": (chosen.groupby(level="study_accession").mean(), held.loc[seen, "auroc"], int(sizes.loc[chosen.index, "without_twins"].sum()),
                               int(held.loc[seen, "records"].sum())),
                "twins_across_studies": (across["all"], across["without_twins"], int(counts.loc[across.index, "all"].sum()),
                                         int(counts.loc[across.index, "without_twins"].sum()))}
        for effect, (first, second, first_records, second_records) in sets.items():
            row = stored.loc[(task, model, effect)]
            if len(first) < least:
                estimates = estimates and row["status"] == "not_estimable" and int(row["studies"]) == len(first)
                continue
            first, second = first.sort_index(), second.sort_index()
            gaps = first.to_numpy() - second.loc[first.index].to_numpy()
            estimates = (estimates and row["status"] == "estimable" and int(row["studies"]) == len(first) and close([row["estimate"]], [gaps.mean()])
                         and close([row["first"], row["second"]], [first.mean(), second.mean()])
                         and (int(row["first_records"]), int(row["second_records"])) == (first_records, second_records))
            generator = np.random.default_rng(M.effect_seed(task, model, effect))
            draws = generator.integers(0, len(gaps), (policy["bootstrap"]["replicates"], len(gaps)))
            intervals = intervals and close(bounds([gaps[draw].mean() for draw in draws]), [row["low"], row["high"]])
            replays += 1
        # The study shortcut: the pooled AUROC of every repeat counted again by ranks, and its bootstrap replayed by repeating the records of every drawn study.
        pooled = predictions[(predictions["task"] == task) & (predictions["model"] == model) & (predictions["scheme"] == "naive_record_split")]
        names = sorted(set(pooled["study_accession"]))
        repeats = sorted(set(pooled["fold"]))
        within = table["all"].dropna()
        lookup = within.to_dict()
        parts = [pooled[pooled["fold"] == repeat] for repeat in repeats]
        cases = [(part["label"].to_numpy(dtype="int64"), part["score"].to_numpy(dtype="float64"), part["study_accession"].map({name: index for index, name in
                                                                                                                            enumerate(names)}).to_numpy())
                 for part in parts]

        def shortcut(counts):
            pooled_values, within_values = [], []
            for repeat, (labels, scores, owner) in zip(repeats, cases):
                taken = np.repeat(np.arange(len(labels)), counts[owner])
                if min(int((labels[taken] == 1).sum()), int((labels[taken] == 0).sum())) >= minimum:
                    pooled_values.append(rank_auroc(labels[taken], scores[taken]))
                drawn = [lookup.get((names[index], repeat)) for index in np.repeat(np.arange(len(names)), counts)]
                drawn = [value for value in drawn if value is not None]
                if drawn:
                    within_values.append(float(np.mean(drawn)))
            return (float(np.mean(pooled_values)) if pooled_values else np.nan) - (float(np.mean(within_values)) if within_values else np.nan), pooled_values, within_values

        row = stored.loc[(task, model, "study_shortcut")]
        value, first, second = shortcut(np.ones(len(names), dtype="int64"))
        estimates = (estimates and row["status"] == "estimable" and int(row["studies"]) == len(names) and close([row["estimate"]], [value])
                     and close([row["first"]], [float(np.mean(first))]) and close([row["first"]], [keyed.loc[(task, model, "naive_record_split"), "auroc_mean"]])
                     and close([row["second"]], [float(np.mean(second))]) and int(row["first_records"]) == len(pooled)
                     and int(row["second_records"]) == int(sizes.loc[within.index, "all"].sum()))
        generator = np.random.default_rng(M.effect_seed(task, model, "study_shortcut"))
        draws = generator.integers(0, len(names), (policy["bootstrap"]["replicates"], len(names)))
        intervals = intervals and close(bounds([shortcut(np.bincount(draw, minlength=len(names)))[0] for draw in draws]), [row["low"], row["high"]])
        replays += 1
        total, base = stored.loc[(task, model, "total_inflation")], keyed.loc[(task, model, reference)]
        estimates = (estimates and total["status"] == "reported" and close([total["estimate"]], [keyed.loc[(task, model, "naive_record_split"), "inflation_auroc"]])
                     and close([total["first"], total["second"]], [keyed.loc[(task, model, "naive_record_split"), "auroc_mean"], base["auroc_mean"]])
                     and int(total["studies"]) == int(base["units_scored"]))
    marked = {name: set(zip(*(studies.loc[studies[name], column] for column in ("task", "model", "scheme", "fold", "study_accession")))) for name in M.FLAGS}
    qc.check("effect_sets_rebuilt_from_the_study_rows", all(marked[name] == flags[name] for name in M.FLAGS[1:]) and len(effects) == len(M.EFFECTS) * len(
        set(zip(effects["task"], effects["model"]))), {name: len(flags[name]) for name in M.FLAGS[1:]})
    qc.check("effect_estimates_and_records_follow_from_the_study_rows_and_predictions", estimates)
    qc.check("effect_intervals_replayed_from_their_seeds", intervals, replays)
    # The additive sensitivity analysis over one matched set of studies S: the set is rebuilt and every part computed again.
    additive_flags, matched = set(), True
    for row in additive.itertuples(index=False):
        own = scored[(scored["task"] == row.task) & (scored["model"] == row.model)]
        split = own[own["scheme"] == "naive_record_split"]
        scores = split.pivot_table(index=["study_accession", "fold"], columns="subset", values="auroc").reindex(columns=list(M.SUBSETS)).dropna()
        counts = split.pivot_table(index=["study_accession", "fold"], columns="subset", values="records", aggfunc="sum").reindex(columns=list(M.SUBSETS))
        safe = own[(own["scheme"] == reference) & (own["subset"] == "all")]
        held = safe.groupby("study_accession").agg(auroc=("auroc", "mean"), records=("records", "sum"))
        members = sorted(set(scores.index.get_level_values("study_accession")) & set(held.index))
        loso = own[own["scheme"] == "naive_leave_one_study_out"]
        across = loso.pivot_table(index="fold", columns="subset", values="auroc").reindex(columns=list(M.SUBSETS)).dropna()
        sizes = loso.pivot_table(index="fold", columns="subset", values="records", aggfunc="sum").reindex(columns=list(M.SUBSETS))
        values, records = [np.nan] * len(M.VALUES), [0] * 6
        if len(members) >= least:
            kept = scores[scores.index.get_level_values("study_accession").isin(members)]
            pooled = predictions[(predictions["task"] == row.task) & (predictions["model"] == row.model) & (predictions["scheme"] == "naive_record_split")
                                 & predictions["study_accession"].isin(members)]
            repeats = [(pairwise(part.loc[part["label"] == 1, "score"], part.loc[part["label"] == 0, "score"]), len(part)) for _, part in pooled.groupby("fold")
                       if min(int((part["label"] == 1).sum()), int((part["label"] == 0).sum())) >= minimum]
            a = float(np.mean([value for value, _ in repeats]))
            b, c = (float(kept[subset].groupby(level="study_accession").mean().mean()) for subset in M.SUBSETS)
            d = float(held.loc[members, "auroc"].mean())
            values[:8] = [a, b, c, d, a - b, b - c, c - d, a - d]
            records[:4] = [sum(size for _, size in repeats), int(counts.loc[kept.index, "all"].sum()), int(counts.loc[kept.index, "without_twins"].sum()),
                           int(held.loc[members, "records"].sum())]
            additive_flags |= {(row.task, row.model, "naive_record_split", fold, study) for study, fold in kept.index}
            additive_flags |= {(row.task, row.model, reference, fold, study) for fold, study in zip(safe["fold"], safe["study_accession"]) if study in members}
        if len(across) >= least:
            values[8:] = [float(across["all"].mean()), float(across["without_twins"].mean()), float(across["all"].mean() - across["without_twins"].mean())]
            records[4:] = [int(sizes.loc[across.index, "all"].sum()), int(sizes.loc[across.index, "without_twins"].sum())]
            additive_flags |= {(row.task, row.model, "naive_leave_one_study_out", fold, fold) for fold in across.index}
        kept_values = [getattr(row, name) for name in M.VALUES]
        matched = (matched and close(values, kept_values) and close([kept_values[4] + kept_values[5] + kept_values[6]], [kept_values[7]])
                   and close([kept_values[7]], [kept_values[0] - kept_values[3]])
                   and [row.a_records, row.b_records, row.c_records, row.d_records, row.e_records, row.f_records] == records
                   and row.matched_studies == len(members) and row.loso_matched_studies == len(across)
                   and (row.status, row.reason) == (("estimable", "") if len(members) >= least else ("not_estimable", "fewer_than_minimum_studies"))
                   and (row.loso_status, row.loso_reason) == (("estimable", "") if len(across) >= least else ("not_estimable", "fewer_than_minimum_studies")))
    qc.check("additive_matched_set_rebuilt_from_the_study_rows", marked["in_additive_set"] == additive_flags, len(additive_flags))
    qc.check("additive_matched_parts_add_up_over_the_same_studies", matched, {f"{row.task}:{row.model}": int(row.matched_studies)
                                                                             for row in additive.itertuples(index=False)})
    receipt = C.read_json(C.receipt_path(2))
    qc.check("effect_receipt_agrees_with_the_tables", all(
        receipt["inflation_effects"][row.task][row.model][row.effect] == dict(
            status=row.status, estimate=C.rounded(row.estimate), low=C.rounded(row.low), high=C.rounded(row.high), studies=int(row.studies),
            first_records=int(row.first_records), second_records=int(row.second_records)) for row in effects.itertuples(index=False))
        and all(receipt["additive_matched"][row.task][row.model] == dict(status=row.status, loso_status=row.loso_status, matched_studies=int(row.matched_studies),
                                                                        loso_matched_studies=int(row.loso_matched_studies),
                                                                        **{name: C.rounded(getattr(row, name)) for name in M.VALUES})
                for row in additive.itertuples(index=False)))
    state["summary"].update(inflation_effects={f"{row.task}:{row.model}:{row.effect}": dict(estimate=C.rounded(row.estimate), low=C.rounded(row.low),
                                                                                            high=C.rounded(row.high), studies=int(row.studies))
                                               for row in effects.itertuples(index=False)})


def shortcut_checks(qc, state):
    populations = state["populations"]
    folds = pd.read_parquet(C.DATA / "shortcut_folds.parquet")
    predictions = pd.read_parquet(C.DATA / "shortcut_predictions.parquet")
    qc.check("shortcut_tables_have_their_columns", list(folds.columns) == S.FOLD_COLUMNS and list(predictions.columns) == S.PREDICTION_COLUMNS)
    from sklearn.metrics import roc_auc_score

    right, within, minimum = True, True, C.CONFIG["splits"]["minimum_per_class"]
    for task in [entry["task"] for entry in state["policy"]["tasks"]]:
        frame = populations[(populations["population"] == "naive") & (populations["task"] == task)].reset_index(drop=True)
        labels, studies = frame["label"].to_numpy(dtype="float64"), frame["study_accession"].to_numpy()
        plans = [("naive_record_split", C.repeat_name(repeat), random_test(frame, repeat)) for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1)]
        plans += [("naive_leave_one_study_out", study, studies == study) for study in sorted(set(studies))]
        for scheme, fold, test in plans:
            train = pd.DataFrame({"study": studies[~test], "label": labels[~test]})
            share, overall = train.groupby("study")["label"].mean().to_dict(), float(train["label"].mean())
            expected = {key: share.get(study, overall) for key, study in zip(zip(frame.loc[test, "study_accession"], frame.loc[test, "sample_accession"]), studies[test])}
            got = predictions[(predictions["task"] == task) & (predictions["scheme"] == scheme) & (predictions["fold"] == fold)]
            row = folds[(folds["task"] == task) & (folds["scheme"] == scheme) & (folds["fold"] == fold)]
            keys = list(zip(got["study_accession"], got["sample_accession"]))
            ready = int((got["label"] == 1).sum()) >= minimum and int((got["label"] == 0).sum()) >= minimum
            right = (right and len(row) == 1 and set(keys) == set(expected) and len(keys) == len(expected) and close([expected[key] for key in keys], got["score"], 1e-12)
                     and (row["status"].iloc[0] == "scored") == ready)
            if scheme == "naive_record_split":
                right = right and close([roc_auc_score(got["label"], got["score"])], row["auroc"]) and close(row["pooled_auroc"], row["auroc"])
            elif ready:
                # Every record of a study unseen in training gets one score, so its AUROC is exactly one half.
                within = within and got["score"].nunique() == 1 and float(row["auroc"].iloc[0]) == 0.5 and float(row["auroc_low"].iloc[0]) == 0.5
            else:
                right = right and row[list(S.METRICS)].isna().all(axis=None) and row["reason"].iloc[0] == "fewer_than_minimum_of_a_class"
        held = predictions[(predictions["task"] == task) & (predictions["scheme"] == "naive_leave_one_study_out")]
        rows = folds[(folds["task"] == task) & (folds["scheme"] == "naive_leave_one_study_out")]
        right = right and len(held) == len(frame) and close(rows["pooled_auroc"], np.full(len(rows), roc_auc_score(held["label"], held["score"])))
    qc.check("study_only_predictions_and_pooled_scores_follow_from_the_populations", right, len(folds))
    held = folds[(folds["scheme"] == "naive_leave_one_study_out") & (folds["status"] == "scored")]
    qc.check("every_held_out_study_scores_exactly_one_half", within and (held["auroc"] == 0.5).all() and len(held) > 0, len(held))
    identity = pd.read_parquet(C.DATA / "study_identification.parquet")
    guessed = pd.read_parquet(C.DATA / "study_identification_predictions.parquet")
    from sklearn.model_selection import StratifiedKFold

    policy = C.CONFIG["shortcuts"]["study_identification"]
    tables = [guessed[guessed["table"] == table].sort_values(["fold", "study_accession", "sample_accession"]).reset_index(drop=True) for table in S.TABLES]
    same = all(len(part) == len(tables[0]) and (part[["fold", "study_accession", "sample_accession", "unit"]].to_numpy() == tables[0][
        ["fold", "study_accession", "sample_accession", "unit"]].to_numpy()).all() for part in tables)
    samples = pd.read_parquet(C.CORRECTION / "corrected" / f"{policy['group']}_samples.parquet", columns=[*C.KEYS, "sample_id"])
    samples["unit"] = samples["study_accession"] + "__" + samples["platform"]
    samples = samples[samples.groupby("unit")["sample_id"].transform("size") >= policy["minimum_unit_samples"]].reset_index(drop=True)
    splitter = StratifiedKFold(n_splits=policy["folds"], shuffle=True, random_state=policy["seed"])
    folds_of = {}
    for number, (_, test) in enumerate(splitter.split(np.zeros(len(samples)), samples["unit"].to_numpy()), start=1):
        folds_of.update({key: f"fold_{number:02d}" for key in zip(samples.loc[test, "study_accession"], samples.loc[test, "sample_accession"])})
    first = tables[0]
    placed = len(first) == len(samples) and [folds_of.get(key) for key in zip(first["study_accession"], first["sample_accession"])] == list(first["fold"])

    def balanced(truth, guess):
        frame = pd.DataFrame({"truth": truth, "right": np.asarray(truth) == np.asarray(guess)})
        return float(frame.groupby("truth")["right"].mean().mean())

    scored = True
    for row in identity.itertuples(index=False):
        part = guessed[(guessed["table"] == row.table) & ((guessed["fold"] == row.fold) | (row.fold == "all"))]
        scored = scored and abs(balanced(part["unit"], part["predicted"]) - row.balanced_accuracy) < 1e-12 and row.test_samples == len(part)
    qc.check("study_identification_follows_from_its_folds_and_predictions", same and placed and scored and len(identity) == len(S.TABLES) * (policy["folds"] + 1),
             dict(samples=len(samples), units=int(samples["unit"].nunique())))
    overall = identity[identity["fold"] == "all"].set_index("table")["balanced_accuracy"]
    state["summary"].update(study_identification={table: round(float(value), 6) for table, value in overall.items()},
                            study_only={task: dict(random_split_auroc=C.rounded(part.loc[part["scheme"] == "naive_record_split", "auroc"].mean()),
                                                   held_out_studies_scored=int(((part["scheme"] == "naive_leave_one_study_out") & (part["status"] == "scored")).sum()))
                                        for task, part in folds.groupby("task")})


def signature_checks(qc, state):
    policy, table = C.CONFIG["signatures"], state["table"]
    found = pd.read_parquet(C.DATA / "signature_genes.parquet")
    coverage = pd.read_parquet(C.DATA / "signature_coverage.parquet")
    tests = pd.read_parquet(C.DATA / "signature_tests.parquet")
    scores = pd.read_parquet(C.DATA / "signature_scores.parquet")
    summary = pd.read_parquet(C.DATA / "signature_summary.parquet")
    qc.check("signature_tables_have_their_columns", list(found.columns) == G.GENE_COLUMNS and list(coverage.columns) == G.COVERAGE_COLUMNS
             and list(tests.columns) == G.TEST_COLUMNS and list(scores.columns) == G.SCORE_COLUMNS and list(summary.columns) == G.SUMMARY_COLUMNS)
    reference = pd.read_parquet(C.ANNOTATION / "genes.parquet", columns=["hgnc_id", "symbol"]).set_index("hgnc_id")["symbol"]
    definitions = {entry["name"]: entry for entry in policy["definitions"]}
    qc.check("every_signature_gene_is_its_hgnc_gene", len(found) == sum(len(entry["symbol_map"]) for entry in policy["definitions"])
             and (found["hgnc_id"].map(reference) == found["symbol"]).all()
             and all(set(found.loc[found["signature"] == name, "symbol"]) == set(entry["genes_up"] + entry["genes_down"]) for name, entry in definitions.items()))
    # The strata and the first reason of every test are rebuilt from the cohort table and the coverage table.
    matrices = pd.read_parquet(C.VALIDATION / "validation_matrices.parquet").set_index(["study_accession", "platform"])
    held = coverage.set_index(["signature", "study_accession", "platform"])
    expected, members, minimum = {}, {}, policy["minimum_group"]
    for comparison, rule in policy["comparisons"].items():
        value, septic = table[rule["variable"]], table["sample_class"] == "sepsis"
        allowed = ~septic | table["sepsis_phase"].isin(policy["sepsis_phases"])
        within = (table["sample_class"] == rule["within_class"]) if "within_class" in rule else True
        positive, negative = (value == rule["positive"]) & within & allowed, (value == rule["negative"]) & within & allowed
        counts = pd.DataFrame({"positives": positive, "negatives": negative}).groupby([table["study_accession"], table["platform"], table["tissue"]]).sum()
        counts = counts[(counts["positives"] > 0) & (counts["negatives"] > 0)]
        for key, part in table[positive | negative].groupby(["study_accession", "platform", "tissue"]):
            members[(comparison, *key)] = (part["sample_id"].to_numpy(), int((~part["is_primary_record"]).sum()))
        for name, entry in definitions.items():
            if comparison not in entry["intended_comparison"]:
                continue
            for (study, platform, tissue), row in counts.iterrows():
                cover = held.loc[(name, study, platform)]
                sided = all(have >= entry["minimum_genes_per_side"] for have, named in ((cover["up_present"], cover["up_published"]),
                                                                                          (cover["down_present"], cover["down_published"])) if named)
                reason = ("definition_not_verified" if not entry["verified"] else "scale_not_comparable" if not matrices.loc[(study, platform), "comparable_scale"]
                          else "fewer_than_minimum_of_a_class" if min(row["positives"], row["negatives"]) < minimum
                          else "gene_coverage_below_threshold" if cover["coverage"] < entry.get("minimum_gene_coverage", policy["minimum_gene_coverage"])
                          else "too_few_genes_on_a_side" if not sided else "")
                for variant, _ in G.variants(entry):
                    expected[(name, variant, comparison, study, platform, tissue)] = (int(row["positives"]), int(row["negatives"]), reason)
    keyed = tests.set_index(["signature", "variant", "comparison", "study_accession", "platform", "tissue"])
    follows = len(keyed) == len(expected) and set(keyed.index) == set(expected) and all(
        (keyed.at[key, "positives"], keyed.at[key, "negatives"]) == value[:2] and (keyed.at[key, "reason"] == value[2]
                                                                                    or not value[2] and keyed.at[key, "reason"] == "fewer_than_minimum_with_a_score")
        for key, value in expected.items())
    qc.check("every_stratum_is_tested_or_skipped_by_its_first_reason", follows and ((tests["status"] == "tested") == (tests["reason"] == "")).all(),
             C.counts(tests.loc[tests["reason"] != "", "reason"]))
    tested = tests[tests["status"] == "tested"]
    grouped = scores.groupby(["signature", "variant", "comparison", "study_accession", "platform", "tissue"])
    counted, worst = True, 0.0
    for row in tested.itertuples(index=False):
        part = grouped.get_group((row.signature, row.variant, row.comparison, row.study_accession, row.platform, row.tissue))
        again = pairwise(part.loc[part["label"] == 1, "score"], part.loc[part["label"] == 0, "score"])
        worst = max(worst, abs(again - row.auroc))
        patients = part.loc[part["label"] == 1, "patient"].nunique(), part.loc[part["label"] == 0, "patient"].nunique()
        counted = (counted and row.scored_positives == int((part["label"] == 1).sum()) and row.scored_negatives == int((part["label"] == 0).sum())
                   and (row.patients_positive, row.patients_negative) == patients
                   and (np.isfinite(row.auroc_low) == (min(patients) >= minimum)) and (not np.isfinite(row.auroc_low) or row.auroc_low <= row.auroc_high))
    qc.check("every_auroc_counted_again_from_the_sample_scores", counted and worst <= C.CONFIG["quality_control"]["metric_tolerance"]
             and grouped.ngroups == len(tested), worst)
    # A sample of the tested strata is scored again by gene symbol from the gene matrices, with the formulas written once more.
    order = sorted(tested.itertuples(index=False), key=lambda row: hashlib.sha256(f"{row.signature}:{row.variant}:{row.comparison}:{row.study_accession}:"
                                                                                f"{row.platform}:{row.tissue}".encode("utf-8")).hexdigest())
    largest, again_right = 0.0, True
    for row in order[:C.CONFIG["quality_control"]["signature_strata"]]:
        entry = definitions[row.signature]
        matrix = matrices.loc[(row.study_accession, row.platform)]
        frame = pd.read_parquet(C.project_path(matrix["gene_matrix"]))
        whole = frame.drop(columns="hgnc_id").set_index("symbol").astype("float64")
        whole = whole if matrix["value_scale"] in C.CONFIG["values"]["log_scales"] else np.log2(whole.clip(lower=0) + C.CONFIG["values"]["linear_offset"])
        part = grouped.get_group((row.signature, row.variant, row.comparison, row.study_accession, row.platform, row.tissue))
        up, down = [gene for gene in entry["genes_up"] if gene in whole.index], [gene for gene in entry["genes_down"] if gene in whole.index]
        chosen = whole[list(part["sample_accession"])]
        if entry["formula"] == "geometric_mean_difference":
            lowest = float(whole.min().min())
            raised = chosen + (1.0 - lowest if lowest < 1 else 0.0)
            weight = dict(G.variants(entry))[row.variant]
            weight = 1.0 if weight == "one" else len(down) / len(up) if weight == "down_over_up" else len(up) / len(down)
            score = (np.exp(np.log(raised.loc[up]).mean(axis=0)) if up else 0.0) - weight * (np.exp(np.log(raised.loc[down]).mean(axis=0)) if down else 0.0)
        else:
            score = chosen.loc[up].sum(axis=0, skipna=False) - chosen.loc[down].sum(axis=0, skipna=False)
        score = score if entry["direction"] == "higher_in_positive" else -score
        distance = np.abs(np.asarray(score, dtype="float64") - part["score"].to_numpy())
        largest = max(largest, float(distance.max()))
        again_right = again_right and abs(pairwise(part.loc[part["label"] == 1, "score"], part.loc[part["label"] == 0, "score"]) - row.auroc) < 1e-12
    qc.check("sample_of_strata_scored_again_by_gene_symbol", again_right and largest <= 1e-9, dict(strata=min(len(order), C.CONFIG["quality_control"]["signature_strata"]),
                                                                                                   largest_difference=largest))
    stage = pd.read_parquet(C.VALIDATION / "sepsis_score_by_study.parquet")
    stage = stage[stage["score_computed"]].set_index(G.STRATUM)
    mine = tested[(tested["signature"] == policy["stage_05_score"]) & (tested["variant"] == "primary") & (tested["comparison"] == "sepsis_vs_healthy_control")]
    mine = mine.set_index(G.STRATUM)
    qc.check("sepsis_metascore_equals_that_of_the_validation_stage", set(mine.index) == set(stage.index)
             and close(mine.loc[stage.index, "auroc"], stage["auroc"]) and close(mine.loc[stage.index, "auroc_low"], stage["auroc_low"])
             and close(mine.loc[stage.index, "auroc_high"], stage["auroc_high"]), len(mine))
    # Independence is decided again by the validation stage's rule: no role in the derivation, no sample shared with a study that had one, no repeat deposit.
    roles = {name: {**{study: "validation" for study in entry["derivation_studies"]["validation"]},
                    **{study: "discovery" for study in entry["derivation_studies"]["discovery"]}} for name, entry in definitions.items()}
    derived = {name: set(table.loc[table["study_accession"].isin(role_of), "sample_id"]) for name, role_of in roles.items()}
    tissues, marked = G.scope_tissues(), True
    for row in tests.itertuples(index=False):
        role = roles[row.signature].get(row.study_accession, "")
        arrays, repeats = members[(row.comparison, row.study_accession, row.platform, row.tissue)]
        shared = 0 if role else int(np.isin(arrays, list(derived[row.signature])).sum())
        marked = (marked and row.derivation_role == role and row.samples_shared_with_derivation == shared and row.repeat_deposits == repeats
                  and row.independent == (not role and shared == 0 and repeats == 0) and row.in_scope_group == (row.tissue in tissues))
    qc.check("derivation_roles_shared_samples_and_repeats_follow_the_validation_rule", marked, int(tests["independent"].sum()))
    alone = stage[(stage["score_derivation"] == "") & (stage["samples_shared_with_derivation"] == 0) & (stage["repeat_deposits"] == 0)]
    apart = mine[mine["independent"]]
    receipt = C.read_json(C.EARLIER["Validation"] / "04.json")["score"]["independent"]
    qc.check("independent_metascore_strata_are_those_of_the_validation_stage", set(apart.index) == set(alone.index) and len(apart) == receipt["strata"]
             and C.rounded(apart["auroc"].median()) == receipt["auroc_median"] and close([apart["auroc"].median()], [alone["auroc"].median()]),
             dict(strata=len(apart), auroc_median=C.rounded(apart["auroc"].median())))
    limit, rows = policy["auroc_report_limit"], []
    for key, part in tests.groupby(["signature", "variant", "comparison"]):
        for scope, inside in ((G.WHOLE_SCOPE, part), (policy["scope_group"], part[part["tissue"].isin(tissues)])):
            for kept, chosen in (("all", inside), ("independent", inside[inside["independent"]])):
                values = chosen.loc[chosen["status"] == "tested", "auroc"].to_numpy()
                rows.append((*key, scope, kept, len(chosen), len(values), np.median(values) if len(values) else np.nan,
                             np.quantile(values, 0.25) if len(values) else np.nan, np.quantile(values, 0.75) if len(values) else np.nan,
                             (values > limit).mean() if len(values) else np.nan))
    index = ["signature", "variant", "comparison", "scope", "strata_set"]
    rebuilt = pd.DataFrame(rows, columns=[*index, "strata", "tested", "auroc_median", "auroc_q1", "auroc_q3", "share_above_limit"]).set_index(index)
    stored = summary.set_index(index)
    qc.check("signature_summary_follows_from_the_tests", len(stored) == len(rebuilt) and set(stored.index) == set(rebuilt.index)
             and all(close(rebuilt.loc[stored.index, name], stored[name]) for name in ("strata", "tested", "auroc_median", "auroc_q1", "auroc_q3", "share_above_limit")))
    primary = summary[summary["primary_variant"] & (summary["strata_set"] == "all") & (summary["scope"] == G.WHOLE_SCOPE)]
    state["summary"].update(signature_median_auroc={f"{row.signature}:{row.comparison}": C.rounded(row.auroc_median) for row in primary.itertuples(index=False)})


def total_checks(qc, state):
    leakage, models, shortcuts, signatures = (C.read_json(C.receipt_path(number)) for number in (1, 2, 3, 4))
    units = pd.read_parquet(C.DATA / "model_units.parquet")
    tests = pd.read_parquet(C.DATA / "signature_tests.parquet")
    qc.check("receipt_totals_agree_with_the_tables", leakage["records"]["total"] == state["summary"]["records"]
             and leakage["patient_safe_folds"]["folds"] == state["summary"]["patient_safe_folds"]
             and models["units"]["total"] == len(units) and models["units"]["scored"] == int((units["status"] == "scored").sum())
             and models["predictions"] == len(pd.read_parquet(C.DATA / "model_predictions.parquet"))
             and shortcuts["study_identification"]["balanced_accuracy"] == {table: C.rounded(value) for table, value in state["summary"]["study_identification"].items()}
             and signatures["tests"]["total"] == len(tests) and signatures["tests"]["tested"] == int((tests["status"] == "tested").sum()))


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/7 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/7 Populations, links and split leakage", lambda: leakage_checks(qc, state))
        qc.phase("3/7 Baseline models, folds and metrics", lambda: model_checks(qc, state))
        qc.phase("4/7 Sources of the inflation", lambda: inflation_checks(qc, state))
        qc.phase("5/7 Study shortcuts", lambda: shortcut_checks(qc, state))
        qc.phase("6/7 Published signatures", lambda: signature_checks(qc, state))
        qc.phase("7/7 Stage receipts and totals", lambda: total_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, max(C.STAGES)) if C.receipt_path(number).exists()]),
                  eligibility="BENCHMARK_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[05 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    Leakage: {summary['share_exposed_random_split']:.1%} of test records exposed under a random split, "
              f"{summary['share_exposed_leave_one_study_out']:.1%} under leave-one-study-out, none in {summary['patient_safe_folds']} patient-safe folds; "
              f"study identification balanced accuracy {summary['study_identification']['aligned']:.3f} before and "
              f"{summary['study_identification']['corrected']:.3f} after the correction", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
