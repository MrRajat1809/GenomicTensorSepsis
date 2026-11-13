"""Count, for every benchmark task and split scheme, the test records that have a deposit of the same array or a sample of the same patient in training."""
import hashlib
import importlib.metadata
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ANNOTATION = ROOT / "data/processed/annotation"
COHORT = ROOT / "data/processed/cohort"
VALIDATION = ROOT / "data/processed/validation"
CORRECTION = ROOT / "data/processed/correction"
EARLIER = {"Annotation": ROOT / "outputs/gene_annotation", "Cohort": ROOT / "outputs/cohort_curation", "Validation": ROOT / "outputs/technical_validation",
           "Correction": ROOT / "outputs/batch_correction"}
DATA = ROOT / "data/processed/benchmark"
OUT = ROOT / "outputs/benchmark"
CONFIG = json.loads((HERE / "benchmark_config.json").read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
STAGES = {1: "01_audit_split_leakage.py", 2: "02_train_baseline_models.py", 3: "03_test_study_shortcuts.py", 4: "04_retest_published_signatures.py",
          5: "05_benchmark_master_qc.py"}
KEYS = ["source", "study_accession", "sample_accession", "platform"]
NAIVE_SCHEMES = ("naive_record_split", "naive_leave_one_study_out")
GTS_SCHEMES = ("leave_one_group_out", "technology_transfer", "temporal")
POPULATIONS = ("all_records", "naive", "gts")
MODELS = ("logistic_regression", "gradient_boosting")
FORMULAS = ("geometric_mean_difference", "sum_difference")
WEIGHTS = ("one", "down_over_up", "up_over_down")
DIRECTIONS = ("higher_in_positive", "lower_in_positive")
TYPES = ("signature", "expression_marker")
POPULATION_COLUMNS = ["population", "task", *KEYS, "sample_id", "patient_key", "patient_identified", "patient_safe_group", "label", "linked_records",
                      "same_array_links", "same_patient_links", "linked_in_other_study"]
LEAKAGE_COLUMNS = ["population", "task", "scheme", "fold", "records", "test_records", "train_records", "exposed", "exposed_same_array", "exposed_same_patient",
                   "share_exposed", "expected_share", "linked_cohort_records"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_path(recorded):
    path = (ROOT / recorded).resolve()
    require(path.is_relative_to(ROOT), f"Path outside project: {recorded}")
    return path


def record(path):
    path = Path(path).resolve()
    return {"sha256": sha(path), "size_bytes": path.stat().st_size}


def bind(paths):
    paths = sorted(set(Path(path).resolve() for path in paths), key=lambda path: path.as_posix())
    return {path.relative_to(ROOT).as_posix(): record(path) for path in paths}


def verify_files(records):
    for name, expected in records.items():
        path = project_path(name)
        require(path.is_file() and record(path) == expected, f"Missing or changed artifact: {name}")


def code_binding(number=max(STAGES)):
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "benchmark_config.json"])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("numpy", "scipy", "pandas", "pyarrow", "scikit-learn"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def receipt_path(number):
    return OUT / f"{number:02d}.json"


def verify_stage(number):
    result = read_json(receipt_path(number))
    require(result["version"] == VERSION and result["status"] == "COMPLETE", f"Stage {number} incomplete")
    require(result["code"] == code_binding(number), f"Benchmark code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
    verify_files(result["outputs"])
    verify_files(result["dependencies"])
    return result


def earlier_stages(number):
    # Every analysis reads the populations of stage 01 only; the analyses do not build on one another.
    return [1] if 1 < number < max(STAGES) else list(range(1, number))


def begin(number):
    validate_config()
    print(f"[{number:02d}] Verifying inputs and provenance v{VERSION}", flush=True)
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if number < max(STAGES) and (OUT / "qc.json").exists():
        write_json(OUT / "qc.json", dict(version=VERSION, status="STALE", reason=f"Stage {number:02d} is being rerun; rerun the master QC after completion"))
    for previous in earlier_stages(number):
        verify_stage(previous)
    if number > 1:
        verify_files(read_json(receipt_path(1))["external"])
    return {"code": code_binding(number), "dependencies": bind([receipt_path(previous) for previous in earlier_stages(number)])}


def validate_config():
    require(VERSION == "1.0.0", "Unsupported benchmark version")
    splits, features, models, metrics = (CONFIG[name] for name in ("splits", "features", "models", "metrics"))
    shortcuts, signatures, values, quality = (CONFIG[name] for name in ("shortcuts", "signatures", "values", "quality_control"))
    require(isinstance(splits["seed"], int) and splits["repeats"] >= 1 and 0 < splits["test_fraction"] < 1 and splits["minimum_per_class"] >= 2,
            "Invalid split policy")
    require(features["group"] and 0 < features["minimum_present_share"] <= 1 and features["genes"] >= 10 and features["chunk_genes"] >= 100,
            "Invalid feature policy")
    linear, boosting = models["logistic_regression"], models["gradient_boosting"]
    require(set(models) == {"method", *MODELS} and linear["C"] > 0 and linear["max_iter"] >= 100 and linear["tol"] > 0 and boosting["learning_rate"] > 0
            and boosting["max_iter"] >= 1 and boosting["max_leaf_nodes"] >= 2 and boosting["min_samples_leaf"] >= 1 and boosting["l2_regularization"] >= 0
            and isinstance(boosting["random_state"], int), "Invalid model policy")
    require(metrics["bootstrap"]["replicates"] >= 200 and 0 < metrics["bootstrap"]["level"] < 1 and isinstance(metrics["bootstrap"]["seed"], int)
            and metrics["inflation_reference"] in GTS_SCHEMES, "Invalid metric policy")
    effects = metrics["inflation_effects"]
    require(isinstance(effects["minimum_studies"], int) and effects["minimum_studies"] >= 2 and effects["bootstrap"]["replicates"] >= 200
            and 0 < effects["bootstrap"]["level"] < 1 and isinstance(effects["bootstrap"]["seed"], int), "Invalid inflation-effect policy")
    identify = shortcuts["study_identification"]
    require(identify["group"] and identify["folds"] >= 2 and identify["genes"] >= 10 and identify["minimum_unit_samples"] >= identify["folds"]
            and identify["C"] > 0 and identify["max_iter"] >= 100 and identify["tol"] > 0, "Invalid study-identification policy")
    comparisons, names = signatures["comparisons"], [entry["name"] for entry in signatures["definitions"]]
    require(signatures["minimum_group"] >= 3 and signatures["sepsis_phases"] and signatures["bootstrap"]["replicates"] >= 200
            and 0 < signatures["bootstrap"]["level"] < 1 and 0 < signatures["minimum_gene_coverage"] <= 1 and 0.5 < signatures["auroc_report_limit"] < 1
            and signatures["scope_group"]
            and all({"variable", "positive", "negative"} <= set(rule) for rule in comparisons.values()), "Invalid signature policy")
    require(names and len(names) == len(set(names)) and signatures["stage_05_score"] in names, "Every signature needs a unique name, and the score of stage 05 one")
    for entry in signatures["definitions"]:
        name, intended, up, down = entry["name"], entry["intended_comparison"], entry["genes_up"], entry["genes_down"]
        require(entry["type"] in TYPES and entry["citation"] and entry["doi"] and entry["pmid"] and entry["source"],
                f"Signature {name} needs a type of {TYPES}, a citation, its DOI, PMID and source")
        require(intended and set(intended) <= set(comparisons) and len(intended) == len(set(intended)),
                f"Signature {name} needs its intended comparisons, each named in signatures.comparisons")
        require(entry["direction"] in DIRECTIONS and entry["formula"] in FORMULAS and isinstance(entry["verified"], bool),
                f"Signature {name} needs a direction of {DIRECTIONS}, a formula of {FORMULAS} and verified true or false; add a new formula to the code first")
        variants = [entry["down_weight"], *(variant["down_weight"] for variant in entry["sensitivity"])]
        labels = ["primary", *(variant["variant"] for variant in entry["sensitivity"])]
        require(set(variants) <= set(WEIGHTS) and len(labels) == len(set(labels)) and (entry["formula"] != "sum_difference" or set(variants) == {"one"}),
                f"Signature {name} needs down weights of {WEIGHTS}, a sum only the weight one, and distinct variant names")
        require(up or down, f"Signature {name} names no genes")
        require(len(up) == len(set(up)) and len(down) == len(set(down)) and not set(up) & set(down), f"Signature {name} names a gene twice")
        require(sorted(row["hgnc_symbol"] for row in entry["symbol_map"]) == sorted(up + down) and all(row["published"] and row["hgnc_id"] for row in entry["symbol_map"]),
                f"The symbol map of signature {name} must give the published symbol and HGNC identifier of every one of its genes once")
        require(entry["minimum_genes_per_side"] >= 1 and 0 < entry.get("minimum_gene_coverage", signatures["minimum_gene_coverage"]) <= 1
                and set(entry["derivation_studies"]) >= {"discovery", "validation"}, f"Invalid gene limits or derivation studies of signature {name}")
    require(values["linear_offset"] > 0 and values["log_scales"] and not set(values["log_scales"]) & set(values["linear_scales"]), "Invalid value policy")
    require(quality["metric_tolerance"] > 0 and quality["refit_tolerance"] > 0 and quality["signature_strata"] >= 1 and quality["refit"]["scheme"] in GTS_SCHEMES,
            "Invalid quality-control policy")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Benchmark code changed during execution")
    verify_files(start["dependencies"])
    result = dict(version=VERSION, status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=runtime(), outputs=bind(outputs), **start, **detail)
    write_json(receipt_path(number), result)
    return result


def save_table(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(path)


def save_review(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(path)


def earlier_binding():
    """Confirm that the annotation, cohort, validation and correction stages passed their quality control unchanged and return the artifacts that this
    stage reads."""
    paths = [ANNOTATION / "genes.parquet"]
    for name, folder in EARLIER.items():
        qc = read_json(folder / "qc.json")
        require(qc["status"] == "PASS" and qc["failures"] == 0, f"{name} quality control has not passed")
        verify_files(qc["code"])
        verify_files(qc["dependencies"])
        paths += [folder / "qc.json", *map(project_path, qc["dependencies"])]
    built = {}
    for name, numbers in (("Cohort", (5, 6)), ("Validation", (1, 4)), ("Correction", (1, 2, 7))):
        for number in numbers:
            built.update(read_json(EARLIER[name] / f"{number:02d}.json")["outputs"])
    matrices = {name: value for name, value in read_json(EARLIER["Annotation"] / "04.json")["outputs"].items()
                if "/gene_matrices/" in name or name.endswith("gene_matrix_summary.parquet")}
    verify_files(built)
    verify_files(matrices)
    return paths + [project_path(name) for name in [*built, *matrices] if name.endswith(".parquet")]


def cohort_policy():
    """The benchmark policy of the cohort stage, as its receipt records it: the tasks, their tissues and the minimum of each class in a scored group."""
    policy = read_json(EARLIER["Cohort"] / "06.json")["policy"]
    require(policy["minimum_per_class"] == CONFIG["splits"]["minimum_per_class"],
            "splits.minimum_per_class must equal the minimum of each class of the cohort stage's benchmark")
    return policy


# Shared arithmetic of the stage. A population is a frame of cohort records with a fresh index; a split is a pair of boolean arrays over it.

def naive_members(table, task, tissues):
    """The records that a naive pooling takes for a task: every record of its tissues and its two classes, without removing repeats, flags or phases."""
    value = table[task["variable"]]
    chosen = table["tissue"].isin(tissues) & value.isin([task["positive"], task["negative"]])
    if "within_class" in task:
        chosen &= table["sample_class"] == task["within_class"]
    return chosen


def task_frames(table, tasks, policy):
    """The naive and the patient-safe population of every task as (population, task, frame); a frame holds the records of the cohort table with a label."""
    frames = []
    for task in policy["tasks"]:
        name = task["task"]
        naive = table[naive_members(table, task, policy["tissues"])]
        naive = naive.assign(label=(naive[task["variable"]] == task["positive"]).astype("int8")).reset_index(drop=True)
        chosen = tasks[(tasks["task"] == name) & tasks["eligible"]][[*KEYS, "label"]]
        safe = table.merge(chosen, on=KEYS, how="inner", validate="one_to_one").assign(label=lambda frame: frame["label"].astype("int8"))
        require(len(safe) == len(chosen), f"Every sample of task {name} must be a record of the cohort table")
        frames += [("naive", name, naive), ("gts", name, safe.reset_index(drop=True))]
    return frames


def twin_pairs(frame):
    """Every ordered pair of linked records of a frame, by position: two deposits of one array (kind 0), or samples of two arrays of one identified
    patient (kind 1). Unidentified samples share an empty patient key, which is never grouped on."""
    import numpy as np
    import pandas as pd

    arrays, patients = frame["sample_id"].to_numpy(), frame["patient_key"].to_numpy()
    require((arrays != "").all() and (frame["patient_identified"].to_numpy() == (patients != "")).all(),
            "Every record needs its array, and exactly the identified patients a patient key")
    first, second, kind = [np.zeros(0, dtype="int64")], [np.zeros(0, dtype="int64")], [np.zeros(0, dtype="int8")]
    for code, keys in ((0, arrays), (1, patients)):
        held = np.flatnonzero(keys != "")
        groups = pd.Series(held).groupby(keys[held]).indices
        for name in sorted(groups):
            positions = held[groups[name]]
            if len(positions) < 2:
                continue
            one, other = np.repeat(positions, len(positions)), np.tile(positions, len(positions))
            keep = (one != other) & ((arrays[one] != arrays[other]) if code == 1 else True)
            first.append(one[keep])
            second.append(other[keep])
            kind.append(np.full(int(keep.sum()), code, dtype="int8"))
    return np.concatenate(first), np.concatenate(second), np.concatenate(kind)


def exposed(pairs, test, train):
    """Mark the test records that have a linked record in training: by any link, by the same array and by the same patient."""
    import numpy as np

    first, second, kind = pairs
    hit = test[first] & train[second]
    marks = {}
    for name, chosen in (("any", hit), ("same_array", hit & (kind == 0)), ("same_patient", hit & (kind == 1))):
        mark = np.zeros(len(test), dtype=bool)
        mark[first[chosen]] = True
        marks[name] = mark
    return marks


def record_split(frame, repeat):
    """The test records of one repeat of the naive random split: the first test_fraction of the records, ordered by SHA-256 of seed, repeat and record."""
    import numpy as np

    policy = CONFIG["splits"]
    digests = np.array([hashlib.sha256(f"{policy['seed']}:{repeat}:{study}:{sample}".encode("utf-8")).hexdigest()
                        for study, sample in zip(frame["study_accession"], frame["sample_accession"])], dtype=object)
    test = np.zeros(len(frame), dtype=bool)
    test[np.argsort(digests, kind="stable")[:math.floor(policy["test_fraction"] * len(frame) + 0.5)]] = True
    return test


def repeat_name(repeat):
    return f"repeat_{repeat:02d}"


def gts_folds(splits, task):
    """The folds of the patient-safe schemes of one task: scheme, fold, training groups, scored test groups and held-out groups that are not scored."""
    part = splits[splits["task"] == task]
    for (scheme, fold), rows in part.groupby(["scheme", "fold"], sort=False):
        groups = {role: sorted(rows.loc[rows["role"] == role, "patient_safe_group"]) for role in ("train", "test", "held_out_not_scored")}
        require(groups["train"] and groups["test"] and not set(groups["train"]) & set(groups["test"] + groups["held_out_not_scored"]),
                f"Fold {fold} of {scheme} for {task} must hold separate training and test groups")
        yield scheme, fold, groups["train"], groups["test"], groups["held_out_not_scored"]


def unit_seed(*parts):
    """A seed for one unit of work from the configured seed and the names of the unit, so that every unit draws the same numbers in every run."""
    text = ":".join(str(part) for part in (CONFIG["metrics"]["bootstrap"]["seed"], *parts))
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:16], 16)


def auroc_rows(scores, cases):
    """AUROC of every row of a score array whose first `cases` columns are positives: the share of positive and negative pairs ordered right, ties as
    half."""
    from scipy.stats import rankdata

    others = scores.shape[1] - cases
    ranks = rankdata(scores, axis=1)
    return (ranks[:, :cases].sum(axis=1) - cases * (cases + 1) / 2) / (cases * others)


def average_precision_rows(scores, truth):
    """Average precision of every row, as scikit-learn computes it: the precision at every distinct score, weighed by the gain in recall there."""
    import numpy as np

    order = np.argsort(-scores, axis=1, kind="stable")
    ordered, labels = np.take_along_axis(scores, order, axis=1), np.take_along_axis(truth, order, axis=1)
    hits = np.cumsum(labels, axis=1)
    last = np.ones(ordered.shape, dtype=bool)
    last[:, :-1] = ordered[:, 1:] != ordered[:, :-1]
    held = np.where(last, hits, 0)
    before = np.zeros(hits.shape, dtype=hits.dtype)
    before[:, 1:] = np.maximum.accumulate(held, axis=1)[:, :-1]
    gain = np.where(last, hits - before, 0)
    return (gain * hits / np.arange(1, ordered.shape[1] + 1)).sum(axis=1) / hits[:, -1]


def score_unit(labels, scores, seed):
    """AUROC and AUPRC of one scored unit by scikit-learn, with percentile intervals from a bootstrap that resamples the records of each class."""
    import numpy as np
    from sklearn.metrics import average_precision_score, roc_auc_score

    policy = CONFIG["metrics"]["bootstrap"]
    labels, scores = np.asarray(labels, dtype="int64"), np.asarray(scores, dtype="float64")
    cases, others = np.flatnonzero(labels == 1), np.flatnonzero(labels == 0)
    auroc, auprc = float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))
    ordered = np.concatenate([scores[cases], scores[others]])[None, :]
    truth = np.concatenate([np.ones(len(cases)), np.zeros(len(others))])[None, :]
    require(abs(auroc_rows(ordered, len(cases))[0] - auroc) < 1e-9 and abs(average_precision_rows(ordered, truth)[0] - auprc) < 1e-9,
            "The bootstrap and scikit-learn disagree on a metric")
    generator = np.random.default_rng(seed)
    drawn = np.concatenate([cases[generator.integers(0, len(cases), (policy["replicates"], len(cases)))],
                            others[generator.integers(0, len(others), (policy["replicates"], len(others)))]], axis=1)
    truth = np.concatenate([np.ones(len(cases)), np.zeros(len(others))])[None, :].repeat(policy["replicates"], axis=0)
    tail = (1 - policy["level"]) / 2
    found = {}
    for name, values in (("auroc", auroc_rows(scores[drawn], len(cases))), ("auprc", average_precision_rows(scores[drawn], truth))):
        found[f"{name}_low"], found[f"{name}_high"] = float(np.quantile(values, tail)), float(np.quantile(values, 1 - tail))
    return dict(auroc=auroc, auprc=auprc, prevalence=float(len(cases) / len(labels)), **found)


def rounded(value, digits=6):
    """A number for a receipt: rounded, or null where it is missing or not finite."""
    return None if value is None or not math.isfinite(float(value)) else round(float(value), digits)


def counts(series):
    return {str(name): int(count) for name, count in series.value_counts().sort_index().items()}


def scorable(labels):
    """Whether a unit holds the minimum of each class to be scored."""
    import numpy as np

    labels = np.asarray(labels)
    minimum = CONFIG["splits"]["minimum_per_class"]
    return int((labels == 1).sum()) >= minimum and int((labels == 0).sum()) >= minimum


def leakage_row(population, task, scheme, fold, pairs, test, train, expected=None, cohort=None):
    marks = exposed(pairs, test, train)
    count = int(test.sum())
    return dict(population=population, task=task, scheme=scheme, fold=fold, records=len(test), test_records=count, train_records=int(train.sum()),
                exposed=int(marks["any"].sum()), exposed_same_array=int(marks["same_array"].sum()), exposed_same_patient=int(marks["same_patient"].sum()),
                share_exposed=float(marks["any"].sum() / count) if count else float("nan"), expected_share=float("nan") if expected is None else expected,
                linked_cohort_records=cohort)


def main():
    import numpy as np
    import pandas as pd

    print(f"[01] Auditing split leakage v{VERSION}", flush=True)
    start = begin(1)
    external = earlier_binding()
    policy = cohort_policy()
    table = pd.read_parquet(COHORT / "cohort_samples.parquet")
    tasks = pd.read_parquet(COHORT / "benchmark_tasks.parquet")
    splits = pd.read_parquet(COHORT / "benchmark_splits.parquet")
    require(not table.duplicated(["study_accession", "sample_accession"]).any(), "Every record of the cohort table must be keyed once")
    whole = table.reset_index(drop=True)
    whole_pairs = twin_pairs(whole)
    where = pd.MultiIndex.from_frame(whole[KEYS])
    share = CONFIG["splits"]["test_fraction"]

    frames = [("all_records", "", whole.assign(label=pd.NA))] + task_frames(whole, tasks, policy)
    linked_pairs = {(population, task): whole_pairs if population == "all_records" else twin_pairs(frame) for population, task, frame in frames}
    rows, population_parts = [], []
    for population, task, frame in frames:
        pairs = linked_pairs[(population, task)]
        first, second, kind = pairs
        studies = frame["study_accession"].to_numpy()
        linked = np.bincount(first, minlength=len(frame))
        population_parts.append(frame[KEYS + ["sample_id", "patient_key", "patient_identified", "patient_safe_group"]].assign(
            population=population, task=task, label=pd.array(frame["label"], dtype="Int8"), linked_records=linked,
            same_array_links=np.bincount(first[kind == 0], minlength=len(frame)), same_patient_links=np.bincount(first[kind == 1], minlength=len(frame)),
            linked_in_other_study=np.bincount(first[studies[first] != studies[second]], minlength=len(frame)) > 0)[POPULATION_COLUMNS])
        if population == "gts":
            continue
        expected = float(np.mean(1 - share ** linked))
        for repeat in range(1, CONFIG["splits"]["repeats"] + 1):
            test = record_split(frame, repeat)
            rows.append(leakage_row(population, task, "naive_record_split", repeat_name(repeat), pairs, test, ~test, expected=expected))
        for study in sorted(set(studies)):
            test = studies == study
            rows.append(leakage_row(population, task, "naive_leave_one_study_out", study, pairs, test, ~test))

    # The patient-safe folds are audited twice: against the samples of the task in the training groups, and against every record of the cohort in them.
    for population, task, frame in frames:
        if population != "gts":
            continue
        pairs, groups = linked_pairs[(population, task)], frame["patient_safe_group"].to_numpy()
        places = where.get_indexer(pd.MultiIndex.from_frame(frame[KEYS]))
        require((places >= 0).all(), f"Every sample of task {task} must be a record of the cohort table")
        for scheme, fold, train_groups, test_groups, unscored in gts_folds(splits, task):
            test, train = np.isin(groups, test_groups + unscored), np.isin(groups, train_groups)
            wide_test = np.zeros(len(whole), dtype=bool)
            wide_test[places[test]] = True
            wide_train = whole["patient_safe_group"].isin(train_groups).to_numpy()
            cohort = int(exposed(whole_pairs, wide_test, wide_train)["any"].sum())
            rows.append(leakage_row(population, task, scheme, fold, pairs, test, train, cohort=cohort))
    leakage = pd.DataFrame(rows, columns=LEAKAGE_COLUMNS)
    leakage["linked_cohort_records"] = leakage["linked_cohort_records"].astype("Int64")
    safe = leakage[leakage["population"] == "gts"]
    require(len(safe) == splits.groupby(["task", "scheme", "fold"]).ngroups, "Every fold of the cohort stage's split table must be audited")
    require((safe["exposed"] == 0).all() and (safe["linked_cohort_records"] == 0).all(),
            "A patient-safe fold of the cohort stage holds a linked record on both sides; rerun the cohort stage and its master QC")
    populations = pd.concat(population_parts, ignore_index=True)
    populations["label"] = populations["label"].astype("Int8")

    summary, review = {}, []
    for (population, task, scheme), part in leakage.groupby(["population", "task", "scheme"], sort=False):
        tested = int(part["test_records"].sum())
        entry = dict(population=population, task=task, scheme=scheme, folds=len(part), test_records=tested, exposed=int(part["exposed"].sum()),
                     exposed_same_array=int(part["exposed_same_array"].sum()), exposed_same_patient=int(part["exposed_same_patient"].sum()),
                     share_exposed=round(float(part["exposed"].sum() / tested), 6) if tested else None,
                     expected_share=round(float(part["expected_share"].iloc[0]), 6) if scheme == "naive_record_split" else None)
        review.append(entry)
        summary.setdefault(task or population, {})[f"{population}:{scheme}" if task else scheme] = {key: value for key, value in entry.items()
                                                                                                   if key not in ("population", "task", "scheme")}
    review = pd.DataFrame(review)
    every = populations[populations["population"] == "all_records"]
    spread = every.groupby("sample_id")["study_accession"].nunique()
    leakage = leakage.sort_values(["population", "task", "scheme", "fold"], kind="stable").reset_index(drop=True)
    populations = populations.sort_values(["population", "task", "study_accession", "sample_accession"], kind="stable").reset_index(drop=True)
    save_table(populations, DATA / "split_populations.parquet")
    save_table(leakage, DATA / "split_leakage.parquet")
    save_review(review, OUT / "split_leakage_review.csv")
    finish(1, start, [DATA / "split_populations.parquet", DATA / "split_leakage.parquet", OUT / "split_leakage_review.csv"], external=bind(external),
           method=CONFIG["splits"]["method"], records=dict(total=len(every), with_a_linked_record=int((every["linked_records"] > 0).sum()),
                                                           share_with_a_linked_record=round(float((every["linked_records"] > 0).mean()), 6),
                                                           arrays_in_more_than_one_study=int((spread > 1).sum()),
                                                           studies_sharing_an_array=int(every.loc[every["sample_id"].isin(spread[spread > 1].index),
                                                                                                  "study_accession"].nunique())),
           populations={f"{population}:{task}" if task else population: int(count) for (population, task), count in
                        populations.groupby(["population", "task"], sort=True).size().items()},
           leakage=summary, patient_safe_folds=dict(folds=len(safe), exposed=0, linked_cohort_records=0), split_policy=CONFIG["splits"],
           protocol=CONFIG["protocol"])
    headline = summary["all_records"]
    print(f"[01] Complete: {headline['naive_record_split']['share_exposed']:.1%} of test records exposed under a naive random split (expected "
          f"{headline['naive_record_split']['expected_share']:.1%}), {headline['naive_leave_one_study_out']['share_exposed']:.1%} under leave-one-study-out; "
          f"{len(safe)} patient-safe folds with no linked record across the split; review sheet: {(OUT / 'split_leakage_review.csv').relative_to(ROOT).as_posix()}",
          flush=True)


if __name__ == "__main__":
    main()
