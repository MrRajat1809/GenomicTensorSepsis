"""Select what the batch correction works on: every sample record with its sample-type group or the reason it is left out, the study units and the genes."""
import hashlib
import importlib.metadata
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ANNOTATION = ROOT / "data/processed/annotation"
COHORT = ROOT / "data/processed/cohort"
VALIDATION = ROOT / "data/processed/validation"
EARLIER = {"Annotation": ROOT / "outputs/gene_annotation", "Cohort": ROOT / "outputs/cohort_curation", "Validation": ROOT / "outputs/technical_validation"}
DATA = ROOT / "data/processed/correction"
ALIGNED = DATA / "aligned"
CORRECTED = DATA / "corrected"
OUT = ROOT / "outputs/batch_correction"
CONFIG = json.loads((HERE / "correction_config.json").read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
STAGES = {1: "01_select_correction_inputs.py", 2: "02_align_value_distributions.py", 3: "03_select_flat_genes.py", 4: "04_correct_anchored_studies.py",
          5: "05_correct_unanchored_studies.py", 6: "06_evaluate_corrections.py", 7: "07_build_corrected_tables.py", 8: "08_correction_master_qc.py"}
KEYS = ["study_accession", "sample_accession", "platform"]
UNIT = ["group", "study_accession", "platform"]
RECORD_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "sample_id", "is_primary_record", "tissue", "sample_class", "sepsis_phase",
                  "age_group", "patient_key", "time_order", "sex", "recorded_sex_unreliable", "expression_sex", "sex_used", "female", "any_flag", "family",
                  "value_scale", "gene_matrix", "group", "unit_role", "estimation", "status", "reason"]
UNIT_COLUMNS = ["group", "source", "study_accession", "platform", "family", "value_scale", "gene_matrix", "role", "records", "primary_records", "estimation_samples",
                "healthy_anchors", "sepsis", "non_septic_control", "other_condition", "neonate_share", "neonatal", "tier", "genes", "genes_in_group_set", "status",
                "reason"]
GENE_COLUMNS = ["group", "hgnc_id", "symbol", "chromosome", "merged_units", "units_with_gene", "anchored_units_with_gene", "check_role", "in_group_set"]
REASONS = {"sample_type_not_corrected": "Its recorded sample type is in no corrected group (swab, plasma, sorted cells or not recorded).",
           "value_scale_not_comparable": "Its gene matrix is not on a log scale.",
           "study_has_too_few_samples": "Its study unit has too few estimation samples (records without a quality flag) to estimate the study's effects.",
           "newborn_study_not_corrected": "Its study unit is of newborns; no study of newborns has healthy controls against which a correction could be tested.",
           "no_class_to_line_up": "Its study unit has no healthy controls and too few samples of a class with a reference (sepsis, healthy or non-septic controls).",
           "no_method_passed_the_checks": "The correction of its study unit did not pass every check that applies to it.",
           "no_check_applicable": "Its study unit has no healthy controls, and the study-mixing check, the one check that tests such a correction, applies "
                                  "to none of its classes (too few samples); it stays in the per-study tables.",
           "repeat_deposit": "Another record of the same sample is its primary record; this record was corrected only for the comparison of deposits."}


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


def accession_order(accession):
    """Order accessions by their letters and then by their number, so that GSM99999 comes before GSM100000, as the cohort and validation stages do."""
    match = re.fullmatch(r"(.*?)([0-9]+)", str(accession))
    return (match.group(1), int(match.group(2)), str(accession)) if match else (str(accession), -1, str(accession))


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
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "correction_config.json"])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("numpy", "scipy", "pandas", "pyarrow"):
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
    require(result["code"] == code_binding(number), f"Correction code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
    verify_files(result["outputs"])
    verify_files(result["dependencies"])
    return result


def begin(number):
    validate_config()
    print(f"[{number:02d}] Verifying inputs and provenance v{VERSION}", flush=True)
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if number < max(STAGES) and (OUT / "qc.json").exists():
        write_json(OUT / "qc.json", dict(version=VERSION, status="STALE", reason=f"Stage {number:02d} is being rerun; rerun the master QC after completion"))
    for previous in range(1, number):
        verify_stage(previous)
    if number > 1:
        verify_files(read_json(receipt_path(1))["external"])
    return {"code": code_binding(number), "dependencies": bind([receipt_path(previous) for previous in range(1, number)])}


def validate_config():
    require(VERSION == "1.0.0", "Unsupported correction version")
    groups, inputs, alignment, flat, anchored, unanchored, checks, quality = (CONFIG[name] for name in (
        "groups", "inputs", "alignment", "flat_genes", "anchored", "unanchored", "checks", "quality_control"))
    tissues = [tissue for group in groups.values() for tissue in group["tissues"]]
    require(groups and len(tissues) == len(set(tissues)), "Every sample type belongs to at most one group")
    require(inputs["value_scales"] and 0 < inputs["minimum_gene_share"] <= 1 and inputs["minimum_unit_samples"] >= 3 and inputs["minimum_healthy_controls"] >= 3
            and inputs["minimum_reference_units"] >= 2 and 0 < inputs["neonatal_share"] <= 1 and "during_sepsis" not in inputs["phases_not_estimated"],
            "Invalid input policy")
    require(alignment["grid_points"] >= 100 and 0 < alignment["minimum_rank_correlation"] <= 1, "Invalid alignment policy")
    require(flat["minimum_group"] >= 3 and flat["count"] >= 50 and 0 <= flat["minimum_expression_quantile"] < 1 and flat["minimum_training_units"] >= 3
            and flat["minimum_shared_genes"] >= max(unanchored["components"]) and flat["minimum_complete_genes"] >= 1, "Invalid flat-gene policy")
    require(anchored["anchor_class"] in anchored["reference_classes"] and 0 < anchored["floor_fraction"] < 1 and anchored["minimum_prior_genes"] >= 3,
            "Invalid anchored policy")
    require(unanchored["methods"] == ["matched_groups"] and set(unanchored["matched_classes"]) <= set(anchored["reference_classes"])
            and unanchored["minimum_matched"] >= 3 and unanchored["components"] == sorted(set(unanchored["components"])) and min(unanchored["components"]) >= 1
            and unanchored["comparison"]["minimum_units"] >= 2, "Invalid unanchored policy")
    require(checks["pca_genes"] >= 100 and 2 <= checks["components"] <= checks["pca_genes"] and checks["neighbours"] >= 3 and checks["minimum_class_samples"] >= 3
            and 0 <= checks["mixing_worsening_limit"] < 1 and 0.5 < checks["sex_agreement"] <= 1 and 0.5 < checks["score_auroc"] <= 1 and checks["minimum_score_genes"] >= 1
            and checks["minimum_sex_calls"] >= 3 and set(checks["classes"]) <= set(anchored["reference_classes"]) and checks["location_report_limit"] > 0
            and checks["scale_report_limit"] > 1, "Invalid check policy")
    require(all(quality[name] > 0 for name in ("alignment_tolerance", "anchor_centre_limit", "rebuild_tolerance", "model_tolerance")), "Invalid quality-control policy")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Correction code changed during execution")
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
    """Confirm that the annotation, cohort and validation stages passed their quality control unchanged and return the artifacts that this stage reads."""
    paths = [ANNOTATION / "genes.parquet"]
    for name, folder in EARLIER.items():
        qc = read_json(folder / "qc.json")
        require(qc["status"] == "PASS" and qc["failures"] == 0, f"{name} quality control has not passed")
        verify_files(qc["code"])
        verify_files(qc["dependencies"])
        paths += [folder / "qc.json", *map(project_path, qc["dependencies"])]
    built = {**read_json(EARLIER["Cohort"] / "05.json")["outputs"], **read_json(EARLIER["Validation"] / "01.json")["outputs"],
             **read_json(EARLIER["Validation"] / "02.json")["outputs"], **read_json(EARLIER["Validation"] / "03.json")["outputs"]}
    matrices = {name: value for name, value in read_json(EARLIER["Annotation"] / "04.json")["outputs"].items()
                if "/gene_matrices/" in name or name.endswith("gene_matrix_summary.parquet")}
    verify_files(built)
    verify_files(matrices)
    return paths + [project_path(name) for name in [*built, *matrices] if name.endswith(".parquet")]


# Shared arithmetic of the stage. Values are genes in rows and samples in columns; a gene that a unit does not measure is a row of empty values.

def unit_name(study, platform):
    return f"{study}__{platform}"


def aligned_path(group, study, platform):
    return ALIGNED / group / f"{unit_name(study, platform)}.parquet"


def group_genes(genes, group):
    """The genes of a group in the order of their HGNC identifiers."""
    chosen = genes.loc[(genes["group"] == group) & genes["in_group_set"], "hgnc_id"]
    return sorted(chosen, key=lambda name: int(name.split(":")[1]))


def unit_records(records, group, study, platform):
    """The records of one study unit, in a fixed order, with a fresh index."""
    own = records[(records["group"] == group) & (records["study_accession"] == study) & (records["platform"] == platform) & (records["status"] == "in_unit")]
    return own.sort_values("sample_accession", kind="stable").reset_index(drop=True)


def read_aligned(group, study, platform, columns):
    """The aligned values of the named records of one unit: genes in rows in the order of the group, records in columns."""
    import pandas as pd

    return pd.read_parquet(aligned_path(group, study, platform), columns=list(columns)).to_numpy(dtype="float64")


def model_arrays(group, genes):
    """The sex effect, anchor mean and anchor standard deviation of every gene of a group, in the order of the group."""
    import pandas as pd

    model = pd.read_parquet(DATA / "correction_model.parquet")
    model = model[model["group"] == group].set_index("hgnc_id").reindex(genes)
    require(model.index.notna().all() and len(model) == len(genes), f"The model of {group} does not hold every gene")
    return {name: model[name].to_numpy(dtype="float64") for name in ("sex_effect", "anchor_mean", "anchor_sd")}


def class_cells(own, study):
    """The positions of the given records of one unit (with a fresh index) in every cell of sample class and age group."""
    require(own["sample_class"].notna().all() and own["age_group"].notna().all(),
            f"An estimation sample of {study} has no sample class or age group; the cohort stage gives every record both")
    return {cell: list(index) for cell, index in own.groupby(["sample_class", "age_group"]).groups.items()}


def unit_moments(group, units, records, sex_effect):
    """The moments of the sex-adjusted aligned values of the estimation samples of every unit of a group, by sample class and age group."""
    result = {}
    for unit in units[(units["group"] == group) & (units["status"] == "in_unit")].itertuples(index=False):
        own = unit_records(records, group, unit.study_accession, unit.platform)
        own = own[own["estimation"]].reset_index(drop=True)
        values = read_aligned(group, unit.study_accession, unit.platform, own["sample_accession"]) - sex_effect[:, None] * own["female"].to_numpy()[None, :]
        result[(unit.study_accession, unit.platform)] = {cell: moments(values[:, index]) for cell, index in class_cells(own, unit.study_accession).items()}
    return result


def gather(cells, classes=None, ages=None):
    """The combined moments of the cells of one unit of the given classes and age groups, or None where it has none."""
    chosen = [part for (name, age), part in cells.items() if (classes is None or name in classes) and (ages is None or age in ages)]
    return combine(chosen) if chosen else None


def class_reference(moments_of, training, mean, deviation, effects):
    """For every reference class and age group, and for all ages, the mean over training units of the corrected mean and variance within the unit."""
    import warnings

    import numpy as np

    policy = CONFIG["unanchored"]
    result = {}
    for name in CONFIG["anchored"]["reference_classes"]:
        for age in [*policy["matched_age_groups"], "all"]:
            centres, spreads = [], []
            for key in training:
                part = gather(moments_of[key], [name], None if age == "all" else [age])
                if part is None:
                    continue
                count, centre, spread = corrected_moments(part, mean, deviation, *effects[key])
                enough = count >= policy["minimum_matched"]
                centres.append(np.where(enough, centre, np.nan))
                spreads.append(np.where(enough, spread, np.nan))
            if not centres:
                continue
            centres, spreads = np.vstack(centres), np.vstack(spreads)
            units = np.isfinite(centres).sum(axis=0)
            usable = units >= CONFIG["inputs"]["minimum_reference_units"]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                result[(name, age)] = (np.where(usable, np.nanmean(centres, axis=0), np.nan), np.where(usable, np.nanmean(spreads, axis=0), np.nan), units)
    return result


def matched_samples(cells):
    """The samples that the matched-group method lines up in one unit, and their degrees of freedom within classes, at the gene that has most."""
    import numpy as np

    if not cells:
        return 0, 0
    held = [np.where(np.isfinite(centre) & np.isfinite(reference), part[0], 0.0) for part, centre, reference in cells]
    return int(np.nanmax(sum(held))), int(np.nanmax(sum(np.maximum(count - 1, 0.0) for count in held)))


def matched_cells(cells, reference, excluded=()):
    """The cells of one unit that the matched-group method lines up, each with the reference of its class and age group, or of all ages where that is missing."""
    import numpy as np

    policy = CONFIG["unanchored"]
    chosen = []
    for (name, age), part in cells.items():
        if name not in policy["matched_classes"] or name in excluded or (name, "all") not in reference:
            continue
        overall = reference[(name, "all")]
        own = reference.get((name, age)) if age in policy["matched_age_groups"] else None
        centre = overall[0] if own is None else np.where(np.isfinite(own[0]), own[0], overall[0])
        spread = overall[1] if own is None else np.where(np.isfinite(own[0]), own[1], overall[1])
        chosen.append((part, centre, spread))
    return chosen


def read_values(path, columns, genes):
    """Read the named sample columns of a matrix as genes in rows, in the order of the given HGNC identifiers, absent genes as empty rows."""
    import pandas as pd

    frame = pd.read_parquet(path, columns=["hgnc_id", *columns]).set_index("hgnc_id")
    require(frame.index.is_unique, f"Repeated gene rows in {path}")
    return frame.reindex(list(genes))[list(columns)].to_numpy(dtype="float64")


def moments(values):
    """Count, mean and sum of squared deviations of every gene over the columns that hold a value."""
    import numpy as np

    present = ~np.isnan(values)
    count = present.sum(axis=1).astype("float64")
    total = np.where(present, values, 0.0).sum(axis=1)
    mean = np.divide(total, count, out=np.full(len(count), np.nan), where=count > 0)
    squares = np.where(present, (values - mean[:, None]) ** 2, 0.0).sum(axis=1)
    return count, mean, np.where(count > 0, squares, np.nan)


def combine(parts):
    """Combine the moments of several groups of samples of the same genes into the moments of all of them."""
    import numpy as np

    parts = list(parts)
    count = sum(part[0] for part in parts)
    total = sum(np.where(part[0] > 0, part[0] * part[1], 0.0) for part in parts)
    mean = np.divide(total, count, out=np.full(len(count), np.nan), where=count > 0)
    squares = sum(np.where(part[0] > 0, part[2] + part[0] * (part[1] - mean) ** 2, 0.0) for part in parts)
    return count, mean, np.where(count > 0, squares, np.nan)


def variance(part):
    import numpy as np

    count, _, squares = part
    return np.divide(squares, count - 1, out=np.full(len(count), np.nan), where=count > 1)


def shrink_variance(estimate, freedom):
    """Draw the variance of every gene toward an inverse-gamma prior fitted across the genes by moments, as ComBat does."""
    import numpy as np

    usable = np.isfinite(estimate) & (freedom >= 1)
    require(usable.sum() >= CONFIG["anchored"]["minimum_prior_genes"], "Too few genes to fit the prior of the scale")
    mean, spread = float(estimate[usable].mean()), float(estimate[usable].var(ddof=1))
    shape, rate = (2 * spread + mean ** 2) / spread, (mean * spread + mean ** 3) / spread
    return np.where(usable, (0.5 * freedom * estimate + rate) / (0.5 * freedom + shape - 1), np.nan)


def shrink_location(estimate, noise):
    """Draw the location of every gene toward a normal prior fitted across the genes, weighing its own estimate by its sampling variance."""
    import numpy as np

    usable = np.isfinite(estimate) & np.isfinite(noise)
    require(usable.sum() >= CONFIG["anchored"]["minimum_prior_genes"], "Too few genes to fit the prior of the location")
    mean, spread = float(estimate[usable].mean()), float(estimate[usable].var(ddof=1))
    return np.where(usable, (spread * estimate + noise * mean) / (spread + noise), np.nan)


def relative_floor(spread, power=1):
    """The floor of a spread over the genes: anchored.floor_fraction of its median over the genes that hold it, for a standard deviation (power 1), or the
    square of that fraction times the median, for a variance (power 2)."""
    import numpy as np

    held = spread[np.isfinite(spread)]
    require(len(held) > 0 and np.median(held) > 0, "No spread over the genes to set a floor from")
    return CONFIG["anchored"]["floor_fraction"] ** power * float(np.median(held))


def anchor_model(parts):
    """The anchor mean and standard deviation of every gene from the moments of the anchor samples of several units, and the number of units with the gene.
    The standard deviation is pooled within units and raised to the relative floor, so that a gene almost constant in the anchor samples (one not expressed,
    in RNA sequencing) does not inflate the standardised values of the units that measure it with noise."""
    import numpy as np

    count = sum(part[0] for part in parts)
    units = sum((part[0] >= 2).astype(int) for part in parts)
    held = sum((part[0] > 0).astype(int) for part in parts)
    total = sum(np.where(part[0] > 0, part[0] * part[1], 0.0) for part in parts)
    squares = sum(np.where(part[0] > 0, part[2], 0.0) for part in parts)
    usable = (units >= CONFIG["inputs"]["minimum_reference_units"]) & (count > 0)
    mean = np.divide(total, count, out=np.full(len(count), np.nan), where=usable)
    # Every unit that holds a gene spends one degree of freedom on its own mean.
    deviation = np.sqrt(np.divide(squares, count - held, out=np.full(len(count), np.nan), where=usable))
    deviation = np.where(np.isfinite(deviation), np.maximum(deviation, relative_floor(deviation)), np.nan)
    return mean, deviation, units


def standardised(part, mean, deviation):
    """The moments of standardised values from the moments of values: the count, the mean and the variance."""
    count, centre, _ = part
    return count, (centre - mean) / deviation, variance(part) / deviation ** 2


def anchored_effects(part, mean, deviation):
    """The location and scale of one unit from the moments of its healthy controls, shrunk by empirical Bayes."""
    import numpy as np

    count, location, spread = standardised(part, mean, deviation)
    scale2 = shrink_variance(spread, count - 1)
    location = shrink_location(np.where(count >= 2, location, np.nan), scale2 / np.maximum(count, 1))
    return location, np.sqrt(scale2)


def corrected_moments(part, mean, deviation, location, scale):
    """The mean and variance, in standardised units, of a group of samples after its unit's correction."""
    count, centre, spread = standardised(part, mean, deviation)
    return count, (centre - location) / scale, spread / scale ** 2


def matched_effects(cells, mean, deviation):
    """The location and scale of one unit that line its class groups up with their references: cells are (moments, reference mean, reference variance).
    The scale makes the unit's spread within its classes that of the references; with it, the location lines the class means up with the references.
    A cell counts for a gene only where its reference holds that gene; a reference variance is raised to the relative floor of that reference."""
    import numpy as np

    valid = [np.isfinite(centre) & np.isfinite(reference) & (part[0] > 0) for part, centre, reference in cells]
    floored = [np.where(ok, np.maximum(reference, relative_floor(reference, 2)), 1.0) if ok.any() else np.ones(len(ok)) for ok, (_, _, reference) in zip(valid, cells)]
    freedom = sum(np.where(ok & (part[0] > 1), part[0] - 1, 0.0) for ok, (part, _, _) in zip(valid, cells))
    pooled = sum(np.where(ok & (part[0] > 1), part[2] / deviation ** 2 / reference, 0.0) for ok, (part, _, _), reference in zip(valid, cells, floored))
    scale2 = shrink_variance(np.divide(pooled, freedom, out=np.full(len(freedom), np.nan), where=freedom >= 1), freedom)
    scale = np.sqrt(scale2)
    count = sum(np.where(ok, part[0], 0.0) for ok, (part, _, _) in zip(valid, cells))
    shift = sum(np.where(ok, part[0] * ((part[1] - mean) / deviation - scale * np.where(ok, centre, 0.0)), 0.0) for ok, (part, centre, _) in zip(valid, cells))
    expected = sum(np.where(ok, part[0] * reference, 0.0) for ok, (part, _, _), reference in zip(valid, cells, floored))
    usable = (count >= 2) & np.isfinite(scale)
    location = np.divide(shift, count, out=np.full(len(count), np.nan), where=usable)
    noise = scale2 * np.divide(expected, count, out=np.full(len(count), np.nan), where=usable) / np.maximum(count, 1)
    return shrink_location(location, noise), np.where(usable, scale, np.nan)


def plain_effects(part, mean, deviation):
    """The location and scale of one unit over all of its samples, shrunk by empirical Bayes: plain ComBat with the study as batch."""
    return anchored_effects(part, mean, deviation)


def correct(values, female, model, location, scale):
    """Corrected values of one unit: the sex effect removed, standardised, the unit's effects removed, returned to the anchor scale, the sex effect restored."""
    sex = model["sex_effect"][:, None] * female[None, :]
    standard = (values - sex - model["anchor_mean"][:, None]) / model["anchor_sd"][:, None]
    return (standard - location[:, None]) / scale[:, None] * model["anchor_sd"][:, None] + model["anchor_mean"][:, None] + sex


def fit_components(features, targets, components):
    """Regress every gene's effect over the training units on the leading principal components of their flat-gene effects; return a predictor."""
    import numpy as np

    centre = features.mean(axis=0)
    _, _, axes = np.linalg.svd(features - centre, full_matrices=False)
    basis = axes[:components].T
    design = np.column_stack([np.ones(len(features)), (features - centre) @ basis])
    coefficients = np.full((design.shape[1], targets.shape[1]), np.nan)
    missing = np.isnan(targets)
    patterns = {}
    for gene, pattern in enumerate(map(bytes, np.packbits(missing, axis=0).T)):
        patterns.setdefault(pattern, []).append(gene)
    for genes in patterns.values():
        rows = ~missing[:, genes[0]]
        if rows.sum() >= design.shape[1] + 1:
            coefficients[:, genes] = np.linalg.lstsq(design[rows], targets[np.ix_(rows, genes)], rcond=None)[0]
    return lambda new: np.concatenate([[1.0], (new - centre) @ basis]) @ coefficients


def choose_components(features, targets):
    """The number of components with the lowest error when every training unit is predicted from the others, over the genes all of them hold."""
    import numpy as np

    policy = CONFIG["flat_genes"]
    complete = ~np.isnan(targets).any(axis=0)
    require(complete.sum() >= policy["minimum_complete_genes"], "Too few genes held by every training unit")
    require(len(features) >= policy["minimum_training_units"], "Too few training units for the flat-gene method")
    # Leaving one unit out must leave more units than parameters (intercept and components) to fit, so at most len(features) - 3 components are tried.
    options = [value for value in CONFIG["unanchored"]["components"] if value <= len(features) - 3]
    if not options:
        # Three training units leave nothing to leave out: the smallest number of components is taken.
        return min(CONFIG["unanchored"]["components"]), {}
    squared = {components: [] for components in options}
    for left in range(len(features)):
        keep = np.arange(len(features)) != left
        train, target = features[keep], targets[np.ix_(keep, complete)]
        centre = train.mean(axis=0)
        _, _, axes = np.linalg.svd(train - centre, full_matrices=False)
        for components in options:
            basis = axes[:components].T
            design = np.column_stack([np.ones(len(train)), (train - centre) @ basis])
            coefficients = np.linalg.lstsq(design, target, rcond=None)[0]
            predicted = np.concatenate([[1.0], (features[left] - centre) @ basis]) @ coefficients
            squared[components].append(float(np.mean((predicted - targets[left, complete]) ** 2)))
    errors = {components: float(np.mean(values)) for components, values in squared.items()}
    best = min(errors, key=lambda value: (errors[value], value))
    return best, errors


def quantile_function(values, grid):
    """The sorted values of one sample read at the shared quantile positions."""
    import numpy as np

    values = np.sort(values[~np.isnan(values)])
    return np.interp(grid, (np.arange(len(values)) + 0.5) / len(values), values)


def align_sample(values, grid, reference):
    """Replace the values of one sample by the reference value at the quantile position of their rank; ties take their mean rank."""
    import numpy as np
    from scipy.stats import rankdata

    aligned = np.full(len(values), np.nan)
    present = ~np.isnan(values)
    aligned[present] = np.interp((rankdata(values[present]) - 0.5) / present.sum(), grid, reference)
    return aligned


def auroc(cases, others):
    """Area under the ROC curve of a score: the share of case and other pairs in which the case has the higher score, ties counted as half."""
    import numpy as np
    import pandas as pd

    cases, others = np.asarray(cases, dtype=float), np.asarray(others, dtype=float)
    ranks = pd.Series(np.concatenate([cases, others])).rank().to_numpy()
    return float((ranks[:len(cases)].sum() - len(cases) * (len(cases) + 1) / 2) / (len(cases) * len(others)))


def two_means(values):
    """Split values into two groups by two-means in one dimension, started at the lowest and the highest value; return the midpoint."""
    import numpy as np

    low, high = float(np.min(values)), float(np.max(values))
    if low == high:
        return low
    for _ in range(100):
        midpoint = (low + high) / 2
        new = (float(values[values <= midpoint].mean()), float(values[values > midpoint].mean()))
        if new == (low, high):
            break
        low, high = new
    return (low + high) / 2


def flat_selection(effects, genes, excluded=()):
    """Select the flat genes of one group from the differences between sepsis and healthy controls of its units, without the excluded units."""
    import numpy as np
    import pandas as pd

    policy = CONFIG["flat_genes"]
    keyed = effects[~(effects["study_accession"] + "|" + effects["platform"]).isin({f"{study}|{platform}" for study, platform in excluded})]
    # Every unit's difference weighs by its precision, n1 n2 / (n1 + n2), so that small units do not decide which genes look flat.
    keyed = keyed.assign(weight=keyed["sepsis"] * keyed["healthy_control"] / (keyed["sepsis"] + keyed["healthy_control"]))
    wide = keyed.pivot(index="hgnc_id", columns=["study_accession", "platform"], values="effect")
    require(wide.shape[1] >= 2, "Too few units to select flat genes")
    weights = keyed.pivot(index="hgnc_id", columns=["study_accession", "platform"], values="weight").reindex(index=wide.index, columns=wide.columns)
    allowed = genes[genes["in_group_set"] & (genes["check_role"] == "") & ~genes["chromosome"].isin(policy["excluded_chromosomes"])]["hgnc_id"]
    keep = wide.index.isin(set(allowed)) & wide.notna().all(axis=1)
    effect, weight = wide[keep].to_numpy(), weights[keep].to_numpy()
    table = pd.DataFrame({"units": wide.shape[1], "flatness": np.sqrt((weight * effect ** 2).sum(axis=1) / weight.sum(axis=1)),
                          "pooled_effect": (weight * effect).sum(axis=1) / weight.sum(axis=1),
                          "healthy_mean": keyed.groupby("hgnc_id")["healthy_mean"].mean().reindex(wide.index[keep]).to_numpy()}, index=wide.index[keep])
    table = table[table["healthy_mean"] >= table["healthy_mean"].quantile(policy["minimum_expression_quantile"])]
    chosen = table.sort_values(["flatness", "healthy_mean"], ascending=[True, False], kind="stable").head(policy["count"]).index
    return table.assign(selected=table.index.isin(chosen)).rename_axis("hgnc_id").reset_index()


def main():
    import numpy as np
    import pandas as pd
    import pyarrow.parquet as pq

    print(f"[01] Selecting correction inputs v{VERSION}", flush=True)
    start = begin(1)
    external = earlier_binding()
    inputs = CONFIG["inputs"]
    cohort = pd.read_parquet(COHORT / "cohort_samples.parquet")
    summary = pd.read_parquet(ANNOTATION / "gene_matrix_summary.parquet", columns=["accession", "platform", "family", "value_scale", "gene_matrix"])
    summary = summary.rename(columns={"accession": "study_accession"})
    sex = pd.read_parquet(VALIDATION / "sex_check_samples.parquet", columns=[*KEYS, "expression_sex", "recorded_sex_unreliable"])
    records = cohort.merge(summary, on=["study_accession", "platform"], how="left", validate="many_to_one").merge(sex, on=KEYS, how="left", validate="one_to_one")
    require(len(records) == len(cohort) and records["gene_matrix"].notna().all(), "Every record of the cohort table needs its gene matrix")
    for column in ("family", "value_scale", "gene_matrix", "expression_sex"):
        records[column] = records[column].fillna("")
    records["recorded_sex_unreliable"] = records["recorded_sex_unreliable"].eq(True)

    groups = {tissue: name for name, group in CONFIG["groups"].items() for tissue in group["tissues"]}
    records["group"] = records["tissue"].map(groups).fillna("")
    records["reason"] = np.select([records["group"] == "", ~records["value_scale"].isin(inputs["value_scales"])],
                                  ["sample_type_not_corrected", "value_scale_not_comparable"], "")
    records.loc[records["reason"] != "", "group"] = ""
    taken = records["reason"] == ""
    has_primary = records.groupby(UNIT)["is_primary_record"].transform("any") & taken
    records["unit_role"] = np.where(taken, np.where(has_primary, "merged", "check_only"), "")
    eligible = taken & ~records["any_flag"] & (records["is_primary_record"] | ~has_primary) & ~records["sepsis_phase"].isin(inputs["phases_not_estimated"])
    # One sample of every patient within a unit estimates its effects: the earliest by its time order, then the lowest accession in natural order; a sample
    # without an identified patient is its own patient.
    patient = records["patient_key"].where(records["patient_key"] != "", records["sample_id"])
    ordered = records[eligible].assign(patient=patient[eligible], no_time=records.loc[eligible, "time_order"].isna()).sort_values(
        [*UNIT, "patient", "no_time", "time_order", "sample_accession"], kind="stable",
        key=lambda column: column.map(accession_order) if column.name == "sample_accession" else column)
    records["estimation"] = records.index.isin(ordered.drop_duplicates([*UNIT, "patient"]).index)

    # The sex of a record: from its own expression, else as recorded (unless the recorded sex of its matrix is unreliable), else the share of females among
    # the estimation samples of its unit.
    recorded = records["sex"].isin(["female", "male"]) & ~records["recorded_sex_unreliable"]
    known = np.select([records["expression_sex"].isin(["female", "male"]), recorded], ["expression", "recorded"], "")
    records["female"] = np.select([known == "expression", known == "recorded"], [records["expression_sex"] == "female", records["sex"] == "female"], np.nan)
    records["sex_used"] = known
    share = records[records["estimation"] & (known != "")].groupby(UNIT)["female"].mean()
    for key, part in records[taken & (known == "")].groupby(UNIT):
        records.loc[part.index, "female"] = share.get(key, 0.5)
        records.loc[part.index, "sex_used"] = "share_of_unit" if key in share.index else "one_half"

    rows, held = [], {}
    for key, part in records[taken].groupby(UNIT, sort=True):
        group, study, platform = key
        first = part.iloc[0]
        estimation = part[part["estimation"]]
        require(set(pq.read_schema(project_path(first["gene_matrix"])).names) >= set(part["sample_accession"]), f"Samples missing from the matrix of {study}")
        frame = pd.read_parquet(project_path(first["gene_matrix"]), columns=["hgnc_id", *part["sample_accession"]])
        values = frame.drop(columns="hgnc_id").to_numpy(dtype="float64")
        held[key] = set(frame.loc[~np.isnan(values).all(axis=1), "hgnc_id"])
        count = lambda name: int((estimation["sample_class"] == name).sum())
        neonates = float((estimation["age_group"] == "neonate").mean()) if len(estimation) else 0.0
        rows.append(dict(group=group, source=first["source"], study_accession=study, platform=platform, family=first["family"], value_scale=first["value_scale"],
                         gene_matrix=first["gene_matrix"], role=first["unit_role"], records=len(part), primary_records=int(part["is_primary_record"].sum()),
                         estimation_samples=len(estimation), healthy_anchors=count(CONFIG["anchored"]["anchor_class"]), sepsis=count("sepsis"),
                         non_septic_control=count("non_septic_control"), other_condition=count("other_condition"), neonate_share=round(neonates, 6),
                         neonatal=neonates >= inputs["neonatal_share"], tier="", genes=len(held[key]), genes_in_group_set=0, status="", reason=""))
    units = pd.DataFrame(rows, columns=UNIT_COLUMNS)
    small = units["estimation_samples"] < inputs["minimum_unit_samples"]
    units["status"] = np.where(small, "left_out", "in_unit")
    units["reason"] = np.where(small, "study_has_too_few_samples", "")
    units["tier"] = np.where(units["healthy_anchors"] >= inputs["minimum_healthy_controls"], "anchored", "unanchored")
    require(not (~small & (units["tier"] == "anchored") & units["neonatal"]).any(),
            "An anchored unit is of newborns; unanchored.basis and the reason newborn_study_not_corrected rest on none being so: revise them first")
    left = set(map(tuple, units.loc[small, UNIT].to_numpy()))
    out = [tuple(key) in left for key in records[UNIT].to_numpy()]
    records.loc[out, "reason"] = "study_has_too_few_samples"
    records.loc[out, ["unit_role", "estimation"]] = ["", False]
    records["status"] = np.where(records["reason"] == "", "in_unit", "left_out")

    # The genes of a group: measured by enough of its units that hold primary records, and by enough anchored units.
    reference = pd.read_parquet(ANNOTATION / "genes.parquet", columns=["hgnc_id", "symbol", "location"]).set_index("hgnc_id")
    chromosome = reference["location"].fillna("").str.extract(r"^(mitochondria|X|Y|\d+)")[0].fillna("")
    checks = pd.read_parquet(VALIDATION / "validation_genes.parquet")
    roles = (checks["test"] + ":" + checks["role"]).groupby(checks["hgnc_id"]).agg(lambda part: "; ".join(sorted(set(part)))).to_dict()
    gene_rows = []
    for group in CONFIG["groups"]:
        active = units[(units["group"] == group) & (units["status"] == "in_unit")]
        merged = active[active["role"] == "merged"]
        require(len(merged) >= inputs["minimum_reference_units"], f"Too few study units in group {group}")
        counts, anchored = {}, {}
        for unit in merged.itertuples(index=False):
            for gene in held[(group, unit.study_accession, unit.platform)]:
                counts[gene] = counts.get(gene, 0) + 1
                anchored[gene] = anchored.get(gene, 0) + (unit.tier == "anchored")
        needed = math.ceil(inputs["minimum_gene_share"] * len(merged) - 1e-9)
        for gene in sorted(counts):
            gene_rows.append(dict(group=group, hgnc_id=gene, symbol=reference.at[gene, "symbol"], chromosome=chromosome.get(gene, ""), merged_units=len(merged),
                                  units_with_gene=counts[gene], anchored_units_with_gene=anchored[gene], check_role=roles.get(gene, ""),
                                  in_group_set=counts[gene] >= needed and anchored[gene] >= inputs["minimum_reference_units"]))
    genes = pd.DataFrame(gene_rows, columns=GENE_COLUMNS)
    for position, unit in units.iterrows():
        chosen = set(genes.loc[(genes["group"] == unit["group"]) & genes["in_group_set"], "hgnc_id"])
        units.at[position, "genes_in_group_set"] = len(held[(unit["group"], unit["study_accession"], unit["platform"])] & chosen)

    records = records[RECORD_COLUMNS].sort_values(["source", "study_accession", "platform", "sample_accession"], kind="stable").reset_index(drop=True)
    units = units.sort_values(UNIT, kind="stable").reset_index(drop=True)
    require(set(records["reason"]) <= set(REASONS) | {""}, "Unknown reason")
    save_table(records, DATA / "correction_records.parquet")
    save_table(units, DATA / "correction_units.parquet")
    save_table(genes, DATA / "correction_genes.parquet")
    save_review(units, OUT / "correction_units_review.csv")
    by_group = {}
    for group in CONFIG["groups"]:
        part = units[(units["group"] == group) & (units["status"] == "in_unit")]
        by_group[group] = dict(units=len(part), merged_units=int((part["role"] == "merged").sum()), check_only_units=int((part["role"] == "check_only").sum()),
                               anchored_units=int((part["tier"] == "anchored").sum()), newborn_units=int(part["neonatal"].sum()),
                               records=int(part["records"].sum()), primary_records=int(part["primary_records"].sum()),
                               estimation_samples=int(part["estimation_samples"].sum()),
                               genes=int(genes.loc[genes["group"] == group, "in_group_set"].sum()),
                               units_left_out=int(((units["group"] == group) & (units["status"] == "left_out")).sum()))
    finish(1, start, [DATA / "correction_records.parquet", DATA / "correction_units.parquet", DATA / "correction_genes.parquet",
                      OUT / "correction_units_review.csv"], external=bind(external), method=inputs["method"], groups=by_group,
           estimation=dict(records=int(records["estimation"].sum()), left_out_by_phase=int((taken & records["sepsis_phase"].isin(inputs["phases_not_estimated"])).sum()),
                           later_samples_of_a_patient_in_the_unit=int((eligible & ~records["estimation"] & (records["status"] == "in_unit")).sum())),
           records=dict(total=len(records), in_a_unit=int((records["status"] == "in_unit").sum()),
                        by_reason={str(name): int(count) for name, count in records.loc[records["reason"] != "", "reason"].value_counts().sort_index().items()},
                        sex_used={str(name or "not_applicable"): int(count) for name, count in records["sex_used"].value_counts().sort_index().items()}),
           reasons=REASONS, protocol=CONFIG["protocol"])
    print("[01] Complete: " + "; ".join(f"{group}: {value['units']} study units ({value['anchored_units']} anchored, {value['check_only_units']} for checks only), "
                                        f"{value['primary_records']:,} primary records, {value['genes']:,} genes" for group, value in by_group.items())
          + f"; review sheet: {(OUT / 'correction_units_review.csv').relative_to(ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
