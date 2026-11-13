"""Verify the figures of the resource: their provenance, the exact set of files, every image, and the plotted numbers against the stage tables."""
import argparse
import hashlib
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

import visualization_common as V

MEANING = ("PASS means that every figure script completed with the current configuration and code, and so passed the layout check of every file it "
           "saved; that the stages whose tables the figures read are unchanged since their quality control; that the figure directory holds exactly "
           "the recorded files, each unchanged, each PDF a PDF and each PNG at 300 dpi with the recorded size, not blank and with nothing drawn on its "
           "outermost pixels; that figure_data.csv documents every figure; that the plotted numbers of the flow, the samples by year, the composition "
           "and its flows, the recorded variables, the genes per table, the sex check of every table and unit, the repeated deposits, the sepsis score, "
           "the cell markers, the quality flags and the batch correction and its ledger equal the stage tables row by row; that every drawn ROC curve "
           "encloses the area the validation stage reports; that the principal components of the stacked view and of the corrected tables agree "
           "with their stated shares of variance and with the shares explained by technology, study and class; and that the benchmark figures count "
           "every test record of every fold with the kind of twin that the folds rebuilt again give, in shares that add up to the leakage table and "
           "its receipt (no twin in a patient-safe fold), draw every scored unit of the baseline models with medians and middle halves that the units "
           "and the summary give, the two AUROCs of every inflation effect with its stored difference and the band of its interval, ROC curves of "
           "the study-only score that the stored predictions give again and that enclose the stored AUROCs (one half within every held-out study), "
           "shares of the own study platform predicted that the stored predictions give again and whose means are the stored balanced accuracies, "
           "exactly the tested strata of every signature with the medians of the summary, every stratum of the signature tests, and ROC curves of "
           "the tested strata that the stored scores give again and that enclose the AUROC of their test. It does not judge whether a figure is well "
           "designed.")
# Two pixels at 300 dpi: the anti-aliasing of a line cut by the edge reaches the outermost two pixels, while the layout keeps a margin of a dozen.
EDGE = 2
# The figure data hold the component scores rounded to four decimals on scales of tens to hundreds, which moves the centring, the correlation and the
# ratio of variances by about a millionth; a thousandth leaves room for that and still fails scores of another component or another table.
COMPONENTS = 1e-3


def explained_by(scores, labels):
    """The R squared of a one-way layout: the share of the variance of the scores that the means of the labels explain."""
    table = pd.DataFrame({"score": np.asarray(scores, dtype=float), "label": np.asarray(labels)})
    total = ((table["score"] - table["score"].mean()) ** 2).sum()
    return float(((table.groupby("label")["score"].transform("mean") - table["score"].mean()) ** 2).sum() / total)


def components_agree(first, second, explained_1, explained_2):
    """Two principal-component scores are centred, uncorrelated and in the ratio of their stated shares of variance."""
    first, second = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    scale = np.sqrt(first.var())
    return bool(np.isfinite(first).all() and np.isfinite(second).all() and abs(first.mean()) <= COMPONENTS * scale and abs(second.mean()) <= COMPONENTS * scale
                and abs(np.corrcoef(first, second)[0, 1]) <= COMPONENTS and abs(first.var() / second.var() / (explained_1 / explained_2) - 1) <= COMPONENTS)


def same_rows(plotted, table, keys, columns, tolerance=1e-4):
    """Whether every plotted row equals the row of the stage table with the same keys, and both hold the same rows; figure_data.csv writes a missing
    value as an empty field."""
    left = plotted.set_index(keys)[columns].replace("", "nan").astype(float).sort_index()
    right = table.set_index(keys)[columns].astype(float).sort_index()
    return len(left) == len(right) and left.index.equals(right.index) and bool(((left - right).abs().fillna(0) <= tolerance).all().all()
                                                                            and (left.isna() == right.isna()).all().all())


def same_labels(plotted, table, keys, columns):
    """Whether the plotted rows carry the labels and flags of the stage table, row for row; figure_data.csv writes a flag as True or False."""
    left = plotted.set_index(keys)[columns].sort_index()
    right = table.set_index(keys)[columns].astype(str).sort_index()
    return len(left) == len(right) and left.index.equals(right.index) and bool((left.to_numpy() == right.to_numpy()).all())


def twin_kinds(frame, test, train):
    """The kind of twin of every test record among the training records, found by joining the test records with the training records on the array
    and on the identified patient."""
    tested = frame[test].reset_index()
    trained = frame[train]
    by_array = set(tested.merge(trained[["sample_id"]].drop_duplicates(), on="sample_id")["index"])
    known = trained.loc[trained["patient_key"] != "", ["patient_key", "sample_id"]].drop_duplicates()
    joined = tested[tested["patient_key"] != ""].merge(known, on="patient_key", suffixes=("", "_other"))
    by_patient = set(joined.loc[joined["sample_id"] != joined["sample_id_other"], "index"])
    array, patient = tested["index"].isin(by_array).to_numpy(), tested["index"].isin(by_patient).to_numpy()
    return tested[["study_accession", "sample_accession"]].assign(twin=np.select([array & patient, array, patient], ["both", "array_only", "patient_only"],
                                                                                 "none"))


def rebuilt_records(tasks):
    """Every test record of every repeat of the naive random split, of the naive leave-one-study-out folds and of the patient-safe folds, with the kind
    of its twin in training, rebuilt from the populations of the benchmark stage and the split table of the cohort stage."""
    populations = pd.read_parquet(V.ROOT / "data/processed/benchmark/split_populations.parquet")
    splits = pd.read_parquet(V.ROOT / "data/processed/cohort/benchmark_splits.parquet")
    settings = V.read_json(V.ROOT / "src/07_benchmark/benchmark_config.json")["splits"]
    rows = []
    for task in tasks:
        naive = populations[(populations["population"] == "naive") & (populations["task"] == task)].reset_index(drop=True)
        folds = []
        for repeat in range(1, settings["repeats"] + 1):
            digests = [hashlib.sha256(f"{settings['seed']}:{repeat}:{study}:{sample}".encode("utf-8")).hexdigest()
                       for study, sample in zip(naive["study_accession"], naive["sample_accession"])]
            chosen = sorted(range(len(naive)), key=digests.__getitem__)[:math.floor(settings["test_fraction"] * len(naive) + 0.5)]
            test = naive.index.isin(chosen)
            folds.append(("naive_record_split", "naive_record_split", f"repeat_{repeat:02d}", naive, test, ~test))
        folds += [("naive_leave_one_study_out", "naive_leave_one_study_out", study, naive, (naive["study_accession"] == study).to_numpy(),
                   (naive["study_accession"] != study).to_numpy()) for study in sorted(set(naive["study_accession"]))]
        safe = populations[(populations["population"] == "gts") & (populations["task"] == task)].reset_index(drop=True)
        for (scheme, fold), part in splits[splits["task"] == task].groupby(["scheme", "fold"]):
            folds.append(("patient_safe", scheme, fold, safe, safe["patient_safe_group"].isin(part.loc[part["role"] != "train", "patient_safe_group"]).to_numpy(),
                          safe["patient_safe_group"].isin(part.loc[part["role"] == "train", "patient_safe_group"]).to_numpy()))
        for split, scheme, fold, frame, test, train in folds:
            rows.append(twin_kinds(frame, test, train).assign(task=task, split=split, scheme=scheme, fold=fold))
    return pd.concat(rows, ignore_index=True)


def roc_points(labels, scores):
    """The ROC curve of a score that is higher in the positive class, from the highest score down, one point for every distinct score."""
    table = pd.DataFrame({"score": np.asarray(scores, dtype=float), "positive": np.asarray(labels) == 1, "negative": np.asarray(labels) == 0})
    counts = table.groupby("score", sort=True)[["positive", "negative"]].sum().iloc[::-1].cumsum()
    return np.r_[0.0, counts["negative"] / counts["negative"].iloc[-1]], np.r_[0.0, counts["positive"] / counts["positive"].iloc[-1]]


def roc_area(fpr, tpr):
    fpr, tpr = np.asarray(fpr, dtype=float), np.asarray(tpr, dtype=float)
    return float(((fpr[1:] - fpr[:-1]) * (tpr[1:] + tpr[:-1]) / 2).sum())


def curves_agree(drawn, scores, keys, stored, tolerance=1e-6):
    """Whether every drawn ROC curve, read in the order of its points, is the curve that the stored scores give again, and encloses the stored AUROC of
    its stratum or fold; and whether the drawn curves are exactly those of the scores."""
    drawn = drawn.assign(point=drawn["point"].astype(float), fpr=drawn["fpr"].astype(float), tpr=drawn["tpr"].astype(float))
    found = {key: part.sort_values("point") for key, part in drawn.groupby(keys, sort=True)}
    expected = {key: roc_points(part["label"], part["score"]) for key, part in scores.groupby(keys, sort=True)}
    if set(found) != set(expected) or set(found) != set(stored):
        return False
    for key, part in found.items():
        fpr, tpr = expected[key]
        if len(part) != len(fpr) or not (np.allclose(part["fpr"], fpr, rtol=0, atol=tolerance) and np.allclose(part["tpr"], tpr, rtol=0, atol=tolerance)
                                         and abs(roc_area(part["fpr"], part["tpr"]) - stored[key]) <= 2 * tolerance
                                         and abs(roc_area(fpr, tpr) - stored[key]) <= 1e-9):
            return False
    return True


def benchmark(qc, data, config):
    """The plotted numbers of Fig7 and SuppFig7 against the tables and receipts of the benchmark stage."""
    policy = config["benchmark"]
    folder = V.ROOT / "data/processed/benchmark"
    models, tasks = policy["model_order"], policy["task_order"]

    def panel(figure, name):
        return data[(data["figure"] == figure) & (data["panel"] == name)]

    # A: every test record of every fold, counted in the rings of its task by the kind of its twin in training.
    records = panel("Fig7", "twins")
    leakage = pd.read_parquet(folder / "split_leakage.parquet")
    stored = leakage[leakage["population"].isin(["naive", "gts"])].set_index(["task", "scheme", "fold"]).sort_index()
    counted = records.assign(exposed=records["twin"] != "none", exposed_same_array=records["twin"].isin(["array_only", "both"]),
                             exposed_same_patient=records["twin"].isin(["patient_only", "both"])).groupby(["task", "scheme", "fold"]).agg(
        test_records=("twin", "size"), exposed=("exposed", "sum"), exposed_same_array=("exposed_same_array", "sum"),
        exposed_same_patient=("exposed_same_patient", "sum")).sort_index()
    columns = ["test_records", "exposed", "exposed_same_array", "exposed_same_patient"]
    qc.check("donut_records_equal_the_leakage_table", counted.index.equals(stored.index) and (counted[columns].to_numpy() == stored[columns].to_numpy()).all()
             and set(records["twin"]) <= set(policy["twin_order"]))
    keys = ["task", "split", "scheme", "fold", "study_accession", "sample_accession"]
    rebuilt = rebuilt_records(tasks).set_index(keys)["twin"].sort_index()
    plotted = records.set_index(keys)["twin"].sort_index()
    qc.check("donut_records_equal_the_fold_rebuild", plotted.index.is_unique and plotted.index.equals(rebuilt.index) and (plotted == rebuilt).all())
    shares = panel("Fig7", "twin_shares")
    recount = records.groupby(["task", "split", "twin"]).size().unstack(fill_value=0).reindex(columns=policy["twin_order"], fill_value=0).stack()
    recount = recount.rename("records").reset_index()
    recount["total"] = recount.groupby(["task", "split"])["records"].transform("sum")
    recount["share"] = recount["records"] / recount["total"]
    receipt = V.read_json(V.ROOT / "outputs/benchmark/01.json")
    exposed = recount[recount["twin"] != "none"].groupby(["task", "split"])["records"].sum()
    totals = recount.groupby(["task", "split"])["total"].first()
    agree = [exposed[(task, "patient_safe")] == 0 and totals[(task, "patient_safe")]
             == sum(entry["test_records"] for key, entry in receipt["leakage"][task].items() if key.startswith("gts:"))
             and all(totals[(task, split)] == receipt["leakage"][task][f"naive:{split}"]["test_records"]
                     and exposed[(task, split)] == receipt["leakage"][task][f"naive:{split}"]["exposed"]
                     and abs(exposed[(task, split)] / totals[(task, split)] - receipt["leakage"][task][f"naive:{split}"]["share_exposed"]) <= 1e-6
                     for split in ("naive_record_split", "naive_leave_one_study_out")) for task in tasks]
    safe = records[records["split"] == "patient_safe"]
    qc.check("donut_shares_equal_the_records_and_the_leakage_receipt", same_rows(shares, recount, ["task", "split", "twin"], ["records", "total", "share"],
                                                                                 tolerance=1e-6)
             and len(agree) == len(tasks) and all(agree) and set(shares["split"]) == set(policy["split_order"])
             and safe[["task", "scheme", "fold"]].drop_duplicates().shape[0] == receipt["patient_safe_folds"]["folds"]
             and receipt["patient_safe_folds"]["exposed"] == 0)

    # B: every scored unit of the two models, and the median and middle half of each cloud.
    units = pd.read_parquet(folder / "model_units.parquet")
    units = units[(units["status"] == "scored") & units["model"].isin(models)]
    plotted = panel("Fig7", "baselines")
    qc.check("raincloud_units_equal_the_scored_units", same_rows(plotted, units, ["task", "model", "scheme", "fold", "unit"], ["records", "auroc"],
                                                                 tolerance=1e-6))
    summary = pd.read_parquet(folder / "model_summary.parquet")
    summary = summary[summary["model"].isin(models)]
    quantiles = units.groupby(["task", "model", "scheme"])["auroc"].quantile([0.25, 0.5, 0.75]).unstack()
    quantiles.columns = ["auroc_q1", "auroc_median", "auroc_q3"]
    quantiles = quantiles.join(summary.set_index(["task", "model", "scheme"])[["units_scored", "auroc_median"]].rename(columns={"auroc_median": "stated"}))
    qc.check("raincloud_medians_and_middle_halves_equal_the_units_and_the_summary",
             same_rows(panel("Fig7", "baseline_quantiles"), quantiles.reset_index(), ["task", "model", "scheme"],
                       ["auroc_q1", "auroc_median", "auroc_q3", "units_scored"], tolerance=1e-6)
             and np.allclose(quantiles["auroc_median"], quantiles["stated"], rtol=0, atol=1e-6) and len(quantiles) == len(summary))

    # C: the two ends of every arrow, the difference beside it and the band of its interval around the second AUROC.
    effects = pd.read_parquet(folder / "model_inflation_effects.parquet")
    effects = effects[effects["model"].isin(models)]
    plotted = panel("Fig7", "inflation")
    qc.check("arrow_ends_equal_the_first_and_second_auroc", same_rows(plotted, effects, ["task", "model", "effect"], ["first", "second", "studies"],
                                                                      tolerance=1e-6) and sorted(set(plotted["effect"])) == sorted(policy["effect_order"]))
    totals = plotted[plotted["effect"] == "total_inflation"]
    banded = plotted[plotted["effect"] != "total_inflation"].astype({name: float for name in ("first", "low", "high", "band_low", "band_high")})
    reference = summary[summary["scheme"] == "naive_record_split"].set_index(["task", "model"])["inflation_auroc"]
    qc.check("arrow_labels_and_bands_equal_the_estimates_and_intervals", same_rows(plotted, effects, ["task", "model", "effect"], ["estimate", "low", "high"],
                                                                                   tolerance=1e-6)
             and np.allclose(plotted["first"].astype(float) - plotted["second"].astype(float), plotted["estimate"].astype(float), rtol=0, atol=2e-6)
             and np.allclose(banded["band_low"], banded["first"] - banded["high"], rtol=0, atol=2e-6)
             and np.allclose(banded["band_high"], banded["first"] - banded["low"], rtol=0, atol=2e-6)
             and (totals[["low", "high", "band_low", "band_high"]] == "").all().all()
             and all(abs(float(row.estimate) - reference[(row.task, row.model)]) <= 1e-6 for row in totals.itertuples(index=False)))

    # D: the ROC curve of the study-only score of every repeat of the naive random split, and one half within every held-out study.
    folds = pd.read_parquet(folder / "shortcut_folds.parquet")
    pooled = folds[(folds["scheme"] == "naive_record_split") & (folds["scoring"] == "pooled")]
    within = folds[(folds["scheme"] == "naive_leave_one_study_out") & (folds["scoring"] == "within_study") & (folds["status"] == "scored")]
    predictions = pd.read_parquet(folder / "shortcut_predictions.parquet")
    stated = V.read_json(V.ROOT / "outputs/benchmark/03.json")["study_only"]
    curves = panel("Fig7", "shortcut_curves")
    qc.check("shortcut_curves_recomputed_from_the_predictions_enclose_the_stored_aurocs",
             curves_agree(curves, predictions[predictions["scheme"] == "naive_record_split"], ["task", "fold"],
                          {(row.task, row.fold): float(row.auroc) for row in pooled.itertuples(index=False)})
             and all(abs(pooled.loc[pooled["task"] == task, "auroc"].mean() - stated[task]["random_split"]["pooled_auroc_mean"]) <= 1e-6
                     and (pooled["task"] == task).sum() == stated[task]["random_split"]["repeats"] for task in tasks))
    plotted = panel("Fig7", "shortcut_scores")
    held = plotted[plotted["scoring"] == "within_study"]
    qc.check("study_only_scores_equal_the_benchmark_stage", same_rows(plotted, pd.concat([pooled, within]), ["task", "scheme", "scoring", "fold", "unit"],
                                                                      ["records", "auroc"], tolerance=1e-6)
             and len(held) == len(within) > 0 and (held["auroc"].astype(float) == 0.5).all()
             and all((within["task"] == task).sum() == stated[task]["held_out_study"]["studies_scored"] for task in tasks))

    # E: the share of every study platform's samples predicted as its own, from the aligned and from the corrected values.
    predictions = pd.read_parquet(folder / "study_identification_predictions.parquet")
    hits = predictions.assign(hit=predictions["unit"] == predictions["predicted"]).groupby(["unit", "table"])["hit"].agg(["mean", "size"]).unstack("table")
    recount = pd.DataFrame({"aligned": hits[("mean", "aligned")], "corrected": hits[("mean", "corrected")], "samples": hits[("size", "aligned")],
                            "samples_corrected": hits[("size", "corrected")]}).reset_index()
    plotted = panel("Fig7", "identification")
    qc.check("identification_shares_recomputed_from_the_predictions", same_rows(plotted, recount, ["unit"], ["aligned", "corrected", "samples"], tolerance=1e-6)
             and (recount["samples"] == recount["samples_corrected"]).all())
    identification = pd.read_parquet(folder / "study_identification.parquet")
    stored = identification[identification["fold"] == "all"].set_index("table")
    receipt = V.read_json(V.ROOT / "outputs/benchmark/03.json")["study_identification"]
    # The bold line of the slopes is the mean of the drawn shares, which must be the stored balanced accuracy; its chance line one in the platforms.
    means = panel("Fig7", "identification_means").set_index("table")
    qc.check("identification_mean_line_is_the_stored_balanced_accuracy", sorted(stored.index) == ["aligned", "corrected"] == sorted(means.index) and all(
        abs(plotted[table].astype(float).mean() - float(stored.loc[table, "balanced_accuracy"])) <= 1e-6
        and abs(float(means.loc[table, "mean"]) - float(stored.loc[table, "balanced_accuracy"])) <= 1e-6
        and abs(float(stored.loc[table, "balanced_accuracy"]) - receipt["balanced_accuracy"][table]) <= 1e-6
        and len(plotted) == int(stored.loc[table, "units"]) == int(float(means.loc[table, "units"]))
        and abs(float(stored.loc[table, "chance"]) - 1 / len(plotted)) <= 1e-6 and abs(float(means.loc[table, "chance"]) - 1 / len(plotted)) <= 1e-6
        for table in stored.index))

    # F and SuppFig7: exactly the tested strata with the medians of the summary, every stratum with its result or reason, and the ROC curve of every
    # tested stratum of a primary formula.
    tests = pd.read_parquet(folder / "signature_tests.parquet")
    keys = ["signature", "variant", "comparison", "study_accession", "platform", "tissue"]
    tested = tests[(tests["status"] == "tested") & (tests["variant"] == "primary")]
    plotted = panel("Fig7", "signatures")
    columns = {tuple(column) for column in policy["signature_columns"]}
    qc.check("signature_strata_are_exactly_the_tested_strata", same_rows(plotted, tested, keys, ["auroc"], tolerance=1e-6)
             and same_labels(plotted, tested, keys, ["independent", "derivation_role"])
             and set(map(tuple, plotted[keys[:3]].drop_duplicates().to_numpy())) == {column for column in columns if column[1] == "primary"})
    summary = pd.read_parquet(folder / "signature_summary.parquet")
    summary = summary[(summary["scope"] == "all_tissues") & (summary["variant"] == "primary") & summary["strata_set"].isin(["all", "independent"])]
    drawn = plotted.assign(auroc=plotted["auroc"].astype(float))
    agree = []
    for row in summary.itertuples(index=False):
        own = drawn[(drawn["signature"] == row.signature) & (drawn["variant"] == row.variant) & (drawn["comparison"] == row.comparison)]
        own = own[own["independent"] == "True"] if row.strata_set == "independent" else own
        agree.append(len(own) == row.tested and abs(float(own["auroc"].median()) - row.auroc_median) <= 1e-6)
    qc.check("signature_medians_equal_the_summary_and_the_strata", same_rows(panel("Fig7", "signature_medians"), summary,
                                                                             ["signature", "variant", "comparison", "strata_set"], ["tested", "auroc_median"],
                                                                             tolerance=1e-6) and len(agree) > 0 and all(agree))
    plotted = panel("SuppFig7", "heatmap")
    layout = plotted[["study_accession", "platform", "tissue", "in_scope_group", "material", "column", "label"]].drop_duplicates()
    qc.check("every_signature_stratum_equals_the_benchmark_stage", same_rows(plotted, tests, keys, ["positives", "negatives", "auroc", "auroc_low",
                                                                                                    "auroc_high"], tolerance=1e-6)
             and same_labels(plotted, tests, keys, ["status", "reason", "derivation_role", "independent", "in_scope_group"])
             and set(map(tuple, tests[keys[:3]].drop_duplicates().to_numpy())) == columns
             and not layout.duplicated(["study_accession", "platform", "tissue"]).any() and layout["label"].is_unique
             and ((layout["material"] == policy["material_order"][0]) == (layout["in_scope_group"] == "True")).all()
             and not layout.duplicated(["in_scope_group", "column"]).any())
    scores = pd.read_parquet(folder / "signature_scores.parquet")
    qc.check("signature_curves_recomputed_from_the_scores_enclose_the_tested_aurocs",
             curves_agree(panel("SuppFig7", "roc"), scores[scores["variant"] == "primary"], keys,
                          {tuple(getattr(row, key) for key in keys): float(row.auroc) for row in tested.itertuples(index=False)}))


class QC:
    def __init__(self):
        self.checks = []

    def check(self, name, passed, detail=None):
        self.checks.append(dict(name=name, passed=bool(passed), detail=detail))
        if not passed:
            raise ValueError(f"QC failed: {name}")


def verify(qc, out, config_path):
    config = V.read_json(config_path)
    manifest = V.read_json(out / "manifest.json")
    qc.check("configuration_is_the_one_of_the_figures", config["version"] == "1.0.0" and config["dpi"] == 300
             and manifest.get("configuration") == V.record(config_path))
    stages = manifest["stages"]
    qc.check("all_figure_scripts_complete", set(stages) == {str(number) for number in V.STAGES} and all(
        entry["status"] == "COMPLETE" and entry["figures"] == V.FIGURES[int(number)] for number, entry in stages.items()))
    for name in V.UPSTREAM:
        report = V.read_json(V.ROOT / "outputs" / name / "qc.json")
        qc.check(f"stage_quality_control_intact_{name}", report.get("status") == "PASS" and all(
            V.record(V.ROOT / path) == expected for path, expected in {**report["code"], **report["dependencies"]}.items()))
    bound = all(V.record(V.ROOT / path) == expected for entry in stages.values() for path, expected in {**entry["inputs"], **entry["upstream"], **entry["code"]}.items())
    qc.check("all_source_and_code_hashes", bound)

    outputs = {name: value for entry in stages.values() for name, value in entry["outputs"].items()}
    found = {path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()} - {"manifest.json", "figure_data.csv", "qc.json", "captions.md"}
    qc.check("exact_set_of_figure_files", found == set(outputs), sorted(found ^ set(outputs))[:10])
    figures = [name for names in V.FIGURES.values() for name in names]
    qc.check("every_figure_as_pdf_and_png", all(f"{name}.{extension}" in outputs for name in figures for extension in ("pdf", "png")))
    qc.check("every_figure_has_a_caption", set(config["captions"]) == set(figures) and all(str(text).strip() for text in config["captions"].values()))
    for name, expected in sorted(outputs.items()):
        path = out / name
        actual = V.record(path)
        fine = actual["sha256"] == expected["sha256"] and actual["size_bytes"] == expected["size_bytes"]
        if name.endswith(".pdf"):
            fine = fine and path.read_bytes()[:5] == b"%PDF-"
        else:
            with Image.open(path) as image:
                width, height = image.size
                dpi = image.info.get("dpi", (0, 0))
                grey = np.asarray(image.convert("L"))
            # Anything drawn on the outermost pixels was cut by the edge of the figure.
            frame = np.concatenate([grey[:EDGE].ravel(), grey[-EDGE:].ravel(), grey[:, :EDGE].ravel(), grey[:, -EDGE:].ravel()])
            fine = (fine and image.format == "PNG" and abs(dpi[0] - 300) < 1 and abs(dpi[1] - 300) < 1 and grey.min() < grey.max() and (frame == 255).all()
                    and abs(width - expected["width_inches"] * 300) <= 1 and abs(height - expected["height_inches"] * 300) <= 1)
        qc.check(f"file_{name}", fine)

    data = pd.read_csv(out / "figure_data.csv", dtype=str, keep_default_na=False)
    qc.check("figure_data_documents_every_figure", set(data["figure"]) == set(figures) and all(
        int((data["figure"].isin(entry["figures"])).sum()) == entry["data_rows"]
        and all(sorted(set(data.loc[data["figure"] == figure, "panel"])) == panels for figure, panels in entry["data_panels"].items()) for entry in stages.values()))

    samples = pd.read_parquet(V.ROOT / "data/processed/cohort/cohort_samples.parquet")
    ledger = pd.read_parquet(V.ROOT / "data/processed/cohort/records_without_values.parquet")
    records = V.read_json(V.ROOT / "outputs/cohort_curation/01.json")["records"]
    families = pd.read_parquet(V.ROOT / "data/processed/annotation/matrix_index.parquet", columns=["accession", "platform", "family"])
    distinct = samples[samples["is_primary_record"]]
    flow = data[(data["figure"] == "Fig1")].set_index("step")["count"].astype(float)
    # Fig1 B: the included records are those with expression values and those of the ledger, by its reasons; from the records with values on, every
    # count follows from the cohort table.
    qc.check("flow_counts_equal_the_cohort_tables", flow["included_sample_records"] == records["included"] == len(samples) + len(ledger)
             and flow["records_without_values"] == len(ledger) == records["without_values"]["total"]
             and {step[len("without_values_"):]: count for step, count in flow.items() if step.startswith("without_values_")}
             == ledger["reason"].value_counts().to_dict() == records["without_values"]["by_reason"]
             and flow["sample_records_with_values"] == len(samples) == records["with_values"]
             and flow["repeat_deposits"] == (~samples["is_primary_record"]).sum() == len(samples) - len(distinct)
             and flow["distinct_samples"] == len(distinct) and flow["one_sample_per_patient"] == samples["one_per_patient"].sum()
             and flow["identified_patients"] == distinct.loc[distinct["patient_identified"], "patient_key"].nunique()
             and flow["samples_without_identifier"] == (~distinct["patient_identified"]).sum()
             and flow["one_sample_per_patient"] == flow["identified_patients"] + flow["samples_without_identifier"]
             and flow["sepsis_before_onset"] == (distinct["sepsis_phase"] == "before_onset").sum()
             and all(flow[f"class_{name}"] == (distinct["sample_class"] == name).sum() for name in config["class_order"])
             and sum(flow[f"class_{name}"] for name in config["class_order"]) == flow["distinct_samples"]
             and flow["study_platforms_with_values"] == samples[["study_accession", "platform"]].drop_duplicates().shape[0]
             and flow["included_study_platforms"] == len(families))
    platforms = pd.read_parquet(V.ROOT / "data/processed/cohort/cohort_summary.parquet", columns=["study_accession", "platform", "release_date"]).merge(
        families.rename(columns={"accession": "study_accession"}), on=["study_accession", "platform"], validate="one_to_one")
    keyed = distinct.merge(platforms.assign(year=platforms["release_date"].str[:4].astype(int)), on=["study_accession", "platform"], how="left",
                           validate="many_to_one")
    by_year = keyed.groupby(["year", "family"]).size().rename("samples").reset_index()
    # figure_data.csv writes the year as a number that other panels made decimal, so the year is compared as a number; the years without samples of
    # a technology are drawn as zero.
    plotted = data[(data["figure"] == "Fig2") & (data["panel"] == "samples_by_year")]
    plotted = plotted.assign(year=plotted["year"].astype(float).astype(int))
    qc.check("samples_by_year_equal_the_cohort_table", same_rows(plotted[plotted["samples"].astype(float) > 0], by_year, ["year", "family"], ["samples"])
             and plotted["samples"].astype(float).sum() == len(distinct))
    material = {tissue: group for group, tissues in config["tissue_groups"].items() for tissue in tissues}
    composition = pd.crosstab(distinct["tissue"].map(material), distinct["sample_class"]).reindex(columns=config["class_order"], fill_value=0)
    plotted = data[(data["figure"] == "Fig2") & (data["panel"] == "composition_material")].set_index("group")[config["class_order"]].astype(float)
    qc.check("composition_equals_the_cohort_table", all(float(plotted.loc[group, name]) == float(composition.loc[group, name]) if group in composition.index
                                                        else float(plotted.loc[group, name]) == 0 for group in plotted.index for name in config["class_order"]))
    flows = data[(data["figure"] == "Fig2") & (data["panel"] == "flows")].assign(value=lambda frame: frame["value"].astype(float))
    first = flows[flows["step"] == "technology_to_material"].set_index(["source", "target"])["value"]
    second = flows[flows["step"] == "material_to_class"].set_index(["source", "target"])["value"]
    expected = keyed.assign(material=keyed["tissue"].map(material)).groupby(["family", "material"]).size().astype(float)
    qc.check("flows_equal_the_cohort_table", first.sort_index().equals(expected.rename_axis(["source", "target"]).sort_index().rename("value"))
             and all(second.get((group, name), 0.0) == float(composition.loc[group, name]) for group in composition.index for name in config["class_order"])
             and first.sum() == second.sum() == len(distinct))
    scored = pd.read_parquet(V.ROOT / "data/processed/cohort/sample_scores.parquet", columns=["study_accession", "sample_accession"]).drop_duplicates()
    scored = set(zip(scored["study_accession"], scored["sample_accession"]))
    recorded = pd.DataFrame({
        "sample_class": samples["sample_class"] != "not_recorded", "tissue": samples["tissue"] != "not_recorded", "sex": samples["sex"].isin(["female", "male"]),
        "age_group": samples["age_group"] != "not_recorded", "age_years": samples["age_years"].notna(), "patient": samples["patient_identified"],
        "time_point": samples["time_order"].notna(), "outcome": samples["outcome"] != "not_recorded", "septic_shock": samples["septic_shock"] != "not_recorded",
        "score": [key in scored for key in zip(samples["study_accession"], samples["sample_accession"])]}).astype(float)
    recorded["label"] = samples["study_accession"] + np.where(samples["platform"] != "", "  " + samples["platform"], "")
    shares = recorded.groupby("label").mean().join(recorded.groupby("label").size().rename("records")).reset_index()
    qc.check("recorded_variables_equal_the_cohort_table", list(config["variables"]) == list(recorded.columns[:-1]) and same_rows(
        data[(data["figure"] == "SuppFig2") & (data["panel"] == "coverage")], shares, ["label"], ["records", *config["variables"]], tolerance=1e-6))
    summary = pd.read_parquet(V.ROOT / "data/processed/annotation/gene_matrix_summary.parquet")
    plotted = data[(data["figure"] == "Fig4") & (data["panel"] == "genes_per_table")]
    qc.check("genes_per_table_equal_the_annotation_stage", same_rows(plotted, summary.rename(columns={"accession": "study_accession"}),
                                                                     ["study_accession", "platform"], ["genes", "protein_coding_genes"]))
    view = config["stacked_view"]
    chosen = summary[summary["value_scale"].isin(view["scales"]) & (summary["genes"] >= view["minimum_genes"])]
    plotted = data[(data["figure"] == "Fig4") & (data["panel"] == "stacked_view")]
    classes = samples.assign(table=samples["study_accession"] + "|" + samples["platform"]).set_index(["table", "sample_accession"])["sample_class"]
    keys = list(zip(plotted["table"], plotted["sample_accession"]))
    per_table = plotted.groupby("table").size()
    first = plotted.iloc[0] if len(plotted) else None
    qc.check("stacked_view_holds_the_samples_of_the_chosen_tables", first is not None and set(per_table.index) == set(chosen["accession"] + "|" + chosen["platform"])
             and all(per_table[f"{row.accession}|{row.platform}"] == row.samples for row in chosen.itertuples()) and len(set(keys)) == len(keys)
             and set(keys) <= set(classes.index) and (classes.loc[keys].to_numpy() == plotted["sample_class"].to_numpy()).all()
             and float(first["genes_used"]) == min(view["genes"], float(first["genes_shared"]))
             and components_agree(plotted["pc1"], plotted["pc2"], float(first["explained_1"]), float(first["explained_2"])))
    sex = pd.read_parquet(V.ROOT / "data/processed/validation/sex_check_matrices.parquet")
    plotted = data[(data["figure"] == "Fig5") & (data["panel"] == "sex_by_table")].assign(marked=lambda frame: (frame["marked"] == "True").astype(float))
    sex = sex.assign(label=[f"{study}  {platform}" if platform else study for study, platform in zip(sex["study_accession"], sex["platform"])],
                     marked=sex["recorded_sex_unreliable"].astype(float))
    sex = sex[(sex["concordant"] + sex["discordant"]) > 0]
    qc.check("sex_check_rows_equal_the_validation_stage", same_rows(plotted, sex, ["label"], ["concordant", "discordant", "marked"]))
    scores = pd.read_parquet(V.ROOT / "data/processed/validation/sepsis_score_by_study.parquet")
    scores = scores[scores["score_computed"]]
    plotted = data[(data["figure"] == "Fig5") & (data["panel"] == "score")]
    # The figure data hold the scores rounded to six decimals; an interval that the validation stage left empty must be empty in the figure data too.
    qc.check("score_rows_equal_the_validation_stage", same_rows(plotted, scores, ["study_accession", "platform", "tissue"],
                                                                ["auroc", "auroc_low", "auroc_high", "auroc_one_per_patient"], tolerance=1e-6))
    # Every drawn ROC curve encloses the area that the validation stage reports for its study platform and tissue.
    curves = data[(data["figure"] == "Fig5") & (data["panel"] == "roc")]
    areas = {stratum: float(np.trapz(part["tpr"].astype(float), part["fpr"].astype(float))) for stratum, part in curves.groupby("stratum")}
    reported = {f"{row.study_accession}|{row.platform}|{row.tissue}": float(row.auroc) for row in scores.itertuples(index=False)}
    qc.check("roc_curves_have_the_areas_of_the_validation_stage", set(areas) == set(reported) and all(abs(areas[key] - reported[key]) <= 1e-6 for key in areas))
    cells = pd.read_parquet(V.ROOT / "data/processed/validation/cell_markers_by_study.parquet")
    plotted = data[(data["figure"] == "Fig5") & (data["panel"] == "markers")]
    own = plotted[plotted["fraction"] == plotted["marker_set"]].set_index(["study_accession", "platform", "fraction"])["mean"].astype(float)
    qc.check("marker_means_equal_the_validation_stage", len(own) == len(cells) and all(
        abs(own.get((row.study_accession, row.platform, row.cell_type), np.nan) - row.mean_in_cell_type) <= 1e-6 for row in cells.itertuples(index=False)))
    units = data[(data["figure"] == "SuppFig3") & (data["panel"] == "splits")]
    matrices = pd.read_parquet(V.ROOT / "data/processed/validation/sex_check_matrices.parquet")
    qc.check("sex_agreement_of_the_units_adds_up_to_the_tables", units["concordant"].astype(float).sum() == matrices["concordant"].sum()
             and units["compared"].astype(float).sum() == (matrices["concordant"] + matrices["discordant"]).sum())
    pairs = pd.read_parquet(V.ROOT / "data/processed/validation/repeat_pairs.parquet")
    plotted = data[(data["figure"] == "Fig5") & (data["panel"] == "repeats")]
    qc.check("repeat_rows_equal_the_validation_stage", same_rows(plotted.assign(centred_correlation=plotted["centred_correlation"].replace("", "nan")), pairs,
                                                                 ["sample_id", "study_a", "platform_a", "study_b", "platform_b"], ["centred_correlation"]))
    quality = pd.read_parquet(V.ROOT / "data/processed/cohort/sample_quality.parquet")
    plotted = data[(data["figure"] == "SuppFig5") & (data["panel"] == "flags")].set_index("flag")
    qc.check("flag_counts_equal_the_cohort_stage", sorted(plotted.index) == sorted(config["flag_labels"])
             and (quality[list(config["flag_labels"])].any(axis=1) == quality["any_flag"]).all()
             and all(float(plotted.loc[name, "records"]) == quality[name].sum() for name in config["flag_labels"]))
    ledger = pd.read_parquet(V.ROOT / "data/processed/correction/correction_ledger.parquet")
    plotted = data[(data["figure"] == "Fig6") & (data["panel"] == "outcomes")]
    qc.check("record_outcomes_equal_the_correction_stage", plotted["records"].astype(float).sum() == len(ledger)
             and plotted.loc[plotted["outcome"].str.startswith("in_table_"), "records"].astype(float).sum() == (ledger["outcome"] == "in_corrected_table").sum())
    flows = data[(data["figure"] == "Fig6") & (data["panel"] == "ledger_flows") & (data["step"] == "group_to_outcome")]
    target = np.where(ledger["outcome"] == "in_corrected_table", "in_corrected_table", np.where(ledger["outcome"] == "corrected_for_checks_only",
                                                                                                "repeat_deposit", ledger["reason"]))
    expected = pd.DataFrame({"source": ledger["group"].to_numpy(), "target": target}).value_counts().astype(float).sort_index()
    drawn = flows.assign(value=flows["value"].astype(float)).set_index(["source", "target"])["value"].sort_index()
    qc.check("ledger_flows_equal_the_correction_stage", drawn.sum() == len(ledger) and drawn.to_dict() == expected.to_dict())
    test = pd.read_parquet(V.ROOT / "data/processed/correction/held_out_test.parquet")
    plotted = data[(data["figure"] == "Fig6") & (data["panel"] == "held_out")]
    qc.check("held_out_rows_equal_the_correction_stage", same_rows(plotted, test, ["group", "study_accession", "platform", "method"],
                                                                   ["location_error", "scale_error", "pool_dispersion", "sampling_floor"]))
    checks = pd.read_parquet(V.ROOT / "data/processed/correction/unit_checks.parquet")
    applied = checks[checks["applicable"]]
    plotted = data[(data["figure"] == "Fig6") & (data["panel"] == "mixing")]
    qc.check("mixing_rows_equal_the_correction_stage", same_rows(plotted, applied[applied["check"] == "study_mixing"],
                                                                 ["group", "study_accession", "platform", "sample_class"], ["before", "after"]))
    corrected = sum(len(pd.read_parquet(V.ROOT / f"data/processed/correction/corrected/{group}_samples.parquet", columns=["sample_id"])) for group in
                    config["correction"]["groups"])
    plotted = data[(data["figure"].isin(["Fig6", "SuppFig6"])) & (data["panel"].str.startswith("views_"))]
    qc.check("views_hold_every_corrected_sample", len(plotted) == corrected and plotted["sample_id"].is_unique)
    # The shares of panel D are read again from the plotted scores of every sample, and the scores against their stated shares of variance.
    factors, agree = data[(data["figure"] == "Fig6") & (data["panel"] == "factors")], []
    for row in factors.itertuples(index=False):
        points = plotted[plotted["panel"] == f"views_{row.group}"]
        own = points[points["unit_method"] == "healthy_anchored"] if row.factor == "sample_class_anchored" else points
        labels = own["sample_class" if row.factor == "sample_class_anchored" else row.factor]
        agree.append(float(row.samples) == len(own) and all(abs(explained_by(own[f"pc{number}_{row.stage}"], labels) - float(getattr(row, f"component_{number}")))
                                                          <= COMPONENTS for number in (1, 2))
                     and components_agree(points[f"pc1_{row.stage}"], points[f"pc2_{row.stage}"], float(row.variance_1), float(row.variance_2)))
    qc.check("factor_shares_equal_the_plotted_components", set(factors["group"]) == set(config["correction"]["groups"])
             and not factors.duplicated(["group", "stage", "factor"]).any() and len(agree) > 0 and all(agree))
    benchmark(qc, data, config)
    return manifest, dict(figures=len(figures), files=len(outputs), panels=sum(name.startswith("panels/") for name in outputs) // 2, figure_data_rows=len(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(V.HERE / "visualization_config.json"), help="Visualization configuration")
    parser.add_argument("--outdir", default="outputs/visualization", help="Figure directory inside outputs/")
    args = parser.parse_args()
    out = Path(args.outdir)
    out = (out if out.is_absolute() else V.ROOT / out).resolve()
    qc, failure, summary, manifest = QC(), None, None, None
    print("[08] Verifying the figures", flush=True)
    try:
        manifest, summary = verify(qc, out, Path(args.config).resolve())
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(name="verification_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version="1.0.0", status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(), required_checks=len(qc.checks),
                  failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks, exception=failure,
                  eligibility="FIGURES_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=summary)
    if passed:
        manifest["status"] = "VERIFIED"
        V.write_json(out / "manifest.json", manifest)
        # The captions belong outside the figures; they are written beside them for the README and the supplementary document.
        captions = V.read_json(Path(args.config).resolve())["captions"]
        figures = [name for names in V.FIGURES.values() for name in names]
        text = "# Figure captions\n\n" + "\n\n".join(f"**{name}.** {captions[name]}" for name in figures) + "\n"
        temporary = out / "captions.md.tmp"
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(out / "captions.md")
    V.write_json(out / "qc.json", report)
    print(f"\n[08 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(out / 'qc.json').relative_to(V.ROOT).as_posix()}", flush=True)
    if passed:
        print(f"    {summary['figures']} figures and {summary['panels']} panels as PDF and PNG; {summary['figure_data_rows']:,} rows of figure data", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
