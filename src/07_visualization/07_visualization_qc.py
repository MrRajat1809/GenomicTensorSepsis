"""Verify the figures of the resource: their provenance, the exact set of files, every image, and the plotted numbers against the stage tables."""
import argparse
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
           "encloses the area the validation stage reports; and that the principal components of the stacked view and of the corrected tables agree "
           "with their stated shares of variance and with the shares explained by technology, study and class. It does not judge whether a figure is "
           "well designed.")
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
    return manifest, dict(figures=len(figures), files=len(outputs), panels=sum(name.startswith("panels/") for name in outputs) // 2, figure_data_rows=len(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(V.HERE / "visualization_config.json"), help="Visualization configuration")
    parser.add_argument("--outdir", default="outputs/visualization", help="Figure directory inside outputs/")
    args = parser.parse_args()
    out = Path(args.outdir)
    out = (out if out.is_absolute() else V.ROOT / out).resolve()
    qc, failure, summary, manifest = QC(), None, None, None
    print("[07] Verifying the figures", flush=True)
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
    print(f"\n[07 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(out / 'qc.json').relative_to(V.ROOT).as_posix()}", flush=True)
    if passed:
        print(f"    {summary['figures']} figures and {summary['panels']} panels as PDF and PNG; {summary['figure_data_rows']:,} rows of figure data", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
