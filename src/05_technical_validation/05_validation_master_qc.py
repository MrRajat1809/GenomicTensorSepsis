"""Verify the index, the sex check, the comparison of repeated deposits and the biology check of the validation stage from their files."""
import importlib
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_index_validation_inputs")

RECORDED = ("female", "male")
MEANING = ("PASS means that every sample of the cohort table was indexed and tested exactly once; that the sex check, recomputed from the gene matrices, gives the "
           "recorded differences and that every call follows the stated rule without the recorded sex, split again whole or by tissue with the same limits, "
           "with every matrix whose recorded sex disagrees with its expression marked as the rule says; that every two records of one sample were compared "
           "once and a recomputed correlation agrees; and that every study platform and tissue with sepsis and healthy controls was tested, with the scores "
           "recomputed by the published formula, the areas under the curve recomputed from them over all samples and over one sample of every patient "
           "chosen again, an ordered interval exactly where both groups hold enough patients, and the role of every study in the derivation of the score. "
           "It verifies that the tests "
           "were carried out as described, not that every sample passed them: the findings are in the tables and the review sheets.")


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


def same_number(value, recorded):
    """Whether a recomputed number agrees with a recorded one: within 1e-6 where both are finite, and equal where either is missing or infinite; the
    receipt writes a missing or infinite number as null."""
    if recorded is None:
        return not np.isfinite(value)
    if isinstance(recorded, bool) or not isinstance(recorded, (int, float)):
        return False
    if np.isnan(value) or np.isnan(recorded):
        return bool(np.isnan(value) and np.isnan(recorded))
    if np.isinf(value) or np.isinf(recorded):
        return bool(value == recorded)
    return abs(value - recorded) <= 1e-6


def split_again(values, male, female, policy):
    """The split of the sex check, written again: two-means from the lowest and the highest value, then the limits on the separation, on the smaller
    group and on the Y-chromosome genes and XIST."""

    def apart(low, high):
        gap = high.mean() - low.mean()
        pooled = np.sqrt((((low - low.mean()) ** 2).sum() + ((high - high.mean()) ** 2).sum()) / max(len(low) + len(high) - 2, 1))
        return float(gap / pooled) if pooled > 0 else (np.inf if gap > 0 else -np.inf)

    result = dict(samples=len(values), status="too_few_samples", lower_group=0, higher_group=0, gap=np.nan, separation=np.nan, male_genes_separation=np.nan,
                  female_gene_separation=np.nan, midpoint=np.nan)
    if len(values) < policy["minimum_samples"]:
        return result, None
    low, high = values.min(), values.max()
    higher = np.zeros(len(values), dtype=bool)
    for _ in range(1000):
        assigned = np.abs(values - high) < np.abs(values - low)
        if not assigned.any() or assigned.all() or (assigned == higher).all():
            higher = assigned
            break
        higher, low, high = assigned, values[~assigned].mean(), values[assigned].mean()
    result["status"] = "one_group"
    if higher.all() or not higher.any():
        return result, None
    result.update(lower_group=int((~higher).sum()), higher_group=int(higher.sum()), gap=float(values[higher].mean() - values[~higher].mean()),
                  separation=apart(values[~higher], values[higher]), midpoint=float((values[higher].mean() + values[~higher].mean()) / 2),
                  male_genes_separation=apart(male[~higher], male[higher]))
    if female is not None:
        result["female_gene_separation"] = apart(-female[~higher], -female[higher])
    if not (result["separation"] >= policy["minimum_separation"]
            and min(result["lower_group"], result["higher_group"]) >= max(policy["minimum_group"], policy["minimum_group_share"] * len(values))):
        return result, None
    if not (result["male_genes_separation"] >= policy["component_separation"]
            and (female is None or result["female_gene_separation"] >= policy["component_separation"])):
        result["status"] = "components_not_separated"
        return result, None
    result["status"] = "tested"
    return result, higher


def provenance_checks(qc, state):
    state["start"] = C.begin(5)
    qc.check("all_stage_receipts_and_code_bindings", True)


def index_checks(qc, state):
    cohort = pd.read_parquet(C.COHORT / "cohort_samples.parquet")
    samples = pd.read_parquet(C.DATA / "validation_samples.parquet")
    matrices = pd.read_parquet(C.DATA / "validation_matrices.parquet")
    genes = pd.read_parquet(C.DATA / "validation_genes.parquet")
    qc.check("every_sample_indexed_once", len(samples) == len(cohort) and not samples.duplicated(C.KEYS).any()
             and set(zip(samples["study_accession"], samples["sample_accession"])) == set(zip(cohort["study_accession"], cohort["sample_accession"])))
    order = list(zip(samples["study_accession"], samples["sample_accession"]))
    source = cohort.set_index(C.KEYS).loc[order]
    qc.check("indexed_labels_agree_with_the_cohort_table", all((samples[name].to_numpy() == source[name].to_numpy()).all() for name in
                                                               ("sample_class", "sepsis_phase", "tissue", "sex", "sample_id", "records_of_sample", "patient_key"))
             and ((samples["time_order"].to_numpy() == source["time_order"].to_numpy()) | (samples["time_order"].isna().to_numpy()
                                                                                        & source["time_order"].isna().to_numpy())).all())
    qc.check("every_gene_matrix_listed_once_with_its_samples", not matrices.duplicated(["study_accession", "platform"]).any()
             and int(matrices["samples"].sum()) == len(samples) and all(C.project_path(path).is_file() for path in matrices["gene_matrix"]))
    reference = pd.read_parquet(C.ANNOTATION / "genes.parquet", columns=["hgnc_id", "symbol"]).set_index("hgnc_id")["symbol"]
    qc.check("every_named_gene_is_one_hgnc_gene", set(genes["hgnc_id"]) <= set(reference.index) and (genes["hgnc_id"].map(reference) == genes["symbol"]).all()
             and len(genes) == len(C.gene_sets()))
    state.update(samples=samples, matrices=matrices, genes=genes,
                 summary=dict(samples=len(samples), gene_matrices=len(matrices), genes_named=int(genes["hgnc_id"].nunique())))


def sex_checks(qc, state):
    policy = C.CONFIG["sex"]
    table = pd.read_parquet(C.DATA / "sex_check_samples.parquet")
    summaries = pd.read_parquet(C.DATA / "sex_check_matrices.parquet")
    samples, genes = state["samples"], state["genes"]
    qc.check("every_sample_has_one_sex_row", len(table) == len(samples) and not table.duplicated(C.KEYS).any())
    order = list(zip(table["study_accession"], table["sample_accession"]))
    qc.check("recorded_sex_is_that_of_the_cohort_table", (table["recorded_sex"].to_numpy() == samples.set_index(C.KEYS).loc[order, "sex"].to_numpy()).all())
    female = genes.loc[(genes["test"] == "sex") & (genes["role"] == "female"), "symbol"].iloc[0]
    male = set(genes.loc[(genes["test"] == "sex") & (genes["role"] == "male"), "symbol"])
    agree, worst, status = True, 0.0, {}
    # The differences are read again from the gene matrices by gene symbol, a second route to the same rows; a missing value must be missing in both.
    for matrix in state["matrices"].itertuples(index=False):
        own = table[(table["study_accession"] == matrix.study_accession) & (table["platform"] == matrix.platform)].set_index("sample_accession")
        if not matrix.comparable_scale:
            status[(matrix.study_accession, matrix.platform)] = "scale_not_comparable"
            agree = agree and own["difference"].isna().all()
            continue
        frame = pd.read_parquet(C.project_path(matrix.gene_matrix), filters=[("symbol", "in", sorted({female, *male}))])
        frame = C.log_values(frame.set_index("symbol").drop(columns="hgnc_id").astype("float64"), matrix.value_scale)
        read = [name for name in frame.index if name in male]
        if not read:
            status[(matrix.study_accession, matrix.platform)] = "genes_missing"
        again = (frame.loc[read].mean(axis=0, skipna=True) - (frame.loc[female] if female in frame.index else 0.0)).reindex(own.index) if read else \
            pd.Series(np.nan, index=own.index)
        both = again.notna() & own["difference"].notna()
        agree = agree and (again.isna() == own["difference"].isna()).all() and (np.abs(again[both] - own.loc[both, "difference"]) <= 1e-6).all()
        if both.any():
            worst = max(worst, float(np.abs(again[both] - own.loc[both, "difference"]).max()))
    qc.check("differences_agree_with_the_gene_matrices", agree, worst)
    # The split is made again from the recorded differences: by tissue where the rule says so, with the limits of the configuration.
    units = {(unit["study_accession"], unit["platform"], unit["tissue"]): unit for unit in C.read_json(C.receipt_path(2))["units"]}
    rows = summaries.set_index(["study_accession", "platform"])
    tissue_of = samples.set_index(C.KEYS)["tissue"]
    expected, found, follows = pd.Series("not_called", index=table.index), set(), True
    for (study, platform), own in table.groupby(["study_accession", "platform"], sort=True):
        if (study, platform) in status:
            continue
        known = own[own["difference"].notna()]
        tissues = tissue_of.loc[list(zip(known["study_accession"], known["sample_accession"]))].to_numpy()
        counts = pd.Series(tissues).value_counts()
        by_tissue = len(counts) > 1 and counts.min() >= policy["minimum_samples_per_tissue"]
        tested = []
        for tissue in (sorted(counts.index) if by_tissue else [""]):
            unit = known[tissues == tissue] if tissue else known
            result, higher = split_again(unit["difference"].to_numpy(), unit["male_genes_value"].to_numpy(),
                                         unit["female_gene_value"].to_numpy() if unit["female_gene_value"].notna().all() and len(unit) else None, policy)
            found.add((study, platform, tissue))
            recorded = units.get((study, platform, tissue), {})
            follows = follows and recorded.get("status") == result["status"] and all(
                recorded.get(name) == result[name] for name in ("samples", "lower_group", "higher_group")) and all(
                same_number(result[name], recorded.get(name, "missing")) for name in ("gap", "separation", "male_genes_separation", "female_gene_separation",
                                                                                     "midpoint"))
            if result["status"] == "tested":
                tested.append(result)
                near = np.abs(unit["difference"] - result["midpoint"]) < policy["ambiguous_zone"] * result["gap"]
                expected.loc[unit.index] = np.where(near, "ambiguous", np.where(higher, "male", "female"))
        row = rows.loc[(study, platform)]
        if by_tissue:
            status[(study, platform)] = "tested_by_tissue" if tested else "no_tissue_shows_two_sexes"
            follows = (follows and row["lower_group"] == sum(unit["lower_group"] for unit in tested) and row["higher_group"] == sum(unit["higher_group"] for unit in tested)
                       and row[["gap", "separation", "midpoint"]].isna().all())
        else:
            status[(study, platform)] = result["status"]
            follows = follows and row["lower_group"] == result["lower_group"] and row["higher_group"] == result["higher_group"] and all(
                same_number(result[name], row[name]) for name in ("gap", "separation", "midpoint"))
    qc.check("calls_follow_the_two_group_rule", follows and found == set(units) and (table["expression_sex"] == expected).all()
             and all(rows.loc[key, "status"] == value for key, value in status.items()) and len(status) == len(rows),
             {name: int(count) for name, count in summaries["status"].value_counts().items()})
    tested = summaries[summaries["status"].isin(["tested", "tested_by_tissue"])]
    keyed = list(zip(table["study_accession"], table["platform"]))
    compared_by_matrix = summaries["concordant"] + summaries["discordant"]
    expected = (compared_by_matrix >= policy["minimum_compared"]) & (summaries["concordant"] / compared_by_matrix.where(compared_by_matrix > 0) < policy["unreliable_below"])
    marks = summaries.set_index(["study_accession", "platform"])["recorded_sex_unreliable"]
    qc.check("unreliable_recorded_sex_follows_the_rule_and_reaches_every_sample", (summaries["recorded_sex_unreliable"] == expected).all()
             and (table["recorded_sex_unreliable"].to_numpy() == marks.loc[keyed].to_numpy()).all(),
             sorted(summaries.loc[summaries["recorded_sex_unreliable"], "study_accession"]))
    stated = table["recorded_sex"].isin(RECORDED)
    decided = table["expression_sex"].isin(RECORDED)
    qc.check("agreement_follows_the_call_and_the_record", (table.loc[~stated, "agreement"] == "sex_not_recorded").all()
             and (table.loc[stated & decided, "agreement"] == np.where(table.loc[stated & decided, "expression_sex"] == table.loc[stated & decided, "recorded_sex"],
                                                                      "concordant", "discordant")).all()
             and (table.loc[stated & ~decided, "agreement"] == table.loc[stated & ~decided, "expression_sex"]).all())
    qc.check("matrix_summaries_agree_with_the_samples", int(summaries["concordant"].sum()) == int((table["agreement"] == "concordant").sum())
             and int(summaries["discordant"].sum()) == int((table["agreement"] == "discordant").sum()) and int(summaries["samples"].sum()) == len(table))
    compared = int((table["agreement"].isin(["concordant", "discordant"])).sum())
    state["summary"].update(matrices_showing_two_sexes=len(tested), samples_compared_for_sex=compared,
                            sex_concordant=int((table["agreement"] == "concordant").sum()), sex_discordant=int((table["agreement"] == "discordant").sum()))


def repeat_checks(qc, state):
    pairs = pd.read_parquet(C.DATA / "repeat_pairs.parquet")
    samples = state["samples"]
    repeated = samples[samples["records_of_sample"] > 1]
    expected = int(sum(count * (count - 1) // 2 for count in repeated.groupby("sample_id").size()))
    of_record = samples.set_index(C.KEYS)["sample_id"]
    qc.check("every_two_records_of_one_sample_compared_once", len(pairs) == expected
             and not pairs.duplicated(["study_a", "sample_a", "study_b", "sample_b"]).any()
             and (of_record.loc[list(zip(pairs["study_a"], pairs["sample_a"]))].to_numpy() == pairs["sample_id"].to_numpy()).all()
             and (of_record.loc[list(zip(pairs["study_b"], pairs["sample_b"]))].to_numpy() == pairs["sample_id"].to_numpy()).all())
    qc.check("correlations_are_valid", pairs["spearman"].between(-1 - 1e-9, 1 + 1e-9).all() and pairs["shared_genes"].gt(0).all()
             and (pairs.loc[pairs["most_similar"], ["rank_among_a", "rank_among_b"]] == 1).all().all())
    matrices = state["matrices"].set_index(["study_accession", "platform"])["gene_matrix"]
    worst, counted = 0.0, True
    # One pair of every two matrices is computed again with the Spearman correlation of pandas; a missing correlation fails.
    for _, pair in pairs.drop_duplicates(["study_a", "platform_a", "study_b", "platform_b"]).iterrows():
        first = pd.read_parquet(C.project_path(matrices.loc[(pair["study_a"], pair["platform_a"])])).set_index("hgnc_id").drop(columns="symbol").dropna()
        second = pd.read_parquet(C.project_path(matrices.loc[(pair["study_b"], pair["platform_b"])])).set_index("hgnc_id").drop(columns="symbol").dropna()
        shared = first.index.intersection(second.index)
        again = first.loc[shared, pair["sample_a"]].astype("float64").corr(second.loc[shared, pair["sample_b"]].astype("float64"), method="spearman")
        distance = abs(float(again) - pair["spearman"])
        counted = counted and len(shared) == pair["shared_genes"] and distance < 1e-6
        worst = max(worst, distance) if np.isfinite(distance) else worst
    qc.check("recomputed_correlations_and_shared_genes_agree", counted, worst)
    state["summary"].update(pairs_of_repeated_records=len(pairs), lowest_spearman_of_a_pair=round(float(pairs["spearman"].min()), 6),
                            median_spearman_of_a_pair=round(float(pairs["spearman"].median()), 6),
                            pairs_whose_linked_record_is_the_most_similar=int(pairs["most_similar"].sum()))


def biology_checks(qc, state):
    policy = C.CONFIG["biology"]
    strata = pd.read_parquet(C.DATA / "sepsis_score_by_study.parquet")
    rows = pd.read_parquet(C.DATA / "sepsis_genes_by_study.parquet")
    scores = pd.read_parquet(C.DATA / "sepsis_score_samples.parquet")
    cells = pd.read_parquet(C.DATA / "cell_markers_by_study.parquet")
    samples = state["samples"].merge(state["matrices"][["study_accession", "platform", "comparable_scale"]], on=["study_accession", "platform"])
    in_group = samples["sample_class"].eq("healthy_control") | (samples["sample_class"].eq("sepsis") & samples["sepsis_phase"].isin(policy["sepsis_phases"]))
    groups = samples[samples["comparable_scale"] & in_group]
    counts = groups.groupby(["study_accession", "platform", "tissue", "sample_class"]).size().unstack(fill_value=0).reindex(columns=["sepsis", "healthy_control"],
                                                                                                                         fill_value=0)
    eligible = counts[counts.min(axis=1) >= policy["minimum_group"]]
    keyed = strata.set_index(["study_accession", "platform", "tissue"])
    qc.check("every_eligible_study_platform_and_tissue_tested", set(keyed.index) == set(eligible.index)
             and (keyed.loc[eligible.index, "sepsis"].to_numpy() == eligible["sepsis"].to_numpy()).all()
             and (keyed.loc[eligible.index, "healthy_control"].to_numpy() == eligible["healthy_control"].to_numpy()).all())
    tested = rows[rows["present"]]
    expected = np.where(tested["role"] == "higher_in_sepsis", tested["difference"] > 0, tested["difference"] < 0)
    qc.check("directions_follow_the_differences", (tested["as_published"].to_numpy() == expected).all() and not rows.loc[~rows["present"], "as_published"].any()
             and len(rows) == len(strata) * int(state["genes"]["test"].isin(["sepsis_score", "single_gene"]).sum())
             and np.allclose(tested["difference"], tested["mean_sepsis"] - tested["mean_healthy_control"]))
    def pairwise(cases, others):
        if not len(cases) or not len(others):
            return np.nan
        return float(((cases[:, None] > others[None, :]).sum() + 0.5 * (cases[:, None] == others[None, :]).sum()) / (len(cases) * len(others)))

    agree, single, worst = True, True, 0.0
    time_of = state["samples"].set_index(C.KEYS)["time_order"]
    # The area under the curve is counted again from the sample scores, pair by pair, over all samples and over the one sample of every patient that
    # is chosen again here: within the stratum, the earliest by its time order, then the lowest accession.
    for key, part in scores.groupby(["study_accession", "platform", "tissue"]):
        row = keyed.loc[key]
        sepsis, healthy = (part["sample_class"] == "sepsis").to_numpy(), (part["sample_class"] == "healthy_control").to_numpy()
        again = pairwise(part.loc[sepsis, "score"].to_numpy(), part.loc[healthy, "score"].to_numpy())
        agree = agree and abs(again - float(row["auroc"])) < 1e-9
        worst = max(worst, abs(again - float(row["auroc"]))) if np.isfinite(again - float(row["auroc"])) else worst
        timing = time_of.loc[list(zip(part["study_accession"], part["sample_accession"]))].to_numpy(dtype="float64")
        ranked = sorted(zip(part["patient"], np.isnan(timing), np.where(np.isnan(timing), 0.0, timing), part["sample_accession"].map(C.accession_order),
                            range(len(part))))
        first = {}
        for patient, *_, position in ranked:
            first.setdefault(patient, position)
        marked = np.isin(np.arange(len(part)), list(first.values()))
        patients = part.loc[sepsis, "patient"].nunique(), part.loc[healthy, "patient"].nunique()
        expected = (pairwise(part.loc[sepsis & marked, "score"].to_numpy(), part.loc[healthy & marked, "score"].to_numpy())
                    if min((sepsis & marked).sum(), (healthy & marked).sum()) >= 2 else np.nan)
        interval = row[["auroc_low", "auroc_high"]].notna().all() if min(patients) >= policy["minimum_group"] else row[["auroc_low", "auroc_high"]].isna().all()
        single = (single and interval and (part["one_per_patient"].to_numpy() == marked).all()
                  and row["patients_sepsis"] == patients[0] == (sepsis & marked).sum() == row["sepsis_one_per_patient"]
                  and row["patients_healthy_control"] == patients[1] == (healthy & marked).sum() == row["healthy_control_one_per_patient"]
                  and np.isnan(expected) == np.isnan(row["auroc_one_per_patient"]) and (np.isnan(expected) or abs(expected - row["auroc_one_per_patient"]) < 1e-9))
    scored = strata[strata["score_computed"]]
    qc.check("areas_under_the_curve_agree_with_the_sample_scores", agree and scores.groupby(["study_accession", "platform", "tissue"]).ngroups == len(scored)
             and scored["auroc"].between(0, 1).all(), worst)
    enough = scored[["patients_sepsis", "patients_healthy_control"]].min(axis=1) >= policy["minimum_group"]
    qc.check("intervals_and_one_per_patient_areas_follow_the_patients", single and (scored["auroc_low"].notna() == enough).all()
             and (scored["auroc_high"].notna() == enough).all() and (scored.loc[enough, "auroc_low"] <= scored.loc[enough, "auroc_high"] + 1e-12).all()
             and scored.loc[enough, ["auroc_low", "auroc_high"]].stack().between(0, 1).all() and (scored["patients_sepsis"] <= scored["sepsis"]).all(),
             dict(with_an_interval=int(enough.sum()), without=int((~enough).sum())))
    derivation = policy["score_derivation"]
    role_of = {**{study: "validation" for study in derivation["validation"]}, **{study: "discovery" for study in derivation["discovery"]}}
    qc.check("derivation_roles_follow_the_configuration", (strata["score_derivation"] == strata["study_accession"].map(lambda study: role_of.get(study, ""))).all())
    # The score is computed again by gene symbol, a second route to the same rows, with the published formula written out once more.
    matrices = state["matrices"].set_index(["study_accession", "platform"])
    genes = state["genes"]
    up = set(genes.loc[(genes["test"] == "sepsis_score") & (genes["role"] == "higher_in_sepsis"), "symbol"])
    down = set(genes.loc[(genes["test"] == "sepsis_score") & (genes["role"] == "lower_in_sepsis"), "symbol"])
    largest, follows = 0.0, True
    for (study, platform, tissue), part in scores.groupby(["study_accession", "platform", "tissue"]):
        matrix = matrices.loc[(study, platform)]
        whole = C.log_values(pd.read_parquet(C.project_path(matrix["gene_matrix"])).set_index("symbol").drop(columns="hgnc_id").astype("float64"), matrix["value_scale"])
        lowest = float(whole.min().min())
        raised = whole[part["sample_accession"]] + (1.0 - lowest if lowest < 1 else 0.0)
        again = (np.exp(np.log(raised.loc[raised.index.isin(up)]).mean(axis=0, skipna=True))
                 - np.exp(np.log(raised.loc[raised.index.isin(down)]).mean(axis=0, skipna=True)))
        distance = np.abs(again.to_numpy() - part["score"].to_numpy())
        follows = follows and len(distance) == len(part) and bool((distance < 1e-6).all())
        largest = max(largest, float(distance.max())) if np.isfinite(distance).all() else largest
    qc.check("scores_follow_the_published_formula", follows, largest)
    qc.check("cell_marker_comparisons_are_valid", cells["auroc"].between(0, 1).all() and (cells["samples"] >= policy["minimum_cells"]).all()
             and (cells["other_samples"] >= policy["minimum_cells"]).all())
    state["summary"].update(study_platforms_and_tissues_with_sepsis_and_healthy_controls=len(strata), gene_comparisons=len(tested),
                            gene_comparisons_as_published=int(tested["as_published"].sum()), median_auroc_of_the_score=round(float(scored["auroc"].median()), 6),
                            lowest_auroc_of_the_score=round(float(scored["auroc"].min()), 6), cell_marker_comparisons=len(cells))


def total_checks(qc, state):
    sex, repeats, biology = (C.read_json(C.receipt_path(number)) for number in (2, 3, 4))
    summary = state["summary"]
    qc.check("receipt_totals_agree_with_the_tables", sex["samples"]["concordant"] == summary["sex_concordant"] and sex["samples"]["discordant"] == summary["sex_discordant"]
             and repeats["pairs"]["total"] == summary["pairs_of_repeated_records"] and biology["strata"]["total"] == len(pd.read_parquet(
                 C.DATA / "sepsis_score_by_study.parquet")) and biology["genes"]["as_published"] == summary["gene_comparisons_as_published"])


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/6 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/6 Index against the cohort table and the gene matrices", lambda: index_checks(qc, state))
        qc.phase("3/6 Sex check against the gene matrices", lambda: sex_checks(qc, state))
        qc.phase("4/6 Repeated deposits", lambda: repeat_checks(qc, state))
        qc.phase("5/6 Known biology", lambda: biology_checks(qc, state))
        qc.phase("6/6 Stage receipts and totals", lambda: total_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, 5) if C.receipt_path(number).exists()]),
                  eligibility="VALIDATION_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[05 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    Sex: {summary['sex_concordant']:,} of {summary['samples_compared_for_sex']:,} compared samples concordant; repeated deposits: "
              f"{summary['pairs_of_repeated_records']:,} pairs, median Spearman {summary['median_spearman_of_a_pair']:.4f}; sepsis score: median AUROC "
              f"{summary['median_auroc_of_the_score']:.3f} in {summary['study_platforms_and_tissues_with_sepsis_and_healthy_controls']} study platforms and tissues",
              flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
