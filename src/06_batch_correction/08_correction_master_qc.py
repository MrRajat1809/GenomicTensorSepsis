"""Verify the selection, the alignment, the corrections, the choice of methods and the corrected tables of the batch correction stage from their files."""
import importlib
import math
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_select_correction_inputs")

MEANING = ("PASS means that every sample record of the cohort has one outcome with its reason; that every record entered the group of its recorded sample type "
           "and scale; that the aligned values of a sample, rebuilt from its gene matrix, equal the stored ones and keep the order of its genes; that the sex "
           "effect, derived again as the precision-weighted difference between the sexes within cells of unit, class and age group, equals the stored one, and "
           "that no anchor standard deviation lies below its floor; that every anchored unit that can be held out was held out by every method that applies "
           "to it, and that every comparison of two methods mirrors its reverse; that every unit was kept only when its correction passed every applicable "
           "check, and a unit without healthy controls only when its study-mixing check applied, and that every pass follows from its stated values and "
           "limits; and that the corrected tables, rebuilt from the aligned values, the study effects and the model, equal the stored ones, with every sample "
           "once and its metadata equal to the cohort table. It verifies that the correction was carried out and checked as described; how well it removed "
           "the differences between studies is in the check, test and comparison tables.")


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
        print(f"[08] {name}", flush=True)
        action()


def largest(difference, worst):
    """The larger of the worst absolute difference so far and the largest one here, infinite where one is not a number: the detail of a check, which
    passes or fails on its own all-values comparison."""
    difference = np.abs(np.asarray(difference, dtype="float64"))
    if not np.isfinite(difference).all():
        return np.inf
    return max(worst, float(difference.max())) if difference.size else worst


def shown(value):
    return value if np.isfinite(value) else "not a finite number"


def provenance_checks(qc, state):
    state["start"] = C.begin(8)
    qc.check("all_stage_receipts_and_code_bindings", True)


def selection_checks(qc, state):
    inputs = C.CONFIG["inputs"]
    cohort = pd.read_parquet(C.COHORT / "cohort_samples.parquet")
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    qc.check("every_cohort_record_once", len(records) == len(cohort) and not records.duplicated(C.KEYS).any()
             and set(map(tuple, records[C.KEYS].to_numpy())) == set(map(tuple, cohort[C.KEYS].to_numpy())))
    groups = {tissue: name for name, group in C.CONFIG["groups"].items() for tissue in group["tissues"]}
    taken = records["tissue"].isin(groups) & records["value_scale"].isin(inputs["value_scales"])
    qc.check("groups_follow_sample_type_and_scale", (records.loc[taken & (records["reason"] != "study_has_too_few_samples"), "group"]
                                                      == records.loc[taken & (records["reason"] != "study_has_too_few_samples"), "tissue"].map(groups)).all()
             and (records.loc[~taken, "status"] == "left_out").all() and (records.loc[~taken, "reason"] != "").all())
    inside = records[records["status"] == "in_unit"]
    counted = inside.groupby(C.UNIT).agg(records=("sample_accession", "size"), estimation=("estimation", "sum")).reset_index()
    healthy = inside[inside["estimation"] & (inside["sample_class"] == C.CONFIG["anchored"]["anchor_class"])].groupby(C.UNIT).size()
    active = units[units["status"] == "in_unit"].set_index(C.UNIT)
    keyed = counted.set_index(C.UNIT)
    qc.check("unit_counts_agree_with_the_records", set(keyed.index) == set(active.index)
             and (keyed.loc[active.index, "records"].to_numpy() == active["records"].to_numpy()).all()
             and (keyed.loc[active.index, "estimation"].to_numpy() == active["estimation_samples"].to_numpy()).all()
             and (healthy.reindex(active.index, fill_value=0).to_numpy() == active["healthy_anchors"].to_numpy()).all())
    qc.check("tiers_and_sizes_follow_the_limits", ((active["healthy_anchors"] >= inputs["minimum_healthy_controls"]) == (active["tier"] == "anchored")).all()
             and (active["estimation_samples"] >= inputs["minimum_unit_samples"]).all()
             and (units.loc[units["status"] == "left_out", "estimation_samples"] < inputs["minimum_unit_samples"]).all())
    qc.check("no_anchored_unit_is_of_newborns", not (active["neonatal"] & (active["tier"] == "anchored")).any())
    # An estimation sample is an unflagged primary record (or any unflagged record of a unit for checks only), of a phase that stands for its class, and the
    # first sample of its patient within the unit by time order and accession, in natural order.
    eligible = (~inside["any_flag"] & (inside["is_primary_record"] | (inside["unit_role"] == "check_only"))
                & ~inside["sepsis_phase"].isin(inputs["phases_not_estimated"]))
    patient = inside["patient_key"].where(inside["patient_key"] != "", inside["sample_id"])
    first = set()
    for _, part in inside[eligible].assign(patient=patient).groupby([*C.UNIT, "patient"]):
        timed = part[part["time_order"].notna()]
        candidates = timed if len(timed) else part
        first.add(min(zip(candidates["time_order"].fillna(0.0), candidates["sample_accession"].map(C.accession_order), candidates.index))[2])
    qc.check("estimation_samples_are_the_first_unflagged_primary_record_of_every_patient", (inside["estimation"] == inside.index.isin(first)).all()
             and not records.loc[records["status"] == "left_out", "estimation"].any())
    qc.check("sex_used_is_stated", records.loc[records["status"] == "in_unit", "female"].between(0, 1).all()
             and records.loc[records["sex_used"] == "expression", "female"].isin([0, 1]).all()
             and not (records["recorded_sex_unreliable"] & (records["sex_used"] == "recorded")).any())
    rule = all(((part["in_group_set"]) == ((part["units_with_gene"] >= math.ceil(inputs["minimum_gene_share"] * part["merged_units"].iloc[0] - 1e-9))
                                            & (part["anchored_units_with_gene"] >= inputs["minimum_reference_units"]))).all() for _, part in genes.groupby("group"))
    qc.check("group_genes_follow_the_share_rule", rule)
    state.update(cohort=cohort, records=records, units=units, genes=genes, summary={})


def alignment_checks(qc, state):
    records, units, genes = state["records"], state["units"], state["genes"]
    reference = pd.read_parquet(C.DATA / "alignment_reference.parquet")
    points = C.CONFIG["alignment"]["grid_points"]
    expected = {C.aligned_path(row.group, row.study_accession, row.platform) for row in units[units["status"] == "in_unit"].itertuples(index=False)}
    found = set(C.ALIGNED.glob("*/*.parquet"))
    qc.check("exactly_one_aligned_file_per_unit", found == expected, sorted(path.name for path in found ^ expected)[:5])
    tolerance, correlation = C.CONFIG["quality_control"]["alignment_tolerance"], C.CONFIG["alignment"]["minimum_rank_correlation"]
    worst, rebuilt = 0.0, True
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        own_reference = reference[reference["group"] == group]
        qc.check(f"reference_of_{group}_is_ordered", len(own_reference) == points and own_reference["value"].is_monotonic_increasing)
        grid = own_reference["position"].to_numpy()
        for unit in units[(units["group"] == group) & (units["status"] == "in_unit")].itertuples(index=False):
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            frame = pd.read_parquet(C.aligned_path(group, unit.study_accession, unit.platform))
            qc.check(f"aligned_layout_{group}_{C.unit_name(unit.study_accession, unit.platform)}", list(frame["hgnc_id"]) == order
                     and list(frame.columns[1:]) == list(own["sample_accession"]))
            # One sample of every unit is aligned again from its gene matrix.
            column = own["sample_accession"].iloc[len(own) // 2]
            raw = C.read_values(C.project_path(unit.gene_matrix), [column], order)[:, 0]
            again = C.align_sample(raw, grid, own_reference["value"].to_numpy())
            stored = frame[column].to_numpy(dtype="float64")
            present = ~np.isnan(raw)
            difference = again[present] - stored[present]
            worst = largest(difference, worst)
            rebuilt = (rebuilt and bool((np.isnan(stored) == ~present).all()) and bool((np.abs(difference) <= tolerance).all())
                       and pd.Series(raw[present]).corr(pd.Series(stored[present]), method="spearman") > correlation)
    qc.check("aligned_values_rebuilt_from_the_gene_matrices", rebuilt, shown(worst))


def sex_slope(group, units, records, order):
    """The sex effect derived again: in every cell of unit, sample class and age group, the difference between the means of the females and of the males of
    the estimation samples of known sex, weighted over the cells by its precision nf nm / (nf + nm); and the number of units with any weight."""
    numerator, denominator, informative = np.zeros(len(order)), np.zeros(len(order)), np.zeros(len(order), dtype=int)
    for unit in units[(units["group"] == group) & (units["status"] == "in_unit") & (units["role"] == "merged")].itertuples(index=False):
        own = C.unit_records(records, group, unit.study_accession, unit.platform)
        own = own[own["estimation"] & own["sex_used"].isin(["expression", "recorded"])].reset_index(drop=True)
        if own.empty:
            continue
        values = C.read_aligned(group, unit.study_accession, unit.platform, own["sample_accession"])
        weight_of_unit = np.zeros(len(order))
        for _, part in own.groupby(["sample_class", "age_group"], dropna=False):
            females, males = (C.moments(values[:, part.index[part["female"] == sex]]) for sex in (1, 0))
            weight = np.divide(females[0] * males[0], females[0] + males[0], out=np.zeros(len(order)), where=(females[0] > 0) & (males[0] > 0))
            numerator += np.where(weight > 0, weight * (females[1] - males[1]), 0.0)
            denominator += weight
            weight_of_unit += weight
        informative += (weight_of_unit > 0).astype(int)
    return np.divide(numerator, denominator, out=np.zeros(len(order)), where=denominator > 0), informative


def held_out_methods(training, moments_of):
    """The methods of the held-out test that apply to every anchored unit with primary records, derived again from the rules of stage 05: none where the unit
    has too few samples besides its healthy controls or too few units remain to learn from."""
    inputs, policy, anchor = C.CONFIG["inputs"], C.CONFIG["unanchored"], C.CONFIG["anchored"]["anchor_class"]
    visible = {key: {cell: part for cell, part in moments_of[key].items() if cell[0] != anchor} for key in training}
    count = {key: int(np.nanmax(C.combine(cells.values())[0])) if cells else 0 for key, cells in visible.items()}
    expected = {}
    for key in training:
        rest = [other for other in training if other != key]
        if count[key] < inputs["minimum_unit_samples"] or len(rest) < inputs["minimum_reference_units"]:
            expected[key] = set()
            continue
        expected[key] = {"none", "plain_combat"}
        if sum(count[other] >= inputs["minimum_unit_samples"] for other in rest) >= C.CONFIG["flat_genes"]["minimum_training_units"]:
            expected[key].add("flat_genes")
        mean, deviation, _ = C.anchor_model([C.gather(moments_of[other], [anchor]) for other in rest])
        effects = {other: C.anchored_effects(C.gather(moments_of[other], [anchor]), mean, deviation) for other in rest}
        matched, freedom = C.matched_samples(C.matched_cells(visible[key], C.class_reference(moments_of, rest, mean, deviation, effects), excluded=[anchor]))
        if matched >= policy["minimum_matched"] and freedom >= policy["minimum_matched"] - 1:
            expected[key].add("matched_groups")
    return expected


def comparison_right(held, comparison, tolerance):
    """Every two methods compared in both orders over the units where both apply, with the count, median and signed-rank test of their differences, and an
    interval that mirrors that of the reverse order."""
    from scipy.stats import wilcoxon

    minimum = C.CONFIG["unanchored"]["comparison"]["minimum_units"]
    keyed = comparison.set_index(["group", "method", "baseline"])
    expected, right = set(), keyed.index.is_unique
    for group, part in held.groupby("group"):
        wide = part.pivot_table(index=["study_accession", "platform"], columns="method", values="location_error")
        for method in wide.columns:
            for baseline in [name for name in wide.columns if name != method]:
                both = wide[[method, baseline]].dropna()
                if len(both) < minimum:
                    continue
                expected.add((group, method, baseline))
                if not right or (group, method, baseline) not in keyed.index or (group, baseline, method) not in keyed.index:
                    right = False
                    continue
                row, reverse = keyed.loc[(group, method, baseline)], keyed.loc[(group, baseline, method)]
                difference = (both[method] - both[baseline]).to_numpy()
                right = (row["units"] == len(both) and row["method_closer"] == int((difference < 0).sum())
                         and abs(row["median_difference"] - float(np.median(difference))) <= tolerance
                         and abs(row["wilcoxon_p"] - float(wilcoxon(difference).pvalue)) <= tolerance
                         and row["difference_low"] == -reverse["difference_high"] and row["difference_high"] == -reverse["difference_low"])
    return bool(right) and set(keyed.index) == expected


def correction_checks(qc, state):
    records, units, genes = state["records"], state["units"], state["genes"]
    policy, anchor = C.CONFIG["flat_genes"], C.CONFIG["anchored"]["anchor_class"]
    flat = pd.read_parquet(C.DATA / "flat_genes.parquet")
    chosen = flat[flat["selected"]].merge(genes[["group", "hgnc_id", "check_role", "in_group_set"]], on=["group", "hgnc_id"])
    qc.check("flat_genes_follow_the_rules", (chosen["check_role"] == "").all() and chosen["in_group_set"].all()
             and not chosen["chromosome"].isin(policy["excluded_chromosomes"]).any()
             and all(len(part) == min(policy["count"], int((flat["group"] == group).sum())) for group, part in chosen.groupby("group")))
    anchored = pd.read_parquet(C.DATA / "anchored_effects.parquet")
    candidates = pd.read_parquet(C.DATA / "candidate_effects.parquet")
    fits = pd.read_parquet(C.DATA / "candidate_fits.parquet")
    qc.check("every_scale_is_positive", (anchored["scale"].dropna() > 0).all() and (candidates["scale"].dropna() > 0).all())
    quality = C.CONFIG["quality_control"]
    stored = pd.read_parquet(C.DATA / "correction_model.parquet")
    worst, settled, sex_worst, sex_right, floored, expected = 0.0, True, 0.0, True, True, {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        model = C.model_arrays(group, order)
        moments_of = C.unit_moments(group, units, records, model["sex_effect"])
        own = units[(units["group"] == group) & (units["status"] == "in_unit") & (units["tier"] == "anchored") & (units["role"] == "merged")]
        centres = []
        for unit in own.itertuples(index=False):
            part = anchored[(anchored["group"] == group) & (anchored["study_accession"] == unit.study_accession) & (anchored["platform"] == unit.platform)]
            part = part.set_index("hgnc_id").reindex(order)
            centres.append(C.corrected_moments(C.gather(moments_of[(unit.study_accession, unit.platform)], [anchor]), model["anchor_mean"], model["anchor_sd"],
                                               part["location"].to_numpy(), part["scale"].to_numpy())[1])
        median = float(np.nanmedian(np.abs(np.vstack(centres))))
        settled = settled and median <= quality["anchor_centre_limit"]
        worst = largest([median], worst)
        slope, informative = sex_slope(group, units, records, order)
        sex_units = stored[stored["group"] == group].set_index("hgnc_id").reindex(order)["sex_units"].to_numpy()
        sex_right = sex_right and bool((np.abs(slope - model["sex_effect"]) <= quality["model_tolerance"]).all()) and bool((informative == sex_units).all())
        sex_worst = largest(slope - model["sex_effect"], sex_worst)
        deviation = model["anchor_sd"][np.isfinite(model["anchor_sd"])]
        floored = floored and deviation.size > 0 and bool((deviation >= C.CONFIG["anchored"]["floor_fraction"] * float(np.median(deviation))
                                                           - quality["model_tolerance"]).all())
        expected.update({(group, *key): found for key, found in held_out_methods([(unit.study_accession, unit.platform) for unit in own.itertuples(index=False)],
                                                                                   moments_of).items()})
    # The anchored correction sets the location to the mean of the healthy controls, so this check is nearly true by construction; it guards the arithmetic.
    qc.check("healthy_controls_sit_on_the_anchor_after_correction", settled, shown(worst))
    qc.check("sex_effect_is_the_weighted_difference_between_the_sexes_within_cells", sex_right, shown(sex_worst))
    qc.check("no_anchor_standard_deviation_below_its_floor", floored)
    active = units[units["status"] == "in_unit"]
    unanchored = active[active["tier"] == "unanchored"]
    keyed = fits.set_index(C.UNIT)
    neonatal = unanchored.set_index(C.UNIT)["neonatal"]
    qc.check("every_unanchored_unit_corrected_by_matched_groups_or_explained", (fits["method"] == "matched_groups").all() and keyed.index.is_unique
             and set(keyed.index) == set(neonatal.index) and not keyed.loc[neonatal.index[neonatal], "applicable"].any()
             and (keyed.loc[neonatal.index[neonatal], "reason"] == "newborn_study_not_corrected").all()
             and set(fits.loc[~fits["applicable"], "reason"]) <= {"newborn_study_not_corrected", "no_class_to_line_up"})
    built = set(map(tuple, candidates[["group", "study_accession", "platform", "method"]].drop_duplicates().to_numpy()))
    qc.check("candidate_effects_match_the_fits", built == set(fits.loc[fits["applicable"], ["group", "study_accession", "platform", "method"]].apply(tuple, axis=1)))
    held = pd.read_parquet(C.DATA / "held_out_test.parquet")
    qc.check("held_out_rows_are_of_anchored_units_and_valid", not held.duplicated([*C.UNIT, "method"]).any() and set(map(tuple, held[C.UNIT].to_numpy())) <= set(expected)
             and set(held["method"]) <= {"none", "plain_combat", "flat_genes", "matched_groups"} and held["location_error"].ge(0).all()
             and held["scale_error"].ge(0).all())
    found = {key: set(part["method"]) for key, part in held.groupby(C.UNIT)}
    missed = sorted(f"{key[1]} {key[2]}" for key, methods in expected.items() if found.get(key, set()) != methods)
    qc.check("every_anchored_unit_held_out_by_every_method_that_applies", not missed, missed[:5])
    comparison = pd.read_parquet(C.DATA / "held_out_comparison.parquet")
    qc.check("every_comparison_of_two_methods_follows_from_the_test_and_mirrors_its_reverse", comparison_right(held, comparison, quality["model_tolerance"]))
    state["summary"].update(held_out_rows=len(held), held_out_units=len(found), median_location_error={str(name): round(float(part.median()), 6) for name, part in
                                                                                                      held.groupby("method")["location_error"]})


def choice_checks(qc, state):
    units, records = state["units"], state["records"]
    checks = pd.read_parquet(C.DATA / "unit_checks.parquet")
    methods = pd.read_parquet(C.DATA / "unit_methods.parquet")
    fits = pd.read_parquet(C.DATA / "candidate_fits.parquet").set_index(C.UNIT)
    active = units[units["status"] == "in_unit"]
    qc.check("every_unit_has_one_method_row", len(methods) == len(active) and set(map(tuple, methods[C.UNIT].to_numpy())) == set(map(tuple, active[C.UNIT].to_numpy())))
    policy = C.CONFIG["checks"]
    applicable = checks[checks["applicable"]]
    recomputed = np.where(applicable["check"] == "study_mixing", applicable["after"] <= applicable["before"] + policy["mixing_worsening_limit"],
                          applicable["after"] >= applicable["limit"])
    qc.check("every_pass_follows_from_its_values", (recomputed == applicable["passed"].to_numpy()).all() and checks.loc[~checks["applicable"], "passed"].all())
    qc.check("mixing_drops_only_unanchored_units", (checks["gate"] == ((checks["check"] != "study_mixing") | (checks["tier"] != "anchored"))).all())
    passing = (checks["passed"] | ~checks["gate"]).groupby([checks[name] for name in [*C.UNIT, "method"]]).all()
    mixing = checks[checks["check"] == "study_mixing"].groupby([*C.UNIT, "method"])["applicable"].any()
    right = True
    for row in methods.itertuples(index=False):
        key = (row.group, row.study_accession, row.platform)
        tried = "healthy_anchored" if row.tier == "anchored" else ("matched_groups" if fits.loc[key, "applicable"] else "")
        passed = bool(tried) and bool(passing.get((*key, tried), False))
        tested = bool(tried) and (row.tier == "anchored" or bool(mixing.get((*key, tried), False)))
        kept = passed and tested
        right = right and row.method_tried == tried and row.method == (tried if kept else "") and (row.status == "corrected") == kept
        right = right and row.reason == ("" if kept else ("no_method_passed_the_checks" if not passed else "no_check_applicable") if tried
                                         else fits.loc[key, "reason"])
    qc.check("every_unit_kept_only_when_its_correction_passed", right)
    kept_unanchored = methods[(methods["status"] == "corrected") & (methods["tier"] == "unanchored")]
    qc.check("units_without_healthy_controls_kept_only_where_their_mixing_check_applied",
             all(bool(mixing.get((row.group, row.study_accession, row.platform, row.method), False)) for row in kept_unanchored.itertuples(index=False)))
    # The counts of effects beyond the report limits decide nothing; they are derived again from the effect tables.
    effect_table = pd.concat([pd.read_parquet(C.DATA / "anchored_effects.parquet"), pd.read_parquet(C.DATA / "candidate_effects.parquet")], ignore_index=True)
    effect_table = effect_table.set_index([*C.UNIT, "method"]).sort_index()
    reported = {group: value["beyond_report_limits"] for group, value in C.read_json(C.receipt_path(6))["groups"].items()}
    counted = True
    for row in methods[methods["method_tried"] != ""].itertuples(index=False):
        part = effect_table.loc[(row.group, row.study_accession, row.platform, row.method_tried)]
        found = dict(locations=int((np.abs(part["location"].to_numpy()) > policy["location_report_limit"]).sum()),
                     scales=int((np.abs(np.log(part["scale"].to_numpy())) > np.log(policy["scale_report_limit"])).sum()))
        counted = counted and reported.get(row.group, {}).get(f"{row.study_accession} {row.platform}") == found
    counted = counted and sum(len(value) for value in reported.values()) == int((methods["method_tried"] != "").sum())
    qc.check("effects_beyond_the_report_limits_are_counted_in_the_receipt", counted)
    repeats = pd.read_parquet(C.DATA / "repeat_test.parquet")
    qc.check("repeat_comparisons_are_valid", repeats["distance_before"].ge(0).all() and repeats["distance_after"].ge(0).all()
             and (repeats["rank_after"] >= 1).all() and (repeats["rank_after"] <= repeats["candidates"]).all() and repeats["shared_genes"].gt(0).all())
    state["methods"] = methods
    state["summary"].update(units=len(methods), corrected_units=int((methods["status"] == "corrected").sum()),
                            by_method={str(name): int(count) for name, count in methods.loc[methods["status"] == "corrected", "method"].value_counts().items()},
                            repeat_pairs=len(repeats), repeat_median_distance_after=round(float(repeats["distance_after"].median()), 6),
                            repeat_pairs_nearest_after=int((repeats["rank_after"] == 1).sum()))


def table_checks(qc, state):
    records, genes, cohort, methods = state["records"], state["genes"], state["cohort"], state["methods"]
    ledger = pd.read_parquet(C.DATA / "correction_ledger.parquet")
    parameters = pd.read_parquet(C.DATA / "correction_parameters.parquet")
    qc.check("ledger_has_every_record_once_with_an_outcome", len(ledger) == len(records) and not ledger.duplicated(C.KEYS).any()
             and ledger["outcome"].isin(["in_corrected_table", "corrected_for_checks_only", "left_out"]).all()
             and (ledger.loc[ledger["outcome"] == "left_out", "reason"] != "").all() and set(ledger["reason"]) <= set(C.REASONS) | {""})
    seen, worst, metadata, rebuilt = set(), 0.0, True, True
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        model = C.model_arrays(group, order)
        frame = pd.read_parquet(C.CORRECTED / f"{group}_expression.parquet")
        samples = pd.read_parquet(C.CORRECTED / f"{group}_samples.parquet")
        inside = ledger[(ledger["group"] == group) & (ledger["outcome"] == "in_corrected_table")]
        qc.check(f"{group}_columns_are_the_samples_in_the_table", list(frame.columns[2:]) == list(samples["sample_id"]) and samples["sample_id"].is_unique
                 and set(samples["sample_id"]) == set(inside["sample_id"]) and samples["is_primary_record"].all())
        qc.check(f"{group}_rows_are_group_genes", set(frame["hgnc_id"]) <= set(order) and frame["hgnc_id"].is_unique
                 and list(frame["hgnc_id"]) == [gene for gene in order if gene in set(frame["hgnc_id"])])
        seen_here = set(samples["sample_id"])
        qc.check(f"{group}_samples_in_no_other_table", not (seen & seen_here))
        seen |= seen_here
        source = cohort.set_index(["source", "study_accession", "sample_accession", "platform"]).loc[
            list(zip(samples["source"], samples["study_accession"], samples["sample_accession"], samples["platform"]))]
        for name in cohort.columns:
            if name in ("source", "study_accession", "sample_accession", "platform"):
                continue
            left, right = samples[name].reset_index(drop=True), source[name].reset_index(drop=True)
            metadata = metadata and bool(((left == right) | (left.isna() & right.isna())).all())
        qc.check(f"{group}_sample_types_are_of_the_group", samples["tissue"].isin(C.CONFIG["groups"][group]["tissues"]).all())
        # Every column is rebuilt from its aligned values, its study's effects and the model.
        values = frame.set_index("hgnc_id").drop(columns="symbol").reindex(order)
        for (study, platform), part in samples.groupby(["study_accession", "platform"], sort=False):
            effect = parameters[(parameters["group"] == group) & (parameters["study_accession"] == study) & (parameters["platform"] == platform)]
            qc.check(f"one_method_for_{group}_{C.unit_name(study, platform)}", effect["method"].nunique() == 1 and effect["method"].iloc[0] == part["unit_method"].iloc[0])
            effect = effect.set_index("hgnc_id").reindex(order)
            raw = C.read_aligned(group, study, platform, part["sample_accession"])
            again = C.correct(raw, part["female"].to_numpy(), model, effect["location"].to_numpy(), effect["scale"].to_numpy())
            stored = values[list(part["sample_id"])].to_numpy(dtype="float64")
            difference = (again - stored)[~np.isnan(again)]
            rebuilt = rebuilt and bool((np.isnan(again) == np.isnan(stored)).all()) and bool((np.abs(difference) <= C.CONFIG["quality_control"]["rebuild_tolerance"]).all())
            worst = largest(difference, worst)
        state["summary"][f"{group}_samples"] = len(samples)
        state["summary"][f"{group}_genes"] = len(frame)
    qc.check("sample_metadata_equals_the_cohort_table", metadata)
    qc.check("corrected_values_rebuilt_from_aligned_values_effects_and_model", rebuilt, shown(worst))
    chosen = methods[methods["status"] == "corrected"]
    qc.check("parameters_are_those_of_the_chosen_methods", set(map(tuple, parameters[["group", "study_accession", "platform", "method"]].drop_duplicates().to_numpy()))
             == set(map(tuple, chosen[["group", "study_accession", "platform", "method"]].to_numpy())))
    state["summary"].update(records=len(ledger), by_outcome={str(name): int(count) for name, count in ledger["outcome"].value_counts().items()})


def total_checks(qc, state):
    selection, tables = C.read_json(C.receipt_path(1)), C.read_json(C.receipt_path(7))
    summary = state["summary"]
    qc.check("receipt_totals_agree_with_the_tables", tables["records"]["by_outcome"] == dict(sorted(summary["by_outcome"].items()))
             and all(tables["groups"][group]["samples"] == summary[f"{group}_samples"] for group in C.CONFIG["groups"])
             and selection["records"]["total"] == summary["records"])


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/7 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/7 Records, units and genes", lambda: selection_checks(qc, state))
        qc.phase("3/7 Aligned values against the gene matrices", lambda: alignment_checks(qc, state))
        qc.phase("4/7 Flat genes, anchored and unanchored corrections", lambda: correction_checks(qc, state))
        qc.phase("5/7 Checks and the choice of methods", lambda: choice_checks(qc, state))
        qc.phase("6/7 Corrected tables against the aligned values", lambda: table_checks(qc, state))
        qc.phase("7/7 Stage receipts and totals", lambda: total_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, 8) if C.receipt_path(number).exists()]),
                  eligibility="CORRECTION_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[08 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    {summary['corrected_units']} of {summary['units']} study units corrected {summary['by_method']}; "
              + "; ".join(f"{group}: {summary[group + '_samples']:,} samples, {summary[group + '_genes']:,} genes" for group in C.CONFIG["groups"])
              + f"; {summary['repeat_pairs_nearest_after']} of {summary['repeat_pairs']} deposits nearest their twin after", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
