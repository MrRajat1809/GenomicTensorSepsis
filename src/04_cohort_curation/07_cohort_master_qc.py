"""Verify the catalogue, the links, the curated labels, the quality flags and the cohort tables of the cohort stage from their files."""
import importlib
import re
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_catalogue_sample_metadata")
L = importlib.import_module("03_curate_sample_labels")
B = importlib.import_module("06_define_benchmark_splits")

KEYS = ["study_accession", "sample_accession"]
SEPSIS_WORDING = re.compile(r"seps|septic", re.IGNORECASE)
HEALTHY_WORDING = re.compile(r"health|heathy|normal|volunteer", re.IGNORECASE)
MEANING = ("PASS means that every included sample record with expression values is catalogued, linked, labelled and flagged exactly once, and that every "
           "included record without values is listed once with its reason and its recorded fields and is in no cohort table; that every repeated record points to the "
           "earliest released record of its sample, and every link that rests on values or on a named re-analysis is one that the configuration declares "
           "and that the evidence supports, with every record that a study marks as a re-analysis linked; that every curated value lies in its vocabulary "
           "and names the field, wording or study sentence it rests on, with every sentence found verbatim in the study record and every numeric age grouped "
           "from its recorded number; that every sepsis label and every healthy-control label rests on such a wording, and every phase agrees with its "
           "class; that technical replicates are flagged as the records or the values show and share the patient of their first record; and that the cohort tables, the "
           "patients and the one-sample-per-patient marks follow from these tables, with no patient both a healthy control and of another class, and that "
           "the patient-safe groups join exactly the cohort groups that an identified patient links; and that every record of every benchmark task is "
           "eligible or excluded by the first reason that applies, with its label equal to its curated value and one sample for every patient, that the "
           "splits rebuilt from these tasks are the stored ones, and that no fold holds a patient-safe group, a cohort group, a sample or an identified "
           "patient on both sides, with every scored group holding the minimum of each class. It does "
           "not assess whether a deposited label, or the curation rule that reads it, is clinically right.")


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
        print(f"[07] {name}", flush=True)
        action()


def provenance_checks(qc, state):
    state["start"] = C.begin(7)
    qc.check("all_stage_receipts_and_code_bindings", True)
    included = pd.read_parquet(C.ANNOTATION / "sample_index.parquet", columns=["source", *KEYS, "platform", "family", "in_matrix", "reason"])
    qc.check("included_samples_unique", len(set(zip(included["study_accession"], included["sample_accession"]))) == len(included)
             and included["in_matrix"].dtype == bool)
    # The cohort tables hold the records with expression values; the others belong in the ledger of records without values.
    state["included"] = included
    state["index"] = included[included["in_matrix"]].reset_index(drop=True)
    state["records"] = set(zip(state["index"]["study_accession"], state["index"]["sample_accession"]))
    state["absent"] = set(zip(included.loc[~included["in_matrix"], "study_accession"], included.loc[~included["in_matrix"], "sample_accession"]))


def same_records(frame, state):
    return len(frame) == len(state["index"]) and set(zip(frame["study_accession"], frame["sample_accession"])) == state["records"]


def keys_of(frame, study="study_accession", sample="sample_accession"):
    return set(zip(frame[study], frame[sample]))


def catalogue_checks(qc, state):
    metadata = pd.read_parquet(C.DATA / "sample_metadata.parquet")
    attributes = pd.read_parquet(C.ACQUIRED / "sample_attributes.parquet", columns=["study_accession", "sample_accession", "key", "value"])
    attributes = attributes[[key in state["records"] for key in zip(attributes["study_accession"], attributes["sample_accession"])]]
    described = metadata[metadata["origin"] != "record"]
    qc.check("every_sample_catalogued", set(zip(metadata["study_accession"], metadata["sample_accession"])) == state["records"])
    qc.check("every_recorded_attribute_catalogued_verbatim", len(described) == len(attributes)
             and sorted(zip(described["sample_accession"], described["key"], described["value"])) == sorted(zip(attributes["sample_accession"], attributes["key"],
                                                                                                                 attributes["value"])))
    fields = pd.read_parquet(C.DATA / "metadata_fields.parquet")
    qc.check("field_list_agrees_with_the_catalogue", len(fields) == metadata.groupby(["study_accession", "platform", "origin", "field"]).ngroups)

    # The ledger is rebuilt here from the sample index and the sample records, apart from the catalogue of stage 01.
    included, ledger = state["included"], pd.read_parquet(C.DATA / "records_without_values.parquet")
    absent = included[~included["in_matrix"]].set_index(KEYS)
    listed = ledger.set_index(KEYS)
    qc.check("every_record_without_values_is_listed_once_with_its_reason", list(ledger.columns) == C.LEDGER_COLUMNS and not ledger.duplicated(KEYS).any()
             and set(listed.index) == state["absent"] and (absent.loc[listed.index, "reason"].fillna("").to_numpy() == listed["reason"].to_numpy()).all()
             and (listed["reason"] != "").all() and (absent.loc[listed.index, "family"].to_numpy() == listed["family"].to_numpy()).all()
             and (absent.loc[listed.index, "platform"].to_numpy() == listed["platform"].to_numpy()).all()
             and (absent.loc[listed.index, "source"].to_numpy() == listed["source"].to_numpy()).all(), len(ledger))
    records = pd.read_parquet(C.ACQUIRED / "samples.parquet", columns=[*KEYS, "title", "source_name", "description"]).set_index(KEYS)
    attributes = pd.read_parquet(C.ACQUIRED / "sample_attributes.parquet", columns=[*KEYS, "channel", "position", "kind", "key", "value"])
    attributes = attributes[[key in state["absent"] for key in zip(attributes["study_accession"], attributes["sample_accession"])]]
    written = {}
    for row in attributes.sort_values([*KEYS, "channel", "position", "key"], kind="stable").itertuples(index=False):
        name = str(row.key) if str(row.key).strip() else "(no key)"
        text = f"{name}: {row.value}" if row.channel <= 1 else f"channel {row.channel} {name}: {row.value}"
        written.setdefault((row.study_accession, row.sample_accession, row.kind), []).append(text)
    kept = all((records.loc[listed.index, name].fillna("").astype(str).to_numpy() == listed[name].to_numpy()).all() for name in ("title", "source_name", "description"))
    kept = kept and set(kind for _, _, kind in written) <= set(C.LEDGER_KINDS) and all(
        (listed[column].to_numpy() == [" | ".join(written.get((study, sample, kind), [])) for study, sample in listed.index]).all()
        for kind, column in C.LEDGER_KINDS.items())
    qc.check("ledger_keeps_the_recorded_fields_that_the_curation_reads", kept and len(ledger) == len(listed))
    receipt = C.read_json(C.receipt_path(1))["records"]
    reasons = {str(name): int(count) for name, count in ledger["reason"].value_counts().sort_index().items()}
    qc.check("included_records_reconcile_with_the_ledger_by_reason", receipt == dict(included=len(included), with_values=len(state["index"]),
                                                                                       without_values=dict(total=len(ledger), by_reason=reasons))
             and len(state["index"]) + len(ledger) == len(included), dict(included=len(included), with_values=len(state["index"]), without_values=reasons))
    state["ledger"] = ledger
    state["summary"] = dict(included_records=len(included), samples=len(state["index"]), records_without_values=len(ledger),
                            records_without_values_by_reason=reasons, recorded_values=len(metadata), fields=int(described["field"].nunique()))


def link_checks(qc, state):
    links = pd.read_parquet(C.DATA / "sample_links.parquet")
    released = dict(pd.read_parquet(C.ACQUIRED / "included_studies.parquet", columns=["accession", "release_date"]).itertuples(index=False, name=None))
    qc.check("every_sample_linked_once", same_records(links, state))
    qc.check("one_primary_record_for_every_sample", (links.groupby("sample_id")["is_primary_record"].sum() == 1).all()
             and (links["sample_id"] == links["primary_study_accession"] + ":" + links["primary_sample_accession"]).all()
             and (links.groupby("sample_id").size() == links.groupby("sample_id")["records_of_sample"].first()).all())
    dates = links["study_accession"].map(released)
    # Of the records of one sample released on one day, the primary is the lowest by the letters and then the number of its study and sample accessions.
    ordered = links.assign(date=dates, study=links["study_accession"].map(C.accession_order), sample=links["sample_accession"].map(C.accession_order))
    earliest = ordered.sort_values(["date", "study", "sample"], kind="stable").groupby("sample_id")[["study_accession", "sample_accession"]].first()
    qc.check("primary_record_is_the_earliest_released", (links["primary_study_accession"].map(released) <= dates).all()
             and (links.loc[links["is_primary_record"], "study_accession"] == links.loc[links["is_primary_record"], "primary_study_accession"]).all()
             and (earliest["study_accession"] + ":" + earliest["sample_accession"] == earliest.index).all())
    qc.check("records_of_one_sample_share_a_cohort_group", (links.groupby("sample_id")["cohort_group"].nunique() == 1).all()
             and (links.groupby(["study_accession", "platform"])["cohort_group"].nunique() == 1).all())
    qc.check("repeated_records_are_linked_by_evidence_and_agree", (links.loc[links["records_of_sample"] > 1, "linked_by"] != "").all()
             and (links.loc[links["records_of_sample"] == 1, "linked_by"] == "").all() and links["values_agree"].all())
    identity, candidates = C.CONFIG["links"]["identity"], pd.read_parquet(C.DATA / "identity_candidates.parquet")
    sample_of = dict(zip(zip(links["study_accession"], links["sample_accession"]), links["sample_id"]))
    decided = candidates[candidates["decision"] == "same_sample"]
    qc.check("every_identity_link_is_a_declared_mutual_best_match", (candidates["correlation"] >= identity["candidate_correlation"]).all()
             and (candidates["margin"] >= identity["candidate_margin"]).all()
             and all(sample_of[(row.study_a, row.sample_a)] == sample_of[(row.study_b, row.sample_b)] for row in decided.itertuples(index=False))
             and all(any({f"{row.study_a}|{row.platform_a}", f"{row.study_b}|{row.platform_b}"} == {entry["a"], entry["b"]}
                         and row.correlation >= entry["minimum_correlation"] for entry in C.REGISTER["same_sample"]) for row in decided.itertuples(index=False)),
             len(decided))
    qc.check("identity_screen_compares_distinct_samples", not candidates["linked_by_files"].any() and not candidates["linked_by_files"].isna().any())
    replicates = pd.read_parquet(C.DATA / "replicate_libraries.parquet")
    paired = set(zip(candidates.loc[candidates["decision"] == "replicate_libraries", "sample_a"], candidates.loc[candidates["decision"] == "replicate_libraries", "sample_b"]))
    qc.check("every_replicate_pair_is_a_candidate_of_a_declared_study", all((row.first, row.second) in paired or (row.second, row.first) in paired
                                                                            for row in replicates.itertuples(index=False))
             and set(replicates["study_accession"] + "|" + replicates["platform"]) <= {entry["study"] for entry in C.REGISTER["replicates"]}, len(replicates))
    named = links["linked_by"].str.contains("record_names_a_re_analysis")
    qc.check("re_analysis_links_join_their_declared_studies", set(links.loc[named, "study_accession"]) <= {side["study"] for rule in C.REGISTER["record_links"]
                                                                                                          for side in (rule["a"], rule["b"])}, int(named.sum()))
    # Every record that its study marks as a re-analysis is one sample with the record of the other study that carries its key.
    recorded = pd.read_parquet(C.ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", "title", "description", "source_name"])
    joined, marked = True, 0
    for rule in C.REGISTER["record_links"]:
        a, b =(recorded[recorded["study_accession"] == rule[side]["study"]] for side in ("a", "b"))
        keys = {}
        for sample, value in zip(a["sample_accession"], a[rule["a"]["field"]].fillna("").astype(str)):
            match = re.search(rule["a"]["pattern"], value)
            if match:
                keys.setdefault(match.group(1), []).append(sample)
        before = marked
        for sample, value in zip(b["sample_accession"], b[rule["b"]["field"]].fillna("").astype(str)):
            if re.search(rule["b"]["required"], value):
                marked += 1
                match = re.search(rule["b"]["pattern"], value)
                partners = keys.get(match.group(1), []) if match else []
                joined = (joined and len(partners) == 1 and (rule["b"]["study"], sample) in sample_of
                          and sample_of[(rule["b"]["study"], sample)] == sample_of.get((rule["a"]["study"], partners[0])))
        joined = joined and marked > before
    qc.check("every_record_marked_as_a_re_analysis_joins_its_counterpart", joined, marked)
    state["links"] = links
    state["summary"].update(distinct_samples=int(links["sample_id"].nunique()), repeated_records=int((~links["is_primary_record"]).sum()),
                            cohort_groups=int(links["cohort_group"].nunique()), study_platforms=links.groupby(["study_accession", "platform"]).ngroups)


def label_checks(qc, state):
    policy = C.CONFIG["curation"]
    labels = pd.read_parquet(C.DATA / "sample_labels.parquet")
    sources = pd.read_parquet(C.DATA / "label_sources.parquet")
    vocabularies = policy["vocabularies"]
    qc.check("every_sample_labelled_once", same_records(labels, state))
    qc.check("rule_sets_name_exactly_the_study_platforms_with_values", set(C.REGISTER["rules"]) == set(state["index"]["study_accession"] + "|"
                                                                                                       + state["index"]["platform"])
             and set(C.REGISTER["sepsis_definitions"]) <= set(state["index"]["study_accession"]), len(C.REGISTER["rules"]))
    qc.check("curated_values_within_their_vocabularies", all(set(labels[name]) <= set(vocabularies[name]) - {"infection_sepsis_not_stated"}
                                                              for name in L.CATEGORICAL) and set(labels["age_group"]) <= set(vocabularies["age_group"])
             and set(labels["time_unit"]) <= set(vocabularies["time_unit"]) | {""} and set(labels["sepsis_phase"]) <= set(vocabularies["sepsis_phase"])
             and set(labels["sepsis_definition"]) <= set(vocabularies["sepsis_definition"]))
    septic = labels["sample_class"] == "sepsis"
    qc.check("every_phase_agrees_with_its_class", labels.loc[septic, "sepsis_phase"].isin(L.SEPSIS_PHASES).all()
             and labels.loc[~septic, "sepsis_phase"].isin(L.OTHER_PHASES).all()
             and (labels.loc[labels["sepsis_phase"] == "recovery_or_follow_up", "sample_class"] == "other_condition").all(),
             {str(name): int(count) for name, count in labels["sepsis_phase"].value_counts().items()})
    qc.check("every_variable_of_every_sample_has_its_source", len(sources) == len(labels) * (len(L.VARIABLES) - 1)
             and not sources.duplicated([*KEYS, "variable"]).any())
    decided = sources[sources["rule"] > 0]
    # Without a rule a value is not recorded; only the phase follows from the class, and a study without a stated definition says nothing of it.
    undecided = sources[sources["rule"] == 0]
    defaults = ((undecided["variable"] == "phase") & (undecided["basis"] == "from_sample_class")) | (
        (undecided["variable"] == "definition") & (undecided["value"] == "not_stated_in_record"))
    qc.check("every_decided_value_names_its_field_or_its_study_sentence", ((decided["fields"] != "") | (decided["study_statement"] != "")).all()
             and (undecided.loc[~defaults, "value"] == L.NOT_RECORDED).all())
    studies = pd.read_parquet(C.ACQUIRED / "studies.parquet", columns=["accession", "title", "summary", "design"]).set_index("accession")
    texts = {study: L.plain(f"{row.title}. {row.summary} {row.design}") for study, row in studies.iterrows() if study in set(labels["study_accession"])}
    quoted = sources[sources["study_statement"] != ""].drop_duplicates(["study_accession", "study_statement"])
    qc.check("every_study_sentence_is_in_its_study_record", all(L.plain(sentence) in texts[study] for study, sentence in
                                                                zip(quoted["study_accession"], quoted["study_statement"])), len(quoted))
    classes = sources[sources["variable"] == "sample_class"]
    wording = classes["fields"] + " " + classes["recorded"] + " " + classes["study_statement"]
    sepsis = classes["value"] == "sepsis"
    healthy = classes["value"] == "healthy_control"
    qc.check("every_sepsis_label_rests_on_a_sepsis_wording", wording[sepsis].map(lambda text: bool(SEPSIS_WORDING.search(text))).all(), int(sepsis.sum()))
    qc.check("every_healthy_control_rests_on_a_healthy_wording", wording[healthy].map(lambda text: bool(HEALTHY_WORDING.search(text))).all(), int(healthy.sum()))
    qc.check("infection_flag_only_on_other_conditions", (labels.loc[labels["infection_sepsis_not_stated"], "sample_class"] == "other_condition").all())
    # The age is read again from the recorded number and the unit of its rule, and grouped in days before any rounding.
    aged = labels[labels["age_years"].notna()].set_index(KEYS)
    limits = policy["age_groups"]
    numbers = sources[(sources["variable"] == "age") & (sources["rule"] > 0)]
    rules = [L.rule_list(C.REGISTER["rules"][f"{study}|{platform}"]["age"])[number - 1] for study, platform, number in
             zip(numbers["study_accession"], numbers["platform"], numbers["rule"])]
    numbers = numbers[["number" in rule for rule in rules]].set_index(KEYS)
    days = numbers["value"].astype(float) * [L.AGE_DAYS[rule["unit"]] for rule in rules if "number" in rule]
    expected = pd.Series(np.where(days <= limits["neonate_up_to_days"], "neonate", np.where(days >= limits["adult_from_years"] * 365.25, "adult", "child")),
                         index=days.index)
    qc.check("age_groups_follow_the_recorded_ages", aged["age_years"].between(0, 120).all() and set(days.index) == set(aged.index) and len(days) == len(aged)
             and (np.abs(aged.loc[days.index, "age_years"] - (days / 365.25).round(4)) <= 1e-9).all()
             and (aged.loc[days.index, "age_group"] == expected).all(), len(aged))
    qc.check("time_points_carry_their_wording_and_unit", (labels["time_order"].notna() == (labels["time_unit"] != "")).all()
             and (labels.loc[labels["time_order"].notna(), "time_point"] != "").all())
    qc.check("recorded_outcomes_carry_their_definition", ((labels["outcome"] != L.NOT_RECORDED) == (labels["outcome_definition"] != "")).all())
    scores = pd.read_parquet(C.DATA / "sample_scores.parquet")
    qc.check("scores_belong_to_included_samples", set(zip(scores["study_accession"], scores["sample_accession"])) <= state["records"]
             and not scores.duplicated([*KEYS, "score"]).any() and np.isfinite(scores["value"]).all())
    state["labels"] = labels
    state["summary"].update(rule_sets=len(C.REGISTER["rules"]), study_sentences=len(quoted), samples_with_a_score=int(scores[KEYS].drop_duplicates().shape[0]))


def quality_checks(qc, state):
    policy = C.CONFIG["quality"]
    quality = pd.read_parquet(C.DATA / "sample_quality.parquet")
    qc.check("every_sample_has_one_quality_row", same_records(quality, state))
    counts =pd.read_parquet(C.EXPRESSION / "count_table_samples.parquet", columns=[*KEYS, "library_size"]).set_index(KEYS)["library_size"]
    measured = quality[quality["library_size"].notna()]
    qc.check("library_flags_follow_the_counts", len(measured) == len(counts)
             and np.allclose(measured["library_size"].to_numpy(), counts.loc[list(zip(measured["study_accession"], measured["sample_accession"]))].to_numpy())
             and (measured["low_library_size"] == (measured["library_size"] < policy["minimum_library_size"])).all()
             and not quality.loc[quality["library_size"].isna(), "low_library_size"].any())
    ranged = quality[quality["dynamic_range"].notna()]
    arrays = pd.read_parquet(C.EXPRESSION / "intensity_table_arrays.parquet", columns=[*KEYS, "signal_median", "signal_percentile_99"]).set_index(KEYS)
    arrays = arrays.loc[[key for key in arrays.index if key in state["records"]]]
    again = (arrays["signal_percentile_99"] / arrays["signal_median"]).where(arrays["signal_median"] > 0)
    stored = quality.set_index(KEYS)["dynamic_range"]
    qc.check("array_flags_follow_their_measures", (ranged["low_dynamic_range"] == (ranged["dynamic_range_to_study_median"]
                                                                                  < policy["dynamic_range_below_share_of_study_median"])).all()
             and np.isclose(ranged.groupby(["study_accession", "platform"])["dynamic_range_to_study_median"].median(), 1).all()
             and set(ranged.set_index(KEYS).index) == set(again.dropna().index) and not quality.loc[quality["dynamic_range"].isna(), "low_dynamic_range"].any()
             and np.allclose(stored.loc[again.dropna().index].to_numpy(), again.dropna().to_numpy(), rtol=1e-12, atol=0),
             dict(measured=len(ranged), without_a_positive_median=int((~(arrays["signal_median"] > 0)).sum())))
    names = quality[list(importlib.import_module("04_flag_sample_quality").FLAGS)]
    qc.check("flag_summary_agrees_with_the_flags", (quality["any_flag"] == names.any(axis=1)).all() and ((quality["record_flag"] != "") == quality["flagged_in_record"]).all())
    flagged = quality[quality["technical_replicate"]]
    first = set(zip(flagged["study_accession"], flagged["replicate_of"]))
    qc.check("technical_replicates_point_to_an_unflagged_first_record_of_their_study", (flagged["replicate_of"] != "").all()
             and (quality.loc[~quality["technical_replicate"], "replicate_of"] == "").all()
             and first <= set(zip(quality.loc[~quality["technical_replicate"], "study_accession"], quality.loc[~quality["technical_replicate"], "sample_accession"]))
             and all(C.accession_order(lead) < C.accession_order(sample) for lead, sample in zip(flagged["replicate_of"], flagged["sample_accession"])),
             len(flagged))
    state["quality"] = quality
    state["summary"].update(flagged_samples=int(quality["any_flag"].sum()))


def table_checks(qc, state):
    table = pd.read_parquet(C.DATA / "cohort_samples.parquet")
    cohorts = pd.read_parquet(C.DATA / "cohort_summary.parquet")
    groups = pd.read_parquet(C.DATA / "cohort_group_summary.parquet")
    qc.check("sample_table_holds_every_record_once", same_records(table, state))
    order = list(zip(table["study_accession"], table["sample_accession"]))
    labels, links, quality = (state[name].set_index(KEYS).loc[order] for name in ("labels", "links", "quality"))
    qc.check("sample_table_agrees_with_links_labels_and_flags", all((table[name].to_numpy() == labels[name].to_numpy()).all() for name in
                                                                    ("sample_class", "sepsis_phase", "sepsis_definition", "tissue", "sex", "age_group", "patient_id",
                                                                     "outcome", "septic_shock"))
             and all((table[name].to_numpy() == links[name].to_numpy()).all() for name in ("sample_id", "is_primary_record", "cohort_group"))
             and (table["any_flag"].to_numpy() == quality["any_flag"].to_numpy()).all())
    known = table[table["patient_identified"]]
    qc.check("patients_are_consistent", (known.groupby("sample_id")["patient_key"].nunique() == 1).all()
             and (table["patient_identified"] == (table["patient_key"] != "")).all() and (table.loc[table["patient_id"] != "", "patient_identified"]).all())
    # The patient-safe groups are rebuilt here by a search over the cohort groups that an identified patient links, apart from the joining of stage 05.
    linked = {}
    for names in known.groupby("patient_key")["cohort_group"].agg(set):
        for name in names:
            linked.setdefault(name, set()).update(names)
    expected = {}
    for first in sorted(table["cohort_group"].unique()):
        if first in expected:
            continue
        seen, frontier = {first}, [first]
        while frontier:
            for name in linked.get(frontier.pop(), set()) - seen:
                seen.add(name)
                frontier.append(name)
        expected.update({name: "+".join(sorted(seen)) for name in seen})
    qc.check("patient_safe_groups_join_exactly_the_cohort_groups_a_patient_links", (table["patient_safe_group"] == table["cohort_group"].map(expected)).all()
             and (known.groupby("patient_key")["patient_safe_group"].nunique() == 1).all() and (table.groupby("sample_id")["patient_safe_group"].nunique() == 1).all(),
             sorted(name for name in set(expected.values()) if "+" in name))
    replicate = state["quality"][state["quality"]["technical_replicate"]]
    patient = dict(zip(zip(table["study_accession"], table["sample_accession"]), table["patient_key"]))
    qc.check("technical_replicates_share_the_patient_of_their_first_record", all(patient[(study, sample)] != "" and patient[(study, sample)] == patient[(study, lead)]
                                                                                 for study, sample, lead in zip(replicate["study_accession"],
                                                                                                                replicate["sample_accession"], replicate["replicate_of"])),
             len(replicate))
    marked = table[table["one_per_patient"]]
    qc.check("one_sample_marked_for_every_patient", marked["is_primary_record"].all() and not marked.loc[marked["patient_identified"], "patient_key"].duplicated().any()
             and set(marked.loc[marked["patient_identified"], "patient_key"]) == set(known["patient_key"])
             and int((~marked["patient_identified"]).sum()) == int((table["is_primary_record"] & ~table["patient_identified"]).sum()))
    clean = known[known["is_primary_record"] & ~known["any_flag"]].groupby("patient_key").size()
    chosen = marked[marked["patient_identified"]].set_index("patient_key")
    qc.check("a_patient_with_an_unflagged_sample_is_marked_by_one", not chosen.loc[chosen.index.isin(clean.index), "any_flag"].any())
    sample_classes = known.groupby(["patient_key", "sample_id"])["sample_class"].agg(lambda values: frozenset(values) - {"not_recorded"})
    mixed = sample_classes.groupby(level="patient_key").agg(
        lambda sets: any("healthy_control" in first and bool(second - first or first - {"healthy_control", "non_septic_control"})
                         for first in sets for second in sets)).astype(bool)
    qc.check("no_patient_is_both_a_healthy_control_and_of_another_class", not mixed.any(), sorted(mixed[mixed].index)[:20])
    differing = table.groupby("sample_id")["sample_class"].nunique() > 1
    qc.check("class_disagreements_are_marked", (table["class_differs_between_records"] == table["sample_id"].map(differing)).all())
    distinct = table[table["is_primary_record"]]
    qc.check("cohort_summaries_agree_with_the_sample_table", len(cohorts) == state["summary"]["study_platforms"] and int(cohorts["records"].sum()) == len(table)
             and int(cohorts["sepsis"].sum()) == int((table["sample_class"] == "sepsis").sum()) and int(cohorts["primary_records"].sum()) == len(distinct)
             and len(groups) == state["summary"]["cohort_groups"] and int(groups["distinct_samples"].sum()) == len(distinct)
             and int(groups["sepsis"].sum()) == int((distinct["sample_class"] == "sepsis").sum()) and int(groups["one_per_patient"].sum()) == len(marked)
             and (groups["patient_safe_group"] == groups["cohort_group"].map(expected)).all()
             and (cohorts["patient_safe_group"] == cohorts["cohort_group"].map(expected)).all())
    receipt = C.read_json(C.receipt_path(5))
    qc.check("receipt_totals_agree_with_the_tables", receipt["records"]["total"] == len(table) and receipt["distinct_samples"]["total"] == len(distinct)
             and receipt["patients"]["one_per_patient"]["total"] == len(marked) and len(distinct) == state["summary"]["distinct_samples"]
             and receipt["cohorts"]["patient_safe_groups"] == len(set(expected.values())))
    # No table of the stage holds a record without values, or a study platform whose records all lack them.
    candidates, replicates = pd.read_parquet(C.DATA / "identity_candidates.parquet"), pd.read_parquet(C.DATA / "replicate_libraries.parquet")
    held = [keys_of(pd.read_parquet(C.DATA / name, columns=KEYS)) for name in ("sample_metadata.parquet", "label_sources.parquet", "sample_scores.parquet")]
    held += [keys_of(frame) for frame in (state["links"], state["quality"], table)]
    held += [keys_of(candidates, "study_a", "sample_a"), keys_of(candidates, "study_b", "sample_b"), keys_of(replicates, sample="first"),
             keys_of(replicates, sample="second")]
    platforms = set(state["index"]["study_accession"] + "|" + state["index"]["platform"])
    named = [set(frame["study_accession"] + "|" + frame["platform"]) for frame in (pd.read_parquet(C.DATA / "metadata_fields.parquet"),
                                                                                   pd.read_parquet(C.DATA / "cohort_groups.parquet"), cohorts)]
    named.append({name for names in groups["study_platforms"] for name in names.split("; ")})
    qc.check("no_cohort_table_holds_a_record_without_values", all(found <= state["records"] for found in held) and all(found == platforms for found in named),
             dict(records_without_values=len(state["absent"]), study_platforms=len(platforms)))
    state["table"] = table
    state["summary"].update(class_of_distinct_samples={name: int((distinct["sample_class"] == name).sum()) for name in
                                                       C.CONFIG["curation"]["vocabularies"]["sample_class"] if name != "infection_sepsis_not_stated"},
                            patients_identified=int(known["patient_key"].nunique()),
                            distinct_samples_without_a_patient_identifier=int((~distinct["patient_identified"]).sum()), one_per_patient=len(marked),
                            patient_safe_groups=len(set(expected.values())),
                            samples_whose_records_differ_in_class=int(distinct["class_differs_between_records"].sum()))


def benchmark_checks(qc, state):
    policy, table = C.CONFIG["benchmark"], state["table"]
    minimum, split_date = policy["minimum_per_class"], policy["temporal_split_date"]
    names = [task["task"] for task in policy["tasks"]]
    tasks = pd.read_parquet(C.DATA / "benchmark_tasks.parquet")
    splits = pd.read_parquet(C.DATA / "benchmark_splits.parquet")
    qc.check("benchmark_lists_every_record_once_for_every_task", list(tasks.columns) == B.TASK_COLUMNS and len(tasks) == len(table) * len(names)
             and list(dict.fromkeys(tasks["task"])) == names and not tasks.duplicated(["task", *KEYS]).any()
             and all(same_records(tasks[tasks["task"] == name], state) for name in names), len(names))
    # Eligibility is replayed with one selection over all reasons at once, apart from the stepwise marking of stage 06.
    order = list(zip(table["study_accession"], table["sample_accession"]))
    chosen, reasons_right, labels_right = {}, True, True
    for task in policy["tasks"]:
        part = tasks[tasks["task"] == task["task"]].set_index(KEYS).loc[order]
        value = table[task["variable"]].to_numpy()
        labelled = np.isin(value, [task["positive"], task["negative"]])
        in_class = (table["sample_class"] == task["within_class"]).to_numpy() if "within_class" in task else labelled
        expected = np.select([~table["is_primary_record"].to_numpy(), ~table["one_per_patient"].to_numpy(),
                              table["any_flag"].to_numpy(), ~np.isin(table["tissue"].to_numpy(), policy["tissues"]),
                              np.isin(table["sepsis_phase"].to_numpy(), policy["phases_left_out"]), ~in_class, ~labelled], list(B.REASONS), default="")
        keep = expected == ""
        reasons_right = (reasons_right and (part["reason"].to_numpy() == expected).all() and (part["eligible"].to_numpy() == keep).all()
                         and (part["patient_safe_group"].to_numpy() == table["patient_safe_group"].to_numpy()).all())
        labels_right = (labels_right and part["label"][~keep].isna().all() and (part["label_value"][~keep] == "").all() and part["label"][keep].notna().all()
                        and (part["label"][keep].astype(int).to_numpy() == (value[keep] == task["positive"]).astype(int)).all()
                        and (part["label_value"][keep].to_numpy() == value[keep]).all())
        chosen[task["task"]] = table[keep]
    qc.check("benchmark_eligibility_follows_the_sample_table_by_the_first_reason", reasons_right, {name: len(frame) for name, frame in chosen.items()})
    qc.check("benchmark_labels_equal_the_curated_values", labels_right)
    qc.check("one_benchmark_sample_for_every_patient_of_a_task", all(not frame["sample_id"].duplicated().any()
                                                                     and not frame.loc[frame["patient_identified"], "patient_key"].duplicated().any()
                                                                     for frame in chosen.values()))
    # The technology and the first release of every group are read again from the study platforms of the cohort summary.
    family = {(study, platform): name for study, platform, name in
              pd.read_parquet(C.ANNOTATION / "matrix_index.parquet", columns=["accession", "platform", "family"]).itertuples(index=False)}
    technologies, first = {}, {}
    for study, platform, group, date in pd.read_parquet(C.DATA / "cohort_summary.parquet", columns=["study_accession", "platform", "patient_safe_group",
                                                                                                   "release_date"]).itertuples(index=False):
        technologies.setdefault(group, set()).add(policy["technologies"][family[(study, platform)]])
        first[group] = min(first.get(group, str(date)), str(date))
    qc.check("no_patient_safe_group_holds_both_technologies", all(len(found) == 1 for found in technologies.values()), len(technologies))
    technology = {group: min(found) for group, found in technologies.items()}
    rebuilt = []
    for task in policy["tasks"]:
        counts = {}
        for group, value in zip(chosen[task["task"]]["patient_safe_group"], chosen[task["task"]][task["variable"]]):
            counts.setdefault(group, [0, 0])[0 if value == task["positive"] else 1] += 1
        groups = sorted(counts)
        scored = {group for group in groups if min(counts[group]) >= minimum}
        held = {group: "test" if group in scored else "held_out_not_scored" for group in groups}
        folds = [("leave_one_group_out", name, {group: "test" if group == name else "train" for group in groups}) for name in sorted(scored)]
        folds += [("technology_transfer", f"{one}_to_{other}", {group: held[group] if technology[group] == other else "train" for group in groups})
                  for one, other in B.TRANSFERS]
        folds.append(("temporal", B.TEMPORAL_FOLD, {group: held[group] if first[group] >= split_date else "train" for group in groups}))
        rebuilt += [(task["task"], scheme, fold, group, roles[group], technology[group], first[group], *counts[group])
                    for scheme, fold, roles in folds for group in groups]
    qc.check("benchmark_splits_rebuilt_from_the_tasks_are_the_stored_ones", list(splits.columns) == B.SPLIT_COLUMNS
             and list(splits.itertuples(index=False, name=None)) == rebuilt, len(rebuilt))
    # Every record of the cohort, not only the benchmark samples, is placed by its group, so a patient's other samples count as well.
    identified, leaks = table["patient_identified"].to_numpy(), []
    for (name, scheme, fold), part in splits.groupby(["task", "scheme", "fold"], sort=False):
        train, held_out = set(part.loc[part["role"] == "train", "patient_safe_group"]), set(part.loc[part["role"] != "train", "patient_safe_group"])
        left, right = table["patient_safe_group"].isin(train).to_numpy(), table["patient_safe_group"].isin(held_out).to_numpy()
        if (part["patient_safe_group"].duplicated().any() or train & held_out or set(table.loc[left, "cohort_group"]) & set(table.loc[right, "cohort_group"])
                or set(table.loc[left, "sample_id"]) & set(table.loc[right, "sample_id"])
                or set(table.loc[left & identified, "patient_key"]) & set(table.loc[right & identified, "patient_key"])):
            leaks.append(f"{name} {scheme} {fold}")
    qc.check("no_fold_holds_a_group_sample_or_patient_on_both_sides", not leaks, dict(folds=splits.groupby(["task", "scheme", "fold"]).ngroups, leaks=leaks[:20]))
    test, unscored = splits[splits["role"] == "test"], splits[splits["role"] == "held_out_not_scored"]
    logo = splits[(splits["scheme"] == "leave_one_group_out") & (splits["role"] != "train")]
    qc.check("every_scored_group_holds_the_minimum_of_each_class", (test[["positives", "negatives"]].min(axis=1) >= minimum).all()
             and (unscored[["positives", "negatives"]].min(axis=1) < minimum).all() and (logo["role"] == "test").all()
             and (logo["fold"] == logo["patient_safe_group"]).all() and not logo.duplicated(["task", "fold"]).any(), len(test))
    transfer, temporal = splits[splits["scheme"] == "technology_transfer"], splits[splits["scheme"] == "temporal"]
    qc.check("transfer_and_temporal_folds_follow_their_rules",
             all((part.loc[part["role"] == "train", "technology"] == fold.split("_to_")[0]).all()
                 and (part.loc[part["role"] != "train", "technology"] == fold.split("_to_")[1]).all() for (_, fold), part in transfer.groupby(["task", "fold"]))
             and (temporal.loc[temporal["role"] == "train", "first_release"] < split_date).all()
             and (temporal.loc[temporal["role"] != "train", "first_release"] >= split_date).all())
    receipt = C.read_json(C.receipt_path(6))
    qc.check("benchmark_receipt_totals_agree_with_the_tables", receipt["split_rows"] == len(splits)
             and all(receipt["tasks"][name]["records"] == len(table) and receipt["tasks"][name]["eligible"] == len(chosen[name])
                     and receipt["tasks"][name]["excluded"] == {code: int((tasks.loc[tasks["task"] == name, "reason"] == code).sum()) for code in B.REASONS}
                     for name in names))
    state["summary"].update(benchmark={name: dict(samples=len(chosen[name]), scored_groups=int((logo["task"] == name).sum())) for name in names})


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/7 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/7 Catalogue and records without values against the sample records", lambda: catalogue_checks(qc, state))
        qc.phase("3/7 Links between repeated records", lambda: link_checks(qc, state))
        qc.phase("4/7 Curated labels and their sources", lambda: label_checks(qc, state))
        qc.phase("5/7 Quality flags against the processing statistics", lambda: quality_checks(qc, state))
        qc.phase("6/7 Sample table, patients and cohort summaries", lambda: table_checks(qc, state))
        qc.phase("7/7 Benchmark tasks and splits", lambda: benchmark_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, max(C.STAGES)) if C.receipt_path(number).exists()]),
                  eligibility="COHORT_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[07 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    {summary['included_records']:,} included records: {summary['records_without_values']} without expression values, listed with their "
              f"reason; {summary['samples']:,} with values are {summary['distinct_samples']:,} distinct samples in {summary['cohort_groups']} cohort groups "
              f"({summary['patient_safe_groups']} patient-safe groups); "
              f"{summary['patients_identified']:,} patients identified; one sample per patient: {summary['one_per_patient']:,}", flush=True)
        print("    Benchmark: " + "; ".join(f"{name} {entry['samples']:,} samples, {entry['scored_groups']} scored groups"
                                          for name, entry in summary["benchmark"].items()), flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
