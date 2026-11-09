"""Build the sample table and the cohort summaries from the links, the curated labels and the quality flags, and mark one sample for every patient."""
import collections
import importlib

C = importlib.import_module("01_catalogue_sample_metadata")

KEYS = ["source", "study_accession", "sample_accession", "platform"]
CLASSES = ("sepsis", "healthy_control", "non_septic_control", "other_condition", "not_recorded")
SAMPLE_COLUMNS = [*KEYS, "cohort_group", "patient_safe_group", "sample_id", "is_primary_record", "records_of_sample", "linked_by", "sample_class", "infection_sepsis_not_stated",
                  "authors_label", "class_basis", "classes_of_sample_records", "class_differs_between_records", "sepsis_phase", "sepsis_definition", "tissue", "sex",
                  "age_years", "age_group", "patient_id", "patient_key", "patient_identified", "samples_of_patient", "one_per_patient", "time_point", "time_order",
                  "time_unit", "outcome", "outcome_definition", "septic_shock", "quality_flags", "any_flag"]


def class_counts(part, suffix=""):
    return {f"{name}{suffix}": int((part["sample_class"] == name).sum()) for name in CLASSES}


def shares(values):
    counts = values.value_counts()
    return "; ".join(f"{name}={count}" for name, count in counts.items())


def main():
    import pandas as pd

    print(f"[05] Building cohort tables v{C.VERSION}", flush=True)
    start = C.begin(5)
    links = pd.read_parquet(C.DATA / "sample_links.parquet")
    labels = pd.read_parquet(C.DATA / "sample_labels.parquet")
    quality = pd.read_parquet(C.DATA / "sample_quality.parquet", columns=[*KEYS, "quality_flags", "any_flag"])
    table = links.merge(labels, on=KEYS, validate="one_to_one").merge(quality, on=KEYS, validate="one_to_one")
    C.require(len(table) == len(links) == len(labels) == len(quality), "The links, labels and quality flags must describe the same samples")

    # A patient identifier holds within its study, or within the studies of one recorded namespace. Two identifiers are the same patient where one sample
    # carries both, in two deposits; the two libraries of a replicate pair are one patient, and so are the replicate arrays that a study records for one sample.
    namespace = {study: entry["name"] for entry in C.REGISTER["patient_namespaces"] for study in entry["studies"]}
    C.require(set(namespace) <= set(table["study_accession"]), "A patient namespace names a study that has no included record with expression values")
    replicates = pd.read_parquet(C.DATA / "replicate_libraries.parquet")
    table["patient_record"] = (table["study_accession"].map(lambda study: namespace.get(study, study)) + ":" + table["patient_id"]).where(table["patient_id"] != "", "")
    for row in replicates.itertuples(index=False):
        place = table.index[(table["study_accession"] == row.study_accession) & table["sample_accession"].isin([row.first, row.second])]
        C.require(len(place) == 2 and (table.loc[place, "patient_id"] == "").all(), f"Replicate libraries {row.first} and {row.second} must be two records without a patient")
        table.loc[place, "patient_record"] = f"{row.study_accession}:replicate pair {row.first}"
    flags = {rule["flag"] for rule in C.REGISTER["replicate_groups"]}
    arrays = pd.read_parquet(C.DATA / "sample_quality.parquet", columns=["study_accession", "sample_accession", "replicate_of", "replicate_basis"])
    for (study, first), part in arrays[arrays["replicate_basis"].isin(flags)].groupby(["study_accession", "replicate_of"], sort=True):
        place = table.index[(table["study_accession"] == study) & table["sample_accession"].isin([first, *part["sample_accession"]])]
        named = set(table.loc[place, "patient_record"]) - {""}
        C.require(len(place) == len(part) + 1 and len(named) <= 1, f"The replicate arrays of {first} in {study} must be records of at most one patient")
        table.loc[place, "patient_record"] = named.pop() if named else f"{study}:replicate group {first}"
    parent = {name: name for name in set(table["patient_record"]) - {""}}

    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for _, names in table[table["patient_record"] != ""].groupby("sample_id")["patient_record"]:
        names = sorted(set(names))
        for other in names[1:]:
            parent[find(other)] = find(names[0])
    members = collections.defaultdict(list)
    for name in parent:
        members[find(name)].append(name)
    key_of = {name: min(names) for names in members.values() for name in names}
    # Every record of a sample receives the patient that any of its records names.
    of_sample = table[table["patient_record"] != ""].groupby("sample_id")["patient_record"].first().map(key_of)
    table["patient_key"] = table["sample_id"].map(of_sample).fillna("")
    table["patient_identified"] = table["patient_key"] != ""
    C.require((table[table["patient_identified"]].groupby("sample_id")["patient_key"].nunique() == 1).all(), "A sample must belong to one patient")

    # A patient may not be a healthy control in one sample and of another class in another. The records of one sample may differ, as deposits word their
    # controls differently, but a record may call a sample a healthy control only where the others call it a healthy or a non-septic control. A conflict
    # of recorded sex is counted and reported.
    known = table[table["patient_identified"]]
    sample_classes = known.groupby(["patient_key", "sample_id"])["sample_class"].agg(lambda values: frozenset(values) - {"not_recorded"})
    mixed = sample_classes.groupby(level="patient_key").agg(
        lambda sets: any("healthy_control" in first and bool(second - first or first - {"healthy_control", "non_septic_control"})
                         for first in sets for second in sets)).astype(bool)
    C.require(not mixed.any(), f"Patients with a healthy-control sample beside a sample of another class: {sorted(mixed[mixed].index)[:20]}")
    sexes = known[known["sex"].isin(["male", "female"])].groupby("patient_key")["sex"].nunique()
    conflicting_sex = sorted(sexes[sexes > 1].index)
    shared = known.groupby("patient_key")["study_accession"].nunique()

    # A cohort group holds every deposit of a sample, but a patient can give samples to two groups: one study on two platforms, or studies that number
    # their patients in one series. Groups that share an identified patient are joined into one patient-safe group, named by its cohort groups, so that a
    # split by these groups keeps every deposit of a sample and every sample of an identified patient on one side.
    safe = {name: {name} for name in table["cohort_group"].unique()}
    for names in known.groupby("patient_key")["cohort_group"].agg(lambda values: sorted(set(values))):
        joined = set().union(*(safe[name] for name in names))
        for name in joined:
            safe[name] = joined
    table["patient_safe_group"] = table["cohort_group"].map({name: "+".join(sorted(members)) for name, members in safe.items()})
    spanning = known.groupby("patient_key")["cohort_group"].nunique()
    C.require((table[table["patient_identified"]].groupby("patient_key")["patient_safe_group"].nunique() == 1).all(),
              "Every sample of a patient must lie in one patient-safe group")

    classes = table.groupby("sample_id")["sample_class"].agg(lambda values: "; ".join(sorted(set(values))))
    table["classes_of_sample_records"] = table["sample_id"].map(classes)
    table["class_differs_between_records"] = table["classes_of_sample_records"].str.contains("; ")

    primary = table[table["is_primary_record"]].copy()
    primary["no_time"] = primary["time_order"].isna()
    primary["other_tissue"] = primary["tissue"] != "whole_blood"
    primary = primary.sort_values(["patient_key", "any_flag", "no_time", "time_order", "other_tissue", "study_accession", "sample_accession"],
                                  kind="stable", key=lambda column: column.map(C.accession_order) if column.name in ("study_accession", "sample_accession")
                                  else column)
    chosen = set(primary[primary["patient_identified"]].groupby("patient_key")["sample_id"].first()) | set(primary.loc[~primary["patient_identified"], "sample_id"])
    table["one_per_patient"] = table["is_primary_record"] & table["sample_id"].isin(chosen)
    per_patient = primary[primary["patient_identified"]].groupby("patient_key")["sample_id"].nunique()
    table["samples_of_patient"] = table["patient_key"].map(per_patient).fillna(0).astype(int)
    table = table[SAMPLE_COLUMNS].sort_values(["source", "study_accession", "platform", "sample_accession"], kind="stable").reset_index(drop=True)
    C.require(int(table["one_per_patient"].sum()) == table.loc[table["patient_identified"], "patient_key"].nunique()
              + int((table["is_primary_record"] & ~table["patient_identified"]).sum()), "One sample must be marked for every patient")

    groups = pd.read_parquet(C.DATA / "cohort_groups.parquet")
    scales = pd.read_parquet(C.ANNOTATION / "matrix_index.parquet", columns=["accession", "platform", "value_scale"]).rename(columns={"accession": "study_accession"})
    rows = []
    for (study, platform), part in table.groupby(["study_accession", "platform"], sort=True):
        first = part[part["is_primary_record"]]
        rows.append(dict(cohort_group=part["cohort_group"].iloc[0], patient_safe_group=part["patient_safe_group"].iloc[0], source=part["source"].iloc[0],
                         study_accession=study, platform=platform, records=len(part), primary_records=len(first), repeated_records=len(part) - len(first),
                         **class_counts(part), **class_counts(first, "_primary"), infection_sepsis_not_stated=int(part["infection_sepsis_not_stated"].sum()),
                         patients_identified=int(part.loc[part["patient_identified"], "patient_key"].nunique()),
                         records_without_patient=int((~part["patient_identified"]).sum()), one_per_patient=int(part["one_per_patient"].sum()),
                         sepsis_before_onset=int((part["sepsis_phase"] == "before_onset").sum()), sepsis_onset_not_recorded=int((part["sepsis_phase"] == "onset_not_recorded").sum()),
                         recovery_or_follow_up=int((part["sepsis_phase"] == "recovery_or_follow_up").sum()), sepsis_definition=part["sepsis_definition"].iloc[0],
                         tissues=shares(part["tissue"]), age_groups=shares(part["age_group"]), sex_recorded=int((part["sex"] != "not_recorded").sum()),
                         time_point_recorded=int(part["time_order"].notna().sum()), outcome_recorded=int((part["outcome"] != "not_recorded").sum()),
                         septic_shock_recorded=int((part["septic_shock"] != "not_recorded").sum()), flagged=int(part["any_flag"].sum()),
                         class_differs_between_records=int(part["class_differs_between_records"].sum())))
    cohorts = pd.DataFrame(rows).merge(groups[["study_accession", "platform", "release_date", "shares_samples_with"]], on=["study_accession", "platform"],
                                       validate="one_to_one").merge(scales, on=["study_accession", "platform"], how="left", validate="one_to_one")
    cohorts["value_scale"] = cohorts["value_scale"].fillna("")
    rows = []
    for name, part in table[table["is_primary_record"]].groupby("cohort_group", sort=True):
        tables = sorted(set(zip(table.loc[table["cohort_group"] == name, "study_accession"], table.loc[table["cohort_group"] == name, "platform"])))
        rows.append(dict(cohort_group=name, patient_safe_group=part["patient_safe_group"].iloc[0],
                         study_platforms="; ".join(f"{study}|{platform}" for study, platform in tables), distinct_samples=len(part), **class_counts(part),
                         patients_identified=int(part.loc[part["patient_identified"], "patient_key"].nunique()),
                         samples_without_patient=int((~part["patient_identified"]).sum()), one_per_patient=int(part["one_per_patient"].sum()),
                         tissues=shares(part["tissue"]), age_groups=shares(part["age_group"]), flagged=int(part["any_flag"].sum())))
    summary = pd.DataFrame(rows)
    C.save_table(table, C.DATA / "cohort_samples.parquet")
    C.save_table(cohorts, C.DATA / "cohort_summary.parquet")
    C.save_table(summary, C.DATA / "cohort_group_summary.parquet")
    C.save_review(cohorts, C.OUT / "cohort_summary_review.csv")

    distinct = table[table["is_primary_record"]]
    unique = table[table["one_per_patient"]]
    sepsis = distinct[distinct["sample_class"] == "sepsis"]
    C.finish(5, start, [C.DATA / "cohort_samples.parquet", C.DATA / "cohort_summary.parquet", C.DATA / "cohort_group_summary.parquet",
                        C.OUT / "cohort_summary_review.csv"], method=C.CONFIG["cohort_tables"]["method"],
             records=dict(total=len(table), **class_counts(table)),
             distinct_samples=dict(total=len(distinct), **class_counts(distinct),
                                   whose_records_differ_in_class=int(distinct["class_differs_between_records"].sum())),
             patients=dict(identified=int(table.loc[table["patient_identified"], "patient_key"].nunique()),
                           distinct_samples_with_an_identified_patient=int(distinct["patient_identified"].sum()),
                           distinct_samples_without_a_patient_identifier=int((~distinct["patient_identified"]).sum()),
                           one_per_patient=dict(total=len(unique), **class_counts(unique)),
                           identified_with_a_sepsis_sample=int(sepsis.loc[sepsis["patient_identified"], "patient_key"].nunique()),
                           sepsis_samples_without_a_patient_identifier=int((~sepsis["patient_identified"]).sum()),
                           most_samples_of_one_patient=int(table["samples_of_patient"].max()),
                           in_several_studies=int((shared > 1).sum()), in_several_cohort_groups=int((spanning > 1).sum()),
                           with_conflicting_recorded_sex=len(conflicting_sex),
                           conflicting_recorded_sex=conflicting_sex[:50], namespaces=C.REGISTER["patient_namespaces"],
                           replicate_pairs=len(replicates)),
             sepsis_by_phase={str(name): int(count) for name, count in sepsis["sepsis_phase"].value_counts().sort_index().items()},
             cohorts=dict(study_platforms=len(cohorts), cohort_groups=len(summary), patient_safe_groups=int(table["patient_safe_group"].nunique()),
                          joined_by_patients=sorted(name for name in set(table["patient_safe_group"]) if "+" in name),
                          groups_with_sepsis=int((summary["sepsis"] > 0).sum()), groups_with_healthy_controls=int((summary["healthy_control"] > 0).sum()),
                          groups_with_sepsis_and_healthy_controls=int(((summary["sepsis"] > 0) & (summary["healthy_control"] > 0)).sum())),
             flagged=dict(records=int(table["any_flag"].sum()), distinct_samples=int(distinct["any_flag"].sum())), protocol=C.CONFIG["protocol"])
    print(f"[05] Complete: {len(table):,} records with expression values, {len(distinct):,} distinct samples; "
          f"{int(table.loc[table['patient_identified'], 'patient_key'].nunique()):,} patients identified and {int((~distinct['patient_identified']).sum()):,} "
          f"distinct samples without a patient identifier; one sample per patient: {len(unique):,}; {len(summary)} cohort groups form "
          f"{table['patient_safe_group'].nunique()} patient-safe groups; "
          + ", ".join(f"{name} {count:,}" for name, count in class_counts(distinct).items())
          + f" among the distinct samples; review sheet: {(C.OUT / 'cohort_summary_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
