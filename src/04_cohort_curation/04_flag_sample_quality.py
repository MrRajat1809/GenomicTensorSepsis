"""Flag the samples that the processing statistics or the records mark as doubtful, without removing any sample."""
import importlib
import re

C = importlib.import_module("01_catalogue_sample_metadata")

FLAGS = ("low_library_size", "few_genes_detected", "low_detection", "low_dynamic_range", "flagged_in_record", "technical_replicate",
         "identity_uncertain")
QUALITY_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "library_size", "genes_detected", "detected_share", "dynamic_range",
                   "dynamic_range_to_study_median", *FLAGS, "record_flag", "replicate_of", "replicate_basis", "identity_basis",
                   "quality_flags", "any_flag"]


def main():
    import numpy as np
    import pandas as pd

    print(f"[04] Flagging sample quality v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["quality"]
    keys = ["study_accession", "sample_accession"]
    table = C.sample_records(["source", "study_accession", "sample_accession", "platform"])
    counts = pd.read_parquet(C.EXPRESSION / "count_table_samples.parquet", columns=[*keys, "platform", "library_size", "genes_detected"])
    arrays = pd.read_parquet(C.EXPRESSION / "intensity_table_arrays.parquet", columns=[*keys, "platform", "detected_share", "signal_median", "signal_percentile_99"])

    # Every measure is compared with the samples of its own study platform, never across studies.
    counts["genes_to_median"] = counts["genes_detected"] / counts.groupby(["study_accession", "platform"])["genes_detected"].transform("median")
    # A ratio to a median at or below zero (a background-subtracted array) is no dynamic range; it stays missing and is not flagged.
    positive = arrays["signal_median"] > 0
    arrays["dynamic_range"] = arrays["signal_percentile_99"] / arrays["signal_median"].where(positive)
    groups = arrays.groupby(["study_accession", "platform"])
    arrays["dynamic_range_to_study_median"] = arrays["dynamic_range"] / groups["dynamic_range"].transform("median")
    arrays["detected_to_median"] = arrays["detected_share"] / groups["detected_share"].transform("median")
    C.require(not counts.duplicated(keys).any() and not arrays.duplicated(keys).any(), "A sample has two rows of processing statistics")
    table = table.merge(counts.drop(columns="platform"), on=keys, how="left", validate="one_to_one").merge(arrays.drop(columns="platform"), on=keys, how="left",
                                                                                                         validate="one_to_one")
    table["low_library_size"] = table["library_size"] < policy["minimum_library_size"]
    table["few_genes_detected"] = table["genes_to_median"] < policy["genes_detected_below_share_of_study_median"]
    table["low_detection"] = table["detected_to_median"] < policy["detected_share_below_share_of_study_median"]
    table["low_dynamic_range"] = table["dynamic_range_to_study_median"] < policy["dynamic_range_below_share_of_study_median"]

    metadata = pd.read_parquet(C.DATA / "sample_metadata.parquet", columns=[*keys, "platform", "field", "value"])
    recorded = {}
    for rule in C.REGISTER["record_flags"]:
        part = metadata[(metadata["study_accession"] == rule["study"]) & (metadata["field"] == rule["field"])]
        if "platform" in rule:
            part = part[part["platform"] == rule["platform"]]
        hits = part[part["value"].map(lambda value: bool(re.search(rule["pattern"], str(value), re.IGNORECASE)))]
        C.require(len(hits), f"A record flag is met by no sample: {rule['study']} {rule['field']} {rule['pattern']}")
        for study, sample in hits[keys].drop_duplicates().itertuples(index=False, name=None):
            recorded[(study, sample)] = "; ".join(filter(None, [recorded.get((study, sample), ""), rule["flag"]]))
    table["record_flag"] = [recorded.get(key, "") for key in zip(table["study_accession"], table["sample_accession"])]
    table["flagged_in_record"] = table["record_flag"] != ""

    # Technical replicates: of the replicate arrays that a study records for one sample, all but the first by accession; and the later library of every
    # replicate pair that the link stage found by its values.
    replicate = {}
    for rule in C.REGISTER["replicate_groups"]:
        part = metadata[(metadata["study_accession"] == rule["study"]) & (metadata["field"] == rule["field"])].drop_duplicates(keys)
        named = part.assign(group=part["value"].str.extract(rule["pattern"], expand=False)).dropna(subset=["group"])
        C.require(len(named) and named["group"].duplicated().any(), f"A replicate rule is met by no group of arrays: {rule['study']} {rule['pattern']}")
        for _, members in named.sort_values("sample_accession", key=lambda column: column.map(C.accession_order)).groupby("group"):
            for sample in members["sample_accession"].iloc[1:]:
                replicate[(rule["study"], sample)] = (members["sample_accession"].iloc[0], rule["flag"])
    libraries = pd.read_parquet(C.DATA / "replicate_libraries.parquet")
    for row in libraries.itertuples(index=False):
        replicate[(row.study_accession, row.second)] = (row.first, policy["replicate_flag"])
    table["replicate_of"] = [replicate.get(key, ("", ""))[0] for key in zip(table["study_accession"], table["sample_accession"])]
    table["replicate_basis"] = [replicate.get(key, ("", ""))[1] for key in zip(table["study_accession"], table["sample_accession"])]
    table["technical_replicate"] = table["replicate_of"] != ""
    uncertain = {(rule["study"], sample): f"{rule['flag']}: {rule['basis']}" for rule in C.REGISTER["identity_flags"] for sample in rule["samples"]}
    C.require(set(uncertain) <= set(zip(table["study_accession"], table["sample_accession"])),
              "An identity flag names a record that is not an included record with expression values")
    table["identity_basis"] = [uncertain.get(key, "") for key in zip(table["study_accession"], table["sample_accession"])]
    table["identity_uncertain"] = table["identity_basis"] != ""
    table["quality_flags"] = table[list(FLAGS)].apply(lambda row: "; ".join(name for name in FLAGS if row[name]), axis=1)
    table["any_flag"] = table["quality_flags"] != ""
    table = table[QUALITY_COLUMNS]
    C.require(len(table) == len(pd.read_parquet(C.DATA / "sample_links.parquet", columns=["sample_accession"])) and not table.duplicated(keys).any(),
              "Every sample needs one quality row")
    C.save_table(table, C.DATA / "sample_quality.parquet")
    flagged = table[table["any_flag"]]
    C.save_review(flagged, C.OUT / "flagged_samples_review.csv")
    by_flag = {name: int(table[name].sum()) for name in FLAGS}
    by_study = {name: {str(study): int(count) for study, count in table[table[name]].groupby("study_accession").size().items()} for name in FLAGS}
    C.finish(4, start, [C.DATA / "sample_quality.parquet", C.OUT / "flagged_samples_review.csv"], method=policy["method"], samples=len(table),
             flagged=dict(samples=len(flagged), by_flag=by_flag, by_flag_and_study=by_study),
             measured=dict(library_size=int(table["library_size"].notna().sum()), detected_share=int(table["detected_share"].notna().sum()),
                           dynamic_range=int(table["dynamic_range"].notna().sum()), dynamic_range_without_a_positive_median=int((~positive).sum())),
             thresholds={name: value for name, value in policy.items() if isinstance(value, (int, float))},
             record_flags=len(C.REGISTER["record_flags"]), replicate_groups=C.REGISTER["replicate_groups"],
             technical_replicates={str(study): int(count) for study, count in table[table["technical_replicate"]].groupby("study_accession").size().items()},
             protocol=C.CONFIG["protocol"])
    print(f"[04] Complete: {len(flagged):,} of {len(table):,} samples carry a flag; "
          + ", ".join(f"{name} {count:,}" for name, count in by_flag.items())
          + f"; review sheet: {(C.OUT / 'flagged_samples_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
