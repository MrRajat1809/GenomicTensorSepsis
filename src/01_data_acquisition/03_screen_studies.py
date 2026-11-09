"""Apply the recorded screening decision to every candidate study and verify each decision against the harvested evidence."""
import importlib

C = importlib.import_module("01_search_repositories")

REVIEW_COLUMNS = ["accession", "source", "decision", "reason", "flag", "basis", "overlaps_with", "title", "assay", "declared_assay", "samples",
                  "human_samples", "single_cell_samples", "release_date", "access", "harvest_status", "url"]
ASSAYS = ("array", "rna-seq")


def sample_technology(samples):
    """Label every sample record array or rna-seq by what it records itself, or leave it empty where it records neither."""
    policy = C.CONFIG["technology"]
    geo = samples["source"] == "geo"
    unknown = sorted(set(samples.loc[geo, "sample_type"]) - set(policy["geo_sample_types"]))
    C.require(not unknown, f"GEO sample types without a technology: {unknown}; map each in technology.geo_sample_types of acquisition_config.json")
    unknown = sorted(set(samples["library_strategy"]) - set(policy["sequencing_library_strategies"]) - {""})
    C.require(not unknown, f"Library strategies without a technology: {unknown}; list those of sequenced samples in "
                           "technology.sequencing_library_strategies of acquisition_config.json")
    sequenced = (samples["library_strategy"] != "") | (samples["instrument_model"] != "")
    designed = (samples["source"] == "arrayexpress") & (samples["platform"] != "")
    typed = samples["sample_type"].map(policy["geo_sample_types"])
    # A GEO sample names its kind in its sample type; a library strategy or instrument, where it records one, must agree with that type.
    conflict = (geo & ((typed == "rna-seq") != sequenced)) | (designed & sequenced)
    found = [f"{study}:{sample}" for study, sample in samples.loc[conflict, ["study_accession", "sample_accession"]].head(10).itertuples(index=False)]
    C.require(not conflict.any(), f"{int(conflict.sum())} sample records whose array and sequencing fields disagree, e.g. {found}; review them and "
                                  "extend the technology rules of acquisition_config.json")
    return sequenced.map({True: "rna-seq", False: ""}).mask(designed, "array").mask(geo, typed)


def study_technology(table, samples):
    """Name the technology of every candidate from its sample records, joined as array+rna-seq for both; keep the declared type where they name none."""
    found = samples.assign(label=sample_technology(samples)).groupby(["source", "study_accession"])["label"].agg(set).to_dict()
    assays, bases = [], []
    for row in table.itertuples(index=False):
        labels = found.get((row.source, row.accession), set())
        C.require(not ("" in labels and len(labels) > 1), f"{row.accession}: only some of its sample records name their technology; review them and "
                                                           "extend the technology rules of acquisition_config.json")
        named = [label for label in ASSAYS if label in labels]
        assays.append("+".join(named) if named else row.declared_assay)
        bases.append("sample_records" if named else "declared")
    return assays, bases


def decision_problems(row, decision, reasons):
    """Return every way one recorded decision contradicts its own fields or the harvested facts."""
    found = []
    reason, included = decision.get("reason"), decision.get("decision") == "include"
    if decision.get("decision") not in ("include", "exclude") or reason not in reasons or not str(decision.get("basis", "")).strip():
        return ["decision must be include or exclude with a known reason and a basis"]
    if decision.get("flag", "uncertain") != "uncertain":
        found.append("flag may only be uncertain")
    if included != (reason == "eligible"):
        found.append("include and the reason eligible must go together")
    if included and not (row.access == "open" and row.source in ("geo", "arrayexpress") and row.harvest_status == "harvested"
                         and row.human_samples > row.single_cell_samples and row.assay and row.assay_basis == "sample_records"
                         and (row.source == "geo" or row.data_files)):
        found.append("an included study needs open access, harvested records, a technology named by its sample records, human bulk samples and a data file")
    if (reason == "controlled_access") != (row.access != "open" and reason not in reasons[:reasons.index("controlled_access")]):
        found.append("a study outside open access needs a scientific reason or controlled_access, and controlled_access needs a study outside open access")
    if reason == "no_public_reads" and row.harvest_status != "harvested_no_public_reads":
        found.append("no_public_reads requires an ENA study without public reads")
    unfiled = row.source == "arrayexpress" and not row.data_files
    # The series matrix of a GEO sequencing series holds no values, so without supplementary files it may hold reads only; an NCBI-generated count
    # table is not in the harvested records, and the download plan stops at an included study without any file of values.
    unfiled_series = row.source == "geo" and row.assay == "rna-seq" and not row.data_files
    if reason == "raw_reads_only" and not ((row.source == "ena" and row.raw_bytes > 0) or (unfiled and row.assay == "rna-seq") or unfiled_series):
        found.append("raw_reads_only requires public reads and no expression file in a repository")
    if reason == "no_expression_files" and not unfiled:
        found.append("no_expression_files requires an ArrayExpress study that declares no data file")
    # Whether deposited files can be read into comparable values is found at the expression stage; here the study must at least hold such files.
    if reason == "no_usable_expression_values" and not (row.access == "open" and row.source in ("geo", "arrayexpress") and row.harvest_status == "harvested"
                                                        and row.human_samples > row.single_cell_samples and row.data_files):
        found.append("no_usable_expression_values requires an open GEO or ArrayExpress study with human bulk samples and deposited expression files")
    if (reason == "third_party_reanalysis") != (row.harvest_status == "reanalysis_without_own_samples"):
        found.append("third_party_reanalysis is for, and only for, series without samples of their own")
    if reason == "not_human" and row.human_samples > 0:
        found.append("not_human contradicts the human samples on record")
    if reason == "not_a_dataset" and row.harvest_status != "catalogue_only":
        found.append("not_a_dataset is for catalogue entries without a study record")
    return found


def main():
    import pandas as pd

    print(f"[03] Screening candidate studies v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy = C.CONFIG["screening"]
    reasons, decisions = list(policy["reasons"]), C.REGISTER
    candidates = pd.read_parquet(C.DATA / "candidates.parquet")
    studies = pd.read_parquet(C.DATA / "studies.parquet")
    samples = pd.read_parquet(C.DATA / "samples.parquet")
    unscreened = sorted(set(candidates["accession"]) - set(decisions))
    C.require(not unscreened, f"{len(unscreened)} candidates without a screening decision: {unscreened[:40]}; add a decision for each to "
                              f"{C.CONFIG['register']}, and move search.census_date only if you intend to extend the resource")
    unknown = sorted(set(decisions) - set(candidates["accession"]))
    C.require(not unknown, f"Decisions recorded for studies that are not candidates: {unknown[:20]}")

    facts = samples.assign(human=samples["organism"].str.contains("Homo sapiens", case=False),
                           single=samples["library_source"].str.contains("single cell", case=False))
    facts = facts.groupby(["source", "study_accession"], as_index=False).agg(
        samples=("sample_accession", "size"), human_samples=("human", "sum"), single_cell_samples=("single", "sum"))
    table = candidates[["source", "accession", "title", "declared_assay", "access", "release_date", "also_in", "url"]].merge(
        studies[["source", "accession", "harvest_status", "recorded_samples", "raw_bytes", "superseries_of", "subseries_of", "reanalysis_of",
                 "supplementary_files"]],
        on=["source", "accession"], validate="one_to_one").merge(
        facts.rename(columns={"study_accession": "accession"}), on=["source", "accession"], how="left", validate="one_to_one")
    for column in ("samples", "human_samples", "single_cell_samples"):
        table[column] = table[column].fillna(0).astype(int)
    C.require(len(table) == len(candidates), "Screening table does not reconcile with the candidates")
    assays, bases = study_technology(table, samples)
    table.insert(3, "assay", assays)
    table.insert(5, "assay_basis", bases)
    # ArrayExpress declares every file of a study; the sample tables alone are not expression data.
    table["data_files"] = [sum(not name.lower().endswith((".idf.txt", ".sdrf.txt")) for name in files.split(";") if name)
                           for files in table.pop("supplementary_files")]
    for column in ("decision", "reason", "basis", "overlaps_with", "flag"):
        table[column] = [str(decisions[accession].get(column, "")) for accession in table["accession"]]

    problems = [f"{row.accession}: {problem}" for row in table.itertuples(index=False)
                for problem in decision_problems(row, decisions[row.accession], reasons)]
    included = set(table.loc[table["decision"] == "include", "accession"])
    members = samples.groupby("study_accession")["sample_accession"].apply(set).to_dict()
    holders = {}
    for accession in sorted(included):
        for sample in members.get(accession, ()):
            holders.setdefault(sample, []).append(accession)
    shared = sorted({tuple(found) for found in holders.values() if len(found) > 1})
    problems += [f"{' and '.join(pair)}: included studies share samples; keep one and exclude the other as overlaps_included_study" for pair in shared]
    overlap = []
    for row in table.itertuples(index=False):
        own = members.get(row.accession, set())
        covered = sum(sample in holders for sample in own)
        if row.reason == "superseries_of_included_subseries" and not (row.superseries_of and own and covered == len(own)):
            problems.append(f"{row.accession}: not every sample of this SuperSeries belongs to an included study")
        if row.reason == "overlaps_included_study" and not (row.overlaps_with in included and own & members.get(row.overlaps_with, set())):
            problems.append(f"{row.accession}: overlaps_with must name an included study that shares samples with it")
        overlap.append(covered if row.decision == "exclude" else 0)
    table["samples_in_included_studies"] = overlap
    C.require(not problems, f"{len(problems)} screening decisions contradict the evidence:\n  " + "\n  ".join(problems[:40]))

    table["reason"] = pd.Categorical(table["reason"], categories=reasons, ordered=True)
    table = table.sort_values(["reason", "source", "accession"]).reset_index(drop=True)
    table["reason"] = table["reason"].astype(str)
    chosen = table[table["decision"] == "include"].reset_index(drop=True)
    C.require(len(chosen), "No study passed screening")
    C.save_table(table, C.DATA / "study_decisions.parquet")
    C.save_table(chosen, C.DATA / "included_studies.parquet")
    review = C.OUT / "screening_review.csv"
    temporary = review.with_name(review.name + ".tmp")
    table[REVIEW_COLUMNS].to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(review)

    attrition, remaining = [dict(step="candidates", removed=0, remaining=len(table))], len(table)
    for reason in reasons:
        if reason != "eligible":
            part = table[table["reason"] == reason]
            remaining -= len(part)
            attrition.append(dict(step=reason, removed=len(part), remaining=remaining, sample_records=int(part["samples"].sum()),
                                  recorded_samples=int(part["recorded_samples"].sum())))
    C.require(remaining == len(chosen), "Screening attrition does not reconcile with the included studies")
    documented = {reason: dict(studies=int((table["reason"] == reason).sum()),
                               recorded_samples=int(table.loc[table["reason"] == reason, "recorded_samples"].sum()),
                               raw_bytes=int(table.loc[table["reason"] == reason, "raw_bytes"].sum()))
                  for reason in ("controlled_access", "no_public_reads", "raw_reads_only", "no_expression_files", "no_usable_expression_values")}
    changed = table[table["assay"] != table["declared_assay"]]
    C.finish(3, start, [C.DATA / "study_decisions.parquet", C.DATA / "included_studies.parquet", review],
             method=policy["method"], attrition=attrition,
             technology=dict(rules=C.CONFIG["technology"], by_basis=C.counts(table, "assay_basis"),
                             differs_from_declared=[dict(accession=row.accession, declared_assay=row.declared_assay, assay=row.assay, decision=row.decision)
                                                    for row in changed.itertuples(index=False)]),
             included=dict(studies=len(chosen), by_source=C.counts(chosen, "source"), by_assay=C.counts(chosen, "assay"),
                           sample_records=int(chosen["samples"].sum()), human_sample_records=int(chosen["human_samples"].sum())),
             eligible_but_unavailable=documented, flagged_uncertain=dict(total=int((table["flag"] == "uncertain").sum()),
                 included=int((chosen["flag"] == "uncertain").sum())), eligibility=C.CONFIG["eligibility"])
    print(f"[03] Complete: {len(chosen):,} of {len(table):,} candidates included; review sheet: {review.relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
