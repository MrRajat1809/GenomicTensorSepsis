"""Verify provenance, replay the census, records and screening, and check every downloaded file of the acquisition stage."""
import gzip
import importlib
import re
import tarfile
import traceback
import zipfile
import zlib
from datetime import date, datetime, timezone
from pathlib import PurePosixPath

import pandas as pd

C = importlib.import_module("01_search_repositories")
H = importlib.import_module("02_harvest_study_metadata")
S = importlib.import_module("03_screen_studies")
D = importlib.import_module("04_download_raw_data")

BLOCK = 4 * 1024 * 1024
MEANING = ("PASS means that the census and the study and sample tables are reproduced from the stored repository responses, that no study released "
           "after the census date is a candidate, that every screening decision agrees with those tables, and that every downloaded file is present at "
           "its recorded size and SHA-256 with readable archives, as the download lock records it. "
           "It does not assess the scientific quality of a study or the content of an expression file.")


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


def offline(url, params=None, post=False):
    raise ValueError(f"Quality control reads stored records only, and this one is missing: {url}")


def checksums(folder):
    """Read a SHA256SUMS.txt as {file: digest}."""
    lines = (folder / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    return {name: digest for digest, name in (line.split("  ", 1) for line in lines)}


def stored_files(folder):
    return {path.relative_to(folder).as_posix() for path in folder.rglob("*") if path.is_file()} - {"SHA256SUMS.txt"}


def day_of(text, formats):
    """Read a release date with the first format that fits; empty where there is none, None where none fits."""
    text = str(text).strip()
    return next((found for found in (attempt(text, form) for form in formats) if found), None) if text else ""


def attempt(text, form):
    try:
        return datetime.strptime(text, form).date()
    except ValueError:
        return None


def replayed(qc, name, built, table):
    pd.testing.assert_frame_equal(built.reset_index(drop=True), pd.read_parquet(C.DATA / table))
    qc.check(name, True)


def provenance_checks(qc, state):
    state["start"] = C.begin(5)
    qc.check("all_stage_receipts_and_code_bindings", True)


def census_checks(qc, state):
    report, priority = C.read_json(C.receipt_path(1)), C.CONFIG["deduplication"]["source_priority"]
    search = C.saved_search()
    qc.check("search_settings_match_the_snapshot", search["search_sha256"] == report["search_sha256"] and search["queries"] == report["queries"])
    qc.check("search_responses_unchanged", C.verify_snapshots() == report["snapshot_files"])
    terms = {term.replace(" ", "_"): term for term in C.CONFIG["search"]["terms"]}
    rows = C.geo_records(terms) + C.arrayexpress_records(terms) + C.ena_records(terms) + C.catalogue_records(terms)
    for row in rows:
        row.update(duplicate_of="", is_superseries=C.CONFIG["search"]["geo"]["superseries_marker"] in row["description"])
    C.apply_cutoff(rows)
    owner = C.deduplicate(rows)
    rows.sort(key=lambda row: (priority.index(row["source"]), row["accession"]))
    census = pd.DataFrame(rows, columns=C.COLUMNS)
    replayed(qc, "census_replayed_from_search_responses", census, "census.parquet")
    candidates = census[census["census_status"] == "candidate"].reset_index(drop=True)
    replayed(qc, "candidates_replayed_from_census", candidates, "candidates.parquet")
    qc.check("candidate_accessions_unique", candidates["accession"].is_unique)
    duplicate = census["census_status"] == "cross_source_duplicate"
    qc.check("every_census_record_has_a_status", (census["census_status"] != "").all())
    qc.check("every_duplicate_names_a_kept_record", census.loc[duplicate, "duplicate_of"].isin(candidates["accession"]).all()
             and (census.loc[~duplicate, "duplicate_of"] == "").all())
    attrition = [dict(source=source, records=int((census["source"] == source).sum()), **C.counts(census[census["source"] == source], "census_status"))
                 for source in priority]
    qc.check("census_attrition_reconciles", attrition == report["attrition"] and report["candidates"]["total"] == len(candidates))
    cutoff, dates = C.CONFIG["search"]["census_date"], C.CONFIG["search"]["release_dates"]
    stored = search.get("census_date") or search["searched_at_utc"][:10]
    qc.check("census_date_is_that_of_the_stored_search", cutoff == stored == report["census_date"], dict(config=cutoff, stored_search=stored))
    days = [day_of(text, dates["formats"]) for text in census["release_date"]]
    limit = date.fromisoformat(cutoff)
    late = census["census_status"] == "released_after_census_date"
    qc.check("every_release_date_read", all(day is not None for day, status in zip(days, census["census_status"])
                                             if status in ("candidate", "released_after_census_date")))
    chosen = census["census_status"] == "candidate"
    qc.check("no_candidate_released_after_the_census_date",
             all((day != "" and day <= limit) or (day == "" and source in dates["undated"])
                 for day, source, keep in zip(days, census["source"], chosen) if keep)
             and all(day != "" and day > limit for day, flagged in zip(days, late) if flagged)
             and int(late.sum()) == report["released_after_census_date"], dict(released_after_census_date=int(late.sum())))
    recall = [accession for accessions in C.CONFIG["recall"]["studies"].values() for accession in accessions]
    missed = [accession for accession in recall if accession not in owner]
    qc.check("recall_studies_reached", not missed and report["recall"]["found"] == len(recall), missed)
    state.update(candidates=candidates, summary=dict(census_records=len(census), candidates=len(candidates), recall_studies=len(recall)))


def record_checks(qc, state):
    report, candidates = C.read_json(C.receipt_path(2)), state["candidates"]
    expected = checksums(H.RAW)
    qc.check("stored_records_complete", set(expected) == stored_files(H.RAW) and len(expected) == report["raw_files"])
    qc.check("stored_records_unchanged", all(C.sha(H.RAW / name) == digest for name, digest in expected.items()))
    studies, samples, attributes = [], [], []
    part = lambda source: candidates[candidates["source"] == source]
    H.harvest_geo(part("geo")["accession"].tolist(), False, studies, samples, attributes)
    H.harvest_arrayexpress(list(part("arrayexpress")[["accession", "title", "description", "platform"]].itertuples(index=False, name=None)),
                           False, studies, samples, attributes)
    H.harvest_ena(part("ena")["accession"].tolist(), False, studies, samples, attributes)
    catalogued = candidates[candidates["access"] != "open"]
    H.harvest_catalogue(list(catalogued[["source", "accession", "title", "description"]].itertuples(index=False, name=None)), False, studies)
    qc.check("stored_records_are_the_records_the_harvest_reads", {path.relative_to(H.RAW).as_posix() for path in H.RETRIEVED} == set(expected))
    studies = candidates[["source", "accession"]].merge(pd.DataFrame(studies, columns=H.STUDY_COLUMNS), on=["source", "accession"],
                                                        how="left", validate="one_to_one")[H.STUDY_COLUMNS].copy()
    for column, origin in (("title", "title"), ("summary", "description")):
        studies[column] = studies[column].where(studies[column] != "", candidates[origin].to_numpy())
    samples = pd.DataFrame(samples, columns=H.SAMPLE_COLUMNS)
    attributes = pd.DataFrame(attributes, columns=H.ATTRIBUTE_COLUMNS)
    replayed(qc, "study_table_replayed_from_stored_records", studies, "studies.parquet")
    replayed(qc, "sample_table_replayed_from_stored_records", samples, "samples.parquet")
    replayed(qc, "attribute_table_replayed_from_stored_records", attributes, "sample_attributes.parquet")
    qc.check("one_study_record_per_candidate", studies["harvest_status"].notna().all() and len(studies) == len(candidates))
    key = ["source", "study_accession", "sample_accession"]
    qc.check("samples_unique_within_a_study", not samples.duplicated(key).any())
    qc.check("every_sample_belongs_to_a_candidate", set(samples["study_accession"]) <= set(candidates["accession"]))
    qc.check("every_attribute_belongs_to_a_sample", set(map(tuple, attributes[key].drop_duplicates().to_numpy())) <= set(map(tuple, samples[key].to_numpy())))
    qc.check("harvest_summary_reconciles", sum(item["studies"] for item in report["summary"]) == len(studies)
             and sum(item["samples"] for item in report["summary"]) == len(samples)
             and sum(item["attributes"] for item in report["summary"]) == len(attributes))
    state["samples"] = samples
    state["summary"].update(stored_records=len(expected), sample_records=len(samples), attribute_values=len(attributes))


def screening_checks(qc, state):
    report, policy, samples = C.read_json(C.receipt_path(3)), C.CONFIG["screening"], state["samples"]
    reasons, decisions = list(policy["reasons"]), C.REGISTER
    table = pd.read_parquet(C.DATA / "study_decisions.parquet")
    included = pd.read_parquet(C.DATA / "included_studies.parquet")
    qc.check("one_decision_per_candidate", sorted(table["accession"]) == sorted(state["candidates"]["accession"]) == sorted(decisions))
    qc.check("decisions_match_the_register", all(str(decisions[row.accession].get(column, "")) == getattr(row, column)
             for row in table.itertuples(index=False) for column in ("decision", "reason", "basis", "overlaps_with", "flag")))
    facts = samples.assign(human=samples["organism"].str.contains("Homo sapiens", case=False),
                           single=samples["library_source"].str.contains("single cell", case=False))
    facts = facts.groupby("study_accession").agg(samples=("sample_accession", "size"), human_samples=("human", "sum"),
                                                 single_cell_samples=("single", "sum")).reindex(table["accession"]).fillna(0).astype(int)
    qc.check("sample_counts_replayed", all((facts[column].to_numpy() == table[column].to_numpy()).all() for column in facts.columns))
    rules, geo = C.CONFIG["technology"], samples["source"] == "geo"
    qc.check("sample_types_and_library_strategies_classified", samples.loc[geo, "sample_type"].isin(list(rules["geo_sample_types"])).all()
             and samples["library_strategy"].isin(rules["sequencing_library_strategies"] + [""]).all())
    sequenced = samples["library_strategy"].ne("") | samples["instrument_model"].ne("")
    kind = samples["platform"].ne("").map({True: "array", False: ""}).mask(sequenced, "rna-seq")
    kind = kind.mask(geo, samples["sample_type"].map(rules["geo_sample_types"]))
    named = kind.groupby(samples["study_accession"]).agg(lambda labels: "+".join(label for label in ("array", "rna-seq") if label in set(labels)))
    named = table["accession"].map(named).fillna("")
    declared = dict(zip(state["candidates"]["accession"], state["candidates"]["declared_assay"]))
    qc.check("declared_assay_is_the_census_declaration", (table["declared_assay"] == table["accession"].map(declared)).all())
    qc.check("technology_replayed_from_the_sample_records", (table["assay"] == named.where(named != "", table["declared_assay"])).all()
             and (table["assay_basis"] == named.ne("").map({True: "sample_records", False: "declared"})).all())
    problems = [f"{row.accession}: {problem}" for row in table.itertuples(index=False)
                for problem in S.decision_problems(row, decisions[row.accession], reasons)]
    qc.check("decisions_agree_with_the_records", not problems, problems)
    chosen = table[table["decision"] == "include"].reset_index(drop=True)
    pd.testing.assert_frame_equal(chosen, included)
    qc.check("included_studies_are_the_include_decisions", (chosen["reason"] == "eligible").all() and len(chosen) == report["included"]["studies"])
    qc.check("included_technology_counts_reconcile", report["included"]["by_assay"] == C.counts(chosen, "assay")
             and set(chosen["assay_basis"]) == {"sample_records"})
    members = samples[samples["study_accession"].isin(chosen["accession"])]
    qc.check("no_sample_in_two_included_studies", not members["sample_accession"].duplicated().any())
    qc.check("included_sample_records_reconcile", len(members) == int(chosen["samples"].sum()) == report["included"]["sample_records"])
    removed = {reason: int((table["reason"] == reason).sum()) for reason in reasons if reason != "eligible"}
    qc.check("screening_attrition_reconciles", {item["step"]: item["removed"] for item in report["attrition"][1:]} == removed
             and report["attrition"][-1]["remaining"] == len(chosen) == len(table) - sum(removed.values()))
    review = pd.read_csv(C.OUT / "screening_review.csv", dtype=str, keep_default_na=False)
    qc.check("review_sheet_matches_the_decisions", review.equals(table[S.REVIEW_COLUMNS].astype(str)))
    state["included"] = included
    state["summary"].update(included_studies=len(chosen), included_sample_records=len(members), exclusions=removed)


def gzip_intact(stream):
    """Read a gzip stream to its end, which verifies the checksum and length of every member."""
    try:
        with gzip.GzipFile(fileobj=stream) as content:
            while content.read(BLOCK):
                pass
        return True
    except (OSError, EOFError, zlib.error):
        return False


def archive_members(path):
    """List a tar archive as {member: size} and name the gzip members that do not decompress."""
    members, damaged = {}, []
    with tarfile.open(path) as archive:
        for member in archive:
            if member.isfile():
                members[member.name] = member.size
                if member.name.lower().endswith(".gz") and not gzip_intact(archive.extractfile(member)):
                    damaged.append(member.name)
    return members, damaged


def listed_members(path, archive):
    """Size and members that the GEO file list records for one archive."""
    size, members, current = None, {}, None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split("\t")
        if fields[0] == "Archive":
            current = fields[1]
            size = int(fields[3]) if current == archive else size
        elif fields[0] == "File" and current == archive:
            members[fields[1]] = int(fields[3])
    return size, members


def study_problems(folder, part, tally):
    """Return the files of one study that differ from their record, do not open, or disagree with the GEO file list."""
    changed, damaged, unlisted = [], [], []
    recorded = dict(zip(part["file"], part["sha256"]))
    if stored_files(folder) != set(recorded) or checksums(folder) != recorded:
        return [f"{folder.name}: folder contents differ from the recorded files"], damaged, unlisted
    for row in part.itertuples(index=False):
        path, name = folder / row.file, row.file.lower()
        label = f"{folder.name}/{row.file}"
        tally["files"] += 1
        tally["bytes"] += int(row.bytes)
        if path.stat().st_size != row.bytes or C.sha(path) != row.sha256:
            changed.append(label)
            continue
        if name.endswith(".gz"):
            tally["compressed_files"] += 1
            with path.open("rb") as stream:
                if not gzip_intact(stream):
                    damaged.append(label)
        if name.endswith((".zip", ".xlsx")):
            try:
                with zipfile.ZipFile(path) as archive:
                    if archive.testzip() is not None:
                        damaged.append(label)
            except (zipfile.BadZipFile, OSError):
                damaged.append(label)
        if name.endswith((".tar", ".tar.gz", ".tgz")):
            try:
                members, broken = archive_members(path)
            except (tarfile.TarError, EOFError, OSError):
                damaged.append(label)
                continue
            tally["archives"] += 1
            tally["archive_members"] += len(members)
            damaged += [f"{label}:{member}" for member in broken]
            # GEO lists the members of an archive in the file list of the same folder.
            listing = (PurePosixPath(row.file).parent / "filelist.txt").as_posix()
            if listing in recorded and (row.bytes, members) != listed_members(folder / listing, path.name):
                unlisted.append(label)
    return changed, damaged, unlisted


def download_checks(qc, state):
    report, policy, included = C.read_json(C.receipt_path(4)), C.CONFIG["download"], state["included"]
    modes = {number: C.read_json(C.receipt_path(number)).get("search_mode") for number in (1, 2, 4)}
    qc.check("search_mode_recorded_by_every_asking_script", set(modes.values()) <= {"stored", "updated"}
             and (report["planned_from"] == "download lock") == (modes[4] == "stored"), modes)
    files = pd.read_parquet(C.DATA / "download_files.parquet")
    fetched, skipped = files[files["status"] == "downloaded"], files[files["status"] == "skipped_by_policy"]
    unwanted = lambda name: any(re.search(pattern, name, re.IGNORECASE) for pattern in policy["skip_patterns"])
    qc.check("every_planned_file_downloaded_or_skipped_by_policy", len(fetched) + len(skipped) == len(files)
             and all(map(unwanted, skipped["file"])) and not any(map(unwanted, fetched["file"])))
    qc.check("files_unique_within_a_study", not files.duplicated(["source", "accession", "file"]).any())
    # A series matrix carries values only for arrays; for a sequencing series it describes the samples.
    sequenced = set(included.loc[~included["assay"].str.split("+").map(lambda labels: "array" in labels), "accession"])
    valued = ~fetched["kind"].isin(["archive_listing", "series_matrix"]) | ((fetched["kind"] == "series_matrix") & ~fetched["accession"].isin(sequenced))
    qc.check("every_included_study_has_expression_files",
             set(fetched.loc[valued, "accession"]) == set(included["accession"]) == set(files["accession"]))
    locked = pd.read_csv(D.LOCK, dtype=str, keep_default_na=False)
    qc.check("download_lock_matches_the_downloaded_files", locked.equals(files.astype(str)))
    qc.check("download_receipt_reconciles", report["files"] == len(fetched) and report["bytes"] == int(fetched["bytes"].sum())
             and report["studies"] == fetched["accession"].nunique() and report["skipped_by_policy"] == skipped["file"].tolist())
    folders = {(source, path.name) for source in sorted(set(included["source"])) for path in (D.RAW / source).iterdir()}
    qc.check("no_folder_outside_the_included_studies", folders == set(zip(included["source"], included["accession"])))
    tally = dict(files=0, bytes=0, compressed_files=0, archives=0, archive_members=0)
    changed, damaged, unlisted = [], [], []
    groups = fetched.groupby(["source", "accession"], sort=False)
    for index, ((source, accession), part) in enumerate(groups, 1):
        found = study_problems(D.RAW / source / accession, part, tally)
        changed, damaged, unlisted = changed + found[0], damaged + found[1], unlisted + found[2]
        if index % 10 == 0 or index == len(groups):
            print(f"[05] verified {index:,}/{len(groups):,} studies, {D.gigabytes(tally['bytes'])}", flush=True)
    state["summary"].update(downloaded=tally, files_differing_from_record=changed, files_not_readable=damaged, archives_differing_from_file_list=unlisted)
    qc.check("downloaded_files_match_recorded_size_and_sha256", not changed, changed)
    qc.check("compressed_files_and_archives_readable", not damaged, damaged)
    qc.check("archive_members_match_the_geo_file_list", not unlisted, unlisted)
    qc.check("every_recorded_file_verified", tally["files"] == len(fetched) and tally["bytes"] == report["bytes"])


def main():
    # A record that is missing must stop the run; it is never fetched again during quality control.
    C.fetch = offline
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/5 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/5 Search responses and census replay", lambda: census_checks(qc, state))
        qc.phase("3/5 Stored records and table replay", lambda: record_checks(qc, state))
        qc.phase("4/5 Screening decisions against the records", lambda: screening_checks(qc, state))
        qc.phase("5/5 Downloaded files, checksums and archives", lambda: download_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
        required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
        runtime=C.runtime(), code=C.code_binding(), exception=failure,
        dependencies=C.bind([C.receipt_path(i) for i in range(1, 5) if C.receipt_path(i).exists()]),
        eligibility="ACQUISITION_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[05 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
