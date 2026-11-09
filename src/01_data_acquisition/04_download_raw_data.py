"""Plan the expression-data downloads of every included study and, on request, fetch them with checksums."""
import argparse
import importlib
import re
import shutil
import time
from datetime import datetime, timezone
from urllib.parse import quote, unquote

import urllib3

C = importlib.import_module("01_search_repositories")

RAW = C.ROOT / "data/raw"
GEO_SERIES = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
NCBI_COUNTS = "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={0}&format=file&file={0}_raw_counts_GRCh38.p13_NCBI.tsv.gz"
BIOSTUDIES_FILES = "https://www.ebi.ac.uk/biostudies/files/"
COLUMNS = ["source", "accession", "file", "kind", "url", "bytes", "status", "sha256"]
LOCK = RAW / "download_lock.csv"
ARRAYEXPRESS_KINDS = {"Raw Data": "raw", "Processed Data": "processed"}
# Servers compress text files in transit, so sizes and resumed ranges only agree for the stored bytes.
STORED = {"Accept-Encoding": "identity"}


def remote_size(url):
    """Size the server announces for a file, or None when the file does not exist."""
    policy, problem = C.CONFIG["http"], None
    for attempt in range(policy["retries"]):
        time.sleep(policy["minimum_interval_seconds"])
        try:
            response = C.SESSION.head(url, headers=STORED, allow_redirects=True, timeout=policy["timeout_seconds"])
            if response.status_code == 404:
                return None
            if response.status_code == 200 and response.headers.get("Content-Length"):
                return int(response.headers["Content-Length"])
            problem = f"HTTP {response.status_code}"
        except (C.requests.ConnectionError, C.requests.Timeout) as error:
            problem = type(error).__name__
        time.sleep(policy["retry_wait_seconds"] * 2 ** attempt)
    raise ValueError(f"No file size after {policy['retries']} attempts ({problem}): {url}")


def folder_names(url):
    """Names of the files in a GEO folder; a series without supplementary files has no folder."""
    try:
        index = C.fetch(url).decode("utf-8", errors="replace")
    except C.requests.HTTPError as error:
        if error.response.status_code != 404:
            raise
        return []
    return sorted(set(re.findall(r'href="([^"/?]+)"', index)))


def geo_files(accession, assay, supplementary, subseries):
    folder = f"{GEO_SERIES}{accession[:-3]}nnn/{accession}/"
    listed = [(name.rsplit("/", 1)[1], "supplementary") for name in supplementary.split(";") if name.startswith("ftp://")]
    if any(name.endswith("_RAW.tar") for name, _ in listed):
        listed.append(("filelist.txt", "archive_listing"))
    rows = [dict(file=name, kind=kind, url=f"{folder}suppl/{quote(name)}") for name, kind in listed]
    # A SuperSeries holds the arrays of its SubSeries, but the tables its authors deposited stay with each SubSeries.
    for member in filter(None, subseries.split(";")):
        base = f"{GEO_SERIES}{member[:-3]}nnn/{member}/suppl/"
        rows += [dict(file=f"{member}/{unquote(name)}", kind="archive_listing" if name == "filelist.txt" else "subseries_supplementary", url=base + name)
                 for name in folder_names(base)]
    index = C.fetch(folder + "matrix/").decode("utf-8", errors="replace")
    rows += [dict(file=name, kind="series_matrix", url=f"{folder}matrix/{quote(name)}")
             for name in sorted(set(re.findall(r'href="([^"/?]+_series_matrix\.txt\.gz)"', index)))]
    for row in rows:
        row["bytes"] = remote_size(row["url"])
        C.require(row["bytes"] is not None, f"{accession}: listed file is missing at the repository: {row['file']}")
    if "rna-seq" in assay:
        url = NCBI_COUNTS.format(accession)
        size = remote_size(url)
        if size is not None:
            rows.append(dict(file=f"{accession}_raw_counts_GRCh38.p13_NCBI.tsv.gz", kind="ncbi_generated_counts", url=url, bytes=size))
    return rows


def declared_files(node, section=""):
    """Yield (section type, file entry) for every file a BioStudies record declares."""
    if isinstance(node, list):
        for item in node:
            yield from declared_files(item, section)
    elif isinstance(node, dict):
        section = node.get("type", section) if "path" not in node else section
        if "path" in node:
            yield section, node
        for value in node.values():
            if isinstance(value, (list, dict)):
                yield from declared_files(value, section)


def arrayexpress_files(accession):
    study = C.read_snapshot("arrayexpress", f"study.{accession}.json")
    return [dict(file=entry["path"], kind=ARRAYEXPRESS_KINDS[section], bytes=int(entry["size"]),
                 url=f"{BIOSTUDIES_FILES}{accession}/{quote(entry['path'])}")
            for section, entry in declared_files(study) if section in ARRAYEXPRESS_KINDS]


def validators(response):
    """The entity tag and the modification time by which a server names the version of the file it sends."""
    return {name: response.headers.get(name, "") for name in ("ETag", "Last-Modified")}


def download(url, path, expected):
    """Fetch one file, resuming a partial copy only from the same version of the file; it is accepted only at exactly the announced size."""
    policy, partial = C.CONFIG["http"], path.with_name(path.name + ".part")
    version = partial.with_name(partial.name + ".json")
    path.parent.mkdir(parents=True, exist_ok=True)
    stalled, problem = 0, None
    furthest = partial.stat().st_size if partial.exists() else 0
    while stalled < policy["retries"]:
        held = partial.stat().st_size if partial.exists() else 0
        known = C.read_json(version) if version.is_file() else {}
        # The server sends the rest of a file only while it still has the version named here; otherwise it sends the whole file again.
        token = known.get("ETag") if known.get("ETag") and not known["ETag"].startswith("W/") else known.get("Last-Modified")
        if held and not (held < expected and known.get("url") == url and known.get("bytes") == expected and token):
            partial.unlink()
            held = 0
        try:
            headers = dict(STORED, **{"Range": f"bytes={held}-", "If-Range": token}) if held else STORED
            with C.SESSION.get(url, headers=headers, stream=True, timeout=policy["timeout_seconds"]) as response:
                problem = f"HTTP {response.status_code}"
                coding = response.headers.get("Content-Encoding", "").strip().lower()
                same = all(value == known[name] for name, value in validators(response).items() if value and known.get(name))
                resumed = response.status_code == 206 and held > 0 and same and response.headers.get("Content-Range") == f"bytes {held}-{expected - 1}/{expected}"
                if coding not in ("", "identity"):
                    problem = f"content coded in transit ({coding})"
                elif response.status_code == 200 and response.headers.get("Content-Length", str(expected)) != str(expected):
                    problem = f"the server sends {response.headers['Content-Length']} bytes where {expected:,} were announced"
                elif response.status_code == 200 or resumed:
                    if not resumed:
                        C.write_json(version, dict(url=url, bytes=expected, **validators(response)))
                    # The stored bytes as they arrive, never decoded, so that their size and any resumed range refer to the file itself.
                    with partial.open("ab" if resumed else "wb") as stream:
                        for block in response.raw.stream(C.CONFIG["download"]["block_megabytes"] * 1024 * 1024, decode_content=False):
                            stream.write(block)
                elif response.status_code == 206:
                    problem = "HTTP 206 with another range or another version of the file"
                    partial.unlink(missing_ok=True)
        except (C.requests.ConnectionError, C.requests.Timeout, C.requests.exceptions.ChunkedEncodingError, urllib3.exceptions.HTTPError) as error:
            problem = type(error).__name__
        size = partial.stat().st_size if partial.exists() else 0
        if size == expected:
            partial.replace(path)
            version.unlink(missing_ok=True)
            return
        if size > expected:
            problem = f"{size:,} bytes received where {expected:,} were announced"
            partial.unlink()
        # Only an attempt that holds more of the file than every attempt before counts as progress, so that a large file survives any number of
        # dropped connections while a file that keeps arriving too long, or that must start again each time, stops after the retries.
        if furthest < size < expected:
            furthest, stalled = size, 0
        else:
            stalled += 1
        time.sleep(policy["retry_wait_seconds"] * 2 ** stalled)
    if not partial.exists():
        version.unlink(missing_ok=True)
    raise ValueError(f"Download incomplete after {policy['retries']} attempts without progress ({problem}): {url}")


def gigabytes(value):
    return f"{value / 1e9:,.1f} GB"


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Fetch the planned files; without it only the plan is written")
    parser.add_argument("--update-search", action="store_true",
                        help="List the files of the included studies at the repositories again and rewrite data/raw/download_lock.csv; without it the "
                             "plan is the stored lock, and a file is accepted only at its recorded SHA-256")
    args = parser.parse_args()
    print(f"[04] Planning expression-data downloads v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["download"]
    included = pd.read_parquet(C.DATA / "included_studies.parquet")
    studies = pd.read_parquet(C.DATA / "studies.parquet").set_index("accession")
    rows, locked, dropped = [], None, []
    if not args.update_search:
        # The files, sizes and checksums of an earlier run, so that a rebuild fetches exactly what that run fetched and asks no repository for a listing.
        C.require(LOCK.is_file(), f"No download lock at {LOCK.relative_to(C.ROOT).as_posix()}; {C.STORED_ANSWERS}")
        locked = pd.read_csv(LOCK, dtype={"source": str, "accession": str, "file": str, "kind": str, "url": str, "bytes": "int64", "status": str, "sha256": str},
                             keep_default_na=False)
        unlisted = sorted(set(included["accession"]) - set(locked["accession"]))
        C.require(list(locked.columns) == COLUMNS and not unlisted,
                  f"The download lock lists no files of the included studies {unlisted[:20]}; run with --update-search to list them at the repositories")
        # A locked study that screening no longer includes leaves the plan, and the lock keeps only the rows of the included studies; quality control
        # then stops until the stored folder of that study is removed.
        dropped = sorted(set(locked["accession"]) - set(included["accession"]))
        locked = locked[locked["accession"].isin(included["accession"])].reset_index(drop=True)
        rows = locked.assign(status=locked["status"].replace({"downloaded": "planned"})).to_dict("records")
        print(f"[04] planned {len(included):,} studies from the download lock; locked studies no longer included: {dropped or 'none'}", flush=True)
    for index, study in enumerate(included.itertuples(index=False), 1):
        if locked is not None:
            break
        if study.source == "geo":
            found = geo_files(study.accession, study.assay, studies.loc[study.accession, "supplementary_files"], study.superseries_of)
        else:
            found = arrayexpress_files(study.accession)
        for row in found:
            skipped = any(re.search(pattern, row["file"], re.IGNORECASE) for pattern in policy["skip_patterns"])
            rows.append(dict(source=study.source, accession=study.accession, status="skipped_by_policy" if skipped else "planned", sha256="", **row))
        if index % 20 == 0 or index == len(included):
            print(f"[04] planned {index:,}/{len(included):,} studies", flush=True)
    plan = pd.DataFrame(rows, columns=COLUMNS)
    expected = dict(zip(zip(plan["accession"], plan["file"]), plan["sha256"])) if locked is not None else {}
    plan["sha256"] = ""
    C.require(not plan.duplicated(["source", "accession", "file"]).any() and (plan["bytes"] > 0).all(), "Plan holds a repeated or empty file")
    wanted = plan[plan["status"] == "planned"]
    # The member list of an archive holds no values, and neither does the series matrix of a sequencing series, which only describes its samples.
    arrays = {study.accession for study in included.itertuples(index=False) if "array" in study.assay.split("+")}
    holding = wanted[(wanted["kind"] != "archive_listing") & ((wanted["kind"] != "series_matrix") | wanted["accession"].isin(arrays))]
    without = sorted(set(included["accession"]) - set(holding["accession"]))
    report = C.OUT / "download_plan.csv"
    temporary = report.with_name(report.name + ".tmp")
    plan.drop(columns="sha256").to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(report)
    RAW.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(RAW).free
    held = sum(row.bytes for row in wanted.itertuples(index=False)
               if (RAW / row.source / row.accession / row.file).is_file() and (RAW / row.source / row.accession / row.file).stat().st_size == row.bytes)
    print(f"[04] Plan: {len(wanted):,} files in {wanted['accession'].nunique()} studies, {gigabytes(wanted['bytes'].sum())}; "
          f"already on disk {gigabytes(held)}; free space {gigabytes(free)}", flush=True)
    for kind, part in wanted.groupby("kind"):
        print(f"[04]   {kind}: {len(part):,} files, {gigabytes(part['bytes'].sum())}", flush=True)
    if without:
        print(f"[04] Studies with no expression file to download: {without}", flush=True)
    if not args.download:
        print(f"[04] Plan written to {report.relative_to(C.ROOT).as_posix()}; rerun with --download to fetch the files", flush=True)
        return
    C.require(not without, f"Included studies without any expression file; revise their screening decision: {without}")
    C.require(free - (wanted["bytes"].sum() - held) >= policy["minimum_free_gigabytes_after_download"] * 1e9,
              "Not enough free disk space for the planned downloads")
    done = 0
    for accession, part in wanted.groupby("accession", sort=False):
        folder = RAW / part["source"].iloc[0] / accession
        for row in part.sort_values("bytes").itertuples(index=False):
            path = folder / row.file
            if not (path.is_file() and path.stat().st_size == row.bytes):
                download(row.url, path, row.bytes)
        sums = {name: C.sha(folder / name) for name in sorted(part["file"])}
        changed = sorted(name for name, digest in sums.items() if expected and expected[(accession, name)] != digest)
        C.require(not changed, f"{accession}: files differ from the download lock: {changed}; the repository changed them since the locked run")
        (folder / "SHA256SUMS.txt").write_text("".join(f"{digest}  {name}\n" for name, digest in sums.items()), encoding="utf-8")
        plan.loc[(plan["accession"] == accession) & (plan["status"] == "planned"), "sha256"] = plan["file"].map(sums)
        plan.loc[(plan["accession"] == accession) & (plan["status"] == "planned"), "status"] = "downloaded"
        done += 1
        print(f"[04] {done:,}/{wanted['accession'].nunique():,} {accession}: {len(part)} files, {gigabytes(part['bytes'].sum())}", flush=True)
    fetched = plan[plan["status"] == "downloaded"]
    C.require(len(fetched) == len(wanted) and fetched["sha256"].str.len().eq(64).all(), "Downloaded files do not reconcile with the plan")
    C.save_table(plan, C.DATA / "download_files.parquet")
    temporary = LOCK.with_name(LOCK.name + ".tmp")
    plan.to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(LOCK)
    outputs = [C.DATA / "download_files.parquet", report, LOCK] + [RAW / source / accession / "SHA256SUMS.txt"
               for source, accession in fetched[["source", "accession"]].drop_duplicates().itertuples(index=False)]
    C.finish(4, start, outputs, downloaded_at_utc=datetime.now(timezone.utc).isoformat(),
             files=len(fetched), bytes=int(fetched["bytes"].sum()), studies=int(fetched["accession"].nunique()),
             by_kind={kind: dict(files=len(part), bytes=int(part["bytes"].sum())) for kind, part in fetched.groupby("kind")},
             skipped_by_policy=plan.loc[plan["status"] == "skipped_by_policy", "file"].tolist(),
             ncbi_generated_counts=int((fetched["kind"] == "ncbi_generated_counts").sum()), search_mode="updated" if args.update_search else "stored",
             planned_from="download lock" if locked is not None else "repositories",
             lock_studies_not_included=dropped, protocol=C.CONFIG["protocol"])
    print(f"[04] Complete: {len(fetched):,} files, {gigabytes(fetched['bytes'].sum())}, in {fetched['accession'].nunique()} studies", flush=True)


if __name__ == "__main__":
    main()
