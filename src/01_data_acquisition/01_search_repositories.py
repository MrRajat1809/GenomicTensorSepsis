"""Search public repositories for human sepsis transcriptome studies and build the deduplicated candidate census."""
import argparse
import hashlib
import importlib.metadata
import json
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

import requests

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CENSUS = ROOT / "data/raw/census"
DATA = ROOT / "data/processed/acquisition"
OUT = ROOT / "outputs/data_acquisition"
CONFIG = json.loads((HERE / "acquisition_config.json").read_text(encoding="utf-8"))
REGISTER = json.loads((HERE / CONFIG["register"]).read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
STAGES = {
    1: "01_search_repositories.py", 2: "02_harvest_study_metadata.py", 3: "03_screen_studies.py",
    4: "04_download_raw_data.py", 5: "05_acquisition_master_qc.py",
}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1/"
ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
OMICSDI = "https://www.omicsdi.org/ws/dataset/search"
COLUMNS = ["source", "accession", "title", "n_samples", "declared_assay", "platform", "release_date", "access", "organism",
           "links", "raw_hint", "matched_terms", "census_status", "duplicate_of", "also_in", "is_superseries",
           "url", "description"]
# A GEO series lists its dataset types whole, separated by semicolons; methylation, SNP or non-coding RNA arrays are not expression arrays.
GEO_ASSAYS = {"Expression profiling by array": "array", "Expression profiling by high throughput sequencing": "rna-seq"}
STORED_ANSWERS = ("unpack snapshots.tar.gz from the GigaDB package into data/raw/, or run with --update-search to ask the repositories again "
                  "(new studies then need screening decisions)")
# The cutoff is checked against the census date stored with the search, not through the digest of the search settings, so that a search stored
# before the cutoff was set keeps its digest.
CUTOFF = ("census_date", "release_dates")
SESSION = requests.Session()
SESSION.headers["User-Agent"] = CONFIG["http"]["user_agent"]
LAST_REQUEST = {}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_path(recorded):
    path = (ROOT / recorded).resolve()
    require(path.is_relative_to(ROOT), f"Path outside project: {recorded}")
    return path


def record(path):
    path = Path(path).resolve()
    return {"sha256": sha(path), "size_bytes": path.stat().st_size}


def bind(paths):
    paths = sorted(set(Path(path).resolve() for path in paths), key=lambda path: path.as_posix())
    return {path.relative_to(ROOT).as_posix(): record(path) for path in paths}


def verify_files(records):
    for name, expected in records.items():
        path = project_path(name)
        require(path.is_file() and record(path) == expected, f"Missing or changed artifact: {name}")


def code_binding(number=5):
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "acquisition_config.json", HERE / CONFIG["register"]])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("requests", "pandas", "pyarrow"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def receipt_path(number):
    return OUT / f"{number:02d}.json"


def verify_stage(number):
    result = read_json(receipt_path(number))
    require(result["version"] == VERSION and result["status"] == "COMPLETE", f"Stage {number} incomplete")
    require(result["code"] == code_binding(number), f"Acquisition code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
    verify_files(result["outputs"])
    verify_files(result["dependencies"])
    return result


def begin(number):
    validate_config()
    print(f"[{number:02d}] Verifying inputs and provenance v{VERSION}", flush=True)
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if number < 5 and (OUT / "qc.json").exists():
        write_json(OUT / "qc.json", dict(version=VERSION, status="STALE", reason=f"Stage {number:02d} is being rerun; rerun 05 after completion"))
    for previous in range(1, number):
        verify_stage(previous)
    dependencies = bind([receipt_path(previous) for previous in range(1, number)])
    return {"code": code_binding(number), "dependencies": dependencies}


def validate_config():
    require(VERSION == "1.0.0", "Unsupported acquisition version")
    search, http = CONFIG["search"], CONFIG["http"]
    terms = search["terms"]
    require(terms and len(set(terms)) == len(terms) and all(re.fullmatch(r"[a-z]+( [a-z]+)*", term) for term in terms),
            "Search terms must be unique lower-case words or phrases")
    require(all(http[key] > 0 for key in ("minimum_interval_seconds", "timeout_seconds", "retries", "retry_wait_seconds")), "Invalid HTTP policy")
    require(sorted(search["geo"]["dataset_types"]) == sorted(GEO_ASSAYS), "Every searched GEO dataset type needs an assay label, and only those")
    sources = ["geo", "arrayexpress", "ena"] + [name.lower() for name in search["catalogue"]["repositories"]]
    require(sorted(CONFIG["deduplication"]["source_priority"]) == sorted(sources), "Source priority must list every searched source exactly once")
    dates = search["release_dates"]
    require(re.fullmatch(r"\d{4}-\d{2}-\d{2}", search["census_date"]) and release_day(search["census_date"]) == search["census_date"]
            and dates["formats"] and set(dates["undated"]) <= set(sources) and all(str(basis).strip() for basis in dates["undated"].values()),
            "search.census_date must be a date (YYYY-MM-DD) that the release-date formats read, and every undated source a searched one with a reason")
    require(set(search["catalogue"]["repositories"].values()) <= {"controlled", "unknown"}, "Catalogued repositories must be controlled or unknown access")
    require(set(search["arrayexpress"]["study_types"].values()) <= {"array", "rna-seq"}, "Unknown assay label for an ArrayExpress study type")
    require(0 < CONFIG["harvest"]["ena_sample_batch"] <= 100 and "sample_accession" in CONFIG["harvest"]["ena_run_fields"], "Invalid harvest batch policy")
    recall = [accession for accessions in CONFIG["recall"]["studies"].values() for accession in accessions]
    require(set(CONFIG["recall"]["studies"]) <= {"geo", "arrayexpress"} and recall and len(set(recall)) == len(recall), "Invalid recall set")
    technology = CONFIG["technology"]
    require(set(technology["geo_sample_types"].values()) <= {"array", "rna-seq"} and technology["sequencing_library_strategies"]
            and len(set(technology["sequencing_library_strategies"])) == len(technology["sequencing_library_strategies"]),
            "GEO sample types must map to array or rna-seq, and the sequencing library strategies must be unique")
    screening = CONFIG["screening"]
    require(list(screening["reasons"])[-1] == "eligible" and {"controlled_access", "no_public_reads", "raw_reads_only"} <= set(screening["reasons"]),
            "Screening reasons must end with eligible and name the three availability reasons")
    require(CONFIG["register"] == "screening_register.json" and list(REGISTER) == sorted(REGISTER)
            and all({"decision", "reason", "basis"} <= set(entry) <= {"decision", "reason", "basis", "overlaps_with", "flag"} for entry in REGISTER.values()),
            "screening_register.json must hold one decision per accession, sorted by accession, with decision, reason, basis and optionally overlaps_with and flag")
    download = CONFIG["download"]
    require(download["minimum_free_gigabytes_after_download"] >= 0 and download["block_megabytes"] > 0 and isinstance(download["skip_patterns"], list),
            "Invalid download policy")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Acquisition code changed during execution")
    verify_files(start["dependencies"])
    result = dict(version=VERSION, status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=runtime(), outputs=bind(outputs), **start, **detail)
    write_json(receipt_path(number), result)
    return result


def save_table(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(path)


def fetch(url, params=None, post=False):
    """Space requests per host and retry transient failures; any other HTTP error stops the run."""
    policy, host, problem = CONFIG["http"], urlsplit(url).netloc, None
    for attempt in range(policy["retries"]):
        wait = LAST_REQUEST.get(host, 0) + policy["minimum_interval_seconds"] - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        LAST_REQUEST[host] = time.monotonic()
        try:
            if post:
                response = SESSION.post(url, data=params, timeout=policy["timeout_seconds"])
            else:
                response = SESSION.get(url, params=params, timeout=policy["timeout_seconds"])
            if response.status_code not in (429, 500, 502, 503, 504):
                response.raise_for_status()
                return response.content
            problem = f"HTTP {response.status_code}"
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as error:
            problem = type(error).__name__
        time.sleep(policy["retry_wait_seconds"] * 2 ** attempt)
    raise ValueError(f"No valid response after {policy['retries']} attempts ({problem}): {url}")


def snapshot(source, name, content):
    """Store a repository response unmodified and return it parsed; ENA answers an empty result with an empty body."""
    path = CENSUS / source / name
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(content)
    temporary.replace(path)
    return json.loads(content) if content.strip() else []


def read_snapshot(source, name):
    content = (CENSUS / source / name).read_bytes()
    return json.loads(content) if content.strip() else []


def search_digest():
    settings = {key: value for key, value in CONFIG["search"].items() if key not in CUTOFF}
    return hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()


def search_geo(terms, queries):
    policy = CONFIG["search"]["geo"]
    types = " OR ".join(f'"{name}"[DataSet Type]' for name in policy["dataset_types"])
    uids = set()
    for name, term in terms.items():
        query = (f'"{term}"[All Fields] AND "{policy["organism"]}"[Organism] AND gse[Entry Type] AND ({types}) '
                 f'AND ("1900/01/01"[PDAT] : "{CONFIG["search"]["census_date"].replace("-", "/")}"[PDAT])')
        parameters = dict(db="gds", term=query, retmax=10000, retmode="json", tool="GenomicTensorSepsis")
        result = snapshot("geo", f"esearch.{name}.json", fetch(EUTILS + "esearch.fcgi", parameters, post=True))["esearchresult"]
        require(int(result["count"]) == len(result["idlist"]), f"GEO returned a truncated id list for: {term}")
        uids.update(result["idlist"])
        queries.append(dict(source="geo", term=term, query=query, hits=len(result["idlist"])))
    uids = sorted(uids, key=int)
    for index in range(0, len(uids), policy["summary_batch"]):
        parameters = dict(db="gds", id=",".join(uids[index:index + policy["summary_batch"]]), retmode="json", tool="GenomicTensorSepsis")
        snapshot("geo", f"esummary.{index // policy['summary_batch']:03d}.json", fetch(EUTILS + "esummary.fcgi", parameters, post=True))


def search_arrayexpress(terms, queries):
    policy = CONFIG["search"]["arrayexpress"]
    accessions = set()
    for name, term in terms.items():
        page, hits = 1, []
        while True:
            parameters = dict(query=f'"{term}"', organism=policy["organism"], pageSize=policy["page_size"], page=page)
            result = snapshot("arrayexpress", f"search.{name}.{page:02d}.json", fetch(BIOSTUDIES + "arrayexpress/search", parameters))
            hits += [hit["accession"] for hit in result["hits"]]
            if not result["hits"] or page * policy["page_size"] >= result["totalHits"]:
                break
            page += 1
        require(len(set(hits)) == result["totalHits"], f"ArrayExpress returned an incomplete result list for: {term}")
        accessions.update(hits)
        queries.append(dict(source="arrayexpress", term=term, query=f'"{term}" organism={policy["organism"]}', hits=len(hits)))
    for accession in sorted(accessions):
        snapshot("arrayexpress", f"study.{accession}.json", fetch(BIOSTUDIES + f"studies/{accession}"))


def ena_query(term):
    clauses = [f'study_title="{term}"', f'study_description="{term}"']
    if " " not in term:
        clauses += [f'study_title="*{term}*"', f'study_description="*{term}*"']
    return " OR ".join(clauses)


def search_ena(terms, queries):
    policy = CONFIG["search"]["ena"]
    accessions = set()
    for name, term in terms.items():
        query = f'({ena_query(term)}) AND first_public<={CONFIG["search"]["census_date"]}'
        parameters = dict(result="study", query=query, format="json", limit=0,
                          fields="study_accession,secondary_study_accession,study_title,study_description,geo_accession,first_public")
        hits = snapshot("ena", f"studies.{name}.json", fetch(ENA, parameters, post=True))
        accessions.update(hit["study_accession"] for hit in hits)
        queries.append(dict(source="ena", term=term, query=query, hits=len(hits)))
    accessions = sorted(accessions)
    for index in range(0, len(accessions), policy["run_batch"]):
        studies = " OR ".join(f'study_accession="{accession}"' for accession in accessions[index:index + policy["run_batch"]])
        query = (f'({studies}) AND tax_eq({policy["tax_id"]}) AND library_source="{policy["library_source"]}" '
                 f'AND library_strategy="{policy["library_strategy"]}"')
        parameters = dict(result="read_run", query=query, format="json", limit=0,
                          fields="study_accession,secondary_study_accession,study_title,sample_accession,instrument_model")
        snapshot("ena", f"runs.{index // policy['run_batch']:03d}.json", fetch(ENA, parameters, post=True))


def search_catalogue(terms, queries):
    policy = CONFIG["search"]["catalogue"]
    for repository in policy["repositories"]:
        for name, term in terms.items():
            query = f'"{term}" AND repository:"{repository}"'
            start, hits = 0, 0
            while True:
                parameters = dict(query=query, start=start, size=policy["page_size"])
                result = snapshot("catalogue", f"{repository}.{name}.{start // policy['page_size']:02d}.json", fetch(OMICSDI, parameters))
                hits += len(result.get("datasets", []))
                start += policy["page_size"]
                if not result.get("datasets") or start >= result.get("count", 0):
                    break
            queries.append(dict(source=repository.lower(), term=term, query=query, hits=hits))


def run_search(terms):
    if CENSUS.exists():
        shutil.rmtree(CENSUS)
    queries = []
    for label, action in (("GEO", search_geo), ("ArrayExpress", search_arrayexpress), ("ENA", search_ena), ("catalogued repositories", search_catalogue)):
        print(f"[01] Searching {label}", flush=True)
        action(terms, queries)
    files = sorted(path for path in CENSUS.rglob("*") if path.is_file())
    sums = "".join(f"{sha(path)}  {path.relative_to(CENSUS).as_posix()}\n" for path in files)
    (CENSUS / "SHA256SUMS.txt").write_text(sums, encoding="utf-8")
    search = dict(version=VERSION, searched_at_utc=datetime.now(timezone.utc).isoformat(), search_sha256=search_digest(),
                  census_date=CONFIG["search"]["census_date"], queries=queries)
    write_json(CENSUS / "search.json", search)
    return search


def census_day(search):
    """The census date of a stored search; a search stored before the cutoff was recorded reaches as far as the day it ran."""
    return search.get("census_date") or search["searched_at_utc"][:10]


def saved_search():
    require((CENSUS / "search.json").is_file() and (CENSUS / "SHA256SUMS.txt").is_file(), f"No stored search under data/raw/census; {STORED_ANSWERS}")
    search = read_json(CENSUS / "search.json")
    require(search["version"] == VERSION and search["search_sha256"] == search_digest(),
            "The search settings of acquisition_config.json differ from those of the stored search; restore them, or run with --update-search to ask "
            "the repositories again (new studies then need screening decisions)")
    require(census_day(search) == CONFIG["search"]["census_date"],
            f"search.census_date ({CONFIG['search']['census_date']}) is not the census date of the stored search ({census_day(search)}); restore it, or run "
            "with --update-search to search up to the new date (new studies then need screening decisions)")
    return search


def verify_snapshots():
    lines = (CENSUS / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines()
    expected = {name: digest for digest, name in (line.split("  ", 1) for line in lines)}
    actual = {path.relative_to(CENSUS).as_posix() for path in CENSUS.rglob("*") if path.is_file()} - {"SHA256SUMS.txt", "search.json"}
    require(set(expected) == actual and all(sha(CENSUS / name) == digest for name, digest in expected.items()),
            f"The stored search answers under data/raw/census are missing or changed; {STORED_ANSWERS}")
    return len(expected)


def joined(values):
    return ";".join(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))


def release_day(text):
    """The release date of a record as YYYY-MM-DD, an empty string where it has none, or None where it cannot be read."""
    text = str(text).strip()
    for form in CONFIG["search"]["release_dates"]["formats"] if text else ():
        try:
            return datetime.strptime(text, form).date().isoformat()
        except ValueError:
            continue
    return None if text else ""


def apply_cutoff(rows):
    """Keep a candidate released after the census date in the census under its own status; stop at a date that cannot be read, and at an undated
    record of a source whose records should carry a date."""
    policy, unread, undated = CONFIG["search"], [], []
    for row in rows:
        if row["census_status"] != "candidate":
            continue
        day = release_day(row["release_date"])
        if day is None:
            unread.append(f"{row['accession']} ({row['release_date']})")
        elif not day and row["source"] not in policy["release_dates"]["undated"]:
            undated.append(f"{row['source']}:{row['accession']}")
        elif day > policy["census_date"]:
            row["census_status"] = "released_after_census_date"
    require(not unread, f"Release dates that no format of search.release_dates.formats reads: {unread[:20]}; add their format there")
    require(not undated, f"Candidates without a release date: {undated[:20]}; name their source in search.release_dates.undated with the reason its "
                         "records carry none, if they belong in the census")


def geo_records(terms):
    matched, rows = {}, []
    for name, term in terms.items():
        for uid in read_snapshot("geo", f"esearch.{name}.json")["esearchresult"]["idlist"]:
            matched.setdefault(uid, []).append(term)
    for path in sorted((CENSUS / "geo").glob("esummary.*.json")):
        result = read_json(path)["result"]
        for uid in result["uids"]:
            doc = result[uid]
            require(doc["entrytype"] == "GSE" and uid in matched, f"Unexpected GEO record: {uid}")
            types = {name.strip() for name in doc["gdstype"].split(";")}
            assay = "+".join(label for name, label in GEO_ASSAYS.items() if name in types)
            linked = [doc.get("bioproject", "")] + [item["targetobject"] for item in doc.get("extrelations", []) if item["relationtype"] == "SRA"]
            rows.append(dict(source="geo", accession=doc["accession"], title=doc["title"], n_samples=int(doc["n_samples"] or 0),
                declared_assay=assay, platform=joined("GPL" + item for item in doc["gpl"].split(";") if item),
                release_date=doc["pdat"].replace("/", "-"), access="open", organism=doc["taxon"], links=joined(linked),
                raw_hint=doc.get("suppfile", ""), matched_terms=joined(matched.pop(uid)), census_status="candidate",
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={doc['accession']}", description=doc["summary"]))
    require(not matched, f"GEO summaries missing for {len(matched)} matched records")
    return rows


def attributes(node):
    return {item["name"]: str(item.get("value", "")) for item in node.get("attributes", []) if "name" in item}


def sections(node):
    """Yield a BioStudies section and every nested subsection; subsections arrive as dicts or lists."""
    if isinstance(node, list):
        for item in node:
            yield from sections(item)
    elif isinstance(node, dict):
        yield node
        yield from sections(node.get("subsections", []))


def structured_links(node):
    """Yield (type, target) for links a submitter declared; accessions mentioned in free text are not identity."""
    if isinstance(node, list):
        for item in node:
            yield from structured_links(item)
    elif isinstance(node, dict) and "url" in node:
        yield attributes(node).get("Type", ""), node["url"]


def arrayexpress_records(terms):
    policy, matched, rows = CONFIG["search"]["arrayexpress"], {}, []
    for name, term in terms.items():
        for path in sorted((CENSUS / "arrayexpress").glob(f"search.{name}.*.json")):
            for hit in read_json(path)["hits"]:
                matched.setdefault(hit["accession"], []).append(term)
    for accession in sorted(matched):
        study = read_snapshot("arrayexpress", f"study.{accession}.json")
        root = study.get("section", {})
        main = attributes(root)
        parts = {part.get("type"): attributes(part) for part in sections(root.get("subsections", []))}
        linked = [pair for part in sections(root) for pair in structured_links(part.get("links", []))]
        imported = re.fullmatch(r"E-GEOD-(\d+)", accession)
        kind = main.get("Study type", "")
        score = parts.get("MIAME Score") or parts.get("MINSEQE Score") or {}
        rows.append(dict(source="arrayexpress", accession=accession, title=main.get("Title", ""),
            n_samples=int(parts.get("Samples", {}).get("Sample count") or 0), declared_assay=policy["study_types"].get(kind, ""),
            platform=joined(target for label, target in linked if label == "Array Design"),
            release_date=attributes(study).get("ReleaseDate", ""), access="open",
            organism=joined(item.get("value", "") for item in root.get("attributes", []) if item.get("name") == "Organism"),
            links=joined(([f"GSE{imported.group(1)}"] if imported else []) + [target for label, target in linked if label == "ENA"]),
            raw_hint="raw" if score.get("Raw") == "*" else "", matched_terms=joined(matched[accession]),
            census_status="candidate" if kind in policy["study_types"] else "study_type_not_expression",
            url=f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}", description=main.get("Description", "")))
    return rows


def ena_records(terms):
    studies, matched, runs, rows = {}, {}, {}, []
    for name, term in terms.items():
        for hit in read_snapshot("ena", f"studies.{name}.json"):
            studies[hit["study_accession"]] = hit
            matched.setdefault(hit["study_accession"], []).append(term)
    for path in sorted((CENSUS / "ena").glob("runs.*.json")):
        for run in read_snapshot("ena", path.name):
            runs.setdefault(run["study_accession"], []).append(run)
    # A matched umbrella project reports its runs under child studies, which appear here without matched terms.
    for accession in sorted(set(studies) | set(runs)):
        study, found = studies.get(accession, {}), runs.get(accession, [])
        first = found[0] if found else {}
        secondary = study.get("secondary_study_accession") or first.get("secondary_study_accession", "")
        rows.append(dict(source="ena", accession=accession, title=study.get("study_title") or first.get("study_title", ""),
            n_samples=len({run["sample_accession"] for run in found}), declared_assay="rna-seq" if found else "",
            platform=joined(sorted({run["instrument_model"] for run in found})), release_date=study.get("first_public", ""),
            access="open", organism="Homo sapiens" if found else "", links=joined(secondary.split(";") + [study.get("geo_accession", "")]),
            raw_hint="FASTQ" if found else "", matched_terms=joined(matched.get(accession, [])),
            census_status="candidate" if found else "no_human_rnaseq_runs",
            url=f"https://www.ebi.ac.uk/ena/browser/view/{accession}", description=study.get("study_description", "")))
    return rows


def catalogue_records(terms):
    found, matched = {}, {}
    for repository, access in CONFIG["search"]["catalogue"]["repositories"].items():
        for name, term in terms.items():
            for path in sorted((CENSUS / "catalogue").glob(f"{repository}.{name}.*.json")):
                for hit in read_json(path).get("datasets", []):
                    found[hit["id"]] = (repository, access, hit)
                    matched.setdefault(hit["id"], []).append(term)
    return [dict(source=repository.lower(), accession=accession, title=hit.get("title") or "", n_samples=0, declared_assay="",
                 platform="", release_date=hit.get("publicationDate") or "", access=access,
                 organism=joined(item.get("name", "") for item in hit.get("organisms") or []), links="", raw_hint="",
                 matched_terms=joined(matched[accession]), census_status="candidate", url="", description=hit.get("description") or "")
            for accession, (repository, access, hit) in sorted(found.items())]


def deduplicate(rows):
    """Mark a candidate as a duplicate when a higher-priority source already holds one of its accessions."""
    owner, also = {}, {}
    for source in CONFIG["deduplication"]["source_priority"]:
        claimed = {}
        for row in rows:
            if row["source"] != source or row["census_status"] != "candidate":
                continue
            aliases = [row["accession"], *filter(None, row["links"].split(";"))]
            match = next((owner[alias] for alias in aliases if alias in owner), None)
            if match:
                row.update(census_status="cross_source_duplicate", duplicate_of=match)
                also.setdefault(match, []).append(row["accession"])
            for alias in aliases:
                claimed.setdefault(alias, match or row["accession"])
        # Applied once per source so that two records of the same source are never merged.
        for alias, accession in claimed.items():
            owner.setdefault(alias, accession)
    for row in rows:
        row["also_in"] = joined(also.get(row["accession"], []))
    return owner


def counts(frame, column):
    return {str(key): int(value) for key, value in frame[column].value_counts().sort_index().items()}


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-search", action="store_true",
                        help="Ask the repositories again, up to search.census_date, instead of reading the stored answers under data/raw/census; "
                             "new studies then need screening decisions")
    args = parser.parse_args()
    print(f"[01] Building the repository census v{VERSION}", flush=True)
    start = begin(1)
    terms = {term.replace(" ", "_"): term for term in CONFIG["search"]["terms"]}
    search = run_search(terms) if args.update_search else saved_search()
    snapshots = verify_snapshots()
    rows = geo_records(terms) + arrayexpress_records(terms) + ena_records(terms) + catalogue_records(terms)
    for row in rows:
        row.update(duplicate_of="", is_superseries=CONFIG["search"]["geo"]["superseries_marker"] in row["description"])
    apply_cutoff(rows)
    owner = deduplicate(rows)
    priority = CONFIG["deduplication"]["source_priority"]
    rows.sort(key=lambda row: (priority.index(row["source"]), row["accession"]))
    census = pd.DataFrame(rows, columns=COLUMNS)
    require(not census.duplicated(["source", "accession"]).any(), "Repeated census record within a source")
    candidates = census[census["census_status"] == "candidate"].reset_index(drop=True)
    require(len(candidates) and candidates["accession"].is_unique and (candidates["title"].str.len() > 0).all(),
            "Candidates must be non-empty, uniquely identified and titled")
    save_table(census, DATA / "census.parquet")
    save_table(candidates, DATA / "candidates.parquet")
    attrition = []
    for source in priority:
        part = census[census["source"] == source]
        attrition.append(dict(source=source, records=len(part), **counts(part, "census_status")))
        print(f"[01] {source}: {len(part):,} records; {int((part['census_status'] == 'candidate').sum()):,} candidates", flush=True)
    recall = [dict(source=source, accession=accession, found=accession in owner, candidate=owner.get(accession, ""))
              for source, accessions in CONFIG["recall"]["studies"].items() for accession in accessions]
    missed = [item["accession"] for item in recall if not item["found"]]
    require(not missed, f"Recall studies not reached by the search: {missed}; revise the search terms or the recall set in acquisition_config.json")
    unscreened = sorted(set(candidates["accession"]) - set(REGISTER))
    if unscreened:
        print(f"[01] {len(unscreened)} candidates without a screening decision, at which script 03 stops: {unscreened[:20]}", flush=True)
    finish(1, start, [DATA / "census.parquet", DATA / "candidates.parquet", CENSUS / "search.json", CENSUS / "SHA256SUMS.txt"],
           searched_at_utc=search["searched_at_utc"], search_sha256=search["search_sha256"], snapshot_files=snapshots,
           search_mode="updated" if args.update_search else "stored", census_date=census_day(search),
           released_after_census_date=int((census["census_status"] == "released_after_census_date").sum()),
           undated_candidates=counts(candidates[candidates["release_date"] == ""], "source"), queries=search["queries"], attrition=attrition,
           candidates=dict(total=len(candidates), by_source=counts(candidates, "source"), by_declared_assay=counts(candidates, "declared_assay"),
                           by_access=counts(candidates, "access"), superseries_flagged=int(candidates["is_superseries"].sum()),
                           ena_umbrella_children=int(((candidates["source"] == "ena") & (candidates["matched_terms"] == "")).sum())),
           recall=dict(basis=CONFIG["recall"]["basis"], checked=len(recall), found=len(recall) - len(missed), studies=recall),
           protocol=CONFIG["protocol"])
    print(f"[01] Complete: {len(candidates):,} candidates from {len(census):,} records; recall {len(recall) - len(missed)}/{len(recall)}", flush=True)


if __name__ == "__main__":
    main()
