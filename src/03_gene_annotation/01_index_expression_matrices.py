"""Index every expression matrix and every included sample, so that any value can be followed back to its record and its source file."""
import hashlib
import importlib.metadata
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

import requests

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RAW = ROOT / "data/raw/annotation"
# A run reads the stored answers of the repositories; it asks them again only with --update-search.
STORED_ANSWERS = "unpack snapshots.tar.gz from the GigaDB package into data/raw/, or run with --update-search to ask the repositories again"
ACQUIRED = ROOT / "data/processed/acquisition"
EXPRESSION = ROOT / "data/processed/expression"
EXPRESSION_OUT = ROOT / "outputs/expression_processing"
DATA = ROOT / "data/processed/annotation"
OUT = ROOT / "outputs/gene_annotation"
WORK = DATA / "work"
CONFIG = json.loads((HERE / "annotation_config.json").read_text(encoding="utf-8"))
REGISTER = json.loads((HERE / CONFIG["register"]).read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
REGISTER_SECTIONS = ("sources", "identifier_replace", "ncbi_check_exceptions")
STAGES = {1: "01_index_expression_matrices.py", 2: "02_collect_annotation_sources.py", 3: "03_map_features_to_genes.py",
          4: "04_build_gene_matrices.py", 5: "05_annotation_master_qc.py"}
KEYWORDS = ("per_sample_files", "ncbi_counts", "series_matrix", "none")
SESSION = requests.Session()
SESSION.headers["User-Agent"] = CONFIG["http"]["user_agent"]
LAST_REQUEST = {}
R_VERSIONS = r"""
suppressPackageStartupMessages(library(jsonlite))
found <- lapply(fromJSON(commandArgs(TRUE)[1]), function(name) if (nzchar(system.file(package = name))) as.character(packageVersion(name)) else NA)
names(found) <- fromJSON(commandArgs(TRUE)[1])
cat(toJSON(c(list(R = paste(R.version$major, R.version$minor, sep = ".")), found), auto_unbox = TRUE, na = "null"))
"""
# For each processing script of the expression stage: its table of study platforms, its table of samples, the column that counts
# the rows of a matrix, and what it did.
PROCESSED = {
    "04_process_affymetrix_arrays.py": ("affymetrix_studies", "affymetrix_arrays", "features", "RMA from CEL files"),
    "05_process_single_channel_arrays.py": ("single_channel_studies", "single_channel_arrays", "probes",
                                            "limma normexp background correction, quantile normalization and log2, from scanner files"),
    "06_process_intensity_tables.py": ("intensity_table_studies", "intensity_table_arrays", "probes", ""),
    "07_process_count_tables.py": ("count_table_studies", "count_table_samples", "genes",
                                   "counts as deposited; edgeR TMM-normalized log2 counts per million"),
    "08_collect_author_matrices.py": ("author_matrix_studies", "author_matrix_samples", "features", ""),
}
CORRECTIONS = {"detection_p_values": "limma neqc with negative controls inferred from detection p-values, from deposited intensity tables",
               "intensity_distribution": "limma normexp background correction (RMA-75 estimate), quantile normalization and log2, from deposited intensity tables",
               "none": "quantile normalization and log2 of deposited processed signals"}
MATRIX_COLUMNS = ["source", "accession", "platform", "platform_title", "family", "tier", "data_source", "script", "method", "value_scale", "samples",
                  "in_matrix", "features", "identifier_shape", "identifier_shape_share", "identifier_example", "matrix", "counts", "record", "reason"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "family", "title", "source_name", "attribute_values", "tier",
                  "data_source", "script", "in_matrix", "reason", "matrix", "counts", "values_from", "values_label", "array_id", "scan_time", "normalization_group"]


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


def code_binding(number=max(STAGES)):
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "annotation_config.json", HERE / CONFIG["register"]])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("requests", "numpy", "pandas", "pyarrow"):
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
    require(result["code"] == code_binding(number), f"Annotation code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
    verify_files(result["outputs"])
    verify_files(result["dependencies"])
    return result


def begin(number):
    validate_config()
    print(f"[{number:02d}] Verifying inputs and provenance v{VERSION}", flush=True)
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if number < max(STAGES) and (OUT / "qc.json").exists():
        write_json(OUT / "qc.json", dict(version=VERSION, status="STALE", reason=f"Stage {number:02d} is being rerun; rerun the master QC after completion"))
    for previous in range(1, number):
        verify_stage(previous)
    if number > 1:
        verify_files(read_json(receipt_path(1))["external"])
    dependencies = bind([receipt_path(previous) for previous in range(1, number)])
    return {"code": code_binding(number), "dependencies": dependencies}


def validate_config():
    require(VERSION == "1.0.0", "Unsupported annotation version")
    http = CONFIG["http"]
    require(all(http[key] > 0 for key in ("minimum_interval_seconds", "timeout_seconds", "retries", "retry_wait_seconds")), "Invalid HTTP policy")
    shapes = CONFIG["index"]["identifier_shapes"]
    require(shapes and all(entry.get("shape") and re.compile(entry["pattern"]) for entry in shapes)
            and len({entry["shape"] for entry in shapes}) == len(shapes), "Invalid identifier shapes")
    sources = CONFIG["sources"]
    require(sources["hgnc"]["url"] and sources["hgnc"]["file"] and {"hgnc_id", "symbol", "status", "entrez_id", "ensembl_gene_id"} <= set(sources["hgnc"]["columns"])
            and sources["geo_platform"]["url"] and sources["package_extras"] is not None, "Invalid annotation sources")
    require(CONFIG["register"] == "annotation_register.json" and tuple(REGISTER) == REGISTER_SECTIONS
            and all(list(REGISTER[name]) == sorted(REGISTER[name]) for name in REGISTER_SECTIONS)
            and not any("|" in key for key in REGISTER["sources"])
            and all(key.count("|") == 1 for name in REGISTER_SECTIONS[1:] for key in REGISTER[name]),
            "annotation_register.json must hold its sections in their order, sorted by key: the source of every platform by its accession, and the "
            "decisions of every study platform by accession|platform")
    for name, entry in REGISTER["sources"].items():
        packaged = entry.get("source") == "bioconductor" and entry.get("package")
        tabled = (entry.get("source") == "geo_platform_table" and entry.get("identifier") and entry.get("genes")
                  and set(entry["genes"]) <= {"entrez", "ensembl", "symbol", "accession"})
        require((packaged or tabled) and str(entry.get("basis", "")).strip(),
                f"{name}: a platform needs an annotation package, or the identifier and gene columns of its platform table, and a basis")
    mapping = CONFIG["mapping"]
    require(0 < mapping["minimum_identifier_share"] <= 1 and all(re.compile(mapping[key]) for key in ("value_separator", "entrez_gene", "ensembl_gene"))
            and re.compile(mapping["ensembl_gene"]).groups >= 1
            and all(str(entry.get("basis", "")).strip() and entry.get("replace") and all(len(pair) == 2 and re.compile(pair[0]) for pair in entry["replace"])
                    for entry in REGISTER["identifier_replace"].values()), "Invalid mapping policy")
    matrices = CONFIG["gene_matrices"]
    require(matrices["rules"] and set(matrices["rules"].values()) <= {"highest_mean", "sum", "sum_counts"} and matrices["counts_tolerance"] > 0,
            "Invalid gene-matrix policy")
    check = matrices["ncbi_check"]
    require(isinstance(check["minimum_samples"], int) and check["minimum_samples"] >= 2 and check["minimum_shared_genes"] > 0 and check["genes_used"] > 0
            and str(check.get("basis", "")).strip()
            and all(entry.get("status") in ("too_few_samples", "too_few_shared_genes") and str(entry.get("basis", "")).strip()
                    for entry in REGISTER["ncbi_check_exceptions"].values()),
            "Invalid NCBI check policy: an exception needs the status it excuses and a basis")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Annotation code changed during execution")
    verify_files(start["dependencies"])
    result = dict(version=VERSION, status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=runtime(), outputs=bind(outputs), **start, **detail)
    write_json(receipt_path(number), result)
    return result


def fetch(url, params=None):
    """Space requests per host and retry transient failures; any other HTTP error stops the run."""
    policy, host, problem = CONFIG["http"], urlsplit(url).netloc, None
    for attempt in range(policy["retries"]):
        wait = LAST_REQUEST.get(host, 0) + policy["minimum_interval_seconds"] - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        LAST_REQUEST[host] = time.monotonic()
        try:
            response = SESSION.get(url, params=params, timeout=policy["timeout_seconds"])
            if response.status_code not in (429, 500, 502, 503, 504):
                response.raise_for_status()
                return response.content
            problem = f"HTTP {response.status_code}"
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as error:
            problem = type(error).__name__
        time.sleep(policy["retry_wait_seconds"] * 2 ** attempt)
    raise ValueError(f"No valid response after {policy['retries']} attempts ({problem}): {url}")


def retrieval_note(path):
    return path.with_name(path.name + ".retrieved.json")


def retrieve(name, url, params=None, update=False):
    """Read one stored annotation source under data/raw/annotation, or ask its source again and store it unmodified beside a note of when it was
    retrieved; return its path and that note."""
    path = RAW / name
    note = retrieval_note(path)
    require(update or path.is_file(), f"No stored annotation source data/raw/annotation/{name}; {STORED_ANSWERS}")
    if update:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(fetch(url, params))
        temporary.replace(path)
        write_json(note, dict(file=name, url=url, parameters=params or {}, retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
                              recorded_from="download", **record(path)))
    elif not note.is_file():
        # A stored copy without a note is dated once by its modification time, and keeps that date however it is copied later.
        write_json(note, dict(file=name, url=url, parameters=params or {}, retrieved_at_utc=datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
                              recorded_from="file modification time", **record(path)))
    retrieved = read_json(note)
    require(retrieved["file"] == name and {key: retrieved[key] for key in ("sha256", "size_bytes")} == record(path),
            f"{name}: the stored source is not the file its retrieval note describes; {STORED_ANSWERS}")
    return path, retrieved


def tool_versions(packages):
    """Versions of R and of the named R packages; a package that is not installed has none."""
    WORK.mkdir(parents=True, exist_ok=True)
    script = WORK / "versions.R"
    script.write_text(R_VERSIONS, encoding="utf-8")
    done = subprocess.run(["Rscript", "--vanilla", script.as_posix(), json.dumps(sorted(packages))], capture_output=True, text=True)
    script.unlink()
    require(done.returncode == 0, f"R could not report its package versions:\n{done.stderr[-2000:]}")
    return json.loads(done.stdout)


def clear_work(folder):
    """Remove the temporary folder of a script, and the shared work folder once it is empty."""
    shutil.rmtree(folder, ignore_errors=True)
    if WORK.is_dir() and not any(WORK.iterdir()):
        WORK.rmdir()


def save_table(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(path)


def save_review(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(path)


def expression_binding():
    """Confirm that the expression stage passed its quality control unchanged; return the artifacts this stage binds and the checksum of every expression output."""
    qc = read_json(EXPRESSION_OUT / "qc.json")
    require(qc["status"] == "PASS" and qc["failures"] == 0, "Expression quality control has not passed")
    verify_files(qc["code"])
    verify_files(qc["dependencies"])
    recorded = {}
    for name in qc["dependencies"]:
        recorded.update(read_json(project_path(name))["outputs"])
    # The tables are verified here; a matrix is large and is verified when it is read.
    tables = {name: expected for name, expected in recorded.items() if "/study_matrices/" not in name and "/study_counts/" not in name}
    verify_files(tables)
    return [EXPRESSION_OUT / "qc.json", *map(project_path, qc["dependencies"]), *map(project_path, tables)], recorded, qc


def expression_matrix(name, recorded):
    """Return the path of a file of the expression stage after checking it against the checksum in its stage receipt."""
    path = project_path(name)
    require(name in recorded and path.is_file() and record(path) == recorded[name], f"Missing or changed expression matrix: {name}")
    return path


def script_of(tier, family):
    """Name the expression script that processes a route, or nothing for a study platform that is not read."""
    if tier == "raw_array_files":
        return "04_process_affymetrix_arrays.py" if family == "affymetrix" else "05_process_single_channel_arrays.py"
    return {"unnormalized_intensities": "06_process_intensity_tables.py", "counts": "07_process_count_tables.py",
            "author_normalized": "08_collect_author_matrices.py"}.get(tier, "")


def identifier_shape(identifiers, shapes):
    """Name the shape that most identifiers of a matrix have, its share and one identifier of that shape."""
    import collections

    named = [next(name for name, pattern in shapes if pattern.fullmatch(value)) for value in identifiers]
    shape, count = collections.Counter(named).most_common(1)[0]
    return shape, round(count / len(named), 4), identifiers[named.index(shape)]


def main():
    import pandas as pd
    import pyarrow.parquet as pq

    print(f"[01] Indexing expression matrices v{VERSION}", flush=True)
    start = begin(1)
    external, recorded, qc = expression_binding()
    shapes = [(entry["shape"], re.compile(entry["pattern"])) for entry in CONFIG["index"]["identifier_shapes"]]
    routes = pd.read_parquet(EXPRESSION / "processing_routes.parquet")
    held = pd.read_parquet(EXPRESSION / "sample_routes.parquet")
    assigned = pd.read_parquet(EXPRESSION / "table_samples.parquet")
    records = pd.read_parquet(ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", "title", "source_name"])
    attributes = pd.read_parquet(ACQUIRED / "sample_attributes.parquet", columns=["study_accession", "sample_accession"])
    external += [ACQUIRED / "samples.parquet", ACQUIRED / "sample_attributes.parquet"]
    routes["script"] = [script_of(tier, family) for tier, family in zip(routes["tier"], routes["family"])]
    studies = {script: pd.read_parquet(EXPRESSION / f"{names[0]}.parquet").set_index(["accession", "platform"]) for script, names in PROCESSED.items()}
    listed = {script: pd.read_parquet(EXPRESSION / f"{names[1]}.parquet") for script, names in PROCESSED.items()}

    matrices = []
    for index, route in enumerate(routes.itertuples(index=False), 1):
        row = dict(source=route.source, accession=route.accession, platform=route.platform, platform_title=route.platform_title, family=route.family,
                   tier=route.tier, data_source=route.data_source if route.data_source in KEYWORDS else "table", script=route.script, method="",
                   value_scale="", samples=route.samples, in_matrix=0, features=0, identifier_shape="", identifier_shape_share=None,
                   identifier_example="", matrix="", counts="", record="", reason="")
        study = studies[route.script].loc[(route.accession, route.platform)] if route.script else None
        if study is None or not str(study.get("matrix") or ""):
            # A study platform without a matrix is kept with the reason its stage recorded.
            row.update(reason="not_processed" if study is None else f"study_platform_{study['reason']}")
            matrices.append(row)
            continue
        counted, method = PROCESSED[route.script][2:]
        if route.script == "06_process_intensity_tables.py":
            method = CORRECTIONS[study["correction"]]
        elif route.script == "08_collect_author_matrices.py":
            method = f"as deposited by the authors: {study['values']}"
        scale = {"07_process_count_tables.py": "log2 counts per million", "08_collect_author_matrices.py": str(study.get("scale"))}.get(route.script, "log2")
        identifiers = pq.read_table(expression_matrix(study["matrix"], recorded), columns=["feature_id"]).column("feature_id").to_pylist()
        require(len(identifiers) == int(study[counted]), f"{route.accession}|{route.platform}: the matrix does not hold the rows its stage recorded")
        shape, share, example = identifier_shape(identifiers, shapes)
        counts = str(study["counts"]) if route.script == "07_process_count_tables.py" else ""
        if counts:
            expression_matrix(counts, recorded)
        row.update(method=method, value_scale=scale, in_matrix=int(study["in_matrix"]) if "in_matrix" in study.index else int(study["samples"]),
                   features=len(identifiers), identifier_shape=shape, identifier_shape_share=share, identifier_example=example, matrix=study["matrix"],
                   counts=counts, record=str(Path(study["matrix"]).with_suffix(".json").as_posix()))
        matrices.append(row)
        if index % 20 == 0 or index == len(routes):
            print(f"[01] indexed {index:,}/{len(routes):,} study platforms", flush=True)
    matrices = pd.DataFrame(matrices, columns=MATRIX_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)

    # Where the values of each sample were read from, as the processing script of its study platform recorded it.
    origin = {}
    for script, table in listed.items():
        for item in table.itertuples(index=False):
            inside = bool(getattr(item, "in_matrix", True))
            source = getattr(item, "table", None) or "/".join(part for part in (getattr(item, "file", ""), getattr(item, "member", "")) if part)
            origin[(item.study_accession, item.sample_accession)] = dict(
                in_matrix=inside, values_from=source, values_label=getattr(item, "label", "") or "", array_id=getattr(item, "array_id", "") or "",
                scan_time=getattr(item, "scan_time", "") or "", normalization_group=getattr(item, "normalization_group", "") or "",
                reason="" if inside else "incomplete_array_file")
    files = matrices.set_index(["accession", "platform"])
    in_table = assigned.set_index(["study_accession", "sample_accession"])["reason"]
    described = records.set_index(["study_accession", "sample_accession"])
    counted = attributes.groupby(["study_accession", "sample_accession"]).size()
    rows = []
    for item in held.merge(routes[["accession", "platform", "script"]], left_on=["study_accession", "platform"], right_on=["accession", "platform"],
                           validate="many_to_one").itertuples(index=False):
        key, study = (item.study_accession, item.sample_accession), files.loc[(item.study_accession, item.platform)]
        found = origin.get(key, dict(in_matrix=False, values_from="", values_label="", array_id="", scan_time="", normalization_group="", reason=""))
        reason = found["reason"]
        if not found["in_matrix"] and not reason:
            # A sample without a column carries the reason of its study platform, of the label map or of the route, in that order.
            reason = study["reason"] if study["reason"].startswith("study_platform") else in_table.get(key, "") or item.reason
        rows.append(dict(source=item.source, study_accession=item.study_accession, sample_accession=item.sample_accession, platform=item.platform,
                         family=item.family, title=described.loc[key, "title"], source_name=described.loc[key, "source_name"],
                         attribute_values=int(counted.get(key, 0)), tier=item.tier, data_source=study["data_source"], script=item.script,
                         in_matrix=found["in_matrix"], reason=reason, matrix=study["matrix"] if found["in_matrix"] else "",
                         counts=study["counts"] if found["in_matrix"] else "", values_from=found["values_from"], values_label=found["values_label"],
                         array_id=found["array_id"], scan_time=found["scan_time"], normalization_group=found["normalization_group"]))
    samples = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    reasons = {str(reason): int(count) for reason, count in samples.loc[~samples["in_matrix"], "reason"].value_counts().sort_index().items()}
    require(len(samples) == len(held) == qc["summary"]["samples"] and not samples.duplicated(["study_accession", "sample_accession"]).any()
            and (samples["in_matrix"] == (samples["reason"] == "")).all(), "The sample index does not hold every included sample once with one outcome")
    require(int(samples["in_matrix"].sum()) == qc["summary"]["in_a_matrix"] == int(matrices["in_matrix"].sum()) and reasons == qc["summary"]["without_a_matrix"]
            and int((matrices["matrix"] != "").sum()) == qc["summary"]["matrices"], "The indexes do not reconcile with the quality control of the expression stage")

    save_table(matrices, DATA / "matrix_index.parquet")
    save_table(samples, DATA / "sample_index.parquet")
    save_review(matrices, OUT / "matrix_index_review.csv")
    by_shape = {str(shape): dict(matrices=len(part), samples=int(part["in_matrix"].sum()))
                for shape, part in matrices[matrices["matrix"] != ""].groupby("identifier_shape")}
    finish(1, start, [DATA / "matrix_index.parquet", DATA / "sample_index.parquet", OUT / "matrix_index_review.csv"], external=bind(external),
           study_platforms=dict(total=len(matrices), with_a_matrix=int((matrices["matrix"] != "").sum())),
           samples=dict(total=len(samples), in_a_matrix=int(samples["in_matrix"].sum()), without_a_matrix=reasons,
                        with_a_source_file_or_column=int((samples["values_from"] != "").sum()), with_an_array_identifier=int((samples["array_id"] != "").sum()),
                        with_a_scan_time=int((samples["scan_time"] != "").sum()), attribute_values=int(samples["attribute_values"].sum())),
           by_identifier_shape=by_shape, identifier_shapes=CONFIG["index"]["identifier_shapes"], protocol=CONFIG["protocol"])
    print(f"[01] Complete: {len(matrices)} study platforms, {int((matrices['matrix'] != '').sum())} with a matrix; {len(samples):,} samples, "
          f"{int(samples['in_matrix'].sum()):,} in a matrix; review sheet: {(OUT / 'matrix_index_review.csv').relative_to(ROOT).as_posix()}", flush=True)
    for shape, item in sorted(by_shape.items(), key=lambda pair: -pair[1]["samples"]):
        print(f"[01]   {shape}: {item['matrices']} matrices, {item['samples']:,} samples", flush=True)
    print(f"[01]   without a matrix: {reasons}", flush=True)


if __name__ == "__main__":
    main()
