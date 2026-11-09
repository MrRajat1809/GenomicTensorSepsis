"""Retrieve the platform records of the included samples and inventory the expression files each study provides."""
import argparse
import collections
import csv
import gzip
import hashlib
import importlib.metadata
import io
import json
import re
import shutil
import struct
import subprocess
import sys
import tarfile
import time
import warnings
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

import requests

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RAW = ROOT / "data/raw"
PLATFORMS = RAW / "platforms"
ACQUIRED = ROOT / "data/processed/acquisition"
ACQUISITION = ROOT / "outputs/data_acquisition"
DATA = ROOT / "data/processed/expression"
OUT = ROOT / "outputs/expression_processing"
MATRICES = DATA / "study_matrices"
COUNTS = DATA / "study_counts"
WORK = DATA / "work"
CONFIG = json.loads((HERE / "expression_config.json").read_text(encoding="utf-8"))
REGISTER = json.loads((HERE / CONFIG["register"]).read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
REGISTER_SECTIONS = ("routes", "label_rewrites", "labels_by_values", "reassigned_files", "incomplete_files", "count_tables", "normalization_groups",
                     "author_matrices", "not_collected", "value_check_exceptions")
STAGES = {1: "01_inventory_expression_files.py", 2: "02_assign_processing_routes.py", 3: "03_map_table_labels.py",
          4: "04_process_affymetrix_arrays.py", 5: "05_process_single_channel_arrays.py", 6: "06_process_intensity_tables.py",
          7: "07_process_count_tables.py", 8: "08_collect_author_matrices.py", 9: "09_expression_master_qc.py"}
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
ARRAY_DESIGN = "https://www.ebi.ac.uk/biostudies/files/{0}/{0}.adf.txt"
# A run reads the stored answers of the repositories; it asks them again only with --update-search.
STORED_ANSWERS = "unpack snapshots.tar.gz from the GigaDB package into data/raw/, or run with --update-search to ask the repositories again"
ARCHIVES = (".tar", ".tar.gz", ".tgz")
SAMPLE = re.compile(r"(GSM\d+)(?!\d)")
NUMBER = re.compile(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?")
SUFFIX = re.compile(CONFIG["inventory"]["label_suffix"], re.IGNORECASE)
SEPARATORS = {"tab": "\t", "comma": ",", "semicolon": ";", "whitespace": r"\s+"}
NOT_RECORDED = ["", "NA", "N/A", "#N/A", "NaN", "nan", "NULL", "null"]
PLATFORM_COLUMNS = ["platform", "family", "title", "technology", "manufacturer", "distribution", "organism", "studies", "samples"]
FILE_COLUMNS = ["source", "accession", "file", "member", "sheet", "kind", "bytes", "format", "compressed", "sample_accession", "delimiter",
                "columns", "labels", "samples_matched", "value_kind", "value_minimum", "value_maximum", "declared_rows", "first_labels"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "family", "own_files", "own_formats", "author_tables",
                  "ncbi_counts", "series_matrix_rows", "data_sources", "author_table_files"]
STUDY_COLUMNS = ["source", "accession", "platform", "family", "platform_title", "samples", "with_own_file", "own_formats", "with_author_table",
                 "with_ncbi_counts", "with_series_matrix_values", "unmatched", "tables"]
BLANK = dict(sheet="", delimiter="", columns=0, labels=0, matched={}, value_kind="none", value_minimum=None, value_maximum=None,
             declared_rows=0, first_labels="")
SESSION = requests.Session()
SESSION.headers["User-Agent"] = CONFIG["http"]["user_agent"]
LAST_REQUEST = {}
R_VERSIONS = r"""
suppressPackageStartupMessages(library(jsonlite))
found <- lapply(fromJSON(commandArgs(TRUE)[1]), function(name) if (nzchar(system.file(package = name))) as.character(packageVersion(name)) else NA)
names(found) <- fromJSON(commandArgs(TRUE)[1])
cat(toJSON(c(list(R = paste(R.version$major, R.version$minor, sep = ".")), found), auto_unbox = TRUE, na = "null"))
"""
csv.field_size_limit(2 ** 31 - 1)


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
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "expression_config.json", HERE / CONFIG["register"]])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("requests", "numpy", "pandas", "pyarrow", "openpyxl", "xlrd"):
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
    require(result["code"] == code_binding(number), f"Expression code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
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
    require(VERSION == "1.0.0", "Unsupported expression version")
    http, platforms, inventory = CONFIG["http"], CONFIG["platforms"], CONFIG["inventory"]
    require(all(http[key] > 0 for key in ("minimum_interval_seconds", "timeout_seconds", "retries", "retry_wait_seconds")), "Invalid HTTP policy")
    require(platforms["families"] and all(rule.get("family") and len(rule) > 1 and set(rule) <= {"family", "technology", "keyword"}
                                           for rule in platforms["families"]), "A platform family rule needs a family and a technology or keyword")
    require(platforms["unmatched_family"] and set(platforms["unrecorded"]) == {"library_strategy", "family"}, "Invalid platform fallbacks")
    assigned = {rule["family"] for rule in platforms["families"]} | {platforms["unmatched_family"], platforms["unrecorded"]["family"]}
    require(sorted(family for families in platforms["assays"].values() for family in families) == sorted(assigned),
            "platforms.assays must name every platform family once, under the study technology of the acquisition stage that it measures")
    require(CONFIG["register"] == "routes_register.json" and tuple(REGISTER) == REGISTER_SECTIONS
            and all(list(REGISTER[name]) == sorted(REGISTER[name]) and all(key.count("|") == 1 for key in REGISTER[name]) for name in REGISTER_SECTIONS),
            "routes_register.json must hold its sections in their order, each with one entry per study platform or table, keyed by accession|platform "
            "or accession|table and sorted by key")
    require(all({"tier", "data_source", "basis"} <= set(entry) <= {"tier", "data_source", "basis", "flag", "raw_not_usable"} for entry in REGISTER["routes"].values()),
            "A route needs a tier, a data source and a basis, and may carry besides only a flag and raw_not_usable")
    require(0 < inventory["member_preview_kilobytes"] <= inventory["file_preview_kilobytes"] <= inventory["maximum_preview_megabytes"] * 1024
            and inventory["preview_rows"] > 0 and inventory["minimum_identifier_characters"] > 0, "Invalid preview policy")
    require(inventory["per_sample_data_formats"] and not set(inventory["per_sample_data_formats"]) & set(inventory["annotation_formats"]),
            "A format cannot be both sample data and platform annotation")
    routes = CONFIG["routes"]
    require({"raw_array_files", "counts"} <= set(routes["tiers"])
            and all(tier["families"] and tier["data_sources"] for tier in routes["tiers"].values()), "Invalid route policy")
    labels = CONFIG["labels"]
    require(labels["identifier_parts"] > 0 and labels["values"]["minimum_shared_rows"] > 0
            and 0 < labels["values"]["minimum_margin"] < labels["values"]["minimum_correlation"] <= 1
            and 0 < labels["values"]["minimum_present_share"] < 1 and labels["values"]["minimum_columns"] >= 3, "Invalid label policy")
    require(all(str(entry.get("basis", "")).strip() for entry in [*REGISTER["label_rewrites"].values(), *REGISTER["labels_by_values"].values()])
            and all(entry.get("replace") and all(len(pair) == 2 and re.compile(pair[0]) for pair in entry["replace"])
                    for entry in REGISTER["label_rewrites"].values()), "A label decision needs a basis, and a rewrite its patterns and replacements")
    affymetrix = CONFIG["affymetrix"]
    require(all(affymetrix[key] > 0 for key in ("workers", "block_megabytes", "minimum_free_gigabytes")) and affymetrix["designs"]
            and all(design.get("package") and set(design) <= {"package", "target"} for design in affymetrix["designs"].values()),
            "Invalid Affymetrix policy")
    single = CONFIG["single_channel"]
    require(single["families"] and single["offset"] >= 0 and single["layouts"]
            and all({"feature", "probe", "foreground", "background", "kind", "measured"} <= set(layout) and layout["feature"] and layout["measured"]
                    for layout in single["layouts"].values())
            and all(entry.get("samples") and str(entry.get("basis", "")).strip() for entry in REGISTER["incomplete_files"].values()),
            "Invalid single-channel policy")
    intensity = CONFIG["intensity_tables"]
    require(intensity["offset"] >= 0 and 0 < intensity["detected_below"] < 1 and intensity["same_median_intensity"] > 0 and intensity["columns"]
            and intensity["estimator_without_detection"] in ("rma75", "saddle", "mle")
            and all(rule.get("signal") and isinstance(rule.get("detection"), list) and rule.get("correction") in ("normexp", "none")
                    for rule in intensity["columns"].values()), "Invalid intensity-table policy")
    counts = CONFIG["count_tables"]
    allowed = {"basis", "identifier_column", "decimal", "identifier_replace", "genes", "repeated_columns", "repeated_identifiers"}
    require(counts["prior_count"] > 0 and counts["count_measures"] and re.compile(counts["summary_rows"])
            and counts["header_identifiers"] and all(plain(word) for word in counts["header_identifiers"])
            and all(str(entry.get("basis", "")).strip() and len(entry) > 1 and set(entry) <= allowed
                    and entry.get("genes", "shared") in ("shared", "absent_is_zero")
                    and entry.get("repeated_columns", "sum") == entry.get("repeated_identifiers", "sum") == "sum"
                    and entry.get("decimal", ",") == "," and all(len(pair) == 2 and re.compile(pair[0]) for pair in entry.get("identifier_replace", []))
                    for entry in REGISTER["count_tables"].values()), "Invalid count-table policy")
    require(CONFIG["author_matrices"]["largest_log2_like_value"] > 0 and not set(REGISTER["author_matrices"]) & set(REGISTER["not_collected"])
            and all(entry.get("values") and entry.get("scale") in ("log2_like", "linear", "other") and str(entry.get("basis", "")).strip()
                    for entry in REGISTER["author_matrices"].values())
            and all(str(entry.get("basis", "")).strip() for entry in REGISTER["not_collected"].values()), "Invalid author-matrix policy")
    groups = REGISTER["normalization_groups"]
    require(all(str(entry.get("basis", "")).strip() and ("field" in entry) != ("attribute" in entry) and entry.get("patterns")
                and all(len(pair) == 2 and re.compile(pair[0]) and pair[1] for pair in entry["patterns"]) for entry in groups.values()),
            "A normalization-group decision needs a basis, one field or attribute, and its patterns with their groups")
    quality, reference, matrix = CONFIG["quality_control"], CONFIG["quality_control"]["rma_reference"], CONFIG["quality_control"]["matrix_check"]
    require(0 < quality["value_check"]["minimum_shared_rows"] <= quality["value_check"]["rows_used"] and quality["value_check"]["minimum_samples"] >= 3
            and reference["tolerance"] > 0 and reference["studies"] and len(set(reference["studies"])) == len(reference["studies"])
            and str(reference.get("basis", "")).strip()
            and all(matrix[key] > 0 for key in ("value_range_relative", "value_range_absolute", "library_size_relative", "log_cpm_absolute"))
            and str(matrix.get("basis", "")).strip()
            and quality["tmm_report_log2"] > 0 and all(str(entry.get("basis", "")).strip() and set(entry) <= {"basis", "samples"}
                                                       for entry in REGISTER["value_check_exceptions"].values())
            and all(str(entry.get("basis", "")).strip() and entry.get("files") for entry in REGISTER["reassigned_files"].values()),
            "Invalid quality-control policy or file reassignment")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Expression code changed during execution")
    verify_files(start["dependencies"])
    result = dict(version=VERSION, status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=runtime(), outputs=bind(outputs), **start, **detail)
    write_json(receipt_path(number), result)
    return result


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode("utf-8")).hexdigest()


def tool_versions(packages):
    """Versions of R and of the named R packages; a package that is not installed has none."""
    WORK.mkdir(parents=True, exist_ok=True)
    script = WORK / "versions.R"
    script.write_text(R_VERSIONS, encoding="utf-8")
    done = subprocess.run(["Rscript", "--vanilla", script.as_posix(), json.dumps(sorted(packages))], capture_output=True, text=True)
    script.unlink()
    require(done.returncode == 0, f"R could not report its package versions:\n{done.stderr[-2000:]}")
    return json.loads(done.stdout)


def kept_record(record_path, matrix_path, identity):
    """Return the stored record of a study platform when what it was computed from is the same and its matrix is unchanged, otherwise None."""
    if not (record_path.is_file() and matrix_path.is_file()):
        return None
    stored_record = read_json(record_path)
    return stored_record if stored_record.get("identity") == identity and stored_record.get("matrix") == record(matrix_path) else None


def normalization_groups(key, samples):
    """The normalization group of every named sample of a study platform by its recorded decision, or None where the study platform is one group."""
    import pandas as pd

    entry = REGISTER["normalization_groups"].get(key)
    if entry is None:
        return None
    study = key.split("|")[0]
    if "field" in entry:
        records = pd.read_parquet(ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", entry["field"]])
        records = records[records["study_accession"] == study]
        values = dict(zip(records["sample_accession"], records[entry["field"]].fillna("").astype(str)))
    else:
        attributes = pd.read_parquet(ACQUIRED / "sample_attributes.parquet", columns=["study_accession", "sample_accession", "key", "value"])
        attributes = attributes[(attributes["study_accession"] == study) & (attributes["key"] == entry["attribute"])]
        values = attributes.groupby("sample_accession")["value"].agg(lambda part: " + ".join(dict.fromkeys(map(str, part)))).to_dict()
    groups = {}
    for sample in samples:
        found = [group for pattern, group in entry["patterns"] if re.search(pattern, values.get(sample, ""))]
        require(len(found) == 1, f"{key}: sample {sample} matches {len(found)} patterns of its normalization groups, not exactly one")
        groups[sample] = found[0]
    require(len(set(groups.values())) > 1, f"{key}: all of its samples fall in one normalization group; remove its decision")
    return groups


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


def acquisition_binding():
    """Confirm that the acquisition stage passed its quality control unchanged and return the artifacts this stage reads."""
    qc = read_json(ACQUISITION / "qc.json")
    require(qc["status"] == "PASS" and qc["failures"] == 0, "Acquisition quality control has not passed")
    verify_files(qc["code"])
    verify_files(qc["dependencies"])
    paths = [ACQUISITION / "qc.json"]
    for name in qc["dependencies"]:
        outputs = read_json(project_path(name))["outputs"]
        verify_files(outputs)
        paths += [project_path(output) for output in outputs]
    return paths


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


def retrieve(name, url, params=None, update=False):
    """Read one stored platform record under data/raw/platforms, or ask the repository again and store its answer unmodified."""
    path = PLATFORMS / name
    if update:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(fetch(url, params))
        temporary.replace(path)
    require(path.is_file(), f"No stored platform record data/raw/platforms/{name}; {STORED_ANSWERS}")
    return path


def geo_platform(path, accession):
    fields = collections.defaultdict(list)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    require(lines and lines[0].strip() == f"^PLATFORM = {accession}", f"{accession}: GEO returned no platform record")
    for line in lines[1:]:
        key, _, value = line[1:].partition(" = ")
        fields[key.strip()].append(value.strip())
    first = lambda key: fields[key][0] if fields[key] else ""
    return dict(title=first("Platform_title"), technology=first("Platform_technology"), manufacturer=first("Platform_manufacturer"),
                distribution=first("Platform_distribution"), organism=";".join(dict.fromkeys(fields["Platform_organism"])))


def array_design(path, accession):
    """Read the header of an ArrayExpress array design file, which ends where the probe table begins."""
    fields = {}
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            cells = [cell.strip() for cell in line.rstrip("\r\n").split("\t")]
            if cells[0].startswith("["):
                break
            if len(cells) > 1 and cells[0]:
                fields.setdefault(cells[0], cells[1])
    require(fields.get("Comment[ArrayExpressAccession]") == accession and fields.get("Array Design Name"),
            f"{accession}: ArrayExpress returned no array design")
    return dict(title=fields["Array Design Name"], technology=fields.get("Technology Type", ""), manufacturer="", distribution="",
                organism=fields.get("Comment[Organism]", ""))


def platform_family(found):
    text = f"{found['manufacturer']} {found['title']}".lower()
    for rule in CONFIG["platforms"]["families"]:
        if rule.get("technology", found["technology"]) == found["technology"] and rule.get("keyword", "") in text:
            return rule["family"]
    return CONFIG["platforms"]["unmatched_family"]


def harvest_platforms(names, update):
    """Read one record per platform: GEO for GPL accessions, the ArrayExpress array design for the others."""
    listing = PLATFORMS / "SHA256SUMS.txt"
    require(update or listing.is_file(), f"No checksum list of the stored platform records at data/raw/platforms/SHA256SUMS.txt; {STORED_ANSWERS}")
    stored = "" if update else listing.read_text(encoding="utf-8")
    rows = []
    for name in names:
        if name.startswith("GPL"):
            found = geo_platform(retrieve(f"{name}.soft", GEO, dict(acc=name, targ="self", form="text", view="brief"), update), name)
        else:
            found = array_design(retrieve(f"{name}.adf.txt", ARRAY_DESIGN.format(name), update=update), name)
        rows.append(dict(platform=name, family=platform_family(found), **found))
    files = sorted(path for path in PLATFORMS.iterdir() if path.name != "SHA256SUMS.txt")
    require(len(files) == len(names), "Platform records do not match the platforms of the included samples; run with --update-search")
    sums = "".join(f"{sha(path)}  {path.name}\n" for path in files)
    require(update or sums == stored, f"The stored platform records under data/raw/platforms differ from their checksum list; {STORED_ANSWERS}")
    listing.write_text(sums, encoding="utf-8")
    return rows


def signature(head):
    """Name a file format from its first bytes; None means text that still has to be read as a table."""
    if head[:4] == b"IDAT":
        return "illumina_idat"
    if head[:5] == b"[CEL]" or head[:8] == struct.pack("<ii", 64, 4):
        return "affymetrix_cel"
    if head[:2] == b"\x3b\x01":
        # Command Console files share one container; the data type inside names what the file holds.
        return "affymetrix_cel" if b"affymetrix-calvin-intensity" in head else "affymetrix_chp" if b"-analysis" in head else "affymetrix_other"
    if head[:4] == struct.pack("<i", 65):
        return "affymetrix_chp"
    if head[:4] == struct.pack("<i", 67) or head[:5] == b"[CDF]":
        return "affymetrix_cdf"
    if head[:2] == b"#%":
        return "affymetrix_library_file"
    if b"[Heading]" in head[:200] and b"Illumina" in head[:200]:
        return "illumina_manifest"
    if head[:5] == b"TYPE\t" and b"FEPARAMS" in head:
        return "agilent_feature_extraction"
    if head[:3] == b"ATF":
        return "genepix_results"
    if head[:12] == b"Begin Header":
        return "imagene_results"
    if head[:14] == b"%%MatrixMarket":
        return "matrix_market"
    if head[:8] == b"!Series_":
        return "geo_series_matrix"
    if head[:4] == b"PK\x03\x04":
        return "zip"
    if head[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1":
        return "ole"
    if head[:5] == b"%PDF-":
        return "pdf"
    if head[257:262] == b"ustar":
        return "tar_archive"
    if b"\x00" in head and head[:2] not in (b"\xff\xfe", b"\xfe\xff"):
        return "unrecognised_binary"
    return None


def stored(opener):
    """Open a file as its stored bytes, looking through gzip compression."""
    stream = opener()
    packed = stream.read(2) == b"\x1f\x8b"
    stream.seek(0)
    return (gzip.GzipFile(fileobj=stream) if packed else stream), packed


def text_rows(content, complete):
    """Split a text preview into rows of cells on the delimiter that most lines share; a cut last line is dropped."""
    lines = content.decode("utf-16" if content[:2] in (b"\xff\xfe", b"\xfe\xff") else "utf-8-sig", errors="replace").splitlines()
    lines = [line for line in (lines if complete else lines[:-1]) if line.strip()]
    delimiter = next((mark for mark in ("\t", ",", ";") if sum(mark in line for line in lines) * 2 > len(lines)), None)
    if delimiter == "\t":
        return [[cell.strip().strip('"') for cell in line.split("\t")] for line in lines], "tab"
    if delimiter:
        return [[cell.strip() for cell in row] for row in csv.reader(lines, delimiter=delimiter)], "comma" if delimiter == "," else "semicolon"
    return [line.split() for line in lines], "whitespace"


def table_header(rows):
    """Return the width of a previewed table, its header labels and the row at which its body starts."""
    widths = collections.Counter(len(row) for row in rows)
    width = max(widths, key=lambda size: (widths[size], size)) if widths else 0
    if width < 2:
        return width, [], 0
    first = next(index for index, row in enumerate(rows) if len(row) == width)
    if first and len(rows[first - 1]) == width - 1:
        # R writes row names without a corner label, so the header is one cell shorter than the rows.
        return width, [""] + rows[first - 1], first
    if not any(NUMBER.fullmatch(cell) for cell in rows[first][1:]):
        return width, rows[first], first + 1
    return width, [], first


def profile(rows, resolve):
    """Describe a previewed table: its columns, header labels, the samples the labels name and the kind of values it holds."""
    width, header, start = table_header(rows)
    if width < 2:
        return dict(columns=width, labels=0, matched={}, value_kind="none", value_minimum=None, value_maximum=None, first_labels="")
    body = [row for row in rows[start:] if len(row) == width][:CONFIG["inventory"]["preview_rows"]]
    matched = {index: sample for index, sample in ((index, resolve(label)) for index, label in enumerate(header)) if sample}
    # Values are profiled in the sample columns when the header names any, otherwise in every column after the first.
    values = [float(row[index]) for row in body for index in (sorted(matched) or range(1, width)) if NUMBER.fullmatch(row[index])]
    kind = "none" if not values else "whole_numbers" if all(value.is_integer() for value in values) else "decimals"
    return dict(columns=width, labels=len(header), matched=matched, value_kind=kind, value_minimum=min(values, default=None),
                value_maximum=max(values, default=None), first_labels=" | ".join(header[:6])[:200])


def series_matrix(stream, resolve):
    """Read the sample line and the declared row counts of a GEO series matrix and profile the start of its value table."""
    accessions, counts, table, begun = [], [], [], False
    for line in io.TextIOWrapper(stream, encoding="utf-8", errors="replace"):
        cells = [cell.strip('"') for cell in line.rstrip("\r\n").split("\t")]
        if cells[0] == "!series_matrix_table_end" or len(table) > CONFIG["inventory"]["preview_rows"]:
            break
        if begun:
            table.append(cells)
        elif cells[0] == "!Sample_geo_accession":
            accessions = cells[1:]
        elif cells[0] == "!Sample_data_row_count":
            counts = [int(cell) for cell in cells[1:]]
        elif cells[0] == "!series_matrix_table_begin":
            begun = True
    require(accessions and len(counts) == len(accessions) and table and table[0][1:] == accessions,
            "A series matrix must name its samples, their row counts and a value table")
    return dict(profile(table, resolve), declared_rows=max(counts))


def sheet_rows(frame):
    """Return the rows of a workbook sheet as cells of text, without trailing empty cells and without empty rows."""
    import pandas as pd

    rows = [["" if pd.isna(cell) else str(cell).strip() for cell in row] for row in frame.itertuples(index=False, name=None)]
    for row in rows:
        while row and not row[-1]:
            row.pop()
    return [row for row in rows if row]


def workbook_sheets(content, engine):
    """Return the first rows of every sheet of a workbook as cells of text, without trailing empty cells."""
    import pandas as pd

    with warnings.catch_warnings():
        # openpyxl leaves out the worksheet extensions it does not know; they hold no cell values.
        warnings.filterwarnings("ignore", message="Unknown extension is not supported and will be removed", module="openpyxl")
        sheets = pd.read_excel(io.BytesIO(content), sheet_name=None, header=None, dtype=str, engine=engine,
                               nrows=CONFIG["inventory"]["preview_rows"] + 20)
    return {str(title): sheet_rows(frame) for title, frame in sheets.items()}


def preview_rows(stream, head, size):
    """Split the start of a text file into rows; a wide matrix holds few rows in its first block, so the preview grows until it has them."""
    policy = CONFIG["inventory"]
    complete = len(head) < size
    while not complete and head.count(b"\n") <= policy["preview_rows"] and len(head) < policy["maximum_preview_megabytes"] * 1024 * 1024:
        more = stream.read(len(head))
        complete = not more
        head += more
    return text_rows(head, complete)


def describe(opener, member, resolve):
    """Return one record per table a file holds: the file itself, or every sheet of a workbook."""
    policy = CONFIG["inventory"]
    stream, packed = stored(opener)
    size = policy["member_preview_kilobytes" if member else "file_preview_kilobytes"] * 1024
    head = stream.read(size)
    name, blank = signature(head), dict(BLANK, compressed=packed)
    if name == "geo_series_matrix":
        return [dict(blank, format=name, delimiter="tab", **series_matrix(stored(opener)[0], resolve))]
    if name in ("zip", "ole"):
        content = stored(opener)[0].read()
        if name == "zip" and "xl/workbook.xml" not in zipfile.ZipFile(io.BytesIO(content)).namelist():
            return [dict(blank, format="zip_archive")]
        sheets = workbook_sheets(content, "openpyxl" if name == "zip" else "xlrd")
        return [dict(blank, format="workbook_sheet", sheet=title, **profile(rows, resolve)) for title, rows in sheets.items()]
    if name:
        return [dict(blank, format=name)]
    rows, delimiter = preview_rows(stream, head, size)
    found = profile(rows, resolve)
    return [dict(blank, format="table" if found["columns"] > 1 else "text", delimiter=delimiter if found["columns"] > 1 else "", **found)]


def opened(item):
    """Open one listed table, which is a file of its study or a member of one of its archives."""
    path = RAW / item.source / item.accession / item.file
    if not item.member:
        return path.open("rb")
    with tarfile.open(path) as archive:
        return io.BytesIO(archive.extractfile(item.member).read())


def table_labels(item):
    """Read the header of one table the way the inventory previewed it: its labels, the row at which its body starts and its delimiter."""
    stream, _ = stored(lambda: opened(item))
    if item.format == "workbook_sheet":
        content = stream.read()
        rows, delimiter = workbook_sheets(content, "openpyxl" if signature(content[:512]) == "zip" else "xlrd")[item.sheet], ""
    else:
        size = CONFIG["inventory"]["member_preview_kilobytes" if item.member else "file_preview_kilobytes"] * 1024
        rows, delimiter = preview_rows(stream, stream.read(size), size)
    _, header, body = table_header(rows)
    return header, body, delimiter


def table_numbers(item, body, delimiter, positions, identifier=0, decimal="."):
    """Read the identifier column of one table as text and chosen columns as numbers; a cell under the header lines that is neither a number nor left empty stops the run."""
    import pandas as pd

    require(identifier not in positions and len(set(positions)) == len(positions), f"{item.accession}: the columns asked of {item.member or item.file} overlap")
    stream, wanted = stored(lambda: opened(item))[0], [identifier, *positions]
    if item.format == "workbook_sheet":
        content = stream.read()
        with warnings.catch_warnings():
            # openpyxl leaves out the worksheet extensions it does not know; they hold no cell values.
            warnings.filterwarnings("ignore", message="Unknown extension is not supported and will be removed", module="openpyxl")
            sheet = pd.read_excel(io.BytesIO(content), sheet_name=item.sheet, header=None, dtype=str,
                                  engine="openpyxl" if signature(content[:512]) == "zip" else "xlrd")
        frame = pd.DataFrame([[row[position] if position < len(row) else "" for position in wanted] for row in sheet_rows(sheet)[body:]],
                             columns=wanted)
    else:
        text, passed = io.TextIOWrapper(stream, encoding="utf-8-sig", errors="replace"), 0
        for line in text if body else ():
            # The body starts after the rows the preview counted, and the preview leaves out blank lines; a table without a header starts at once.
            passed += bool(line.strip())
            if passed == body:
                break
        frame = pd.read_csv(text, sep=SEPARATORS[delimiter], header=None, usecols=wanted, dtype={identifier: str}, keep_default_na=False,
                            na_values={position: NOT_RECORDED for position in positions}, low_memory=False)[wanted]
    frame.columns = ["identifier", *positions]
    frame["identifier"] = frame["identifier"].fillna("").str.strip()
    read, first = {}, 0
    for position in positions:
        if frame[position].dtype.kind not in "iuf":
            cells = frame[position].fillna("").astype(str).str.strip()
            # A table written with a decimal comma is read so only where the caller says that it is one.
            numbers = pd.to_numeric((cells if decimal == "." else cells.str.replace(decimal, ".", regex=False)).mask(cells.isin(NOT_RECORDED)), errors="coerce")
            read[position] = (cells, numbers)
    if len(frame) and len(read) == len(positions):
        # Rows directly under the header whose chosen cells all hold text are further header lines, not rows of the table; a row of missing values
        # is a row of the table.
        text = pd.concat([numbers.isna() & ~cells.isin(NOT_RECORDED) for cells, numbers in read.values()], axis=1).all(axis=1)
        first = int((~text).to_numpy().argmax())
    for position, (cells, numbers) in read.items():
        odd = cells[first:][numbers[first:].isna() & ~cells[first:].isin(NOT_RECORDED)]
        require(odd.empty, f"{item.accession}: column {position} of {item.member or item.file} holds cells that are not numbers: {odd.unique()[:5].tolist()}")
        frame[position] = numbers
    return frame.iloc[first:].reset_index(drop=True)


def plain(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def study_lookup(own, values):
    """Return how one study names its samples: a resolver for table labels and one for per-sample file names."""
    accessions, minimum = set(own["sample_accession"]), CONFIG["inventory"]["minimum_identifier_characters"]
    names, declared = collections.defaultdict(set), collections.defaultdict(set)
    for row in own.itertuples(index=False):
        files = [name.rsplit("/", 1)[-1] for name in row.raw_files.split(";") if name]
        for name in files:
            declared[name].add(row.sample_accession)
        for value in [row.sample_accession, row.title, row.description, *values.get(row.sample_accession, ()), *files,
                      *(name.split(".")[0] for name in files)]:
            if len(plain(value)) >= minimum:
                names[plain(value)].add(row.sample_accession)
    # An identifier shared by several samples names none of them.
    names = {key: next(iter(found)) for key, found in names.items() if len(found) == 1}
    declared = {key: next(iter(found)) for key, found in declared.items() if len(found) == 1}

    def resolve(label):
        found = SAMPLE.search(label)
        if found and found.group(1) in accessions:
            return found.group(1)
        return names.get(plain(SUFFIX.sub("", label.strip()))) or names.get(plain(label))

    def owner(name):
        name = name.rsplit("/", 1)[-1]
        found = SAMPLE.match(name)
        return declared.get(name) or (found.group(1) if found and found.group(1) in accessions else "")

    return resolve, owner


def list_study(files, resolve, owner):
    """Describe every file of one study and every member of its archives."""
    rows = []
    for item in files.itertuples(index=False):
        path = RAW / item.source / item.accession / item.file
        require(path.is_file() and path.stat().st_size == item.bytes, f"{item.accession}: {item.file} is missing or differs from its recorded size")
        base = dict(source=item.source, accession=item.accession, file=item.file, kind=item.kind)
        if item.file.lower().endswith(ARCHIVES):
            rows.append(dict(base, member="", bytes=item.bytes, sample_accession="", **dict(BLANK, format="tar_archive", compressed=False)))
            with tarfile.open(path) as archive:
                for member in archive:
                    if member.isfile():
                        for found in describe(lambda: archive.extractfile(member), True, resolve):
                            rows.append(dict(base, member=member.name, bytes=member.size, sample_accession=owner(member.name), **found))
        else:
            for found in describe(lambda: path.open("rb"), False, resolve):
                rows.append(dict(base, member="", bytes=item.bytes, sample_accession=owner(item.file), **found))
    return rows


def sample_evidence(samples, records):
    """Collect, for every sample, the files that are its own and the tables that name it."""
    policy = CONFIG["inventory"]
    found = {key: dict(own_files=0, own_formats=set(), tables=[], ncbi=False, matrix=0)
             for key in zip(samples["study_accession"], samples["sample_accession"])}
    for item in records:
        if item["sample_accession"]:
            if item["format"] in policy["per_sample_data_formats"]:
                entry = found[(item["accession"], item["sample_accession"])]
                entry["own_files"] += 1
                entry["own_formats"].add(item["format"])
            continue
        for sample in sorted(set(item["matched"].values())):
            entry = found[(item["accession"], sample)]
            if item["format"] == "geo_series_matrix":
                entry["matrix"] = max(entry["matrix"], item["declared_rows"])
            elif item["kind"] == "ncbi_generated_counts":
                entry["ncbi"] = True
            else:
                entry["tables"].append("/".join(part for part in (item["file"], item["member"], item["sheet"]) if part))
    return found


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-search", action="store_true",
                        help="Ask the repositories again for every platform record instead of reading the stored records under data/raw/platforms")
    args = parser.parse_args()
    print(f"[01] Inventorying expression files v{VERSION}", flush=True)
    start = begin(1)
    external = acquisition_binding()
    policy = CONFIG["platforms"]
    included = pd.read_parquet(ACQUIRED / "included_studies.parquet")
    samples = pd.read_parquet(ACQUIRED / "samples.parquet")
    samples = samples[samples["study_accession"].isin(included["accession"])].reset_index(drop=True)
    attributes = pd.read_parquet(ACQUIRED / "sample_attributes.parquet")
    attributes = attributes[attributes["study_accession"].isin(included["accession"])]
    files = pd.read_parquet(ACQUIRED / "download_files.parquet")
    files = files[files["status"] == "downloaded"]
    require(set(files["accession"]) == set(included["accession"]) == set(samples["study_accession"]) and len(samples) == included["samples"].sum(),
            "Included studies, their samples and their downloaded files do not reconcile")

    names = sorted(set(samples["platform"]) - {""})
    require(all(re.fullmatch(r"GPL\d+|A-[A-Z]+-\d+", name) for name in names), "A sample names a platform that is not one GEO or ArrayExpress accession")
    if PLATFORMS.exists() and args.update_search:
        shutil.rmtree(PLATFORMS)
    print(f"[01] Retrieving {len(names)} platform records", flush=True)
    platforms = pd.DataFrame(harvest_platforms(names, args.update_search))
    family = dict(zip(platforms["platform"], platforms["family"]))
    unrecorded = (samples["platform"] == "") & (samples["library_strategy"] == policy["unrecorded"]["library_strategy"])
    samples["family"] = samples["platform"].map(family).mask(unrecorded, policy["unrecorded"]["family"])
    require(samples["family"].notna().all(), "A sample has neither a platform record nor a recorded sequencing strategy")
    counted = samples.groupby("platform").agg(studies=("study_accession", "nunique"), samples=("sample_accession", "size"))
    platforms = platforms.join(counted, on="platform")[PLATFORM_COLUMNS]

    values = collections.defaultdict(lambda: collections.defaultdict(list))
    for study, sample, value in attributes[["study_accession", "sample_accession", "value"]].itertuples(index=False, name=None):
        values[study][sample].append(value)
    records = []
    for index, accession in enumerate(included["accession"], 1):
        resolve, owner = study_lookup(samples[samples["study_accession"] == accession], values[accession])
        records += list_study(files[files["accession"] == accession], resolve, owner)
        if index % 10 == 0 or index == len(included):
            print(f"[01] described {index:,}/{len(included):,} studies, {len(records):,} files and archive members", flush=True)

    evidence = sample_evidence(samples, records)
    for item in records:
        item["samples_matched"] = len(set(item["matched"].values()))
    listing = pd.DataFrame(records)[FILE_COLUMNS]
    rows = []
    for row in samples.itertuples(index=False):
        entry = evidence[(row.study_accession, row.sample_accession)]
        sources = [label for label, held in (("per_sample_file", entry["own_files"]), ("author_table", entry["tables"]),
                                             ("ncbi_counts", entry["ncbi"]), ("series_matrix", entry["matrix"])) if held]
        rows.append(dict(source=row.source, study_accession=row.study_accession, sample_accession=row.sample_accession, platform=row.platform,
                         family=row.family, own_files=entry["own_files"], own_formats=";".join(sorted(entry["own_formats"])),
                         author_tables=len(entry["tables"]), ncbi_counts=entry["ncbi"], series_matrix_rows=entry["matrix"],
                         data_sources=";".join(sources), author_table_files=";".join(entry["tables"])))
    inventory = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    require(len(inventory) == len(samples) and not inventory.duplicated(["study_accession", "sample_accession"]).any(),
            "Sample inventory does not hold every included sample exactly once")

    titles, summary = dict(zip(platforms["platform"], platforms["title"])), []
    tables = listing[listing["format"].isin(["table", "workbook_sheet"]) & (listing["sample_accession"] == "")
                     & (listing["kind"] != "ncbi_generated_counts")]
    for (source, accession, platform), part in inventory.groupby(["source", "study_accession", "platform"], sort=False):
        formats = collections.Counter(name for names in part["own_formats"] for name in names.split(";") if name)
        described = [f"{'/'.join(name for name in (row.file, row.member, row.sheet) if name)} [{row.value_kind}, {row.samples_matched}/{row.labels} labels]"
                     for row in tables[tables["accession"] == accession].itertuples(index=False)]
        summary.append(dict(source=source, accession=accession, platform=platform, family=part["family"].iloc[0],
                            platform_title=titles.get(platform, ""), samples=len(part), with_own_file=int((part["own_files"] > 0).sum()),
                            own_formats=";".join(f"{name}:{count}" for name, count in sorted(formats.items())),
                            with_author_table=int((part["author_tables"] > 0).sum()), with_ncbi_counts=int(part["ncbi_counts"].sum()),
                            with_series_matrix_values=int((part["series_matrix_rows"] > 0).sum()),
                            unmatched=int((part["data_sources"] == "").sum()), tables="; ".join(described)))
    summary = pd.DataFrame(summary, columns=STUDY_COLUMNS)
    require(summary["samples"].sum() == len(samples), "Study summary does not reconcile with the included samples")

    save_table(platforms, DATA / "platforms.parquet")
    save_table(listing, DATA / "expression_files.parquet")
    save_table(inventory, DATA / "sample_expression.parquet")
    save_table(summary, DATA / "study_platforms.parquet")
    save_review(summary, OUT / "inventory_review.csv")
    counts = lambda frame, column: {str(key): int(value) for key, value in frame[column].value_counts().sort_index().items()}
    finish(1, start, [DATA / "platforms.parquet", DATA / "expression_files.parquet", DATA / "sample_expression.parquet",
                      DATA / "study_platforms.parquet", OUT / "inventory_review.csv", PLATFORMS / "SHA256SUMS.txt"],
           external=bind(external), retrieved_at_utc=datetime.now(timezone.utc).isoformat(), search_mode="updated" if args.update_search else "stored",
           platforms=dict(records=len(platforms), by_family=counts(platforms, "family")),
           samples=dict(total=len(inventory), by_family=counts(inventory, "family"), by_data_sources=counts(inventory, "data_sources"),
                        unmatched=int((inventory["data_sources"] == "").sum())),
           files=dict(listed=len(listing), by_format=counts(listing, "format"),
                      unrecognised=["/".join(name for name in (item["accession"], item["file"], item["member"]) if name)
                                    for item in records if item["format"] in ("unrecognised_binary", "affymetrix_other")]),
           protocol=CONFIG["protocol"])
    print(f"[01] Complete: {len(platforms)} platforms, {len(listing):,} files and archive members, {len(inventory):,} samples; "
          f"{int((inventory['data_sources'] == '').sum()):,} samples that no file or table label names yet; "
          f"review sheet: {(OUT / 'inventory_review.csv').relative_to(ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
