"""Normalise the single-channel array files of every Agilent, ImaGene and CodeLink study platform with limma, each study platform on its own."""
import csv
import importlib
import io
import re
import shutil
import subprocess
import tarfile
import time
from datetime import datetime

C = importlib.import_module("01_inventory_expression_files")

WORK = C.WORK / "single_channel"
SCRIPT = C.HERE / C.STAGES[5]
TOOLS = ["limma"]
STARTS = {"TYPE\t": "agilent_feature_extraction", "Begin Header": "imagene", "CodeLink Expression Analysis": "codelink_export"}
AGILENT_BARCODE = re.compile(r"\d+_\d+_\d+")
IMAGENE_FIELD = re.compile(r"^\t([^\t\r\n]+)\t([^\t\r\n]*)\r?$", re.MULTILINE)
CODELINK_RULE = re.compile(r"^-{10,}\r?\n", re.MULTILINE)
CODELINK_SLIDE = re.compile(r"Slide \((.+?)\)")
SCAN_FORMAT = "%m-%d-%Y %H:%M:%S"
ARRAY_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "layout", "file", "member", "in_matrix", "array_id", "scan_time",
                 "date_recorded", "software", "design", "features", "missing_features", "foreground_median", "background_median"]
STUDY_COLUMNS = ["source", "accession", "platform", "samples", "in_matrix", "layout", "software", "design", "features", "measured_features",
                 "probes", "missing_values", "value_minimum", "value_maximum", "first_scan", "last_scan", "matrix"]
LIMMA = r"""
suppressPackageStartupMessages({library(jsonlite); library(limma)})
job <- fromJSON(commandArgs(TRUE)[1])
rows <- job$features
count <- job$arrays
stored <- function(name) matrix(readBin(file.path(job$work, name), "double", rows * count, size = 8), rows, count)
probes <- readLines(file.path(job$work, "probes.txt"))
measured <- readLines(file.path(job$work, "measured.txt")) == "1"
if (length(probes) != rows || length(measured) != rows) stop("FEATURES ", length(probes), " ", length(measured))

# Step 1: normexp background correction of each array from its feature and local background intensities.
corrected <- backgroundCorrect.matrix(stored("foreground.bin"), stored("background.bin"), method = "normexp", normexp.method = "saddle",
                                      offset = job$offset, verbose = FALSE)
# Step 2: quantile normalization across the arrays of the study platform, control features included, then log2,
# in the order limma applies to single-channel intensities.
normalised <- log2(normalizeBetweenArrays(corrected, method = "quantile"))
# Step 3: control features are left out, and the features of one probe are averaged.
result <- avereps(normalised[measured, , drop = FALSE], ID = probes[measured])

writeBin(as.vector(result), file.path(job$work, "values.bin"), size = 8)
writeLines(rownames(result), file.path(job$work, "features.txt"))
"""


def scan_time(recorded):
    """Write the scan date of an Agilent file in ISO form; the other layouts do not record a scan date."""
    try:
        return datetime.strptime(recorded.strip(), SCAN_FORMAT).strftime("%Y-%m-%dT%H:%M:%S")
    except ValueError:
        return ""


def split_file(text):
    """Name the layout of an array file and return what its header records and the text of its feature table."""
    name = next((layout for start, layout in STARTS.items() if text.startswith(start)), None)
    C.require(name is not None, "the file is in none of the known single-channel layouts")
    if name == "agilent_feature_extraction":
        cut = text.index("\nFEATURES\t") + 1
        rows = [line.split("\t") for line in text[:cut].splitlines()]
        position = next(index for index, row in enumerate(rows) if row[0] == "FEPARAMS")
        # One deposited header lacks two cell breaks in its middle, so the scan and grid fields are read from the left and the version from the right.
        found, ending = dict(zip(rows[position], rows[position + 1])), dict(zip(rows[position][::-1], rows[position + 1][::-1]))
        barcode = AGILENT_BARCODE.match(found.get("FeatureExtractor_Barcode", ""))
        facts = dict(array_id=barcode.group(0) if barcode else found.get("FeatureExtractor_Barcode", ""), scan_time=scan_time(found.get("Scan_Date", "")),
                     date_recorded=found.get("Scan_Date", ""), software=f"Feature Extraction {ending.get('FeatureExtractor_Version', '')}".strip(),
                     design=found.get("Grid_Name", ""))
        return name, facts, text[cut:]
    if name == "imagene":
        head, body = text.split("Begin Raw Data", 1)
        found = dict(IMAGENE_FIELD.findall(head))
        image = re.split(r"[\\/]", found.get("Image File", ""))[-1].rsplit(".", 1)[0]
        facts = dict(array_id=image, scan_time="", date_recorded=found.get("Date", ""), software=f"ImaGene {found.get('version', '')}".strip(), design="")
        return name, facts, body.split("End Raw Data", 1)[0].lstrip("\r\n")
    head, body = CODELINK_RULE.split(text, maxsplit=1)
    lines = [line.strip() for line in head.splitlines()]
    slide = CODELINK_SLIDE.search(head)
    design = next((line.split("\t", 1)[1].strip() for line in lines if line.startswith("LAYOUT\t")), "")
    return name, dict(array_id=slide.group(1) if slide else "", scan_time="", date_recorded="", software=lines[0], design=design), body


def read_array(content, layouts):
    """Return the layout of one array file, what its header records, and its features with their probe, kind and intensities."""
    import numpy as np
    import pandas as pd

    name, facts, body = split_file(content.decode("latin-1"))
    C.require(name in layouts, f"no columns are recorded for the layout {name}")
    layout = layouts[name]
    wanted = [*layout["feature"], layout["probe"], layout["kind"], layout["foreground"], layout["background"]]
    # A file that was cut before its end may stop inside a line; a last line with fewer cells than the header is not a feature.
    body = body.rstrip("\r\n")
    final = body.rfind("\n")
    if final > 0 and body[final + 1:].count("\t") < body[:body.index("\n")].count("\t"):
        body = body[:final]
    table = pd.read_csv(io.StringIO(body), sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE,
                        usecols=lambda label: label.strip() in wanted).fillna("")
    table.columns = [label.strip() for label in table.columns]
    C.require(set(wanted) <= set(table.columns), f"the feature table lacks {sorted(set(wanted) - set(table.columns))}")
    # A line without a feature position is padding after the last feature.
    table = table[table[layout["feature"][0]].str.strip() != ""]
    position = table[layout["feature"][0]].str.strip()
    for column in layout["feature"][1:]:
        position = position + "|" + table[column].str.strip()
    features = pd.DataFrame(dict(key=position, probe=table[layout["probe"]].str.strip(), kind=table[layout["kind"]].str.strip(),
                                 foreground=pd.to_numeric(table[layout["foreground"]]), background=pd.to_numeric(table[layout["background"]])))
    if "masked_value" in layout:
        masked = (features["foreground"] == layout["masked_value"]) | (features["background"] == layout["masked_value"])
        features.loc[masked, ["foreground", "background"]] = np.nan
    C.require(len(features) > 0 and features["key"].is_unique, "the features of the file are not distinct")
    return name, facts, features.reset_index(drop=True)


def read_study(key, members, layouts):
    """Read every array of one study platform: the features of the array design and, for each file, what it records and its intensities."""
    archives, reference, layout, facts, vectors = {}, None, None, [], []
    try:
        for item in members:
            path, label = C.RAW / item.source / item.accession / item.file, item.member or item.file
            if item.member and path not in archives:
                archives[path] = tarfile.open(path)
            opener = (lambda: archives[path].extractfile(item.member)) if item.member else (lambda: path.open("rb"))
            try:
                name, found, features = read_array(C.stored(opener)[0].read(), layouts)
            except (ValueError, KeyError, StopIteration) as problem:
                raise ValueError(f"{key}: {label} cannot be read: {problem}") from problem
            layout = layout or name
            C.require(name == layout, f"{key}: {label} is in layout {name}, the first array in {layout}")
            described = features[["key", "probe", "kind"]]
            if reference is not None and len(described) == len(reference) and not described["key"].equals(reference["key"]):
                C.require(set(described["key"]) == set(reference["key"]), f"{key}: {label} does not hold the features of the other arrays")
                features = features.set_index("key").loc[reference["key"]].reset_index()
                described = features[["key", "probe", "kind"]]
            # A file may end before the last feature of the array; the features it holds must still be the first ones of the design.
            shorter, longer = sorted((described, described if reference is None else reference), key=len)
            C.require(shorter.equals(longer.iloc[:len(shorter)]), f"{key}: {label} does not hold the features of the other arrays")
            reference = longer
            vectors.append((features["foreground"].to_numpy(dtype="float64"), features["background"].to_numpy(dtype="float64")))
            facts.append(dict(found, features=len(features)))
    finally:
        for archive in archives.values():
            archive.close()
    return layout, reference, facts, vectors


def run_limma(name, reference, measured, foreground, background, policy):
    """Normalise the arrays of one study platform with limma and return the matrix and its probes."""
    import numpy as np

    work = WORK / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        # R reads a matrix column by column, one array after the other.
        np.ascontiguousarray(foreground.T).tofile(work / "foreground.bin")
        np.ascontiguousarray(background.T).tofile(work / "background.bin")
        (work / "probes.txt").write_text("".join(f"{probe}\n" for probe in reference["probe"]), encoding="utf-8", newline="\n")
        (work / "measured.txt").write_text("".join(f"{int(flag)}\n" for flag in measured), encoding="utf-8", newline="\n")
        C.write_json(work / "job.json", dict(work=work.as_posix(), features=foreground.shape[0], arrays=foreground.shape[1], offset=policy["offset"]))
        (work / "limma.R").write_text(LIMMA, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (work / "limma.R").as_posix(), (work / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"{name}: limma failed:\n{done.stderr[-3000:]}")
        probes = (work / "features.txt").read_text(encoding="utf-8").splitlines()
        values = np.fromfile(work / "values.bin", dtype="<f8")
        C.require(len(probes) == len(set(probes)) == reference.loc[measured, "probe"].nunique() and values.size == len(probes) * foreground.shape[1],
                  f"{name}: the limma output is incomplete")
        return values.reshape(foreground.shape[1], len(probes)).T, probes
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    import numpy as np
    import pandas as pd

    print(f"[05] Processing single-channel arrays v{C.VERSION}", flush=True)
    start = C.begin(5)
    policy, formats = C.CONFIG["single_channel"], C.CONFIG["routes"]["tiers"]["raw_array_files"]["formats"]
    incomplete_files = C.REGISTER["incomplete_files"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    routes = routes[routes["family"].isin(policy["families"]) & (routes["tier"] == "raw_array_files")].sort_values(["samples", "accession", "platform"])
    unknown = sorted(set(incomplete_files) - set(routes["accession"] + "|" + routes["platform"]))
    C.require(not unknown, f"Incomplete files recorded for study platforms that are not single-channel raw-array routes: {unknown}")
    held = pd.read_parquet(C.DATA / "sample_routes.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing = listing[listing["sample_accession"] != ""]
    recorded = pd.read_parquet(C.ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", "channel_count"])
    archives = pd.read_parquet(C.ACQUIRED / "download_files.parquet")
    archives = {(row.accession, row.file): row.sha256 for row in archives[archives["status"] == "downloaded"].itertuples(index=False)}

    versions = C.tool_versions(TOOLS)
    missing = sorted(name for name in TOOLS if not versions.get(name))
    C.require(not missing, f"The image lacks R packages; rebuild it from docker/Dockerfile: {missing}")
    C.MATRICES.mkdir(parents=True, exist_ok=True)
    method = dict(script=C.sha(SCRIPT), reader=C.sha(C.HERE / C.STAGES[1]), inventory=C.CONFIG["inventory"], layouts=policy["layouts"],
                  offset=policy["offset"], tools=versions, runtime=C.runtime())

    studies, arrays, outputs, computed = [], [], [], 0
    for index, route in enumerate(routes.itertuples(index=False), 1):
        key = f"{route.accession}|{route.platform}"
        group = held[(held["study_accession"] == route.accession) & (held["platform"] == route.platform) & held["in_source"]]
        own = listing[(listing["accession"] == route.accession) & listing["sample_accession"].isin(group["sample_accession"])
                      & listing["format"].isin(formats[route.family])].sort_values("sample_accession")
        C.require(len(own) == len(group) == route.samples and own["sample_accession"].is_unique,
                  f"{key}: the route does not hold exactly one array file for every sample")
        channels = recorded[(recorded["study_accession"] == route.accession) & recorded["sample_accession"].isin(group["sample_accession"])]
        C.require(len(channels) == len(group) and not (pd.to_numeric(channels["channel_count"], errors="coerce") > 1).any(),
                  f"{key}: a sample record declares more than one channel")
        members = list(own.itertuples(index=False))
        expected = sorted(incomplete_files.get(key, {}).get("samples", []))
        # A stored matrix stays valid while this script, the inventory script and its inventory settings, the layouts, the software and the files behind
        # it are the same.
        identity = C.digest(dict(method=method, incomplete=expected, arrays=[[item.sample_accession, item.file, item.member, archives[(item.accession, item.file)]]
                                                                             for item in members]))
        name = f"{route.accession}__{route.platform}"
        record_path, matrix_path, began = C.MATRICES / f"{name}.json", C.MATRICES / f"{name}.parquet", time.monotonic()
        record = C.kept_record(record_path, matrix_path, identity)
        fresh = record is None
        if fresh:
            layout, reference, facts, vectors = read_study(key, members, policy["layouts"])
            whole = [found["features"] == len(reference) for found in facts]
            incomplete = sorted(item.sample_accession for item, inside in zip(members, whole) if not inside)
            C.require(incomplete == expected, f"{key}: the array files that end before the last feature are {incomplete}; "
                                              f"incomplete_files of {C.CONFIG['register']} records {expected}")
            measured = reference["kind"].isin(policy["layouts"][layout]["measured"]).to_numpy()
            C.require(measured.any() and (reference.loc[measured, "probe"] != "").all(), f"{key}: measured features without a probe name")
            foreground = np.column_stack([vectors[position][0] for position, inside in enumerate(whole) if inside])
            background = np.column_stack([vectors[position][1] for position, inside in enumerate(whole) if inside])
            absent = np.isnan(foreground) | np.isnan(background)
            C.require("masked_value" in policy["layouts"][layout] or not absent.any(), f"{key}: intensities are missing in a layout that masks none")
            values, probes = run_limma(name, reference, measured, foreground, background, policy)
            C.require(bool((np.isfinite(values) | np.isnan(values)).all()) and int(np.isnan(values).sum()) <= int(absent[measured].sum()),
                      f"{key}: limma returned values that are not finite where the arrays hold intensities")
            frame = pd.DataFrame(values.astype("float32"), columns=[item.sample_accession for item, inside in zip(members, whole) if inside])
            frame.insert(0, "feature_id", probes)
            C.save_table(frame, matrix_path)
            record = dict(identity=identity, accession=route.accession, platform=route.platform, layout=layout, offset=policy["offset"],
                          columns=policy["layouts"][layout], features=len(reference), measured_features=int(measured.sum()), probes=len(probes),
                          missing_values=int(np.isnan(values).sum()), value_minimum=float(np.nanmin(values)), value_maximum=float(np.nanmax(values)),
                          arrays=[dict(sample_accession=item.sample_accession, file=item.file, member=item.member, in_matrix=inside, **found,
                                       missing_features=int((np.isnan(first) | np.isnan(second)).sum()),
                                       foreground_median=float(np.nanmedian(first)), background_median=float(np.nanmedian(second)))
                                  for item, inside, found, (first, second) in zip(members, whole, facts, vectors)],
                          tools=versions, matrix=C.record(matrix_path))
            C.write_json(record_path, record)
            computed += 1
        C.require([item["sample_accession"] for item in record["arrays"]] == [item.sample_accession for item in members]
                  and sorted(item["sample_accession"] for item in record["arrays"] if not item["in_matrix"]) == expected,
                  f"{key}: the stored record does not list the arrays of the route")
        inside = [item for item in record["arrays"] if item["in_matrix"]]
        dates = sorted(item["scan_time"] for item in inside if item["scan_time"])
        studies.append(dict(source=route.source, accession=route.accession, platform=route.platform, samples=route.samples, in_matrix=len(inside),
                            layout=record["layout"], software="; ".join(sorted({item["software"] for item in inside})),
                            design="; ".join(sorted({item["design"] for item in inside} - {""})), features=record["features"],
                            measured_features=record["measured_features"], probes=record["probes"], missing_values=record["missing_values"],
                            value_minimum=record["value_minimum"], value_maximum=record["value_maximum"],
                            first_scan=dates[0] if dates else "", last_scan=dates[-1] if dates else "",
                            matrix=matrix_path.relative_to(C.ROOT).as_posix()))
        arrays += [dict(source=route.source, study_accession=route.accession, platform=route.platform, layout=record["layout"], **item)
                   for item in record["arrays"]]
        outputs += [matrix_path, record_path]
        print(f"[05] {index}/{len(routes)} {route.accession} {route.platform}: {len(inside):,} of {len(members):,} arrays, {record['layout']}, "
              f"{record['probes']:,} probes, {'computed' if fresh else 'kept'} in {time.monotonic() - began:,.0f} s", flush=True)

    C.clear_work(WORK)
    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)
    arrays = pd.DataFrame(arrays, columns=ARRAY_COLUMNS)
    C.require(len(arrays) == studies["samples"].sum() == routes["samples"].sum() and len(studies) == len(routes)
              and arrays["in_matrix"].sum() == studies["in_matrix"].sum() and not arrays.duplicated(["study_accession", "sample_accession"]).any(),
              "Processed arrays do not reconcile with the single-channel routes")
    # A slide identifier is typed by hand at some scanners, so an array counts as deposited twice only when its intensities agree as well.
    repeated = arrays[(arrays["array_id"] != "") & arrays.duplicated(["layout", "array_id", "foreground_median", "background_median"], keep=False)]
    C.save_table(studies, C.DATA / "single_channel_studies.parquet")
    C.save_table(arrays, C.DATA / "single_channel_arrays.parquet")
    C.save_review(studies, C.OUT / "single_channel_review.csv")
    by_layout = {str(layout): dict(study_platforms=len(part), arrays=int(part["samples"].sum()), in_matrix=int(part["in_matrix"].sum()))
                 for layout, part in studies.groupby("layout")}
    C.finish(5, start, [C.DATA / "single_channel_studies.parquet", C.DATA / "single_channel_arrays.parquet", C.OUT / "single_channel_review.csv"] + outputs,
             method=policy["method"], study_platforms=dict(processed=len(studies), computed_in_this_run=computed),
             arrays=dict(routed=len(arrays), in_matrix=int(arrays["in_matrix"].sum()), incomplete_file=int((~arrays["in_matrix"]).sum())),
             by_layout=by_layout, arrays_deposited_more_than_once=dict(arrays=len(repeated), identifiers=int(repeated["array_id"].nunique()),
                                                                       studies=sorted(set(repeated["study_accession"])),
                                                                       basis="same slide identifier and the same median feature and background intensity"),
             layouts=policy["layouts"], offset=policy["offset"], incomplete_files=incomplete_files, tools=versions, protocol=C.CONFIG["protocol"])
    print(f"[05] Complete: {len(studies)} study platforms normalised with limma ({computed} computed in this run); {int(arrays['in_matrix'].sum()):,} of "
          f"{len(arrays):,} arrays are in a matrix; review sheet: {(C.OUT / 'single_channel_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for layout, item in by_layout.items():
        print(f"[05]   {layout}: {item['study_platforms']} study platforms, {item['in_matrix']:,} of {item['arrays']:,} arrays", flush=True)
    print(f"[05]   {len(repeated):,} array records carry the slide identifier and the median intensities of another record "
          f"({', '.join(sorted(set(repeated['study_accession']))) or 'none'})", flush=True)


if __name__ == "__main__":
    main()
