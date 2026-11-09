"""Summarise the CEL files of every Affymetrix study platform by RMA, each study platform on its own."""
import importlib
import shutil
import subprocess
import tarfile
import time
from datetime import datetime

C = importlib.import_module("01_inventory_expression_files")

WORK = C.WORK / "affymetrix"
SCRIPT = C.HERE / C.STAGES[4]
TOOLS = ["oligo", "oligoClasses", "affxparser", "affyio", "preprocessCore"]
SCAN_FORMATS = ("%m/%d/%y %H:%M:%S", "%Y-%m-%dT%H:%M:%SZ")
ARRAY_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "chip_type", "normalization_group", "file", "member", "scan_time",
                 "scan_date_recorded", "raw_pm_median"]
STUDY_COLUMNS = ["source", "accession", "platform", "samples", "chip_type", "package", "package_version", "target", "pm_probes", "features",
                 "value_minimum", "value_maximum", "first_scan", "last_scan", "matrix"]
RMA = r"""
suppressPackageStartupMessages({
  library(jsonlite); library(parallel); library(oligo); library(affxparser); library(affyio); library(preprocessCore)
})
job <- fromJSON(commandArgs(TRUE)[1])
files <- job$files
count <- length(files)
store <- function(index) file.path(job$work, sprintf("%05d.pm", index))

# Every array of a study platform must be the same chip, and the chip must have a recorded design.
headers <- lapply(files, read.celfile.header, info = "full")
chip <- unique(vapply(headers, function(header) header$cdfName, ""))
if (length(chip) != 1) stop("SEVERAL_CHIPS ", paste(chip, collapse = ", "))
design <- job$designs[[chip]]
if (is.null(design)) stop("NO_DESIGN ", chip)
first <- read.celfiles(files[1], verbose = FALSE)
if (annotation(first) != design$package) stop("DESIGN_PACKAGE ", annotation(first))

# The perfect-match probes and their probe sets come from the accessors that oligo::rma itself uses.
if (is(first, "ExpressionFeatureSet")) {
  if (!is.null(design$target)) stop("TARGET_ON_EXPRESSION_ARRAY ", chip)
  probes <- data.frame(fid = pmindex(first), set = probeNames(first), stringsAsFactors = FALSE)
} else {
  if (is.null(design$target)) stop("NO_TARGET ", chip)
  found <- oligo:::stArrayPmInfo(first, target = design$target)
  probes <- data.frame(fid = found$fid, set = as.character(found$fsetid), stringsAsFactors = FALSE)
}
probes <- probes[order(probes$set, probes$fid, method = "radix"), ]
cells <- sort(unique(probes$fid))
place <- match(probes$fid, cells)
size <- nrow(probes)

# Step 1: background correction of each array; the sorted values add up to the normalization target.
corrected <- function(index) {
  raw <- readCel(files[index], indices = cells, readHeader = FALSE, readXY = FALSE, readIntensities = TRUE, readStdvs = FALSE,
                 readPixels = FALSE, readOutliers = FALSE, readMasked = FALSE)$intensities[place]
  values <- rma.background.correct(matrix(as.double(raw), ncol = 1), copy = FALSE)[, 1]
  writeBin(values, store(index), size = 8)
  list(sorted = sort(values), median = median(raw))
}
target <- numeric(size)
medians <- numeric(count)
for (batch in split(seq_len(count), ceiling(seq_len(count) / job$workers))) {
  done <- mclapply(batch, corrected, mc.cores = job$workers, mc.preschedule = FALSE)
  for (position in seq_along(batch)) {
    if (!is.list(done[[position]])) stop("BACKGROUND ", files[batch[position]], " ", as.character(done[[position]]))
    target <- target + done[[position]]$sorted / count
    medians[batch[position]] <- done[[position]]$median
  }
}

# Step 2: quantile normalization of each array to the target of the study platform.
normalised <- function(index) {
  values <- readBin(store(index), "double", size, size = 8)
  values <- normalize.quantiles.use.target(matrix(values, ncol = 1), target, copy = FALSE)[, 1]
  writeBin(values, store(index), size = 8)
  TRUE
}
done <- mclapply(seq_len(count), normalised, mc.cores = job$workers, mc.preschedule = FALSE)
if (!all(vapply(done, isTRUE, TRUE))) stop("NORMALIZATION ", paste(unlist(done[!vapply(done, isTRUE, TRUE)]), collapse = " "))

# Step 3: median polish of the log2 intensities of each probe set, in blocks of whole probe sets.
sets <- unique(probes$set)
starts <- match(sets, probes$set)
ends <- c(starts[-1] - 1, size)
block <- ceiling(ends / max(1, floor(job$block_megabytes * 2^20 / (8 * count))))
result <- matrix(0, length(sets), count)
for (current in unique(block)) {
  inside <- which(block == current)
  from <- starts[inside[1]]
  rows <- ends[inside[length(inside)]] - from + 1
  values <- matrix(0, rows, count)
  for (index in seq_len(count)) {
    connection <- file(store(index), "rb")
    seek(connection, (from - 1) * 8)
    values[, index] <- readBin(connection, "double", rows, size = 8)
    close(connection)
  }
  summary <- subColSummarizeMedianpolishLog(values, probes$set[from:(from + rows - 1)])
  result[inside, ] <- summary[sets[inside], , drop = FALSE]
}

writeBin(as.vector(result), file.path(job$work, "values.bin"), size = 8)
writeLines(sets, file.path(job$work, "features.txt"))
dates <- vapply(headers, function(header) if (length(header$ScanDate) == 1) as.character(header$ScanDate) else "", "")
report <- list(chip_type = chip, package = design$package, package_version = as.character(packageVersion(design$package)),
               feature_set = class(first)[1], pm_probes = size, features = length(sets), scan_dates = dates, raw_pm_medians = medians)
writeLines(toJSON(report, auto_unbox = TRUE, digits = NA), file.path(job$work, "report.json"))
"""


def scan_time(recorded):
    """Write the scan date of a CEL header in ISO form; the CEL generations record it differently, and some do not record it."""
    for pattern in SCAN_FORMATS:
        try:
            return datetime.strptime(recorded.strip(), pattern).strftime("%Y-%m-%dT%H:%M:%S")
        except ValueError:
            continue
    return ""


def stage_arrays(arrays, folder):
    """Give R one plain CEL file per array: the stored file when it is uncompressed and outside an archive, otherwise an unpacked copy."""
    paths, archives = {}, {}
    folder.mkdir(parents=True, exist_ok=True)
    try:
        for position, item in enumerate(arrays):
            path = C.RAW / item.source / item.accession / item.file
            if not item.member and not item.compressed:
                paths[position] = path
                continue
            if item.member and path not in archives:
                archives[path] = tarfile.open(path)
            opener = (lambda: archives[path].extractfile(item.member)) if item.member else (lambda: path.open("rb"))
            paths[position] = folder / f"{position:05d}.CEL"
            with paths[position].open("wb") as target:
                shutil.copyfileobj(C.stored(opener)[0], target, 4 * 1024 * 1024)
    finally:
        for archive in archives.values():
            archive.close()
    return [paths[position] for position in range(len(arrays))]


def run_rma(name, arrays, policy):
    """Run RMA over the arrays of one study platform and return the matrix, its probe sets and what R reported."""
    import numpy as np

    work = WORK / name
    shutil.rmtree(work, ignore_errors=True)
    try:
        files = stage_arrays(arrays, work / "cel")
        C.write_json(work / "job.json", dict(files=[path.as_posix() for path in files], work=work.as_posix(), designs=policy["designs"],
                                             workers=policy["workers"], block_megabytes=policy["block_megabytes"]))
        (work / "rma.R").write_text(RMA, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (work / "rma.R").as_posix(), (work / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"{name}: RMA failed:\n{done.stderr[-3000:]}")
        report = C.read_json(work / "report.json")
        features = (work / "features.txt").read_text(encoding="utf-8").splitlines()
        values = np.fromfile(work / "values.bin", dtype="<f8")
        C.require(len(features) == report["features"] == len(set(features)) and values.size == len(features) * len(arrays)
                  and len(report["scan_dates"]) == len(report["raw_pm_medians"]) == len(arrays), f"{name}: RMA output is incomplete")
        # R writes a matrix column by column, one array after the other.
        values = values.reshape(len(arrays), len(features)).T
        C.require(bool(np.isfinite(values).all()), f"{name}: RMA returned values that are not finite")
        return values, features, report
    finally:
        shutil.rmtree(work, ignore_errors=True)


def run_groups(name, arrays, groups, policy):
    """Run RMA over the arrays of one study platform, or over the arrays of each of its normalization groups on their own, in the order of the arrays."""
    import numpy as np

    if not groups:
        return run_rma(name, arrays, policy)
    parts = []
    for group in sorted(set(groups.values())):
        chosen = [position for position, item in enumerate(arrays) if groups[item.sample_accession] == group]
        parts.append((chosen, *run_rma(f"{name}__{group}", [arrays[position] for position in chosen], policy)))
    features, report = parts[0][2], dict(parts[0][3], scan_dates=[""] * len(arrays), raw_pm_medians=[0.0] * len(arrays))
    values = np.full((len(features), len(arrays)), np.nan)
    for chosen, part, names, found in parts:
        C.require(names == features and found["chip_type"] == report["chip_type"], f"{name}: its normalization groups do not hold the same probe sets")
        values[:, chosen] = part
        for index, position in enumerate(chosen):
            report["scan_dates"][position], report["raw_pm_medians"][position] = found["scan_dates"][index], found["raw_pm_medians"][index]
    return values, features, report


def stored_result(record_path, matrix_path, identity, policy, versions):
    """Return the stored record of a study platform when nothing it depends on has changed, otherwise None."""
    record = C.kept_record(record_path, matrix_path, identity)
    design = policy["designs"].get(record.get("chip_type"), {}) if record else {}
    current = record and record.get("design") == design and record.get("package_version") == versions.get(design.get("package"))
    return record if current else None


def main():
    import pandas as pd

    print(f"[04] Processing Affymetrix arrays v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["affymetrix"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    routes = routes[(routes["family"] == "affymetrix") & (routes["tier"] == "raw_array_files")].sort_values(["samples", "accession", "platform"])
    held = pd.read_parquet(C.DATA / "sample_routes.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing = listing[(listing["format"] == "affymetrix_cel") & (listing["sample_accession"] != "")].copy()
    # A file that the repository links to the wrong record is read as the sample that a recorded decision names; the files of a study only change places.
    for key, entry in C.REGISTER["reassigned_files"].items():
        study = key.split("|")[0]
        moved = (listing["accession"] == study) & listing["member"].isin(entry["files"])
        C.require(int(moved.sum()) == len(entry["files"]) and sorted(listing.loc[moved, "sample_accession"]) == sorted(entry["files"].values()),
                  f"{key}: the reassigned files must be files of the study that only change places between its records")
        listing.loc[moved, "sample_accession"] = listing.loc[moved, "member"].map(entry["files"])
    archives = pd.read_parquet(C.ACQUIRED / "download_files.parquet")
    archives = {(row.accession, row.file): row.sha256 for row in archives[archives["status"] == "downloaded"].itertuples(index=False)}

    packages = sorted({design["package"] for design in policy["designs"].values()})
    versions = C.tool_versions(TOOLS + packages)
    missing = sorted(name for name in TOOLS + packages if not versions.get(name))
    C.require(not missing, f"The image lacks R packages; rebuild it from docker/Dockerfile: {missing}")
    C.MATRICES.mkdir(parents=True, exist_ok=True)
    tools = {name: versions[name] for name in ["R"] + TOOLS}
    method = dict(script=C.sha(SCRIPT), reader=C.sha(C.HERE / C.STAGES[1]), inventory=C.CONFIG["inventory"], runtime=C.runtime(), tools=tools)

    studies, arrays, outputs, computed = [], [], [], 0
    for index, route in enumerate(routes.itertuples(index=False), 1):
        key = f"{route.accession}|{route.platform}"
        group = held[(held["study_accession"] == route.accession) & (held["platform"] == route.platform) & held["in_source"]]
        own = listing[(listing["accession"] == route.accession) & listing["sample_accession"].isin(group["sample_accession"])]
        own = own.sort_values("sample_accession")
        C.require(len(own) == len(group) == route.samples and own["sample_accession"].is_unique,
                  f"{key}: the route does not hold exactly one CEL file for every sample")
        members = list(own.itertuples(index=False))
        groups = C.normalization_groups(key, [item.sample_accession for item in members])
        # A stored matrix stays valid while this script, the inventory script that stages the files and assigns the groups, its inventory settings,
        # the software, the arrays behind it and their normalization groups are the same.
        content = dict(method=method, arrays=[[item.sample_accession, item.file, item.member, archives[(item.accession, item.file)]] for item in members])
        if groups:
            content["groups"] = [[item.sample_accession, groups[item.sample_accession]] for item in members]
        identity = C.digest(content)
        name = f"{route.accession}__{route.platform}"
        record_path, matrix_path, began = C.MATRICES / f"{name}.json", C.MATRICES / f"{name}.parquet", time.monotonic()
        record = stored_result(record_path, matrix_path, identity, policy, versions)
        fresh = record is None
        if fresh:
            C.require(shutil.disk_usage(C.MATRICES).free >= policy["minimum_free_gigabytes"] * 1e9, f"{key}: not enough free disk space to unpack the CEL files")
            values, features, report = run_groups(name, members, groups, policy)
            frame = pd.DataFrame(values.astype("float32"), columns=[item.sample_accession for item in members])
            frame.insert(0, "feature_id", features)
            C.save_table(frame, matrix_path)
            record = dict(identity=identity, accession=route.accession, platform=route.platform, chip_type=report["chip_type"],
                          design=policy["designs"][report["chip_type"]], package_version=report["package_version"],
                          feature_set=report["feature_set"], pm_probes=report["pm_probes"], features=len(features),
                          value_minimum=float(frame.iloc[:, 1:].min().min()), value_maximum=float(frame.iloc[:, 1:].max().max()),
                          arrays=[dict(sample_accession=item.sample_accession, normalization_group=groups[item.sample_accession] if groups else "",
                                       file=item.file, member=item.member, scan_time=scan_time(date), scan_date_recorded=date, raw_pm_median=median)
                                  for item, date, median in zip(members, report["scan_dates"], report["raw_pm_medians"])],
                          tools=tools, matrix=C.record(matrix_path))
            C.write_json(record_path, record)
            computed += 1
        C.require([item["sample_accession"] for item in record["arrays"]] == [item.sample_accession for item in members],
                  f"{key}: the stored record does not list the arrays of the route")
        dates = sorted(item["scan_time"] for item in record["arrays"] if item["scan_time"])
        studies.append(dict(source=route.source, accession=route.accession, platform=route.platform, samples=route.samples,
                            chip_type=record["chip_type"], package=record["design"]["package"], package_version=record["package_version"],
                            target=record["design"].get("target", ""), pm_probes=record["pm_probes"], features=record["features"],
                            value_minimum=record["value_minimum"], value_maximum=record["value_maximum"],
                            first_scan=dates[0] if dates else "", last_scan=dates[-1] if dates else "",
                            matrix=matrix_path.relative_to(C.ROOT).as_posix()))
        # A stored record that names no normalization group is of a study platform normalized as one group.
        arrays += [dict(source=route.source, study_accession=route.accession, platform=route.platform, chip_type=record["chip_type"],
                        **{"normalization_group": "", **item}) for item in record["arrays"]]
        outputs += [matrix_path, record_path]
        print(f"[04] {index}/{len(routes)} {route.accession} {route.platform}: {len(members):,} arrays, {record['chip_type']}, "
              f"{record['features']:,} probe sets, {'computed' if fresh else 'kept'} in {time.monotonic() - began:,.0f} s", flush=True)

    C.clear_work(WORK)
    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)
    arrays = pd.DataFrame(arrays, columns=ARRAY_COLUMNS)
    C.require(len(arrays) == studies["samples"].sum() == routes["samples"].sum() and len(studies) == len(routes)
              and not arrays.duplicated(["study_accession", "sample_accession"]).any(), "Processed arrays do not reconcile with the Affymetrix routes")
    C.save_table(studies, C.DATA / "affymetrix_studies.parquet")
    C.save_table(arrays, C.DATA / "affymetrix_arrays.parquet")
    C.save_review(studies, C.OUT / "affymetrix_review.csv")
    by_chip = {str(chip): dict(study_platforms=len(part), arrays=int(part["samples"].sum()), features=sorted(int(value) for value in set(part["features"])))
               for chip, part in studies.groupby("chip_type")}
    C.finish(4, start, [C.DATA / "affymetrix_studies.parquet", C.DATA / "affymetrix_arrays.parquet", C.OUT / "affymetrix_review.csv"] + outputs,
             method=policy["method"], study_platforms=dict(processed=len(studies), computed_in_this_run=computed), arrays=len(arrays),
             normalization_groups={key: entry for key, entry in C.REGISTER["normalization_groups"].items() if key in set(routes["accession"] + "|" + routes["platform"])},
             reassigned_files=C.REGISTER["reassigned_files"],
             by_chip_type=by_chip, designs=policy["designs"], tools=versions, protocol=C.CONFIG["protocol"])
    print(f"[04] Complete: {len(studies)} study platforms and {len(arrays):,} arrays summarised by RMA ({computed} computed in this run); "
          f"review sheet: {(C.OUT / 'affymetrix_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for chip, item in by_chip.items():
        print(f"[04]   {chip}: {item['study_platforms']} study platforms, {item['arrays']:,} arrays, {item['features']} probe sets", flush=True)


if __name__ == "__main__":
    main()
