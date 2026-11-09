"""Normalise the probe intensities that studies deposit before normalization, in tables or in one file per array, each study platform on its own."""
import importlib
import re
import shutil
import subprocess
import time

C = importlib.import_module("01_inventory_expression_files")

WORK = C.WORK / "intensity_tables"
SCRIPT = C.HERE / C.STAGES[6]
TOOLS = ["limma"]
TABLES = ("table", "workbook_sheet")
BEADCHIP = re.compile(r"(?<!\d)(\d{10,12}_[A-L])(?![A-Za-z0-9])")
ARRAY_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "table", "label", "array_id", "correction", "signal_median",
                 "signal_percentile_99", "negative_share", "detected_share", "background_mean", "background_sd", "signal_mean"]
STUDY_COLUMNS = ["source", "accession", "platform", "family", "samples", "in_matrix", "tables", "correction", "probes", "value_minimum",
                 "value_maximum", "matrix"]
LIMMA = r"""
suppressPackageStartupMessages({library(jsonlite); library(limma)})
job <- fromJSON(commandArgs(TRUE)[1])
rows <- job$features
count <- job$arrays
stored <- function(name) matrix(readBin(file.path(job$work, name), "double", rows * count, size = 8), rows, count)
signal <- stored("signal.bin")
fitted <- matrix(NaN, count, 3)

if (job$correction == "detection_p_values") {
  # The negative control probes are not deposited, so limma infers their mean and spread on each array from its detection p-values;
  # neqc then corrects the background by the normexp model, adds the offset, quantile-normalises the arrays and takes log2.
  detection <- stored("detection.bin")
  fitted <- normexp.fit.detection.p(signal, detection)
  result <- neqc(signal, detection.p = detection, offset = job$offset)
} else if (job$correction == "intensity_distribution") {
  # Without detection p-values the background mean and spread of each array are estimated from the distribution of its intensities,
  # followed by the same normexp correction, offset, quantile normalization and log2.
  for (index in seq_len(count)) fitted[index, ] <- normexp.fit(signal[, index], method = job$estimator)$par
  corrected <- backgroundCorrect.matrix(signal, method = "normexp", normexp.method = job$estimator, offset = job$offset, verbose = FALSE)
  result <- log2(normalizeBetweenArrays(corrected, method = "quantile"))
} else {
  # Signals that the scanner software has already background-corrected are quantile-normalised and put on the log2 scale.
  result <- log2(normalizeBetweenArrays(signal, method = "quantile"))
}

writeBin(as.vector(result), file.path(job$work, "values.bin"), size = 8)
writeBin(as.vector(fitted), file.path(job$work, "fitted.bin"), size = 8)
"""


def measure(label):
    """Return the measure a column label ends with, without case and punctuation; a label that is a name alone has none."""
    return "".join(re.findall(r"[a-z0-9]+", label[len(C.SUFFIX.sub("", label)):].lower()))


def table_columns(key, route, labels, listing, names):
    """Find, for every sample that a table of the route names, its intensity column and its detection column in the label map."""
    mine = labels[(labels["accession"] == route.accession) & labels["table"].isin(route.data_source.split(";"))
                  & (labels["platform"] == route.platform) & (labels["sample_accession"] != "")]
    tables, columns = {}, []
    for sample, part in mine.groupby("sample_accession"):
        signal, detection = part[part["measure"].isin(names["signal"])], part[part["measure"].isin(names["detection"])]
        C.require(len(signal) == 1 and len(detection) <= 1 and set(detection["table"]) <= set(signal["table"]),
                  f"{key}: sample {sample} does not have one intensity column and at most one detection column in one table")
        table = signal["table"].iloc[0]
        if table not in tables:
            found = listing[(listing["accession"] == route.accession) & (listing["table"] == table) & listing["format"].isin(TABLES)]
            C.require(len(found) == 1, f"{key}: the inventory does not list the table {table} exactly once")
            item = next(found.itertuples(index=False))
            tables[table] = (item, *C.table_labels(item))
        item, header, body, delimiter = tables[table]
        C.require(all(row.position < len(header) and header[row.position] == row.label for row in [*signal.itertuples(), *detection.itertuples()]),
                  f"{key}: the columns of {table} are not where the label map recorded them")
        columns.append(dict(sample_accession=sample, table=table, item=item, body=body, delimiter=delimiter, label=signal["label"].iloc[0],
                            signal=int(signal["position"].iloc[0]), detection=int(detection["position"].iloc[0]) if len(detection) else None))
    return columns


def file_columns(key, route, samples, listing, formats, names):
    """Find, in the one table of every array, its intensity column and its detection column by the measure their labels end with."""
    columns = []
    for sample in sorted(samples):
        found = listing[(listing["accession"] == route.accession) & (listing["sample_accession"] == sample) & listing["format"].isin(formats)
                        & (listing["value_kind"] != "none")]
        C.require(len(found) == 1, f"{key}: sample {sample} does not have exactly one table of its own")
        item = next(found.itertuples(index=False))
        header, body, delimiter = C.table_labels(item)
        where = {kind: [position for position, label in enumerate(header) if position and measure(label) and measure(label) in names[kind]]
                 for kind in ("signal", "detection")}
        C.require(len(where["signal"]) == 1 and len(where["detection"]) <= 1,
                  f"{key}: {item.table} does not hold one intensity column and at most one detection column")
        columns.append(dict(sample_accession=sample, table=item.table, item=item, body=body, delimiter=delimiter, label=header[where["signal"][0]],
                            signal=where["signal"][0], detection=where["detection"][0] if where["detection"] else None))
    return columns


def read_columns(key, columns):
    """Read the intensities and detection p-values of every array of one study platform, with the probes of the first table as the rows."""
    import numpy as np

    tables, probes, signal, detection = {}, None, {}, {}
    for column in columns:
        tables.setdefault(column["table"], []).append(column)
    for table, group in tables.items():
        first = group[0]
        positions = sorted({column[kind] for column in group for kind in ("signal", "detection") if column[kind] is not None})
        frame = C.table_numbers(first["item"], first["body"], first["delimiter"], positions)
        C.require(len(frame) > 0 and (frame["identifier"] != "").all() and frame["identifier"].is_unique,
                  f"{key}: a probe identifier of {table} is empty or repeated")
        if probes is None:
            probes = frame["identifier"].tolist()
        elif frame["identifier"].tolist() != probes:
            # Tables of one study platform may list the probes in another order, but they must hold the same probes.
            C.require(set(frame["identifier"]) == set(probes), f"{key}: {table} does not hold the probes of the other tables of the study platform")
            frame = frame.set_index("identifier", drop=False).loc[probes]
        for column in group:
            signal[column["sample_accession"]] = frame[column["signal"]].to_numpy(dtype="float64")
            if column["detection"] is not None:
                detection[column["sample_accession"]] = frame[column["detection"]].to_numpy(dtype="float64")
    order = [column["sample_accession"] for column in columns]
    C.require(len(detection) in (0, len(order)), f"{key}: only some arrays of the study platform have detection p-values")
    return probes, np.column_stack([signal[sample] for sample in order]), np.column_stack([detection[sample] for sample in order]) if detection else None


def run_limma(name, signal, detection, correction, policy):
    """Normalise the arrays of one study platform with limma and return the matrix and the normexp parameters fitted to each array."""
    import numpy as np

    work = WORK / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        # R reads a matrix column by column, one array after the other.
        np.ascontiguousarray(signal.T).tofile(work / "signal.bin")
        if detection is not None:
            np.ascontiguousarray(detection.T).tofile(work / "detection.bin")
        C.write_json(work / "job.json", dict(work=work.as_posix(), features=signal.shape[0], arrays=signal.shape[1], correction=correction,
                                             offset=policy["offset"], estimator=policy["estimator_without_detection"]))
        (work / "limma.R").write_text(LIMMA, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (work / "limma.R").as_posix(), (work / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"{name}: limma failed:\n{done.stderr[-3000:]}")
        values, fitted = np.fromfile(work / "values.bin", dtype="<f8"), np.fromfile(work / "fitted.bin", dtype="<f8")
        C.require(values.size == signal.size and fitted.size == signal.shape[1] * 3, f"{name}: the limma output is incomplete")
        return values.reshape(signal.shape[1], signal.shape[0]).T, fitted.reshape(3, signal.shape[1]).T
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    import numpy as np
    import pandas as pd

    print(f"[06] Processing intensity tables v{C.VERSION}", flush=True)
    start = C.begin(6)
    policy, formats = C.CONFIG["intensity_tables"], C.CONFIG["routes"]["tiers"]["unnormalized_intensities"]["formats"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    routes = routes[routes["tier"] == "unnormalized_intensities"].sort_values(["samples", "accession", "platform"])
    C.require(set(routes["family"]) <= set(policy["columns"]), f"No intensity columns are recorded for {sorted(set(routes['family']) - set(policy['columns']))}")
    held = pd.read_parquet(C.DATA / "sample_routes.parquet")
    assigned = pd.read_parquet(C.DATA / "table_samples.parquet")
    labels = pd.read_parquet(C.DATA / "label_map.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing["table"] = ["/".join(name for name in names if name) for names in zip(listing["file"], listing["member"], listing["sheet"])]
    archives = pd.read_parquet(C.ACQUIRED / "download_files.parquet")
    archives = {(row.accession, row.file): row.sha256 for row in archives[archives["status"] == "downloaded"].itertuples(index=False)}

    versions = C.tool_versions(TOOLS)
    missing = sorted(name for name in TOOLS if not versions.get(name))
    C.require(not missing, f"The image lacks R packages; rebuild it from docker/Dockerfile: {missing}")
    C.MATRICES.mkdir(parents=True, exist_ok=True)
    method = dict(script=C.sha(SCRIPT), reader=C.sha(C.HERE / C.STAGES[1]), inventory=C.CONFIG["inventory"], offset=policy["offset"],
                  estimator=policy["estimator_without_detection"], detected_below=policy["detected_below"], tools=versions, runtime=C.runtime())

    studies, arrays, outputs, computed = [], [], [], 0
    for index, route in enumerate(routes.itertuples(index=False), 1):
        key, names = f"{route.accession}|{route.platform}", policy["columns"][route.family]
        if route.data_source == "per_sample_files":
            group = held[(held["study_accession"] == route.accession) & (held["platform"] == route.platform) & held["in_source"]]
            columns = file_columns(key, route, group["sample_accession"], listing[listing["sample_accession"] != ""], formats.get(route.family, []), names)
        else:
            group = assigned[(assigned["study_accession"] == route.accession) & (assigned["platform"] == route.platform) & assigned["in_table"]]
            columns = table_columns(key, route, labels, listing[listing["sample_accession"] == ""], names)
        C.require(columns and [column["sample_accession"] for column in columns] == sorted(group["sample_accession"]),
                  f"{key}: the samples with an intensity column are not the samples that the data source of the route holds")
        # A stored matrix stays valid while this script, the inventory script that reads the tables and its inventory settings, the software, the files
        # and the columns behind it are the same.
        identity = C.digest(dict(method=method, columns=names, arrays=[[column["sample_accession"], column["table"], column["signal"], column["detection"],
                                                                       archives[(route.accession, column["item"].file)]] for column in columns]))
        name = f"{route.accession}__{route.platform}"
        record_path, matrix_path, began = C.MATRICES / f"{name}.json", C.MATRICES / f"{name}.parquet", time.monotonic()
        record = C.kept_record(record_path, matrix_path, identity)
        fresh = record is None
        if fresh:
            probes, signal, detection = read_columns(key, columns)
            C.require(bool(np.isfinite(signal).all()) and (detection is None or bool(np.isfinite(detection).all())),
                      f"{key}: an intensity or a detection p-value is missing")
            correction = "none" if names["correction"] == "none" else "detection_p_values" if detection is not None else "intensity_distribution"
            C.require(correction != "none" or bool((signal > 0).all()), f"{key}: background-corrected signals must be above zero")
            if detection is not None:
                # limma decides from the first array whether small or large detection values mark bright probes, so every array must agree.
                across = np.arange(signal.shape[1])
                bright, dark = detection[signal.argmax(axis=0), across], detection[signal.argmin(axis=0), across]
                C.require(bool(((detection >= 0) & (detection <= 1)).all()) and bool((bright < dark).all()),
                          f"{key}: the detection values are not p-values that fall as the intensity rises")
            values, fitted = run_limma(name, signal, detection, correction, policy)
            C.require(bool(np.isfinite(values).all()), f"{key}: limma returned values that are not finite")
            frame = pd.DataFrame(values.astype("float32"), columns=[column["sample_accession"] for column in columns])
            frame.insert(0, "feature_id", probes)
            C.save_table(frame, matrix_path)
            number = lambda value: None if np.isnan(value) else float(value)
            record = dict(identity=identity, accession=route.accession, platform=route.platform, family=route.family, correction=correction,
                          offset=None if correction == "none" else policy["offset"],
                          estimator=policy["estimator_without_detection"] if correction == "intensity_distribution" else None,
                          columns=names, tables=list(dict.fromkeys(column["table"] for column in columns)),
                          probes=len(probes), value_minimum=float(values.min()), value_maximum=float(values.max()),
                          arrays=[dict(sample_accession=column["sample_accession"], table=column["table"], label=column["label"],
                                       array_id=(BEADCHIP.search(column["label"]) or [""])[0], signal_median=float(np.median(signal[:, position])),
                                       signal_percentile_99=float(np.quantile(signal[:, position], 0.99)),
                                       negative_share=float(np.mean(signal[:, position] < 0)),
                                       detected_share=None if detection is None else float(np.mean(detection[:, position] < policy["detected_below"])),
                                       background_mean=number(fitted[position, 0]), background_sd=number(np.exp(fitted[position, 1])),
                                       signal_mean=number(np.exp(fitted[position, 2])))
                                  for position, column in enumerate(columns)],
                          tools=versions, matrix=C.record(matrix_path))
            C.write_json(record_path, record)
            computed += 1
        C.require([item["sample_accession"] for item in record["arrays"]] == [column["sample_accession"] for column in columns],
                  f"{key}: the stored record does not list the arrays of the route")
        studies.append(dict(source=route.source, accession=route.accession, platform=route.platform, family=route.family, samples=route.samples,
                            in_matrix=len(record["arrays"]), tables=len(record["tables"]), correction=record["correction"], probes=record["probes"],
                            value_minimum=record["value_minimum"], value_maximum=record["value_maximum"],
                            matrix=matrix_path.relative_to(C.ROOT).as_posix()))
        arrays += [dict(source=route.source, study_accession=route.accession, platform=route.platform, correction=record["correction"], **item)
                   for item in record["arrays"]]
        outputs += [matrix_path, record_path]
        print(f"[06] {index}/{len(routes)} {route.accession} {route.platform}: {len(record['arrays']):,} of {route.samples:,} arrays, "
              f"{record['correction']}, {record['probes']:,} probes, {'computed' if fresh else 'kept'} in {time.monotonic() - began:,.0f} s", flush=True)

    C.clear_work(WORK)
    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)
    arrays = pd.DataFrame(arrays, columns=ARRAY_COLUMNS)
    C.require(len(studies) == len(routes) and len(arrays) == studies["in_matrix"].sum() and (studies["in_matrix"] <= studies["samples"]).all()
              and not arrays.duplicated(["study_accession", "sample_accession"]).any(), "Processed arrays do not reconcile with the intensity routes")
    # A BeadChip position is read by the scanner, so one identifier is one array; the median intensity shows whether both deposits hold the same
    # export, and deposits round their intensities to different numbers of decimals.
    repeated = arrays[(arrays["array_id"] != "") & arrays.duplicated("array_id", keep=False)]
    same_export = repeated[repeated.groupby("array_id")["signal_median"].transform(lambda medians: medians.max() - medians.min())
                           <= policy["same_median_intensity"]]
    C.save_table(studies, C.DATA / "intensity_table_studies.parquet")
    C.save_table(arrays, C.DATA / "intensity_table_arrays.parquet")
    C.save_review(studies, C.OUT / "intensity_tables_review.csv")
    by_correction = {str(correction): dict(study_platforms=len(part), arrays=int(part["in_matrix"].sum())) for correction, part in studies.groupby("correction")}
    C.finish(6, start, [C.DATA / "intensity_table_studies.parquet", C.DATA / "intensity_table_arrays.parquet", C.OUT / "intensity_tables_review.csv"] + outputs,
             method=policy["method"], study_platforms=dict(processed=len(studies), computed_in_this_run=computed),
             arrays=dict(routed=int(studies["samples"].sum()), in_matrix=len(arrays)), by_correction=by_correction,
             arrays_deposited_more_than_once=dict(arrays=len(repeated), identifiers=int(repeated["array_id"].nunique()),
                                                  with_the_same_median_intensity=len(same_export), studies=sorted(set(repeated["study_accession"])),
                                                  basis=f"same BeadChip barcode and position in the column label; medians within "
                                                        f"{policy['same_median_intensity']} (intensity_tables.same_median_intensity) count as the same"),
             offset=policy["offset"], estimator_without_detection=policy["estimator_without_detection"], detected_below=policy["detected_below"],
             columns=policy["columns"], tools=versions, protocol=C.CONFIG["protocol"])
    print(f"[06] Complete: {len(studies)} study platforms normalised with limma ({computed} computed in this run); {len(arrays):,} of "
          f"{int(studies['samples'].sum()):,} arrays are in a matrix; review sheet: {(C.OUT / 'intensity_tables_review.csv').relative_to(C.ROOT).as_posix()}",
          flush=True)
    for correction, item in by_correction.items():
        print(f"[06]   {correction}: {item['study_platforms']} study platforms, {item['arrays']:,} arrays", flush=True)
    print(f"[06]   {len(repeated):,} array records carry the BeadChip identifier of another record "
          f"({', '.join(sorted(set(repeated['study_accession']))) or 'none'}); {len(same_export):,} of them with the same median intensity", flush=True)


if __name__ == "__main__":
    main()
