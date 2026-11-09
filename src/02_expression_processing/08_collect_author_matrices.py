"""Collect, as deposited, the expression values of every study platform whose deposit holds only values that its authors had already normalised."""
import importlib
import io
import time

C = importlib.import_module("01_inventory_expression_files")

SCRIPT = C.HERE / C.STAGES[8]
TABLES = ("table", "workbook_sheet")
CLOSING = "!series_matrix_table_end"
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "table", "label", "missing_values", "value_median"]
STUDY_COLUMNS = ["source", "accession", "platform", "family", "samples", "in_matrix", "data_source", "features", "values", "scale", "missing_values",
                 "value_minimum", "value_median", "value_maximum", "reason", "matrix"]


def matrix_header(item):
    """Return the column labels of the value table of a GEO series matrix and the number of lines before its body."""
    text, passed = io.TextIOWrapper(C.stored(lambda: C.opened(item))[0], encoding="utf-8", errors="replace"), 0
    for line in text:
        passed += bool(line.strip())
        if line.startswith("!series_matrix_table_begin"):
            return [cell.strip('"') for cell in next(text).rstrip("\r\n").split("\t")], passed + 1
    raise ValueError(f"{item.accession}: {item.file} holds no value table")


def leftover_files(collected):
    """The files of the matrix folders that no current route produces: neither an output of stages 04 to 07 nor the matrix or record of a collected
    study platform of this script."""
    produced = {C.project_path(name) for number in range(4, 8) for name in C.read_json(C.receipt_path(number))["outputs"]}
    produced |= {(C.MATRICES / f"{name}{suffix}").resolve() for name in collected for suffix in (".json", ".parquet")}
    return sorted((path for folder in (C.MATRICES, C.COUNTS) if folder.is_dir() for path in folder.iterdir() if path.resolve() not in produced),
                  key=lambda path: path.as_posix())


def source_columns(key, route, samples, labels, listing):
    """Find the one table that holds the values of a study platform and the column of every sample that it holds."""
    files = listing[(listing["accession"] == route.accession) & (listing["sample_accession"] == "")]
    if route.data_source == "series_matrix":
        # A study with several platforms has one series matrix for each; the one that names the samples of this platform is read.
        found = [(item, *matrix_header(item)) for item in files[files["format"] == "geo_series_matrix"].itertuples(index=False)]
        found = [(item, header, body) for item, header, body in found if set(samples) <= set(header)]
        C.require(len(found) == 1 and len(set(found[0][1])) == len(found[0][1]), f"{key}: not exactly one series matrix names every sample of the study platform")
        item, header, body = found[0]
        return [dict(sample_accession=sample, table=item.table, item=item, body=body, delimiter="tab", position=header.index(sample), label=sample)
                for sample in sorted(samples)]
    mine = labels[(labels["accession"] == route.accession) & labels["table"].isin(route.data_source.split(";")) & (labels["platform"] == route.platform)
                  & (labels["sample_accession"] != "")]
    C.require(mine["table"].nunique() == 1 and mine["sample_accession"].is_unique,
              f"{key}: its samples are not in one table with one column each; a study platform that cannot be read so is recorded in not_collected of "
              f"{C.CONFIG['register']}")
    found = files[(files["table"] == mine["table"].iloc[0]) & files["format"].isin(TABLES)]
    C.require(len(found) == 1, f"{key}: the inventory does not list the table {mine['table'].iloc[0]} exactly once")
    item = next(found.itertuples(index=False))
    header, body, delimiter = C.table_labels(item)
    C.require(all(row.position < len(header) and header[row.position] == row.label for row in mine.itertuples()),
              f"{key}: the columns of {item.table} are not where the label map recorded them")
    return [dict(sample_accession=row.sample_accession, table=item.table, item=item, body=body, delimiter=delimiter, position=int(row.position), label=row.label)
            for row in mine.sort_values("sample_accession").itertuples()]


def main():
    import numpy as np
    import pandas as pd

    print(f"[08] Collecting author matrices v{C.VERSION}", flush=True)
    start = C.begin(8)
    policy, recorded, skipped = C.CONFIG["author_matrices"], C.REGISTER["author_matrices"], C.REGISTER["not_collected"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    routes = routes[routes["tier"] == "author_normalized"].sort_values(["samples", "accession", "platform"])
    keys = list(routes["accession"] + "|" + routes["platform"])
    C.require(sorted(keys) == sorted([*recorded, *skipped]),
              f"{C.CONFIG['register']} must record every author-normalized study platform once, under author_matrices or under not_collected: "
              f"{sorted(set(keys) ^ set([*recorded, *skipped]))}")
    held = pd.read_parquet(C.DATA / "sample_routes.parquet")
    assigned = pd.read_parquet(C.DATA / "table_samples.parquet")
    labels = pd.read_parquet(C.DATA / "label_map.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing["table"] = ["/".join(name for name in names if name) for names in zip(listing["file"], listing["member"], listing["sheet"])]
    archives = pd.read_parquet(C.ACQUIRED / "download_files.parquet")
    archives = {(row.accession, row.file): row.sha256 for row in archives[archives["status"] == "downloaded"].itertuples(index=False)}
    C.MATRICES.mkdir(parents=True, exist_ok=True)
    method = dict(script=C.sha(SCRIPT), reader=C.sha(C.HERE / C.STAGES[1]), inventory=C.CONFIG["inventory"],
                  largest_log2_like_value=policy["largest_log2_like_value"], runtime=C.runtime())
    # The matrices of a route that was removed or renamed would stay beside the current ones, so every file that no current route produces is removed.
    removed = leftover_files([f"{route.accession}__{route.platform}" for route, key in zip(routes.itertuples(index=False), keys) if key in recorded])
    for path in removed:
        path.unlink()
    if removed:
        print(f"[08] removed {len(removed)} files of the matrix folders that no current route produces", flush=True)

    studies, samples, outputs, computed = [], [], [], 0
    for index, (route, key) in enumerate(zip(routes.itertuples(index=False), keys), 1):
        name = f"{route.accession}__{route.platform}"
        record_path, matrix_path, began = C.MATRICES / f"{name}.json", C.MATRICES / f"{name}.parquet", time.monotonic()
        base = dict(source=route.source, accession=route.accession, platform=route.platform, family=route.family, samples=route.samples,
                    data_source="series_matrix" if route.data_source == "series_matrix" else "table")
        if key in skipped:
            # The deposit cannot be read into one matrix of comparable values; the study platform is kept with the recorded reason and no matrix.
            C.require(not matrix_path.exists() and not record_path.exists(), f"{key}: a matrix is stored for a study platform recorded as not collected")
            studies.append(dict(base, in_matrix=0, reason="not_collected", matrix=""))
            print(f"[08] {index}/{len(routes)} {route.accession} {route.platform}: not collected, as recorded in not_collected of {C.CONFIG['register']}",
                  flush=True)
            continue
        if route.data_source == "series_matrix":
            group = held[(held["study_accession"] == route.accession) & (held["platform"] == route.platform) & held["in_source"]]
        else:
            group = assigned[(assigned["study_accession"] == route.accession) & (assigned["platform"] == route.platform) & assigned["in_table"]]
        columns, entry = source_columns(key, route, sorted(group["sample_accession"]), labels, listing), recorded[key]
        C.require(columns and [column["sample_accession"] for column in columns] == sorted(group["sample_accession"]),
                  f"{key}: the samples with a value column are not the samples that the data source of the route holds")
        # A stored matrix stays valid while this script, the inventory script that reads the tables and its inventory settings, the Python packages, the
        # recorded scale, the file and the columns behind it are the same.
        identity = C.digest(dict(method=method, values=entry, table=[columns[0]["table"], archives[(route.accession, columns[0]["item"].file)]],
                                 columns=[[column["sample_accession"], column["position"]] for column in columns]))
        record = C.kept_record(record_path, matrix_path, identity)
        fresh = record is None
        if fresh:
            first, positions = columns[0], [column["position"] for column in columns]
            frame = C.table_numbers(first["item"], first["body"], first["delimiter"], sorted(positions))
            frame = frame[frame["identifier"] != CLOSING]
            C.require(len(frame) > 0 and (frame["identifier"] != "").all() and frame["identifier"].is_unique,
                      f"{key}: a feature identifier of {first['table']} is empty or repeated")
            values = frame[positions].to_numpy(dtype="float64")
            C.require(bool(np.isfinite(values).any(axis=0).all()), f"{key}: a sample has no value at all")
            # The recorded scale is checked against the values: a log2-like matrix holds no large value, a linear one no negative value.
            C.require(entry["scale"] != "log2_like" or float(np.nanmax(values)) <= policy["largest_log2_like_value"],
                      f"{key}: recorded as log2_like but holds values up to {np.nanmax(values):.1f}")
            C.require(entry["scale"] != "linear" or float(np.nanmin(values)) >= 0, f"{key}: recorded as linear but holds values below zero")
            stored = pd.DataFrame(values, columns=[column["sample_accession"] for column in columns])
            stored.insert(0, "feature_id", frame["identifier"].tolist())
            C.save_table(stored, matrix_path)
            record = dict(identity=identity, accession=route.accession, platform=route.platform, data_source=base["data_source"], table=first["table"],
                          values=entry["values"], scale=entry["scale"], basis=entry["basis"], features=len(frame),
                          missing_values=int(np.isnan(values).sum()), value_minimum=float(np.nanmin(values)), value_median=float(np.nanmedian(values)),
                          value_maximum=float(np.nanmax(values)),
                          samples=[dict(sample_accession=column["sample_accession"], table=column["table"], label=column["label"],
                                        missing_values=int(np.isnan(values[:, position]).sum()), value_median=float(np.nanmedian(values[:, position])))
                                   for position, column in enumerate(columns)],
                          matrix=C.record(matrix_path))
            C.write_json(record_path, record)
            computed += 1
        C.require([item["sample_accession"] for item in record["samples"]] == [column["sample_accession"] for column in columns],
                  f"{key}: the stored record does not list the samples of the route")
        studies.append(dict(base, in_matrix=len(record["samples"]), reason="", matrix=matrix_path.relative_to(C.ROOT).as_posix(),
                            **{field: record[field] for field in ("features", "values", "scale", "missing_values", "value_minimum", "value_median", "value_maximum")}))
        samples += [dict(source=route.source, study_accession=route.accession, platform=route.platform, **item) for item in record["samples"]]
        outputs += [matrix_path, record_path]
        print(f"[08] {index}/{len(routes)} {route.accession} {route.platform}: {len(record['samples']):,} of {route.samples:,} samples, {record['features']:,} features, "
              f"{record['values']} ({record['scale']}), {'collected' if fresh else 'kept'} in {time.monotonic() - began:,.0f} s", flush=True)

    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)
    samples = pd.DataFrame(samples, columns=SAMPLE_COLUMNS)
    C.require(len(studies) == len(routes) and len(samples) == studies["in_matrix"].sum() and (studies["in_matrix"] <= studies["samples"]).all()
              and not samples.duplicated(["study_accession", "sample_accession"]).any(), "Collected samples do not reconcile with the author-normalized routes")
    C.save_table(studies, C.DATA / "author_matrix_studies.parquet")
    C.save_table(samples, C.DATA / "author_matrix_samples.parquet")
    C.save_review(studies, C.OUT / "author_matrices_review.csv")
    collected = studies[studies["reason"] == ""]
    by_scale = {str(scale): dict(study_platforms=len(part), samples=int(part["in_matrix"].sum())) for scale, part in collected.groupby("scale")}
    C.finish(8, start, [C.DATA / "author_matrix_studies.parquet", C.DATA / "author_matrix_samples.parquet", C.OUT / "author_matrices_review.csv"] + outputs,
             method=policy["method"], study_platforms=dict(processed=len(studies), collected=len(collected), not_collected=len(studies) - len(collected),
                                                           computed_in_this_run=computed),
             samples=dict(routed=int(studies["samples"].sum()), in_matrix=len(samples)), by_scale=by_scale,
             largest_log2_like_value=policy["largest_log2_like_value"], values=recorded, not_collected=skipped,
             removed_leftover_files=[path.relative_to(C.ROOT).as_posix() for path in removed], protocol=C.CONFIG["protocol"])
    print(f"[08] Complete: {len(collected)} of {len(studies)} study platforms collected as deposited ({computed} read in this run); {len(samples):,} of "
          f"{int(studies['samples'].sum()):,} samples are in a matrix; review sheet: {(C.OUT / 'author_matrices_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for scale, item in by_scale.items():
        print(f"[08]   {scale}: {item['study_platforms']} study platforms, {item['samples']:,} samples", flush=True)
    print(f"[08]   not collected: {', '.join(sorted(skipped)) or 'none'}", flush=True)


if __name__ == "__main__":
    main()
