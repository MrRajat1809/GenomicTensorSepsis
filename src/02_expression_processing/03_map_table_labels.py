"""Map the columns of every routed table to the samples of its study platforms and record the rule behind each match."""
import collections
import importlib
import io
import re

C = importlib.import_module("01_inventory_expression_files")

TABLES = ("table", "workbook_sheet")
KEYWORDS = ("per_sample_files", "ncbi_counts", "series_matrix", "none")
FIELDS = ("sample_accession", "title", "description", "source_name")
MAP_COLUMNS = ["source", "accession", "table", "position", "label", "measure", "sample_accession", "platform", "rule", "field", "rewritten",
               "repeated", "correlation", "margin"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "tier", "table", "columns", "rule", "in_table", "reason"]
REVIEW_COLUMNS = ["source", "accession", "table", "platforms", "tier", "columns", "sample_columns", "samples_named", "rules", "measures", "rewritten",
                  "repeated_samples", "unnamed_labels", "first_unnamed", "ambiguous_labels", "route_samples", "route_samples_unnamed", "flag", "basis"]


def parts(value):
    """Split an identifier into its runs of letters and digits, without case."""
    return re.findall(r"[a-z0-9]+", str(value).lower())


def digit_runs(found):
    """The runs of digits of an identifier split into parts: P1_1 has two, P11 has one."""
    return tuple(re.findall(r"\d+", " ".join(found)))


def identifiers(own, values, limit):
    """Index the candidate samples by every whole field of their records and by every run of up to `limit` parts of a field, with the runs of digits
    of each and the first field that holds them."""
    whole, part = collections.defaultdict(dict), collections.defaultdict(dict)
    for row in own.itertuples(index=False):
        fields = [(name, getattr(row, name)) for name in FIELDS] + [("characteristic", value) for value in values.get(row.sample_accession, ())]
        fields += [("file_name", name.rsplit("/", 1)[-1]) for name in row.raw_files.split(";") if name]
        for field, value in fields:
            found = parts(value)
            if found:
                whole["".join(found)].setdefault(row.sample_accession, {}).setdefault(digit_runs(found), field)
            for size in range(1, min(limit, len(found) - 1) + 1):
                for first in range(len(found) - size + 1):
                    window = found[first:first + size]
                    part["".join(window)].setdefault(row.sample_accession, {}).setdefault(digit_runs(window), field)
    return whole, part


def name_sample(text, accessions, whole, part):
    """Name the sample a label carries, with the rule and field behind it; a label that would name several samples names none."""
    base = C.SUFFIX.sub("", text)
    measure, found = "".join(parts(text[len(base):])), C.SAMPLE.search(text)
    if found and found.group(1) in accessions:
        return dict(sample=found.group(1), rule="accession", field="sample_accession", measure=measure)
    if not parts(base):
        return dict(sample="", rule="measure_only" if measure else "", field="", measure=measure)
    for rule, index in (("whole_field", whole), ("part_of_field", part)):
        for pieces, kind in ((parts(base), measure), (parts(text), "")):
            key = "".join(pieces)
            # A bare number is never taken as an identifier, because ages and day counts are numbers too.
            samples = index.get(key, {}) if re.search("[a-z]", key) else {}
            if len(samples) == 1:
                (sample, fields), = samples.items()
                # Separators are ignored except where they part two runs of digits, so that P1_1 does not name P11.
                if digit_runs(pieces) in fields:
                    return dict(sample=sample, rule=rule, field=fields[digit_runs(pieces)], measure=kind)
            elif samples:
                return dict(sample="", rule="ambiguous", field="", measure="")
    return dict(sample="", rule="", field="", measure="")


def map_columns(labels, own, values, replace, limit):
    """Assign every column of one table to the sample its label names, or to none; a column labelled only by a measure follows its sample."""
    accessions, platform = set(own["sample_accession"]), dict(zip(own["sample_accession"], own["platform"]))
    whole, part = identifiers(own, values, limit)
    rows, previous = [], None
    for position, label in enumerate(labels):
        text = label.strip()
        for pattern, replacement in replace:
            text = re.sub(pattern, replacement, text)
        found = name_sample(text, accessions, whole, part)
        if found["rule"] == "measure_only":
            found = dict(previous, rule="preceding_column", measure=found["measure"]) if previous else dict(found, rule="")
        else:
            previous = found if found["sample"] else None
        rows.append(dict(position=position, label=label, measure=found["measure"], sample_accession=found["sample"],
                         platform=platform.get(found["sample"], ""), rule=found["rule"], field=found["field"],
                         rewritten=bool(found["sample"]) and text != label.strip(), repeated=False, correlation=None, margin=None))
    return rows


def table_values(item, body, delimiter, positions):
    """Read chosen columns of a text table as numbers, keyed by its first column."""
    import pandas as pd

    text, passed = io.TextIOWrapper(C.stored(lambda: C.opened(item))[0], encoding="utf-8-sig", errors="replace"), 0
    for line in text:
        # The body starts after the rows the preview counted, and the preview leaves out blank lines.
        passed += bool(line.strip())
        if passed == body:
            break
    frame = pd.read_csv(text, sep=C.SEPARATORS[delimiter], header=None, usecols=[0, *positions], index_col=0, dtype={0: str}, low_memory=False)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame[~frame.index.duplicated()]


def matrix_values(path):
    """Read the value table of a GEO series matrix; its columns are sample accessions."""
    import pandas as pd

    text = io.TextIOWrapper(C.stored(lambda: path.open("rb"))[0], encoding="utf-8", errors="replace")
    for line in text:
        if line.startswith("!series_matrix_table_begin"):
            break
    frame = pd.read_csv(text, sep="\t", index_col=0, comment="!", dtype={"ID_REF": str}, low_memory=False)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame[~frame.index.duplicated()]


def centred_ranks(frame):
    """Rank the rows within each column, then remove the mean rank of each row, so that what all arrays share does not count."""
    ranks = frame.rank().to_numpy(dtype="float64")
    ranks -= ranks.mean(axis=1, keepdims=True)
    return (ranks - ranks.mean(axis=0)) / ranks.std(axis=0)


def match_by_values(item, body, delimiter, rows, own, matrices, policy, name):
    """Name the unnamed columns of a table by the series-matrix sample that holds the same values."""
    import numpy as np

    positions = [row["position"] for row in rows[1:] if not row["sample_accession"]]
    wanted = sorted(set(own["sample_accession"]) - {row["sample_accession"] for row in rows})
    C.require(delimiter and positions and len(wanted) >= policy["minimum_columns"], f"{name}: no unnamed columns or samples are left to match by values")
    table = table_values(item, body, delimiter, positions)
    table = table.loc[:, table.notna().mean() > policy["minimum_present_share"]]
    reference = next((frame for frame in map(matrix_values, matrices) if set(wanted) <= set(frame.columns)), None)
    C.require(reference is not None, f"{name}: no series matrix holds every sample that is still unnamed")
    reference = reference[wanted]
    shared = table.index.intersection(reference.index)
    shared = shared[(table.loc[shared].notna().all(axis=1) & reference.loc[shared].notna().all(axis=1)).to_numpy()]
    C.require(len(shared) >= policy["minimum_shared_rows"] and table.shape[1] >= policy["minimum_columns"],
              f"{name}: only {len(shared):,} complete rows are shared with the series matrix")
    scores = centred_ranks(table.loc[shared]).T @ centred_ranks(reference.loc[shared]) / len(shared)
    forward, backward, platform = scores.argmax(axis=1), scores.argmax(axis=0), dict(zip(own["sample_accession"], own["platform"]))
    for index, position in enumerate(table.columns):
        best, runner = np.sort(scores[index])[-1], np.sort(scores[index])[-2]
        if backward[forward[index]] == index and best >= policy["minimum_correlation"] and best - runner >= policy["minimum_margin"]:
            sample = wanted[forward[index]]
            rows[position].update(sample_accession=sample, platform=platform[sample], rule="values", field="series_matrix",
                                  correlation=round(float(best), 6), margin=round(float(best - runner), 6))


def check_columns(rows, name):
    """Mark the samples that one label names in several columns; different labels may not name the same sample."""
    named = collections.defaultdict(list)
    for row in rows:
        if row["sample_accession"]:
            named[(row["sample_accession"], row["measure"])].append(row)
    for (sample, _), group in named.items():
        C.require(len({row["label"] for row in group}) == 1, f"{name}: different labels name sample {sample}: {[row['label'] for row in group]}")
        for row in group:
            row["repeated"] = len(group) > 1


def summarise(rows, decided):
    """Count what the columns of one table were assigned to, for the review sheet."""
    named = [row for row in rows if row["sample_accession"]]
    unnamed = [row["label"] for row in rows[1:] if not row["sample_accession"] and row["rule"] != "ambiguous"]
    rules = collections.Counter(":".join(name for name in (row["rule"], row["field"]) if name) for row in named if row["rule"] != "preceding_column")
    ambiguous = sum(row["rule"] == "ambiguous" for row in rows)
    repeated = len({row["sample_accession"] for row in named if row["repeated"]})
    return dict(columns=len(rows), sample_columns=len(named), samples_named=len({row["sample_accession"] for row in named}),
                rules="; ".join(f"{rule}={count}" for rule, count in sorted(rules.items())),
                measures=";".join(sorted({row["measure"] or "value" for row in named})), rewritten=sum(row["rewritten"] for row in named),
                repeated_samples=repeated, unnamed_labels=len(unnamed), first_unnamed=" | ".join(unnamed[:8])[:200], ambiguous_labels=ambiguous,
                flag="uncertain" if decided or ambiguous or repeated or len(rules) > 1 else "")


def main():
    import pandas as pd

    print(f"[03] Mapping table labels v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy, rewrites, by_values = C.CONFIG["labels"], C.REGISTER["label_rewrites"], C.REGISTER["labels_by_values"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    before = pd.read_parquet(C.DATA / "sample_routes.parquet")
    samples = pd.read_parquet(C.ACQUIRED / "samples.parquet")
    attributes = pd.read_parquet(C.ACQUIRED / "sample_attributes.parquet")
    routed = routes[~routes["data_source"].isin(KEYWORDS)]
    readers = collections.defaultdict(list)
    for route in routed.itertuples(index=False):
        for name in route.data_source.split(";"):
            readers[(route.source, route.accession, name)].append(route)
    decisions = {**rewrites, **by_values}
    unknown = sorted(set(decisions) - {f"{accession}|{name}" for _, accession, name in readers})
    C.require(not unknown, f"Label decisions recorded for tables that no route reads: {unknown}")
    C.require(not set(rewrites) & set(by_values), "A table is either rewritten or matched by values, not both")

    listing = listing[listing["format"].isin(TABLES + ("geo_series_matrix",)) & (listing["sample_accession"] == "")].copy()
    listing["table"] = ["/".join(name for name in names if name) for names in zip(listing["file"], listing["member"], listing["sheet"])]
    values = collections.defaultdict(lambda: collections.defaultdict(list))
    for study, sample, value in attributes[["study_accession", "sample_accession", "value"]].itertuples(index=False, name=None):
        values[study][sample].append(value)

    columns, tables = [], []
    for index, ((source, accession, name), users) in enumerate(sorted(readers.items()), 1):
        key, files = f"{accession}|{name}", listing[listing["accession"] == accession]
        found = files[(files["table"] == name) & files["format"].isin(TABLES)]
        C.require(len(found) == 1, f"{key}: the route names a table that the inventory does not list exactly once")
        item = next(found.itertuples(index=False))
        own = samples[(samples["study_accession"] == accession) & samples["platform"].isin([route.platform for route in users])]
        header, body, delimiter = C.table_labels(item)
        C.require(header, f"{key}: the table has no header row")
        rows = map_columns(header, own, values[accession], rewrites.get(key, {}).get("replace", []), policy["identifier_parts"])
        if key in by_values:
            matrices = [C.RAW / source / accession / file for file in files.loc[files["format"] == "geo_series_matrix", "file"]]
            match_by_values(item, body, delimiter, rows, own, matrices, policy["values"], key)
        check_columns(rows, key)
        columns += [dict(source=source, accession=accession, table=name, **row) for row in rows]
        tables.append(dict(source=source, accession=accession, table=name, platforms=";".join(route.platform for route in users),
                           tier=users[0].tier, basis=decisions.get(key, {}).get("basis", ""), **summarise(rows, key in decisions)))
        if index % 10 == 0 or index == len(readers):
            print(f"[03] mapped {index:,}/{len(readers):,} tables, {len(columns):,} columns", flush=True)

    named = collections.defaultdict(list)
    for row in columns:
        if row["sample_accession"]:
            named[(row["accession"], row["sample_accession"])].append(row)
    rows = []
    for route in routed.itertuples(index=False):
        group = samples[(samples["study_accession"] == route.accession) & (samples["platform"] == route.platform)]
        for sample in group["sample_accession"]:
            mine = [row for row in named[(route.accession, sample)] if row["table"] in route.data_source.split(";")]
            holding = sorted({row["table"] for row in mine})
            C.require(len(holding) < 2, f"{route.accession}: sample {sample} is named in two tables of its route: {holding}")
            rows.append(dict(source=route.source, study_accession=route.accession, sample_accession=sample, platform=route.platform, tier=route.tier,
                             table="".join(holding), columns=len(mine), rule=";".join(sorted({row["rule"] for row in mine} - {"preceding_column"})),
                             in_table=bool(mine), reason="" if mine else "absent_from_table"))
    assigned = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    C.require(len(assigned) == routed["samples"].sum() and not assigned.duplicated(["study_accession", "sample_accession"]).any(),
              "The label map does not hold every sample of the table routes exactly once")
    held = before[before["in_source"] & ~before["data_source"].isin(KEYWORDS)]
    lost = sorted(set(zip(held["study_accession"], held["sample_accession"]))
                  - set(zip(assigned.loc[assigned["in_table"], "study_accession"], assigned.loc[assigned["in_table"], "sample_accession"])))
    C.require(not lost, f"{len(lost)} samples that the route check found in their table are no longer named: {lost[:20]}")

    size = assigned.groupby(["study_accession", "platform"])["in_table"].agg(["size", "sum"])
    for table in tables:
        part = size.loc[[(table["accession"], platform) for platform in table["platforms"].split(";")]]
        table.update(route_samples=int(part["size"].sum()), route_samples_unnamed=int(part["size"].sum() - part["sum"].sum()))
    label_map = pd.DataFrame(columns, columns=MAP_COLUMNS).astype({"correlation": "float64", "margin": "float64"})
    review = pd.DataFrame(tables, columns=REVIEW_COLUMNS)
    C.save_table(label_map, C.DATA / "label_map.parquet")
    C.save_table(assigned, C.DATA / "table_samples.parquet")
    C.save_review(review, C.OUT / "labels_review.csv")
    by_rule = {str(rule): int(count) for rule, count in assigned.loc[assigned["in_table"], "rule"].value_counts().sort_index().items()}
    absent = {str(study): int(count) for study, count in assigned.loc[~assigned["in_table"], "study_accession"].value_counts().sort_index().items()}
    C.finish(3, start, [C.DATA / "label_map.parquet", C.DATA / "table_samples.parquet", C.OUT / "labels_review.csv"],
             method=policy["method"], tables=len(tables), study_platforms=len(routed),
             samples=dict(total=len(assigned), named=int(assigned["in_table"].sum()), by_rule=by_rule, absent_from_table=absent),
             columns=dict(total=len(label_map), naming_a_sample=int((label_map["sample_accession"] != "").sum()),
                          rewritten=int(label_map["rewritten"].sum()), repeated=int(label_map["repeated"].sum()),
                          ambiguous=int((label_map["rule"] == "ambiguous").sum())),
             flagged_uncertain=int((review["flag"] == "uncertain").sum()), identifier_parts=policy["identifier_parts"], values=policy["values"],
             rewrites=rewrites, by_values=by_values, protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: {len(tables)} tables of {len(routed)} study platforms; {int(assigned['in_table'].sum()):,} of {len(assigned):,} samples "
          f"named by a column; review sheet: {(C.OUT / 'labels_review.csv').relative_to(C.ROOT).as_posix()}",
          flush=True)
    for rule, count in by_rule.items():
        print(f"[03]   {rule}: {count:,} samples", flush=True)
    print(f"[03]   absent_from_table: {sum(absent.values()):,} samples in {len(absent)} studies: {absent}", flush=True)


if __name__ == "__main__":
    main()
