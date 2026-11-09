"""Collect the RNA-seq read counts of every study platform as deposited and compute TMM-normalised log2 counts per million, each study platform on its own."""
import hashlib
import importlib
import re
import shutil
import subprocess
import sys
import time

C = importlib.import_module("01_inventory_expression_files")

WORK = C.WORK / "count_tables"
SCRIPT = C.HERE / C.STAGES[7]
TOOLS = ["edgeR", "limma"]
TABLES = ("table", "workbook_sheet")
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "normalization_group", "table", "label", "columns", "library_size",
                  "genes_detected", "norm_factor", "log2_cpm_of_zero", "counts_digest"]
STUDY_COLUMNS = ["source", "accession", "platform", "samples", "in_matrix", "data_source", "tables", "genes", "identifiers_in_any_table", "gene_rule",
                 "summary_rows", "count_kind", "library_size_minimum", "library_size_median", "library_size_maximum", "counts", "matrix"]
EDGER = r"""
suppressPackageStartupMessages({library(jsonlite); library(edgeR)})
job <- fromJSON(commandArgs(TRUE)[1])
counts <- matrix(readBin(file.path(job$work, "counts.bin"), "double", job$genes * job$samples, size = 8), job$genes, job$samples)

# Step 1: TMM normalization factors from the counts of all genes; the library size of a sample is the sum of its counts.
y <- calcNormFactors(DGEList(counts = counts), method = "TMM")
# Step 2: log2 counts per million with the library sizes scaled by their TMM factors and the prior count of the job.
result <- cpm(y, log = TRUE, prior.count = job$prior_count)

writeBin(as.vector(result), file.path(job$work, "values.bin"), size = 8)
writeBin(as.vector(y$samples$norm.factors), file.path(job$work, "factors.bin"), size = 8)
writeBin(as.vector(y$samples$lib.size), file.path(job$work, "sizes.bin"), size = 8)
"""


def measure(label):
    """Return the measure a column label ends with, without case and punctuation; a label that is a name alone has none."""
    return "".join(re.findall(r"[a-z0-9]+", label[len(C.SUFFIX.sub("", label)):].lower()))


def source_columns(key, route, samples, labels, listing):
    """Name, for every sample that the data source of the route holds, the table its counts are in and, where the label map knows them, its columns."""
    if route.data_source == "per_sample_files":
        own = listing[(listing["accession"] == route.accession) & listing["sample_accession"].isin(samples) & listing["format"].isin(TABLES)
                      & (listing["value_kind"] != "none")]
        C.require(sorted(own["sample_accession"]) == sorted(samples), f"{key}: a sample does not have exactly one table of its own")
        return [dict(sample_accession=item.sample_accession, table=item.table, item=item, positions=None)
                for item in own.sort_values("sample_accession").itertuples(index=False)]
    tables = listing[(listing["accession"] == route.accession) & (listing["sample_accession"] == "") & listing["format"].isin(TABLES)]
    if route.data_source == "ncbi_counts":
        found = tables[tables["kind"] == "ncbi_generated_counts"]
        C.require(len(found) == 1, f"{key}: the study does not hold exactly one table of NCBI-generated counts")
        item = next(found.itertuples(index=False))
        header = C.table_labels(item)[0]
        C.require(len(set(header)) == len(header), f"{key}: a column label of {item.table} is repeated")
        return [dict(sample_accession=sample, table=item.table, item=item, positions=[header.index(sample)]) for sample in sorted(samples) if sample in header]
    mine = labels[(labels["accession"] == route.accession) & labels["table"].isin(route.data_source.split(";")) & (labels["platform"] == route.platform)
                  & (labels["sample_accession"] != "")]
    columns = []
    for sample, part in mine.groupby("sample_accession"):
        C.require(part["table"].nunique() == 1, f"{key}: sample {sample} is named in two tables")
        # Several columns of one sample hold one quantity measured again only where the label map marks them as repeats of one label and one measure.
        C.require(part["label"].nunique() == part["measure"].nunique() == 1 and (len(part) == 1 or bool(part["repeated"].all())),
                  f"{key}: sample {sample} is named by columns of different labels or measures in {part['table'].iloc[0]}: {part['label'].tolist()}")
        found = tables[tables["table"] == part["table"].iloc[0]]
        C.require(len(found) == 1, f"{key}: the inventory does not list the table {part['table'].iloc[0]} exactly once")
        columns.append(dict(sample_accession=sample, table=part["table"].iloc[0], item=next(found.itertuples(index=False)),
                            positions=sorted(int(position) for position in part["position"]), label=part["label"].iloc[0]))
    return columns


def settle_columns(key, columns, decision, policy):
    """Read the header of every table of one study platform once and settle its identifier column and the count columns of each of its samples."""
    headers = {}
    for column in columns:
        item = column["item"]
        if column["table"] not in headers:
            header, body, delimiter = C.table_labels(item)
            identifier = header.index(decision["identifier_column"]) if "identifier_column" in decision else 0
            headers[column["table"]] = dict(header=header, body=body, delimiter=delimiter, identifier=identifier)
        column.update(headers[column["table"]])
        header, identifier = column["header"], column["identifier"]
        if column["positions"] is None:
            # A file of one sample holds its counts in its second column, or in the one column that its label names as counts.
            named = [position for position, label in enumerate(header) if position != identifier and measure(label) and measure(label) in policy["count_measures"]]
            C.require(item.columns == 2 or len(named) == 1, f"{key}: {item.table} does not hold one count column")
            column["positions"] = [1] if item.columns == 2 else named
        C.require(not header or all(position < len(header) for position in column["positions"]),
                  f"{key}: the columns of {item.table} are not where the label map recorded them")
        column.setdefault("label", header[column["positions"][0]] if header else "")


def read_table(key, group, decision, policy):
    """Read the count columns of one table: its gene identifiers, the counts of each of its samples and the summary rows that were left out."""
    first = group[0]
    item = first["item"]
    positions = sorted({position for column in group for position in column["positions"]})
    frame = C.table_numbers(item, first["body"], first["delimiter"], positions, identifier=first["identifier"], decimal=decision.get("decimal", "."))
    # A table without a header row starts with its first gene; a first identifier that is a column label is a header whose labels are numbers.
    opening = frame["identifier"].iloc[0] if len(frame) else ""
    C.require(first["header"] or C.plain(opening) not in {C.plain(word) for word in policy["header_identifiers"]},
              f"{key}: {item.table} has no header row, yet its first identifier {opening} is a column label (count_tables.header_identifiers)")
    for pattern, replacement in decision.get("identifier_replace", []):
        frame["identifier"] = frame["identifier"].str.replace(pattern, replacement, regex=True)
    summary = frame["identifier"].str.fullmatch(policy["summary_rows"])
    left_out, frame = frame.loc[summary, "identifier"].tolist(), frame[~summary]
    # A missing count stops the run before rows are summed, because a sum would pass over it.
    C.require(bool(frame[positions].notna().all().all()), f"{key}: {item.table} holds a count that is missing")
    if decision.get("repeated_identifiers") == "sum":
        # Counts add up, so rows that carry one identifier are summed where a decision says so; the genes keep the order of their first rows.
        frame = frame.groupby("identifier", sort=False, as_index=False)[positions].sum()
    values, genes = frame[positions].to_numpy(dtype="float64"), frame["identifier"]
    C.require(len(genes) > 0 and (genes != "").all() and genes.is_unique, f"{key}: a gene identifier of {item.table} is empty or repeated")
    C.require(bool((values >= 0).all()), f"{key}: {item.table} holds a count that is missing or below zero")
    place = {position: index for index, position in enumerate(positions)}
    counts = {column["sample_accession"]: values[:, [place[position] for position in column["positions"]]].sum(axis=1) for column in group}
    # Many files of one study repeat the same identifiers, so each distinct identifier is kept in memory once.
    return [sys.intern(gene) for gene in genes], counts, left_out


def read_counts(key, columns, decision, policy):
    """Read the counts of every sample of one study platform and put them on one list of gene identifiers."""
    import numpy as np
    import pandas as pd

    tables, read, summary = {}, [], set()
    for column in columns:
        tables.setdefault(column["table"], []).append(column)
    for group in tables.values():
        C.require(all(len(column["positions"]) == 1 for column in group) or decision.get("repeated_columns") == "sum",
                  f"{key}: a sample has several count columns in {group[0]['table']} and count_tables of {C.CONFIG['register']} does not say how to read them")
        identifiers, counts, left_out = read_table(key, group, decision, policy)
        read.append((read[-1][0] if read and identifiers == read[-1][0] else identifiers, counts))
        summary.update(left_out)
    genes, rule, union = read[0][0], "same_genes", set().union(*(identifiers for identifiers, _ in read))
    if any(identifiers is not genes and set(identifiers) != set(genes) for identifiers, _ in read):
        # Tables of one study platform may list other genes only where a recorded decision says how to read that.
        rule = decision.get("genes")
        C.require(rule in ("shared", "absent_is_zero"),
                  f"{key}: its tables do not hold the same genes and count_tables of {C.CONFIG['register']} does not say how to read them")
        if rule == "shared":
            shared = set(genes).intersection(*(identifiers for identifiers, _ in read))
            genes = [gene for gene in genes if gene in shared]
        else:
            # A gene that a file does not list is zero only if the files leave out genes without reads, so none of them may list a zero.
            C.require(all(bool((values > 0).all()) for _, counts in read for values in counts.values()),
                      f"{key}: a table lists a count of zero, so a gene that another table does not list is not known to be zero")
            genes = list(dict.fromkeys(gene for identifiers, _ in read for gene in identifiers))
    C.require(len(genes) > 0, f"{key}: no gene is left")
    vectors = {}
    for identifiers, counts in read:
        for sample, values in counts.items():
            # A gene that a table does not list has no reads there; under the rule "shared" every gene of the matrix is in every table.
            same = identifiers is genes or identifiers == genes
            vectors[sample] = values if same else pd.Series(values, index=identifiers).reindex(genes).fillna(0).to_numpy()
    return genes, np.column_stack([vectors[column["sample_accession"]] for column in columns]), sorted(summary), len(union), rule


def run_edger(name, counts, prior_count):
    """Compute TMM factors and log2 counts per million for one study platform with edgeR."""
    import numpy as np

    work = WORK / name
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        # R reads a matrix column by column, one sample after the other.
        np.ascontiguousarray(counts.T).tofile(work / "counts.bin")
        C.write_json(work / "job.json", dict(work=work.as_posix(), genes=counts.shape[0], samples=counts.shape[1], prior_count=prior_count))
        (work / "edger.R").write_text(EDGER, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (work / "edger.R").as_posix(), (work / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"{name}: edgeR failed:\n{done.stderr[-3000:]}")
        values, factors, sizes = (np.fromfile(work / f"{part}.bin", dtype="<f8") for part in ("values", "factors", "sizes"))
        C.require(values.size == counts.size and factors.size == sizes.size == counts.shape[1], f"{name}: the edgeR output is incomplete")
        return values.reshape(counts.shape[1], counts.shape[0]).T, factors, sizes
    finally:
        shutil.rmtree(work, ignore_errors=True)


def normalise(name, counts, prior_count, groups):
    """TMM factors and log2 counts per million of one study platform, computed over all of its samples or within each of its normalization groups."""
    import numpy as np

    if groups is None:
        return run_edger(name, counts, prior_count)
    values, factors, sizes = np.full(counts.shape, np.nan), np.full(counts.shape[1], np.nan), np.full(counts.shape[1], np.nan)
    for group in sorted(set(groups)):
        chosen = np.flatnonzero(np.array(groups) == group)
        values[:, chosen], factors[chosen], sizes[chosen] = run_edger(f"{name}__{group}", counts[:, chosen], prior_count)
    return values, factors, sizes


def zero_floor(sizes, factors, groups, prior_count):
    """The log2 count per million that a gene without reads receives in each sample: edgeR scales the prior count to the mean TMM-scaled library size of
    the samples normalized together."""
    import numpy as np

    effective = np.asarray(sizes) * np.asarray(factors)
    labels = np.array(groups if groups is not None else [""] * len(effective))
    mean = np.array([effective[labels == label].mean() for label in labels])
    scaled = prior_count * effective / mean
    return np.log2(scaled / (effective + 2 * scaled) * 1e6)


def main():
    import numpy as np
    import pandas as pd

    print(f"[07] Processing count tables v{C.VERSION}", flush=True)
    start = C.begin(7)
    policy, decisions = C.CONFIG["count_tables"], C.REGISTER["count_tables"]
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    routes = routes[routes["tier"] == "counts"].sort_values(["samples", "accession", "platform"])
    unknown = sorted(set(decisions) - set(routes["accession"] + "|" + routes["platform"]))
    C.require(not unknown, f"Count decisions recorded for study platforms that are not count routes: {unknown}")
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
    C.COUNTS.mkdir(parents=True, exist_ok=True)
    method = dict(script=C.sha(SCRIPT), reader=C.sha(C.HERE / C.STAGES[1]), inventory=C.CONFIG["inventory"], prior_count=policy["prior_count"],
                  count_measures=policy["count_measures"], summary_rows=policy["summary_rows"], header_identifiers=policy["header_identifiers"],
                  tools=versions, runtime=C.runtime())

    studies, samples, outputs, computed = [], [], [], 0
    for index, route in enumerate(routes.itertuples(index=False), 1):
        key, decision = f"{route.accession}|{route.platform}", decisions.get(f"{route.accession}|{route.platform}", {})
        if route.data_source in ("per_sample_files", "ncbi_counts"):
            group = held[(held["study_accession"] == route.accession) & (held["platform"] == route.platform) & held["in_source"]]
        else:
            group = assigned[(assigned["study_accession"] == route.accession) & (assigned["platform"] == route.platform) & assigned["in_table"]]
        columns = source_columns(key, route, sorted(group["sample_accession"]), labels, listing)
        C.require(columns and [column["sample_accession"] for column in columns] == sorted(group["sample_accession"]),
                  f"{key}: the samples with a count column are not the samples that the data source of the route holds")
        settle_columns(key, columns, decision, policy)
        found = C.normalization_groups(key, [column["sample_accession"] for column in columns])
        groups = [found[column["sample_accession"]] for column in columns] if found else None
        # Stored matrices stay valid while this script, the inventory script that reads the tables and its inventory settings, the software, the recorded
        # decision, the files, the columns read from them and the normalization groups are the same.
        identity = C.digest(dict(method=method, decision=decision, groups=groups, samples=[[column["sample_accession"], column["table"], column["positions"],
                                                                                            archives[(route.accession, column["item"].file)]] for column in columns]))
        name = f"{route.accession}__{route.platform}"
        record_path, matrix_path, counts_path = C.MATRICES / f"{name}.json", C.MATRICES / f"{name}.parquet", C.COUNTS / f"{name}.parquet"
        began, record = time.monotonic(), C.kept_record(record_path, matrix_path, identity)
        if record is not None and not (counts_path.is_file() and record.get("counts") == C.record(counts_path)):
            record = None
        fresh = record is None
        if fresh:
            genes, counts, summary, union, rule = read_counts(key, columns, decision, policy)
            sizes = counts.sum(axis=0)
            C.require(bool((sizes > 0).all()), f"{key}: a sample has no counts at all")
            values, factors, used = normalise(name, counts, policy["prior_count"], groups)
            floors = zero_floor(sizes, factors, groups, policy["prior_count"])
            C.require(bool(np.isfinite(values).all()) and bool(np.isfinite(factors).all()) and bool(np.allclose(used, sizes)),
                      f"{key}: edgeR returned values that are not finite, or other library sizes")
            whole = bool((counts == np.round(counts)).all())
            titles = [column["sample_accession"] for column in columns]
            stored = pd.DataFrame(counts.astype("int32") if whole and counts.max() < 2 ** 31 else counts, columns=titles)
            stored.insert(0, "feature_id", genes)
            C.save_table(stored, counts_path)
            frame = pd.DataFrame(values.astype("float32"), columns=titles)
            frame.insert(0, "feature_id", genes)
            C.save_table(frame, matrix_path)
            # Counts are hashed in the order of their sorted identifiers, so that one sample deposited in two studies has one digest.
            order = np.argsort(np.array(genes, dtype=object), kind="stable")
            record = dict(identity=identity, accession=route.accession, platform=route.platform, data_source=route.data_source, decision=decision,
                          tables=list(dict.fromkeys(column["table"] for column in columns)), genes=len(genes), identifiers_in_any_table=union,
                          gene_rule=rule, summary_rows=summary,
                          count_kind="whole_numbers" if whole else "decimals", prior_count=policy["prior_count"],
                          value_minimum=float(values.min()), value_maximum=float(values.max()),
                          samples=[dict(sample_accession=column["sample_accession"], normalization_group=groups[position] if groups else "",
                                        table=column["table"], label=column["label"], columns=len(column["positions"]), library_size=float(sizes[position]),
                                        genes_detected=int((counts[:, position] > 0).sum()), norm_factor=float(factors[position]),
                                        log2_cpm_of_zero=float(floors[position]),
                                        counts_digest=hashlib.sha256(np.ascontiguousarray(counts[order, position]).tobytes()).hexdigest())
                                   for position, column in enumerate(columns)],
                          tools=versions, counts=C.record(counts_path), matrix=C.record(matrix_path))
            C.write_json(record_path, record)
            computed += 1
        C.require([item["sample_accession"] for item in record["samples"]] == [column["sample_accession"] for column in columns],
                  f"{key}: the stored record does not list the samples of the route")
        sizes = [item["library_size"] for item in record["samples"]]
        studies.append(dict(source=route.source, accession=route.accession, platform=route.platform, samples=route.samples, in_matrix=len(record["samples"]),
                            data_source=route.data_source if route.data_source in ("per_sample_files", "ncbi_counts") else "table",
                            tables=len(record["tables"]), genes=record["genes"], identifiers_in_any_table=record["identifiers_in_any_table"],
                            gene_rule=record["gene_rule"], summary_rows=len(record["summary_rows"]), count_kind=record["count_kind"],
                            library_size_minimum=min(sizes), library_size_median=float(np.median(sizes)), library_size_maximum=max(sizes),
                            counts=counts_path.relative_to(C.ROOT).as_posix(), matrix=matrix_path.relative_to(C.ROOT).as_posix()))
        samples += [dict(source=route.source, study_accession=route.accession, platform=route.platform, **item) for item in record["samples"]]
        outputs += [matrix_path, counts_path, record_path]
        print(f"[07] {index}/{len(routes)} {route.accession} {route.platform or '(no platform)'}: {len(record['samples']):,} of {route.samples:,} samples, "
              f"{record['genes']:,} genes, {record['count_kind']}, {'computed' if fresh else 'kept'} in {time.monotonic() - began:,.0f} s", flush=True)

    C.clear_work(WORK)
    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS).sort_values(["source", "accession", "platform"]).reset_index(drop=True)
    samples = pd.DataFrame(samples, columns=SAMPLE_COLUMNS)
    C.require(len(studies) == len(routes) and len(samples) == studies["in_matrix"].sum() and (studies["in_matrix"] <= studies["samples"]).all()
              and not samples.duplicated(["study_accession", "sample_accession"]).any(), "Processed samples do not reconcile with the count routes")
    repeated = samples[samples.duplicated("counts_digest", keep=False)]
    C.save_table(studies, C.DATA / "count_table_studies.parquet")
    C.save_table(samples, C.DATA / "count_table_samples.parquet")
    C.save_review(studies, C.OUT / "count_tables_review.csv")
    by_source = {str(source): dict(study_platforms=len(part), samples=int(part["samples"].sum()), in_matrix=int(part["in_matrix"].sum()))
                 for source, part in studies.groupby("data_source")}
    C.finish(7, start, [C.DATA / "count_table_studies.parquet", C.DATA / "count_table_samples.parquet", C.OUT / "count_tables_review.csv"] + outputs,
             method=policy["method"], study_platforms=dict(processed=len(studies), computed_in_this_run=computed),
             samples=dict(routed=int(studies["samples"].sum()), in_matrix=len(samples)), by_data_source=by_source,
             counts_that_are_not_whole_numbers=sorted(studies.loc[studies["count_kind"] == "decimals", "accession"] + "|"
                                                      + studies.loc[studies["count_kind"] == "decimals", "platform"]),
             samples_deposited_more_than_once=dict(samples=len(repeated), digests=int(repeated["counts_digest"].nunique()),
                                                   studies=sorted(set(repeated["study_accession"])), basis="the same counts for the same gene identifiers"),
             prior_count=policy["prior_count"], count_measures=policy["count_measures"], summary_rows=policy["summary_rows"],
             decisions=decisions, tools=versions, protocol=C.CONFIG["protocol"],
             normalization_groups={key: entry for key, entry in C.REGISTER["normalization_groups"].items() if key in set(routes["accession"] + "|" + routes["platform"])},
             log2_cpm_of_zero=dict(lowest=round(float(samples["log2_cpm_of_zero"].min()), 4), median=round(float(samples["log2_cpm_of_zero"].median()), 4),
                                   highest=round(float(samples["log2_cpm_of_zero"].max()), 4)))
    print(f"[07] Complete: {len(studies)} study platforms ({computed} computed in this run); {len(samples):,} of {int(studies['samples'].sum()):,} samples "
          f"have counts and TMM-normalised log2 counts per million; review sheet: {(C.OUT / 'count_tables_review.csv').relative_to(C.ROOT).as_posix()}",
          flush=True)
    for source, item in by_source.items():
        print(f"[07]   {source}: {item['study_platforms']} study platforms, {item['in_matrix']:,} of {item['samples']:,} samples", flush=True)
    print(f"[07]   {len(repeated):,} sample records carry the counts of another record ({', '.join(sorted(set(repeated['study_accession']))) or 'none'})",
          flush=True)


if __name__ == "__main__":
    main()
