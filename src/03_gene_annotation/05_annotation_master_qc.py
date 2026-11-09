"""Verify the indexes, the annotation sources, the feature maps and the gene matrices of the annotation stage from their files."""
import importlib
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_index_expression_matrices")
G = importlib.import_module("04_build_gene_matrices")

STATUSES = ("one_gene", "several_genes", "no_gene_in_annotation", "gene_not_in_hgnc", "accession_not_resolved", "not_in_annotation")
NCBI_COLUMNS = ["accession", "platform", "ncbi_table", "status", "samples_compared", "own_column_most_similar", "shared_genes", "genes_used",
                "median_own_correlation", "smallest_margin", "not_most_similar"]
MEANING = ("PASS means that every study platform and every included sample is indexed once and agrees with the expression stage; that every "
           "annotation source is stored unchanged; that every row of every expression matrix is in its map with one status; that every gene-level "
           "value is the kept row of its matrix (the highest-mean row among the rows of the preferred probe quality where the gene has one), or the sum "
           "of its distinct rows, as the stage records, with the log2 counts per million and the expressed flag of every count study recomputed with "
           "its normalization groups; and that every count study read from an author table or from files of single samples, whose study also deposits "
           "NCBI-generated counts, has every sample most similar to the NCBI counts of its own accession, unless the configuration records it as an "
           "exception that cannot be compared. It does not assess whether an annotation "
           "assigns a probe to the right gene.")


class QC:
    def __init__(self):
        self.checks = []
        self.section = ""

    def check(self, name, passed, detail=None):
        self.checks.append(dict(section=self.section, name=name, passed=bool(passed), detail=detail))
        if not passed:
            raise ValueError(f"QC failed: {name}")

    def phase(self, name, action):
        self.section = name
        print(f"[05] {name}", flush=True)
        action()


def provenance_checks(qc, state):
    state["start"] = C.begin(5)
    qc.check("all_stage_receipts_and_code_bindings", True)


def index_checks(qc, state):
    recorded, expression = C.expression_binding()[1:]
    matrices = pd.read_parquet(C.DATA / "matrix_index.parquet")
    samples = pd.read_parquet(C.DATA / "sample_index.parquet")
    held = matrices[matrices["matrix"] != ""]
    qc.check("every_study_platform_indexed_once", not matrices.duplicated(["accession", "platform"]).any()
             and len(matrices) == expression["summary"]["study_platforms"] and len(held) == expression["summary"]["matrices"])
    qc.check("every_included_sample_indexed_once", not samples.duplicated(["study_accession", "sample_accession"]).any()
             and len(samples) == expression["summary"]["samples"])
    reasons = {str(reason): int(count) for reason, count in samples.loc[~samples["in_matrix"], "reason"].value_counts().sort_index().items()}
    qc.check("sample_outcomes_agree_with_the_expression_stage", int(samples["in_matrix"].sum()) == expression["summary"]["in_a_matrix"]
             and reasons == expression["summary"]["without_a_matrix"] and (samples["in_matrix"] == (samples["reason"] == "")).all())
    inside = samples[samples["in_matrix"]]
    qc.check("every_sample_in_a_matrix_names_its_matrix_and_its_source", (inside["matrix"] != "").all() and (inside["values_from"] != "").all()
             and set(inside["matrix"]) == set(held["matrix"]) and (inside["title"] != "").all())
    qc.check("every_indexed_matrix_is_an_output_of_the_expression_stage", set(held["matrix"]) <= set(recorded)
             and set(held.loc[held["counts"] != "", "counts"]) <= set(recorded))
    state.update(matrices=held.sort_values(["source", "accession", "platform"]), recorded=recorded,
                 summary=dict(study_platforms=len(matrices), matrices=len(held), samples=len(samples), samples_in_a_matrix=int(samples["in_matrix"].sum())))


def source_checks(qc, state):
    sources = pd.read_parquet(C.DATA / "annotation_sources.parquet")
    listed = dict(reversed(line.split("  ", 1)) for line in (C.RAW / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines())
    stored = {path.name for path in C.RAW.iterdir()} - {"SHA256SUMS.txt"}
    qc.check("stored_sources_complete_and_unchanged", stored == set(listed) and all(C.sha(C.RAW / name) == digest for name, digest in listed.items()))
    dated = sources[sources["raw_file"] != ""]
    names = [C.project_path(name).name for name in dated["raw_file"]]
    qc.check("stored_sources_are_the_recorded_sources_and_their_retrieval_notes",
             set(listed) == set(names) | {C.retrieval_note(C.RAW / name).name for name in names})
    notes = [C.read_json(C.retrieval_note(C.RAW / name)) for name in names]
    qc.check("every_stored_source_carries_the_date_of_its_retrieval_note", all(
        note["file"] == name and note["sha256"] == listed[name] and version == f"retrieved {note['retrieved_at_utc'][:10]}"
        for name, version, note in zip(names, dated["version"], notes)))
    genes = pd.read_parquet(C.DATA / "genes.parquet")
    qc.check("gene_reference_identifiers_unique", genes["hgnc_id"].is_unique and genes["symbol"].is_unique and (genes["hgnc_id"] != "").all())
    arrays = state["matrices"][state["matrices"]["family"] != "rna_seq"]
    named = {platform: entry.get("package", platform) for platform, entry in C.REGISTER["sources"].items()}
    qc.check("every_array_platform_has_a_source", set(arrays["platform"]) == set(named) and set(named.values()) <= set(sources["name"])
             and all(C.project_path(table).is_file() for table in sources["table"]))
    state["genes"] = genes.set_index("hgnc_id")
    state["summary"].update(genes=len(genes), annotation_sources=len(sources))


def matrix_checks(qc, state):
    genes, policy = state["genes"], C.CONFIG["gene_matrices"]
    mapped = pd.read_parquet(C.DATA / "feature_map_summary.parquet").set_index(["accession", "platform"])
    built = pd.read_parquet(C.DATA / "gene_matrix_summary.parquet").set_index(["accession", "platform"])
    libraries = pd.read_parquet(C.EXPRESSION / "count_table_samples.parquet").set_index(["study_accession", "platform", "sample_accession"])
    prior = C.read_json(C.EXPRESSION_OUT / "07.json")["prior_count"]
    qc.check("every_matrix_has_a_map_and_a_gene_matrix", len(mapped) == len(built) == len(state["matrices"])
             and set(mapped.index) == set(built.index) == set(zip(state["matrices"]["accession"], state["matrices"]["platform"])))
    maps, numbered, shapes, kept, sums, tally = [], [], [], [], [], dict(rows=0, rows_with_one_gene=0, gene_values=0)
    for number, matrix in enumerate(state["matrices"].itertuples(index=False), 1):
        key, name = (matrix.accession, matrix.platform), f"{matrix.accession}|{matrix.platform}"
        table = pd.read_parquet(C.expression_matrix(matrix.matrix, state["recorded"]))
        feature_map = pd.read_parquet(C.project_path(mapped.loc[key, "map"]))
        single = (feature_map["status"] == "one_gene").to_numpy()
        named = feature_map.loc[single, "hgnc_id"]
        counts = feature_map["status"].value_counts()
        if (feature_map["feature_id"].tolist() != table["feature_id"].tolist() or not set(feature_map["status"]) <= set(STATUSES)
                or not named.isin(genes.index).all() or (feature_map.loc[~single, "hgnc_id"] != "").any()
                or (feature_map.loc[single, "symbol"].to_numpy() != genes.loc[named, "symbol"].to_numpy()).any()
                or any(int(counts.get(status, 0)) != int(mapped.loc[key, status]) for status in STATUSES)):
            maps.append(name)
        if matrix.family == "rna_seq":
            # Row numbers of a table, read as Entrez genes, would name real genes; they run from 0 or 1 without a gap, and human Entrez genes have no 0
            # and none from 4 to 8.
            digits = feature_map["feature_id"][feature_map["feature_id"].str.fullmatch(C.CONFIG["mapping"]["entrez_gene"])]
            numbers = np.unique(digits.astype("int64").to_numpy())
            if len(numbers) and numbers[0] <= 1 and np.array_equal(numbers, np.arange(numbers[0], numbers[0] + len(numbers))):
                numbered.append(name)
        gene_matrix = pd.read_parquet(C.project_path(built.loc[key, "gene_matrix"]))
        rows = pd.read_parquet(C.project_path(built.loc[key, "gene_rows"]))
        order = sorted(set(named), key=lambda gene: int(gene.split(":")[1]))
        if (gene_matrix["hgnc_id"].tolist() != order or rows["hgnc_id"].tolist() != order or list(gene_matrix.columns[2:]) != list(table.columns[1:])
                or (gene_matrix["symbol"].to_numpy() != genes.loc[order, "symbol"].to_numpy()).any() or int(rows["rows"].sum()) != int(single.sum())
                or len(order) != int(built.loc[key, "genes"])):
            shapes.append(name)
            continue
        values, result = table.iloc[:, 1:].to_numpy(), gene_matrix.iloc[:, 2:].to_numpy()
        place = {feature: position for position, feature in enumerate(table["feature_id"])}
        code = {gene: index for index, gene in enumerate(order)}
        codes, positions = np.fromiter((code[gene] for gene in named), dtype="int64", count=len(named)), np.flatnonzero(single)
        rule = built.loc[key, "rule"]
        if rule == "highest_mean":
            # The gene-level row must be the kept row of the matrix, of the preferred probe quality where the gene has such a row, and no other such
            # row of the gene may have a higher mean.
            chosen = np.fromiter((place[feature] for feature in rows["kept_feature_id"]), dtype="int64", count=len(rows))
            block = values[positions]
            held = np.isfinite(block).sum(axis=1)
            means = np.where(held > 0, np.nansum(block, axis=1) / np.maximum(held, 1), -np.inf)
            quality = feature_map["probe_quality"].fillna("").astype(str).to_numpy()
            good = np.array([any(value.startswith(prefix) for prefix in policy["preferred_probe_quality"]) for value in quality[positions]], dtype=bool)
            has_good = np.zeros(len(order), dtype=bool)
            np.logical_or.at(has_good, codes, good)
            eligible = good | ~has_good[codes]
            highest = np.full(len(order), -np.inf)
            np.maximum.at(highest, codes[eligible], means[eligible])
            own = dict(zip(positions.tolist(), means.tolist()))
            if (not np.array_equal(result, values[chosen], equal_nan=True) or not single[chosen].all()
                    or (feature_map["hgnc_id"].to_numpy()[chosen] != np.array(order, dtype=object)).any()
                    or not np.allclose([own[position] for position in chosen.tolist()], highest, rtol=0, atol=1e-12, equal_nan=True)
                    or (rows["kept_probe_quality"].to_numpy() != quality[chosen]).any() or (rows["rows_counted"] != 1).any()):
                kept.append(name)
        else:
            raw = values
            if rule == "sum_counts":
                stored = pd.read_parquet(C.expression_matrix(matrix.counts, state["recorded"]))
                raw = stored.iloc[:, 1:].to_numpy(dtype="float64")
            # Rows of one gene with the same values in every sample are one measurement, counted once; a row of zeros adds nothing and is not a copy.
            distinct = ~(pd.DataFrame(raw[positions]).assign(gene=codes).duplicated().to_numpy() & (raw[positions] != 0).any(axis=1))
            expected = G.summed(raw[positions][distinct].astype("float64"), codes[distinct], len(order))
            counted = np.bincount(codes[distinct], minlength=len(order))
            if (rows["rows_counted"].to_numpy() != counted).any() or (rows["identical_rows"].to_numpy() != rows["rows"].to_numpy() - counted).any():
                sums.append(name)
            elif rule == "sum_counts":
                gene_counts = pd.read_parquet(C.project_path(built.loc[key, "gene_counts"]))
                own = libraries.loc[[(*key, sample) for sample in table.columns[1:]]]
                sizes, factors, groups = own["library_size"].to_numpy(dtype="float64"), own["norm_factor"].to_numpy(dtype="float64"), own["normalization_group"].to_numpy()
                effective = sizes * factors
                scaled = prior * effective / np.array([effective[groups == group].mean() for group in groups])
                again = np.log2((expected + scaled) / (effective + 2 * scaled) * 1e6)
                smallest = float(min((groups == group).sum() for group in set(groups.tolist())))
                rule_of = policy["expressed"]
                needed = smallest if smallest <= rule_of["large_n"] else rule_of["large_n"] + (smallest - rule_of["large_n"]) * rule_of["minimum_proportion"]
                flag = (((expected / effective * 1e6) >= rule_of["minimum_count"] / np.median(effective) * 1e6).sum(axis=1) >= needed - 1e-14) & (
                    expected.sum(axis=1) >= rule_of["minimum_total_count"] - 1e-14)
                if (gene_counts["hgnc_id"].tolist() != order or list(gene_counts.columns) != list(gene_matrix.columns)
                        or not np.allclose(gene_counts.iloc[:, 2:].to_numpy(dtype="float64"), expected, rtol=1e-12, atol=0)
                        or not (np.abs(again - result) <= policy["counts_tolerance"]).all() or (rows["expressed"].to_numpy(dtype=bool) != flag).any()):
                    sums.append(name)
            elif not np.allclose(result, expected, rtol=1e-9, atol=0) or rows["expressed"].notna().any():
                sums.append(name)
        tally["rows"] += len(feature_map)
        tally["rows_with_one_gene"] += int(single.sum())
        tally["gene_values"] += int(result.size)
        if number % 15 == 0 or number == len(state["matrices"]):
            print(f"[05] verified {number:,}/{len(state['matrices']):,} matrices, {tally['rows']:,} rows, {tally['gene_values']:,} gene-level values", flush=True)
    qc.check("every_map_holds_the_rows_of_its_matrix_with_one_status", not maps, maps)
    qc.check("no_rna_seq_matrix_is_identified_by_row_numbers", not numbered, numbered)
    qc.check("every_gene_matrix_holds_the_genes_and_samples_of_its_map", not shapes, shapes)
    qc.check("kept_rows_are_the_highest_mean_rows_of_their_matrices", not kept, kept)
    qc.check("added_rows_and_counts_equal_the_sums_of_their_distinct_rows", not sums, sums)
    state.update(mapped=mapped, built=built)
    state["summary"].update(tally)


def centred_ranks(values):
    """Rank the rows within each column, remove the mean rank of each row, and scale every column to unit variance."""
    ranks = pd.DataFrame(values).rank().to_numpy(dtype="float64")
    ranks -= ranks.mean(axis=1, keepdims=True)
    ranks -= ranks.mean(axis=0)
    spread = ranks.std(axis=0)
    return ranks / np.where(spread > 0, spread, 1.0)


def ncbi_checks(qc, state):
    """Compare the gene-level counts of every count study read from author tables or single-sample files with the NCBI-generated counts of its study."""
    policy, built = C.CONFIG["gene_matrices"]["ncbi_check"], state["built"]
    routes = pd.read_parquet(C.EXPRESSION / "processing_routes.parquet")
    listing = pd.read_parquet(C.EXPRESSION / "expression_files.parquet")
    listing = listing[listing["kind"] == "ncbi_generated_counts"]
    genes = state["genes"]
    hgnc_of = {str(entrez): gene for gene, entrez in genes["entrez_id"].items() if str(entrez).strip()}
    rows, failed = [], []
    for route in routes[(routes["tier"] == "counts") & (routes["data_source"] != "ncbi_counts") & routes["accession"].isin(listing["accession"])].itertuples(index=False):
        key = (route.accession, route.platform)
        item = next(listing[listing["accession"] == route.accession].itertuples(index=False))
        path = C.ROOT / "data/raw" / item.source / item.accession / item.file
        ours = pd.read_parquet(C.project_path(built.loc[key, "gene_counts"])).set_index("hgnc_id").drop(columns="symbol")
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        shared = [sample for sample in ours.columns if sample in header]
        row = dict(accession=route.accession, platform=route.platform, ncbi_table=item.file, status="too_few_samples", samples_compared=len(shared),
                   own_column_most_similar=0, shared_genes=0, genes_used=0, median_own_correlation=np.nan, smallest_margin=np.nan, not_most_similar="")
        if len(shared) >= policy["minimum_samples"]:
            theirs = pd.read_csv(path, sep="\t", usecols=[header[0], *shared], dtype={header[0]: str}).set_index(header[0])
            theirs = theirs[theirs.index.isin(hgnc_of)]
            theirs = theirs.groupby(theirs.index.map(hgnc_of)).sum()
            common = sorted(set(theirs.index) & set(ours.index), key=lambda gene: int(gene.split(":")[1]))
            common = [gene for gene, a, b in zip(common, ours.loc[common, shared].sum(axis=1) > 0, theirs.loc[common, shared].sum(axis=1) > 0) if a and b]
            row.update(shared_genes=len(common), status="too_few_shared_genes")
            if len(common) >= policy["minimum_shared_genes"]:
                chosen = [common[position] for position in np.unique(np.linspace(0, len(common) - 1, min(len(common), policy["genes_used"])).astype(int))]
                scores = centred_ranks(ours.loc[chosen, shared].to_numpy(dtype="float64")).T @ centred_ranks(theirs.loc[chosen, shared].to_numpy(dtype="float64"))
                scores /= len(chosen)
                hits = scores.argmax(axis=1) == np.arange(len(shared))
                others = np.where(np.eye(len(shared), dtype=bool), -np.inf, scores).max(axis=1)
                row.update(status="compared", own_column_most_similar=int(hits.sum()), genes_used=len(chosen),
                           median_own_correlation=round(float(np.median(np.diag(scores))), 6), smallest_margin=round(float((np.diag(scores) - others).min()), 6),
                           not_most_similar="; ".join(sample for sample, hit in zip(shared, hits) if not hit)[:500])
                if not hits.all():
                    failed.append(f"{route.accession}|{route.platform}")
        rows.append(row)
        print(f"[05] NCBI counts of {route.accession}|{route.platform}: {row['status']}, {row['own_column_most_similar']} of {row['samples_compared']} "
              f"most similar to their own column", flush=True)
    review = pd.DataFrame(rows, columns=NCBI_COLUMNS)
    C.save_review(review, C.OUT / "ncbi_check_review.csv")
    qc.check("every_compared_library_is_most_similar_to_its_own_ncbi_counts", not failed, failed)
    # A study platform that cannot be compared passes only as a recorded exception, so that the check cannot pass by comparing nothing.
    uncompared = {f"{row.accession}|{row.platform}": row.status for row in review.itertuples(index=False) if row.status != "compared"}
    excused = {name: entry["status"] for name, entry in C.REGISTER["ncbi_check_exceptions"].items()}
    qc.check("every_study_platform_with_ncbi_counts_compared_or_recorded_as_an_exception", uncompared == excused,
             dict(not_compared=uncompared, exceptions=excused))
    compared = review[review["status"] == "compared"]
    state["summary"].update(ncbi_check=dict(study_platforms_compared=len(compared), samples_compared=int(compared["samples_compared"].sum()),
                                            own_column_most_similar=int(compared["own_column_most_similar"].sum()), not_compared=uncompared))


def total_checks(qc, state):
    mapped, built, summary = state["mapped"], state["built"], state["summary"]
    third, fourth = C.read_json(C.receipt_path(3)), C.read_json(C.receipt_path(4))
    qc.check("map_receipt_reconciles", third["features"]["total"] == int(mapped["features"].sum()) == summary["rows"]
             and third["features"]["one_gene"] == int(mapped["one_gene"].sum()) == summary["rows_with_one_gene"] and third["matrices"] == len(mapped))
    qc.check("gene_matrix_receipt_reconciles", fourth["rows"]["used"] == int(built["rows_used"].sum()) == summary["rows_with_one_gene"]
             and fourth["matrices"] == len(built) and fourth["samples"] == int(built["samples"].sum()) == summary["samples_in_a_matrix"])
    files = {folder.name: len(list(folder.glob("*.parquet"))) for folder in (C.DATA / "feature_maps", G.GENES, G.ROWS, G.COUNTS)}
    counted = int((built["gene_counts"] != "").sum())
    qc.check("one_file_for_every_map_and_gene_matrix", files == {"feature_maps": len(mapped), "gene_matrices": len(built), "gene_rows": len(built), "gene_counts": counted})
    summary.update(gene_matrices=len(built), gene_count_matrices=counted, fewest_genes_in_a_matrix=int(built["genes"].min()),
                   median_genes_in_a_matrix=int(built["genes"].median()), most_genes_in_a_matrix=int(built["genes"].max()),
                   median_protein_coding_genes_in_a_matrix=int(built["protein_coding_genes"].median()))


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/6 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/6 Indexes against the expression stage", lambda: index_checks(qc, state))
        qc.phase("3/6 Annotation sources and gene reference", lambda: source_checks(qc, state))
        qc.phase("4/6 Maps and gene matrices against the expression matrices", lambda: matrix_checks(qc, state))
        qc.phase("5/6 Gene-level counts against the NCBI-generated counts", lambda: ncbi_checks(qc, state))
        qc.phase("6/6 Stage receipts and totals", lambda: total_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, 5) if C.receipt_path(number).exists()]),
                  eligibility="ANNOTATION_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[05 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    {summary['rows']:,} rows of {summary['matrices']} matrices are mapped; {summary['gene_matrices']} gene matrices hold "
              f"{summary['fewest_genes_in_a_matrix']:,} to {summary['most_genes_in_a_matrix']:,} genes (median {summary['median_genes_in_a_matrix']:,})", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
