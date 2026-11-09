"""Build one gene-level matrix for every expression matrix from the rows that map to one gene, each study platform on its own."""
import importlib
import shutil

C = importlib.import_module("01_index_expression_matrices")

GENES = C.DATA / "gene_matrices"
COUNTS = C.DATA / "gene_counts"
ROWS = C.DATA / "gene_rows"
ROW_COLUMNS = ["hgnc_id", "symbol", "locus_group", "rows", "feature_ids", "matched_by", "rows_counted", "identical_rows", "kept_feature_id",
               "kept_probe_quality", "rule", "expressed"]
SUMMARY_COLUMNS = ["source", "accession", "platform", "family", "value_scale", "rule", "samples", "features", "rows_used", "rows_left_out", "genes",
                   "protein_coding_genes", "genes_with_several_rows", "identical_rows", "kept_rows_of_preferred_quality", "expressed_genes", "gene_matrix",
                   "gene_counts", "gene_rows"]


def log_cpm(counts, sizes, factors, prior, groups=None):
    """Log2 counts per million as edgeR computes them: library sizes scaled by their TMM factors, and a prior count scaled to each library by the mean
    scaled library size of the samples that were normalized together (its normalization group)."""
    import numpy as np

    effective = sizes * factors
    labels = np.asarray(groups if groups is not None else [""] * len(effective))
    scaled = prior * effective / np.array([effective[labels == label].mean() for label in labels])
    return np.log2((counts + scaled) / (effective + 2 * scaled) * 1e6)


def identical_copies(values, codes):
    """Mark every row that holds the same values in every sample as an earlier row of the same gene; such copies count once in a sum. A row of zeros
    adds nothing to a sum and is not called a copy."""
    import numpy as np

    seen, copies = set(), np.zeros(len(codes), dtype=bool)
    for place, (code, row) in enumerate(zip(codes, values)):
        if not np.any(row):
            continue
        signature = (int(code), np.ascontiguousarray(row, dtype="float64").tobytes())
        copies[place] = signature in seen
        seen.add(signature)
    return copies


def expressed(counts, sizes, factors, groups, policy):
    """The genes that edgeR's filterByExpr keeps with its default limits, with the normalization groups as groups and the TMM-scaled library sizes."""
    import numpy as np

    effective = sizes * factors
    labels = np.asarray(groups if groups is not None else [""] * len(effective))
    smallest = float(min((labels == label).sum() for label in set(labels.tolist())))
    if smallest > policy["large_n"]:
        smallest = policy["large_n"] + (smallest - policy["large_n"]) * policy["minimum_proportion"]
    cutoff = policy["minimum_count"] / np.median(effective) * 1e6
    above = (counts / effective * 1e6 >= cutoff).sum(axis=1) >= smallest - 1e-14
    return above & (counts.sum(axis=1) >= policy["minimum_total_count"] - 1e-14)


def preferred(qualities, prefixes):
    import numpy as np

    return np.array([bool(quality) and any(str(quality).startswith(prefix) for prefix in prefixes) for quality in qualities], dtype=bool)


def gene_groups(feature_map, genes):
    """Group the rows that map to one gene: for each gene, in the order of the HGNC identifiers, the positions of its rows in the matrix."""
    import numpy as np

    single = feature_map["status"].to_numpy() == "one_gene"
    positions, named = np.flatnonzero(single), feature_map.loc[single, "hgnc_id"].to_numpy()
    order = sorted(set(named), key=lambda gene: int(gene.split(":")[1]))
    code = {gene: index for index, gene in enumerate(order)}
    codes = np.fromiter((code[gene] for gene in named), dtype="int64", count=len(named))
    C.require(all(gene in genes for gene in order), "A map names a gene that the gene reference does not hold")
    return order, positions, codes


def summed(values, codes, count):
    import numpy as np

    result = np.zeros((count, values.shape[1]), dtype="float64")
    np.add.at(result, codes, values)
    return result


def main():
    import numpy as np
    import pandas as pd

    print(f"[04] Building gene matrices v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["gene_matrices"]
    recorded = C.expression_binding()[1]
    matrices = pd.read_parquet(C.DATA / "matrix_index.parquet")
    matrices = matrices[matrices["matrix"] != ""].sort_values(["source", "accession", "platform"])
    C.require(set(matrices["value_scale"]) <= set(policy["rules"]), f"No rule is recorded for the value scales {sorted(set(matrices['value_scale']) - set(policy['rules']))}")
    reference = pd.read_parquet(C.DATA / "genes.parquet", columns=["hgnc_id", "symbol", "locus_group"]).set_index("hgnc_id")
    genes = reference.to_dict("index")
    libraries = pd.read_parquet(C.EXPRESSION / "count_table_samples.parquet").set_index(["study_accession", "platform", "sample_accession"])
    prior = C.read_json(C.EXPRESSION_OUT / "07.json")["prior_count"]
    for folder in (GENES, COUNTS, ROWS):
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir(parents=True)

    summary, outputs = [], []
    for number, matrix in enumerate(matrices.itertuples(index=False), 1):
        key, name, rule = f"{matrix.accession}|{matrix.platform}", f"{matrix.accession}__{matrix.platform}", policy["rules"][matrix.value_scale]
        table = pd.read_parquet(C.expression_matrix(matrix.matrix, recorded))
        feature_map = pd.read_parquet(C.DATA / "feature_maps" / f"{name}.parquet")
        C.require(table["feature_id"].tolist() == feature_map["feature_id"].tolist(), f"{key}: the map does not hold the rows of the matrix in their order")
        samples, values = list(table.columns[1:]), table.iloc[:, 1:].to_numpy()
        order, positions, codes = gene_groups(feature_map, genes)
        C.require(len(order) > 0, f"{key}: no row maps to one gene")
        features = feature_map["feature_id"].to_numpy()[positions]
        qualities = feature_map["probe_quality"].fillna("").to_numpy()[positions]
        rows = pd.DataFrame(dict(code=codes, feature=features, matched=feature_map["matched_by"].fillna("").to_numpy()[positions])).groupby("code").agg(
            size=("feature", "size"), join=("feature", "|".join), matched=("matched", "|".join))
        kept, kept_quality, gene_counts, flags = [""] * len(order), [""] * len(order), None, pd.Series(pd.NA, index=range(len(order)), dtype="boolean")
        copies = np.zeros(len(codes), dtype=bool)
        if rule == "highest_mean":
            # Of the rows of one gene, the one with the highest mean over the samples of the study is kept as it is; a row without any value comes last.
            # Where the annotation rates its probes, a probe of a preferred quality goes first whenever the gene has one.
            block = values[positions]
            held = np.isfinite(block).sum(axis=1)
            means = np.where(held > 0, np.nansum(block, axis=1) / np.maximum(held, 1), -np.inf)
            good = preferred(qualities, policy["preferred_probe_quality"])
            ranked = pd.DataFrame(dict(code=codes, good=good, mean=means, place=np.arange(len(codes)))).sort_values(
                ["code", "good", "mean", "place"], ascending=[True, False, False, True])
            chosen = ranked.drop_duplicates("code")["place"].to_numpy()
            result, kept, kept_quality = block[chosen], features[chosen].tolist(), qualities[chosen].tolist()
            counted = np.ones(len(order), dtype=int)
        elif rule == "sum":
            # Values on a linear scale add up over the rows of one gene; identical copies of a row count once.
            C.require(bool(np.isfinite(values[positions]).all()), f"{key}: a value is missing in a matrix whose rows are added up")
            copies = identical_copies(values[positions], codes)
            result = summed(values[positions][~copies], codes[~copies], len(order))
        else:
            # Counts add up over the rows of one gene, identical copies once; the normalized values are computed again with the library size, the TMM
            # factor and the normalization group that each sample has in the expression stage, so that they stay on the scale of the matrix they come from.
            counts = pd.read_parquet(C.expression_matrix(matrix.counts, recorded))
            C.require(list(counts.columns) == list(table.columns) and counts["feature_id"].tolist() == table["feature_id"].tolist(),
                      f"{key}: the count matrix does not hold the rows and samples of the normalized matrix")
            own = libraries.loc[[(matrix.accession, matrix.platform, sample) for sample in samples]]
            sizes, factors = own["library_size"].to_numpy(dtype="float64"), own["norm_factor"].to_numpy(dtype="float64")
            raw = counts.iloc[:, 1:].to_numpy(dtype="float64")
            C.require(bool(np.allclose(raw.sum(axis=0), sizes, rtol=1e-9)), f"{key}: the counts do not add up to the library sizes of the expression stage")
            groups = own["normalization_group"].tolist()
            C.require(float(np.abs(log_cpm(raw, sizes, factors, prior, groups) - values).max()) <= policy["counts_tolerance"],
                      f"{key}: the stored normalized values are not reproduced with the normalization groups of the samples")
            copies = identical_copies(raw[positions], codes)
            gene_counts = summed(raw[positions][~copies], codes[~copies], len(order))
            result = log_cpm(gene_counts, sizes, factors, prior, groups).astype("float32")
            flags = pd.Series(expressed(gene_counts, sizes, factors, groups, policy["expressed"]), dtype="boolean")
        if rule != "highest_mean":
            counted = np.bincount(codes[~copies], minlength=len(order))
        labels = dict(hgnc_id=order, symbol=[genes[gene]["symbol"] for gene in order])
        path = GENES / f"{name}.parquet"
        C.save_table(pd.concat([pd.DataFrame(labels), pd.DataFrame(result, columns=samples)], axis=1), path)
        rows_path = ROWS / f"{name}.parquet"
        C.save_table(pd.DataFrame(dict(labels, locus_group=[genes[gene]["locus_group"] for gene in order], rows=rows["size"].to_numpy(),
                                       feature_ids=rows["join"].to_numpy(), matched_by=rows["matched"].to_numpy(), rows_counted=counted,
                                       identical_rows=rows["size"].to_numpy() - counted if rule != "highest_mean" else 0, kept_feature_id=kept,
                                       kept_probe_quality=kept_quality, rule=rule, expressed=flags.array), columns=ROW_COLUMNS), rows_path)
        outputs += [path, rows_path]
        counts_path = ""
        if gene_counts is not None:
            whole = bool((gene_counts == np.round(gene_counts)).all()) and gene_counts.max() < 2 ** 31
            counts_path = COUNTS / f"{name}.parquet"
            C.save_table(pd.concat([pd.DataFrame(labels), pd.DataFrame(gene_counts.astype("int32") if whole else gene_counts, columns=samples)], axis=1), counts_path)
            outputs.append(counts_path)
        coding = sum(genes[gene]["locus_group"] == "protein-coding gene" for gene in order)
        summary.append(dict(source=matrix.source, accession=matrix.accession, platform=matrix.platform, family=matrix.family, value_scale=matrix.value_scale,
                            rule=rule, samples=len(samples), features=len(table), rows_used=len(positions), rows_left_out=len(table) - len(positions),
                            genes=len(order), protein_coding_genes=coding, genes_with_several_rows=int((rows["size"] > 1).sum()),
                            identical_rows=int(copies.sum()), kept_rows_of_preferred_quality=int(preferred(kept_quality, policy["preferred_probe_quality"]).sum()),
                            expressed_genes=int(flags.sum()) if flags.notna().any() else 0,
                            gene_matrix=path.relative_to(C.ROOT).as_posix(), gene_counts=counts_path.relative_to(C.ROOT).as_posix() if counts_path else "",
                            gene_rows=rows_path.relative_to(C.ROOT).as_posix()))
        if number % 15 == 0 or number == len(matrices):
            print(f"[04] built {number:,}/{len(matrices):,} gene matrices", flush=True)

    summary = pd.DataFrame(summary, columns=SUMMARY_COLUMNS)
    C.save_table(summary, C.DATA / "gene_matrix_summary.parquet")
    C.save_review(summary, C.OUT / "gene_matrices_review.csv")
    by_rule = {str(rule): dict(matrices=len(part), samples=int(part["samples"].sum()), median_genes=int(part["genes"].median()))
               for rule, part in summary.groupby("rule")}
    C.finish(4, start, [C.DATA / "gene_matrix_summary.parquet", C.OUT / "gene_matrices_review.csv", *outputs], method=policy["method"],
             matrices=len(summary), samples=int(summary["samples"].sum()), by_rule=by_rule,
             rows=dict(total=int(summary["features"].sum()), used=int(summary["rows_used"].sum()), left_out=int(summary["rows_left_out"].sum())),
             genes=dict(fewest=int(summary["genes"].min()), median=int(summary["genes"].median()), most=int(summary["genes"].max()),
                        fewest_protein_coding=int(summary["protein_coding_genes"].min()), median_protein_coding=int(summary["protein_coding_genes"].median())),
             identical_rows=dict(total=int(summary["identical_rows"].sum()),
                                 by_study_platform={f"{row.accession}|{row.platform}": int(row.identical_rows) for row in summary.itertuples(index=False) if row.identical_rows}),
             rules=policy["rules"], counts_tolerance=policy["counts_tolerance"], prior_count=prior, preferred_probe_quality=policy["preferred_probe_quality"],
             expressed=policy["expressed"], refinements=policy["refinements"], protocol=C.CONFIG["protocol"])
    print(f"[04] Complete: {len(summary)} gene matrices for {int(summary['samples'].sum()):,} samples; {int(summary['genes'].min()):,} to "
          f"{int(summary['genes'].max()):,} genes in a matrix (median {int(summary['genes'].median()):,}); review sheet: "
          f"{(C.OUT / 'gene_matrices_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for rule, item in by_rule.items():
        print(f"[04]   {rule}: {item['matrices']} matrices, {item['samples']:,} samples", flush=True)


if __name__ == "__main__":
    main()
