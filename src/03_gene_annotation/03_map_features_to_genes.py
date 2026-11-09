"""Relate every row of every expression matrix to HGNC genes, keeping every row with the gene it maps to or the reason it maps to none."""
import collections
import importlib
import re
import shutil

C = importlib.import_module("01_index_expression_matrices")

MAPS = C.DATA / "feature_maps"
GENE_REFERENCE = "HGNC complete gene set"
STATUSES = ("one_gene", "several_genes", "no_gene_in_annotation", "gene_not_in_hgnc", "accession_not_resolved", "not_in_annotation")
DIALECTS = ("probe_id", "array_address", "nuid")
MAP_COLUMNS = ["feature_id", "status", "hgnc_id", "symbol", "entrez_id", "ensembl_gene_id", "locus_group", "candidate_genes", "annotation_source",
               "annotation_key", "annotation_gene", "matched_by", "probe_quality"]
SUMMARY_COLUMNS = ["source", "accession", "platform", "family", "annotation_source", "annotation_key", "features", "in_annotation", *STATUSES, "genes",
                   "protein_coding_genes", "map"]


def gene_reference(genes):
    """Index the HGNC genes by Entrez gene, Ensembl gene, approved symbol, previous symbol and alias."""
    index = {name: collections.defaultdict(set) for name in ("entrez", "ensembl", "approved", "previous", "alias")}
    for gene in genes.itertuples(index=False):
        for name, values in (("entrez", [gene.entrez_id]), ("ensembl", [gene.ensembl_gene_id]), ("approved", [gene.symbol]),
                             ("previous", gene.prev_symbol.split("|")), ("alias", gene.alias_symbol.split("|"))):
            for value in values:
                if value.strip():
                    index[name][value.strip()].add(gene.hgnc_id)
    return index


def genes_of(kind, values, index, accessions):
    """Return the HGNC genes that identifiers of one kind name, and how the first of them was matched."""
    found, matched = set(), ""
    for value in values:
        if kind == "accession":
            hits, how = set().union(*(index["entrez"].get(gene, ()) for gene in accessions.get(value.split(".")[0], ()))), "accession"
        elif kind == "symbol":
            # A symbol names the gene that carries it as its approved symbol; failing that, the genes that carried it before or list it as an alias.
            hits, how = next(((index[name][value], f"{name}_symbol") for name in ("approved", "previous", "alias") if value in index[name]), ((), ""))
        else:
            hits, how = index[kind].get(value, ()), f"{kind}_gene_id"
        if hits:
            found |= set(hits)
            matched = matched or how
    return found, matched


def package_map(table, index):
    """For every probe of an annotation package: the Entrez genes it names, the HGNC genes they are, and its other identifiers."""
    probes = {}
    for row in table.itertuples(index=False):
        entry = probes.setdefault(row.probe_id, dict(named=[], array_address=row.array_address, nuid=row.nuid, probe_quality=row.probe_quality))
        if row.entrez_id:
            entry["named"].append(row.entrez_id)
    for entry in probes.values():
        entry["genes"], entry["matched_by"] = genes_of("entrez", entry["named"], index, {})
    return probes


def table_map(table, entry, index, accessions, separator):
    """For every identifier of a platform table: what its gene columns name and the HGNC genes that is, by the first column that leads to a gene."""
    probes = {}
    for row in table.to_dict("records"):
        found = probes.setdefault(row[entry["identifier"]].strip(), dict(named=[], genes=set(), matched_by="", probe_quality="", only_accessions=True))
        for kind in ("entrez", "accession", "symbol"):
            values = [value for value in separator.split(row.get(entry["genes"].get(kind, ""), "").strip()) if value and value != "---"]
            if not values:
                continue
            # A row that names only accessions, none of which org.Hs.eg.db relates to an Entrez gene, is unresolved rather than outside HGNC.
            found["only_accessions"] = found["only_accessions"] and kind == "accession" and not any(accessions.get(value.split(".")[0]) for value in values)
            found["named"] += [value for value in values if value not in found["named"]]
            genes, matched = genes_of(kind, values, index, accessions)
            if genes:
                found["genes"] |= genes
                found["matched_by"] = found["matched_by"] or matched
                break
    return probes


def row_numbers(features, entrez):
    """Whether the numeric identifiers of a matrix run from 0 or 1 without a gap, as row numbers do; human Entrez genes have no 0 and none from 4 to 8."""
    numbers = sorted({int(feature) for feature in features if entrez.fullmatch(feature)})
    return bool(numbers) and numbers[0] in (0, 1) and numbers[-1] - numbers[0] + 1 == len(numbers)


def rows_of(features, lookup, source, key, genes):
    """Write one row for every feature of a matrix: its gene, or the reason it has none."""
    rows = []
    for feature in features:
        found = lookup(feature)
        row = dict(feature_id=feature, hgnc_id="", symbol="", entrez_id="", ensembl_gene_id="", locus_group="", candidate_genes="", annotation_source=source,
                   annotation_key=key, annotation_gene="", matched_by="", probe_quality="")
        if found is None:
            rows.append(dict(row, status="not_in_annotation"))
            continue
        candidates = sorted(found["genes"])
        status = ("no_gene_in_annotation" if not found["named"] else "accession_not_resolved" if not candidates and found.get("only_accessions")
                  else "gene_not_in_hgnc" if not candidates else "one_gene" if len(candidates) == 1 else "several_genes")
        row.update(status=status, annotation_gene="|".join(found["named"]), matched_by=found["matched_by"], probe_quality=found.get("probe_quality", ""),
                   candidate_genes="|".join(candidates) if len(candidates) > 1 else "")
        if status == "one_gene":
            row.update(genes[candidates[0]])
        rows.append(row)
    return rows


def main():
    import pandas as pd
    import pyarrow.parquet as pq

    print(f"[03] Mapping features to genes v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy, platforms, rewrites = C.CONFIG["mapping"], C.REGISTER["sources"], C.REGISTER["identifier_replace"]
    recorded = C.expression_binding()[1]
    matrices = pd.read_parquet(C.DATA / "matrix_index.parquet")
    matrices = matrices[matrices["matrix"] != ""].sort_values(["source", "accession", "platform"])
    unknown = sorted(set(rewrites) - set(matrices["accession"] + "|" + matrices["platform"]))
    C.require(not unknown, f"Identifier rewrites recorded for study platforms without a matrix: {unknown}")
    reference = pd.read_parquet(C.DATA / "genes.parquet")
    index = gene_reference(reference)
    genes = {gene.hgnc_id: dict(hgnc_id=gene.hgnc_id, symbol=gene.symbol, entrez_id=gene.entrez_id, ensembl_gene_id=gene.ensembl_gene_id,
                                locus_group=gene.locus_group) for gene in reference.itertuples(index=False)}
    coding = {gene for gene, item in genes.items() if item["locus_group"] == "protein-coding gene"}
    accessions = collections.defaultdict(set)
    for row in pd.read_parquet(C.DATA / "accession_genes.parquet").itertuples(index=False):
        accessions[row.accession].add(row.entrez_id)
    separator, entrez, ensembl = re.compile(policy["value_separator"]), re.compile(policy["entrez_gene"]), re.compile(policy["ensembl_gene"])
    shutil.rmtree(MAPS, ignore_errors=True)
    MAPS.mkdir(parents=True)

    annotations, summary, outputs = {}, [], []
    for number, matrix in enumerate(matrices.itertuples(index=False), 1):
        key, name = f"{matrix.accession}|{matrix.platform}", f"{matrix.accession}__{matrix.platform}"
        features = pq.read_table(C.expression_matrix(matrix.matrix, recorded), columns=["feature_id"]).column("feature_id").to_pylist()
        if matrix.family == "rna_seq":
            # The rows of an RNA-seq matrix are gene identifiers themselves: an Entrez gene, an Ensembl gene or a symbol, named by its own shape.
            replace = rewrites.get(key, {}).get("replace", [])
            C.require(not row_numbers(features, entrez), f"{key}: its numeric row identifiers run from 0 or 1 without a gap, as the row numbers of a table do; "
                                                         "they cannot be read as Entrez genes")

            def lookup(feature):
                gene = ensembl.fullmatch(feature)
                kind, value = ("entrez", feature) if entrez.fullmatch(feature) else ("ensembl", gene.group(1)) if gene else ("symbol", feature)
                found, matched = genes_of(kind, [value], index, accessions)
                if kind == "symbol" and not found and replace:
                    # A recorded rewrite is tried only for a symbol that names no gene as it is written.
                    rewritten = value
                    for pattern, replacement in replace:
                        rewritten = re.sub(pattern, replacement, rewritten)
                    found, matched = genes_of(kind, [rewritten], index, accessions)
                    value = rewritten if found else value
                return dict(named=[value], genes=found, matched_by=matched)

            source, dialect = GENE_REFERENCE, "gene_identifier"
        else:
            entry = platforms[matrix.platform]
            source = entry.get("package", matrix.platform)
            if source not in annotations:
                if entry["source"] == "bioconductor":
                    probes = package_map(pd.read_parquet(C.DATA / "package_annotations" / f"{source}.parquet"), index)
                    # An Illumina matrix may name its probes by probe identifier, by array address or by nuID; each is a way into the same annotation.
                    annotations[source] = {dialect: {(probe if dialect == "probe_id" else item[dialect]): item for probe, item in probes.items()
                                                     if dialect == "probe_id" or item[dialect]} for dialect in DIALECTS}
                else:
                    table = pd.read_parquet(C.DATA / "platform_tables" / f"{source}.parquet")
                    annotations[source] = {entry["identifier"]: table_map(table, entry, index, accessions, separator)}
            shared = {dialect: sum(feature in probes for feature in features) for dialect, probes in annotations[source].items()}
            dialect = max(shared, key=shared.get)
            C.require(shared[dialect] >= policy["minimum_identifier_share"] * len(features),
                      f"{key}: only {shared[dialect]:,} of its {len(features):,} rows are identifiers of {source}; the recorded annotation does not fit the matrix")
            lookup = annotations[source][dialect].get
        table = pd.DataFrame(rows_of(features, lookup, source, dialect, genes), columns=MAP_COLUMNS)
        C.require(len(table) == matrix.features and table["feature_id"].tolist() == features, f"{key}: the map does not hold every row of the matrix")
        path = MAPS / f"{name}.parquet"
        C.save_table(table, path)
        outputs.append(path)
        counts, mapped = table["status"].value_counts(), set(table.loc[table["status"] == "one_gene", "hgnc_id"])
        summary.append(dict(source=matrix.source, accession=matrix.accession, platform=matrix.platform, family=matrix.family, annotation_source=source,
                            annotation_key=dialect, features=len(table), in_annotation=int((table["status"] != "not_in_annotation").sum()),
                            **{status: int(counts.get(status, 0)) for status in STATUSES}, genes=len(mapped), protein_coding_genes=len(mapped & coding),
                            map=path.relative_to(C.ROOT).as_posix()))
        if number % 15 == 0 or number == len(matrices):
            print(f"[03] mapped {number:,}/{len(matrices):,} matrices", flush=True)

    summary = pd.DataFrame(summary, columns=SUMMARY_COLUMNS)
    C.save_table(summary, C.DATA / "feature_map_summary.parquet")
    C.save_review(summary, C.OUT / "feature_maps_review.csv")
    totals = {status: int(summary[status].sum()) for status in STATUSES}
    by_family = {str(family): dict(matrices=len(part), features=int(part["features"].sum()), one_gene=int(part["one_gene"].sum()),
                                   median_protein_coding_genes=int(part["protein_coding_genes"].median())) for family, part in summary.groupby("family")}
    C.finish(3, start, [C.DATA / "feature_map_summary.parquet", C.OUT / "feature_maps_review.csv", *outputs], method=policy["method"],
             matrices=len(summary), features=dict(total=int(summary["features"].sum()), **totals), by_family=by_family,
             genes=dict(reference=len(genes), protein_coding=len(coding), fewest_protein_coding_in_a_matrix=int(summary["protein_coding_genes"].min()),
                        most_protein_coding_in_a_matrix=int(summary["protein_coding_genes"].max())),
             minimum_identifier_share=policy["minimum_identifier_share"], identifier_replace=rewrites, protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: {int(summary['features'].sum()):,} rows of {len(summary)} matrices kept; {totals['one_gene']:,} map to one gene; "
          f"review sheet: {(C.OUT / 'feature_maps_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for status in STATUSES[1:]:
        print(f"[03]   {status}: {totals[status]:,} rows", flush=True)
    print(f"[03]   protein-coding genes reached in a matrix: {int(summary['protein_coding_genes'].min()):,} to {int(summary['protein_coding_genes'].max()):,} "
          f"of {len(coding):,} (median {int(summary['protein_coding_genes'].median()):,})", flush=True)


if __name__ == "__main__":
    main()
