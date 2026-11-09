"""Collect the gene reference and the annotation source of every array platform, and record the version and the checksum of each."""
import argparse
import csv
import importlib
import io
import re
import shutil
import subprocess
from datetime import datetime, timezone

C = importlib.import_module("01_index_expression_matrices")

WORK = C.WORK / "sources"
PACKAGES = C.DATA / "package_annotations"
PLATFORMS = C.DATA / "platform_tables"
SOURCE_COLUMNS = ["kind", "name", "version", "built_from", "platforms", "rows", "identifiers", "with_a_gene", "table", "raw_file"]
EXPORT = r"""
suppressPackageStartupMessages({library(jsonlite); library(AnnotationDbi)})
job <- fromJSON(commandArgs(TRUE)[1])
described <- list()
for (name in job$packages) {
  suppressPackageStartupMessages(library(name, character.only = TRUE))
  db <- get(name)
  probes <- keys(db, keytype = "PROBEID")
  # One row for every pair of a probe and an Entrez gene; a probe without a gene keeps one row with no gene.
  pairs <- suppressMessages(select(db, keys = probes, columns = c("ENTREZID", "SYMBOL"), keytype = "PROBEID"))
  names(pairs) <- c("probe_id", "entrez_id", "symbol")
  prefix <- sub("[.]db$", "", name)
  for (extra in names(job$extras)) {
    map <- paste0(prefix, job$extras[[extra]])
    pairs[[extra]] <- ""
    if (exists(map)) {
      values <- as.list(get(map))
      pairs[[extra]] <- vapply(values[pairs$probe_id], function(value) if (length(value) && !is.na(value[1])) as.character(value[1]) else "", "")
    }
  }
  pairs[is.na(pairs)] <- ""
  write.table(pairs, file.path(job$work, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  info <- metadata(db)
  described[[name]] <- list(version = as.character(packageVersion(name)), probes = length(probes), metadata = as.list(setNames(info$value, info$name)))
}

# GenBank and RefSeq accessions of the platform tables, related to Entrez genes by the organism package.
organism <- job$organism
suppressPackageStartupMessages(library(organism, character.only = TRUE))
accessions <- readLines(file.path(job$work, "accessions.txt"))
known <- accessions[accessions %in% keys(get(organism), keytype = "ACCNUM")]
found <- data.frame(ACCNUM = character(), ENTREZID = character())
if (length(known)) found <- suppressMessages(select(get(organism), keys = known, columns = "ENTREZID", keytype = "ACCNUM"))
names(found) <- c("accession", "entrez_id")
write.table(found[!is.na(found$entrez_id), ], file.path(job$work, "accessions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
info <- metadata(get(organism))
described[[organism]] <- list(version = as.character(packageVersion(organism)), probes = length(accessions), metadata = as.list(setNames(info$value, info$name)))
write(toJSON(described, auto_unbox = TRUE, null = "null"), file.path(job$work, "packages.json"))
"""


def platform_table(path, accession):
    """Read the probe table of a full GEO platform record with every column as text, and the number of rows the record declares."""
    import pandas as pd

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    C.require(lines and lines[0].strip() == f"^PLATFORM = {accession}", f"{accession}: GEO returned no platform record")
    declared = next((int(line.partition(" = ")[2]) for line in lines if line.startswith("!Platform_data_row_count")), None)
    first, last = lines.index("!platform_table_begin"), lines.index("!platform_table_end")
    table = pd.read_csv(io.StringIO("\n".join(lines[first + 1:last])), sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    C.require(len(table) == declared and len(table) > 0, f"{accession}: the platform table holds {len(table):,} rows and the record declares {declared}")
    return table


def export_packages(packages, accessions, policy):
    """Export the probe-to-gene pairs of the annotation packages, and the Entrez genes of the named accessions, with R."""
    import pandas as pd

    shutil.rmtree(WORK, ignore_errors=True)
    WORK.mkdir(parents=True)
    try:
        (WORK / "accessions.txt").write_text("".join(f"{accession}\n" for accession in accessions), encoding="utf-8", newline="\n")
        C.write_json(WORK / "job.json", dict(work=WORK.as_posix(), packages=packages, extras=policy["package_extras"], organism=policy["organism_package"]))
        (WORK / "export.R").write_text(EXPORT, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (WORK / "export.R").as_posix(), (WORK / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"The annotation packages could not be exported:\n{done.stderr[-3000:]}")
        read = lambda name: pd.read_csv(WORK / f"{name}.tsv", sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
        return {name: read(name) for name in packages}, read("accessions"), C.read_json(WORK / "packages.json")
    finally:
        C.clear_work(WORK)


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-search", action="store_true",
                        help="Retrieve every source again instead of reading the stored sources under data/raw/annotation")
    args = parser.parse_args()
    print(f"[02] Collecting annotation sources v{C.VERSION}", flush=True)
    start = C.begin(2)
    policy, assigned = C.CONFIG["sources"], C.REGISTER["sources"]
    matrices = pd.read_parquet(C.DATA / "matrix_index.parquet")
    arrays = matrices[(matrices["matrix"] != "") & (matrices["family"] != "rna_seq")]
    needed = sorted(set(arrays["platform"]))
    C.require(needed == sorted(assigned), f"The sources of {C.CONFIG['register']} must record every array platform with a matrix once: "
                                          f"{sorted(set(needed) ^ set(assigned))}")
    listing = C.RAW / "SHA256SUMS.txt"
    C.require(args.update_search or listing.is_file(), f"No checksum list of the stored sources at data/raw/annotation/SHA256SUMS.txt; {C.STORED_ANSWERS}")
    listed = "" if args.update_search else listing.read_text(encoding="utf-8")
    if C.RAW.exists() and args.update_search:
        shutil.rmtree(C.RAW)
    for folder in (PACKAGES, PLATFORMS):
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir(parents=True)
    users = lambda test: {name: sorted(platform for platform, entry in assigned.items() if test(platform, entry) == name)
                          for name in {test(platform, entry) for platform, entry in assigned.items()} - {None}}
    sources, outputs = [], []

    # 1. The gene reference: the complete HGNC gene set, every column kept.
    print("[02] Retrieving the HGNC gene set", flush=True)
    path, retrieved = C.retrieve(policy["hgnc"]["file"], policy["hgnc"]["url"], update=args.update_search)
    stored = [path, C.retrieval_note(path)]
    genes = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)
    genes = genes.apply(lambda column: column.str.strip('"'))
    C.require(set(policy["hgnc"]["columns"]) <= set(genes.columns) and genes["hgnc_id"].is_unique and genes["symbol"].is_unique and (genes["hgnc_id"] != "").all(),
              "The HGNC table lacks a recorded column, or an HGNC identifier or symbol is repeated")
    C.save_table(genes, C.DATA / "genes.parquet")
    outputs.append(C.DATA / "genes.parquet")
    sources.append(dict(kind="gene_reference", name="HGNC complete gene set", version=f"retrieved {retrieved['retrieved_at_utc'][:10]}", built_from=policy["hgnc"]["url"], platforms="",
                        rows=len(genes), identifiers=len(genes), with_a_gene=len(genes), table=(C.DATA / "genes.parquet").relative_to(C.ROOT).as_posix(),
                        raw_file=path.relative_to(C.ROOT).as_posix()))

    # 2. The platform tables of GEO for the chips without an annotation package, every column kept.
    accessions, separator = set(), re.compile(C.CONFIG["mapping"]["value_separator"])
    for name, platforms in sorted(users(lambda platform, entry: platform if entry["source"] == "geo_platform_table" else None).items()):
        entry = assigned[name]
        print(f"[02] Retrieving the platform table of {name}", flush=True)
        path, retrieved = C.retrieve(f"{name}.soft", policy["geo_platform"]["url"], dict(policy["geo_platform"]["parameters"], acc=name), args.update_search)
        stored += [path, C.retrieval_note(path)]
        table = platform_table(path, name)
        wanted = [entry["identifier"], *entry["genes"].values()]
        C.require(set(wanted) <= set(table.columns), f"{name}: the platform table lacks {sorted(set(wanted) - set(table.columns))}")
        if "accession" in entry["genes"]:
            # A cell may name several accessions; they are split and read without their version as the mapping stage reads them.
            accessions.update(value.split(".")[0] for cell in table[entry["genes"]["accession"]] for value in separator.split(cell.strip()) if value and value != "---")
        described = (table[list(entry["genes"].values())].apply(lambda column: column.str.strip()) != "").any(axis=1)
        C.save_table(table, PLATFORMS / f"{name}.parquet")
        outputs.append(PLATFORMS / f"{name}.parquet")
        sources.append(dict(kind="geo_platform_table", name=name, version=f"retrieved {retrieved['retrieved_at_utc'][:10]}",
                            built_from="columns " + ", ".join(f"{kind}: {column}" for kind, column in entry["genes"].items()), platforms=";".join(platforms),
                            rows=len(table), identifiers=int(table[entry["identifier"]].nunique()), with_a_gene=int(described.sum()),
                            table=(PLATFORMS / f"{name}.parquet").relative_to(C.ROOT).as_posix(), raw_file=path.relative_to(C.ROOT).as_posix()))

    # 3. The annotation packages of the image, and the Entrez genes of the accessions that the platform tables name.
    packaged = users(lambda platform, entry: entry.get("package"))
    versions = C.tool_versions([*packaged, policy["organism_package"], "AnnotationDbi"])
    missing = sorted(name for name, version in versions.items() if not version)
    C.require(not missing, f"The image lacks R packages; rebuild it from docker/Dockerfile: {missing}")
    print(f"[02] Exporting {len(packaged)} annotation packages and {len(accessions):,} accessions", flush=True)
    exported, related, described = export_packages(sorted(packaged), sorted(accessions), policy)
    for name, table in exported.items():
        info = described[name]
        C.require(table["probe_id"].nunique() == info["probes"] and info["version"] == versions[name], f"{name}: the export does not hold every probe of the package")
        C.save_table(table, PACKAGES / f"{name}.parquet")
        outputs.append(PACKAGES / f"{name}.parquet")
        sources.append(dict(kind="bioconductor_package", name=name, version=info["version"],
                            built_from=f"{info['metadata'].get('EGSOURCENAME', 'Entrez Gene')} {info['metadata'].get('EGSOURCEDATE', '')}".strip(),
                            platforms=";".join(packaged[name]), rows=len(table), identifiers=info["probes"],
                            with_a_gene=int(table.loc[table["entrez_id"] != "", "probe_id"].nunique()),
                            table=(PACKAGES / f"{name}.parquet").relative_to(C.ROOT).as_posix(), raw_file=""))
    info = described[policy["organism_package"]]
    C.save_table(related, C.DATA / "accession_genes.parquet")
    outputs.append(C.DATA / "accession_genes.parquet")
    sources.append(dict(kind="accession_reference", name=policy["organism_package"], version=info["version"],
                        built_from=f"{info['metadata'].get('EGSOURCENAME', 'Entrez Gene')} {info['metadata'].get('EGSOURCEDATE', '')}".strip(),
                        platforms=";".join(sorted(platform for platform, entry in assigned.items() if "accession" in entry.get("genes", {}))),
                        rows=len(related), identifiers=len(accessions), with_a_gene=int(related["accession"].nunique()),
                        table=(C.DATA / "accession_genes.parquet").relative_to(C.ROOT).as_posix(), raw_file=""))

    files = sorted(stored)
    stray = sorted(path.name for path in C.RAW.iterdir() if path.name != "SHA256SUMS.txt" and path not in stored)
    C.require(not stray, f"data/raw/annotation holds files that are not recorded sources: {stray}; remove them, or run with --update-search")
    sums = "".join(f"{C.sha(path)}  {path.name}\n" for path in files)
    C.require(args.update_search or sums == listed, f"The stored sources under data/raw/annotation differ from their checksum list; {C.STORED_ANSWERS}")
    listing.write_text(sums, encoding="utf-8")
    sources = pd.DataFrame(sources, columns=SOURCE_COLUMNS)
    C.save_table(sources, C.DATA / "annotation_sources.parquet")
    C.save_review(sources, C.OUT / "annotation_sources_review.csv")
    C.finish(2, start, [C.DATA / "annotation_sources.parquet", C.OUT / "annotation_sources_review.csv", C.RAW / "SHA256SUMS.txt", *files, *outputs],
             method=policy["method"], retrieved_at_utc=datetime.now(timezone.utc).isoformat(), search_mode="updated" if args.update_search else "stored",
             genes=dict(total=len(genes), by_status={str(key): int(value) for key, value in genes["status"].value_counts().items()},
                        with_an_entrez_gene=int((genes["entrez_id"] != "").sum()), with_an_ensembl_gene=int((genes["ensembl_gene_id"] != "").sum())),
             platforms=dict(total=len(needed), with_a_package=sum(map(len, packaged.values())), with_a_platform_table=len(needed) - sum(map(len, packaged.values()))),
             sources={str(kind): len(part) for kind, part in sources.groupby("kind")}, tools=versions, recorded_platforms=assigned,
             protocol=C.CONFIG["protocol"])
    print(f"[02] Complete: {len(genes):,} HGNC genes; {len(needed)} array platforms, {sum(map(len, packaged.values()))} annotated by {len(packaged)} packages and "
          f"{len(needed) - sum(map(len, packaged.values()))} by their GEO platform tables; review sheet: "
          f"{(C.OUT / 'annotation_sources_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for row in sources.itertuples(index=False):
        print(f"[02]   {row.name} ({row.version}): {row.with_a_gene:,} of {row.identifiers:,} identifiers with a gene", flush=True)


if __name__ == "__main__":
    main()
