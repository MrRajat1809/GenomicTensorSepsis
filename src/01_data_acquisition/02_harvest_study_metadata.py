"""Retrieve the study and sample records of every census candidate and normalise them into study, sample and attribute tables."""
import argparse
import csv
import importlib
import io
import json
import re
import shutil
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from urllib.parse import quote

C = importlib.import_module("01_search_repositories")

RAW = C.ROOT / "data/raw/metadata"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
BIOSTUDIES_FILES = "https://www.ebi.ac.uk/biostudies/files/"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
EGA = "https://metadata.ega-archive.org/"
STUDY_COLUMNS = ["source", "accession", "title", "summary", "design", "study_type", "platforms", "pubmed_ids",
                 "superseries_of", "subseries_of", "reanalysis_of", "supplementary_files", "submission_date", "last_update_date",
                 "recorded_samples", "harvested_samples", "raw_bytes", "harvest_status"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "title", "sample_type", "source_name", "organism",
                  "taxid", "molecule", "platform", "channel_count", "library_strategy", "library_source",
                  "library_selection", "instrument_model", "description", "data_processing", "raw_files", "relations",
                  "last_update_date"]
ATTRIBUTE_COLUMNS = ["source", "study_accession", "sample_accession", "channel", "position", "kind", "key", "value"]
SDRF_ATTRIBUTE = re.compile(r"(characteristics|factor ?value)\s*\[(.+)\]", re.IGNORECASE)
SDRF_FIELDS = {
    "source_name": ("source name",), "organism": ("characteristics[organism]",), "platform": ("array design ref",), "description": ("description",),
    "library_strategy": ("comment[library_strategy]",), "library_source": ("comment[library_source]",),
    "library_selection": ("comment[library_selection]",), "instrument_model": ("comment[instrument_model]",),
    "raw_files": ("array data file", "comment[fastq_uri]", "comment[submitted_file_name]"),
    "relations": ("comment[ena_sample]", "comment[biosd_sample]", "comment[ena_run]"),
}
# Every stored record that the harvest reads, so that the checksum list covers exactly these.
RETRIEVED = set()


def retrieve(relative, url, params=None, update=False):
    """Read one stored repository record under data/raw/metadata, or ask the repository again and store its answer unmodified."""
    path = RAW / relative
    if update:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(C.fetch(url, params))
        temporary.replace(path)
    C.require(path.is_file(), f"No stored record data/raw/metadata/{relative}; {C.STORED_ANSWERS}")
    RETRIEVED.add(path)
    return path


def stored_checksums():
    """The checksum list of the stored records, which a run without --update-search reads and leaves as it is."""
    listing = RAW / "SHA256SUMS.txt"
    C.require(listing.is_file(), f"No checksum list of the stored records at data/raw/metadata/SHA256SUMS.txt; {C.STORED_ANSWERS}")
    return {name: digest for digest, name in (line.split("  ", 1) for line in listing.read_text(encoding="utf-8").splitlines())}


def soft_records(path, entity):
    """Parse GEO SOFT text into {accession: {field: [values]}}; fields repeat, so every value is kept in order."""
    records, current = {}, None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("^"):
            kind, _, accession = line[1:].partition(" = ")
            current = records.setdefault(accession.strip(), {}) if kind.strip() == entity else None
        elif line.startswith("!") and current is not None:
            key, _, value = line[1:].partition(" = ")
            current.setdefault(key.strip(), []).append(value.strip())
    return records


def first(record, key):
    return record.get(key, [""])[0]


def related(record, key, label):
    return C.joined(value.split(":", 1)[1] for value in record.get(key, []) if value.startswith(label + ":"))


def harvest_geo(accessions, update, studies, samples, attributes):
    for index, accession in enumerate(accessions, 1):
        query = dict(acc=accession, form="text", view="brief")
        series = soft_records(retrieve(f"geo/{accession}.series.soft", GEO, dict(query, targ="self"), update), "SERIES")
        records = soft_records(retrieve(f"geo/{accession}.samples.soft", GEO, dict(query, targ="gsm"), update), "SAMPLE")
        C.require(list(series) == [accession], f"{accession}: GEO returned an unexpected series record")
        series = series[accession]
        listed = series.get("Series_sample_id", [])
        # A third-party reanalysis series owns no samples; it points at the samples of the series it reanalyses.
        C.require(set(listed) == set(records) and (records or related(series, "Series_relation", "Reanalysis of")),
                  f"{accession}: sample records do not match the series sample list")
        studies.append(dict(source="geo", accession=accession, title=first(series, "Series_title"),
            summary=" ".join(series.get("Series_summary", [])), design=" ".join(series.get("Series_overall_design", [])),
            study_type=C.joined(series.get("Series_type", [])), platforms=C.joined(series.get("Series_platform_id", [])),
            pubmed_ids=C.joined(series.get("Series_pubmed_id", [])), superseries_of=related(series, "Series_relation", "SuperSeries of"),
            subseries_of=related(series, "Series_relation", "SubSeries of"), reanalysis_of=related(series, "Series_relation", "Reanalysis of"),
            supplementary_files=C.joined(series.get("Series_supplementary_file", [])),
            submission_date=first(series, "Series_submission_date"), last_update_date=first(series, "Series_last_update_date"),
            recorded_samples=len(listed), harvested_samples=len(records), raw_bytes=0,
            harvest_status="harvested" if records else "reanalysis_without_own_samples"))
        for sample, item in records.items():
            channels = int(first(item, "Sample_channel_count") or 1)
            across = lambda field: C.joined(value for channel in range(1, channels + 1) for value in item.get(f"{field}_ch{channel}", []))
            samples.append(dict(source="geo", study_accession=accession, sample_accession=sample, title=first(item, "Sample_title"),
                sample_type=first(item, "Sample_type"), source_name=across("Sample_source_name"), organism=across("Sample_organism"),
                taxid=across("Sample_taxid"), molecule=across("Sample_molecule"), platform=first(item, "Sample_platform_id"),
                channel_count=channels, library_strategy=first(item, "Sample_library_strategy"),
                library_source=first(item, "Sample_library_source"), library_selection=first(item, "Sample_library_selection"),
                instrument_model=first(item, "Sample_instrument_model"), description=" ".join(item.get("Sample_description", [])),
                data_processing=" ".join(item.get("Sample_data_processing", [])),
                raw_files=C.joined(value for key, values in item.items() if key.startswith("Sample_supplementary_file") for value in values),
                relations=C.joined(item.get("Sample_relation", [])), last_update_date=first(item, "Sample_last_update_date")))
            for channel in range(1, channels + 1):
                for position, entry in enumerate(item.get(f"Sample_characteristics_ch{channel}", [])):
                    key, separator, value = entry.partition(":")
                    attributes.append(dict(source="geo", study_accession=accession, sample_accession=sample, channel=channel,
                        position=position, kind="characteristic", key=key.strip() if separator else "",
                        value=value.strip() if separator else entry.strip()))
        if index % 50 == 0 or index == len(accessions):
            print(f"[02] geo: {index:,}/{len(accessions):,} studies", flush=True)


def declared_files(node):
    """Yield every file entry a BioStudies record declares, wherever it is nested."""
    if isinstance(node, list):
        for item in node:
            yield from declared_files(item)
    elif isinstance(node, dict):
        if "path" in node:
            yield node
        for value in node.values():
            if isinstance(value, (list, dict)):
                yield from declared_files(value)


def harvest_arrayexpress(candidates, update, studies, samples, attributes):
    for accession, title, description, platform in candidates:
        declared = [item["path"] for item in declared_files(C.read_snapshot("arrayexpress", f"study.{accession}.json"))]
        tables = sorted(name for name in declared if name.lower().endswith(".sdrf.txt"))
        # The sample and description tables are stored by their names alone, which stay apart only while each begins with its own accession.
        foreign = sorted(name for name in declared if name.lower().endswith((".idf.txt", ".sdrf.txt")) and not name.startswith(accession + "."))
        C.require(not foreign, f"{accession}: description or sample tables not named after the study: {foreign}")
        for name in sorted(name for name in declared if name.lower().endswith(".idf.txt")):
            retrieve(f"arrayexpress/{name}", BIOSTUDIES_FILES + f"{accession}/{quote(name)}", update=update)
        sources = {}
        for name in tables:
            content = retrieve(f"arrayexpress/{name}", BIOSTUDIES_FILES + f"{accession}/{quote(name)}", update=update).read_text(encoding="utf-8", errors="replace")
            rows = list(csv.reader(io.StringIO(content), delimiter="\t"))
            header = [column.strip() for column in rows[0]]
            count = lambda column: len({row[header.index(column)].strip() for row in rows[1:] if len(row) > header.index(column)} - {""}) if column in header else 0
            assay = next((column for column in ("Assay Name", "Hybridization Name") if column in header), None)
            origin = next((column for column in ("Source Name", "Extract Name", "Labeled Extract Name") if column in header), None)
            C.require(assay or origin, f"{accession}: {name} has no assay, hybridization or source column")
            # One record per assay, the unit that carries one expression profile; a source can be assayed more than once.
            # Fall back to the source when a submitter gave every sample of a data matrix the same assay name.
            unit = header.index(assay if assay and count(assay) * max(count("Label"), 1) >= count(origin) else origin)
            for row in rows[1:]:
                row = row + [""] * (len(header) - len(row))
                if not row[unit].strip():
                    continue
                source = sources.setdefault(row[unit].strip(), {})
                for position, (column, value) in enumerate(zip(header, row)):
                    if value.strip():
                        source.setdefault((position, column), []).append(value.strip())
        studies.append(dict(source="arrayexpress", accession=accession, title=title, summary=description, design="", study_type="",
            platforms=platform, pubmed_ids="", superseries_of="", subseries_of="", reanalysis_of="", supplementary_files=C.joined(declared),
            submission_date="", last_update_date="", recorded_samples=len(sources), harvested_samples=len(sources), raw_bytes=0,
            harvest_status="harvested" if sources else "no_sample_table"))
        for sample, columns in sources.items():
            def field(name):
                return C.joined(value for (_, column), values in columns.items()
                                if re.sub(r"\s*\[", "[", column.lower()) in SDRF_FIELDS[name] for value in values)
            samples.append(dict(source="arrayexpress", study_accession=accession, sample_accession=sample, title=sample, sample_type="",
                source_name=field("source_name"), organism=field("organism"), taxid="", molecule="", platform=field("platform"), channel_count=1,
                library_strategy=field("library_strategy"), library_source=field("library_source"),
                library_selection=field("library_selection"), instrument_model=field("instrument_model"),
                description=field("description"), data_processing="", raw_files=field("raw_files"), relations=field("relations"),
                last_update_date=""))
            for (position, column), values in columns.items():
                matched = SDRF_ATTRIBUTE.fullmatch(column)
                if matched:
                    kind = "characteristic" if matched.group(1).lower() == "characteristics" else "factor"
                    attributes.append(dict(source="arrayexpress", study_accession=accession, sample_accession=sample, channel=1,
                        position=position, kind=kind, key=matched.group(2).strip(), value=C.joined(values)))
    print(f"[02] arrayexpress: {len(candidates):,} studies", flush=True)


def harvest_ena(accessions, update, studies, samples, attributes):
    policy, search, census = C.CONFIG["harvest"], C.CONFIG["search"]["ena"], {}
    for path in sorted((C.CENSUS / "ena").glob("runs.*.json")):
        for run in C.read_snapshot("ena", path.name):
            census.setdefault(run["study_accession"], set()).add(run["sample_accession"])
    for index, accession in enumerate(accessions, 1):
        query = (f'study_accession="{accession}" AND tax_eq({search["tax_id"]}) AND library_source="{search["library_source"]}" '
                 f'AND library_strategy="{search["library_strategy"]}"')
        parameters = dict(result="read_run", query=query, fields=policy["ena_run_fields"], format="json", limit=0)
        content = retrieve(f"ena/{accession}.runs.json", C.ENA, parameters, update).read_bytes()
        runs = json.loads(content) if content.strip() else []
        wanted = sorted({run["sample_accession"] for run in runs if run["sample_accession"]})
        C.require(set(wanted) == census.get(accession, set()), f"{accession}: runs changed since the census; rerun stage 01")
        by_sample, found = {}, {}
        for run in runs:
            by_sample.setdefault(run["sample_accession"], []).append(run)
        for start in range(0, len(wanted), policy["ena_sample_batch"]):
            batch = wanted[start:start + policy["ena_sample_batch"]]
            path = retrieve(f"ena/{accession}.samples.{start // policy['ena_sample_batch']:03d}.xml", ENA_XML + ",".join(batch), update=update)
            for node in ET.fromstring(path.read_bytes()).iter("SAMPLE"):
                for key in {node.get("accession"), *(item.text for item in node.iter("PRIMARY_ID")), *(item.text for item in node.iter("EXTERNAL_ID"))}:
                    found[key] = node
        size = sum(int(part) for run in runs for part in (run.get("fastq_bytes") or "").split(";") if part)
        studies.append(dict(source="ena", accession=accession, title="", summary="", design="", study_type="", platforms="",
            pubmed_ids="", superseries_of="", subseries_of="", reanalysis_of="", supplementary_files="", submission_date="", last_update_date="",
            recorded_samples=len(wanted), harvested_samples=len(set(wanted) & set(found)), raw_bytes=size,
            harvest_status="harvested" if size else "harvested_no_public_reads"))
        for sample in wanted:
            node, sample_runs = found.get(sample), by_sample[sample]
            if node is None:
                continue
            samples.append(dict(source="ena", study_accession=accession, sample_accession=sample, title=node.findtext("TITLE") or "",
                sample_type="SRA", source_name="", organism=node.findtext("SAMPLE_NAME/SCIENTIFIC_NAME") or "",
                taxid=node.findtext("SAMPLE_NAME/TAXON_ID") or "", molecule="", platform="", channel_count=1,
                library_strategy=search["library_strategy"], library_source=search["library_source"],
                library_selection=C.joined(run.get("library_selection", "") for run in sample_runs),
                instrument_model=C.joined(run.get("instrument_model", "") for run in sample_runs),
                description=node.findtext("DESCRIPTION") or "", data_processing="",
                raw_files=C.joined(run["run_accession"] for run in sample_runs), relations=node.get("alias") or "", last_update_date=""))
            for position, item in enumerate(node.iter("SAMPLE_ATTRIBUTE")):
                attributes.append(dict(source="ena", study_accession=accession, sample_accession=sample, channel=1, position=position,
                    kind="characteristic", key=(item.findtext("TAG") or "").strip(), value=(item.findtext("VALUE") or "").strip()))
        if index % 20 == 0 or index == len(accessions):
            print(f"[02] ena: {index:,}/{len(accessions):,} studies", flush=True)


def harvest_catalogue(candidates, update, studies):
    """Controlled-access holdings expose study records only; EGA also reports dataset sizes."""
    for source, accession, title, description in candidates:
        row = dict(source=source, accession=accession, title=title, summary=description, design="", study_type="", platforms="",
                   pubmed_ids="", superseries_of="", subseries_of="", reanalysis_of="", supplementary_files="", submission_date="", last_update_date="",
                   recorded_samples=0, harvested_samples=0, raw_bytes=0, harvest_status="catalogue_only")
        endpoint = {"EGAS": f"studies/{accession}/datasets", "EGAD": f"datasets/{accession}"}.get(accession[:4])
        if source == "ega" and endpoint:
            content = json.loads(retrieve(f"ega/{accession}.json", EGA + endpoint, update=update).read_bytes())
            datasets = content if isinstance(content, list) else [content]
            row.update(recorded_samples=sum(int(item.get("num_samples") or 0) for item in datasets),
                       platforms=C.joined(name for item in datasets for name in item.get("technologies") or []),
                       study_type=C.joined(name for item in datasets for name in item.get("dataset_types") or []),
                       harvest_status="study_record_only")
        studies.append(row)
    print(f"[02] controlled-access catalogue: {len(candidates):,} records", flush=True)


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--update-search", action="store_true",
                        help="Ask the repositories again for every record instead of reading the stored records under data/raw/metadata")
    args = parser.parse_args()
    print(f"[02] Harvesting study and sample records v{C.VERSION}", flush=True)
    start = C.begin(2)
    candidates = pd.read_parquet(C.DATA / "candidates.parquet")
    listed = None if args.update_search else stored_checksums()
    if RAW.exists() and args.update_search:
        shutil.rmtree(RAW)
    studies, samples, attributes = [], [], []
    part = lambda source: candidates[candidates["source"] == source]
    harvest_geo(part("geo")["accession"].tolist(), args.update_search, studies, samples, attributes)
    harvest_arrayexpress(list(part("arrayexpress")[["accession", "title", "description", "platform"]].itertuples(index=False, name=None)),
                         args.update_search, studies, samples, attributes)
    harvest_ena(part("ena")["accession"].tolist(), args.update_search, studies, samples, attributes)
    catalogued = candidates[candidates["access"] != "open"]
    harvest_catalogue(list(catalogued[["source", "accession", "title", "description"]].itertuples(index=False, name=None)),
                      args.update_search, studies)
    studies = pd.DataFrame(studies, columns=STUDY_COLUMNS)
    samples = pd.DataFrame(samples, columns=SAMPLE_COLUMNS)
    attributes = pd.DataFrame(attributes, columns=ATTRIBUTE_COLUMNS)
    order = candidates[["source", "accession"]].merge(studies, on=["source", "accession"], how="left", validate="one_to_one")
    C.require(len(studies) == len(candidates) and order["harvest_status"].notna().all(), "Every candidate needs exactly one study record")
    studies = order[STUDY_COLUMNS].copy()
    # ENA and catalogue records carry their title and description in the census; fill them from there.
    for column, origin in (("title", "title"), ("summary", "description")):
        studies[column] = studies[column].where(studies[column] != "", candidates[origin].to_numpy())
    C.require(not samples.duplicated(["source", "study_accession", "sample_accession"]).any(), "Repeated sample within a study")
    keys = set(zip(samples["source"], samples["study_accession"], samples["sample_accession"]))
    C.require(set(zip(attributes["source"], attributes["study_accession"], attributes["sample_accession"])) <= keys,
              "Attribute rows without a sample record")
    C.save_table(studies, C.DATA / "studies.parquet")
    C.save_table(samples, C.DATA / "samples.parquet")
    C.save_table(attributes, C.DATA / "sample_attributes.parquet")
    stray = sorted(path.relative_to(RAW).as_posix() for path in RAW.rglob("*") if path.is_file() and path != RAW / "SHA256SUMS.txt" and path not in RETRIEVED)
    C.require(not stray, f"data/raw/metadata holds {len(stray)} files that the harvest does not read: {stray[:20]}; remove them, or run with --update-search")
    files = sorted(RETRIEVED)
    sums = {path.relative_to(RAW).as_posix(): C.sha(path) for path in files}
    C.require(listed is None or sums == listed, f"The stored records under data/raw/metadata differ from their checksum list; {C.STORED_ANSWERS}")
    (RAW / "SHA256SUMS.txt").write_text("".join(f"{digest}  {name}\n" for name, digest in sums.items()), encoding="utf-8")
    summary = []
    for source in C.CONFIG["deduplication"]["source_priority"]:
        chosen = studies[studies["source"] == source]
        summary.append(dict(source=source, studies=len(chosen), **C.counts(chosen, "harvest_status"),
            samples=int((samples["source"] == source).sum()), attributes=int((attributes["source"] == source).sum()),
            unique_samples=int(samples.loc[samples["source"] == source, "sample_accession"].nunique())))
    counted = candidates[["source", "accession", "n_samples"]].merge(studies[["source", "accession", "harvested_samples"]], on=["source", "accession"])
    differing = counted[(counted["n_samples"] != counted["harvested_samples"]) & counted["source"].isin(["geo", "arrayexpress", "ena"])]
    C.finish(2, start, [C.DATA / "studies.parquet", C.DATA / "samples.parquet", C.DATA / "sample_attributes.parquet", RAW / "SHA256SUMS.txt"],
             harvested_at_utc=datetime.now(timezone.utc).isoformat(), search_mode="updated" if args.update_search else "stored", raw_files=len(files),
             summary=summary, sample_count_differs_from_census=json.loads(differing.to_json(orient="records")),
             open_fastq_bytes=int(studies["raw_bytes"].sum()), protocol=C.CONFIG["protocol"])
    print(f"[02] Complete: {len(studies):,} studies; {len(samples):,} sample records; {len(attributes):,} attribute values", flush=True)


if __name__ == "__main__":
    main()
