"""Catalogue every recorded field of every included sample with expression values, verbatim, list the fields of each study platform, and list the records without values with their reason."""
import hashlib
import importlib.metadata
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ACQUIRED = ROOT / "data/processed/acquisition"
EXPRESSION = ROOT / "data/processed/expression"
ANNOTATION = ROOT / "data/processed/annotation"
ANNOTATION_OUT = ROOT / "outputs/gene_annotation"
DATA = ROOT / "data/processed/cohort"
OUT = ROOT / "outputs/cohort_curation"
CONFIG = json.loads((HERE / "cohort_config.json").read_text(encoding="utf-8"))
REGISTER = json.loads((HERE / CONFIG["register"]).read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
STAGES = {1: "01_catalogue_sample_metadata.py", 2: "02_link_repeated_samples.py", 3: "03_curate_sample_labels.py", 4: "04_flag_sample_quality.py",
          5: "05_build_cohort_tables.py", 6: "06_define_benchmark_splits.py", 7: "07_cohort_master_qc.py"}
METADATA_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "origin", "channel", "position", "key", "field", "value"]
FIELD_COLUMNS = ["source", "study_accession", "platform", "origin", "field", "keys", "samples", "samples_with_value", "distinct_values", "most_frequent_values"]
LEDGER_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "family", "reason", "title", "source_name", "description", "characteristics",
                  "factor_values"]
# The kinds of recorded attribute that the ledger lists, each in its own column.
LEDGER_KINDS = {"characteristic": "characteristics", "factor": "factor_values"}
# The register lists, each sorted by the accession that this path reads from its entries; rules and sepsis_definitions are keyed by accession.
REGISTER_LISTS = {"same_sample": ("a",), "replicates": ("study",), "record_links": ("a", "study"), "patient_namespaces": ("studies", 0),
                  "record_flags": ("study",), "replicate_groups": ("study",), "identity_flags": ("study",)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_path(recorded):
    path = (ROOT / recorded).resolve()
    require(path.is_relative_to(ROOT), f"Path outside project: {recorded}")
    return path


def record(path):
    path = Path(path).resolve()
    return {"sha256": sha(path), "size_bytes": path.stat().st_size}


def bind(paths):
    paths = sorted(set(Path(path).resolve() for path in paths), key=lambda path: path.as_posix())
    return {path.relative_to(ROOT).as_posix(): record(path) for path in paths}


def verify_files(records):
    for name, expected in records.items():
        path = project_path(name)
        require(path.is_file() and record(path) == expected, f"Missing or changed artifact: {name}")


def code_binding(number=max(STAGES)):
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "cohort_config.json", HERE / CONFIG["register"]])


def runtime():
    result = {"python": sys.version.split()[0]}
    for name in ("numpy", "pandas", "pyarrow"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def receipt_path(number):
    return OUT / f"{number:02d}.json"


def verify_stage(number):
    result = read_json(receipt_path(number))
    require(result["version"] == VERSION and result["status"] == "COMPLETE", f"Stage {number} incomplete")
    require(result["code"] == code_binding(number), f"Cohort code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
    verify_files(result["outputs"])
    verify_files(result["dependencies"])
    return result


def begin(number):
    validate_config()
    print(f"[{number:02d}] Verifying inputs and provenance v{VERSION}", flush=True)
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if number < max(STAGES) and (OUT / "qc.json").exists():
        write_json(OUT / "qc.json", dict(version=VERSION, status="STALE", reason=f"Stage {number:02d} is being rerun; rerun the master QC after completion"))
    for previous in range(1, number):
        verify_stage(previous)
    if number > 1:
        verify_files(read_json(receipt_path(1))["external"])
    dependencies = bind([receipt_path(previous) for previous in range(1, number)])
    return {"code": code_binding(number), "dependencies": dependencies}


def accession_in(entry, path):
    for step in path:
        entry = entry[step]
    return entry


def validate_config():
    require(VERSION == "1.0.0", "Unsupported cohort version")
    require(CONFIG["register"] == "curation_register.json" and set(REGISTER) == {"rules", "sepsis_definitions", *REGISTER_LISTS}
            and all(isinstance(REGISTER[name], dict) and list(REGISTER[name]) == sorted(REGISTER[name]) for name in ("rules", "sepsis_definitions"))
            and all(isinstance(REGISTER[name], list) and [accession_in(entry, path) for entry in REGISTER[name]]
                    == sorted(accession_in(entry, path) for entry in REGISTER[name]) for name, path in REGISTER_LISTS.items()),
            "curation_register.json must hold the rule sets and the sepsis definitions keyed by accession, and the lists "
            + ", ".join(REGISTER_LISTS) + ", each sorted by accession")
    catalogue, links = CONFIG["catalogue"], CONFIG["links"]
    require(catalogue["record_fields"] and len(set(catalogue["record_fields"])) == len(catalogue["record_fields"]) and catalogue["examples"] > 0,
            "Invalid catalogue policy")
    require(links["median_tolerance"] >= 0 and links["primary_record"] == "earliest_release", "Invalid link policy")
    identity = links["identity"]
    require(0 < identity["candidate_correlation"] <= 1 and identity["candidate_margin"] >= 0 and 0 < identity["minimum_common_genes"] <= identity["genes_used"]
            and all(str(entry.get("basis", "")).strip() and 0.9 < entry["minimum_correlation"] <= 1 for entry in REGISTER["same_sample"] + REGISTER["replicates"])
            and all(str(rule.get("basis", "")).strip() for rule in REGISTER["record_links"] + REGISTER["patient_namespaces"])
            and all(re.compile(side["pattern"]).groups == 1 for rule in REGISTER["record_links"] for side in (rule["a"], rule["b"]))
            and all(str(rule["b"].get("required", "")).strip() and re.compile(rule["b"]["required"]) for rule in REGISTER["record_links"]),
            "Invalid identity, record-link or patient-namespace entries")
    require(all(re.compile(rule["pattern"]).groups == 1 and rule["flag"] for rule in REGISTER["replicate_groups"]) and CONFIG["quality"]["replicate_flag"]
            and all({"study", "field", "pattern", "flag"} <= set(rule) for rule in REGISTER["record_flags"])
            and all({"study", "samples", "flag", "basis"} <= set(rule) for rule in REGISTER["identity_flags"]), "Invalid replicate or quality-flag entries")
    curation = CONFIG["curation"]
    require(REGISTER["rules"] and {"sample_class", "tissue", "sex", "outcome", "septic_shock", "age_group", "time_unit"} <= set(curation["vocabularies"])
            and 0 < curation["age_groups"]["neonate_up_to_days"] < 365 and curation["age_groups"]["adult_from_years"] > 1, "Invalid curation policy")
    benchmark, vocabularies = CONFIG["benchmark"], curation["vocabularies"]
    require(benchmark["tissues"] and set(benchmark["tissues"]) <= set(vocabularies["tissue"]) - {"not_recorded"}
            and set(benchmark["phases_left_out"]) <= set(vocabularies["sepsis_phase"]) and benchmark["minimum_per_class"] >= 1
            and re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", benchmark["temporal_split_date"]) and set(benchmark["technologies"].values()) == {"array", "rna_seq"}
            and benchmark["tasks"] and len({task["task"] for task in benchmark["tasks"]}) == len(benchmark["tasks"])
            and all(task["variable"] in ("sample_class", "outcome", "septic_shock") and task["positive"] != task["negative"]
                    and {task["positive"], task["negative"]} <= set(vocabularies[task["variable"]]) - {"not_recorded", "infection_sepsis_not_stated"}
                    and task.get("within_class", "sepsis") in vocabularies["sample_class"] for task in benchmark["tasks"]),
            "Invalid benchmark policy")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Cohort code changed during execution")
    verify_files(start["dependencies"])
    result = dict(version=VERSION, status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  runtime=runtime(), outputs=bind(outputs), **start, **detail)
    write_json(receipt_path(number), result)
    return result


def save_table(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_parquet(temporary, index=False, compression="zstd")
    temporary.replace(path)


def save_review(frame, path):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
    temporary.replace(path)


def annotation_binding():
    """Confirm that the annotation stage passed its quality control unchanged and return the artifacts of earlier stages that this stage reads."""
    qc = read_json(ANNOTATION_OUT / "qc.json")
    require(qc["status"] == "PASS" and qc["failures"] == 0, "Annotation quality control has not passed")
    verify_files(qc["code"])
    verify_files(qc["dependencies"])
    # The first annotation script bound the tables of the expression stage and the sample tables of the acquisition stage; they must be unchanged.
    first = read_json(ANNOTATION_OUT / "01.json")
    verify_files(first["external"])
    verify_files(first["outputs"])
    tables = [name for name in first["external"] if name.endswith(".parquet")] + [name for name in first["outputs"] if name.endswith(".parquet")]
    return [ANNOTATION_OUT / "qc.json", *map(project_path, qc["dependencies"]), *map(project_path, tables), ACQUIRED / "included_studies.parquet",
            ACQUIRED / "studies.parquet"]


def sample_records(columns):
    """Return the included sample records that have expression values, in the order of the sample index of the annotation stage; stage 01 lists the
    others in records_without_values.parquet."""
    import pandas as pd

    index = pd.read_parquet(ANNOTATION / "sample_index.parquet", columns=[*columns, "in_matrix"])
    require(index["in_matrix"].dtype == bool, "The sample index must mark every record with or without expression values")
    return index[index["in_matrix"]].drop(columns="in_matrix").reset_index(drop=True)


def accession_order(accession):
    """Order accessions by their letters and then by their number, so that GSM99999 comes before GSM100000."""
    match = re.fullmatch(r"(.*?)([0-9]+)", str(accession))
    return (match.group(1), int(match.group(2)), str(accession)) if match else (str(accession), -1, str(accession))


def field_name(key):
    """Write a key in lower case with single spaces, so that one field keeps one name across the spellings of a study; a value without a key is named so."""
    return re.sub(r"\s+", " ", str(key).replace("_", " ")).strip().lower() or "(no key)"


def main():
    import pandas as pd

    print(f"[01] Cataloguing sample metadata v{VERSION}", flush=True)
    start = begin(1)
    external = annotation_binding()
    policy = CONFIG["catalogue"]
    keys = ["source", "study_accession", "sample_accession", "platform"]
    included = pd.read_parquet(ANNOTATION / "sample_index.parquet", columns=[*keys, "family", "in_matrix", "reason"])
    included["reason"] = included["reason"].fillna("")
    require(included["in_matrix"].dtype == bool and (included["in_matrix"] == (included["reason"] == "")).all(),
            "Every included record without expression values needs the reason that the expression stage recorded, and a record with values none")
    index = sample_records(keys)
    records = pd.read_parquet(ACQUIRED / "samples.parquet").merge(included[keys], on=keys, validate="one_to_one")
    require(len(records) == len(included) and set(policy["record_fields"]) <= set(records.columns),
            "The sample records do not hold every included sample and every recorded field")
    attributes = pd.read_parquet(ACQUIRED / "sample_attributes.parquet").merge(included[keys], on=["source", "study_accession", "sample_accession"],
                                                                               validate="many_to_one")

    parts = []
    for position, name in enumerate(policy["record_fields"]):
        part = records[["source", "study_accession", "sample_accession", "platform", name]].rename(columns={name: "value"})
        part = part[part["value"].fillna("").astype(str).str.strip() != ""]
        parts.append(part.assign(origin="record", channel=0, position=position, key=name, value=part["value"].astype(str)))
    parts.append(attributes.rename(columns={"kind": "origin"})[["source", "study_accession", "sample_accession", "platform", "origin", "channel", "position", "key", "value"]])
    metadata = pd.concat(parts, ignore_index=True)
    metadata["field"] = metadata["key"].map(field_name)
    require(set(zip(metadata["study_accession"], metadata["sample_accession"])) == set(zip(included["study_accession"], included["sample_accession"])),
            "A sample has no recorded field at all")
    metadata = metadata.merge(included[["study_accession", "sample_accession", "in_matrix"]], on=["study_accession", "sample_accession"], validate="many_to_one")
    metadata = metadata.sort_values(["source", "study_accession", "sample_accession", "origin", "channel", "position", "key"], kind="stable")
    absent, metadata = metadata[~metadata["in_matrix"]], metadata[metadata["in_matrix"]][METADATA_COLUMNS].reset_index(drop=True)

    # A record without values leaves the cohort tables; the ledger keeps its reason and every field that the curation reads, in their recorded order.
    unknown = sorted(set(absent["origin"]) - {"record", *LEDGER_KINDS})
    require(not unknown, f"A record without expression values carries attributes of a kind that the ledger does not list: {unknown}; name the kind in "
                         "LEDGER_KINDS of 01_catalogue_sample_metadata.py")
    ledger = included.loc[~included["in_matrix"], [*keys, "family", "reason"]].merge(
        records[["study_accession", "sample_accession", "title", "source_name", "description"]], on=["study_accession", "sample_accession"], validate="one_to_one")
    for name in ("title", "source_name", "description"):
        ledger[name] = ledger[name].fillna("").astype(str)
    for kind, column in LEDGER_KINDS.items():
        part = absent[absent["origin"] == kind]
        named = part["key"].where(part["key"].astype(str).str.strip() != "", "(no key)")
        written = named.where(part["channel"] <= 1, "channel " + part["channel"].astype(str) + " " + named) + ": " + part["value"].astype(str)
        joined = written.groupby([part["study_accession"], part["sample_accession"]], sort=False).agg(" | ".join)
        ledger[column] = [joined.get(key, "") for key in zip(ledger["study_accession"], ledger["sample_accession"])]
    ledger = ledger.sort_values(["source", "study_accession", "platform", "sample_accession"], kind="stable")[LEDGER_COLUMNS].reset_index(drop=True)
    reasons = {str(name): int(count) for name, count in ledger["reason"].value_counts().sort_index().items()}
    require(len(index) + len(ledger) == len(included) == len(records) and sum(reasons.values()) == len(ledger),
            "The included records must be those with expression values and those without")
    platforms = set(zip(index["study_accession"], index["platform"]))
    left = sorted(f"{study}|{platform}" for study, platform in set(zip(ledger["study_accession"], ledger["platform"])) - platforms)

    sizes = index.groupby(["study_accession", "platform"]).size()
    fields = []
    for (source, study, platform, origin, field), part in metadata.groupby(["source", "study_accession", "platform", "origin", "field"], sort=True):
        counts = part["value"].value_counts()
        fields.append(dict(source=source, study_accession=study, platform=platform, origin=origin, field=field, keys=" | ".join(sorted(set(part["key"]))),
                           samples=int(sizes[(study, platform)]), samples_with_value=int(part["sample_accession"].nunique()), distinct_values=len(counts),
                           most_frequent_values="; ".join(f"{str(value)[:80]} = {count}" for value, count in counts.head(policy["examples"]).items())))
    fields = pd.DataFrame(fields, columns=FIELD_COLUMNS)
    save_table(metadata, DATA / "sample_metadata.parquet")
    save_table(fields, DATA / "metadata_fields.parquet")
    save_table(ledger, DATA / "records_without_values.parquet")
    save_review(fields, OUT / "metadata_fields_review.csv")
    described = fields[fields["origin"] != "record"]
    finish(1, start, [DATA / "sample_metadata.parquet", DATA / "metadata_fields.parquet", DATA / "records_without_values.parquet",
                      OUT / "metadata_fields_review.csv"], external=bind(external),
           records=dict(included=len(included), with_values=len(index), without_values=dict(total=len(ledger), by_reason=reasons)),
           samples=len(index), study_platforms=len(sizes), study_platforms_without_any_values=left,
           values=dict(total=len(metadata), by_origin={str(key): int(value) for key, value in metadata["origin"].value_counts().sort_index().items()}),
           fields=dict(distinct=int(metadata["field"].nunique()), distinct_among_characteristics_and_factors=int(described["field"].nunique()),
                       study_platforms_with_two_or_fewer=int((described.groupby(["study_accession", "platform"])["field"].nunique() <= 2).sum())),
           record_fields=policy["record_fields"], protocol=CONFIG["protocol"])
    print(f"[01] Complete: {len(metadata):,} recorded values of {len(index):,} sample records with expression values in {len(sizes)} study platforms; "
          f"{int(described['field'].nunique()):,} distinct characteristic and factor fields; {len(ledger)} of {len(included):,} included records have no "
          f"values and are listed with their reason ({', '.join(f'{name} {count}' for name, count in reasons.items())}); review sheet: "
          f"{(OUT / 'metadata_fields_review.csv').relative_to(ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
