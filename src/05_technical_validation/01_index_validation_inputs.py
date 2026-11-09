"""Index what the technical validation can test: the gene matrix and column of every sample, and the genes that the tests read."""
import hashlib
import importlib.metadata
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ANNOTATION = ROOT / "data/processed/annotation"
COHORT = ROOT / "data/processed/cohort"
ANNOTATION_OUT = ROOT / "outputs/gene_annotation"
COHORT_OUT = ROOT / "outputs/cohort_curation"
DATA = ROOT / "data/processed/validation"
OUT = ROOT / "outputs/technical_validation"
CONFIG = json.loads((HERE / "validation_config.json").read_text(encoding="utf-8"))
VERSION = CONFIG["version"]
STAGES = {1: "01_index_validation_inputs.py", 2: "02_check_recorded_sex.py", 3: "03_compare_repeated_deposits.py", 4: "04_check_known_biology.py",
          5: "05_validation_master_qc.py"}
KEYS = ["study_accession", "sample_accession"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "cohort_group", "sample_id", "is_primary_record", "records_of_sample", "linked_by",
                  "sample_class", "sepsis_phase", "tissue", "sex", "patient_key", "one_per_patient", "time_order", "any_flag", "gene_matrix", "value_scale"]
MATRIX_COLUMNS = ["source", "study_accession", "platform", "gene_matrix", "value_scale", "comparable_scale", "samples", "genes", "sex_recorded", "female", "male",
                  "repeated_records", "sepsis", "healthy_control", "tissues"]
GENE_COLUMNS = ["test", "role", "published_symbol", "hgnc_id", "symbol", "resolved_by", "matrices_with_gene"]


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
    return bind([HERE / name for stage, name in STAGES.items() if stage <= number] + [HERE / "validation_config.json"])


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
    require(result["code"] == code_binding(number), f"Validation code/config changed; rerun the earliest changed stage and its dependants (stage {number})")
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
    # Every test reads the index of stage 01 only; the tests do not build on one another.
    for previous in ([1] if 1 < number < max(STAGES) else range(1, number)):
        verify_stage(previous)
    if number > 1:
        verify_files(read_json(receipt_path(1))["external"])
    dependencies = bind([receipt_path(previous) for previous in ([1] if 1 < number < max(STAGES) else range(1, number))])
    return {"code": code_binding(number), "dependencies": dependencies}


def validate_config():
    require(VERSION == "1.0.0", "Unsupported validation version")
    values, sex, biology = CONFIG["values"], CONFIG["sex"], CONFIG["biology"]
    require(values["linear_offset"] > 0 and values["log_scales"] and not set(values["log_scales"]) & set(values["linear_scales"]), "Invalid value policy")
    require(sex["female_gene"] and sex["male_genes"] and sex["minimum_samples"] >= 4 and sex["minimum_separation"] > 0 and sex["minimum_group"] >= 1
            and 0 <= sex["minimum_group_share"] < 0.5 and 0 <= sex["ambiguous_zone"] < 1 and sex["minimum_compared"] >= 1 and 0 < sex["unreliable_below"] <= 1
            and 0 < sex["component_separation"] <= sex["minimum_separation"] and sex["minimum_samples_per_tissue"] >= sex["minimum_samples"],
            "Invalid sex policy")
    require(biology["minimum_group"] >= 3 and biology["score"]["up"] and biology["score"]["down"] and biology["minimum_score_genes"] >= 1
            and biology["minimum_cells"] >= 2 and set(biology["cell_markers"]) <= set(biology["cell_fractions"]) and biology["sepsis_phases"]
            and biology["bootstrap"]["replicates"] >= 200 and 0 < biology["bootstrap"]["level"] < 1
            and not set(biology["score_derivation"]["discovery"]) & set(biology["score_derivation"]["validation"]), "Invalid biology policy")
    require(all(0 < limit < 1 for limits in (*biology["report"].values(), *CONFIG["repeats"]["report"].values()) for limit in limits),
            "Invalid report limits")


def finish(number, start, outputs, **detail):
    require(start["code"] == code_binding(number), "Validation code changed during execution")
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


def earlier_binding():
    """Confirm that the cohort and annotation stages passed their quality control unchanged and return the artifacts of earlier stages that this stage reads."""
    paths = [ANNOTATION / "genes.parquet"]
    for folder, name in ((COHORT_OUT, "Cohort"), (ANNOTATION_OUT, "Annotation")):
        qc = read_json(folder / "qc.json")
        require(qc["status"] == "PASS" and qc["failures"] == 0, f"{name} quality control has not passed")
        verify_files(qc["code"])
        verify_files(qc["dependencies"])
        paths += [folder / "qc.json", *map(project_path, qc["dependencies"])]
    built = read_json(COHORT_OUT / "05.json")["outputs"]
    matrices = {name: value for name, value in read_json(ANNOTATION_OUT / "04.json")["outputs"].items()
                if "/gene_matrices/" in name or name.endswith("gene_matrix_summary.parquet")}
    verify_files(built)
    verify_files(matrices)
    return paths + [project_path(name) for name in [*built, *matrices] if name.endswith(".parquet")]


def resolve_genes(symbols, genes):
    """Find the HGNC gene of every published symbol: the approved symbol of that name, or else the one gene that carried the name before."""
    approved = dict(zip(genes["symbol"], genes["hgnc_id"]))
    current = dict(zip(genes["hgnc_id"], genes["symbol"]))
    earlier = {}
    for identifier, names in zip(genes["hgnc_id"], genes["prev_symbol"].fillna("")):
        for name in filter(None, names.split("|")):
            earlier.setdefault(name, []).append(identifier)
    found = {}
    for symbol in symbols:
        if symbol in approved:
            found[symbol] = (approved[symbol], symbol, "approved symbol")
        else:
            require(len(earlier.get(symbol, [])) == 1, f"A published gene symbol does not name one HGNC gene: {symbol}")
            found[symbol] = (earlier[symbol][0], current[earlier[symbol][0]], "previous symbol")
    return found


def read_genes(path, identifiers):
    """Read the rows of the named genes from one gene matrix: genes in rows, by HGNC identifier, and samples in columns."""
    import pandas as pd

    frame = pd.read_parquet(path, filters=[("hgnc_id", "in", sorted(identifiers))])
    return frame.drop(columns="symbol").set_index("hgnc_id").astype("float64")


def log_values(frame, scale):
    """Return the values on a log2 scale: as they are where the matrix holds logarithms, and log2 of the value plus an offset where it is linear."""
    import numpy as np

    policy = CONFIG["values"]
    if scale in policy["log_scales"]:
        return frame
    require(scale in policy["linear_scales"], f"Values of this scale cannot be compared: {scale}")
    return np.log2(frame.clip(lower=0) + policy["linear_offset"])


def accession_order(accession):
    """Order accessions by their letters and then by their number, so that GSM99999 comes before GSM100000."""
    match = re.fullmatch(r"(.*?)([0-9]+)", str(accession))
    return (match.group(1), int(match.group(2)), str(accession)) if match else (str(accession), -1, str(accession))


def at_least(name, limit):
    """The name under which a receipt counts the values at or beyond a report limit, such as at_least_0_95."""
    return f"{name}_{str(limit).replace('.', '_')}"


def comparable(scale):
    return scale in CONFIG["values"]["log_scales"] or scale in CONFIG["values"]["linear_scales"]


def auroc(cases, others):
    """Area under the ROC curve of a score: the share of case and other pairs in which the case has the higher score, ties counted as half."""
    import numpy as np
    import pandas as pd

    cases, others = np.asarray(cases, dtype=float), np.asarray(others, dtype=float)
    ranks = pd.Series(np.concatenate([cases, others])).rank().to_numpy()
    return float((ranks[:len(cases)].sum() - len(cases) * (len(cases) + 1) / 2) / (len(cases) * len(others)))


def gene_sets():
    """Return every published symbol that the tests read, with the test and the role it has in it."""
    sex, biology = CONFIG["sex"], CONFIG["biology"]
    rows = [("sex", "female", sex["female_gene"]), *(("sex", "male", symbol) for symbol in sex["male_genes"]),
            *(("sepsis_score", "higher_in_sepsis", symbol) for symbol in biology["score"]["up"]),
            *(("sepsis_score", "lower_in_sepsis", symbol) for symbol in biology["score"]["down"]),
            *(("single_gene", "lower_in_sepsis", symbol) for symbol in biology["lower_in_sepsis"]),
            *(("cell_marker", cell, symbol) for cell, symbols in biology["cell_markers"].items() for symbol in symbols)]
    return rows


def main():
    import pandas as pd
    import pyarrow.parquet as pq

    print(f"[01] Indexing validation inputs v{VERSION}", flush=True)
    start = begin(1)
    external = earlier_binding()
    samples = pd.read_parquet(COHORT / "cohort_samples.parquet")
    summary = pd.read_parquet(ANNOTATION / "gene_matrix_summary.parquet", columns=["accession", "platform", "value_scale", "genes", "gene_matrix"])
    summary = summary.rename(columns={"accession": "study_accession"})
    held = samples.merge(summary, on=["study_accession", "platform"], how="left", validate="many_to_one")
    require(held["gene_matrix"].notna().all() and (held["gene_matrix"] != "").all(), "A sample of the cohort table has no gene matrix")

    present, rows = {}, []
    for (source, study, platform, path, scale, genes), part in held.groupby(["source", "study_accession", "platform", "gene_matrix", "value_scale", "genes"], sort=True):
        schema = pq.read_schema(project_path(path))
        columns = set(schema.names) - {"hgnc_id", "symbol"}
        require(columns == set(part["sample_accession"]), f"The gene matrix of {study} {platform} does not hold exactly its samples")
        present[(study, platform)] = set(pd.read_parquet(project_path(path), columns=["hgnc_id"])["hgnc_id"])
        require(len(present[(study, platform)]) == genes, f"The gene matrix of {study} {platform} does not hold the recorded number of genes")
        rows.append(dict(source=source, study_accession=study, platform=platform, gene_matrix=path, value_scale=scale, comparable_scale=comparable(scale),
                         samples=len(part), genes=int(genes), sex_recorded=int(part["sex"].isin(["female", "male"]).sum()), female=int((part["sex"] == "female").sum()),
                         male=int((part["sex"] == "male").sum()), repeated_records=int((part["records_of_sample"] > 1).sum()),
                         sepsis=int((part["sample_class"] == "sepsis").sum()), healthy_control=int((part["sample_class"] == "healthy_control").sum()),
                         tissues="; ".join(f"{name}={count}" for name, count in part["tissue"].value_counts().items())))
    matrices = pd.DataFrame(rows, columns=MATRIX_COLUMNS)

    reference = pd.read_parquet(ANNOTATION / "genes.parquet", columns=["hgnc_id", "symbol", "prev_symbol"])
    named = gene_sets()
    found = resolve_genes(sorted({symbol for _, _, symbol in named}), reference)
    genes = pd.DataFrame([dict(test=test, role=role, published_symbol=symbol, hgnc_id=found[symbol][0], symbol=found[symbol][1], resolved_by=found[symbol][2],
                               matrices_with_gene=sum(found[symbol][0] in held_genes for held_genes in present.values())) for test, role, symbol in named],
                         columns=GENE_COLUMNS)
    require(not genes.duplicated(["test", "role", "hgnc_id"]).any(), "A gene is named twice in one role")
    held = held[SAMPLE_COLUMNS].sort_values(["source", "study_accession", "platform", "sample_accession"], kind="stable").reset_index(drop=True)
    save_table(held, DATA / "validation_samples.parquet")
    save_table(matrices, DATA / "validation_matrices.parquet")
    save_table(genes, DATA / "validation_genes.parquet")
    save_review(genes, OUT / "validation_genes_review.csv")
    finish(1, start, [DATA / "validation_samples.parquet", DATA / "validation_matrices.parquet", DATA / "validation_genes.parquet",
                      OUT / "validation_genes_review.csv"], external=bind(external),
           samples=dict(records=len(held), with_recorded_sex=int(held["sex"].isin(["female", "male"]).sum()),
                        repeated_records=int((held["records_of_sample"] > 1).sum())),
           matrices=dict(total=len(matrices), of_comparable_scale=int(matrices["comparable_scale"].sum()),
                         by_scale={str(name): int(count) for name, count in matrices["value_scale"].value_counts().sort_index().items()}),
           genes=dict(named=int(genes["hgnc_id"].nunique()), resolved_by_previous_symbol=sorted(set(genes.loc[genes["resolved_by"] == "previous symbol",
                                                                                                       "published_symbol"]))),
           protocol=CONFIG["protocol"])
    print(f"[01] Complete: {len(held):,} samples are columns of {len(matrices)} gene matrices ({int(matrices['comparable_scale'].sum())} of a comparable scale); "
          f"{int(genes['hgnc_id'].nunique())} genes named by the tests; review sheet: {(OUT / 'validation_genes_review.csv').relative_to(ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
