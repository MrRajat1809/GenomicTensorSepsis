"""Verify that every included sample has one recorded outcome and that every expression matrix agrees with the tables that describe it."""
import gzip
import importlib
import shutil
import subprocess
import traceback
from datetime import datetime, timezone

import numpy as np
import pandas as pd

C = importlib.import_module("01_inventory_expression_files")

KEYWORDS = ("per_sample_files", "ncbi_counts", "series_matrix", "none")
TABLES = ("table", "workbook_sheet")
# For each processing stage: its table of study platforms, its table of samples, and the column that counts the rows of a matrix.
STAGES = {4: ("affymetrix_studies", "affymetrix_arrays", "features"), 5: ("single_channel_studies", "single_channel_arrays", "probes"),
          6: ("intensity_table_studies", "intensity_table_arrays", "probes"), 7: ("count_table_studies", "count_table_samples", "genes"),
          8: ("author_matrix_studies", "author_matrix_samples", "features")}
VALUE_COLUMNS = ["source", "accession", "platform", "series_matrix", "status", "samples_compared", "own_column_most_similar", "shared_rows", "rows_used",
                 "median_own_correlation", "lowest_own_correlation", "smallest_margin", "exception", "not_most_similar"]
RMA_CHECK = r"""
suppressPackageStartupMessages({library(jsonlite); library(oligo)})
job <- fromJSON(commandArgs(TRUE)[1])
raw <- read.celfiles(job$files, verbose = FALSE)
summarised <- if (is.null(job$target)) rma(raw) else rma(raw, target = job$target)
writeBin(as.vector(exprs(summarised)), file.path(job$work, "values.bin"), size = 8)
writeLines(rownames(exprs(summarised)), file.path(job$work, "features.txt"))
"""
MEANING = ("PASS means that every included sample has exactly one recorded outcome, a column in one matrix of its study platform or a reason for having "
           "none that a second table confirms; that the platform family of every study platform agrees with the technology that the sample records of "
           "its study name; that every matrix holds the rows, the samples and the value range that its stage recorded, and the matrix "
           "folders hold exactly the files of the stage receipts; that the counts of every stage reconcile with the routes; that every sample whose study "
           "deposits a comparable GEO series matrix is most similar to the column of its own accession there, except where a recorded exception that is "
           "still needed explains why (a matrix read from the series matrix itself is not compared with it); that the streamed RMA agrees with "
           "oligo::rma on the reference studies, a 3' array study and a study summarised at a target; and that the log2 counts per million of every count "
           "matrix follow from its counts, library sizes, TMM factors and normalization groups. It does not assess the biological quality of a sample or "
           "whether a normalization suits an analysis.")


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
        print(f"[09] {name}", flush=True)
        action()


def stage_of(route):
    """Name the stage that processes a route, or None for a study platform that is not read."""
    if route.tier == "raw_array_files":
        return 4 if route.family == "affymetrix" else 5
    return {"unnormalized_intensities": 6, "counts": 7, "author_normalized": 8}.get(route.tier)


def provenance_checks(qc, state):
    state["start"] = C.begin(9)
    qc.check("all_stage_receipts_and_code_bindings", True)


def reason_evidence(studies, samples):
    """The tables by which the reason that a sample has no matrix is confirmed, other than the tables that gave the reason."""
    inventory = pd.read_parquet(C.DATA / "sample_expression.parquet").set_index(["study_accession", "sample_accession"])
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing = listing[(listing["sample_accession"] != "") & (~listing["format"].isin(TABLES) | (listing["value_kind"] != "none"))]
    labels = pd.read_parquet(C.DATA / "label_map.parquet")
    labels = labels[labels["sample_accession"] != ""]
    arrays, collected = samples[5][~samples[5]["in_matrix"]], studies[8][studies[8]["reason"] == "not_collected"]
    return dict(inventory=inventory, files=set(zip(listing["accession"], listing["sample_accession"], listing["format"])),
                named=set(zip(labels["accession"], labels["sample_accession"], labels["table"])),
                incomplete=set(zip(arrays["study_accession"], arrays["sample_accession"])), not_collected=set(zip(collected["accession"], collected["platform"])))


def reason_confirmed(row, evidence):
    """Whether the inventory, the label map or the table of a later script shows the reason that a sample has no matrix."""
    key = (row.study_accession, row.sample_accession)
    if row.reason == "not_processed":
        return row.tier == "not_processed"
    if row.reason == "no_own_data_file":
        formats = C.CONFIG["routes"]["tiers"][row.tier].get("formats", {}).get(row.family, [])
        return not any((*key, name) in evidence["files"] for name in formats)
    if row.reason == "not_in_ncbi_counts":
        return not bool(evidence["inventory"].loc[key, "ncbi_counts"])
    if row.reason == "no_series_matrix_values":
        return int(evidence["inventory"].loc[key, "series_matrix_rows"]) == 0
    if row.reason == "absent_from_table":
        return not any((*key, table) in evidence["named"] for table in row.data_source.split(";"))
    if row.reason == "incomplete_array_file":
        return key in evidence["incomplete"]
    return row.reason == "study_platform_not_collected" and (row.study_accession, row.platform) in evidence["not_collected"]


def technology_disagreements(routes, included):
    """The study platforms whose platform family is not one that measures the technology of their study, as the acquisition stage named it from the
    sample records."""
    assays = C.CONFIG["platforms"]["assays"]
    technology = dict(zip(included["accession"], included["assay"]))
    unknown = sorted({label for value in technology.values() for label in value.split("+")} - set(assays))
    C.require(not unknown, f"Study technologies without platform families: {unknown}; record them in platforms.assays of expression_config.json")
    return [f"{route.accession}|{route.platform}: platform family {route.family}, study technology {technology[route.accession]}"
            for route in routes.itertuples(index=False)
            if route.family not in {family for label in technology[route.accession].split("+") for family in assays[label]}]


def outcome_checks(qc, state):
    routes = pd.read_parquet(C.DATA / "processing_routes.parquet")
    held = pd.read_parquet(C.DATA / "sample_routes.parquet")
    assigned = pd.read_parquet(C.DATA / "table_samples.parquet")
    included = pd.read_parquet(C.ACQUIRED / "included_studies.parquet")
    recorded = pd.read_parquet(C.ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", "platform"])
    recorded = recorded[recorded["study_accession"].isin(included["accession"])]
    pair = lambda frame: list(zip(frame["study_accession"], frame["sample_accession"]))
    qc.check("every_included_sample_routed_once", sorted(pair(held)) == sorted(pair(recorded)) and len(set(pair(held))) == len(held))
    platforms = pd.read_parquet(C.DATA / "study_platforms.parquet")
    qc.check("every_study_platform_has_one_route", sorted(zip(routes["accession"], routes["platform"])) == sorted(zip(platforms["accession"], platforms["platform"]))
             and not routes.duplicated(["accession", "platform"]).any() and routes["samples"].sum() == len(held))
    disagreeing = technology_disagreements(routes, included)
    qc.check("every_platform_family_agrees_with_the_technology_of_its_study", not disagreeing, disagreeing[:50])
    routes["stage"] = [stage_of(route) for route in routes.itertuples(index=False)]
    studies = {number: pd.read_parquet(C.DATA / f"{names[0]}.parquet") for number, names in STAGES.items()}
    samples = {number: pd.read_parquet(C.DATA / f"{names[1]}.parquet") for number, names in STAGES.items()}
    wrong = [number for number, table in studies.items()
             if sorted(zip(table["accession"], table["platform"])) != sorted(zip(*[routes.loc[routes["stage"] == number, column] for column in ("accession", "platform")]))]
    qc.check("every_routed_study_platform_processed_by_its_stage", not wrong, wrong)

    # A sample is expected in a matrix when its data source holds it, its array file is complete and its study platform is collected.
    in_table = assigned.set_index(["study_accession", "sample_accession"])["in_table"]
    outcome = held.merge(routes[["accession", "platform", "stage"]], left_on=["study_accession", "platform"], right_on=["accession", "platform"],
                         validate="many_to_one")
    # For a route that reads a table, the label map of stage 03 says which samples the table holds.
    by_table = ~outcome["data_source"].isin(KEYWORDS)
    outcome.loc[by_table, "in_source"] = [bool(in_table[key]) for key in pair(outcome[by_table])]
    outcome.loc[by_table, "reason"] = ["" if held_there else "absent_from_table" for held_there in outcome.loc[by_table, "in_source"]]
    incomplete = {(key.split("|")[0], sample) for key, entry in C.REGISTER["incomplete_files"].items() for sample in entry["samples"]}
    outcome.loc[[key in incomplete for key in pair(outcome)], "reason"] = "incomplete_array_file"
    skipped = {tuple(key.split("|")) for key in C.REGISTER["not_collected"]}
    outcome.loc[[key in skipped for key in zip(outcome["study_accession"], outcome["platform"])], "reason"] = "study_platform_not_collected"
    expected = outcome["reason"] == ""
    evidence = reason_evidence(studies, samples)
    unconfirmed = sorted(f"{row.study_accession}|{row.sample_accession}: {row.reason}" for row in outcome[~expected].itertuples(index=False)
                         if not reason_confirmed(row, evidence))
    qc.check("every_sample_without_a_matrix_has_a_reason_that_another_table_confirms", not unconfirmed, unconfirmed[:50])
    present = set()
    for number, table in samples.items():
        inside = table[table["in_matrix"]] if "in_matrix" in table.columns else table
        qc.check(f"stage_{number:02d}_samples_listed_once", not table.duplicated(["study_accession", "sample_accession"]).any())
        present |= set(pair(inside))
    missing = sorted(set(pair(outcome[expected])) - present)
    extra = sorted(present - set(pair(outcome[expected])))
    qc.check("every_expected_sample_is_in_a_matrix", not missing, missing[:50])
    qc.check("no_sample_in_a_matrix_without_a_data_source", not extra, extra[:50])
    outcome["in_matrix"] = expected
    by_tier = [dict(tier=tier, study_platforms=len(part[["study_accession", "platform"]].drop_duplicates()), samples=len(part),
                    in_matrix=int(part["in_matrix"].sum())) for tier, part in outcome.groupby("tier", sort=False)]
    reasons = {str(reason): int(count) for reason, count in outcome.loc[~expected, "reason"].value_counts().sort_index().items()}
    # The merge with the routes keeps only the samples whose study platform has a route.
    qc.check("every_routed_sample_has_an_outcome", len(outcome) == len(held) == len(recorded))
    state.update(routes=routes, studies=studies, samples=samples,
                 summary=dict(samples=len(outcome), in_a_matrix=int(outcome["in_matrix"].sum()), without_a_matrix=reasons, by_tier=by_tier))


def within(found, expected, relative, absolute):
    """Whether every value agrees with its expected value within the tolerances; a value that is not finite never agrees."""
    found, expected = np.asarray(found, dtype="float64"), np.asarray(expected, dtype="float64")
    return found.shape == expected.shape and bool((np.abs(found - expected) <= absolute + relative * np.abs(expected)).all())


def matrix_checks(qc, state):
    shape, values, counts, tally = [], [], [], dict(matrices=0, count_matrices=0, values=0)
    tolerance = C.CONFIG["quality_control"]["matrix_check"]
    for number, (_, _, rows) in STAGES.items():
        table, listed = state["studies"][number], state["samples"][number]
        listed = listed[listed["in_matrix"]] if "in_matrix" in listed.columns else listed
        for study in table[table["matrix"].fillna("") != ""].itertuples(index=False):
            name = f"{study.accession}|{study.platform}"
            frame = pd.read_parquet(C.project_path(study.matrix))
            own = listed[(listed["study_accession"] == study.accession) & (listed["platform"] == study.platform)]
            block = frame.iloc[:, 1:].to_numpy()
            if (frame.columns[0] != "feature_id" or not frame["feature_id"].is_unique or (frame["feature_id"] == "").any() or len(frame) != getattr(study, rows)
                    or list(frame.columns[1:]) != own["sample_accession"].tolist() or block.dtype.kind != "f"):
                shape.append(name)
            absent = int(np.isnan(block).sum())
            allowed = int(getattr(study, "missing_values", 0))
            low, high = (getattr(study, field, None) for field in ("value_minimum", "value_maximum"))
            if number == 7:
                # The table of count studies does not carry the range of the log2 counts per million; the record beside the matrix does.
                stored = C.read_json(C.project_path(study.matrix).with_suffix(".json"))
                low, high = stored["value_minimum"], stored["value_maximum"]
            if (absent != allowed or np.isinf(block).any() or low is None or high is None
                    or not within([np.nanmin(block), np.nanmax(block)], [low, high], tolerance["value_range_relative"], tolerance["value_range_absolute"])):
                values.append(name)
            tally["matrices"] += 1
            tally["values"] += int(block.size)
            if number == 7:
                # A count study has a second matrix, the counts themselves, on the same genes and samples; their sums are the recorded library sizes.
                raw = pd.read_parquet(C.project_path(study.counts))
                totals = raw.iloc[:, 1:].to_numpy(dtype="float64")
                if (list(raw.columns) != list(frame.columns) or raw["feature_id"].tolist() != frame["feature_id"].tolist() or not np.isfinite(totals).all()
                        or (totals < 0).any() or not within(totals.sum(axis=0), own["library_size"].to_numpy(), tolerance["library_size_relative"], 0.0)):
                    counts.append(name)
                tally["count_matrices"] += 1
            if tally["matrices"] % 20 == 0:
                print(f"[09] verified {tally['matrices']:,} matrices, {tally['values']:,} values", flush=True)
    print(f"[09] verified {tally['matrices']:,} matrices, {tally['values']:,} values", flush=True)
    qc.check("matrix_rows_and_samples_agree_with_the_stage_tables", not shape, shape)
    qc.check("matrix_values_agree_with_their_records", not values, values)
    qc.check("count_matrices_agree_with_the_library_sizes", not counts, counts)
    state["summary"].update(tally)


def group_checks(qc, state):
    """Every study platform with recorded normalization groups is processed in them, and every count matrix follows from its counts and factors."""
    decisions = C.REGISTER["normalization_groups"]
    grouped = pd.concat([state["samples"][4], state["samples"][7]], ignore_index=True)
    wrong = []
    for (study, platform), part in grouped.groupby(["study_accession", "platform"]):
        key = f"{study}|{platform}"
        expected = C.normalization_groups(key, part["sample_accession"].tolist()) if key in decisions else None
        found = part["normalization_group"].tolist()
        if (expected is None and set(found) != {""}) or (expected is not None and found != [expected[sample] for sample in part["sample_accession"]]):
            wrong.append(key)
    processed = set(zip(grouped["study_accession"], grouped["platform"]))
    qc.check("every_recorded_normalization_group_is_a_processed_study_platform", all(tuple(key.split("|")) in processed for key in decisions),
             sorted(key for key in decisions if tuple(key.split("|")) not in processed))
    qc.check("samples_are_normalized_in_their_recorded_groups", not wrong, wrong)
    prior, samples, largest, differing = C.CONFIG["count_tables"]["prior_count"], state["samples"][7], 0.0, []
    tolerance = C.CONFIG["quality_control"]["matrix_check"]["log_cpm_absolute"]
    for study in state["studies"][7].itertuples(index=False):
        own = samples[(samples["study_accession"] == study.accession) & (samples["platform"] == study.platform)]
        counts = pd.read_parquet(C.project_path(study.counts)).iloc[:, 1:].to_numpy(dtype="float64")
        values = pd.read_parquet(C.project_path(study.matrix)).iloc[:, 1:].to_numpy(dtype="float64")
        effective = own["library_size"].to_numpy() * own["norm_factor"].to_numpy()
        labels = own["normalization_group"].to_numpy()
        # edgeR scales the prior count to the mean TMM-scaled library size of the samples that were normalized together.
        scaled = prior * effective / np.array([effective[labels == label].mean() for label in labels])
        expected = np.log2((counts + scaled) / (effective + 2 * scaled) * 1e6)
        floor = np.log2(scaled / (effective + 2 * scaled) * 1e6)
        gaps = np.concatenate([np.abs(expected - values).ravel(), np.abs(floor - own["log2_cpm_of_zero"].to_numpy())])
        if not (gaps <= tolerance).all():
            differing.append(f"{study.accession}|{study.platform}")
        finite = gaps[np.isfinite(gaps)]
        largest = max(largest, float(finite.max()) if finite.size else 0.0)
    qc.check("count_log_cpm_follows_from_counts_factors_and_groups", not differing, dict(largest_finite_difference=round(largest, 8), study_platforms=differing))
    extreme = samples[np.abs(np.log2(samples["norm_factor"])) > C.CONFIG["quality_control"]["tmm_report_log2"]]
    state["summary"].update(
        normalization_groups={key: {str(name): int(count) for name, count in grouped.loc[(grouped["study_accession"] + "|" + grouped["platform"]) == key,
                                                                                         "normalization_group"].value_counts().sort_index().items()}
                              for key in decisions},
        largest_log_cpm_difference=round(largest, 8),
        tmm_factors_beyond_report_limit={f"{study}|{platform}": int(count) for (study, platform), count in extreme.groupby(["study_accession", "platform"]).size().items()},
        log2_cpm_of_zero=dict(lowest=round(float(samples["log2_cpm_of_zero"].min()), 4), highest=round(float(samples["log2_cpm_of_zero"].max()), 4)))


def opened_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(encoding="utf-8", errors="replace")


def series_columns(path):
    """The column labels of the value table of a GEO series matrix, and the number of lines before them."""
    with opened_text(path) as text:
        for number, line in enumerate(text, 1):
            if line.startswith("!series_matrix_table_begin"):
                return [label.strip().strip('"') for label in next(text).rstrip("\n").split("\t")], number
    return [], 0


def series_values(path, skip, columns):
    """Read chosen columns of the value table of a GEO series matrix as numbers, indexed by its row identifiers; repeated identifiers are left out."""
    with opened_text(path) as text:
        frame = pd.read_csv(text, sep="\t", skiprows=skip, usecols=columns, index_col=0, dtype={columns[0]: str}, comment="!", low_memory=False)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return frame[~frame.index.duplicated(keep=False)]


def centred_ranks(values):
    """Rank the rows within each column, remove the mean rank of each row, and scale every column to unit variance."""
    ranks = pd.DataFrame(values).rank().to_numpy(dtype="float64")
    ranks -= ranks.mean(axis=1, keepdims=True)
    ranks -= ranks.mean(axis=0)
    spread = ranks.std(axis=0)
    return ranks / np.where(spread > 0, spread, 1.0)


def value_checks(qc, state):
    """Compare every matrix with the GEO series matrix of its study, sample by sample, where both hold the same row identifiers."""
    policy, exceptions = C.CONFIG["quality_control"]["value_check"], C.REGISTER["value_check_exceptions"]
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing = listing[(listing["format"] == "geo_series_matrix") & (listing["value_kind"] != "none") & (listing["sample_accession"] == "")]
    rows, failed, failing = [], [], {}
    for number in STAGES:
        table, listed = state["studies"][number], state["samples"][number]
        listed = listed[listed["in_matrix"]] if "in_matrix" in listed.columns else listed
        for study in table[table["matrix"].fillna("") != ""].itertuples(index=False):
            key = f"{study.accession}|{study.platform}"
            mine = listed[(listed["study_accession"] == study.accession) & (listed["platform"] == study.platform)]
            own = mine["sample_accession"].tolist()
            row = dict(source=study.source, accession=study.accession, platform=study.platform, series_matrix="", status="no_series_matrix_values",
                       samples_compared=0, own_column_most_similar=0, shared_rows=0, rows_used=0, median_own_correlation=np.nan, lowest_own_correlation=np.nan,
                       smallest_margin=np.nan, exception=key in exceptions, not_most_similar="")
            best = None
            if number == 8 and study.data_source == "series_matrix":
                # A matrix that stage 08 read from the series matrix itself would only be compared with itself.
                row.update(series_matrix=mine["table"].iloc[0], status="same_source")
            else:
                for item in listing[listing["accession"] == study.accession].itertuples(index=False):
                    path = C.RAW / item.source / item.accession / item.file
                    labels, skip = series_columns(path)
                    shared = [sample for sample in own if sample in labels]
                    if labels and (best is None or len(shared) > len(best[2])):
                        best = (path, skip, shared, labels[0])
            if best is not None and len(best[2]) < policy["minimum_samples"]:
                row.update(series_matrix=best[0].name, status="too_few_shared_samples")
            elif best is not None:
                path, skip, shared, identifier = best
                matrix = pd.read_parquet(C.project_path(study.matrix), columns=["feature_id", *shared]).set_index("feature_id")
                deposit = series_values(path, skip, [identifier, *shared])
                matrix.index, deposit.index = matrix.index.astype(str), deposit.index.astype(str)
                common = sorted(set(matrix.index) & set(deposit.index))
                complete = [name for name, ours, theirs in zip(common, np.isfinite(matrix.loc[common].to_numpy(dtype="float64")).all(axis=1),
                                                                np.isfinite(deposit.loc[common].to_numpy(dtype="float64")).all(axis=1)) if ours and theirs]
                row.update(series_matrix=path.name, shared_rows=len(complete), status="too_few_shared_rows")
                if len(complete) >= policy["minimum_shared_rows"]:
                    # Rows spread evenly over the sorted identifiers, so that the choice depends on nothing but the identifiers.
                    chosen = [complete[position] for position in np.unique(np.linspace(0, len(complete) - 1, min(len(complete), policy["rows_used"])).astype(int))]
                    scores = centred_ranks(matrix.loc[chosen, shared].to_numpy(dtype="float64")).T @ centred_ranks(deposit.loc[chosen, shared].to_numpy(dtype="float64"))
                    scores /= len(chosen)
                    own_score = np.diag(scores)
                    others = np.where(np.eye(len(shared), dtype=bool), -np.inf, scores).max(axis=1)
                    hits = scores.argmax(axis=1) == np.arange(len(shared))
                    row.update(status="compared", samples_compared=len(shared), own_column_most_similar=int(hits.sum()), rows_used=len(chosen),
                               median_own_correlation=round(float(np.median(own_score)), 6), lowest_own_correlation=round(float(own_score.min()), 6),
                               smallest_margin=round(float((own_score - others).min()), 6),
                               not_most_similar="; ".join(sample for sample, hit in zip(shared, hits) if not hit)[:500])
                    # An exception covers a whole study platform, or only the samples that it names.
                    failing[key] = {sample for sample, hit in zip(shared, hits) if not hit}
                    allowed = set(exceptions[key].get("samples", shared)) if key in exceptions else set()
                    if failing[key] - allowed:
                        failed.append(key)
            rows.append(row)
            print(f"[09] values of {key}: {row['status']}, {row['own_column_most_similar']} of {row['samples_compared']} most similar to their own column",
                  flush=True)
    review = pd.DataFrame(rows, columns=VALUE_COLUMNS)
    C.save_review(review, C.OUT / "value_check_review.csv")
    compared = review[review["status"] == "compared"]
    qc.check("every_compared_sample_is_most_similar_to_its_own_series_matrix_column", not failed, failed)
    qc.check("every_value_check_exception_names_a_compared_study_platform",
             set(exceptions) <= set(compared["accession"] + "|" + compared["platform"]), sorted(exceptions))
    # An exception that no longer applies must be removed: one of a whole study platform where every sample is most similar to its own column, or one that
    # names a sample that is now most similar to its own column.
    needless = sorted(key for key, entry in exceptions.items()
                      if not failing.get(key) or ("samples" in entry and set(entry["samples"]) - failing[key]))
    qc.check("every_value_check_exception_is_still_needed", not needless, needless)
    state["summary"].update(value_check=dict(study_platforms_compared=len(compared), samples_compared=int(compared["samples_compared"].sum()),
                                             own_column_most_similar=int(compared["own_column_most_similar"].sum()), exceptions=sorted(exceptions),
                                             not_compared={str(name): int(count) for name, count in review.loc[review["status"] != "compared", "status"].value_counts().items()}))


def oligo_rma(name, members, target):
    """Summarise the CEL files of one study platform with oligo::rma itself, at the target of its design where it has one."""
    A = importlib.import_module("04_process_affymetrix_arrays")
    work = C.WORK / "rma_check"
    shutil.rmtree(work, ignore_errors=True)
    try:
        files = A.stage_arrays(members, work / "cel")
        C.write_json(work / "job.json", dict(files=[path.as_posix() for path in files], work=work.as_posix(), target=target))
        (work / "rma_check.R").write_text(RMA_CHECK, encoding="utf-8")
        done = subprocess.run(["Rscript", "--vanilla", (work / "rma_check.R").as_posix(), (work / "job.json").as_posix()], capture_output=True, text=True)
        C.require(done.returncode == 0, f"oligo::rma failed on the reference study {name}:\n{done.stderr[-3000:]}")
        features = (work / "features.txt").read_text(encoding="utf-8").splitlines()
        return features, np.fromfile(work / "values.bin", dtype="<f8").reshape(len(members), len(features)).T
    finally:
        C.clear_work(work)


def rma_check(qc, state):
    """Summarise every reference Affymetrix study again with oligo::rma itself and compare it with the stored matrix: a 3' array and an array summarised
    at a target, so that both ways in which stage 04 takes the probes from the design package are replayed."""
    policy, designs = C.CONFIG["quality_control"]["rma_reference"], C.CONFIG["affymetrix"]["designs"]
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    listing = listing[(listing["format"] == "affymetrix_cel") & (listing["sample_accession"] != "")].copy()
    # A file that the repository links to the wrong record is read as the sample that the register names, as stage 04 reads it.
    for key, entry in C.REGISTER["reassigned_files"].items():
        moved = (listing["accession"] == key.split("|")[0]) & listing["member"].isin(entry["files"])
        listing.loc[moved, "sample_accession"] = listing.loc[moved, "member"].map(entry["files"])
    studies = state["studies"][4].set_index(["accession", "platform"])
    targets = {name: designs[studies.loc[tuple(name.split("|")), "chip_type"]].get("target") for name in policy["studies"]}
    qc.check("rma_references_replay_both_ways_of_taking_the_probes", {target is not None for target in targets.values()} == {False, True}, targets)
    results = []
    for name, target in targets.items():
        accession, platform = name.split("|")
        study = studies.loc[(accession, platform)]
        arrays = state["samples"][4]
        arrays = arrays[(arrays["study_accession"] == accession) & (arrays["platform"] == platform)]
        C.require(set(arrays["normalization_group"]) == {""}, f"The RMA reference study {name} must be one normalization group")
        own = listing[(listing["accession"] == accession) & listing["sample_accession"].isin(arrays["sample_accession"])]
        members = list(own.sort_values("sample_accession").itertuples(index=False))
        C.require([item.sample_accession for item in members] == arrays["sample_accession"].tolist(), f"The RMA reference study {name} does not hold one CEL file per array")
        features, reference = oligo_rma(name, members, target)
        stored = pd.read_parquet(C.project_path(study["matrix"])).set_index("feature_id")
        same_rows = len(features) == len(set(features)) and sorted(features) == sorted(stored.index)
        gaps = np.abs(stored.loc[features, arrays["sample_accession"]].to_numpy(dtype="float64") - reference) if same_rows else np.array([np.inf])
        results.append(dict(study=name, chip_type=study["chip_type"], target=target or "", arrays=len(members), probe_sets=len(features),
                            largest_difference=float(gaps.max()) if np.isfinite(gaps).all() else None,
                            agrees=bool(same_rows and (gaps <= policy["tolerance"]).all())))
    qc.check("streamed_rma_agrees_with_oligo_rma", all(result["agrees"] for result in results), results)
    state["summary"].update(rma_check=results)


def total_checks(qc, state):
    receipts = {number: C.read_json(C.receipt_path(number)) for number in range(2, 9)}
    studies, samples, routes = state["studies"], state["samples"], state["routes"]
    inside = {number: int(table["in_matrix"].sum()) if "in_matrix" in table.columns else len(table) for number, table in samples.items()}
    qc.check("route_receipt_reconciles", receipts[2]["samples"]["total"] == state["summary"]["samples"]
             and {row["tier"]: row["study_platforms"] for row in receipts[2]["by_tier"]} == routes.groupby("tier").size().to_dict())
    qc.check("affymetrix_receipt_reconciles", receipts[4]["arrays"] == inside[4] and receipts[4]["study_platforms"]["processed"] == len(studies[4]))
    qc.check("single_channel_receipt_reconciles", receipts[5]["arrays"]["in_matrix"] == inside[5] and receipts[5]["study_platforms"]["processed"] == len(studies[5]))
    qc.check("intensity_table_receipt_reconciles", receipts[6]["arrays"]["in_matrix"] == inside[6] and receipts[6]["study_platforms"]["processed"] == len(studies[6]))
    qc.check("count_table_receipt_reconciles", receipts[7]["samples"]["in_matrix"] == inside[7] and receipts[7]["study_platforms"]["processed"] == len(studies[7]))
    qc.check("author_matrix_receipt_reconciles", receipts[8]["samples"]["in_matrix"] == inside[8] and receipts[8]["study_platforms"]["processed"] == len(studies[8]))
    qc.check("samples_in_matrices_reconcile", sum(inside.values()) == state["summary"]["in_a_matrix"])
    with_matrix = sum(int((table["matrix"].fillna("") != "").sum()) for table in studies.values())
    qc.check("one_matrix_for_every_processed_study_platform", with_matrix == state["summary"]["matrices"]
             and with_matrix == len(routes) - int(routes["stage"].isna().sum()) - len(C.REGISTER["not_collected"]))
    # The matrix folders hold the files of the current routes only; a file of a removed or renamed route would be released with them.
    folders = (C.MATRICES.resolve(), C.COUNTS.resolve())
    bound = {name for number in range(4, 9) for name in receipts[number]["outputs"] if C.project_path(name).parent in folders}
    held = {path.resolve().relative_to(C.ROOT).as_posix() for folder in folders if folder.is_dir() for path in folder.iterdir()}
    qc.check("matrix_folders_hold_exactly_the_files_of_the_stage_receipts", held == bound,
             dict(not_bound=sorted(held - bound)[:50], missing=sorted(bound - held)[:50]))
    state["summary"].update(study_platforms=len(routes), study_platforms_with_a_matrix=with_matrix,
                            samples_by_stage={C.STAGES[number]: count for number, count in inside.items()})


def main():
    qc, state = QC(), {}
    failure = None
    try:
        qc.phase("1/7 Receipts, code bindings and stage outputs", lambda: provenance_checks(qc, state))
        qc.phase("2/7 Routes, study platforms and sample outcomes", lambda: outcome_checks(qc, state))
        qc.phase("3/7 Matrices against the stage tables", lambda: matrix_checks(qc, state))
        qc.phase("4/7 Normalization groups and log2 counts per million", lambda: group_checks(qc, state))
        qc.phase("5/7 Matrices against the series matrices of their studies", lambda: value_checks(qc, state))
        qc.phase("6/7 Streamed RMA against oligo::rma", lambda: rma_check(qc, state))
        qc.phase("7/7 Stage receipts and totals", lambda: total_checks(qc, state))
    except Exception:
        failure = traceback.format_exc()
        if not qc.checks or qc.checks[-1]["passed"]:
            qc.checks.append(dict(section=qc.section, name="phase_completed", passed=False, detail=failure))
    passed = failure is None and bool(qc.checks)
    report = dict(version=C.VERSION, status="PASS" if passed else "FAIL", generated_at_utc=datetime.now(timezone.utc).isoformat(),
                  required_checks=len(qc.checks), failures=sum(not row["passed"] for row in qc.checks), checks=qc.checks,
                  runtime=C.runtime(), code=C.code_binding(), exception=failure,
                  dependencies=C.bind([C.receipt_path(number) for number in range(1, 9) if C.receipt_path(number).exists()]),
                  eligibility="EXPRESSION_VERIFIED" if passed else "NOT_ELIGIBLE", interpretation=MEANING, summary=state.get("summary"))
    C.write_json(C.OUT / "qc.json", report)
    print(f"\n[09 {'PASS' if passed else 'FAIL'}] {report['eligibility']}", flush=True)
    print(f"    Checks: {report['required_checks']}; failures: {report['failures']}; report: {(C.OUT / 'qc.json').relative_to(C.ROOT).as_posix()}", flush=True)
    if passed:
        summary = report["summary"]
        print(f"    {summary['in_a_matrix']:,} of {summary['samples']:,} samples are in one of {summary['matrices']:,} matrices; without a matrix: "
              f"{summary['without_a_matrix']}", flush=True)
    if failure:
        print(failure, flush=True)
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
