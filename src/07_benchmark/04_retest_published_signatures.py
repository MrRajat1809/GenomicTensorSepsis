"""Score published sepsis signatures by their published formulas within every study platform and tissue, and summarise their AUROC across studies."""
import importlib

C = importlib.import_module("01_audit_split_leakage")

REASONS = {"definition_not_verified": "The definition of the signature could not be verified against its publication.",
           "scale_not_comparable": "The gene matrix of the study platform is on a scale that the validation stage does not take as log2.",
           "fewer_than_minimum_of_a_class": "A class of the comparison has fewer than signatures.minimum_group samples in the stratum.",
           "gene_coverage_below_threshold": "Fewer than the signature's minimum_gene_coverage of its genes are measured on the platform.",
           "too_few_genes_on_a_side": "A side of the signature that names genes has fewer than its minimum_genes_per_side measured on the platform.",
           "fewer_than_minimum_with_a_score": "A class has fewer than signatures.minimum_group samples with a score (every gene of the score missing)."}
GENE_COLUMNS = ["signature", "side", "published_symbol", "hgnc_id", "symbol", "resolved_by"]
COVERAGE_COLUMNS = ["signature", "source", "study_accession", "platform", "value_scale", "comparable_scale", "genes_published", "up_published", "up_present",
                    "down_published", "down_present", "coverage", "missing"]
TEST_COLUMNS = ["signature", "label", "type", "variant", "primary_variant", "comparison", "primary_comparison", "source", "study_accession", "platform", "tissue",
                "value_scale", "positives", "negatives", "scored_positives", "scored_negatives", "patients_positive", "patients_negative", "genes_up_used",
                "genes_down_used", "coverage", "down_weight", "score_shift", "status", "reason", "auroc", "auroc_low", "auroc_high", "derivation_role",
                "samples_shared_with_derivation", "repeat_deposits", "in_scope_group", "independent", "doi", "pmid"]
SCORE_COLUMNS = ["signature", "variant", "comparison", *C.KEYS, "tissue", "label", "patient", "score"]
SUMMARY_COLUMNS = ["signature", "label", "type", "variant", "primary_variant", "comparison", "primary_comparison", "scope", "strata_set", "strata", "tested",
                   *(f"skipped_{reason}" for reason in REASONS), "auroc_median", "auroc_q1", "auroc_q3", "share_above_limit", "citation", "doi", "pmid"]
STRATUM = ["study_accession", "platform", "tissue"]
WHOLE_SCOPE = "all_tissues"
STRATA_SETS = ("all", "independent")


def log_values(frame, scale):
    """Return the values on a log2 scale as the validation stage does: as they are on a log scale, log2 of the value plus an offset on a linear one."""
    import numpy as np

    policy = C.CONFIG["values"]
    if scale in policy["log_scales"]:
        return frame
    C.require(scale in policy["linear_scales"], f"Values of this scale cannot be compared: {scale}")
    return np.log2(frame.clip(lower=0) + policy["linear_offset"])


def read_genes(path, identifiers):
    """Read the rows of the named genes from one gene matrix: genes in rows, by HGNC identifier, and samples in columns."""
    import pandas as pd

    frame = pd.read_parquet(path, filters=[("hgnc_id", "in", sorted(identifiers))])
    return frame.drop(columns="symbol").set_index("hgnc_id").astype("float64")


def matrix_minimum(path, scale):
    """The lowest value of a gene matrix on the log2 scale of the scores."""
    import pandas as pd

    frame = pd.read_parquet(path).drop(columns=["hgnc_id", "symbol"]).astype("float64")
    return float(log_values(frame, scale).min().min())


def geometric_mean(values):
    """The geometric mean of every column over the rows that hold a value, as MetaIntegrator and the validation stage compute it."""
    import numpy as np

    return np.exp(np.log(values).mean(axis=0, skipna=True))


def signature_score(values, up, down, formula, weight, shift):
    """The score of every sample (column) by a published formula: geometric means of the raised values of each side, the lower side times its weight,
    or sums of the log values of each side; a side without genes counts as zero."""
    if formula == "geometric_mean_difference":
        high = geometric_mean(values.loc[up] + shift) if up else 0.0
        low = geometric_mean(values.loc[down] + shift) if down else 0.0
        return high - weight * low
    return values.loc[up].sum(axis=0, skipna=False) - values.loc[down].sum(axis=0, skipna=False)


def down_weight(rule, up, down):
    """The weight of the lower side: one, or the ratio of the counts of the genes present on the two sides."""
    if rule == "one":
        return 1.0
    C.require(up and down, "A weight from the counts of genes needs genes on both sides")
    return len(down) / len(up) if rule == "down_over_up" else len(up) / len(down)


def auroc(cases, others):
    """Area under the ROC curve of a score: the share of case and other pairs in which the case has the higher score, ties counted as half."""
    import numpy as np
    import pandas as pd

    cases, others = np.asarray(cases, dtype=float), np.asarray(others, dtype=float)
    ranks = pd.Series(np.concatenate([cases, others])).rank().to_numpy()
    return float((ranks[:len(cases)].sum() - len(cases) * (len(cases) + 1) / 2) / (len(cases) * len(others)))


def cluster_interval(cases, others, case_groups, other_groups, policy, seed):
    """A percentile interval of the area under the curve from a bootstrap that resamples patients, with all of their samples, within each group, as the
    validation stage draws it."""
    import numpy as np

    generator = np.random.default_rng(seed)
    case_sets = [np.flatnonzero(case_groups == name) for name in np.unique(case_groups)]
    other_sets = [np.flatnonzero(other_groups == name) for name in np.unique(other_groups)]
    found = []
    for _ in range(policy["replicates"]):
        first = np.concatenate([case_sets[index] for index in generator.integers(0, len(case_sets), len(case_sets))])
        second = np.concatenate([other_sets[index] for index in generator.integers(0, len(other_sets), len(other_sets))])
        found.append(auroc(cases[first], others[second]))
    tail = (1 - policy["level"]) / 2
    return float(np.quantile(found, tail)), float(np.quantile(found, 1 - tail))


def resolve(definition, reference):
    """The HGNC gene of every gene of a signature, from its symbol map, checked against the HGNC table: the identifier carries the approved symbol, and the
    published symbol is that symbol or the previous symbol of that one gene alone."""
    approved = dict(zip(reference["symbol"], reference["hgnc_id"]))
    current = dict(zip(reference["hgnc_id"], reference["symbol"]))
    earlier = {}
    for identifier, names in zip(reference["hgnc_id"], reference["prev_symbol"].fillna("")):
        for name in filter(None, names.split("|")):
            earlier.setdefault(name, []).append(identifier)
    rows = []
    for row in definition["symbol_map"]:
        published, symbol, identifier = row["published"], row["hgnc_symbol"], row["hgnc_id"]
        C.require(current.get(identifier) == symbol, f"{identifier} is not the HGNC gene {symbol} (signature {definition['name']})")
        by = "approved symbol" if approved.get(published) == identifier else "previous symbol" if earlier.get(published) == [identifier] else ""
        C.require(by, f"The published symbol {published} of signature {definition['name']} does not name {identifier} alone")
        rows.append(dict(signature=definition["name"], side="up" if symbol in definition["genes_up"] else "down", published_symbol=published, hgnc_id=identifier,
                         symbol=symbol, resolved_by=by))
    return rows


def comparison_members(samples, rule):
    """The positive and negative samples of a comparison; a sepsis sample counts only in the phases of signatures.sepsis_phases."""
    value = samples[rule["variable"]]
    allowed = (samples["sample_class"] != "sepsis") | samples["sepsis_phase"].isin(C.CONFIG["signatures"]["sepsis_phases"])
    within = samples["sample_class"] == rule["within_class"] if "within_class" in rule else True
    return (value == rule["positive"]) & within & allowed, (value == rule["negative"]) & within & allowed


def variants(definition):
    return [("primary", definition["down_weight"]), *((variant["variant"], variant["down_weight"]) for variant in definition["sensitivity"])]


def stage_05_checks(found, tests, scores):
    """Require that the Sepsis MetaScore of this stage gives the validation stage's AUROC, interval and sample scores in every stratum it tested."""
    import numpy as np
    import pandas as pd

    policy = C.CONFIG["signatures"]
    name = policy["stage_05_score"]
    definition = next(entry for entry in policy["definitions"] if entry["name"] == name)
    C.require("sepsis_vs_healthy_control" in definition["intended_comparison"] and definition["down_weight"] == "one"
              and definition["formula"] == "geometric_mean_difference", f"Signature {name} must be scored as the validation stage scores it, against healthy controls")
    named = pd.read_parquet(C.VALIDATION / "validation_genes.parquet")
    named = named[named["test"] == "sepsis_score"]
    sides, found = {"higher_in_sepsis": "up", "lower_in_sepsis": "down"}, found[found["signature"] == name]
    C.require(sorted(zip(named["published_symbol"], named["hgnc_id"], named["role"].map(sides))) == sorted(
        zip(found["published_symbol"], found["hgnc_id"], found["side"])), f"Signature {name} must name the genes of the validation stage's score")
    stage = pd.read_parquet(C.VALIDATION / "sepsis_score_by_study.parquet")
    stage = stage[stage["score_computed"]].set_index(STRATUM)
    mine = tests[(tests["signature"] == name) & (tests["variant"] == "primary") & (tests["comparison"] == "sepsis_vs_healthy_control")
                 & (tests["status"] == "tested")].set_index(STRATUM)
    C.require(set(mine.index) == set(stage.index), f"Signature {name} must be tested in exactly the strata of the validation stage")
    mine = mine.loc[stage.index]
    gap = np.abs(mine["auroc"].to_numpy() - stage["auroc"].to_numpy())
    low, high = (np.where(mine[side].isna() & stage[side].isna(), 0.0, np.abs(mine[side].to_numpy() - stage[side].to_numpy()))
                 for side in ("auroc_low", "auroc_high"))
    tolerance = C.CONFIG["quality_control"]["metric_tolerance"]
    C.require((gap <= tolerance).all() and (np.nan_to_num(low, nan=np.inf) <= tolerance).all() and (np.nan_to_num(high, nan=np.inf) <= tolerance).all(),
              f"Signature {name} disagrees with the AUROC or interval of the validation stage")
    theirs = pd.read_parquet(C.VALIDATION / "sepsis_score_samples.parquet").set_index(["study_accession", "sample_accession", "platform", "tissue"])["score"]
    own = scores[(scores["signature"] == name) & (scores["variant"] == "primary") & (scores["comparison"] == "sepsis_vs_healthy_control")]
    own = own.set_index(["study_accession", "sample_accession", "platform", "tissue"])["score"]
    C.require(set(own.index) == set(theirs.index), f"Signature {name} must score exactly the samples of the validation stage")
    distance = np.abs(own.loc[theirs.index].to_numpy() - theirs.to_numpy())
    C.require((distance <= tolerance).all(), f"The sample scores of signature {name} disagree with those of the validation stage")
    # The validation stage's rule of independence: no role in the derivation of the score, no sample shared with a study that had one, no repeat deposit.
    alone = stage[(stage["score_derivation"] == "") & (stage["samples_shared_with_derivation"] == 0) & (stage["repeat_deposits"] == 0)]
    apart = mine[mine["independent"]]
    C.require(set(apart.index) == set(alone.index) and (apart.loc[alone.index, "repeat_deposits"].to_numpy() == alone["repeat_deposits"].to_numpy()).all()
              and (mine.loc[stage.index, "samples_shared_with_derivation"].to_numpy() == stage["samples_shared_with_derivation"].to_numpy()).all()
              and (mine.loc[stage.index, "derivation_role"].to_numpy() == stage["score_derivation"].to_numpy()).all(),
              f"The independent strata of signature {name} must be those of the validation stage")
    return dict(strata=len(mine), samples=len(own), largest_auroc_difference=float(gap.max()), largest_interval_difference=float(np.nanmax([*low, *high])),
                largest_score_difference=float(distance.max()), independent_strata=len(apart), independent_auroc_median=float(apart["auroc"].median()))


def summarise(tests, definitions, tissues):
    """The AUROC of every signature, variant and comparison across the strata: over every tissue and over the tissues of the scope group, each over all
    strata and over those independent of the signature's derivation."""
    import numpy as np

    policy = C.CONFIG["signatures"]
    rows, limit = [], policy["auroc_report_limit"]
    for (name, variant, comparison), part in tests.groupby(["signature", "variant", "comparison"], sort=False):
        definition = definitions[name]
        for scope, within in ((WHOLE_SCOPE, part), (policy["scope_group"], part[part["tissue"].isin(tissues)])):
            for kept, chosen in zip(STRATA_SETS, (within, within[within["independent"]])):
                tested = chosen.loc[chosen["status"] == "tested", "auroc"].to_numpy(dtype="float64")
                rows.append(dict(signature=name, label=definition["label"], type=definition["type"], variant=variant, primary_variant=variant == "primary",
                                 comparison=comparison, primary_comparison=definition["intended_comparison"][0] == comparison, scope=scope, strata_set=kept,
                                 strata=len(chosen), tested=len(tested), **{f"skipped_{reason}": int((chosen["reason"] == reason).sum()) for reason in REASONS},
                                 auroc_median=float(np.median(tested)) if len(tested) else np.nan,
                                 auroc_q1=float(np.quantile(tested, 0.25)) if len(tested) else np.nan,
                                 auroc_q3=float(np.quantile(tested, 0.75)) if len(tested) else np.nan,
                                 share_above_limit=float((tested > limit).mean()) if len(tested) else np.nan, citation=definition["citation"],
                                 doi=definition["doi"], pmid=definition["pmid"]))
    return rows


def scope_tissues():
    """The tissues of the scope group, as the correction stage's configuration, bound by its quality control, defines them."""
    qc = C.read_json(C.EARLIER["Correction"] / "qc.json")
    C.verify_files(qc["code"])
    groups = C.read_json(C.project_path("src/06_batch_correction/correction_config.json"))["groups"]
    name = C.CONFIG["signatures"]["scope_group"]
    C.require(name in groups, f"signatures.scope_group must name a group of the correction stage: {sorted(groups)}")
    return list(groups[name]["tissues"])


def main():
    import numpy as np
    import pandas as pd

    print(f"[04] Retesting published signatures v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["signatures"]
    stage = C.read_json(C.EARLIER["Validation"] / "04.json")
    C.require(stage["limits"]["minimum_group"] == policy["minimum_group"] and stage["bootstrap"] == policy["bootstrap"]
              and stage["sepsis_phases"] == policy["sepsis_phases"], "signatures.minimum_group, bootstrap and sepsis_phases must equal those of the validation stage")
    qc = C.read_json(C.EARLIER["Validation"] / "qc.json")
    C.verify_files(qc["code"])
    C.require(C.read_json(C.project_path("src/05_technical_validation/validation_config.json"))["values"] == C.CONFIG["values"],
              "The value policy must equal that of the validation stage")
    samples = pd.read_parquet(C.COHORT / "cohort_samples.parquet")
    matrices = pd.read_parquet(C.VALIDATION / "validation_matrices.parquet")
    samples = samples.merge(matrices[["study_accession", "platform", "gene_matrix", "value_scale", "comparable_scale"]], on=["study_accession", "platform"],
                            how="left", validate="many_to_one")
    C.require(samples["gene_matrix"].notna().all(), "Every record of the cohort table needs its gene matrix")
    samples["patient"] = samples["patient_key"].where(samples["patient_key"] != "", samples["sample_id"])
    reference = pd.read_parquet(C.ANNOTATION / "genes.parquet", columns=["hgnc_id", "symbol", "prev_symbol"])
    definitions = {entry["name"]: entry for entry in policy["definitions"]}
    found = pd.DataFrame([row for entry in policy["definitions"] for row in resolve(entry, reference)], columns=GENE_COLUMNS)
    wanted = sorted(set(found["hgnc_id"]))
    sides = {name: (list(part.loc[part["side"] == "up", "hgnc_id"]), list(part.loc[part["side"] == "down", "hgnc_id"])) for name, part in found.groupby("signature")}
    published = dict(zip(found["hgnc_id"], found["published_symbol"]))

    coverage, held = [], {}
    for matrix in matrices.itertuples(index=False):
        present = set(pd.read_parquet(C.project_path(matrix.gene_matrix), columns=["hgnc_id"])["hgnc_id"])
        held[(matrix.study_accession, matrix.platform)] = present
        for name, (up, down) in sides.items():
            have_up, have_down = [gene for gene in up if gene in present], [gene for gene in down if gene in present]
            coverage.append(dict(signature=name, source=matrix.source, study_accession=matrix.study_accession, platform=matrix.platform,
                                 value_scale=matrix.value_scale, comparable_scale=bool(matrix.comparable_scale), genes_published=len(up) + len(down),
                                 up_published=len(up), up_present=len(have_up), down_published=len(down), down_present=len(have_down),
                                 coverage=(len(have_up) + len(have_down)) / (len(up) + len(down)),
                                 missing="; ".join(published[gene] for gene in up + down if gene not in present)))
    coverage = pd.DataFrame(coverage, columns=COVERAGE_COLUMNS)

    # The role of every study in the derivation of a signature, and the samples of those studies, which other studies may hold as repeat deposits; a
    # discovery role outranks a validation role, as in the validation stage.
    roles = {name: {**{study: "validation" for study in entry["derivation_studies"]["validation"]},
                    **{study: "discovery" for study in entry["derivation_studies"]["discovery"]}} for name, entry in definitions.items()}
    derived = {name: set(samples.loc[samples["study_accession"].isin(role_of), "sample_id"]) for name, role_of in roles.items()}
    tissues = scope_tissues()
    tests, scores, read, lowest = [], [], {}, {}
    minimum = policy["minimum_group"]
    for comparison, rule in policy["comparisons"].items():
        positive, negative = comparison_members(samples, rule)
        chosen = [entry for entry in policy["definitions"] if comparison in entry["intended_comparison"]]
        if not chosen:
            continue
        members = samples[positive | negative].assign(label=positive[positive | negative].astype("int8"))
        for (study, platform, tissue), part in members.groupby(STRATUM, sort=True):
            if part["label"].nunique() < 2:
                continue
            first, path = part.iloc[0], C.project_path(part.iloc[0]["gene_matrix"])
            cases = part[part["label"] == 1].sort_values("sample_accession", kind="stable")
            others = part[part["label"] == 0].sort_values("sample_accession", kind="stable")
            for entry in chosen:
                name = entry["name"]
                up, down = sides[name]
                present = held[(study, platform)]
                have_up, have_down = [gene for gene in up if gene in present], [gene for gene in down if gene in present]
                share = (len(have_up) + len(have_down)) / (len(up) + len(down))
                sided = all(len(have) >= entry["minimum_genes_per_side"] for have, named in ((have_up, up), (have_down, down)) if named)
                role = roles[name].get(study, "")
                shared = 0 if role else int(part["sample_id"].isin(derived[name]).sum())
                repeats = int((~part["is_primary_record"]).sum())
                for variant, rule_of_weight in variants(entry):
                    reason = ("definition_not_verified" if not entry["verified"] else "scale_not_comparable" if not first["comparable_scale"]
                              else "fewer_than_minimum_of_a_class" if min(len(cases), len(others)) < minimum
                              else "gene_coverage_below_threshold" if share < entry.get("minimum_gene_coverage", policy["minimum_gene_coverage"])
                              else "too_few_genes_on_a_side" if not sided else "")
                    row = dict(signature=name, label=entry["label"], type=entry["type"], variant=variant, primary_variant=variant == "primary",
                               comparison=comparison, primary_comparison=entry["intended_comparison"][0] == comparison, source=first["source"],
                               study_accession=study, platform=platform, tissue=tissue, value_scale=first["value_scale"], positives=len(cases), negatives=len(others),
                               scored_positives=0, scored_negatives=0, patients_positive=0, patients_negative=0, genes_up_used=len(have_up),
                               genes_down_used=len(have_down), coverage=share, down_weight=np.nan, score_shift=np.nan, status="skipped", reason=reason,
                               auroc=np.nan, auroc_low=np.nan, auroc_high=np.nan, derivation_role=role, samples_shared_with_derivation=shared,
                               repeat_deposits=repeats, in_scope_group=tissue in tissues, independent=not role and shared == 0 and repeats == 0,
                               doi=entry["doi"], pmid=entry["pmid"])
                    if not reason:
                        if path not in read:
                            read[path] = log_values(read_genes(path, wanted), first["value_scale"])
                        values, shift = read[path], 0.0
                        if entry["formula"] == "geometric_mean_difference":
                            if path not in lowest:
                                lowest[path] = matrix_minimum(path, first["value_scale"])
                            # As in the validation stage: values are raised so that the lowest is one where any lies below one, since the logarithm of
                            # a value at or below zero is undefined and that of a value below one negative.
                            shift = 1.0 - lowest[path] if lowest[path] < 1 else 0.0
                        weight = down_weight(rule_of_weight, have_up, have_down)
                        score = signature_score(values[list(cases["sample_accession"]) + list(others["sample_accession"])], have_up, have_down,
                                                entry["formula"], weight, shift)
                        score = (score if entry["direction"] == "higher_in_positive" else -score).dropna()
                        own = part.set_index("sample_accession").loc[score.index].rename_axis("sample_accession")
                        hits, misses = own["label"].eq(1).to_numpy(), own["label"].eq(0).to_numpy()
                        found_scores = score.to_numpy(dtype="float64")
                        row.update(scored_positives=int(hits.sum()), scored_negatives=int(misses.sum()), down_weight=weight,
                                   score_shift=shift if entry["formula"] == "geometric_mean_difference" else np.nan)
                        if min(hits.sum(), misses.sum()) < minimum:
                            row["reason"] = "fewer_than_minimum_with_a_score"
                        else:
                            patients = int(own.loc[hits, "patient"].nunique()), int(own.loc[misses, "patient"].nunique())
                            low, high = np.nan, np.nan
                            if min(patients) >= minimum:
                                low, high = cluster_interval(found_scores[hits], found_scores[misses], own.loc[hits, "patient"].to_numpy(),
                                                             own.loc[misses, "patient"].to_numpy(), policy["bootstrap"], policy["bootstrap"]["seed"])
                            row.update(status="tested", patients_positive=patients[0], patients_negative=patients[1],
                                       auroc=auroc(found_scores[hits], found_scores[misses]), auroc_low=low, auroc_high=high)
                            scores.append(own.reset_index()[[*C.KEYS, "tissue", "label", "patient"]].assign(
                                signature=name, variant=variant, comparison=comparison, score=found_scores))
                    tests.append(row)
        done = pd.DataFrame([row for row in tests if row["comparison"] == comparison and row["variant"] == "primary"])
        print(f"[04] {comparison}: " + "; ".join(f"{name} {int((part['status'] == 'tested').sum())} of {len(part)} strata tested"
                                                 for name, part in done.groupby("signature", sort=False)), flush=True)
    tests = pd.DataFrame(tests, columns=TEST_COLUMNS)
    scores = pd.concat(scores, ignore_index=True)[SCORE_COLUMNS]
    C.require(set(tests["reason"]) <= set(REASONS) | {""}, "Unknown reason")
    agreement = stage_05_checks(found, tests, scores)
    summary = pd.DataFrame(summarise(tests, definitions, tissues), columns=SUMMARY_COLUMNS)

    order = ["signature", "variant", "comparison", *STRATUM]
    tests = tests.sort_values(order, kind="stable").reset_index(drop=True)
    scores = scores.sort_values([*order[:3], "study_accession", "platform", "tissue", "sample_accession"], kind="stable").reset_index(drop=True)
    coverage = coverage.sort_values(["signature", "study_accession", "platform"], kind="stable").reset_index(drop=True)
    summary = summary.sort_values(["signature", "variant", "comparison", "scope", "strata_set"], kind="stable").reset_index(drop=True)
    outputs = {name: C.DATA / f"{name}.parquet" for name in ("signature_genes", "signature_coverage", "signature_tests", "signature_scores", "signature_summary")}
    for name, frame in zip(outputs, (found, coverage, tests, scores, summary)):
        C.save_table(frame, outputs[name])
    C.save_review(summary, C.OUT / "signature_summary_review.csv")
    headline = {}
    for row in summary[summary["primary_variant"]].itertuples(index=False):
        headline.setdefault(row.signature, {}).setdefault(row.comparison, {}).setdefault(row.scope, {})[row.strata_set] = dict(
            strata=int(row.strata), tested=int(row.tested), auroc_median=C.rounded(row.auroc_median), auroc_q1=C.rounded(row.auroc_q1),
            auroc_q3=C.rounded(row.auroc_q3), share_above_limit=C.rounded(row.share_above_limit),
            skipped={reason: int(getattr(row, f"skipped_{reason}")) for reason in REASONS if getattr(row, f"skipped_{reason}")})
    C.finish(4, start, [*outputs.values(), C.OUT / "signature_summary_review.csv"], method=policy["method"], reasons=REASONS,
             signatures={name: dict(label=entry["label"], type=entry["type"], citation=entry["citation"], doi=entry["doi"], pmid=entry["pmid"],
                                    intended_comparison=entry["intended_comparison"], formula=entry["formula"],
                                    variants={variant: weight for variant, weight in variants(entry)}, genes=len(entry["genes_up"]) + len(entry["genes_down"]),
                                    derivation_studies=entry["derivation_studies"]) for name, entry in definitions.items()},
             tests=dict(total=len(tests), tested=int((tests["status"] == "tested").sum()), skipped_by_reason=C.counts(tests.loc[tests["reason"] != "", "reason"])),
             summary=headline, auroc_report_limit=policy["auroc_report_limit"], stage_05_agreement=agreement, scope=dict(group=policy["scope_group"], tissues=tissues),
             limits={name: policy[name] for name in ("minimum_group", "minimum_gene_coverage", "sepsis_phases", "bootstrap")}, protocol=C.CONFIG["protocol"])
    print("[04] Complete: " + "; ".join(f"{name} {comparison} median AUROC {entry[WHOLE_SCOPE]['all']['auroc_median']} over "
                                        f"{entry[WHOLE_SCOPE]['all']['tested']} strata" for name, by in headline.items() for comparison, entry in by.items())
          + f"; MetaScore agrees with the validation stage in {agreement['strata']} strata, {agreement['independent_strata']} of them independent; review sheet: "
          f"{(C.OUT / 'signature_summary_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
