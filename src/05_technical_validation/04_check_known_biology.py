"""Check inside single studies the direction of published sepsis genes between sepsis and healthy controls, and the marker genes of sorted blood cells."""
import importlib

C = importlib.import_module("01_index_validation_inputs")

GENE_COLUMNS = ["source", "study_accession", "platform", "tissue", "test", "role", "published_symbol", "hgnc_id", "symbol", "present", "sepsis", "healthy_control",
                "mean_sepsis", "mean_healthy_control", "difference", "standardised_difference", "as_published"]
STRATUM_COLUMNS = ["source", "study_accession", "platform", "tissue", "value_scale", "sepsis", "healthy_control", "sepsis_left_out_by_phase", "genes_present",
                   "genes_as_published", "score_genes_higher", "score_genes_lower", "score_computed", "score_shift", "score_difference", "auroc", "auroc_low",
                   "auroc_high", "patients_sepsis", "patients_healthy_control", "auroc_one_per_patient", "sepsis_one_per_patient", "healthy_control_one_per_patient",
                   "repeat_deposits", "score_derivation", "samples_shared_with_derivation"]
SCORE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "tissue", "sample_class", "patient", "one_per_patient", "score"]
CELL_COLUMNS = ["source", "study_accession", "platform", "cell_type", "markers_read", "samples", "other_samples", "other_fractions", "mean_in_cell_type",
                "mean_in_others", "auroc"]
GROUPS = ("sepsis", "healthy_control")


def geometric_mean(values):
    """The geometric mean of every column over the rows that hold a value, as MetaIntegrator computes it."""
    import numpy as np

    return np.exp(np.log(values).mean(axis=0, skipna=True))


def metascore(values, up, down, shift):
    """The Sepsis MetaScore of every sample: the geometric mean of the genes published as higher in sepsis minus that of the genes published as lower,
    on values raised so that the lowest value of the matrix is one where any is below zero."""
    return geometric_mean(values.loc[up] + shift) - geometric_mean(values.loc[down] + shift)


def matrix_minimum(path, scale):
    """The lowest value of a gene matrix on the log2 scale of the tests."""
    import pandas as pd

    frame = pd.read_parquet(path).drop(columns=["hgnc_id", "symbol"]).astype("float64")
    return float(C.log_values(frame, scale).min().min())


def cluster_interval(cases, others, case_groups, other_groups, policy, seed):
    """A percentile interval of the area under the curve from a bootstrap that resamples patients, with all of their samples, within each group."""
    import numpy as np

    generator = np.random.default_rng(seed)
    case_sets = [np.flatnonzero(case_groups == name) for name in np.unique(case_groups)]
    other_sets = [np.flatnonzero(other_groups == name) for name in np.unique(other_groups)]
    found = []
    for _ in range(policy["replicates"]):
        first = np.concatenate([case_sets[index] for index in generator.integers(0, len(case_sets), len(case_sets))])
        second = np.concatenate([other_sets[index] for index in generator.integers(0, len(other_sets), len(other_sets))])
        found.append(C.auroc(cases[first], others[second]))
    tail = (1 - policy["level"]) / 2
    return float(np.quantile(found, tail)), float(np.quantile(found, 1 - tail))


def one_per_patient(own):
    """Mark, among the samples of one stratum, the one of every patient that is earliest by its time order, then lowest in accession."""
    timed = own["time_order"].notna()
    ranked = sorted(zip(own["patient"], ~timed, own["time_order"].where(timed, 0.0), map(C.accession_order, own.index), own.index))
    chosen = {}
    for patient, *_, sample in ranked:
        chosen.setdefault(patient, sample)
    return own.index.isin(list(chosen.values()))


def main():
    import numpy as np
    import pandas as pd

    print(f"[04] Checking known biology v{C.VERSION}", flush=True)
    start = C.begin(4)
    policy = C.CONFIG["biology"]
    samples = pd.read_parquet(C.DATA / "validation_samples.parquet")
    matrices = pd.read_parquet(C.DATA / "validation_matrices.parquet")
    genes = pd.read_parquet(C.DATA / "validation_genes.parquet")
    # Patients: an identified patient, or else the distinct sample, so that the repeated samples of one patient are resampled together.
    samples["patient"] = samples["patient_key"].where(samples["patient_key"] != "", samples["sample_id"])
    # The studies from which the score was derived, and the distinct samples they hold, which other studies may share.
    derivation = policy["score_derivation"]
    role_of = {**{study: "validation" for study in derivation["validation"]}, **{study: "discovery" for study in derivation["discovery"]}}
    derived = set(samples.loc[samples["study_accession"].isin(role_of), "sample_id"])
    samples["in_group"] = samples["sample_class"].eq("healthy_control") | (samples["sample_class"].eq("sepsis") & samples["sepsis_phase"].isin(policy["sepsis_phases"]))
    named = genes[genes["test"].isin(["sepsis_score", "single_gene"])]
    higher = list(genes.loc[(genes["test"] == "sepsis_score") & (genes["role"] == "higher_in_sepsis"), "hgnc_id"])
    lower = list(genes.loc[(genes["test"] == "sepsis_score") & (genes["role"] == "lower_in_sepsis"), "hgnc_id"])
    markers = {cell: list(part["hgnc_id"]) for cell, part in genes[genes["test"] == "cell_marker"].groupby("role")}
    gene_rows, strata, scores, cells = [], [], [], []
    for matrix in matrices[matrices["comparable_scale"]].itertuples(index=False):
        part = samples[(samples["study_accession"] == matrix.study_accession) & (samples["platform"] == matrix.platform)]
        counts = part[part["in_group"]].groupby(["tissue", "sample_class"]).size().unstack(fill_value=0).reindex(columns=list(GROUPS), fill_value=0)
        eligible = [tissue for tissue, row in counts.iterrows() if min(row["sepsis"], row["healthy_control"]) >= policy["minimum_group"]]
        fractions = part[part["tissue"].isin(policy["cell_fractions"])].groupby("tissue").size()
        fractions = fractions[fractions >= policy["minimum_cells"]]
        tested_cells = [cell for cell in fractions.index if cell in markers] if len(fractions) >= 2 else []
        if not eligible and not tested_cells:
            continue
        wanted = set(named["hgnc_id"]) | {gene for cell in tested_cells for gene in markers[cell]}
        values = C.log_values(C.read_genes(C.project_path(matrix.gene_matrix), wanted), matrix.value_scale)
        base = dict(source=matrix.source, study_accession=matrix.study_accession, platform=matrix.platform)
        shift = None
        for tissue in eligible:
            stratum = part[(part["tissue"] == tissue) & part["in_group"]]
            left_out = int(((part["tissue"] == tissue) & part["sample_class"].eq("sepsis") & ~part["in_group"]).sum())
            sepsis = list(stratum.loc[stratum["sample_class"] == "sepsis", "sample_accession"])
            healthy = list(stratum.loc[stratum["sample_class"] == "healthy_control", "sample_accession"])
            agreed = 0
            for gene in named.itertuples(index=False):
                row = dict(**base, tissue=tissue, test=gene.test, role=gene.role, published_symbol=gene.published_symbol, hgnc_id=gene.hgnc_id, symbol=gene.symbol,
                           present=gene.hgnc_id in values.index, sepsis=0, healthy_control=0, mean_sepsis=np.nan, mean_healthy_control=np.nan, difference=np.nan,
                           standardised_difference=np.nan, as_published=False)
                if row["present"]:
                    first, second = values.loc[gene.hgnc_id, sepsis].dropna(), values.loc[gene.hgnc_id, healthy].dropna()
                    if len(first) >= 2 and len(second) >= 2:
                        difference = float(first.mean() - second.mean())
                        pooled = float(np.sqrt(((len(first) - 1) * first.var() + (len(second) - 1) * second.var()) / (len(first) + len(second) - 2)))
                        row.update(sepsis=len(first), healthy_control=len(second), mean_sepsis=float(first.mean()), mean_healthy_control=float(second.mean()),
                                   difference=difference, standardised_difference=difference / pooled if pooled > 0 else np.nan,
                                   as_published=difference > 0 if gene.role == "higher_in_sepsis" else difference < 0)
                    else:
                        row["present"] = False
                agreed += row["as_published"]
                gene_rows.append(row)
            up, down = [gene for gene in higher if gene in values.index], [gene for gene in lower if gene in values.index]
            shared = stratum["sample_id"].isin(derived) & ~stratum["study_accession"].isin(role_of)
            summary = dict(**base, tissue=tissue, value_scale=matrix.value_scale, sepsis=len(sepsis), healthy_control=len(healthy), sepsis_left_out_by_phase=left_out,
                           genes_present=sum(row["present"] for row in gene_rows[-len(named):]), genes_as_published=int(agreed), score_genes_higher=len(up),
                           score_genes_lower=len(down), score_computed=min(len(up), len(down)) >= policy["minimum_score_genes"], score_shift=np.nan,
                           score_difference=np.nan, auroc=np.nan, auroc_low=np.nan, auroc_high=np.nan, patients_sepsis=0, patients_healthy_control=0,
                           auroc_one_per_patient=np.nan, sepsis_one_per_patient=0, healthy_control_one_per_patient=0,
                           repeat_deposits=int((~stratum["is_primary_record"]).sum()), score_derivation=role_of.get(matrix.study_accession, ""),
                           samples_shared_with_derivation=int(shared.sum()))
            if summary["score_computed"]:
                if shift is None:
                    lowest = matrix_minimum(C.project_path(matrix.gene_matrix), matrix.value_scale)
                    # MetaIntegrator raises negative values so that the lowest is one; a lowest value between zero and one is raised the same way, since its
                    # logarithm would otherwise be negative or undefined.
                    shift = 1.0 - lowest if lowest < 1 else 0.0
                score = metascore(values[sepsis + healthy], up, down, shift).dropna()
                own = stratum.set_index("sample_accession").loc[score.index]
                cases, others = own["sample_class"].eq("sepsis").to_numpy(), own["sample_class"].eq("healthy_control").to_numpy()
                values_of = score.to_numpy()
                patients = int(own.loc[cases, "patient"].nunique()), int(own.loc[others, "patient"].nunique())
                # A percentile interval over one or two patients collapses, so it is computed only where both groups hold enough patients.
                low, high = np.nan, np.nan
                if min(patients) >= policy["minimum_group"]:
                    low, high = cluster_interval(values_of[cases], values_of[others], own.loc[cases, "patient"].to_numpy(),
                                                 own.loc[others, "patient"].to_numpy(), policy["bootstrap"], policy["bootstrap"]["seed"])
                single = one_per_patient(own)
                summary.update(score_shift=shift, score_difference=float(values_of[cases].mean() - values_of[others].mean()),
                               auroc=C.auroc(values_of[cases], values_of[others]), auroc_low=low, auroc_high=high,
                               patients_sepsis=patients[0], patients_healthy_control=patients[1],
                               sepsis_one_per_patient=int((cases & single).sum()), healthy_control_one_per_patient=int((others & single).sum()))
                if (cases & single).sum() >= 2 and (others & single).sum() >= 2:
                    summary["auroc_one_per_patient"] = C.auroc(values_of[cases & single], values_of[others & single])
                scores.extend(dict(**base, sample_accession=sample, tissue=tissue, sample_class=row.sample_class, patient=row.patient,
                                   one_per_patient=bool(marked), score=float(value))
                              for (sample, row), marked, value in zip(own.iterrows(), single, values_of))
            strata.append(summary)
        for cell in tested_cells:
            read = [gene for gene in markers[cell] if gene in values.index]
            if not read:
                continue
            score = values.loc[read].mean(axis=0, skipna=True)
            own = score[list(part.loc[part["tissue"] == cell, "sample_accession"])].dropna()
            rest = part[part["tissue"].isin(fractions.index) & (part["tissue"] != cell)]
            other = score[list(rest["sample_accession"])].dropna()
            cells.append(dict(**base, cell_type=cell, markers_read=len(read), samples=len(own), other_samples=len(other),
                              other_fractions="; ".join(sorted(set(rest["tissue"]))), mean_in_cell_type=float(own.mean()), mean_in_others=float(other.mean()),
                              auroc=C.auroc(own, other)))
    gene_rows = pd.DataFrame(gene_rows, columns=GENE_COLUMNS)
    strata = pd.DataFrame(strata, columns=STRATUM_COLUMNS)
    scores = pd.DataFrame(scores, columns=SCORE_COLUMNS)
    cells = pd.DataFrame(cells, columns=CELL_COLUMNS)
    C.require(len(strata) > 0 and len(gene_rows) == len(strata) * len(named), "Every study platform and tissue with both groups needs one row for every gene")
    C.save_table(gene_rows, C.DATA / "sepsis_genes_by_study.parquet")
    C.save_table(strata, C.DATA / "sepsis_score_by_study.parquet")
    C.save_table(scores, C.DATA / "sepsis_score_samples.parquet")
    C.save_table(cells, C.DATA / "cell_markers_by_study.parquet")
    C.save_review(strata, C.OUT / "sepsis_score_review.csv")
    C.save_review(cells, C.OUT / "cell_markers_review.csv")
    tested = gene_rows[gene_rows["present"]]
    by_gene = {str(symbol): dict(role=str(part["role"].iloc[0]), symbol_now=str(part["symbol"].iloc[0]), tested=len(part), as_published=int(part["as_published"].sum()))
               for symbol, part in tested.groupby("published_symbol", sort=False)}
    scored = strata[strata["score_computed"]]
    # A stratum is independent of the derivation of the score when its study took no part in it and it shares no sample with such a study or another study.
    independent = scored[(scored["score_derivation"] == "") & (scored["samples_shared_with_derivation"] == 0) & (scored["repeat_deposits"] == 0)]
    report = policy["report"]
    without = scored[scored["auroc_low"].isna()]
    change = (scored["auroc_one_per_patient"] - scored["auroc"]).abs()
    C.finish(4, start, [C.DATA / "sepsis_genes_by_study.parquet", C.DATA / "sepsis_score_by_study.parquet", C.DATA / "sepsis_score_samples.parquet",
                        C.DATA / "cell_markers_by_study.parquet", C.OUT / "sepsis_score_review.csv", C.OUT / "cell_markers_review.csv"], method=policy["method"],
             strata=dict(total=len(strata), study_platforms=int(strata.groupby(["study_accession", "platform"]).ngroups), sepsis=int(strata["sepsis"].sum()),
                         healthy_control=int(strata["healthy_control"].sum()),
                         by_tissue={str(name): int(count) for name, count in strata["tissue"].value_counts().sort_index().items()}),
             genes=dict(comparisons=len(tested), as_published=int(tested["as_published"].sum()),
                        share_as_published=round(float(tested["as_published"].mean()), 6), by_gene=by_gene),
             score=dict(computed=len(scored), auroc_lowest=round(float(scored["auroc"].min()), 6), auroc_median=round(float(scored["auroc"].median()), 6),
                        **{C.at_least("auroc_at_least", limit): int((scored["auroc"] >= limit).sum()) for limit in report["auroc_at_least"]},
                        **{C.at_least("auroc_below", limit): int((scored["auroc"] < limit).sum()) for limit in report["auroc_below"]},
                        **{C.at_least("interval_lower_bound_at_least", limit): int((scored["auroc_low"] >= limit).sum())
                           for limit in report["interval_lower_bound_at_least"]},
                        intervals=dict(computed=int(scored["auroc_low"].notna().sum()), not_computed=len(without),
                                       reason_not_computed=f"fewer than biology.minimum_group ({policy['minimum_group']}) patients in a group",
                                       strata_without_an_interval=sorted(without["study_accession"] + "|" + without["platform"] + "|" + without["tissue"])),
                        sepsis_left_out_by_phase=int(strata["sepsis_left_out_by_phase"].sum()),
                        independent=dict(strata=len(independent), auroc_median=round(float(independent["auroc"].median()), 6) if len(independent) else None,
                                         auroc_lowest=round(float(independent["auroc"].min()), 6) if len(independent) else None),
                        by_derivation={str(name or "not_used"): int(count) for name, count in scored["score_derivation"].value_counts().items()},
                        largest_change_with_one_sample_per_patient=round(float(change.max()), 6) if change.notna().any() else None,
                        one_per_patient_not_computed=int(change.isna().sum())),
             score_derivation=derivation, bootstrap=policy["bootstrap"], sepsis_phases=policy["sepsis_phases"],
             cell_markers=[dict(study=row.study_accession, cell_type=row.cell_type, samples=int(row.samples), other_samples=int(row.other_samples),
                                auroc=round(float(row.auroc), 6)) for row in cells.itertuples(index=False)],
             limits={name: policy[name] for name in ("minimum_group", "minimum_score_genes", "minimum_cells")}, protocol=C.CONFIG["protocol"])
    print(f"[04] Complete: {len(strata)} study platforms and tissues with sepsis and healthy controls; {int(tested['as_published'].sum()):,} of {len(tested):,} "
          f"gene comparisons have the published direction; score AUROC median {scored['auroc'].median():.3f}, lowest {scored['auroc'].min():.3f}; "
          f"{len(cells)} cell-marker comparisons; review sheets: {(C.OUT / 'sepsis_score_review.csv').relative_to(C.ROOT).as_posix()}, "
          f"{(C.OUT / 'cell_markers_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
