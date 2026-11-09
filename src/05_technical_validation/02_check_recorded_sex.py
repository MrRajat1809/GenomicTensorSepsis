"""Read the sex of every sample from its own expression, inside its gene matrix, and compare it with the recorded sex."""
import importlib

C = importlib.import_module("01_index_validation_inputs")

SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "recorded_sex", "female_gene_value", "male_genes_value", "male_genes_read",
                  "difference", "expression_sex", "agreement", "recorded_sex_unreliable"]
MATRIX_COLUMNS = ["source", "study_accession", "platform", "samples", "status", "female_gene_read", "male_genes_read",
                  "lower_group", "higher_group", "gap", "separation", "midpoint",
                  "called_female", "called_male", "ambiguous", "with_recorded_sex", "concordant", "discordant", "concordance", "recorded_sex_unreliable"]
UNIT_FIELDS = ["study_accession", "platform", "tissue", "samples", "status", "lower_group", "higher_group", "gap", "separation", "male_genes_separation",
               "female_gene_separation", "midpoint"]
RECORDED = ("female", "male")


def two_groups(values):
    """Split values into a lower and a higher group by two-means in one dimension, started at the lowest and the highest value."""
    import numpy as np

    low, high = float(values.min()), float(values.max())
    higher = np.zeros(len(values), dtype=bool)
    for _ in range(1000):
        assigned = np.abs(values - high) < np.abs(values - low)
        if not assigned.any() or assigned.all() or (assigned == higher).all():
            return assigned
        higher, low, high = assigned, float(values[~assigned].mean()), float(values[assigned].mean())
    return higher


def separation(low, high):
    """The distance from the mean of the lower group to that of the higher, in pooled standard deviations of the two groups."""
    import numpy as np

    gap = float(high.mean() - low.mean())
    pooled = float(np.sqrt((((low - low.mean()) ** 2).sum() + ((high - high.mean()) ** 2).sum()) / max(len(low) + len(high) - 2, 1)))
    return gap / pooled if pooled > 0 else (np.inf if gap > 0 else -np.inf)


def split_sexes(known, male, female, policy):
    """Split the samples of one matrix, or of one tissue of it, into two groups by their difference, and return the summary of the split and the higher
    group. The split stands only where the groups lie apart and are large enough, and where the Y-chromosome genes are higher and XIST (where it is read)
    lower in the higher group, each by sex.component_separation pooled standard deviations, so that a few outliers of a one-sex cohort cannot pass."""
    import numpy as np

    result = dict(samples=len(known), status="too_few_samples", lower_group=0, higher_group=0, gap=np.nan, separation=np.nan, male_genes_separation=np.nan,
                  female_gene_separation=np.nan, midpoint=np.nan)
    if len(known) < policy["minimum_samples"]:
        return result, None
    # The recorded sex takes no part: the samples are split by their own expression alone.
    higher = two_groups(known.to_numpy())
    low, high = known[~higher], known[higher]
    result["status"] = "one_group"
    if not len(low) or not len(high):
        return result, None
    result.update(lower_group=len(low), higher_group=len(high), gap=float(high.mean() - low.mean()), separation=separation(low, high),
                  midpoint=float((high.mean() + low.mean()) / 2), male_genes_separation=separation(male[~higher], male[higher]))
    if female is not None:
        result["female_gene_separation"] = separation(-female[~higher], -female[higher])
    smaller = min(len(low), len(high))
    if not (result["separation"] >= policy["minimum_separation"] and smaller >= max(policy["minimum_group"], policy["minimum_group_share"] * len(known))):
        return result, None
    if not (result["male_genes_separation"] >= policy["component_separation"]
            and (female is None or result["female_gene_separation"] >= policy["component_separation"])):
        result["status"] = "components_not_separated"
        return result, None
    result["status"] = "tested"
    return result, higher


def main():
    import numpy as np
    import pandas as pd

    print(f"[02] Checking recorded sex v{C.VERSION}", flush=True)
    start = C.begin(2)
    policy = C.CONFIG["sex"]
    samples = pd.read_parquet(C.DATA / "validation_samples.parquet")
    matrices = pd.read_parquet(C.DATA / "validation_matrices.parquet")
    genes = pd.read_parquet(C.DATA / "validation_genes.parquet")
    female = genes.loc[(genes["test"] == "sex") & (genes["role"] == "female"), "hgnc_id"].iloc[0]
    male = list(genes.loc[(genes["test"] == "sex") & (genes["role"] == "male"), "hgnc_id"])
    rows, summaries, units = [], [], []
    for matrix in matrices.itertuples(index=False):
        part = samples[(samples["study_accession"] == matrix.study_accession) & (samples["platform"] == matrix.platform)]
        recorded = part.set_index("sample_accession")["sex"]
        tissue_of = part.set_index("sample_accession")["tissue"]
        summary = dict(source=matrix.source, study_accession=matrix.study_accession, platform=matrix.platform, samples=len(part), status="tested",
                       female_gene_read=False, male_genes_read=0,
                       lower_group=0, higher_group=0, gap=np.nan, separation=np.nan, midpoint=np.nan, called_female=0, called_male=0, ambiguous=0,
                       with_recorded_sex=int(recorded.isin(RECORDED).sum()), concordant=0, discordant=0, concordance=np.nan, recorded_sex_unreliable=False)
        table = pd.DataFrame(dict(female_gene_value=np.nan, male_genes_value=np.nan, difference=np.nan), index=recorded.index)
        called = pd.Series("not_called", index=recorded.index)
        if not matrix.comparable_scale:
            summary["status"] = "scale_not_comparable"
        else:
            values = C.log_values(C.read_genes(C.project_path(matrix.gene_matrix), [female, *male]), matrix.value_scale)
            read = [gene for gene in male if gene in values.index]
            summary["male_genes_read"] = len(read)
            if not read:
                summary["status"] = "genes_missing"
            else:
                summary["female_gene_read"] = female in values.index
                if summary["female_gene_read"]:
                    table["female_gene_value"] = values.loc[female].reindex(table.index)
                table["male_genes_value"] = values.loc[read].mean(axis=0, skipna=True).reindex(table.index)
                table["difference"] = table["male_genes_value"] - (table["female_gene_value"] if summary["female_gene_read"] else 0.0)
                known = table["difference"].dropna()
                # The difference of the two sexes is not the same in every tissue, so a matrix of several tissues is split within each tissue where every
                # tissue holds enough samples, and as a whole otherwise.
                tissues = tissue_of.loc[known.index].value_counts()
                by_tissue = len(tissues) > 1 and tissues.min() >= policy["minimum_samples_per_tissue"]
                for tissue in (sorted(tissues.index) if by_tissue else [""]):
                    values_of = known[tissue_of.loc[known.index] == tissue] if tissue else known
                    result, higher = split_sexes(values_of, table.loc[values_of.index, "male_genes_value"],
                                                 table.loc[values_of.index, "female_gene_value"] if summary["female_gene_read"] else None, policy)
                    units.append(dict(study_accession=matrix.study_accession, platform=matrix.platform, tissue=tissue, **result))
                    if higher is not None:
                        near = (values_of - result["midpoint"]).abs() < policy["ambiguous_zone"] * result["gap"]
                        called.loc[values_of.index] = np.where(near, "ambiguous", np.where(higher, "male", "female"))
                if by_tissue:
                    tested = [unit for unit in units[-len(tissues):] if unit["status"] == "tested"]
                    summary.update(status="tested_by_tissue" if tested else "no_tissue_shows_two_sexes", lower_group=sum(unit["lower_group"] for unit in tested),
                                   higher_group=sum(unit["higher_group"] for unit in tested))
                else:
                    summary.update({name: units[-1][name] for name in ("status", "lower_group", "higher_group", "gap", "separation", "midpoint")})
        agreement = pd.Series("sex_not_recorded", index=recorded.index)
        stated = recorded.isin(RECORDED)
        agreement[stated] = np.where(called[stated] == recorded[stated], "concordant", np.where(called[stated].isin(RECORDED), "discordant", called[stated]))
        summary.update(called_female=int((called == "female").sum()), called_male=int((called == "male").sum()), ambiguous=int((called == "ambiguous").sum()),
                       concordant=int((agreement == "concordant").sum()), discordant=int((agreement == "discordant").sum()))
        if summary["concordant"] + summary["discordant"]:
            summary["concordance"] = summary["concordant"] / (summary["concordant"] + summary["discordant"])
        # A recorded sex that disagrees with the expression of most of a matrix is a property of the deposit, not of single samples.
        summary["recorded_sex_unreliable"] = bool(summary["concordant"] + summary["discordant"] >= policy["minimum_compared"]
                                                  and summary["concordance"] < policy["unreliable_below"])
        summaries.append(summary)
        rows.append(table.assign(source=matrix.source, study_accession=matrix.study_accession, platform=matrix.platform, recorded_sex=recorded,
                                 male_genes_read=summary["male_genes_read"], expression_sex=called, agreement=agreement,
                                 recorded_sex_unreliable=summary["recorded_sex_unreliable"]).rename_axis("sample_accession").reset_index())
    table = pd.concat(rows, ignore_index=True)[SAMPLE_COLUMNS]
    summaries = pd.DataFrame(summaries, columns=MATRIX_COLUMNS)
    C.require(len(table) == len(samples) and not table.duplicated(C.KEYS).any(), "Every sample needs one row of the sex check")
    C.save_table(table, C.DATA / "sex_check_samples.parquet")
    C.save_table(summaries, C.DATA / "sex_check_matrices.parquet")
    C.save_review(summaries, C.OUT / "sex_check_matrices_review.csv")
    C.save_review(table[table["agreement"] == "discordant"], C.OUT / "sex_discordant_samples_review.csv")
    compared = table[table["agreement"].isin(["concordant", "discordant"])]
    reliable = compared[~compared["recorded_sex_unreliable"]]
    tested = summaries[summaries["status"].isin(["tested", "tested_by_tissue"])]
    discordant = table[table["agreement"] == "discordant"]
    # A separation is infinite where both groups hold one value each; the receipt writes it, like a missing one, as null.
    units = [{name: (round(float(value), 6) if np.isfinite(value) else None) if isinstance(value, float) else value for name, value in unit.items()}
             for unit in units]
    C.finish(2, start, [C.DATA / "sex_check_samples.parquet", C.DATA / "sex_check_matrices.parquet", C.OUT / "sex_check_matrices_review.csv",
                        C.OUT / "sex_discordant_samples_review.csv"], method=policy["method"],
             matrices=dict(total=len(summaries), by_status={str(name): int(count) for name, count in summaries["status"].value_counts().sort_index().items()},
                           tested_with_recorded_sex=int((tested["with_recorded_sex"] > 0).sum()),
                           tested_without_the_female_gene=int((~tested["female_gene_read"]).sum())),
             samples=dict(total=len(table), called=int(table["expression_sex"].isin(RECORDED).sum()), ambiguous=int((table["expression_sex"] == "ambiguous").sum()),
                          not_called=int((table["expression_sex"] == "not_called").sum()), with_recorded_sex=int(table["recorded_sex"].isin(RECORDED).sum()),
                          compared=len(compared), concordant=int((compared["agreement"] == "concordant").sum()), discordant=len(discordant),
                          concordance=round(float((compared["agreement"] == "concordant").mean()), 6) if len(compared) else None,
                          called_without_a_recorded_sex=int((table["expression_sex"].isin(RECORDED) & ~table["recorded_sex"].isin(RECORDED)).sum()),
                          compared_outside_unreliable_matrices=len(reliable), discordant_outside_unreliable_matrices=int((reliable["agreement"] == "discordant").sum()),
                          concordance_outside_unreliable_matrices=round(float((reliable["agreement"] == "concordant").mean()), 6) if len(reliable) else None),
             unreliable_matrices=sorted(summaries.loc[summaries["recorded_sex_unreliable"], "study_accession"] + "|"
                                        + summaries.loc[summaries["recorded_sex_unreliable"], "platform"]),
             concordance_by_matrix=dict(lowest=round(float(tested["concordance"].min()), 6), median=round(float(tested["concordance"].median()), 6)),
             discordant_by_study={str(study): int(count) for study, count in discordant.groupby("study_accession").size().items()},
             split_by_tissue=sorted({f"{unit['study_accession']}|{unit['platform']}" for unit in units if unit["tissue"]}),
             units=[{name: unit[name] for name in UNIT_FIELDS} for unit in units],
             limits={name: policy[name] for name in ("minimum_samples", "minimum_separation", "minimum_group", "minimum_group_share", "component_separation",
                                                     "minimum_samples_per_tissue", "ambiguous_zone", "minimum_compared", "unreliable_below")},
             basis=policy["basis"], genes=dict(female=policy["female_gene"], male=policy["male_genes"]), protocol=C.CONFIG["protocol"])
    print(f"[02] Complete: {len(tested)} of {len(summaries)} matrices show two sexes; {len(compared):,} samples with a recorded sex were compared: "
          f"{int((compared['agreement'] == 'concordant').sum()):,} concordant, {len(discordant):,} discordant; review sheets: "
          f"{(C.OUT / 'sex_check_matrices_review.csv').relative_to(C.ROOT).as_posix()}, {(C.OUT / 'sex_discordant_samples_review.csv').relative_to(C.ROOT).as_posix()}",
          flush=True)


if __name__ == "__main__":
    main()
