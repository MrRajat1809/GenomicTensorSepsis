"""Find the genes that sepsis does not change: the standardised difference between sepsis and healthy controls of every gene in every anchored unit."""
import importlib

C = importlib.import_module("01_select_correction_inputs")

EFFECT_COLUMNS = ["group", "study_accession", "platform", "hgnc_id", "sepsis", "healthy_control", "effect", "healthy_mean"]
FLAT_COLUMNS = ["group", "hgnc_id", "symbol", "chromosome", "units", "flatness", "pooled_effect", "healthy_mean", "selected"]


def main():
    import numpy as np
    import pandas as pd

    print(f"[03] Selecting flat genes v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy = C.CONFIG["flat_genes"]
    anchor = C.CONFIG["anchored"]["anchor_class"]
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    effects, flat, summary = [], [], {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        eligible = units[(units["group"] == group) & (units["status"] == "in_unit") & (units["role"] == "merged") & (units["tier"] == "anchored")
                         & (units["sepsis"] >= policy["minimum_group"]) & (units["healthy_anchors"] >= policy["minimum_group"])]
        for unit in eligible.itertuples(index=False):
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            own = own[own["estimation"]]
            sepsis = C.moments(C.read_aligned(group, unit.study_accession, unit.platform, own.loc[own["sample_class"] == "sepsis", "sample_accession"]))
            healthy = C.moments(C.read_aligned(group, unit.study_accession, unit.platform, own.loc[own["sample_class"] == anchor, "sample_accession"]))
            freedom = sepsis[0] + healthy[0] - 2
            pooled = np.sqrt(np.divide(sepsis[2] + healthy[2], freedom, out=np.full(len(order), np.nan), where=freedom > 0))
            correction = 1 - 3 / (4 * (sepsis[0] + healthy[0]) - 9)
            usable = (sepsis[0] >= 2) & (healthy[0] >= 2) & (pooled > 0)
            effect = np.where(usable, (sepsis[1] - healthy[1]) / np.where(usable, pooled, 1) * correction, np.nan)
            effects.append(pd.DataFrame({"group": group, "study_accession": unit.study_accession, "platform": unit.platform, "hgnc_id": order,
                                         "sepsis": sepsis[0].astype(int), "healthy_control": healthy[0].astype(int), "effect": effect, "healthy_mean": healthy[1]}))
        table = pd.concat([part for part in effects if part["group"].iloc[0] == group], ignore_index=True)
        chosen = C.flat_selection(table, genes[genes["group"] == group])
        named = genes[genes["group"] == group].set_index("hgnc_id")
        chosen = chosen.assign(group=group, symbol=chosen["hgnc_id"].map(named["symbol"]), chromosome=chosen["hgnc_id"].map(named["chromosome"]))
        C.require(int(chosen["selected"].sum()) == min(policy["count"], len(chosen)), f"Flat genes of {group} not selected as configured")
        flat.append(chosen[FLAT_COLUMNS])
        selected = chosen[chosen["selected"]]
        summary[group] = dict(units=len(eligible), candidates=len(chosen), selected=len(selected), flatness_highest_selected=round(float(selected["flatness"].max()), 6),
                              flatness_median_candidate=round(float(chosen["flatness"].median()), 6))
        print(f"[03] {group}: {len(selected):,} flat genes of {len(chosen):,} candidates from {len(eligible)} units; their precision-weighted root mean "
              f"square difference between sepsis and healthy controls is at most {selected['flatness'].max():.3f} standard deviations", flush=True)
    effects = pd.concat(effects, ignore_index=True)[EFFECT_COLUMNS]
    flat = pd.concat(flat, ignore_index=True)
    C.save_table(effects, C.DATA / "flat_gene_effects.parquet")
    C.save_table(flat, C.DATA / "flat_genes.parquet")
    C.save_review(flat[flat["selected"]], C.OUT / "flat_genes_review.csv")
    C.finish(3, start, [C.DATA / "flat_gene_effects.parquet", C.DATA / "flat_genes.parquet", C.OUT / "flat_genes_review.csv"], method=policy["method"],
             groups=summary, limits={name: policy[name] for name in ("minimum_group", "count", "minimum_expression_quantile", "excluded_chromosomes")},
             protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: review sheet: {(C.OUT / 'flat_genes_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
