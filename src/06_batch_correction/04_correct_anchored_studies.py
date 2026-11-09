"""Correct the units with healthy controls on those controls, as COCONUT does, after removing the effect of sex, and build the references of every sample class."""
import importlib

C = importlib.import_module("01_select_correction_inputs")

MODEL_COLUMNS = ["group", "hgnc_id", "symbol", "sex_effect", "sex_units", "anchor_mean", "anchor_sd", "anchor_units"]
EFFECT_COLUMNS = ["group", "study_accession", "platform", "method", "hgnc_id", "location", "scale"]
REFERENCE_COLUMNS = ["group", "sample_class", "age_group", "hgnc_id", "mean", "variance", "units"]


def sex_effect(group, units, records, order):
    """The pooled slope of every gene on the female indicator within the cells of unit, sample class and age group that hold both sexes, over the estimation
    samples of known sex of the units that hold primary records, and the number of units that inform it. Centring within cells keeps a difference between
    the shares of females of the classes of a unit out of the slope, where it would stand for the difference between the classes."""
    import numpy as np

    products, squares, informative = np.zeros(len(order)), np.zeros(len(order)), np.zeros(len(order), dtype=int)
    for unit in units[(units["group"] == group) & (units["status"] == "in_unit") & (units["role"] == "merged")].itertuples(index=False):
        own = C.unit_records(records, group, unit.study_accession, unit.platform)
        own = own[own["estimation"] & own["sex_used"].isin(["expression", "recorded"])].reset_index(drop=True)
        cells = [index for index in C.class_cells(own, unit.study_accession).values() if own.loc[index, "female"].nunique() == 2]
        if not cells:
            continue
        values = C.read_aligned(group, unit.study_accession, unit.platform, own["sample_accession"])
        spread = np.zeros(len(order))
        for index in cells:
            part = values[:, index]
            present = ~np.isnan(part)
            held = present.any(axis=1)[:, None]
            female = np.where(present, own.loc[index, "female"].to_numpy()[None, :], np.nan)
            centred_sex = female - np.nanmean(np.where(held, female, 0.0), axis=1)[:, None]
            centred = part - np.nanmean(np.where(held, part, 0.0), axis=1)[:, None]
            products += np.nansum(centred_sex * centred, axis=1)
            spread += np.nansum(centred_sex ** 2, axis=1)
        squares += spread
        informative += (spread > 0).astype(int)
    return np.divide(products, squares, out=np.zeros(len(order)), where=squares > 0), informative


def main():
    import numpy as np
    import pandas as pd

    print(f"[04] Correcting anchored studies v{C.VERSION}", flush=True)
    start = C.begin(4)
    anchor = C.CONFIG["anchored"]["anchor_class"]
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    models, effects, references, summary = [], [], [], {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        symbols = genes[genes["group"] == group].set_index("hgnc_id")["symbol"]
        sex, sex_units = sex_effect(group, units, records, order)
        moments_of = C.unit_moments(group, units, records, sex)
        active = units[(units["group"] == group) & (units["status"] == "in_unit")]
        anchored = active[active["tier"] == "anchored"]
        training = [(unit.study_accession, unit.platform) for unit in anchored.itertuples(index=False) if unit.role == "merged"]
        C.require(len(training) >= C.CONFIG["inputs"]["minimum_reference_units"], f"Too few anchored units in {group}")
        mean, deviation, anchor_units = C.anchor_model([C.gather(moments_of[key], [anchor]) for key in training])
        found = {}
        for unit in anchored.itertuples(index=False):
            key = (unit.study_accession, unit.platform)
            found[key] = C.anchored_effects(C.gather(moments_of[key], [anchor]), mean, deviation)
            effects.append(pd.DataFrame({"group": group, "study_accession": unit.study_accession, "platform": unit.platform, "method": "healthy_anchored",
                                         "hgnc_id": order, "location": found[key][0], "scale": found[key][1]}))
        reference = C.class_reference(moments_of, training, mean, deviation, found)
        for (name, age), (centre, spread, count) in reference.items():
            references.append(pd.DataFrame({"group": group, "sample_class": name, "age_group": age, "hgnc_id": order, "mean": centre, "variance": spread,
                                            "units": count}))
        models.append(pd.DataFrame({"group": group, "hgnc_id": order, "symbol": [symbols[gene] for gene in order], "sex_effect": sex, "sex_units": sex_units,
                                    "anchor_mean": mean, "anchor_sd": deviation, "anchor_units": anchor_units}))
        # After the correction the healthy controls of every training unit sit on the anchor: mean near zero and variance near one in standard units.
        settled = [C.corrected_moments(C.gather(moments_of[key], [anchor]), mean, deviation, *found[key]) for key in training]
        centres = np.vstack([part[1] for part in settled])
        scales = np.vstack([found[key][1] for key in found])
        healthy = reference.get((anchor, "all"))
        floor = C.relative_floor(deviation)
        summary[group] = dict(anchored_units=len(anchored), training_units=len(training), check_only_units=int((anchored["role"] == "check_only").sum()),
                              anchor_samples=int(sum(np.nanmax(C.gather(moments_of[key], [anchor])[0]) for key in training)),
                              genes=len(order), genes_with_anchor_model=int(np.isfinite(deviation).sum()),
                              anchor_sd_floor=round(floor, 6), genes_at_the_anchor_sd_floor=int((deviation <= floor).sum()),
                              genes_with_sex_effect_over_one=int((np.abs(sex) > 1).sum()),
                              largest_sex_effects={symbols[order[index]]: round(float(sex[index]), 4) for index in np.argsort(-np.abs(sex))[:5]},
                              median_absolute_healthy_mean_after=round(float(np.nanmedian(np.abs(centres))), 6),
                              median_scale=round(float(np.nanmedian(scales)), 6),
                              reference_cells=sorted(f"{name}|{age}" for name, age in reference),
                              healthy_reference_genes=int(np.isfinite(healthy[0]).sum()) if healthy else 0)
        print(f"[04] {group}: {len(anchored)} anchored units corrected on their healthy controls ({len(training)} build the model); "
              f"{int(np.isfinite(deviation).sum()):,} of {len(order):,} genes have an anchor model, {summary[group]['genes_at_the_anchor_sd_floor']:,} of them "
              f"at the floor of the standard deviation ({floor:.4f}); median scale {summary[group]['median_scale']:.3f}", flush=True)
    models = pd.concat(models, ignore_index=True)[MODEL_COLUMNS]
    effects = pd.concat(effects, ignore_index=True)[EFFECT_COLUMNS]
    references = pd.concat(references, ignore_index=True)[REFERENCE_COLUMNS]
    C.require(np.isfinite(models["sex_effect"]).all() and (models["anchor_sd"].dropna() > 0).all(), "Invalid anchor model")
    C.require((effects["scale"].dropna() > 0).all(), "Every scale must be positive")
    C.save_table(models, C.DATA / "correction_model.parquet")
    C.save_table(effects, C.DATA / "anchored_effects.parquet")
    C.save_table(references, C.DATA / "class_references.parquet")
    C.finish(4, start, [C.DATA / "correction_model.parquet", C.DATA / "anchored_effects.parquet", C.DATA / "class_references.parquet"],
             method=C.CONFIG["anchored"]["method"], groups=summary, limits={name: C.CONFIG["anchored"][name] for name in ("floor_fraction", "minimum_prior_genes")},
             protocol=C.CONFIG["protocol"])
    print("[04] Complete", flush=True)


if __name__ == "__main__":
    main()
