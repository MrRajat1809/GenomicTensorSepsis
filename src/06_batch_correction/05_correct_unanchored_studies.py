"""Correct the units without healthy controls by matched groups, and test it, with flat genes, plain ComBat and no correction, where the answer is known."""
import importlib

C = importlib.import_module("01_select_correction_inputs")

EFFECT_COLUMNS = ["group", "study_accession", "platform", "method", "hgnc_id", "location", "scale"]
FIT_COLUMNS = ["group", "study_accession", "platform", "role", "family", "neonatal", "method", "applicable", "reason", "samples", "genes_corrected"]
HELD_COLUMNS = ["group", "study_accession", "platform", "family", "method", "hidden_healthy", "estimation_samples", "genes", "location_error", "scale_error",
                "pool_dispersion", "sampling_floor", "location_components", "scale_components"]
COMPARISON_COLUMNS = ["group", "method", "baseline", "units", "method_closer", "median_difference", "difference_low", "difference_high", "wilcoxon_p"]
METHODS = ("none", "plain_combat", "flat_genes", "matched_groups")


def unit_features(part, mean, deviation):
    """The location and log scale of every gene of a unit over all of the given samples, in standard units of the anchor model."""
    import numpy as np

    _, location, spread = C.standardised(part, mean, deviation)
    with np.errstate(divide="ignore", invalid="ignore"):
        return location, 0.5 * np.log(spread)


def flat_effects(features, targets, new, flat):
    """Predict a unit's location and scale of every gene from its flat genes, with components chosen by leaving out one training unit at a time."""
    import numpy as np

    locations, log_scales = features
    new_location, new_log_scale = new
    usable = flat & np.isfinite(locations).all(axis=0) & np.isfinite(log_scales).all(axis=0) & np.isfinite(new_location) & np.isfinite(new_log_scale)
    C.require(usable.sum() >= C.CONFIG["flat_genes"]["minimum_shared_genes"], "Too few flat genes held by the unit and every training unit")
    chosen = []
    for train, target, own in ((locations, targets[0], new_location), (log_scales, targets[1], new_log_scale)):
        components, _ = C.choose_components(train[:, usable], target)
        chosen.append((components, C.fit_components(train[:, usable], target, components)(own[usable])))
    (location_components, location), (scale_components, log_scale) = chosen
    return location, np.exp(log_scale), location_components, scale_components


def held_out_error(hidden, pool):
    """How far the hidden healthy controls of a unit land from the healthy controls of the other units after the same correction."""
    import warnings

    import numpy as np

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        centre = np.nanmean(np.vstack([part[1] for part in pool]), axis=0)
        deviation = np.sqrt(np.nanmean(np.vstack([part[2] for part in pool]), axis=0))
        distance = (hidden[1] - centre) / deviation
        ratio = np.log(np.sqrt(hidden[2]) / deviation)
    usable = np.isfinite(distance) & np.isfinite(ratio) & (hidden[0] >= 2)
    return float(np.sqrt(np.mean(distance[usable] ** 2))), float(np.sqrt(np.mean(ratio[usable] ** 2))), int(usable.sum())


def pool_dispersion(pool):
    """How far the healthy controls of each unit of the pool lie from those of the other units after the same correction: the median over the units of
    the root mean square distance in pooled standard deviations within units, the error that the frame of the method has before any unit is held out."""
    import warnings

    import numpy as np

    found = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        deviation = np.sqrt(np.nanmean(np.vstack([part[2] for part in pool]), axis=0))
        for left in range(len(pool)):
            centre = np.nanmean(np.vstack([part[1] for index, part in enumerate(pool) if index != left]), axis=0)
            distance = (pool[left][1] - centre) / deviation
            usable = np.isfinite(distance) & (pool[left][0] >= 2)
            found.append(float(np.sqrt(np.mean(distance[usable] ** 2))))
    return float(np.median(found))


def sampling_floor(values):
    """The error that sampling the hidden healthy controls alone gives: half the root mean square distance between the means of two halves of them, in
    their pooled standard deviation, since the halves differ by twice the standard error of the whole."""
    import warnings

    import numpy as np

    first, second = values[:, 0::2], values[:, 1::2]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        spread = np.sqrt((np.nansum((first - np.nanmean(first, axis=1)[:, None]) ** 2, axis=1) + np.nansum((second - np.nanmean(second, axis=1)[:, None]) ** 2, axis=1))
                         / (np.isfinite(values).sum(axis=1) - 2))
        distance = (np.nanmean(first, axis=1) - np.nanmean(second, axis=1)) / spread
    usable = np.isfinite(distance) & (np.isfinite(first).sum(axis=1) >= 2) & (np.isfinite(second).sum(axis=1) >= 2)
    return float(np.sqrt(np.mean(distance[usable] ** 2)) / 2) if usable.any() else float("nan")


def paired(tested, policy):
    """Compare the location errors of every two methods in the units where both were applied: how often each is closer, the median paired difference with
    a bootstrap interval over units, and the signed-rank test. Each pair is computed once; its row in the reverse order mirrors the same draws."""
    import numpy as np
    from scipy.stats import wilcoxon

    rows = []
    wide = tested.pivot_table(index=["study_accession", "platform"], columns="method", values="location_error")
    generator = np.random.default_rng(policy["seed"])
    tail = (1 - policy["level"]) / 2
    present = [name for name in METHODS if name in wide.columns]
    for position, method in enumerate(present):
        for baseline in present[position + 1:]:
            both = wide[[method, baseline]].dropna()
            if len(both) < policy["minimum_units"]:
                continue
            difference = (both[method] - both[baseline]).to_numpy()
            draws = [float(np.median(difference[generator.integers(0, len(difference), len(difference))])) for _ in range(policy["replicates"])]
            median, low, high = float(np.median(difference)), float(np.quantile(draws, tail)), float(np.quantile(draws, 1 - tail))
            p = float(wilcoxon(difference).pvalue)
            rows.append(dict(method=method, baseline=baseline, units=len(both), method_closer=int((difference < 0).sum()), median_difference=median,
                             difference_low=low, difference_high=high, wilcoxon_p=p))
            rows.append(dict(method=baseline, baseline=method, units=len(both), method_closer=int((difference > 0).sum()), median_difference=-median,
                             difference_low=-high, difference_high=-low, wilcoxon_p=p))
    return rows


def main():
    import numpy as np
    import pandas as pd

    print(f"[05] Correcting unanchored studies v{C.VERSION}", flush=True)
    start = C.begin(5)
    policy, flat_policy = C.CONFIG["unanchored"], C.CONFIG["flat_genes"]
    anchor = C.CONFIG["anchored"]["anchor_class"]
    minimum = C.CONFIG["inputs"]["minimum_unit_samples"]
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    flat_effect_table = pd.read_parquet(C.DATA / "flat_gene_effects.parquet")
    reference_table = pd.read_parquet(C.DATA / "class_references.parquet")
    effects, fits, held, comparisons, summary = [], [], [], [], {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        position = {gene: index for index, gene in enumerate(order)}
        model = C.model_arrays(group, order)
        mean, deviation = model["anchor_mean"], model["anchor_sd"]
        moments_of = C.unit_moments(group, units, records, model["sex_effect"])
        active = units[(units["group"] == group) & (units["status"] == "in_unit")].set_index(["study_accession", "platform"])
        training = [key for key, unit in active.iterrows() if unit["tier"] == "anchored" and unit["role"] == "merged"]
        reference = {}
        for (name, age), part in reference_table[reference_table["group"] == group].groupby(["sample_class", "age_group"]):
            part = part.set_index("hgnc_id").reindex(order)
            reference[(name, age)] = (part["mean"].to_numpy(), part["variance"].to_numpy(), part["units"].to_numpy())
        group_genes = genes[genes["group"] == group]
        group_effects = flat_effect_table[flat_effect_table["group"] == group]

        # Every unit without healthy controls: its class means lined up with the references, unless it is of newborns or has no class to line up.
        for key, unit in active[active["tier"] == "unanchored"].iterrows():
            base = dict(group=group, study_accession=key[0], platform=key[1], role=unit["role"], family=unit["family"], neonatal=bool(unit["neonatal"]),
                        method="matched_groups")
            cells = C.matched_cells(moments_of[key], reference)
            matched, freedom = C.matched_samples(cells)
            if unit["neonatal"] or matched < policy["minimum_matched"] or freedom < policy["minimum_matched"] - 1:
                fits.append(dict(**base, applicable=False, samples=matched, genes_corrected=0,
                                 reason="newborn_study_not_corrected" if unit["neonatal"] else "no_class_to_line_up"))
                continue
            location, scale = C.matched_effects(cells, mean, deviation)
            effects.append(pd.DataFrame({"group": group, "study_accession": key[0], "platform": key[1], "method": "matched_groups", "hgnc_id": order,
                                         "location": location, "scale": scale}))
            fits.append(dict(**base, applicable=True, reason="", samples=matched, genes_corrected=int((np.isfinite(location) & np.isfinite(scale)).sum())))

        # The test where the answer is known: every anchored unit in turn, its healthy controls hidden, corrected from the rest by every method.
        visible_of = {key: {cell: part for cell, part in moments_of[key].items() if cell[0] != anchor} for key in training}
        visible_count = {key: int(np.nanmax(C.combine(cells.values())[0])) if cells else 0 for key, cells in visible_of.items()}
        for held_key in training:
            unit = active.loc[held_key]
            rest = [key for key in training if key != held_key]
            hidden = C.gather(moments_of[held_key], [anchor])
            visible_cells = visible_of[held_key]
            visible = C.combine(visible_cells.values()) if visible_cells else None
            base = dict(group=group, study_accession=held_key[0], platform=held_key[1], family=unit["family"], hidden_healthy=int(np.nanmax(hidden[0])),
                        estimation_samples=visible_count[held_key])
            if base["estimation_samples"] < minimum or len(rest) < C.CONFIG["inputs"]["minimum_reference_units"]:
                continue
            rest_mean, rest_deviation, _ = C.anchor_model([C.gather(moments_of[key], [anchor]) for key in rest])
            rest_effects = {key: C.anchored_effects(C.gather(moments_of[key], [anchor]), rest_mean, rest_deviation) for key in rest}
            pool = [C.corrected_moments(C.gather(moments_of[key], [anchor]), rest_mean, rest_deviation, *rest_effects[key]) for key in rest]
            rows, spread = {}, {}
            raw = [(part[0], part[1], C.variance(part)) for part in (C.gather(moments_of[key], [anchor]) for key in rest)]
            rows["none"] = (*held_out_error((hidden[0], hidden[1], C.variance(hidden)), raw), 0, 0)
            spread["none"] = pool_dispersion(raw)
            spread["flat_genes"] = spread["matched_groups"] = pool_dispersion(pool)
            own = C.unit_records(records, group, *held_key)
            own = own[own["estimation"] & (own["sample_class"] == anchor)]
            controls = C.read_aligned(group, held_key[0], held_key[1], own["sample_accession"]) - model["sex_effect"][:, None] * own["female"].to_numpy()[None, :]
            floor = sampling_floor(controls)
            everything = {key: C.gather(moments_of[key]) for key in rest} | {held_key: visible}
            plain_mean, plain_deviation, _ = C.anchor_model(list(everything.values()))
            plain = {key: C.plain_effects(part, plain_mean, plain_deviation) for key, part in everything.items()}
            plain_pool = [C.corrected_moments(C.gather(moments_of[key], [anchor]), plain_mean, plain_deviation, *plain[key]) for key in rest]
            rows["plain_combat"] = (*held_out_error(C.corrected_moments(hidden, plain_mean, plain_deviation, *plain[held_key]), plain_pool), 0, 0)
            spread["plain_combat"] = pool_dispersion(plain_pool)
            # The flat-gene method learns from the samples other than healthy controls of the remaining units that have enough of them, as it predicts
            # the held-out unit from its own.
            learners = [key for key in rest if visible_count[key] >= minimum]
            if len(learners) >= flat_policy["minimum_training_units"]:
                chosen = C.flat_selection(group_effects, group_genes, excluded=[held_key])
                held_flat = np.zeros(len(order), dtype=bool)
                held_flat[[position[gene] for gene in chosen.loc[chosen["selected"], "hgnc_id"]]] = True
                rest_train = [unit_features(C.combine(visible_of[key].values()), rest_mean, rest_deviation) for key in learners]
                location, scale, first, second = flat_effects(
                    (np.vstack([part[0] for part in rest_train]), np.vstack([part[1] for part in rest_train])),
                    (np.vstack([rest_effects[key][0] for key in learners]), np.log(np.vstack([rest_effects[key][1] for key in learners]))),
                    unit_features(visible, rest_mean, rest_deviation), held_flat)
                rows["flat_genes"] = (*held_out_error(C.corrected_moments(hidden, rest_mean, rest_deviation, location, scale), pool), first, second)
            rest_reference = C.class_reference(moments_of, rest, rest_mean, rest_deviation, rest_effects)
            cells = C.matched_cells(visible_cells, rest_reference, excluded=[anchor])
            matched, freedom = C.matched_samples(cells)
            if matched >= policy["minimum_matched"] and freedom >= policy["minimum_matched"] - 1:
                location, scale = C.matched_effects(cells, rest_mean, rest_deviation)
                rows["matched_groups"] = (*held_out_error(C.corrected_moments(hidden, rest_mean, rest_deviation, location, scale), pool), 0, 0)
            for method, (location_error, scale_error, count, first, second) in rows.items():
                held.append(dict(**base, method=method, genes=count, location_error=location_error, scale_error=scale_error, pool_dispersion=spread[method],
                                 sampling_floor=floor, location_components=first, scale_components=second))

        tested = pd.DataFrame([row for row in held if row["group"] == group], columns=HELD_COLUMNS)
        wide = tested.pivot_table(index=["study_accession", "platform"], columns="method", values="location_error") if len(tested) else pd.DataFrame()
        both = {name: wide.dropna(subset=[name, "matched_groups"]) if {name, "matched_groups"} <= set(wide.columns) else wide.iloc[0:0] for name in ("flat_genes", "none")}
        compared = paired(tested, policy["comparison"])
        comparisons += [dict(group=group, **row) for row in compared]
        own = [row for row in fits if row["group"] == group]
        summary[group] = dict(unanchored_units=len(own), corrected_by_matched_groups=sum(row["applicable"] for row in own),
                              median_pool_dispersion={method: round(float(tested.loc[tested["method"] == method, "pool_dispersion"].median()), 6)
                                                      for method in METHODS if (tested["method"] == method).any()},
                              median_sampling_floor=round(float(tested.drop_duplicates(["study_accession", "platform"])["sampling_floor"].median()), 6) if len(tested) else None,
                              paired_comparison=[{name: (round(value, 6) if isinstance(value, float) else value) for name, value in row.items()} for row in compared],
                              not_corrected={reason: sum(row["reason"] == reason for row in own) for reason in ("newborn_study_not_corrected", "no_class_to_line_up")},
                              held_out_units=len(wide),
                              median_location_error={method: round(float(tested.loc[tested["method"] == method, "location_error"].median()), 6)
                                                     for method in METHODS if (tested["method"] == method).any()},
                              **{f"matched_groups_closer_than_{name}": f"{int((part['matched_groups'] < part[name]).sum()) if len(part) else 0} of {len(part)}"
                                 for name, part in both.items()})
        print(f"[05] {group}: {summary[group]['corrected_by_matched_groups']} of {len(own)} unanchored units corrected by matched groups; hidden-control test "
              f"over {len(wide)} units, median location error " + ", ".join(f"{name} {value:.3f}" for name, value in summary[group]["median_location_error"].items())
              + f"; matched groups closer than flat genes in {summary[group]['matched_groups_closer_than_flat_genes']}", flush=True)
    effects = pd.concat(effects, ignore_index=True)[EFFECT_COLUMNS]
    fits = pd.DataFrame(fits, columns=FIT_COLUMNS)
    held = pd.DataFrame(held, columns=HELD_COLUMNS)
    comparisons = pd.DataFrame(comparisons, columns=COMPARISON_COLUMNS)
    C.require((effects["scale"].dropna() > 0).all() and set(fits.loc[~fits["applicable"], "reason"]) <= set(C.REASONS), "Invalid corrections")
    C.save_table(effects, C.DATA / "candidate_effects.parquet")
    C.save_table(fits, C.DATA / "candidate_fits.parquet")
    C.save_table(held, C.DATA / "held_out_test.parquet")
    C.save_table(comparisons, C.DATA / "held_out_comparison.parquet")
    C.save_review(fits, C.OUT / "candidate_fits_review.csv")
    C.save_review(held, C.OUT / "held_out_test_review.csv")
    C.save_review(comparisons, C.OUT / "held_out_comparison_review.csv")
    C.finish(5, start, [C.DATA / "candidate_effects.parquet", C.DATA / "candidate_fits.parquet", C.DATA / "held_out_test.parquet",
                        C.DATA / "held_out_comparison.parquet", C.OUT / "candidate_fits_review.csv", C.OUT / "held_out_test_review.csv",
                        C.OUT / "held_out_comparison_review.csv"], method=policy["method"], groups=summary,
             limits={**{name: policy[name] for name in ("matched_classes", "matched_age_groups", "minimum_matched", "components", "comparison")},
                     **{name: flat_policy[name] for name in ("minimum_training_units", "minimum_shared_genes", "minimum_complete_genes")}},
             protocol=C.CONFIG["protocol"])
    print(f"[05] Complete: review sheets: {(C.OUT / 'candidate_fits_review.csv').relative_to(C.ROOT).as_posix()}, "
          f"{(C.OUT / 'held_out_test_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
