"""Check every corrected unit for study mixing, sex genes and the Sepsis MetaScore, keep those that pass, and compare the two deposits of a sample."""
import importlib

C = importlib.import_module("01_select_correction_inputs")

CHECK_COLUMNS = ["group", "study_accession", "platform", "role", "tier", "method", "check", "sample_class", "gate", "applicable", "samples", "before", "after",
                 "limit", "passed"]
METHOD_COLUMNS = ["group", "study_accession", "platform", "role", "tier", "family", "neonatal", "method_tried", "method", "basis", "status", "reason"]
REPORTED_COLUMNS = ["locations_beyond_limit", "scales_beyond_limit"]
REPEAT_COLUMNS = ["group", "sample_id", "study_a", "sample_a", "platform_a", "study_b", "sample_b", "platform_b", "method_a", "method_b", "shared_genes",
                  "distance_before", "distance_after", "rank_before", "rank_after", "candidates"]
POOLED_COLUMNS = ["group", "scope", "check", "samples", "before", "after"]


class Values:
    """The aligned values of the estimation samples of every unit on the genes that the checks read, before and after a given correction."""

    def __init__(self, group, units, records, order, wanted, model):
        import numpy as np

        self.index = np.array(wanted)
        self.model = {name: values[self.index] for name, values in model.items()}
        self.samples, self.raw = {}, {}
        for unit in units[(units["group"] == group) & (units["status"] == "in_unit")].itertuples(index=False):
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            own = own[own["estimation"]].reset_index(drop=True)
            values = C.read_aligned(group, unit.study_accession, unit.platform, own["sample_accession"])
            self.samples[(unit.study_accession, unit.platform)] = own
            self.raw[(unit.study_accession, unit.platform)] = values[self.index]

    def corrected(self, key, effect):
        location, scale = effect
        return C.correct(self.raw[key], self.samples[key]["female"].to_numpy(), self.model, location[self.index], scale[self.index])


def projection(pool, components):
    """The mean and principal axes of the pooled samples; missing values count as the mean."""
    import numpy as np

    with np.errstate(all="ignore"):
        centre = np.nanmean(pool, axis=1)
    centred = np.nan_to_num(pool - centre[:, None])
    _, _, axes = np.linalg.svd(centred.T, full_matrices=False)
    return lambda values: np.nan_to_num(values - centre[:, None]).T @ axes[:components].T


def neighbour_excess(points, own, neighbours):
    """The mean share of the nearest neighbours of a unit's samples that come from the unit, less the share expected by chance."""
    import numpy as np

    targets = points[own]
    distances = (targets ** 2).sum(axis=1)[:, None] + (points ** 2).sum(axis=1)[None, :] - 2 * targets @ points.T
    distances[np.arange(len(targets)), np.flatnonzero(own)] = np.inf
    nearest = np.argpartition(distances, neighbours, axis=1)[:, :neighbours]
    return float(own[nearest].mean() - (own.sum() - 1) / (len(own) - 1))


def sex_difference(values, female, male):
    """The mean of the Y-chromosome genes less XIST for every sample, or the mean of the Y-chromosome genes where XIST is not read."""
    import numpy as np

    with np.errstate(all="ignore"):
        result = np.nanmean(values[male], axis=0)
    return result - values[female[0]] if female else result


def metascore(values, up, down):
    import numpy as np

    with np.errstate(all="ignore"):
        return np.nanmean(values[up], axis=0) - np.nanmean(values[down], axis=0)


def main():
    import numpy as np
    import pandas as pd

    print(f"[06] Evaluating corrections v{C.VERSION}", flush=True)
    start = C.begin(6)
    policy = C.CONFIG["checks"]
    anchor = C.CONFIG["anchored"]["anchor_class"]
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    named = pd.read_parquet(C.VALIDATION / "validation_genes.parquet")
    pairs = pd.read_parquet(C.VALIDATION / "repeat_pairs.parquet")
    effect_table = pd.concat([pd.read_parquet(C.DATA / "anchored_effects.parquet"), pd.read_parquet(C.DATA / "candidate_effects.parquet")], ignore_index=True)
    fits = pd.read_parquet(C.DATA / "candidate_fits.parquet")
    checks, methods, repeats, pooled, summary = [], [], [], [], {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        position = {gene: index for index, gene in enumerate(order)}
        model = C.model_arrays(group, order)
        active = units[(units["group"] == group) & (units["status"] == "in_unit")].set_index(["study_accession", "platform"])
        pool_keys = [key for key, unit in active.iterrows() if unit["tier"] == "anchored" and unit["role"] == "merged"]
        effects = {}
        for (study, platform, method), part in effect_table[effect_table["group"] == group].groupby(["study_accession", "platform", "method"]):
            part = part.set_index("hgnc_id").reindex(order)
            effects[((study, platform), method)] = (part["location"].to_numpy(), part["scale"].to_numpy())

        # The genes that the checks read: the most variable genes of the corrected pool, the sex genes and the score genes.
        moments_of = C.unit_moments(group, units, records, model["sex_effect"])
        complete = np.all([C.gather(moments_of[key])[0] > 0 for key in pool_keys], axis=0)
        after = C.combine([(count, centre, spread * (count - 1)) for count, centre, spread in
                           (C.corrected_moments(C.gather(moments_of[key]), model["anchor_mean"], model["anchor_sd"], *effects[(key, "healthy_anchored")])
                            for key in pool_keys)])
        spread_after = np.where(complete, C.variance(after) * model["anchor_sd"] ** 2, -np.inf)
        top_after = list(np.argsort(-np.nan_to_num(spread_after, nan=-np.inf))[:policy["pca_genes"]])
        role_genes = {role: [position[gene] for gene in part["hgnc_id"] if gene in position] for role, part in named.groupby(["test", "role"])}
        female, male = role_genes.get(("sex", "female"), []), role_genes.get(("sex", "male"), [])
        up, down = role_genes.get(("sepsis_score", "higher_in_sepsis"), []), role_genes.get(("sepsis_score", "lower_in_sepsis"), [])
        wanted = sorted(set(top_after) | set(female) | set(male) | set(up) | set(down))
        local = {gene: index for index, gene in enumerate(wanted)}
        values = Values(group, units, records, order, wanted, model)
        pick = lambda indices: [local[index] for index in indices]
        pool_after = {key: values.corrected(key, effects[(key, "healthy_anchored")]) for key in pool_keys}
        # One space for before and after: the principal components of the corrected pool, into which the aligned values are projected as well.
        used = pick(top_after)
        project = projection(np.hstack([pool_after[key][used] for key in pool_keys]), policy["components"])
        pool_values = {"after": pool_after, "before": {key: values.raw[key] for key in pool_keys}}
        projected = {stage: {key: project(pool_values[stage][key][used]) for key in pool_keys} for stage in ("before", "after")}
        classes_of = {key: own["sample_class"].to_numpy() for key, own in values.samples.items()}
        splits = {}

        def midpoint(stage, female_read, male_read):
            """The one split of the sex difference over the whole pool, on the sex genes that the unit reads."""
            if (stage, tuple(female_read), tuple(male_read)) not in splits:
                difference = np.concatenate([sex_difference(pool_values[stage][key], female_read, male_read) for key in pool_keys])
                splits[(stage, tuple(female_read), tuple(male_read))] = C.two_means(difference[np.isfinite(difference)])
            return splits[(stage, tuple(female_read), tuple(male_read))]

        def evaluate(key, method, corrected):
            """Every applicable check of one unit under one correction: before (aligned values) and after."""
            unit = active.loc[key]
            own = values.samples[key]
            rows, passed = [], True
            base = dict(group=group, study_accession=key[0], platform=key[1], role=unit["role"], tier=unit["tier"], method=method)
            mixing_gate = unit["tier"] != "anchored"
            own_values = {"before": values.raw[key], "after": corrected}
            members = [other for other in pool_keys if other != key]
            for name in policy["classes"]:
                mine = classes_of[key] == name
                others = [other for other in members if (classes_of[other] == name).sum() >= policy["minimum_class_samples"]]
                applicable = mine.sum() >= policy["minimum_class_samples"] and len(others) >= 2
                row = dict(**base, check="study_mixing", sample_class=name, gate=mixing_gate, applicable=bool(applicable), samples=int(mine.sum()), before=np.nan,
                           after=np.nan, limit=policy["mixing_worsening_limit"], passed=True)
                if applicable:
                    found = {}
                    for stage in ("before", "after"):
                        points = [project(own_values[stage][used][:, mine])] + [projected[stage][other][classes_of[other] == name] for other in members]
                        flags = np.concatenate([np.ones(int(mine.sum()), dtype=bool)] + [np.zeros(len(part), dtype=bool) for part in points[1:]])
                        found[stage] = neighbour_excess(np.vstack(points), flags, policy["neighbours"])
                    row.update(before=found["before"], after=found["after"], passed=found["after"] <= found["before"] + policy["mixing_worsening_limit"])
                rows.append(row)
            # The sex genes that the unit reads; the pool is split on the same genes.
            female_read = [gene for gene in pick(female) if np.isfinite(values.raw[key][gene]).any()]
            male_read = [gene for gene in pick(male) if np.isfinite(values.raw[key][gene]).any()]
            called = own["expression_sex"].isin(["female", "male"]).to_numpy()
            row = dict(**base, check="sex_genes", sample_class="", gate=True, applicable=bool(called.sum() >= policy["minimum_sex_calls"] and male_read),
                       samples=int(called.sum()), before=np.nan, after=np.nan, limit=policy["sex_agreement"], passed=True)
            if row["applicable"]:
                truth = (own.loc[called, "expression_sex"] == "male").to_numpy()
                found = {}
                for stage in ("before", "after"):
                    difference = sex_difference(own_values[stage], female_read, male_read)[called]
                    found[stage] = float(np.mean((difference > midpoint(stage, female_read, male_read)) == truth))
                row.update(before=found["before"], after=found["after"], passed=found["after"] >= policy["sex_agreement"])
            rows.append(row)
            sepsis = classes_of[key] == "sepsis"
            up_read = [gene for gene in pick(up) if np.isfinite(values.raw[key][gene]).any()]
            down_read = [gene for gene in pick(down) if np.isfinite(values.raw[key][gene]).any()]
            applicable = (sepsis.sum() >= policy["minimum_class_samples"] and len(up_read) >= policy["minimum_score_genes"]
                          and len(down_read) >= policy["minimum_score_genes"])
            row = dict(**base, check="sepsis_score", sample_class="sepsis", gate=True, applicable=bool(applicable), samples=int(sepsis.sum()), before=np.nan,
                       after=np.nan, limit=policy["score_auroc"], passed=True)
            if applicable:
                found = {}
                for stage in ("before", "after"):
                    healthy = np.concatenate([metascore(pool_values[stage][other], up_read, down_read)[classes_of[other] == anchor] for other in members])
                    cases = metascore(own_values[stage], up_read, down_read)[sepsis]
                    found[stage] = C.auroc(cases[np.isfinite(cases)], healthy[np.isfinite(healthy)])
                row.update(before=found["before"], after=found["after"], passed=found["after"] >= policy["score_auroc"])
            rows.append(row)
            for row in rows:
                passed = passed and (bool(row["passed"]) or not row["gate"])
            return rows, passed

        chosen, beyond = {}, {}
        for key, unit in active.iterrows():
            if unit["tier"] == "anchored":
                method, reason = "healthy_anchored", ""
            else:
                fit = fits[(fits["group"] == group) & (fits["study_accession"] == key[0]) & (fits["platform"] == key[1])].iloc[0]
                method, reason = ("matched_groups", "") if fit["applicable"] else ("", fit["reason"])
            if method:
                rows, passed = evaluate(key, method, pool_after[key] if key in pool_after else values.corrected(key, effects[(key, method)]))
                checks.extend(rows)
                # Without healthy controls only the study-mixing check tests the correction; a unit where it applies to no class is not kept.
                tested = unit["tier"] == "anchored" or any(row["check"] == "study_mixing" and row["applicable"] for row in rows)
                reason = "" if passed and tested else ("no_method_passed_the_checks" if not passed else "no_check_applicable")
                location, scale = effects[(key, method)]
                beyond[key] = (int((np.abs(location) > policy["location_report_limit"]).sum()), int((np.abs(np.log(scale)) > np.log(policy["scale_report_limit"])).sum()))
            if method and not reason:
                chosen[key] = method
            methods.append(dict(group=group, study_accession=key[0], platform=key[1], role=unit["role"], tier=unit["tier"], family=unit["family"],
                                neonatal=bool(unit["neonatal"]), method_tried=method, method=method if key in chosen else "",
                                basis={"healthy_anchored": "healthy_controls", "matched_groups": "class_references"}.get(method, "") if key in chosen else "",
                                status="corrected" if key in chosen else "left_out", reason=reason,
                                locations_beyond_limit=beyond.get(key, (0, 0))[0], scales_beyond_limit=beyond.get(key, (0, 0))[1]))

        # The final test: two deposits of one sample, held by two corrected units, before and after; never used for choosing.
        lookup = records[(records["group"] == group) & (records["status"] == "in_unit")].set_index(C.KEYS)
        cache = {}

        def unit_values(key):
            if key not in cache:
                if len(cache) >= 3:
                    cache.pop(next(iter(cache)))
                own = C.unit_records(records, group, *key)
                raw = C.read_aligned(group, key[0], key[1], own["sample_accession"])
                cache[key] = (own, raw, C.correct(raw, own["female"].to_numpy(), model, *effects[(key, chosen[key])]))
            return cache[key]

        for pair in pairs.sort_values(["study_a", "platform_a", "study_b", "platform_b"], kind="stable").itertuples(index=False):
            first, second = (pair.study_a, pair.sample_a, pair.platform_a), (pair.study_b, pair.sample_b, pair.platform_b)
            if first not in lookup.index or second not in lookup.index:
                continue
            first_key, second_key = (pair.study_a, pair.platform_a), (pair.study_b, pair.platform_b)
            if first_key not in chosen or second_key not in chosen:
                continue
            row = dict(group=group, sample_id=pair.sample_id, study_a=pair.study_a, sample_a=pair.sample_a, platform_a=pair.platform_a, study_b=pair.study_b,
                       sample_b=pair.sample_b, platform_b=pair.platform_b, method_a=chosen[first_key], method_b=chosen[second_key])
            own_a, raw_a, fixed_a = unit_values(first_key)
            own_b, raw_b, fixed_b = unit_values(second_key)
            column_a, column_b = int(np.flatnonzero(own_a["sample_accession"] == pair.sample_a)[0]), int(np.flatnonzero(own_b["sample_accession"] == pair.sample_b)[0])
            for stage, values_a, values_b in (("before", raw_a, raw_b), ("after", fixed_a, fixed_b)):
                shared = np.isfinite(values_a[:, column_a]) & np.isfinite(values_b).all(axis=1)
                distances = np.sqrt(np.mean((values_b[shared] - values_a[shared, column_a][:, None]) ** 2, axis=0))
                row[f"distance_{stage}"] = float(distances[column_b])
                row[f"rank_{stage}"] = int((distances < distances[column_b]).sum() + 1)
                row["candidates"] = len(distances)
                row["shared_genes"] = int(shared.sum())
            repeats.append(row)

        # The whole group: the Sepsis MetaScore of sepsis against healthy controls, and the one split of the sex genes, over every corrected unit, and over the
        # units corrected on their healthy controls only, whose sepsis samples took no part in their correction.
        merged = [key for key in chosen if active.loc[key, "role"] == "merged"]
        for scope, keys in (("all_units", merged), ("units_with_healthy_controls", [key for key in merged if chosen[key] == "healthy_anchored"])):
            found = {}
            for stage in ("before", "after"):
                fixed = {key: values.raw[key] if stage == "before" else (pool_after[key] if key in pool_after else values.corrected(key, effects[(key, chosen[key])]))
                         for key in keys}
                scores = {name: np.concatenate([metascore(fixed[key], pick(up), pick(down))[classes_of[key] == name] for key in keys]) for name in ("sepsis", anchor)}
                cases, others = scores["sepsis"][np.isfinite(scores["sepsis"])], scores[anchor][np.isfinite(scores[anchor])]
                found[("sepsis_score_auroc", stage)] = (len(cases) + len(others), C.auroc(cases, others))
                difference = np.concatenate([sex_difference(fixed[key], pick(female), pick(male)) for key in keys])
                truth = np.concatenate([values.samples[key]["expression_sex"].to_numpy() for key in keys])
                called = np.isin(truth, ["female", "male"]) & np.isfinite(difference)
                split = C.two_means(difference[called])
                found[("sex_agreement_one_split", stage)] = (int(called.sum()), float(np.mean((difference[called] > split) == (truth[called] == "male"))))
            for name in ("sepsis_score_auroc", "sex_agreement_one_split"):
                pooled.append(dict(group=group, scope=scope, check=name, samples=found[(name, "after")][0], before=found[(name, "before")][1],
                                   after=found[(name, "after")][1]))
        own_methods = [row for row in methods if row["group"] == group]
        # The genes of every unit whose effects lie beyond checks.location_report_limit and checks.scale_report_limit are reported and decide nothing.
        summary[group] = dict(units=len(own_methods), corrected=sum(row["status"] == "corrected" for row in own_methods),
                              left_out=sorted(f"{row['study_accession']} {row['platform']}" for row in own_methods if row["status"] == "left_out"),
                              left_out_by_reason={reason: sum(row["reason"] == reason for row in own_methods)
                                                  for reason in sorted({row["reason"] for row in own_methods if row["status"] == "left_out"})},
                              by_method={name: sum(row["method"] == name for row in own_methods) for name in ("healthy_anchored", "matched_groups")},
                              repeat_pairs=sum(row["group"] == group for row in repeats),
                              beyond_report_limits={f"{study} {platform}": dict(locations=counts[0], scales=counts[1])
                                                    for (study, platform), counts in sorted(beyond.items())})
        print(f"[06] {group}: {summary[group]['corrected']} of {len(own_methods)} units corrected ({summary[group]['by_method']}); "
              f"{len(summary[group]['left_out'])} left out; {summary[group]['repeat_pairs']} pairs of deposits compared", flush=True)
    checks = pd.DataFrame(checks, columns=CHECK_COLUMNS)
    methods = pd.DataFrame(methods, columns=[*METHOD_COLUMNS, *REPORTED_COLUMNS])
    repeats = pd.DataFrame(repeats, columns=REPEAT_COLUMNS)
    pooled = pd.DataFrame(pooled, columns=POOLED_COLUMNS)
    C.save_table(checks, C.DATA / "unit_checks.parquet")
    C.save_table(methods[METHOD_COLUMNS], C.DATA / "unit_methods.parquet")
    C.save_table(repeats, C.DATA / "repeat_test.parquet")
    C.save_table(pooled, C.DATA / "pooled_checks.parquet")
    C.save_review(methods, C.OUT / "unit_methods_review.csv")
    C.save_review(checks[checks["applicable"]], C.OUT / "unit_checks_review.csv")
    C.save_review(repeats, C.OUT / "repeat_test_review.csv")
    # Both deposits hold the same sample, so the distance after the correction is the noise that estimating two corrections adds.
    nearest = int((repeats["rank_after"] == 1).sum())
    C.require(len(repeats) > 0, "No two deposits of one sample were compared")
    C.finish(6, start, [C.DATA / "unit_checks.parquet", C.DATA / "unit_methods.parquet", C.DATA / "repeat_test.parquet", C.DATA / "pooled_checks.parquet",
                        C.OUT / "unit_methods_review.csv", C.OUT / "unit_checks_review.csv", C.OUT / "repeat_test_review.csv"], method=policy["method"],
             groups=summary, repeats=dict(pairs=len(repeats), median_distance_before=round(float(repeats["distance_before"].median()), 6),
                                          median_distance_after=round(float(repeats["distance_after"].median()), 6), nearest_after=nearest,
                                          nearest_before=int((repeats["rank_before"] == 1).sum())),
             pooled=[{name: (value if isinstance(value, str) else round(float(value), 6)) for name, value in row.items()} for row in pooled.to_dict(orient="records")],
             limits={name: policy[name] for name in policy if name != "method"}, protocol=C.CONFIG["protocol"])
    print(f"[06] Complete: {len(repeats)} pairs of deposits of one sample; each stays nearest its twin in {nearest} after the correction (median distance "
          f"{repeats['distance_before'].median():.3f} before, {repeats['distance_after'].median():.3f} after); review sheets: "
          f"{(C.OUT / 'unit_methods_review.csv').relative_to(C.ROOT).as_posix()}, {(C.OUT / 'repeat_test_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
