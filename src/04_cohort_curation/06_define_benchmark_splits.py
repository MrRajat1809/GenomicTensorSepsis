"""Define the benchmark tasks and their split schemes, each of which assigns whole patient-safe groups to one side of every fold."""
import importlib

C = importlib.import_module("01_catalogue_sample_metadata")

KEYS = ["source", "study_accession", "sample_accession", "platform"]
REASONS = ("not_primary_record", "not_one_per_patient", "quality_flag", "tissue_outside_task", "phase_outside_task", "class_outside_task",
           "label_not_recorded")
SCHEMES = ("leave_one_group_out", "technology_transfer", "temporal")
TRANSFERS = (("array", "rna_seq"), ("rna_seq", "array"))
TEMPORAL_FOLD = "earlier_to_later"
TASK_COLUMNS = [*KEYS, "task", "patient_safe_group", "eligible", "reason", "label", "label_value"]
SPLIT_COLUMNS = ["task", "scheme", "fold", "patient_safe_group", "role", "technology", "first_release", "positives", "negatives"]


def exclusions(table, task):
    """Return the test of every reason of REASONS for one task, in that order; a record is excluded by the first test that it meets."""
    policy = C.CONFIG["benchmark"]
    labelled = table[task["variable"]].isin([task["positive"], task["negative"]])
    in_class = table["sample_class"] == task["within_class"] if "within_class" in task else labelled
    return [~table["is_primary_record"], ~table["one_per_patient"], table["any_flag"], ~table["tissue"].isin(policy["tissues"]),
            table["sepsis_phase"].isin(policy["phases_left_out"]), ~in_class, ~labelled]


def fold_roles(groups):
    """Return the scheme, the name and the role of every group for every fold; groups holds the technology, first release and class counts of the groups."""
    policy = C.CONFIG["benchmark"]
    scored = (groups["positives"] >= policy["minimum_per_class"]) & (groups["negatives"] >= policy["minimum_per_class"])
    held_out = scored.map({True: "test", False: "held_out_not_scored"})
    folds = [("leave_one_group_out", name, held_out.where(groups.index == name, "train")) for name in sorted(groups.index[scored])]
    folds += [("technology_transfer", f"{first}_to_{second}", held_out.where(groups["technology"] == second, "train")) for first, second in TRANSFERS]
    folds.append(("temporal", TEMPORAL_FOLD, held_out.where(groups["first_release"] >= policy["temporal_split_date"], "train")))
    return folds


def main():
    import pandas as pd

    print(f"[06] Defining benchmark tasks and splits v{C.VERSION}", flush=True)
    start = C.begin(6)
    policy = C.CONFIG["benchmark"]
    table = pd.read_parquet(C.DATA / "cohort_samples.parquet")
    platforms = pd.read_parquet(C.DATA / "cohort_summary.parquet", columns=["study_accession", "platform", "release_date"]).merge(
        pd.read_parquet(C.ANNOTATION / "matrix_index.parquet", columns=["accession", "platform", "family"]).rename(columns={"accession": "study_accession"}),
        on=["study_accession", "platform"], how="left", validate="one_to_one")
    unknown = sorted(set(platforms["family"].fillna("(none)")) - set(policy["technologies"]))
    C.require(not unknown, f"Technology families without a benchmark technology: {unknown}; classify them in benchmark.technologies of cohort_config.json")
    platforms["release_date"] = platforms["release_date"].astype(str)
    C.require(platforms["release_date"].str.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}").all(), "Every study platform needs its release date")
    platforms["technology"] = platforms["family"].map(policy["technologies"])
    table = table.merge(platforms[["study_accession", "platform", "technology", "release_date"]], on=["study_accession", "platform"], validate="many_to_one")
    mixed = table.groupby("patient_safe_group")["technology"].nunique()
    C.require((mixed == 1).all(), f"Patient-safe groups that hold arrays and RNA-seq: {sorted(mixed[mixed > 1].index)}")
    dated = table.groupby("patient_safe_group").agg(technology=("technology", "first"), first_release=("release_date", "min"))

    parts, splits, review, summary = [], [], [], {}
    for task in policy["tasks"]:
        name = task["task"]
        reason = pd.Series("", index=table.index)
        for code, test in zip(REASONS, exclusions(table, task)):
            reason = reason.mask((reason == "") & test, code)
        eligible = reason == ""
        value = table[task["variable"]]
        part = table[KEYS + ["patient_safe_group"]].assign(task=name, eligible=eligible, reason=reason, label=(value == task["positive"]).astype("Int8").where(eligible),
                                                           label_value=value.where(eligible, ""))[TASK_COLUMNS]
        chosen = part[eligible]
        known = table.loc[eligible & table["patient_identified"], "patient_key"]
        C.require(not known.duplicated().any() and not table.loc[eligible, "sample_id"].duplicated().any(),
                  f"Task {name} must hold one sample of every patient")
        groups = pd.DataFrame({"positives": (chosen["label"] == 1).groupby(chosen["patient_safe_group"]).sum().astype(int),
                               "negatives": (chosen["label"] == 0).groupby(chosen["patient_safe_group"]).sum().astype(int)}).join(dated, validate="one_to_one")
        folds = {}
        for scheme, fold, roles in fold_roles(groups):
            side = groups.assign(role=roles)
            train, test, unscored = (side[side["role"] == role] for role in ("train", "test", "held_out_not_scored"))
            C.require(int(train["positives"].sum()) > 0 and int(train["negatives"].sum()) > 0 and len(test) > 0,
                      f"Fold {fold} of {scheme} for {name} needs both classes among its training groups and at least one test group")
            splits.append(side.rename_axis("patient_safe_group").reset_index().assign(task=name, scheme=scheme, fold=fold)[SPLIT_COLUMNS])
            review.append(dict(task=name, scheme=scheme, fold=fold, train_groups=len(train), train_positives=int(train["positives"].sum()),
                               train_negatives=int(train["negatives"].sum()), test_groups=len(test), test_positives=int(test["positives"].sum()),
                               test_negatives=int(test["negatives"].sum()), held_out_not_scored_groups=len(unscored),
                               held_out_not_scored_samples=int((unscored["positives"] + unscored["negatives"]).sum())))
            folds[scheme] = folds.get(scheme, 0) + 1
        excluded = {code: int((reason == code).sum()) for code in REASONS}
        C.require(sum(excluded.values()) + int(eligible.sum()) == len(table), f"The records of task {name} must reconcile")
        scored = (groups["positives"] >= policy["minimum_per_class"]) & (groups["negatives"] >= policy["minimum_per_class"])
        summary[name] = dict(definition=task, records=len(part), eligible=int(eligible.sum()), positives=int((chosen["label"] == 1).sum()),
                             negatives=int((chosen["label"] == 0).sum()), patient_safe_groups=len(groups), groups_with_the_minimum_of_each_class=int(scored.sum()),
                             samples_in_those_groups=int((groups.loc[scored, "positives"] + groups.loc[scored, "negatives"]).sum()),
                             groups_by_technology={str(key): int(count) for key, count in groups["technology"].value_counts().sort_index().items()},
                             groups_released_before_the_split_date=int((groups["first_release"] < policy["temporal_split_date"]).sum()),
                             excluded=excluded, folds=folds)
        parts.append(part)

    tasks = pd.concat(parts, ignore_index=True)
    splits = pd.concat(splits, ignore_index=True)
    C.save_table(tasks, C.DATA / "benchmark_tasks.parquet")
    C.save_table(splits, C.DATA / "benchmark_splits.parquet")
    C.save_review(pd.DataFrame(review), C.OUT / "benchmark_review.csv")
    C.finish(6, start, [C.DATA / "benchmark_tasks.parquet", C.DATA / "benchmark_splits.parquet", C.OUT / "benchmark_review.csv"],
             method=policy["method"], policy={key: value for key, value in policy.items() if key != "method"}, reasons=list(REASONS), schemes=list(SCHEMES),
             tasks=summary, split_rows=len(splits), protocol=C.CONFIG["protocol"])
    print("[06] Complete: " + "; ".join(f"{name} {entry['eligible']:,} samples ({entry['positives']:,} positive) in {entry['patient_safe_groups']} groups, "
                                        f"{entry['groups_with_the_minimum_of_each_class']} of them scored" for name, entry in summary.items())
          + f"; {len(review)} folds; review sheet: {(C.OUT / 'benchmark_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
