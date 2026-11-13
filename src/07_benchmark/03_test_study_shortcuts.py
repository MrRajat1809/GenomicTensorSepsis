"""Measure how far the study alone predicts the class of a record, and how far expression identifies the study before and after the correction."""
import importlib
import warnings

C = importlib.import_module("01_audit_split_leakage")

FOLD_COLUMNS = ["task", "scheme", "scoring", "fold", "unit", "train_records", "records", "positives", "negatives", "unseen_study_records", "status", "reason",
                "auroc", "auroc_low", "auroc_high", "auprc", "auprc_low", "auprc_high", "prevalence", "pooled_auroc"]
PREDICTION_COLUMNS = ["task", "scheme", "fold", *C.KEYS, "label", "score"]
IDENTITY_COLUMNS = ["table", "fold", "train_samples", "test_samples", "units", "genes", "balanced_accuracy", "accuracy", "chance"]
IDENTITY_PREDICTION_COLUMNS = ["table", "fold", *C.KEYS, "sample_id", "unit", "predicted"]
TABLES = ("aligned", "corrected")
METRICS = ("auroc", "auroc_low", "auroc_high", "auprc", "auprc_low", "auprc_high", "prevalence")
UNIT_REASONS = {"fewer_than_minimum_of_a_class": "The held-out study holds fewer than splits.minimum_per_class records of one of the two classes."}


def study_shares(train_studies, train_labels, test_studies):
    """The share of the positive class among the training records of the study of every test record, or among all training records for an unseen study."""
    import pandas as pd

    shares = pd.Series(train_labels, dtype="float64").groupby(train_studies).mean()
    found = pd.Series(test_studies).map(shares)
    return found.fillna(float(train_labels.mean())).to_numpy(dtype="float64"), int(found.isna().sum())


def shortcut_rows(frame, task):
    """The study-only predictor of one task: over the pooled test records of every naive random split, and within every held-out study of the naive
    leave-one-study-out split that holds the minimum of each class; the AUROC over all held-out records pooled is kept beside the latter."""
    import numpy as np
    from sklearn.metrics import roc_auc_score

    labels, studies = frame["label"].to_numpy(dtype="int64"), frame["study_accession"].to_numpy()
    rows, predictions = [], []

    def row(scheme, scoring, fold, unit, test, missing, status, reason, metrics):
        return dict(task=task, scheme=scheme, scoring=scoring, fold=fold, unit=unit, train_records=int((~test).sum()), records=int(test.sum()),
                    positives=int(labels[test].sum()), negatives=int((labels[test] == 0).sum()), unseen_study_records=missing, status=status, reason=reason,
                    **(metrics or {name: np.nan for name in METRICS}), pooled_auroc=np.nan)

    for repeat in range(1, C.CONFIG["splits"]["repeats"] + 1):
        test, fold = C.record_split(frame, repeat), C.repeat_name(repeat)
        score, missing = study_shares(studies[~test], labels[~test], studies[test])
        C.require(C.scorable(labels[test]), f"The test records of {fold} of the random split for {task} need the minimum of each class")
        found = row("naive_record_split", "pooled", fold, "pooled", test, missing, "scored", "",
                    C.score_unit(labels[test], score, C.unit_seed("study_only", task, "naive_record_split", fold)))
        rows.append(dict(found, pooled_auroc=found["auroc"]))
        predictions.append(frame.loc[test, C.KEYS].assign(task=task, scheme="naive_record_split", fold=fold, label=labels[test], score=score))
    # A held-out study is unseen in training, so all of its records get one score and its own AUROC is 0.5 by construction.
    held, start = np.full(len(frame), np.nan), len(rows)
    for study in sorted(set(studies)):
        test = studies == study
        score, missing = study_shares(studies[~test], labels[~test], studies[test])
        held[test] = score
        ready = C.scorable(labels[test])
        rows.append(row("naive_leave_one_study_out", "within_study", study, study, test, missing, "scored" if ready else "not_scored",
                        "" if ready else "fewer_than_minimum_of_a_class",
                        C.score_unit(labels[test], score, C.unit_seed("study_only", task, "naive_leave_one_study_out", study)) if ready else None))
        predictions.append(frame.loc[test, C.KEYS].assign(task=task, scheme="naive_leave_one_study_out", fold=study, label=labels[test], score=score))
    pooled = float(roc_auc_score(labels, held))
    for found in rows[start:]:
        found["pooled_auroc"] = pooled
    return rows, predictions


def identity_inputs():
    """The samples of the corrected table of the study-identification group whose study platform holds enough of them, and their aligned and corrected
    values over the genes with a value in every one of them in both tables: samples in rows, as float32."""
    import numpy as np
    import pandas as pd

    policy = C.CONFIG["shortcuts"]["study_identification"]
    group = policy["group"]
    samples = pd.read_parquet(C.CORRECTION / "corrected" / f"{group}_samples.parquet", columns=[*C.KEYS, "sample_id"])
    samples["unit"] = samples["study_accession"] + "__" + samples["platform"]
    sizes = samples.groupby("unit")["sample_id"].transform("size")
    samples = samples[sizes >= policy["minimum_unit_samples"]].reset_index(drop=True)
    frame = pd.read_parquet(C.CORRECTION / "corrected" / f"{group}_expression.parquet", columns=["hgnc_id", *samples["sample_id"]])
    corrected = frame[list(samples["sample_id"])].to_numpy(dtype="float32").T
    genes = list(frame["hgnc_id"])
    del frame
    aligned = np.full(corrected.shape, np.nan, dtype="float32")
    units = samples.groupby(["study_accession", "platform"]).indices
    for study, platform in sorted(units):
        positions = units[(study, platform)]
        part = pd.read_parquet(C.CORRECTION / "aligned" / group / f"{study}__{platform}.parquet", columns=["hgnc_id", *samples.loc[positions, "sample_accession"]])
        part = part.set_index("hgnc_id").reindex(genes)
        aligned[positions] = part[list(samples.loc[positions, "sample_accession"])].to_numpy(dtype="float32").T
    complete = np.isfinite(corrected).all(axis=0) & np.isfinite(aligned).all(axis=0)
    C.require(int(complete.sum()) >= policy["genes"], "Too few genes with values in every sample of the study-identification set")
    return samples, {"aligned": aligned[:, complete], "corrected": corrected[:, complete]}, [gene for gene, keep in zip(genes, complete) if keep]


def identify_studies(values, units, train, test):
    """Fit multinomial logistic regression of the study platform on the most variable genes of the training samples and predict the test samples."""
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    policy = C.CONFIG["shortcuts"]["study_identification"]
    spread = values[train].astype("float64").var(axis=0)
    columns = np.sort(np.argsort(-spread, kind="stable")[:policy["genes"]])
    model = make_pipeline(StandardScaler(), LogisticRegression(C=policy["C"], max_iter=policy["max_iter"], tol=policy["tol"]))
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        try:
            model.fit(values[np.ix_(train, columns)], units[train])
        except ConvergenceWarning as error:
            raise ValueError(f"Study identification did not converge ({error}); raise shortcuts.study_identification.max_iter") from error
    return model.predict(values[np.ix_(test, columns)]), len(columns)


def main():
    import numpy as np
    import pandas as pd
    from sklearn.metrics import accuracy_score, balanced_accuracy_score
    from sklearn.model_selection import StratifiedKFold

    print(f"[03] Testing study shortcuts v{C.VERSION}", flush=True)
    start = C.begin(3)
    policy = C.cohort_policy()
    populations = pd.read_parquet(C.DATA / "split_populations.parquet")
    rows, predictions = [], []
    for task in [entry["task"] for entry in policy["tasks"]]:
        frame = populations[(populations["population"] == "naive") & (populations["task"] == task)].reset_index(drop=True)
        found, made = shortcut_rows(frame, task)
        rows += found
        predictions += made
    folds = pd.DataFrame(rows, columns=FOLD_COLUMNS)
    predictions = pd.concat(predictions, ignore_index=True)[PREDICTION_COLUMNS]
    shortcut = {}
    for task, by_task in folds.groupby("task", sort=False):
        split, held = by_task[by_task["scheme"] == "naive_record_split"], by_task[(by_task["scheme"] == "naive_leave_one_study_out") & (by_task["status"] == "scored")]
        shortcut[task] = dict(random_split=dict(repeats=len(split), pooled_auroc_mean=C.rounded(split["auroc"].mean()), auprc_mean=C.rounded(split["auprc"].mean()),
                                                prevalence_mean=C.rounded(split["prevalence"].mean())),
                              held_out_study=dict(studies=int((by_task["scheme"] == "naive_leave_one_study_out").sum()), studies_scored=len(held),
                                                  within_study_auroc_mean=C.rounded(held["auroc"].mean()), within_study_auroc_lowest=C.rounded(held["auroc"].min()),
                                                  within_study_auroc_highest=C.rounded(held["auroc"].max())))
    print("[03] Study only: " + "; ".join(f"{task} {entry['random_split']['pooled_auroc_mean']:.3f} under the random split, "
                                          f"{entry['held_out_study']['within_study_auroc_mean']:.3f} within {entry['held_out_study']['studies_scored']} held-out studies"
                                          for task, entry in shortcut.items()), flush=True)

    identify = C.CONFIG["shortcuts"]["study_identification"]
    samples, values, genes = identity_inputs()
    units = samples["unit"].to_numpy()
    splitter = StratifiedKFold(n_splits=identify["folds"], shuffle=True, random_state=identify["seed"])
    identity, guessed = [], []
    for number, (train, test) in enumerate(splitter.split(np.zeros(len(samples)), units), start=1):
        fold = f"fold_{number:02d}"
        for table in TABLES:
            predicted, used = identify_studies(values[table], units, train, test)
            identity.append(dict(table=table, fold=fold, train_samples=len(train), test_samples=len(test), units=int(len(set(units))), genes=used,
                                 balanced_accuracy=float(balanced_accuracy_score(units[test], predicted)), accuracy=float(accuracy_score(units[test], predicted)),
                                 chance=1 / len(set(units))))
            guessed.append(samples.loc[test, [*C.KEYS, "sample_id", "unit"]].assign(table=table, fold=fold, predicted=predicted))
        print(f"[03] Study identification {fold}: " + "; ".join(f"{row['table']} {row['balanced_accuracy']:.3f}" for row in identity[-len(TABLES):]), flush=True)
    guessed = pd.concat(guessed, ignore_index=True)[IDENTITY_PREDICTION_COLUMNS]
    for table in TABLES:
        part = guessed[guessed["table"] == table]
        identity.append(dict(table=table, fold="all", train_samples=0, test_samples=len(part), units=int(part["unit"].nunique()), genes=identify["genes"],
                             balanced_accuracy=float(balanced_accuracy_score(part["unit"], part["predicted"])),
                             accuracy=float(accuracy_score(part["unit"], part["predicted"])), chance=1 / part["unit"].nunique()))
    identity = pd.DataFrame(identity, columns=IDENTITY_COLUMNS)

    folds = folds.sort_values(["task", "scheme", "fold"], kind="stable").reset_index(drop=True)
    predictions = predictions.sort_values(["task", "scheme", "fold", "study_accession", "sample_accession"], kind="stable").reset_index(drop=True)
    identity = identity.sort_values(["table", "fold"], kind="stable").reset_index(drop=True)
    guessed = guessed.sort_values(["table", "fold", "study_accession", "sample_accession"], kind="stable").reset_index(drop=True)
    outputs = {name: C.DATA / f"{name}.parquet" for name in ("shortcut_folds", "shortcut_predictions", "study_identification", "study_identification_predictions")}
    for name, frame in zip(outputs, (folds, predictions, identity, guessed)):
        C.save_table(frame, outputs[name])
    C.save_review(folds, C.OUT / "shortcut_folds_review.csv")
    C.save_review(identity, C.OUT / "study_identification_review.csv")
    overall = identity[identity["fold"] == "all"].set_index("table")
    C.finish(3, start, [*outputs.values(), C.OUT / "shortcut_folds_review.csv", C.OUT / "study_identification_review.csv"], method=C.CONFIG["shortcuts"]["method"],
             study_only=shortcut, unit_reasons=UNIT_REASONS, study_identification=dict(samples=len(samples), units=int(len(set(units))), genes_complete=len(genes),
                                                            genes_used=identify["genes"], folds=identify["folds"],
                                                            balanced_accuracy={table: C.rounded(overall.at[table, "balanced_accuracy"]) for table in TABLES},
                                                            chance=C.rounded(1 / len(set(units)))),
             policy=identify, protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: study platform identified from expression with balanced accuracy {overall.at['aligned', 'balanced_accuracy']:.3f} before the "
          f"correction and {overall.at['corrected', 'balanced_accuracy']:.3f} after it (chance {1 / len(set(units)):.3f}, {len(set(units))} study platforms); "
          f"review sheets: {(C.OUT / 'shortcut_folds_review.csv').relative_to(C.ROOT).as_posix()}, "
          f"{(C.OUT / 'study_identification_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
