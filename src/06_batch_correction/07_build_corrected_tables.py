"""Write the corrected table of every group with a table of its samples and their full metadata, the effects and model that rebuild it, and a ledger of every record."""
import importlib

C = importlib.import_module("01_select_correction_inputs")

EXTRA_COLUMNS = ["group", "family", "value_scale", "expression_sex", "sex_used", "female", "estimation", "unit_tier", "unit_method", "method_basis", "checks_applied"]
PARAMETER_COLUMNS = ["group", "study_accession", "platform", "role", "method", "hgnc_id", "location", "scale"]
LEDGER_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "sample_id", "is_primary_record", "tissue", "sample_class", "group",
                  "unit_method", "outcome", "reason"]


def main():
    import numpy as np
    import pandas as pd

    print(f"[07] Building corrected tables v{C.VERSION}", flush=True)
    start = C.begin(7)
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    methods = pd.read_parquet(C.DATA / "unit_methods.parquet")
    checks = pd.read_parquet(C.DATA / "unit_checks.parquet")
    cohort = pd.read_parquet(C.COHORT / "cohort_samples.parquet")
    effect_table = pd.concat([pd.read_parquet(C.DATA / "anchored_effects.parquet"), pd.read_parquet(C.DATA / "candidate_effects.parquet")], ignore_index=True)
    applied = checks[checks["applicable"]].assign(name=lambda frame: np.where(frame["sample_class"] != "", frame["check"] + ":" + frame["sample_class"], frame["check"]))
    applied = applied.groupby(["group", "study_accession", "platform", "method"])["name"].agg(lambda part: ";".join(part)).to_dict()
    outputs, parameters, summary = [], [], {}
    records["unit_method"] = ""
    records["outcome"] = np.where(records["status"] == "in_unit", "", "left_out")
    where = records.set_index(C.KEYS).index
    C.require(where.is_unique and (records.index == range(len(records))).all(), "Records must be keyed once")
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        model = C.model_arrays(group, order)
        symbols = genes[genes["group"] == group].set_index("hgnc_id")["symbol"]
        chosen = methods[(methods["group"] == group) & (methods["status"] == "corrected")].sort_values(["study_accession", "platform"], kind="stable")
        merged = chosen[chosen["role"] == "merged"]
        effects = {key: part for key, part in effect_table[effect_table["group"] == group].groupby(["study_accession", "platform", "method"])}
        sizes = {(unit.study_accession, unit.platform): int(C.unit_records(records, group, unit.study_accession, unit.platform)["is_primary_record"].sum())
                 for unit in merged.itertuples(index=False)}
        matrix = np.full((len(order), sum(sizes.values())), np.nan, dtype="float32")
        columns, rows, filled = [], [], 0
        for unit in chosen.itertuples(index=False):
            part = effects[(unit.study_accession, unit.platform, unit.method)].set_index("hgnc_id").reindex(order)
            parameters.append(pd.DataFrame({"group": group, "study_accession": unit.study_accession, "platform": unit.platform, "role": unit.role,
                                            "method": unit.method, "hgnc_id": order, "location": part["location"].to_numpy(), "scale": part["scale"].to_numpy()}))
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            keys = list(zip(own["study_accession"], own["sample_accession"], own["platform"]))
            index = where.get_indexer(keys)
            C.require((index >= 0).all(), "Every record of a unit must be in the records")
            records.loc[index, "unit_method"] = unit.method
            primary = own["is_primary_record"].to_numpy() & (unit.role == "merged")
            records.loc[index[primary], "outcome"] = "in_corrected_table"
            records.loc[index[~primary], ["outcome", "reason"]] = ["corrected_for_checks_only", "repeat_deposit"]
            if unit.role != "merged":
                continue
            own = own[own["is_primary_record"]].reset_index(drop=True)
            raw = C.read_aligned(group, unit.study_accession, unit.platform, own["sample_accession"])
            fixed = C.correct(raw, own["female"].to_numpy(), model, part["location"].to_numpy(), part["scale"].to_numpy())
            matrix[:, filled:filled + len(own)] = fixed
            filled += len(own)
            columns += list(own["sample_id"])
            rows.append(own.assign(unit_tier=unit.tier, unit_method=unit.method, method_basis=unit.basis,
                                   checks_applied=applied.get((group, unit.study_accession, unit.platform, unit.method), "")))
        C.require(filled == matrix.shape[1] and len(set(columns)) == len(columns), f"Every sample of {group} needs exactly one column")
        held = ~np.isnan(matrix).all(axis=1)
        frame = pd.DataFrame(matrix[held], columns=columns)
        frame.insert(0, "symbol", [symbols[gene] for gene, keep in zip(order, held) if keep])
        frame.insert(0, "hgnc_id", [gene for gene, keep in zip(order, held) if keep])
        samples = pd.concat(rows, ignore_index=True)
        samples = samples[["source", "study_accession", "sample_accession", "platform", *EXTRA_COLUMNS]].merge(
            cohort, on=["source", "study_accession", "sample_accession", "platform"], how="left", validate="one_to_one")
        samples = samples[[*cohort.columns, *EXTRA_COLUMNS]]
        C.require((samples["sample_id"].to_numpy() == np.array(columns)).all(), f"The sample table of {group} must follow the columns")
        expression, sample_path = C.CORRECTED / f"{group}_expression.parquet", C.CORRECTED / f"{group}_samples.parquet"
        C.save_table(frame, expression)
        C.save_table(samples, sample_path)
        outputs += [expression, sample_path]
        summary[group] = dict(samples=len(columns), genes=int(held.sum()), genes_without_any_value=int((~held).sum()), studies=int(merged["study_accession"].nunique()),
                              units=len(merged), by_method={str(name): int(count) for name, count in samples["unit_method"].value_counts().sort_index().items()},
                              by_tissue={str(name): int(count) for name, count in samples["tissue"].value_counts().sort_index().items()},
                              by_class={str(name): int(count) for name, count in samples["sample_class"].value_counts().sort_index().items()},
                              empty_share=round(float(np.isnan(matrix[held]).mean()), 6))
        print(f"[07] {group}: {len(columns):,} samples of {len(merged)} study units over {int(held.sum()):,} genes "
              f"({summary[group]['empty_share']:.1%} of cells empty: genes a study does not measure)", flush=True)
        del matrix, frame
    left = records["status"] == "in_unit"
    unchosen = left & (records["outcome"] == "")
    because = methods[methods["status"] == "left_out"].set_index(C.UNIT)["reason"]
    records.loc[unchosen, "outcome"] = "left_out"
    records.loc[unchosen, "reason"] = [because[key] for key in map(tuple, records.loc[unchosen, C.UNIT].to_numpy())]
    C.require((records["outcome"] != "").all() and set(records["reason"]) <= set(C.REASONS) | {""}, "Every record needs an outcome")
    ledger = records[LEDGER_COLUMNS]
    parameters = pd.concat(parameters, ignore_index=True)[PARAMETER_COLUMNS]
    C.save_table(ledger, C.DATA / "correction_ledger.parquet")
    C.save_table(parameters, C.DATA / "correction_parameters.parquet")
    outcomes = {str(name): int(count) for name, count in ledger["outcome"].value_counts().sort_index().items()}
    reasons = {str(name): int(count) for name, count in ledger.loc[ledger["outcome"] == "left_out", "reason"].value_counts().sort_index().items()}
    C.finish(7, start, [*outputs, C.DATA / "correction_ledger.parquet", C.DATA / "correction_parameters.parquet"], groups=summary, records=dict(
        total=len(ledger), by_outcome=outcomes, left_out_by_reason=reasons), reasons=C.REASONS, protocol=C.CONFIG["protocol"])
    print("[07] Complete: " + "; ".join(f"{name} {count:,}" for name, count in outcomes.items()), flush=True)


if __name__ == "__main__":
    main()
