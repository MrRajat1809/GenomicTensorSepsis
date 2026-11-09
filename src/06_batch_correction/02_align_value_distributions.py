"""Give every sample of a group the same distribution of values, so that arrays and RNA sequencing become comparable; the order of the genes in a sample is kept."""
import importlib

C = importlib.import_module("01_select_correction_inputs")


def main():
    import numpy as np
    import pandas as pd

    print(f"[02] Aligning value distributions v{C.VERSION}", flush=True)
    start = C.begin(2)
    records = pd.read_parquet(C.DATA / "correction_records.parquet")
    units = pd.read_parquet(C.DATA / "correction_units.parquet")
    genes = pd.read_parquet(C.DATA / "correction_genes.parquet")
    points = C.CONFIG["alignment"]["grid_points"]
    grid = (np.arange(points) + 0.5) / points
    outputs, references, summary = [], [], {}
    for group in C.CONFIG["groups"]:
        order = C.group_genes(genes, group)
        active = units[(units["group"] == group) & (units["status"] == "in_unit")]
        # The reference weighs every unit that holds primary records the same, whatever its size or technology.
        functions = []
        for unit in active[active["role"] == "merged"].itertuples(index=False):
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            values = C.read_values(C.project_path(unit.gene_matrix), list(own.loc[own["estimation"], "sample_accession"]), order)
            functions.append(np.mean([C.quantile_function(values[:, column], grid) for column in range(values.shape[1])], axis=0))
        reference = np.mean(functions, axis=0)
        C.require(np.all(np.diff(reference) >= -1e-9), f"The reference distribution of {group} is not ordered")
        references.append(pd.DataFrame({"group": group, "position": grid, "value": reference}))
        folder = C.ALIGNED / group
        folder.mkdir(parents=True, exist_ok=True)
        expected = {C.aligned_path(group, unit.study_accession, unit.platform) for unit in active.itertuples(index=False)}
        for stale in folder.glob("*.parquet"):
            if stale not in expected:
                stale.unlink()
        aligned_records, kept = 0, []
        for unit in active.itertuples(index=False):
            own = C.unit_records(records, group, unit.study_accession, unit.platform)
            columns = list(own["sample_accession"])
            values = C.read_values(C.project_path(unit.gene_matrix), columns, order)
            aligned = np.column_stack([C.align_sample(values[:, column], grid, reference) for column in range(len(columns))])
            # Within every sample the order of the genes is that of its own values.
            kept.append(min(float(pd.Series(values[:, column]).corr(pd.Series(aligned[:, column]), method="spearman")) for column in range(len(columns))))
            frame = pd.DataFrame(aligned.astype("float32"), columns=columns)
            frame.insert(0, "hgnc_id", order)
            path = C.aligned_path(group, unit.study_accession, unit.platform)
            C.save_table(frame, path)
            outputs.append(path)
            aligned_records += len(columns)
        C.require(all(value > C.CONFIG["alignment"]["minimum_rank_correlation"] for value in kept), f"Alignment changed the order of genes within a sample of {group}")
        summary[group] = dict(units=len(active), reference_units=len(functions), records=aligned_records, genes=len(order),
                              reference_lowest=round(float(reference[0]), 6), reference_median=round(float(np.median(reference)), 6),
                              reference_highest=round(float(reference[-1]), 6), lowest_rank_correlation_within_a_sample=round(min(kept), 8))
        print(f"[02] {group}: {aligned_records:,} records of {len(active)} units aligned to a reference of {len(functions)} units over {len(order):,} genes", flush=True)
    references = pd.concat(references, ignore_index=True)
    C.save_table(references, C.DATA / "alignment_reference.parquet")
    C.finish(2, start, [C.DATA / "alignment_reference.parquet", *outputs], method=C.CONFIG["alignment"]["method"], groups=summary,
             grid_points=points, protocol=C.CONFIG["protocol"])
    print("[02] Complete: " + "; ".join(f"{group}: {value['records']:,} records" for group, value in summary.items()), flush=True)


if __name__ == "__main__":
    main()
