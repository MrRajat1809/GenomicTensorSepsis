"""Compare the gene-level values of every two records of one sample that is deposited in several study platforms."""
import collections
import importlib
import itertools

C = importlib.import_module("01_index_validation_inputs")

PAIR_COLUMNS = ["cohort_group", "sample_id", "linked_by", "study_a", "sample_a", "platform_a", "study_b", "sample_b", "platform_b", "shared_genes", "spearman",
                "centred_correlation", "rank_among_b", "candidates_b", "rank_among_a", "candidates_a", "best_other_correlation", "most_similar"]
TABLE_COLUMNS = ["cohort_group", "study_a", "platform_a", "study_b", "platform_b", "pairs", "shared_genes", "lowest_spearman", "median_spearman",
                 "lowest_centred_correlation", "median_centred_correlation", "most_similar", "not_most_similar"]


def unit_columns(frame):
    """Centre every column and scale it to unit length, so that the product of two such tables holds the correlations of their columns."""
    import numpy as np

    values = frame.to_numpy(dtype="float64")
    values = values - values.mean(axis=0)
    return values / np.sqrt((values ** 2).sum(axis=0))


def main():
    import numpy as np
    import pandas as pd

    print(f"[03] Comparing repeated deposits v{C.VERSION}", flush=True)
    start = C.begin(3)
    samples = pd.read_parquet(C.DATA / "validation_samples.parquet")
    matrices = pd.read_parquet(C.DATA / "validation_matrices.parquet").set_index(["study_accession", "platform"])
    repeated = samples[samples["records_of_sample"] > 1]
    by_tables = collections.defaultdict(list)
    for sample_id, part in repeated.groupby("sample_id", sort=True):
        records = sorted(part[["study_accession", "platform", "sample_accession"]].itertuples(index=False, name=None))
        for first, second in itertools.combinations(records, 2):
            by_tables[(part["cohort_group"].iloc[0], first[:2], second[:2])].append((sample_id, part["linked_by"].iloc[0], first[2], second[2]))
    expected = int(sum(count * (count - 1) // 2 for count in repeated.groupby("sample_id").size()))
    C.require(sum(map(len, by_tables.values())) == expected, "Every two records of one sample must be compared once")

    rows, loaded, group_loaded = [], {}, None
    for (group, table_a, table_b), pairs in sorted(by_tables.items()):
        if group != group_loaded:
            loaded, group_loaded = {}, group
        for table in (table_a, table_b):
            if table not in loaded:
                entry = matrices.loc[table]
                frame = pd.read_parquet(C.project_path(entry["gene_matrix"])).drop(columns="symbol").set_index("hgnc_id").astype("float64")
                frame = C.log_values(frame, entry["value_scale"]) if entry["comparable_scale"] else frame
                loaded[table] = (frame.dropna(axis=0, how="any"), bool(entry["comparable_scale"]))
        (first, first_log), (second, second_log) = loaded[table_a], loaded[table_b]
        shared = first.index.intersection(second.index).sort_values()
        C.require(len(shared) > 0, f"Two study platforms that share a sample share no gene: {table_a} {table_b}")
        first, second, same = first.loc[shared], second.loc[shared], table_a == table_b
        ranks_a = unit_columns(first[sorted({pair[2] for pair in pairs})].rank())
        ranks_b = unit_columns(second[sorted({pair[3] for pair in pairs})].rank())
        place_a = {name: position for position, name in enumerate(sorted({pair[2] for pair in pairs}))}
        place_b = {name: position for position, name in enumerate(sorted({pair[3] for pair in pairs}))}
        # Every gene is centred on its median within its own matrix, so that what two samples share beyond the average blood profile decides.
        similar = None
        if first_log and second_log:
            similar = unit_columns(first.sub(first.median(axis=1), axis=0)).T @ unit_columns(second.sub(second.median(axis=1), axis=0))
        column_a, column_b = {name: position for position, name in enumerate(first.columns)}, {name: position for position, name in enumerate(second.columns)}
        for sample_id, linked_by, sample_a, sample_b in pairs:
            row = dict(cohort_group=group, sample_id=sample_id, linked_by=linked_by, study_a=table_a[0], sample_a=sample_a, platform_a=table_a[1], study_b=table_b[0],
                       sample_b=sample_b, platform_b=table_b[1], shared_genes=len(shared),
                       spearman=float(ranks_a[:, place_a[sample_a]] @ ranks_b[:, place_b[sample_b]]), centred_correlation=np.nan, rank_among_b=0,
                       candidates_b=0, rank_among_a=0, candidates_a=0, best_other_correlation=np.nan, most_similar=False)
            if similar is not None:
                i, j = column_a[sample_a], column_b[sample_b]
                value = float(similar[i, j])
                others_b = np.delete(similar[i, :], [j, i] if same else [j])
                others_a = np.delete(similar[:, j], [i, j] if same else [i])
                row.update(centred_correlation=value, rank_among_b=1 + int((others_b > value).sum()), candidates_b=len(others_b) + 1,
                           rank_among_a=1 + int((others_a > value).sum()), candidates_a=len(others_a) + 1,
                           best_other_correlation=float(others_b.max()) if len(others_b) else np.nan)
                row["most_similar"] = row["rank_among_b"] == 1 and row["rank_among_a"] == 1
            rows.append(row)
        print(f"[03]   {group} {table_a[0]}|{table_a[1]} and {table_b[0]}|{table_b[1]}: {len(pairs):,} pairs over {len(shared):,} genes", flush=True)
    table = pd.DataFrame(rows, columns=PAIR_COLUMNS)
    compared = table[table["centred_correlation"].notna()]
    keys = ["cohort_group", "study_a", "platform_a", "study_b", "platform_b"]
    summary = table.groupby(keys, sort=True).agg(pairs=("sample_id", "size"), shared_genes=("shared_genes", "first"), lowest_spearman=("spearman", "min"),
                                                 median_spearman=("spearman", "median"), lowest_centred_correlation=("centred_correlation", "min"),
                                                 median_centred_correlation=("centred_correlation", "median"), most_similar=("most_similar", "sum")).reset_index()
    summary["not_most_similar"] = summary["pairs"] - summary["most_similar"]
    summary = summary[TABLE_COLUMNS]
    C.save_table(table, C.DATA / "repeat_pairs.parquet")
    C.save_table(summary, C.DATA / "repeat_pairs_by_matrices.parquet")
    C.save_review(summary, C.OUT / "repeat_pairs_review.csv")
    C.save_review(table[~table["most_similar"]], C.OUT / "repeat_pairs_not_most_similar_review.csv")
    C.finish(3, start, [C.DATA / "repeat_pairs.parquet", C.DATA / "repeat_pairs_by_matrices.parquet", C.OUT / "repeat_pairs_review.csv",
                        C.OUT / "repeat_pairs_not_most_similar_review.csv"], method=C.CONFIG["repeats"]["method"],
             samples=dict(deposited_more_than_once=int(repeated["sample_id"].nunique()), records=len(repeated)),
             pairs=dict(total=len(table), compared_after_centring=len(compared), pairs_of_matrices=len(summary),
                        by_link={str(name): int(count) for name, count in table["linked_by"].value_counts().sort_index().items()}),
             spearman=dict(lowest=round(float(table["spearman"].min()), 6), first_percentile=round(float(table["spearman"].quantile(0.01)), 6),
                           median=round(float(table["spearman"].median()), 6),
                           **{C.at_least("at_least", limit): int((table["spearman"] >= limit).sum()) for limit in C.CONFIG["repeats"]["report"]["spearman_at_least"]}),
             centred_correlation=dict(lowest=round(float(compared["centred_correlation"].min()), 6), median=round(float(compared["centred_correlation"].median()), 6)),
             most_similar=dict(pairs=int(table["most_similar"].sum()), not_most_similar=int((~table["most_similar"]).sum()),
                               share=round(float(table["most_similar"].mean()), 6)),
             protocol=C.CONFIG["protocol"])
    print(f"[03] Complete: {len(table):,} pairs of records of {int(repeated['sample_id'].nunique()):,} samples; Spearman median "
          f"{table['spearman'].median():.4f}, lowest {table['spearman'].min():.4f}; the linked record is the most similar one in "
          f"{int(table['most_similar'].sum()):,} pairs; review sheets: {(C.OUT / 'repeat_pairs_review.csv').relative_to(C.ROOT).as_posix()}, "
          f"{(C.OUT / 'repeat_pairs_not_most_similar_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)


if __name__ == "__main__":
    main()
