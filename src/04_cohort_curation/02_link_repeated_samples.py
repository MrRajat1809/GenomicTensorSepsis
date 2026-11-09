"""Link the records of one array or library that is deposited in several study platforms, and group the study platforms that share samples."""
import collections
import importlib

C = importlib.import_module("01_catalogue_sample_metadata")

LINK_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "sample_id", "primary_study_accession", "primary_sample_accession",
                "is_primary_record", "records_of_sample", "linked_by", "values_agree", "cohort_group"]
CANDIDATE_COLUMNS = ["family", "study_a", "platform_a", "sample_a", "study_b", "platform_b", "sample_b", "correlation", "margin", "same_study", "linked_by_files",
                     "decision"]
GROUP_COLUMNS = ["cohort_group", "source", "study_accession", "platform", "study_platforms_in_group", "release_date", "records", "primary_records",
                 "repeated_records", "shares_samples_with"]


def evidence(expression):
    """Return, for every sample that a matrix holds, what its files record about the physical array or library, and the kind of evidence it is."""
    import pandas as pd

    found = {}
    for row in pd.read_parquet(expression / "affymetrix_arrays.parquet").itertuples(index=False):
        if row.scan_time:
            found[(row.study_accession, row.sample_accession)] = ("chip_type_scan_time_and_median_intensity", (row.chip_type, row.scan_time, row.raw_pm_median), None)
    for row in pd.read_parquet(expression / "single_channel_arrays.parquet").itertuples(index=False):
        if row.in_matrix and row.array_id:
            found[(row.study_accession, row.sample_accession)] = ("slide_identifier_and_median_intensities",
                                                                 (row.layout, row.array_id, row.foreground_median, row.background_median), None)
    for row in pd.read_parquet(expression / "intensity_table_arrays.parquet").itertuples(index=False):
        if row.array_id:
            # A BeadChip position is read by the scanner, so the identifier alone names the array; the median shows whether the exports agree.
            found[(row.study_accession, row.sample_accession)] = ("beadchip_barcode_and_position", (row.array_id,), row.signal_median)
    for row in pd.read_parquet(expression / "count_table_samples.parquet").itertuples(index=False):
        found[(row.study_accession, row.sample_accession)] = ("identical_counts", (row.counts_digest,), None)
    return found


def centred_ranks(frame):
    """Rank the genes within each sample, remove the mean rank of each sample, and scale each sample to unit variance; missing genes count as the mean."""
    import numpy as np

    ranks = frame.rank().to_numpy(dtype="float64")
    ranks = np.nan_to_num(ranks - np.nanmean(ranks, axis=0), nan=0.0)
    spread = ranks.std(axis=0)
    return ranks / np.where(spread > 0, spread, 1.0)


def identity_candidates(policy, group_of, order):
    """Every pair of samples of one technology family that are each other's most similar sample above the limits, by rank correlation of their gene
    profiles over genes that every table of the family holds. Of the records that the files tie to one sample only the first by order is scored, so
    that a further copy of a sample cannot stand beside it as the next most similar sample and hide the pair."""
    import numpy as np
    import pandas as pd

    tables = pd.read_parquet(C.ANNOTATION / "gene_matrix_summary.parquet")
    tables = tables[tables["family"].isin(policy["families"]) & tables["value_scale"].isin(policy["value_scales"])].sort_values(["family", "accession", "platform"])
    found, screened = [], {}
    for family, part in tables.groupby("family", sort=True):
        held = [set(pd.read_parquet(C.project_path(row.gene_matrix), columns=["hgnc_id"])["hgnc_id"]) for row in part.itertuples(index=False)]
        common = sorted(set.intersection(*held), key=lambda gene: int(gene.split(":")[1]))
        C.require(len(common) >= policy["minimum_common_genes"], f"The tables of {family} share only {len(common):,} genes")
        # Genes spread evenly over the sorted identifiers, so that the choice depends on nothing but the identifiers.
        chosen = [common[position] for position in np.unique(np.linspace(0, len(common) - 1, min(len(common), policy["genes_used"])).astype(int))]
        blocks, owners = [], []
        for row in part.itertuples(index=False):
            frame = pd.read_parquet(C.project_path(row.gene_matrix)).set_index("hgnc_id").drop(columns="symbol").loc[chosen]
            blocks.append(centred_ranks(frame))
            owners += [(row.accession, row.platform, sample) for sample in frame.columns]
        C.require(all((study, sample) in group_of for study, _, sample in owners), f"A column of a gene matrix of {family} is not an included sample")
        earliest = {}
        for position, (study, _, sample) in enumerate(owners):
            group = group_of[(study, sample)]
            if group not in earliest or order(study, sample) < order(owners[earliest[group]][0], owners[earliest[group]][2]):
                earliest[group] = position
        kept = sorted(earliest.values())
        screened[family] = dict(tables=len(part), records=len(owners), distinct_samples=len(kept), genes=len(chosen))
        if len(kept) < 2:
            print(f"[02] identity screen of {family}: {len(kept)} distinct sample, no pair to compare", flush=True)
            continue
        data = np.hstack(blocks)[:, kept]
        owners = [owners[position] for position in kept]
        scores = data.T @ data / data.shape[0]
        np.fill_diagonal(scores, -np.inf)
        best, top = scores.argmax(axis=1), scores.max(axis=1)
        second = np.partition(scores, -2, axis=1)[:, -2]
        for first in range(len(owners)):
            other = best[first]
            margin = top[first] - max(second[first], second[other])
            if first < other and best[other] == first and top[first] >= policy["candidate_correlation"] and margin >= policy["candidate_margin"]:
                found.append(dict(family=family, study_a=owners[first][0], platform_a=owners[first][1], sample_a=owners[first][2], study_b=owners[other][0],
                                  platform_b=owners[other][1], sample_b=owners[other][2], correlation=round(float(top[first]), 6), margin=round(float(margin), 6)))
        print(f"[02] identity screen of {family}: {len(part)} tables, {screened[family]['records']:,} records of {len(owners):,} distinct samples, "
              f"{len(chosen):,} genes", flush=True)
        del data, scores
    return pd.DataFrame(found, columns=CANDIDATE_COLUMNS[:9]), screened


def record_pairs(rules):
    """Pairs of records that a study names as re-analyses of records of another study, joined by the key that a pattern reads from a field of each."""
    import re

    import pandas as pd

    records = pd.read_parquet(C.ACQUIRED / "samples.parquet", columns=["study_accession", "sample_accession", "title", "description", "source_name"])
    pairs = []
    for rule in rules:
        keyed = []
        for side in (rule["a"], rule["b"]):
            own = records[records["study_accession"] == side["study"]]
            keys = {}
            for sample, value in zip(own["sample_accession"], own[side["field"]].fillna("")):
                match = re.search(side["pattern"], str(value))
                if match:
                    C.require(match.group(1) not in keys, f"Two records of {side['study']} give the key {match.group(1)}")
                    keys[match.group(1)] = sample
            keyed.append(keys)
        # The key pattern must read every record that its study marks as a re-analysis, so that a pattern too narrow cannot pass unnoticed.
        own = records[records["study_accession"] == rule["b"]["study"]]
        marked = {sample for sample, value in zip(own["sample_accession"], own[rule["b"]["field"]].fillna("")) if re.search(rule["b"]["required"], str(value))}
        C.require(marked and marked == set(keyed[1].values()), f"The key pattern of {rule['b']['study']} must read exactly the records that it marks as "
                                                               f"re-analyses; not read: {sorted(marked - set(keyed[1].values()))[:10]}")
        shared = sorted(set(keyed[0]) & set(keyed[1]))
        C.require(shared and len(shared) == len(keyed[1]), f"Every record of {rule['b']['study']} that the rule names must have its counterpart in {rule['a']['study']}")
        pairs += [((rule["a"]["study"], keyed[0][key]), (rule["b"]["study"], keyed[1][key])) for key in shared]
    return pairs


def main():
    import pandas as pd

    print(f"[02] Linking repeated samples v{C.VERSION}", flush=True)
    start = C.begin(2)
    policy = C.CONFIG["links"]
    index = C.sample_records(["source", "study_accession", "sample_accession", "platform"])
    released = dict(pd.read_parquet(C.ACQUIRED / "included_studies.parquet", columns=["accession", "release_date"]).itertuples(index=False, name=None))
    C.require(set(index["study_accession"]) <= set(released) and all(released[study] for study in set(index["study_accession"])),
              "An included study has no release date")
    found = evidence(C.EXPRESSION)
    inside = set(zip(index["study_accession"], index["sample_accession"]))
    C.require(set(found) <= inside, "The expression stage records an array or library for a record without expression values")

    # Records with the same evidence are one sample; a sample without evidence is its own.
    groups = collections.defaultdict(list)
    for row in index.itertuples(index=False):
        key = (row.study_accession, row.sample_accession)
        kind, value, _ = found.get(key, ("", key, None))
        groups[(kind, value)].append((row.study_accession, row.platform, row.sample_accession))

    def order(study, sample):
        # The study released first, then the lower accession by its letters and number: the order that names the primary record of a sample.
        return released[study], C.accession_order(study), C.accession_order(sample)

    # Near-identical profiles are candidates; they join records only where links.identity names the study platforms, and mark replicate libraries.
    identity = policy["identity"]
    group_of = {(study, sample): group for group, members in groups.items() for study, _, sample in members}
    candidates, screened = identity_candidates(identity, group_of, order)
    platform_of = dict(zip(zip(index["study_accession"], index["sample_accession"]), index["platform"]))
    candidates["same_study"] = candidates["study_a"] == candidates["study_b"]
    candidates["linked_by_files"] = [group_of[(a, x)] == group_of[(b, y)] and group_of[(a, x)][0] != ""
                                     for a, x, b, y in zip(candidates["study_a"], candidates["sample_a"], candidates["study_b"], candidates["sample_b"])]
    candidates["decision"] = ""
    extra, replicates, found_by_entry = [], [], []

    def found_by(entry, chosen):
        return dict(**{name: entry[name] for name in ("a", "b", "study") if name in entry}, pairs=len(chosen),
                    lowest_correlation=float(chosen["correlation"].min()), highest_correlation=float(chosen["correlation"].max()),
                    lowest_margin=float(chosen["margin"].min()))

    for entry in C.REGISTER["same_sample"]:
        sides = {entry["a"], entry["b"]}
        named = pd.Series([{f"{a}|{p}", f"{b}|{q}"} == sides or (len(sides) == 1 and f"{a}|{p}" in sides and f"{b}|{q}" in sides)
                           for a, p, b, q in zip(candidates["study_a"], candidates["platform_a"], candidates["study_b"], candidates["platform_b"])],
                          index=candidates.index, dtype=bool)
        chosen = candidates[named & (candidates["correlation"] >= entry["minimum_correlation"])]
        C.require(len(chosen), f"No identity candidate links {entry['a']} and {entry['b']} above {entry['minimum_correlation']}")
        candidates.loc[chosen.index, "decision"] = "same_sample"
        found_by_entry.append(found_by(entry, chosen))
        extra += [((row.study_a, row.sample_a), (row.study_b, row.sample_b), "expression_identity") for row in chosen.itertuples(index=False)]
    for entry in C.REGISTER["replicates"]:
        study, platform = entry["study"].split("|")
        chosen = candidates[(candidates["study_a"] == study) & (candidates["study_b"] == study) & (candidates["platform_a"] == platform)
                            & (candidates["platform_b"] == platform) & (candidates["correlation"] >= entry["minimum_correlation"])]
        C.require(len(chosen), f"No replicate libraries found in {entry['study']} above {entry['minimum_correlation']}")
        candidates.loc[chosen.index, "decision"] = "replicate_libraries"
        found_by_entry.append(found_by(entry, chosen))
        replicates += [dict(study_accession=study, platform=platform, first=min(row.sample_a, row.sample_b, key=C.accession_order),
                            second=max(row.sample_a, row.sample_b, key=C.accession_order), correlation=row.correlation, margin=row.margin)
                       for row in chosen.itertuples(index=False)]
    record_links = record_pairs(C.REGISTER["record_links"])
    extra += [(first, second, "record_names_a_re_analysis") for first, second in record_links]
    outside = sorted({key for pair in extra for key in pair[:2] if key not in platform_of})
    C.require(not outside, f"A declared link names records that are not included records with expression values: {outside[:10]}; a record without values "
                           f"is listed in records_without_values.parquet and joins no sample, so revise the entry in {C.CONFIG['register']}")

    record_parent = {}

    def root(key):
        record_parent.setdefault(key, key)
        while record_parent[key] != key:
            record_parent[key] = record_parent[record_parent[key]]
            key = record_parent[key]
        return key

    kinds = collections.defaultdict(set)
    for (kind, _), members in groups.items():
        for study, _, sample in members[1:]:
            record_parent[root((study, sample))] = root((members[0][0], members[0][2]))
    for first, second, kind in extra:
        record_parent[root(second)] = root(first)
    for (kind, _), members in groups.items():
        if kind and len(members) > 1:
            kinds[root((members[0][0], members[0][2]))].add(kind)
    for first, _, kind in extra:
        kinds[root(first)].add(kind)
    merged = collections.defaultdict(list)
    for row in index.itertuples(index=False):
        merged[root((row.study_accession, row.sample_accession))].append((row.study_accession, row.platform, row.sample_accession))
    medians = {(study, sample): item[2] for (study, sample), item in found.items()}
    links, parent = {}, {key: key for key in set(zip(index["study_accession"], index["platform"]))}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    for key, members in merged.items():
        primary = min(members, key=lambda member: order(member[0], member[2]))
        # The median intensities of two exports of one array must agree; deposits joined by their values or their records have no such check.
        values = [medians[(study, sample)] for study, _, sample in members if medians.get((study, sample)) is not None]
        agree = len(members) == 1 or len(values) < 2 or max(values) - min(values) <= policy["median_tolerance"]
        for study, platform, sample in members:
            links[(study, sample)] = dict(sample_id=f"{primary[0]}:{primary[2]}", primary_study_accession=primary[0], primary_sample_accession=primary[2],
                                          is_primary_record=(study, sample) == (primary[0], primary[2]), records_of_sample=len(members),
                                          linked_by="; ".join(sorted(kinds[key])) if len(members) > 1 else "", values_agree=bool(agree))
            parent[find((study, platform))] = find((primary[0], primary[1]))
    members = collections.defaultdict(list)
    for key in sorted(parent):
        members[find(key)].append(key)
    names = {root: f"G{number:03d}" for number, root in enumerate(sorted(members, key=lambda root: min(members[root])), 1)}
    table = pd.DataFrame([dict(source=row.source, study_accession=row.study_accession, sample_accession=row.sample_accession, platform=row.platform,
                               cohort_group=names[find((row.study_accession, row.platform))], **links[(row.study_accession, row.sample_accession)])
                          for row in index.itertuples(index=False)], columns=LINK_COLUMNS)
    C.require(len(table) == len(index) and int(table["is_primary_record"].sum()) == table["sample_id"].nunique()
              and (table.groupby("sample_id")["is_primary_record"].sum() == 1).all() and (table.groupby("sample_id")["cohort_group"].nunique() == 1).all(),
              "Every sample must have exactly one primary record, in one cohort group")

    sources = dict(zip(zip(index["study_accession"], index["platform"]), index["source"]))
    rows = []
    for root, keys in members.items():
        for study, platform in keys:
            part = table[(table["study_accession"] == study) & (table["platform"] == platform)]
            rows.append(dict(cohort_group=names[root], source=sources[(study, platform)], study_accession=study, platform=platform, study_platforms_in_group=len(keys),
                             release_date=released[study], records=len(part), primary_records=int(part["is_primary_record"].sum()),
                             repeated_records=int((~part["is_primary_record"]).sum()),
                             shares_samples_with="; ".join(f"{other}|{chip}" for other, chip in keys if (other, chip) != (study, platform))))
    grouped = pd.DataFrame(rows, columns=GROUP_COLUMNS).sort_values(["cohort_group", "release_date", "study_accession", "platform"]).reset_index(drop=True)
    replicate_table = pd.DataFrame(replicates, columns=["study_accession", "platform", "first", "second", "correlation", "margin"])
    C.require(not set(replicate_table["second"]) & set(replicate_table["first"]), "A replicate library may not be the first of another pair")
    C.save_table(table, C.DATA / "sample_links.parquet")
    C.save_table(grouped, C.DATA / "cohort_groups.parquet")
    C.save_table(candidates[CANDIDATE_COLUMNS], C.DATA / "identity_candidates.parquet")
    C.save_table(replicate_table, C.DATA / "replicate_libraries.parquet")
    C.save_review(grouped, C.OUT / "cohort_groups_review.csv")
    C.save_review(candidates[CANDIDATE_COLUMNS], C.OUT / "identity_candidates_review.csv")
    repeated = table[table["records_of_sample"] > 1]
    C.finish(2, start, [C.DATA / "sample_links.parquet", C.DATA / "cohort_groups.parquet", C.DATA / "identity_candidates.parquet", C.DATA / "replicate_libraries.parquet",
                        C.OUT / "cohort_groups_review.csv", C.OUT / "identity_candidates_review.csv"], method=policy["method"], evidence_added=policy["evidence_added"],
             identity=dict(candidates=len(candidates), already_linked_by_files=int(candidates["linked_by_files"].sum()),
                           within_one_study=int(candidates["same_study"].sum()),
                           by_decision={str(name or "candidate_only"): int(count) for name, count in candidates["decision"].value_counts().sort_index().items()},
                           screened=screened, found_by_entry=found_by_entry,
                           policy=identity, same_sample=C.REGISTER["same_sample"], replicates=C.REGISTER["replicates"]),
             record_links=dict(pairs=len(record_links), rules=C.REGISTER["record_links"]),
             records=dict(total=len(table)),
             samples=dict(distinct=int(table["sample_id"].nunique()),
                          deposited_more_than_once=int(repeated["sample_id"].nunique()), repeated_records=int((~table["is_primary_record"]).sum()),
                          most_records_of_one_sample=int(table["records_of_sample"].max())),
             linked_by={str(kind): int(count) for kind, count in repeated["linked_by"].value_counts().sort_index().items()},
             records_whose_values_differ_between_deposits=int((~table["values_agree"]).sum()),
             cohort_groups=dict(total=len(members), of_several_study_platforms=sum(len(keys) > 1 for keys in members.values()),
                                study_platforms=len(parent)),
             median_tolerance=policy["median_tolerance"], primary_record=policy["primary_record"], protocol=C.CONFIG["protocol"])
    print(f"[02] Complete: {len(table):,} sample records with expression values are {int(table['sample_id'].nunique()):,} distinct samples; "
          f"{int((~table['is_primary_record']).sum()):,} records repeat another; "
          f"{len(parent)} study platforms form {len(members)} cohort groups; review sheet: {(C.OUT / 'cohort_groups_review.csv').relative_to(C.ROOT).as_posix()}",
          flush=True)
    for kind, count in repeated["linked_by"].value_counts().sort_index().items():
        print(f"[02]   {kind}: {count:,} records", flush=True)
    for name, part in grouped[grouped["study_platforms_in_group"] > 1].groupby("cohort_group"):
        print(f"[02]   {name}: {'; '.join(part['study_accession'] + '|' + part['platform'])}: {int(part['records'].sum()):,} records, "
              f"{int(part['primary_records'].sum()):,} distinct samples", flush=True)


if __name__ == "__main__":
    main()
