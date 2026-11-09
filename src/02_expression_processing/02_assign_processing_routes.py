"""Apply the recorded processing route to every study platform and verify each route against the inventory."""
import collections
import importlib

C = importlib.import_module("01_inventory_expression_files")

TABLES = ("table", "workbook_sheet")
KEYWORDS = ("per_sample_files", "ncbi_counts", "series_matrix", "none")
REASONS = {"per_sample_files": "no_own_data_file", "ncbi_counts": "not_in_ncbi_counts", "series_matrix": "no_series_matrix_values",
           "table": "not_named_in_table", "none": "not_processed"}
ROUTE_COLUMNS = ["source", "accession", "platform", "family", "platform_title", "samples", "tier", "data_source", "covered", "uncovered",
                 "value_kind", "needs_label_mapping", "flag", "basis"]
SAMPLE_COLUMNS = ["source", "study_accession", "sample_accession", "platform", "family", "tier", "data_source", "in_source", "reason"]


def route_evidence(group, decision, own, tables, policy):
    """Find which samples of one study platform the recorded data source holds, the kind of values there, and every contradiction."""
    tier, source, family = decision.get("tier"), str(decision.get("data_source", "")), group["family"].iloc[0]
    rule, kind = policy["tiers"].get(tier), source if source in KEYWORDS else "table"
    if rule is None or not source or not str(decision.get("basis", "")).strip() or decision.get("flag", "uncertain") != "uncertain":
        return None, set(), ["a route needs a known tier, a data source, a basis and no flag other than uncertain"]
    problems, kinds = [], set()
    if family not in rule["families"] or kind not in rule["data_sources"]:
        problems.append(f"tier {tier} does not allow family {family} with data source {kind}")
    if kind == "per_sample_files":
        formats, held = rule.get("formats", {}).get(family, []), []
        for sample in group["sample_accession"]:
            files = [(name, values) for name, values in own.get(sample, ()) if name in formats and (name not in TABLES or values != "none")]
            held.append(bool(files))
            kinds.update(values for name, values in files if name in TABLES)
    elif kind == "ncbi_counts":
        held, kinds = group["ncbi_counts"].tolist(), {"whole_numbers"}
    elif kind == "series_matrix":
        held = (group["series_matrix_rows"] > 0).tolist()
    elif kind == "table":
        names = source.split(";")
        missing = [name for name in names if name not in tables]
        if missing:
            problems.append(f"data source is not a table of this study: {missing}")
        kinds = {tables[name] for name in names if name in tables}
        held = [bool(set(names) & set(files.split(";"))) for files in group["author_table_files"]]
    else:
        held = [False] * len(group)
    flagged = "flag" in decision
    if kind != "none" and not any(held) and not flagged:
        problems.append("the data source holds no sample of this study platform; flag it uncertain while its labels await mapping")
    if tier == "counts" and kinds - {"whole_numbers"} and not flagged:
        problems.append("counts that are not whole numbers must be flagged uncertain")
    # A route may not sit above a lower level that every sample of the study platform has, unless it records why that level cannot be used.
    raw = [name for name in policy["tiers"]["raw_array_files"]["formats"].get(family, []) if name not in TABLES]
    usable = not str(decision.get("raw_not_usable", "")).strip()
    if tier != "raw_array_files" and raw and usable and all(any(name in raw for name, _ in own.get(sample, ())) for sample in group["sample_accession"]):
        problems.append("every sample has a raw array file, so the tier must be raw_array_files, or the route must record raw_not_usable")
    if tier == "raw_array_files" and not usable:
        problems.append("a route that reads the raw array files cannot record raw_not_usable")
    if tier != "counts" and family in policy["tiers"]["counts"]["families"] and group["ncbi_counts"].all():
        problems.append("NCBI counts cover every sample, so the tier must be counts")
    return held, kinds, problems


def main():
    import pandas as pd

    print(f"[02] Assigning processing routes v{C.VERSION}", flush=True)
    start = C.begin(2)
    policy = C.CONFIG["routes"]
    decisions, order = C.REGISTER["routes"], list(policy["tiers"])
    inventory = pd.read_parquet(C.DATA / "sample_expression.parquet")
    listing = pd.read_parquet(C.DATA / "expression_files.parquet")
    platforms = pd.read_parquet(C.DATA / "study_platforms.parquet")
    keys = [f"{accession}|{platform}" for accession, platform in zip(platforms["accession"], platforms["platform"])]
    unrouted = sorted(set(keys) - set(decisions))
    C.require(not unrouted, f"{len(unrouted)} study platforms without a route; record one in {C.CONFIG['register']}: {unrouted[:20]}")
    unknown = sorted(set(decisions) - set(keys))
    C.require(not unknown, f"Routes recorded for study platforms that are not in the inventory: {unknown[:20]}")

    own, tables = collections.defaultdict(lambda: collections.defaultdict(list)), collections.defaultdict(dict)
    for row in listing[listing["format"] != "tar_archive"].itertuples(index=False):
        if row.sample_accession:
            own[row.accession][row.sample_accession].append((row.format, row.value_kind))
        elif row.format in TABLES and row.kind != "ncbi_generated_counts":
            tables[row.accession]["/".join(part for part in (row.file, row.member, row.sheet) if part)] = row.value_kind

    routes, rows, problems = [], [], []
    for item, key in zip(platforms.itertuples(index=False), keys):
        group, decision = inventory[(inventory["study_accession"] == item.accession) & (inventory["platform"] == item.platform)], decisions[key]
        held, kinds, found = route_evidence(group, decision, own[item.accession], tables[item.accession], policy)
        problems += [f"{key}: {problem}" for problem in found]
        if held is None:
            continue
        kind = decision["data_source"] if decision["data_source"] in KEYWORDS else "table"
        routes.append(dict(source=item.source, accession=item.accession, platform=item.platform, family=item.family,
                           platform_title=item.platform_title, samples=len(group), tier=decision["tier"], data_source=decision["data_source"],
                           covered=sum(held), uncovered=len(held) - sum(held), value_kind=";".join(sorted(kinds)),
                           needs_label_mapping=kind == "table" and sum(held) < len(held), flag=decision.get("flag", ""), basis=decision["basis"]))
        for row, inside in zip(group.itertuples(index=False), held):
            rows.append(dict(source=row.source, study_accession=row.study_accession, sample_accession=row.sample_accession, platform=row.platform,
                             family=row.family, tier=decision["tier"], data_source=decision["data_source"], in_source=bool(inside),
                             reason="" if inside else REASONS[kind]))
    C.require(not problems, f"{len(problems)} routes contradict the inventory:\n  " + "\n  ".join(problems[:40]))

    routes = pd.DataFrame(routes, columns=ROUTE_COLUMNS)
    routes["tier"] = pd.Categorical(routes["tier"], categories=order, ordered=True)
    routes = routes.sort_values(["tier", "source", "accession", "platform"]).reset_index(drop=True)
    routes["tier"] = routes["tier"].astype(str)
    samples = pd.DataFrame(rows, columns=SAMPLE_COLUMNS)
    C.require(len(samples) == len(inventory) and routes["samples"].sum() == len(inventory) and routes["covered"].sum() == samples["in_source"].sum(),
              "Routes do not reconcile with the inventoried samples")
    C.save_table(routes, C.DATA / "processing_routes.parquet")
    C.save_table(samples, C.DATA / "sample_routes.parquet")
    C.save_review(routes, C.OUT / "routes_review.csv")
    summary = [dict(tier=tier, study_platforms=len(part), samples=int(part["samples"].sum()), covered=int(part["covered"].sum()))
               for tier in order for part in [routes[routes["tier"] == tier]] if len(part)]
    reasons = samples.loc[~samples["in_source"], "reason"].value_counts().sort_index()
    C.finish(2, start, [C.DATA / "processing_routes.parquet", C.DATA / "sample_routes.parquet", C.OUT / "routes_review.csv"],
             method=policy["method"], by_tier=summary,
             samples=dict(total=len(samples), in_source=int(samples["in_source"].sum()), outside_source={str(key): int(value) for key, value in reasons.items()}),
             flagged_uncertain=int((routes["flag"] == "uncertain").sum()),
             needs_label_mapping=[f"{row.accession}|{row.platform}" for row in routes[routes["needs_label_mapping"]].itertuples(index=False)],
             tiers=policy["tiers"])
    print(f"[02] Complete: {len(routes)} study platforms routed; {int(samples['in_source'].sum()):,} of {len(samples):,} samples held by their data source; "
          f"review sheet: {(C.OUT / 'routes_review.csv').relative_to(C.ROOT).as_posix()}", flush=True)
    for item in summary:
        print(f"[02]   {item['tier']}: {item['study_platforms']} study platforms, {item['covered']:,} of {item['samples']:,} samples", flush=True)


if __name__ == "__main__":
    main()
