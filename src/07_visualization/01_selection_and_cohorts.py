"""Draw how the studies and samples were selected and what the resource holds: Fig1, Fig2 and SuppFig1."""
import itertools

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import to_rgb
from matplotlib.patches import Arc, FancyBboxPatch, PathPatch, Rectangle, Wedge
from matplotlib.path import Path

import visualization_common as V

REPOSITORIES = {"ena": "ENA", "geo": "GEO", "arrayexpress": "ArrayExpress", "ega": "EGA", "dbgap": "dbGaP", "node": "NODE"}
# The rings of the atlas that give the share of the records of a study platform that record a curated variable, from the outside in.
RINGS = {"sex": "sex", "age": "age", "patient": "patient identifier", "time_point": "time point", "outcome": "outcome", "septic_shock": "septic shock"}
MATERIAL_COLOR = "#9AA3B5"


def words(code, start=False):
    """Write a reason code in words, with the usual spelling of its abbreviations; start=True begins it with a capital."""
    text = code.replace("_", " ")
    for plain, spelled in (("rnaseq", "RNA-seq"), ("rna seq", "RNA-seq"), ("rna-seq", "RNA-seq"), ("ncbi", "NCBI"), ("mrna", "mRNA"), ("superseries", "SuperSeries"),
                           ("subseries", "SubSeries"), ("transcriptome wide", "transcriptome-wide"), ("cross source", "cross-source")):
        text = text.replace(plain, spelled)
    return text[0].upper() + text[1:] if start else text


def box(ax, x, y, width, height, text, face="#FFFFFF", edge="#BBBBBB", size=9, left=False):
    ax.add_patch(FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.003,rounding_size=0.012", facecolor=face, edgecolor=edge, linewidth=1.0))
    ax.text(x + 0.014 if left else x + width / 2, y + height / 2, text, ha="left" if left else "center", va="center", fontsize=size, color=V.INK, linespacing=1.4)


def arrow(ax, start, end):
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="-|>", color=V.MUTED, linewidth=1.0, shrinkA=0, shrinkB=0, mutation_scale=10))


def reasons(counts, limit=None):
    lines = [f"{V.number(count)}  {words(name)}" for name, count in counts.items()]
    return "\n".join(lines if limit is None else lines[:limit])


def flow_studies(canvas, census, decisions, studies, platforms, with_values):
    """Left column of the flow: from the search results to the included studies and their study platforms, with every reason beside its step."""
    ax = canvas.subplots()
    V.blank(ax)
    V.panel(canvas, ax, "A", "From the repositories to the studies")
    found = census["source"].map(REPOSITORIES).value_counts()
    dropped = census.loc[census["census_status"] != "candidate", "census_status"].value_counts()
    excluded = decisions.loc[decisions["decision"] != "include", "reason"].value_counts()
    assays = studies["assay"].value_counts()
    sources = studies["source"].map(REPOSITORIES).value_counts()
    blue = V.tint("#4A6FE3", 0.86)
    listed = [f"{name} {V.number(count)}" for name, count in found.items()]
    box(ax, 0.04, 0.83, 0.50, 0.14, f"Records found by the searches\n{V.number(len(census))}\n" + ", ".join(listed[:3]) + "\n" + ", ".join(listed[3:]),
        blue, "#4A6FE3", 8.5)
    box(ax, 0.61, 0.70, 0.37, 0.13, f"Not kept as candidates: {V.number(dropped.sum())}\n" + reasons(dropped), "#F7F7F7", size=8, left=True)
    box(ax, 0.04, 0.56, 0.50, 0.10, f"Candidate studies screened\n{V.number(len(decisions))}", blue, "#4A6FE3")
    box(ax, 0.61, 0.05, 0.37, 0.45, f"Excluded at screening: {V.number(excluded.sum())}\n(first reason that applies)\n" + reasons(excluded), "#F7F7F7", size=8, left=True)
    box(ax, 0.04, 0.05, 0.50, 0.17, f"Included studies\n{V.number(len(studies))}\n" + ", ".join(f"{name} {count}" for name, count in sources.items()) + "\n"
        + ", ".join(f"{count} {words(name)}" for name, count in assays.items()) + f"\n{platforms} study platforms, {with_values} with values", V.tint("#55A868", 0.82),
        "#55A868")
    arrow(ax, (0.29, 0.83), (0.29, 0.66))
    arrow(ax, (0.29, 0.765), (0.61, 0.765))
    arrow(ax, (0.29, 0.56), (0.29, 0.22))
    arrow(ax, (0.29, 0.39), (0.61, 0.39))


def flow_samples(canvas, included, ledger, samples, config):
    """Right column of the flow: from the sample records of the included studies to the samples. The records without expression values are listed in
    their ledger with the reason; from the sample records with values on, nothing is removed, so that each box follows from the one above it."""
    ax = canvas.subplots()
    V.blank(ax)
    V.panel(canvas, ax, "B", "From the sample records to the samples")
    distinct = samples[samples["is_primary_record"]]
    without = ledger["reason"].value_counts()
    tables = samples[["study_accession", "platform"]].drop_duplicates().shape[0]
    marked = int(samples["one_per_patient"].sum())
    identified = distinct.loc[distinct["patient_identified"], "patient_key"].nunique()
    V.require(marked == identified + int((~distinct["patient_identified"]).sum()),
              "The samples marked one per patient must be the identified patients and the samples without an identifier")
    green, grey = V.tint("#55A868", 0.82), "#F7F7F7"
    box(ax, 0.04, 0.87, 0.50, 0.10, f"Sample records of the included studies\n{V.number(included)}", green, "#55A868")
    box(ax, 0.61, 0.73, 0.37, 0.16, f"Without expression values: {V.number(len(ledger))}\n(listed with the reason)\n" + reasons(without), grey, size=8, left=True)
    box(ax, 0.04, 0.66, 0.50, 0.10, f"Sample records with expression values\n{V.number(len(samples))} in {tables} tables", green, "#55A868")
    box(ax, 0.61, 0.565, 0.37, 0.08, f"Repeat deposits: {V.number(len(samples) - len(distinct))} records\n(kept and marked)", grey, size=8, left=True)
    box(ax, 0.04, 0.45, 0.50, 0.10, f"Distinct samples\n{V.number(len(distinct))}", green, "#55A868")
    box(ax, 0.61, 0.355, 0.37, 0.08, f"Further samples of a patient: {V.number(len(distinct) - marked)}\n(kept and marked)", grey, size=8, left=True)
    box(ax, 0.04, 0.20, 0.50, 0.14, f"One sample per patient\n{V.number(marked)}\n{V.number(identified)} identified patients and\n"
        f"{V.number((~distinct['patient_identified']).sum())} samples without an identifier", green, "#55A868", 8.5)
    for start, end in (((0.29, 0.87), (0.29, 0.76)), ((0.29, 0.815), (0.61, 0.815)), ((0.29, 0.66), (0.29, 0.55)), ((0.29, 0.605), (0.61, 0.605)),
                       ((0.29, 0.45), (0.29, 0.34)), ((0.29, 0.395), (0.61, 0.395))):
        arrow(ax, start, end)
    before = int((distinct["sepsis_phase"] == "before_onset").sum())
    ax.text(0.04, 0.155, f"Distinct samples by curated class: {V.number(len(distinct))}\n(sepsis includes {V.number(before)} drawn before the onset, "
            "marked by their phase)", fontsize=8.5, color=V.INK, va="center", linespacing=1.4)
    width = 0.94 / len(config["class_order"])
    for position, name in enumerate(config["class_order"]):
        count = int((distinct["sample_class"] == name).sum())
        box(ax, 0.04 + position * width, 0.03, width - 0.012, 0.085, f"{config['class_labels'][name]}\n{V.number(count)}", V.tint(config["class_colors"][name], 0.72),
            config["class_colors"][name], 8)


def growth(canvas, table, config):
    """Distinct samples by the year in which their study platform was released, added up over the years and stacked by technology."""
    ax = canvas.subplots()
    years = table.index.to_numpy()
    families = [name for name in config["family_order"] if name in table.columns]
    running = table[families].cumsum()
    ax.stackplot(years, *(running[name].to_numpy() for name in families), colors=[config["family_colors"][name] for name in families],
                 labels=[config["family_labels"][name] for name in families], alpha=0.88, edgecolor="white", linewidth=0.6)
    total = running.sum(axis=1)
    ax.annotate(f"{V.number(total.iloc[-1])} distinct samples", (years[-1], total.iloc[-1]), xytext=(-4, 6), textcoords="offset points", ha="right",
                va="bottom", fontsize=8, color=V.INK)
    # The searches were run in the last year shown, so it holds only the studies released before them.
    ticks = [year for year in years if (years[-1] - year) % 4 == 0]
    ax.set_xticks(ticks, [f"{year}\n(part)" if year == years[-1] else str(year) for year in ticks])
    ax.set_xlim(years[0], years[-1])
    ax.set_ylim(0, total.max() * 1.16)
    ax.set_ylabel("Distinct samples, added up")
    V.clean_axis(ax, grid="y")
    V.panel(canvas, ax, "A", "Samples by year of release")
    ax.legend(loc="upper left", reverse=True)


def material_flow(canvas, flows, config, material_order):
    """Distinct samples from their technology to the material sampled to their curated class."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_xlim(-0.68, 1.6)
    ax.set_ylim(-0.01, 1.01)
    first, second = flows[0], flows[1]
    families = [name for name in config["family_order"] if (first["source"] == name).any()]
    classes = [name for name in config["class_order"] if (second["target"] == name).any()]
    totals = pd.concat([first.groupby("source")["value"].sum(), second.groupby("target")["value"].sum(), first.groupby("target")["value"].sum()])
    colors = {**{name: config["family_colors"][name] for name in families}, **{name: MATERIAL_COLOR for name in material_order},
              **{name: config["class_colors"][name] for name in classes}}
    labels = {**{name: f"{config['family_labels'][name]}  {V.number(totals[name])}" for name in families},
              **{name: f"{name}\n{V.number(totals[name])}" for name in material_order},
              **{name: f"{config['class_labels'][name]}  {V.number(totals[name])}" for name in classes}}
    # The materials stand right of the middle, their labels to their left over the ribbons from the technologies, clear of both other columns.
    V.alluvial(ax, [families, material_order, classes], [(first, "source"), (second, "target")], colors, labels, xs=[0.0, 0.66, 1.0], gap=0.045,
               inner="left")
    V.panel(canvas, ax, "B", "From technology to material to class")


def atlas(canvas, table, shared, config, top):
    """Every study platform as one sector of a circle, ordered by technology and year of release. From the outside in: the technology, the share of
    every curated class, the sample records on a logarithmic scale, and the share of the records that record each of six curated variables; bands join
    the study platforms that hold the same samples."""
    ax = canvas.subplots()
    ax.set_axis_off()
    # The view is set through the data limits, not fixed limits, so that the equal aspect widens the shorter side without a notice.
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0)
    ax.update_datalim([(-1.5, -1.5), (1.5, 1.5)])
    V.panel(canvas, ax, "C", "Every study platform at a glance")
    family_gap, top_gap = 2.4, 12.0
    changes = int((table["family"] != table["family"].shift()).sum()) - 1
    span = (360 - top_gap - family_gap * changes) / len(table)
    shares = V.cmap()
    groups = list(dict.fromkeys(shared["cohort_group"]))
    # Spectral without its pale middle, so that every band stands out on white.
    tones = np.concatenate([np.linspace(0.0, 0.34, (len(groups) + 1) // 2), np.linspace(0.66, 1.0, len(groups) // 2)])
    group_colors = {name: colormaps["Spectral"](value) for name, value in zip(groups, tones)}
    angle, centres, family = 90 - top_gap / 2, {}, None
    largest = np.log10(top)
    for row in table.itertuples(index=False):
        if family is not None and row.family != family:
            angle -= family_gap
        family = row.family
        start, end = angle, angle - span
        centre = np.deg2rad((start + end) / 2)
        centres[row.key] = centre
        ax.add_patch(Wedge((0, 0), 1.0, end, start, width=0.045, facecolor=config["family_colors"][row.family], edgecolor="white", linewidth=0.3))
        # The classes of the study platform share its sector in the proportion of their sample records.
        cursor = start
        for name in config["class_order"]:
            share = getattr(row, name) / row.records
            if share > 0:
                ax.add_patch(Wedge((0, 0), 0.945, cursor - span * share, cursor, width=0.045, facecolor=config["class_colors"][name], edgecolor="none"))
                cursor -= span * share
        length = 0.24 * np.log10(row.records) / largest
        ax.add_patch(Wedge((0, 0), 0.63 + length, end + 0.12 * span, start - 0.12 * span, width=length, facecolor="#AEB5C4", edgecolor="none"))
        for ring, name in enumerate(RINGS):
            value = getattr(row, f"share_{name}")
            ax.add_patch(Wedge((0, 0), 0.6 - ring * 0.036, end, start, width=0.03, facecolor=shares(value) if value > 0 else V.MISSING, edgecolor="white",
                               linewidth=0.2))
        turn = np.rad2deg(centre) % 360
        flipped = 90 < turn < 270
        ax.text(1.03 * np.cos(centre), 1.03 * np.sin(centre), row.label, rotation=turn + 180 if flipped else turn, rotation_mode="anchor",
                ha="right" if flipped else "left", va="center", fontsize=8, color=V.INK)
        angle = end
    # The logarithmic scale of the sample records, written in the gap at the top.
    for count in (10, 100):
        radius = 0.63 + 0.24 * np.log10(count) / largest
        ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=90 + top_gap / 2, theta2=450 - top_gap / 2, edgecolor="#BBBBBB", linewidth=0.6,
                         linestyle=":"))
        ax.text(0, radius, str(count), ha="center", va="center", fontsize=8, color=V.MUTED)
    most = shared["samples"].max()
    inner = 0.6 - (len(RINGS) - 1) * 0.036 - 0.045
    for row in shared.itertuples(index=False):
        first, second = centres[row.first], centres[row.second]
        path = Path([(inner * np.cos(first), inner * np.sin(first)), (0.15 * np.cos((first + second) / 2), 0.15 * np.sin((first + second) / 2)),
                     (inner * np.cos(second), inner * np.sin(second))], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        ax.add_patch(PathPatch(path, facecolor="none", edgecolor=group_colors[row.cohort_group], linewidth=1.0 + 5.5 * row.samples / most, alpha=0.8,
                               capstyle="round"))
    bar = ax.inset_axes([0.74, 0.035, 0.24, 0.022])
    bar.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=shares, extent=(0, 1, 0, 1))
    bar.set_yticks([])
    bar.set_xticks([0, 0.5, 1], ["0%", "50%", "100%"])
    bar.tick_params(length=2, pad=1.5)
    for side in bar.spines.values():
        side.set_visible(False)
    bar.set_title("Share recorded", fontsize=8, color=V.INK, pad=3)
    handles = [Rectangle((0, 0), 1, 1, color=config["family_colors"][name]) for name in config["family_order"]]
    handles += [Rectangle((0, 0), 1, 1, color=config["class_colors"][name]) for name in config["class_order"]]
    canvas.legend(handles, [config["family_labels"][name] for name in config["family_order"]] + [config["class_labels"][name] for name in config["class_order"]],
                  loc="outside lower center", ncols=5, columnspacing=1.2, handlelength=1.2)


def class_table(canvas, table, config, letter, rows, ceiling):
    """One half of the study platforms: the sample records of every curated class as a shaded cell with its count; both halves share their scale."""
    ax = canvas.subplots()
    classes = config["class_order"]
    image = np.ones((rows, len(classes) + 1, 3))
    for column, name in enumerate([*classes, "records"]):
        shade = V.ramp(config["class_colors"][name] if name != "records" else "#8C96AA")
        for position, value in enumerate(table[name].to_numpy()):
            level = np.log1p(value) / np.log1p(ceiling)
            image[position, column] = to_rgb(shade(level)) if value > 0 else (1.0, 1.0, 1.0)
            if value > 0:
                red, green, blue = image[position, column]
                ink = "white" if 0.299 * red + 0.587 * green + 0.114 * blue < 0.55 else V.INK
                ax.text(column, position, V.number(value), ha="center", va="center", fontsize=8, color=ink)
    ax.imshow(image, aspect="auto", interpolation="nearest", extent=(-0.5, len(classes) + 0.5, rows - 0.5, -0.5))
    for position, row in enumerate(table.itertuples(index=False)):
        ax.add_patch(Rectangle((-1.05, position - 0.5), 0.4, 1, color=config["family_colors"][row.family], linewidth=0))
    ax.set_xlim(-1.1, len(classes) + 0.5)
    ax.set_ylim(rows - 0.5, -0.5)
    ax.set_yticks(range(len(table)), table["label"])
    ax.set_xticks(range(len(classes) + 1), [config["class_labels"][name] for name in classes] + ["All records"], rotation=35, ha="left",
                  rotation_mode="anchor")
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    for column in range(len(classes) + 1):
        ax.axvline(column + 0.5, color="white", linewidth=1.2)
    V.panel(canvas, ax, letter, "")
    if letter == "A":
        handles = [Rectangle((0, 0), 1, 1, color=config["family_colors"][name]) for name in config["family_order"]]
        canvas.legend(handles, [config["family_labels"][name] for name in config["family_order"]], loc="outside lower center", ncols=3,
                      title="Technology", title_fontsize=8)


def main():
    run = V.Run(1)
    config = run.config
    census = run.table("acquisition/census.parquet", ["source", "census_status"])
    decisions = run.table("acquisition/study_decisions.parquet", ["source", "accession", "decision", "reason"])
    studies = run.table("acquisition/included_studies.parquet", ["source", "accession", "assay"])
    samples = run.table("cohort/cohort_samples.parquet")
    ledger = run.table("cohort/records_without_values.parquet", ["study_accession", "platform", "reason"])
    records = run.receipt("outputs/cohort_curation/01.json")["records"]
    summary = run.table("cohort/cohort_summary.parquet")
    families = run.table("annotation/matrix_index.parquet", ["accession", "platform", "family"]).rename(columns={"accession": "study_accession"})
    keys = lambda frame: set(zip(frame["study_accession"], frame["platform"]))
    V.require(records["included"] == len(samples) + len(ledger) and records["with_values"] == len(samples)
              and records["without_values"]["by_reason"] == {str(name): int(count) for name, count in ledger["reason"].value_counts().items()},
              "The sample records with and without expression values must add up to the included records of the cohort receipt")
    V.require(keys(families) == keys(summary) | keys(ledger) and keys(summary) == keys(samples),
              "Every included study platform must hold sample records, and the cohort summary those with expression values")
    summary = summary.merge(families, on=["study_accession", "platform"], validate="one_to_one")
    V.require(set(summary["family"]) <= set(config["family_order"]), "A technology family has no colour in the configuration")
    V.require(set(census["source"]) | set(studies["source"]) | set(decisions["source"]) <= set(REPOSITORIES), "A repository has no name in the figure code")
    V.require(set(samples["sample_class"]) <= set(config["class_order"]), "A curated class has no place in the configuration")
    distinct = samples[samples["is_primary_record"]]

    # The flow of panel B: the included records, those without values by reason, the records with values and their repeat deposits, the distinct
    # samples, one sample per patient and the classes.
    flow = pd.DataFrame([("records_found", len(census)), ("candidate_studies", len(decisions)), ("included_studies", len(studies)),
                         ("included_study_platforms", len(families)), ("study_platforms_with_values", len(summary)),
                         ("included_sample_records", records["included"]), ("records_without_values", len(ledger)),
                         ("sample_records_with_values", len(samples)), ("repeat_deposits", len(samples) - len(distinct)), ("distinct_samples", len(distinct)),
                         ("one_sample_per_patient", int(samples["one_per_patient"].sum())),
                         ("identified_patients", int(distinct.loc[distinct["patient_identified"], "patient_key"].nunique())),
                         ("samples_without_identifier", int((~distinct["patient_identified"]).sum())),
                         ("sepsis_before_onset", int((distinct["sepsis_phase"] == "before_onset").sum())),
                         *((f"class_{name}", int((distinct["sample_class"] == name).sum())) for name in config["class_order"]),
                         *((f"not_candidate_{name}", count) for name, count in census.loc[census["census_status"] != "candidate", "census_status"].value_counts().items()),
                         *((f"excluded_{name}", count) for name, count in decisions.loc[decisions["decision"] != "include", "reason"].value_counts().items()),
                         *((f"without_values_{name}", count) for name, count in ledger["reason"].value_counts().items())], columns=["step", "count"])
    run.add("Fig1", "flow", flow)
    run.figure("Fig1", (14, 8.4), dict(nrows=1, ncols=2), [
        ("studies", (0, 0), (7, 8.4), lambda canvas: flow_studies(canvas, census, decisions, studies, len(families), len(summary))),
        ("samples", (0, 1), (7, 8.4), lambda canvas: flow_samples(canvas, records["included"], ledger, samples, config))])

    summary["year"] = summary["release_date"].str[:4].astype(int)
    summary["key"] = summary["study_accession"] + "|" + summary["platform"]
    keyed = distinct.merge(summary[["study_accession", "platform", "year", "family"]], on=["study_accession", "platform"], how="left", validate="many_to_one")
    V.require(keyed["family"].notna().all(), "Every distinct sample must belong to a study platform of the cohort summary")
    years = list(range(int(summary["year"].min()), int(summary["year"].max()) + 1))
    by_year = keyed.groupby(["year", "family"]).size().unstack(fill_value=0).reindex(index=years, fill_value=0)
    run.add("Fig2", "samples_by_year", by_year.stack().rename("samples").reset_index().assign(
        added_up=lambda frame: frame.groupby("family")["samples"].cumsum()))
    material = {tissue: group for group, tissues in config["tissue_groups"].items() for tissue in tissues}
    V.require(set(distinct["tissue"]) <= set(material), "A material has no group in the configuration")
    keyed["material"] = keyed["tissue"].map(material)
    material_order = [name for name in config["tissue_groups"] if (keyed["material"] == name).any()]
    first = keyed.groupby(["family", "material"]).size().rename("value").reset_index().rename(columns={"family": "source", "material": "target"})
    second = keyed.groupby(["material", "sample_class"]).size().rename("value").reset_index().rename(columns={"material": "source", "sample_class": "target"})
    run.add("Fig2", "flows", pd.concat([first.assign(step="technology_to_material"), second.assign(step="material_to_class")], ignore_index=True))
    composition = pd.crosstab(keyed["material"], keyed["sample_class"]).reindex(index=list(config["tissue_groups"]), columns=config["class_order"], fill_value=0)
    run.add("Fig2", "composition_material", composition.rename_axis("group").reset_index())

    variables = {"sex": samples["sex"].isin(["female", "male"]), "age": samples["age_group"] != "not_recorded", "patient": samples["patient_identified"],
                "time_point": samples["time_order"].notna(), "outcome": samples["outcome"] != "not_recorded", "septic_shock": samples["septic_shock"] != "not_recorded"}
    V.require(list(variables) == list(RINGS), "Every ring of the atlas needs its variable")
    recorded = pd.DataFrame(variables).groupby([samples["study_accession"], samples["platform"]]).mean().add_prefix("share_").reset_index()
    order = {name: position for position, name in enumerate(config["family_order"])}
    platforms = summary.merge(recorded, on=["study_accession", "platform"], validate="one_to_one")
    platforms = platforms.assign(order=platforms["family"].map(order)).sort_values(["order", "release_date", "study_accession", "platform"]).reset_index(drop=True)
    # A study with several platforms is named with its platform, so that its sectors can be told apart.
    platforms["label"] = platforms["study_accession"] + np.where(platforms["study_accession"].duplicated(keep=False), " " + platforms["platform"], "")
    V.require((platforms[config["class_order"]].sum(axis=1) == platforms["records"]).all(), "The classes of a study platform must add up to its records")
    repeated = samples[samples["records_of_sample"] > 1].assign(key=lambda frame: frame["study_accession"] + "|" + frame["platform"])
    pairs = pd.DataFrame([pair for keys in repeated.groupby("sample_id")["key"].agg(lambda values: sorted(set(values)))
                          for pair in itertools.combinations(keys, 2)], columns=["first", "second"])
    shared = pairs.groupby(["first", "second"]).size().rename("samples").reset_index()
    groups = summary.set_index("key")["cohort_group"]
    V.require((shared["first"].map(groups) == shared["second"].map(groups)).all(), "Two study platforms that share samples must be in one cohort group")
    shared["cohort_group"] = shared["first"].map(groups)
    shared = shared.sort_values(["cohort_group", "first", "second"]).reset_index(drop=True)
    run.add("Fig2", "atlas", platforms[["label", "family", "year", "records", *config["class_order"], *(f"share_{name}" for name in RINGS)]].round(6))
    run.add("Fig2", "shared_samples", shared[["first", "second", "cohort_group", "samples"]])
    run.figure("Fig2", (14, 10.4), dict(nrows=2, ncols=2, width_ratios=[5.4, 8.6], height_ratios=[3.9, 6.5]), [
        ("samples_by_year", (0, 0), (5.4, 3.9), lambda canvas: growth(canvas, by_year, config)),
        ("flows", (1, 0), (5.4, 6.5), lambda canvas: material_flow(canvas, [first, second], config, material_order)),
        ("atlas", (slice(None), 1), (8.6, 10.4), lambda canvas: atlas(canvas, platforms, shared, config, int(platforms["records"].max())))])

    table = summary.sort_values(["study_accession", "platform"]).reset_index(drop=True)
    table["label"] = table["study_accession"] + np.where(table["platform"] != "", "  " + table["platform"], "")
    run.add("SuppFig1", "sample_records", table[["label", *config["class_order"], "records"]])
    half = (len(table) + 1) // 2
    ceiling = float(table["records"].max())
    run.figure("SuppFig1", (13, 11), dict(nrows=1, ncols=2), [
        ("first_half", (0, 0), (6.5, 11), lambda canvas: class_table(canvas, table.iloc[:half], config, "A", half, ceiling)),
        ("second_half", (0, 1), (6.5, 11), lambda canvas: class_table(canvas, table.iloc[half:], config, "B", half, ceiling))])
    run.finish()


if __name__ == "__main__":
    main()
