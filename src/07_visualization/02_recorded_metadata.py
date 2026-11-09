"""Draw what the records say about the samples: what the class labels rest on, the ages, the recorded shock and outcome of the sepsis samples, when the
samples of the longitudinal studies were drawn, and how many sample records record each curated variable: Fig3, SuppFig2."""
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, LogNorm, Normalize
from matplotlib.patches import Rectangle, Wedge

import visualization_common as V

# The coverage panels: the colour of a share of one, and the gap between two technologies, in columns.
SHARE_COLOR, GAP = "#2F6B94", 0.6
BASIS_TINT = {"sample_record": 0.0, "sample_record_and_study_description": 0.4, "study_description": 0.7, "not_recorded": 0.88}
BASIS_HATCH = {"sample_record_and_study_description": "///", "study_description": "..."}
SHOCK = {"yes": ("Septic shock", "#E97233"), "no": ("No septic shock", "#FCE082"), "not_recorded": ("Shock not recorded", "#D0D0D0")}
OUTCOME = {"survivor": ("Survivor", "#4A6FE3"), "non_survivor": ("Non-survivor", "#DB4325"), "not_recorded": ("Outcome not recorded", "#D0D0D0")}
DAY_UNITS = {"day": 1.0, "hour": 1 / 24}


def columns(families):
    """The place of every study platform in a coverage panel, with a gap where the technology changes."""
    families = np.asarray(families)
    return np.arange(len(families)) + GAP * np.r_[0, np.cumsum(families[1:] != families[:-1])]


def coverage(canvas, table, shares, width, config, letter, title, legend=False):
    """The study platforms of some technologies as columns and the curated variables as rows: a cell shaded by the share of the sample records of its
    study platform that record its variable, and white where none does. Bars above give the sample records of every study platform in the colour of
    its technology, and the figure at the end of a row the share of all the sample records of the panel that record the variable."""
    top, ax = canvas.subplots(2, 1, height_ratios=[1.8, shares.shape[1]], sharex=True)
    xs = columns(table["family"])
    shades = ListedColormap(V.ramp(SHARE_COLOR)(np.linspace(0.15, 1, 256)))
    for (column, row), value in np.ndenumerate(shares.to_numpy()):
        if value > 0:
            ax.add_patch(Rectangle((xs[column] - 0.5, row - 0.5), 1, 1, facecolor=shades(value), edgecolor="white", linewidth=0.7))
    overall = shares.mul(table["records"], axis=0).sum() / table["records"].sum()
    for row, name in enumerate(shares.columns):
        ax.text(xs[-1] + 0.9, row, f"{overall[name]:.0%}", va="center", ha="left", fontsize=8)
    ax.set_xlim(-0.5, width - 0.5)
    ax.set_ylim(shares.shape[1] - 0.5, -0.5)
    ax.set_yticks(range(shares.shape[1]), [config["variables"][name] for name in shares.columns])
    ax.set_xticks(xs, table["label"], rotation=90)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    top.bar(xs, table["records"], width=0.78, color=[config["family_colors"][name] for name in table["family"]])
    top.set_yscale("log")
    top.set_ylim(1, 1000)
    top.set_yticks([10, 1000], ["10", "1,000"])
    top.minorticks_off()
    top.set_ylabel("Records", rotation=0, ha="right", va="center")
    top.tick_params(axis="x", length=0)
    for side in ("top", "right", "bottom"):
        top.spines[side].set_visible(False)
    V.panel(canvas, top, letter, title)
    if legend:
        # The space that the fewer study platforms of this panel leave at the right holds the legend of both panels.
        start = (xs[-1] + 4.5) / width
        handles = [Rectangle((0, 0), 1, 1, color=config["family_colors"][name]) for name in config["family_order"]]
        ax.legend(handles, [config["family_labels"][name] for name in config["family_order"]], loc="upper left", bbox_to_anchor=(start, 1.0))
        bar = canvas.colorbar(ScalarMappable(Normalize(0, 1), shades), cax=ax.inset_axes([start + 0.01, 0.1, 0.15, 0.07]), orientation="horizontal")
        bar.outline.set_visible(False)
        bar.set_ticks([0, 0.5, 1], labels=["0%", "50%", "100%"])
        bar.set_label("Share recorded")


def label_basis(canvas, counts, config):
    """Two rings: the curated class of every sample record inside, and outside what that class rests on."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0)
    ax.update_datalim([(-1.08, -1.08), (1.08, 1.08)])
    V.panel(canvas, ax, "A", "What a class rests on")
    total = float(counts.to_numpy().sum())
    angle = 90.0
    for name, row in counts.iterrows():
        span = 360 * row.sum() / total
        if span <= 0:
            continue
        ax.add_patch(Wedge((0, 0), 0.68, angle - span, angle, width=0.24, facecolor=config["class_colors"][name], edgecolor="white", linewidth=1.2))
        inner = angle
        for basis, value in row.items():
            part = 360 * value / total
            if part > 0:
                ax.add_patch(Wedge((0, 0), 1.0, inner - part, inner, width=0.27, facecolor=V.tint(config["class_colors"][name], BASIS_TINT[basis]), edgecolor="white",
                                   linewidth=0.8, hatch=BASIS_HATCH.get(basis, "")))
                inner -= part
        angle -= span
    ax.text(0, 0, V.number(total), ha="center", va="center", fontsize=12, color=V.INK)
    handles = [Rectangle((0, 0), 1, 1, color=config["class_colors"][name]) for name in counts.index]
    handles += [Rectangle((0, 0), 1, 1, facecolor=V.tint(V.MUTED, BASIS_TINT[basis]), edgecolor="white", hatch=BASIS_HATCH.get(basis, "")) for basis in counts.columns]
    labels = [f"{config['class_labels'][name]}  {V.number(counts.loc[name].sum())}" for name in counts.index]
    labels += [config["basis_labels"][basis] for basis in counts.columns]
    canvas.legend(handles, labels, loc="outside lower center", ncols=2, columnspacing=1.0, handlelength=1.2)


def ages(canvas, table, config, grid):
    """The recorded age in years of the distinct samples of every class, as a density reflected at zero, so that the newborns keep their weight."""
    ax = canvas.subplots()
    ax.axvspan(0, 18, color="#F5F5F5", linewidth=0, zorder=0)
    top = 0.0
    for name in config["class_order"]:
        values = table.loc[table["sample_class"] == name, "age_years"].to_numpy()
        density = V.kde(np.concatenate([values, -values]), grid)
        if density is None:
            continue
        density = 2 * density
        top = max(top, density.max())
        color = config["class_colors"][name]
        ax.fill_between(grid, density, color=color, alpha=0.14, linewidth=0)
        ax.plot(grid, density, color=color, linewidth=1.5, label=f"{config['class_labels'][name]}  ({V.number(len(values))})")
    ax.set_xlim(grid[0], grid[-1])
    # Room above the curves for the legend.
    ax.set_ylim(0, top * 1.55)
    ax.text(9, top * 1.48, "Under 18", ha="center", va="center", fontsize=8, color=V.INK)
    ax.set_xlabel("Recorded age (years)")
    ax.set_ylabel("Density")
    V.clean_axis(ax, grid="y")
    V.panel(canvas, ax, "B", "Ages of the samples")
    ax.legend(loc="upper right", title="Distinct samples with an age", title_fontsize=8)


def clinical_flow(canvas, frame, total):
    """The distinct sepsis samples, from their recorded septic shock to their recorded outcome."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_xlim(-0.5, 1.72)
    ax.set_ylim(-0.01, 1.01)
    first = frame.groupby("shock").size().rename("value").reset_index().assign(source="sepsis").rename(columns={"shock": "target"})
    second = frame.groupby(["shock", "outcome"]).size().rename("value").reset_index().rename(columns={"shock": "source", "outcome": "target"})
    second = second.assign(target="outcome:" + second["target"])
    shocks = [name for name in SHOCK if (frame["shock"] == name).any()]
    outcomes = [f"outcome:{name}" for name in OUTCOME if (frame["outcome"] == name).any()]
    colors = {"sepsis": "#DB4325", **{name: SHOCK[name][1] for name in SHOCK}, **{f"outcome:{name}": OUTCOME[name][1] for name in OUTCOME}}
    counts = {**frame["shock"].value_counts().to_dict(), **{f"outcome:{name}": count for name, count in frame["outcome"].value_counts().items()}}
    labels = {"sepsis": f"Sepsis\n{V.number(total)}", **{name: f"{SHOCK[name][0]}\n{V.number(counts[name])}" for name in shocks},
              **{key: f"{OUTCOME[key.split(':')[1]][0]}  {V.number(counts[key])}" for key in outcomes}}
    V.alluvial(ax, [["sepsis"], shocks, outcomes], [(first, "target"), (second, "target")], colors, labels, xs=[0.0, 0.5, 1.0], gap=0.05, inner="left")
    V.panel(canvas, ax, "C", "Shock and outcome of the sepsis samples")


def calendar(canvas, blocks, config):
    """When the distinct samples of every longitudinal study platform were drawn: one row for every study platform, one column for every day (left), or
    for every time point in its order where the records give no unit or only an order (right); the shade and the number give the samples."""
    axes = canvas.subplots(1, 2, width_ratios=[blocks[0].shape[1] + 4, blocks[1].shape[1] + 4])
    top = max(float(block.to_numpy().max()) for block in blocks)
    shades = V.cmap()
    norm = LogNorm(1, top)
    for ax, block, label in zip(axes, blocks, ("Day of sampling (hours converted to days)", "Time point in the order of the study")):
        values = block.to_numpy(dtype=float)
        image = ax.imshow(np.ma.masked_equal(values, 0), aspect="auto", cmap=shades, norm=norm, interpolation="nearest", rasterized=True)
        for (row, column), value in np.ndenumerate(values):
            if value > 0:
                red, green, blue, _ = shades(norm(value))
                ink = "white" if 0.299 * red + 0.587 * green + 0.114 * blue < 0.55 else V.INK
                ax.text(column, row, V.number(value), ha="center", va="center", fontsize=8, color=ink)
        ax.set_xticks(range(values.shape[1]), [str(name) for name in block.columns])
        ax.set_yticks(range(values.shape[0]), block.index)
        ax.set_xlim(-0.5, values.shape[1] - 0.5)
        ax.set_ylim(values.shape[0] - 0.5, -0.5)
        for row in range(values.shape[0] - 1):
            ax.axhline(row + 0.5, color="white", linewidth=1.0)
        ax.set_facecolor("#F7F7F7")
        ax.set_xlabel(label)
        ax.tick_params(length=0)
        for side in ax.spines.values():
            side.set_visible(False)
    V.panel(canvas, axes[0], "D", "When the samples of the longitudinal studies were drawn")
    bar = canvas.colorbar(image, ax=axes, location="right", shrink=0.8, pad=0.01, aspect=25)
    bar.outline.set_visible(False)
    bar.set_label("Distinct samples")


def main():
    run = V.Run(2)
    config = run.config
    display = config["display"]
    samples = run.table("cohort/cohort_samples.parquet")
    scores = run.table("cohort/sample_scores.parquet", ["study_accession", "sample_accession"]).drop_duplicates()
    families = run.table("annotation/matrix_index.parquet", ["accession", "platform", "family"]).rename(columns={"accession": "study_accession"})
    scored = set(zip(scores["study_accession"], scores["sample_accession"]))
    recorded = pd.DataFrame({
        "sample_class": samples["sample_class"] != "not_recorded", "tissue": samples["tissue"] != "not_recorded", "sex": samples["sex"].isin(["female", "male"]),
        "age_group": samples["age_group"] != "not_recorded", "age_years": samples["age_years"].notna(), "patient": samples["patient_identified"],
        "time_point": samples["time_order"].notna(), "outcome": samples["outcome"] != "not_recorded", "septic_shock": samples["septic_shock"] != "not_recorded",
        "score": [key in scored for key in zip(samples["study_accession"], samples["sample_accession"])]})
    V.require(list(recorded.columns) == list(config["variables"]), "The variables of the coverage panel must be those of the configuration")
    V.require(set(samples["sample_class"]) <= set(config["class_order"]) and set(samples["class_basis"]) <= set(config["basis_labels"])
              and set(families["family"]) <= set(config["family_order"]), "A class, a class basis or a technology has no place in the configuration")
    keys = [samples["study_accession"], samples["platform"]]
    shares = recorded.groupby(keys).mean()
    table = samples.groupby(["study_accession", "platform"]).size().rename("records").reset_index().merge(families, on=["study_accession", "platform"],
                                                                                                       validate="one_to_one")
    shares = shares.loc[list(zip(table["study_accession"], table["platform"]))].reset_index(drop=True)
    table["label"] = table["study_accession"] + np.where(table["platform"] != "", "  " + table["platform"], "")
    # The variables from the most recorded down, and within every technology the study platforms sorted by their shares row by row from the top, so
    # that the recorded cells fall into steps.
    rows = list(recorded.mean().sort_values(ascending=False, kind="stable").index)
    keys = pd.concat([table["family"].map({name: rank for rank, name in enumerate(config["family_order"])}), -shares[rows].round(2), -table["records"]], axis=1)
    order = keys.sort_values(list(keys.columns), kind="stable").index
    table, shares = table.loc[order].reset_index(drop=True), shares.loc[order].reset_index(drop=True)
    run.add("SuppFig2", "coverage", pd.concat([table[["label", "family", "records"]], shares.round(6)], axis=1))
    arrays = table["family"] != "rna_seq"
    width = max(float(columns(table.loc[part, "family"])[-1]) + 1 for part in (arrays, ~arrays))

    basis = pd.crosstab(samples["sample_class"], samples["class_basis"]).reindex(index=config["class_order"], columns=list(config["basis_labels"]), fill_value=0)
    run.add("Fig3", "label_basis", basis.rename_axis("sample_class").reset_index())
    distinct = samples[samples["is_primary_record"]]
    aged = distinct[distinct["age_years"].notna()][["sample_class", "age_years"]]
    V.require((aged["age_years"] >= 0).all(), "A recorded age must not be negative")
    grid = np.linspace(0, max(100.0, float(aged["age_years"].max())), 401)
    run.add("Fig3", "ages", aged.groupby("sample_class")["age_years"].describe().round(6).rename_axis("sample_class").reset_index())
    sepsis = distinct[distinct["sample_class"] == "sepsis"]
    V.require(set(sepsis["septic_shock"]) <= set(SHOCK) and set(sepsis["outcome"]) <= set(OUTCOME), "A recorded shock or outcome has no place in the figure code")
    clinical = pd.DataFrame({"shock": sepsis["septic_shock"].to_numpy(), "outcome": sepsis["outcome"].to_numpy()})
    run.add("Fig3", "clinical_flow", clinical.value_counts().rename("samples").reset_index())

    timed = distinct[distinct["time_order"].notna()].copy()
    timed["label"] = timed["study_accession"] + np.where(timed.groupby("study_accession")["platform"].transform("nunique") > 1, " " + timed["platform"], "")
    in_days = timed[timed["time_unit"].isin(DAY_UNITS)].copy()
    in_days["day"] = (in_days["time_order"] * in_days["time_unit"].map(DAY_UNITS)).round().astype(int)
    V.require(in_days["day"].between(0, display["calendar_days"]).all(), "A sampling day lies beyond the days of the calendar")
    ordered = timed[~timed["time_unit"].isin(DAY_UNITS)].copy()
    ordered["step"] = ordered.groupby("label")["time_order"].rank(method="dense").astype(int)
    days = pd.crosstab(in_days["label"], in_days["day"]).reindex(columns=range(display["calendar_days"] + 1), fill_value=0)
    steps = pd.crosstab(ordered["label"], ordered["step"])
    # Rows from the earliest and densest schedules down.
    days = days.loc[sorted(days.index, key=lambda name: (int(np.flatnonzero(days.loc[name].to_numpy())[0]), -int((days.loc[name] > 0).sum()), name))]
    steps = steps.loc[sorted(steps.index, key=lambda name: (-int((steps.loc[name] > 0).sum()), name))]
    run.add("Fig3", "calendar", pd.concat([days.stack().rename("samples").reset_index().rename(columns={"label": "study", "day": "position"}).assign(scale="day"),
                                           steps.stack().rename("samples").reset_index().rename(columns={"label": "study", "step": "position"}).assign(scale="order")],
                                          ignore_index=True).query("samples > 0"))
    run.figure("Fig3", (14, 10.2), dict(nrows=2, ncols=3, width_ratios=[4.2, 4.9, 4.9], height_ratios=[4.7, 5.5]), [
        ("label_basis", (0, 0), (4.2, 4.7), lambda canvas: label_basis(canvas, basis, config)),
        ("ages", (0, 1), (4.9, 4.7), lambda canvas: ages(canvas, aged, config, grid)),
        ("clinical_flow", (0, 2), (4.9, 4.7), lambda canvas: clinical_flow(canvas, clinical, len(sepsis))),
        ("calendar", (1, slice(None)), (14, 5.5), lambda canvas: calendar(canvas, [days, steps], config))])
    run.figure("SuppFig2", (14, 8.4), dict(nrows=2, ncols=1), [
        ("arrays", (0, 0), (14, 4.2), lambda canvas: coverage(canvas, table[arrays].reset_index(drop=True), shares.loc[arrays, rows].reset_index(drop=True),
                                                               width, config, "A", "Microarray study platforms")),
        ("rna_seq", (1, 0), (14, 4.2), lambda canvas: coverage(canvas, table[~arrays].reset_index(drop=True), shares.loc[~arrays, rows].reset_index(drop=True),
                                                                width, config, "B", "RNA-seq study platforms", legend=True))])
    run.finish()


if __name__ == "__main__":
    main()
