"""Draw the quality measures and flags of the samples: library sizes, dynamic range of the arrays and the sample records that carry each flag: SuppFig5."""
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter1d

import visualization_common as V

BELOW, ABOVE = "#DB4325", "#457B9D"
STRIP_BINS = 300


def strips(canvas, table, column, limit, flag, letter, title, label, ticks):
    """One strip for every study platform, ordered by its median: the density of the measure over its samples on a logarithmic axis, shaded from white at
    no samples to full colour at the densest bin of the strip, red below the limit of the cohort stage; a tick marks the median, and the samples that
    the limit flags are counted at the right."""
    ax = canvas.subplots()
    order = list(table.groupby("label")[column].median().sort_values(ascending=False).index)
    logged = np.log10(table[column])
    edges = np.linspace(float(logged.min()) - 0.15, float(logged.max()) + 0.15, STRIP_BINS + 1)
    centres = (edges[:-1] + edges[1:]) / 2
    cut = np.log10(limit)
    shades = {True: V.ramp(BELOW), False: V.ramp(ABOVE)}
    image = np.ones((len(order), STRIP_BINS, 3))
    medians = []
    for row, name in enumerate(order):
        values = logged[table["label"] == name].to_numpy()
        counts, _ = np.histogram(values, bins=edges)
        smooth = gaussian_filter1d(counts.astype(float), 3.0, mode="constant")
        # A root of the density, so that the tails of a study stay visible beside its peak.
        level = np.sqrt(smooth / smooth.max())
        for position, (centre, value) in enumerate(zip(centres, level)):
            image[row, position] = to_rgb(shades[bool(centre < cut)](value))
        medians.append(float(np.median(values)))
        flagged = int(table.loc[table["label"] == name, flag].sum())
        if flagged:
            ax.text(1.01, row, f"{flagged} flagged", transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=BELOW)
    ax.imshow(image, aspect="auto", interpolation="nearest", extent=(edges[0], edges[-1], len(order) - 0.5, -0.5))
    for row in range(1, len(order)):
        ax.axhline(row - 0.5, color="white", linewidth=0.8)
    ax.vlines(medians, np.arange(len(order)) - 0.38, np.arange(len(order)) + 0.38, color=V.INK, linewidth=1.1)
    ax.axvline(cut, color=V.INK, linestyle="--", linewidth=0.8)
    ax.set_yticks(range(len(order)), order)
    ax.set_xticks(np.log10(ticks), [f"{value:g}" if value < 1e3 else f"10$^{{{int(np.log10(value))}}}$" for value in ticks])
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlabel(label)
    ax.tick_params(axis="y", length=0)
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    V.panel(canvas, ax, letter, title)
    canvas.legend([Patch(color=ABOVE), Patch(color=BELOW), Line2D([], [], color=V.INK, linestyle="--", linewidth=0.8),
                   Line2D([], [], color=V.INK, marker="|", linestyle="", markersize=8, markeredgewidth=1.1)],
                  ["Above the limit", "Below the limit (flagged)", "Limit", "Median"], loc="outside lower center", ncols=4)


def flag_map(canvas, counts, labels, config):
    """The sample records that carry each flag, by technology: a shaded cell with its count, and the records of each flag in all."""
    ax = canvas.subplots()
    families = [name for name in config["family_order"] if name in counts.columns]
    values = counts[families].to_numpy(dtype=float)
    shades = V.cmap()
    norm = LogNorm(1, max(float(values.max()), 2.0))
    ax.imshow(np.ma.masked_equal(values, 0), aspect="auto", cmap=shades, norm=norm, interpolation="nearest")
    for (row, column), value in np.ndenumerate(values):
        if value > 0:
            red, green, blue, _ = shades(norm(value))
            ax.text(column, row, V.number(value), ha="center", va="center", fontsize=8, color="white" if 0.299 * red + 0.587 * green + 0.114 * blue < 0.55 else V.INK)
    for row, total in enumerate(counts[families].sum(axis=1)):
        ax.text(len(families) - 0.35, row, V.number(total), ha="left", va="center", fontsize=8, color=V.INK)
    for column in range(1, len(families)):
        ax.axvline(column - 0.5, color="white", linewidth=1.5)
    for row in range(1, len(counts)):
        ax.axhline(row - 0.5, color="white", linewidth=1.5)
    ax.set_facecolor("#F7F7F7")
    ax.set_xlim(-0.5, len(families) + 0.3)
    ax.set_yticks(range(len(counts)), labels)
    ax.set_xticks(range(len(families)), [config["family_labels"][name].replace(" ", "\n", 1) for name in families])
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    V.panel(canvas, ax, "C", "Sample records that carry each flag")


def main():
    run = V.Run(5)
    config = run.config
    quality = run.table("cohort/sample_quality.parquet")
    families = run.table("annotation/matrix_index.parquet", ["accession", "platform", "family"]).rename(columns={"accession": "study_accession"})
    limits = V.read_json(V.ROOT / "src/04_cohort_curation/cohort_config.json")["quality"]
    quality["label"] = quality["study_accession"] + np.where(quality["platform"] != "", "  " + quality["platform"], "")
    libraries = quality[quality["library_size"].notna()]
    arrays = quality[quality["dynamic_range_to_study_median"].notna()]
    V.require((libraries["library_size"] > 0).all() and (arrays["dynamic_range_to_study_median"] > 0).all(), "A measure on a logarithmic axis must be positive")
    run.add("SuppFig5", "library_sizes", libraries[["label", "sample_accession", "library_size", "low_library_size"]])
    run.add("SuppFig5", "dynamic_range", arrays[["label", "sample_accession", "dynamic_range_to_study_median", "low_dynamic_range"]].round(6))
    counts = pd.DataFrame({"records": {name: int(quality[name].sum()) for name in config["flag_labels"]}})
    V.require((quality[list(config["flag_labels"])].any(axis=1) == quality["any_flag"]).all(), "The flags of the configuration must be exactly those that flag a sample")
    # A label may name a limit of the cohort stage, which is read from its configuration rather than written twice.
    labels = [config["flag_labels"][name].format(**limits) for name in counts.index]
    run.add("SuppFig5", "flags", counts.rename_axis("flag").reset_index())
    keyed = quality.merge(families, on=["study_accession", "platform"], how="left", validate="many_to_one")
    V.require(keyed["family"].notna().all(), "Every sample record must belong to a study platform of the matrix index")
    by_family = pd.DataFrame({name: keyed.loc[keyed[name], "family"].value_counts() for name in config["flag_labels"]}).T.fillna(0).astype(int)
    by_family = by_family.reindex(index=list(config["flag_labels"]), columns=[name for name in config["family_order"] if name in by_family.columns], fill_value=0)
    run.add("SuppFig5", "flags_by_technology", by_family.rename_axis("flag").reset_index())
    run.figure("SuppFig5", (14, 7.6), dict(nrows=2, ncols=2, width_ratios=[1.0, 1.0], height_ratios=[1.0, 1.0]), [
        ("library_sizes", (slice(None), 0), (7, 7.6), lambda canvas: strips(canvas, libraries, "library_size", limits["minimum_library_size"], "low_library_size",
                                                                            "A", "Reads assigned to genes, every RNA-seq library", "Assigned reads",
                                                                            [1e3, 1e4, 1e5, 1e6, 1e7, 1e8])),
        ("dynamic_range", (0, 1), (7, 3.8), lambda canvas: strips(canvas, arrays, "dynamic_range_to_study_median",
                                                                  limits["dynamic_range_below_share_of_study_median"], "low_dynamic_range", "B",
                                                                  "Dynamic range of the arrays read from intensities",
                                                                  "Dynamic range relative to the median of the study platform", [0.03, 0.1, 0.3, 1, 3])),
        ("flags", (1, 1), (7, 3.8), lambda canvas: flag_map(canvas, by_family, labels, config))])
    run.finish()


if __name__ == "__main__":
    main()
