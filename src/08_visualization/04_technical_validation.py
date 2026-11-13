"""Draw the technical validation: sex read from expression, agreement of repeated deposits, the Sepsis MetaScore, the published sepsis genes and the
markers of sorted cells: Fig5, SuppFig3, SuppFig4."""
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy.cluster import hierarchy

import visualization_common as V

# Sex has colours of its own, so that red and blue keep meaning sepsis and healthy controls in every figure.
SEX_COLORS = {"female": "#8E5BA8", "male": "#118AB2", "not_recorded": "#A0A0A0", "ambiguous": "#F3B032"}
SEX_LABELS = {"female": "Recorded female", "male": "Recorded male", "not_recorded": "Sex not recorded"}
MARK, CLOUD, SCORE, SCORE_DARK = "#E97233", "#457B9D", "#DB4325", "#8B1E0F"
TISSUE_LABELS = {"whole_blood": "whole blood", "blood_not_further_specified": "blood", "leukocytes": "leukocytes", "pbmc": "PBMC", "monocytes": "monocytes"}
TISSUE_COLORS = {"whole_blood": "#355070", "blood_not_further_specified": "#6D597A", "leukocytes": "#B56576", "pbmc": "#E56B6F", "monocytes": "#EAAC8B"}
CELL_LABELS = {"neutrophils": "Neutrophils", "monocytes": "Monocytes", "lymphocytes": "Lymphocytes", "b_cells": "B cells", "cd4_t_cells": "CD4 T cells",
               "cd8_t_cells": "CD8 T cells", "nk_cells": "NK cells", "pbmc": "PBMC"}
PLANE = (-0.7, 1.7)


def name(study, platform):
    return f"{study}  {platform}" if platform else study


def sex_plane(canvas, points, agree, compared):
    """Every called sample of the tables that measure XIST, placed by its own table: XIST from the mean of the samples read as male (0) to that of the
    samples read as female (1), and the Y-chromosome genes the other way; marginal densities by recorded sex, in proportion to their samples."""
    axes = canvas.subplots(2, 2, width_ratios=[4.2, 1.0], height_ratios=[1.0, 4.2])
    top, corner, main, side = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    corner.set_axis_off()
    shown = points.iloc[np.random.default_rng(20261013).permutation(len(points))]
    V.bubbles(main, shown["x"], shown["y"], shown["recorded"].map(SEX_COLORS), size=11, alpha=0.5, rasterized=True)
    against = points[points["disagrees"]]
    main.scatter(against["x"], against["y"], s=20, c=against["recorded"].map(SEX_COLORS), edgecolors="#111111", linewidths=0.6, alpha=0.95, zorder=3)
    grid = np.linspace(*PLANE, 400)
    for group in ("not_recorded", "male", "female"):
        part = points[points["recorded"] == group]
        for ax, column in ((top, "x"), (side, "y")):
            density = V.kde(part[column], grid)
            if density is None:
                continue
            density = density * len(part)
            if ax is top:
                ax.fill_between(grid, density, color=SEX_COLORS[group], alpha=0.22, linewidth=0)
                ax.plot(grid, density, color=SEX_COLORS[group], linewidth=1.2)
            else:
                ax.fill_betweenx(grid, density, color=SEX_COLORS[group], alpha=0.22, linewidth=0)
                ax.plot(density, grid, color=SEX_COLORS[group], linewidth=1.2)
    for ax in (main, top):
        ax.set_xlim(*PLANE)
    for ax in (main, side):
        ax.set_ylim(*PLANE)
    for value in (0, 1):
        main.axvline(value, color=V.GRID, linewidth=0.8, linestyle=":", zorder=0)
        main.axhline(value, color=V.GRID, linewidth=0.8, linestyle=":", zorder=0)
    top.set_ylim(bottom=0)
    side.set_xlim(left=0)
    top.set_yticks([])
    top.tick_params(axis="x", labelbottom=False)
    side.set_xticks([])
    side.tick_params(axis="y", labelleft=False)
    top.spines["left"].set_visible(False)
    side.spines["bottom"].set_visible(False)
    main.set_xlabel("XIST, from the male (0) to the female (1)\ngroup mean of its own table")
    main.set_ylabel("Y-chromosome genes, from the female (0)\nto the male (1) group mean of its own table")
    main.text(0.97, 0.97, f"{V.number(agree)} of {V.number(compared)} recorded\nsexes agree with expression", transform=main.transAxes, ha="right",
              va="top", fontsize=8, color=V.INK, linespacing=1.3)
    V.panel(canvas, top, "A", "Sex read from expression")
    counts = points["recorded"].value_counts()
    handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor=SEX_COLORS[group], markeredgecolor="white", markersize=7) for group in SEX_LABELS]
    handles.append(Line2D([], [], marker="o", linestyle="", markerfacecolor=SEX_COLORS["female"], markeredgecolor="#111111", markersize=6))
    canvas.legend(handles, [f"{SEX_LABELS[group]} ({V.number(counts.get(group, 0))})" for group in SEX_LABELS]
                  + [f"Record disagrees with expression ({V.number(len(against))})"], loc="outside lower center", ncols=2, columnspacing=1.0)


def repeats(canvas, pairs, kinds, config):
    """The two records of one sample, with every gene centred on its median in its own matrix, so that only what the sample adds to the average blood
    profile counts; a raw correlation would be near one for any two blood samples. One raincloud for every kind of evidence that links the records:
    the density above, the median and middle half in the middle, the pairs below."""
    ax = canvas.subplots()
    order = [kind for kind in [*config["link_labels"], "several"] if (kinds == kind).any()]
    labels = {**config["link_labels"], "several": "Several kinds\nof evidence"}
    jitter = np.random.default_rng(config["seed"])
    low = float(np.floor((pairs["centred_correlation"].min() - 0.05) * 10) / 10)
    grid = np.linspace(low, 1.0, 400)
    for position, kind in enumerate(order):
        part = pairs[kinds == kind]
        values = part["centred_correlation"].to_numpy()
        hit = part["most_similar"].to_numpy(dtype=bool)
        density = V.kde(values, grid)
        if density is not None:
            density = density / density.max() * 0.42
            # The cloud is drawn where it rises, not as a line along its whole base.
            density[density < 0.004] = np.nan
            ax.fill_between(grid, position - density, position, color=CLOUD, alpha=0.3, linewidth=0)
            ax.plot(grid, position - density, color=CLOUD, linewidth=1.0)
        first, middle, third = np.quantile(values, [0.25, 0.5, 0.75])
        ax.plot([first, third], [position + 0.07] * 2, color=V.INK, linewidth=2.4, solid_capstyle="butt")
        ax.plot([middle], [position + 0.07], marker="o", markersize=4.5, markerfacecolor="white", markeredgecolor=V.INK, markeredgewidth=1.0)
        spread = position + 0.16 + jitter.uniform(0, 0.24, len(values))
        V.bubbles(ax, values[hit], spread[hit], CLOUD, size=9, alpha=0.55, rasterized=True)
        ax.scatter(values[~hit], spread[~hit], s=18, marker="x", color=MARK, linewidth=1.0, zorder=3)
        ax.text(1.02, position, f"n={V.number(len(values))}", transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=V.INK)
    ax.set_yticks(range(len(order)), [labels[kind] for kind in order])
    ax.set_ylim(len(order) - 0.5, -0.62)
    ax.set_xlim(low, 1.01)
    ax.set_xlabel("Correlation of the two records (every gene\ncentred on its median in its own table)")
    V.clean_axis(ax, grid="x", categorical="y")
    V.panel(canvas, ax, "B", "Two deposits of one sample agree")
    canvas.legend([Line2D([], [], color=V.INK, linewidth=2.4), Line2D([], [], marker="x", linestyle="", color=MARK, markersize=5)],
                  ["Median and middle half", "Not the most similar record in the other table"], loc="outside lower center", ncols=1)


def roc_curve(cases, others):
    """The points of the ROC curve of a score that is higher in the cases, one point for every distinct value, so that ties move both rates at once."""
    values = np.concatenate([cases, others])
    labels = np.concatenate([np.ones(len(cases)), np.zeros(len(others))])
    order = np.argsort(-values, kind="mergesort")
    values, labels = values[order], labels[order]
    ends = np.r_[np.flatnonzero(np.diff(values)), len(values) - 1]
    hits = np.cumsum(labels)[ends]
    return np.r_[0.0, (ends + 1 - hits) / len(others)], np.r_[0.0, hits / len(cases)]


def roc_curves(canvas, curves, table):
    """The ROC curve of the Sepsis MetaScore between sepsis and healthy controls in every study platform and tissue, with the median curve and the middle
    half of the curves at every false-positive rate."""
    ax = canvas.subplots()
    grid = np.linspace(0, 1, 201)
    rows = []
    for _, part in curves.groupby("stratum", sort=False):
        ax.plot(part["fpr"], part["tpr"], color=SCORE, alpha=0.32, linewidth=0.9)
        rows.append(np.interp(grid, part["fpr"], part["tpr"]))
    rows = np.array(rows)
    middle = np.median(rows, axis=0)
    middle[0] = 0.0
    first, third = np.quantile(rows, [0.25, 0.75], axis=0)
    ax.fill_between(grid, first, third, color=SCORE, alpha=0.18, linewidth=0)
    ax.plot(grid, middle, color=SCORE_DARK, linewidth=2.4)
    ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1.1, zorder=0)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("False-positive rate (1 - specificity)")
    ax.set_ylabel("True-positive rate (sensitivity)")
    V.clean_axis(ax, grid="both")
    V.panel(canvas, ax, "C", "Sepsis MetaScore in single studies")
    handles = [Line2D([], [], color=SCORE, alpha=0.5, linewidth=0.9), Line2D([], [], color=SCORE_DARK, linewidth=2.4), Patch(color=SCORE, alpha=0.25)]
    ax.legend(handles, [f"One study platform and tissue ({len(table)})", "Median curve", "Middle half of the curves"], loc="lower right",
              title=f"AUROC median {table['auroc'].median():.3f}, lowest {table['auroc'].min():.2f}", title_fontsize=8)


def directions(canvas, matrix, roles, tissues, limit):
    """The standardised difference between sepsis and healthy controls of every published gene (rows) in every study platform and tissue (columns), the
    columns ordered by a clustering tree; left strip: the published direction; top strip: the tissue."""
    tree_ax, strip, ax = canvas.subplots(3, 1, height_ratios=[0.8, 0.16, 4.0])
    values = matrix.to_numpy(dtype=float)
    link = hierarchy.linkage(np.nan_to_num(values.T), "average", optimal_ordering=True)
    leaves = V.tree(tree_ax, link)
    values, columns = values[:, leaves], [matrix.columns[index] for index in leaves]
    colors = V.cmap(diverging=True).copy()
    colors.set_bad(V.MISSING)
    image = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=colors, vmin=-limit, vmax=limit, interpolation="nearest")
    strip.imshow(np.array([[to_rgb(TISSUE_COLORS[tissues[column]]) for column in columns]]), aspect="auto", interpolation="nearest")
    for position, role in enumerate(roles):
        ax.add_patch(Rectangle((-1.15, position - 0.5), 0.45, 1, color=SCORE if role == "higher_in_sepsis" else "#4A6FE3", linewidth=0, clip_on=False))
    for target in (tree_ax, strip, ax):
        target.set_xlim(-0.5, len(columns) - 0.5)
    # The tree and the strip keep their axes, unseen, so that they count as parts of the heatmap they stand over rather than as diagrams of their own.
    for target in (tree_ax, strip):
        target.set_xticks([])
        target.set_yticks([])
        for side in target.spines.values():
            side.set_visible(False)
    ax.set_yticks(range(len(matrix)), matrix.index)
    ax.set_xticks(range(len(columns)), [column[0] for column in columns], rotation=90)
    ax.tick_params(length=0)
    ax.tick_params(axis="y", pad=12)
    for side in ax.spines.values():
        side.set_visible(False)
    V.panel(canvas, tree_ax, "D", "Sepsis genes keep their direction across studies")
    bar = canvas.colorbar(image, ax=[tree_ax, strip, ax], location="right", shrink=0.6, pad=0.01, aspect=22, extend="both")
    bar.outline.set_visible(False)
    bar.set_label("Standardised difference,\nsepsis against healthy controls")
    present = [tissue for tissue in TISSUE_COLORS if tissue in set(tissues.values())]
    canvas.legend([Patch(color=TISSUE_COLORS[tissue]) for tissue in present] + [Patch(color=SCORE), Patch(color="#4A6FE3")],
                  [TISSUE_LABELS[tissue][0].upper() + TISSUE_LABELS[tissue][1:] for tissue in present] + ["Up in sepsis", "Down in sepsis"],
                  loc="outside lower center", ncols=4)


def markers(canvas, table):
    """The marker genes of every cell type in every sorted fraction: the shade is the mean of a marker set in a fraction against the other fractions of
    the same study (standardised over its fractions); the number is that mean on the log2 scale. A frame marks the fraction that a marker set should mark."""
    ax = canvas.subplots()
    rows = list(dict.fromkeys(table["row"]))
    sets = [cell for cell in CELL_LABELS if (table["marker_set"] == cell).any()]
    shade = table.pivot(index="row", columns="marker_set", values="standardised").reindex(index=rows, columns=sets)
    means = table.pivot(index="row", columns="marker_set", values="mean").reindex(index=rows, columns=sets)
    image = ax.imshow(shade.to_numpy(dtype=float), aspect="auto", cmap=V.cmap(diverging=True), vmin=-2, vmax=2, interpolation="nearest")
    norm, colors = Normalize(-2, 2), V.cmap(diverging=True)
    own = table.drop_duplicates("row").set_index("row")
    for (row, column), value in np.ndenumerate(means.to_numpy(dtype=float)):
        red, green, blue, _ = colors(norm(shade.iat[row, column]))
        ax.text(column, row, f"{value:.1f}", ha="center", va="center", fontsize=8, color="white" if 0.299 * red + 0.587 * green + 0.114 * blue < 0.5 else V.INK)
        if own.loc[rows[row], "fraction"] == sets[column]:
            ax.add_patch(Rectangle((column - 0.5, row - 0.5), 1, 1, facecolor="none", edgecolor="#111111", linewidth=1.6))
    studies = list(own["study_accession"])
    for row in range(1, len(rows)):
        if studies[row] != studies[row - 1]:
            ax.axhline(row - 0.5, color="white", linewidth=3)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_xticks(range(len(sets)), [CELL_LABELS[cell] for cell in sets], rotation=40, ha="right", rotation_mode="anchor")
    ax.set_xlabel("Marker genes of")
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    V.panel(canvas, ax, "E", "Sorted cells carry their markers")
    bar = canvas.colorbar(image, ax=ax, location="right", shrink=0.7, pad=0.02, aspect=20)
    bar.outline.set_visible(False)
    bar.set_label("Standardised mean log2")


def splits(canvas, samples, units, columns):
    """One small density plot for every tested unit (a table, or one tissue of a table that the validation stage split by tissue): the difference that
    splits its samples, by the sex read from it, every curve in proportion to its samples; below, how many recorded sexes agree."""
    axes = canvas.subplots(int(np.ceil(len(units) / columns)), columns).ravel()
    for ax, row in zip(axes, units.itertuples(index=False)):
        part = samples[samples["unit"] == row.unit]
        values = part["difference"].to_numpy()
        pad = (values.max() - values.min()) * 0.08
        grid = np.linspace(values.min() - pad, values.max() + pad, 300)
        top = 0.0
        for sex in ("female", "male", "ambiguous"):
            own = part.loc[part["expression_sex"] == sex, "difference"].to_numpy()
            density = V.kde(own, grid)
            if density is None:
                if len(own):
                    ax.vlines(own, 0, 0.25, transform=ax.get_xaxis_transform(), color=SEX_COLORS[sex], linewidth=1.0)
                continue
            density = density * len(own)
            top = max(top, density.max())
            ax.fill_between(grid, density, color=SEX_COLORS[sex], alpha=0.3, linewidth=0)
            ax.plot(grid, density, color=SEX_COLORS[sex], linewidth=1.0)
        ax.axvline(row.midpoint, color=V.INK, linewidth=0.6, linestyle=":")
        ax.set_ylim(0, top * 1.08 if top > 0 else 1)
        tissue = TISSUE_LABELS.get(row.tissue, CELL_LABELS.get(row.tissue, row.tissue.replace("_", " ")))
        title = "\n".join(line for line in (row.study_accession, row.platform, tissue) if line)
        ax.set_title(title + ("" if row.female_gene_read else " *"), fontsize=8, pad=2, color=MARK if row.marked else V.INK)
        ax.set_xlabel(f"{V.number(row.concordant)}/{V.number(row.compared)}" if row.compared else "no sex recorded", fontsize=8, labelpad=1, color=V.INK)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
    for ax in axes[len(units):]:
        ax.set_axis_off()
    handles = [Patch(color=SEX_COLORS[sex], alpha=0.6) for sex in ("female", "male", "ambiguous")]
    canvas.legend(handles, ["Read as female", "Read as male", "Ambiguous"], loc="outside lower center", ncols=3)


def pair_table(canvas, table):
    """Every two tables that share samples, with the agreement of the two records of their samples: one shaded column for every measure."""
    ax = canvas.subplots()
    columns = {"median_centred_correlation": "Median,\ncentred", "lowest_centred_correlation": "Lowest,\ncentred", "median_spearman": "Median,\nSpearman",
               "lowest_spearman": "Lowest,\nSpearman", "similar_share": "Most similar\nrecord"}
    shades, norm = V.cmap(), Normalize(0.4, 1.0)
    image = np.ones((len(table), len(columns), 3))
    for column, key in enumerate(columns):
        for row, value in enumerate(table[key].to_numpy(dtype=float)):
            image[row, column] = to_rgb(shades(norm(value)))
            text = f"{table['most_similar'].iat[row]} of {table['pairs'].iat[row]}" if key == "similar_share" else f"{value:.2f}"
            ax.text(column, row, text, ha="center", va="center", fontsize=8, color=V.INK)
    ax.imshow(image, aspect="auto", interpolation="nearest")
    groups = list(table["cohort_group"])
    for row in range(1, len(table)):
        if groups[row] != groups[row - 1]:
            ax.axhline(row - 0.5, color="white", linewidth=3)
    for column in range(1, len(columns)):
        ax.axvline(column - 0.5, color="white", linewidth=1.5)
    ax.set_yticks(range(len(table)), table["label"])
    ax.set_xticks(range(len(columns)), list(columns.values()))
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    bar = canvas.colorbar(ScalarMappable(norm=norm, cmap=shades), ax=ax, location="right", shrink=0.5, pad=0.02, aspect=20)
    bar.outline.set_visible(False)
    bar.set_label("Correlation, or share of the most similar record")


def main():
    run = V.Run(4)
    config = run.config
    sex = run.table("validation/sex_check_samples.parquet")
    by_table = run.table("validation/sex_check_matrices.parquet")
    pairs = run.table("validation/repeat_pairs.parquet")
    by_pair = run.table("validation/repeat_pairs_by_matrices.parquet")
    scores = run.table("validation/sepsis_score_by_study.parquet")
    score_samples = run.table("validation/sepsis_score_samples.parquet")
    genes = run.table("validation/sepsis_genes_by_study.parquet")
    cells = run.table("validation/cell_markers_by_study.parquet")
    policy = run.table("validation/validation_genes.parquet")
    matrices = run.table("annotation/gene_matrix_summary.parquet", ["accession", "platform", "value_scale", "gene_matrix"])
    tissues = run.table("cohort/cohort_samples.parquet", ["study_accession", "sample_accession", "platform", "tissue"])
    display = config["display"]

    by_table["label"] = [name(study, platform) for study, platform in zip(by_table["study_accession"], by_table["platform"])]
    by_table["marked"] = by_table["recorded_sex_unreliable"]
    by_table["compared"] = by_table["concordant"] + by_table["discordant"]
    compared = by_table[by_table["compared"] > 0].sort_values(["concordance", "compared"], ascending=[True, False]).reset_index(drop=True)
    run.add("Fig5", "sex_by_table", compared[["label", "concordant", "discordant", "compared", "concordance", "marked"]])
    # The validation stage splits a matrix of several tissues within each tissue; the midpoint and gap of every split, or of the whole matrix, are the
    # units of its receipt, so every called sample is placed by the unit that called it.
    receipt = V.read_json(V.ROOT / "outputs/technical_validation/02.json")
    V.require("units" in receipt, "The receipt of the sex check lists no units; rerun the validation stage")
    units = pd.DataFrame(receipt["units"])
    split = set(zip(units.loc[units["tissue"] != "", "study_accession"], units.loc[units["tissue"] != "", "platform"]))
    units = units[units["status"] == "tested"].merge(by_table[["study_accession", "platform", "female_gene_read", "marked"]], on=["study_accession", "platform"],
                                                     validate="many_to_one")
    units = units.sort_values(["study_accession", "platform", "tissue"]).reset_index(drop=True)
    units["unit"] = units["study_accession"] + "|" + units["platform"] + "|" + units["tissue"]
    sex = sex.merge(tissues, on=["study_accession", "sample_accession", "platform"], how="left", validate="one_to_one")
    V.require(sex["tissue"].notna().all(), "Every sample of the sex check must be a sample of the cohort table")
    sex["unit"] = sex["study_accession"] + "|" + sex["platform"] + "|" + np.where([key in split for key in zip(sex["study_accession"], sex["platform"])],
                                                                                    sex["tissue"], "")
    placed = sex.merge(units[["unit", "midpoint", "gap"]], on="unit", how="left")
    V.require(placed.loc[placed["expression_sex"].isin(["female", "male", "ambiguous"]), "gap"].notna().all()
              and set(placed.loc[placed["gap"].notna() & placed["difference"].notna(), "expression_sex"]) <= {"female", "male", "ambiguous"},
              "Every sample read as female, male or ambiguous must be placed by the tested unit that read it, and only those")
    called = placed[placed["gap"].notna() & placed["difference"].notna()].copy()
    called["recorded"] = called["recorded_sex"].where(called["recorded_sex"].isin(["female", "male"]), "not_recorded")
    # Each axis of the plane runs from the mean of one group read in the same unit to that of the other, so that units of every scale meet.
    reads_xist = set(units.loc[units["female_gene_read"].astype(bool), "unit"])
    plane = []
    for unit, part in called[called["unit"].isin(reads_xist)].groupby("unit"):
        female, male = part[part["expression_sex"] == "female"], part[part["expression_sex"] == "male"]
        if female.empty or male.empty:
            continue
        xist = female["female_gene_value"].mean() - male["female_gene_value"].mean()
        y_genes = male["male_genes_value"].mean() - female["male_genes_value"].mean()
        V.require(xist != 0 and y_genes != 0, f"The two sex groups of {unit} must differ in both axes")
        plane.append(part.assign(x=(part["female_gene_value"] - male["female_gene_value"].mean()) / xist,
                                 y=(part["male_genes_value"] - female["male_genes_value"].mean()) / y_genes))
    plane = pd.concat(plane, ignore_index=True)
    plane["disagrees"] = plane["agreement"] == "discordant"
    outside = int((~plane["x"].between(*PLANE) | ~plane["y"].between(*PLANE)).sum())
    run.add("Fig5", "sex_plane", plane[["study_accession", "platform", "sample_accession", "recorded", "expression_sex", "agreement"]].assign(
        x=plane["x"].round(6), y=plane["y"].round(6), outside_view=outside))
    unit_agreement = called.groupby("unit")["agreement"].value_counts().unstack(fill_value=0)
    for column in ("concordant", "discordant"):
        units[column] = units["unit"].map(unit_agreement[column] if column in unit_agreement else pd.Series(dtype=int)).fillna(0).astype(int)
    units["compared"] = units["concordant"] + units["discordant"]

    kinds = pairs["linked_by"].str.split("; ")
    V.require(kinds.map(lambda parts: set(parts) <= set(config["link_labels"])).all(), "A kind of evidence that links two records has no label in the configuration")
    kinds = kinds.map(lambda parts: parts[0] if len(parts) == 1 else "several")
    run.add("Fig5", "repeats", pairs[["sample_id", "study_a", "platform_a", "study_b", "platform_b", "linked_by", "spearman", "centred_correlation",
                                      "most_similar"]].round({"spearman": 6, "centred_correlation": 6}))

    scores = scores[scores["score_computed"]].sort_values("auroc", ascending=False).reset_index(drop=True)
    scores["role"] = np.where(scores["score_derivation"] != "", scores["score_derivation"], np.where(scores["samples_shared_with_derivation"] > 0, "shared", ""))
    run.add("Fig5", "score", scores[["study_accession", "platform", "tissue", "sepsis", "healthy_control", "auroc", "auroc_low", "auroc_high",
                                     "auroc_one_per_patient", "role", "repeat_deposits"]].round({"auroc": 6, "auroc_low": 6, "auroc_high": 6,
                                                                                                  "auroc_one_per_patient": 6}))
    keys = ["study_accession", "platform", "tissue"]
    curves = []
    for row in scores.itertuples(index=False):
        part = score_samples[(score_samples["study_accession"] == row.study_accession) & (score_samples["platform"] == row.platform)
                             & (score_samples["tissue"] == row.tissue)]
        cases = part.loc[part["sample_class"] == "sepsis", "score"].to_numpy()
        others = part.loc[part["sample_class"] == "healthy_control", "score"].to_numpy()
        V.require(len(cases) == row.sepsis and len(others) == row.healthy_control, f"The scored samples of {row.study_accession} {row.tissue} must be its groups")
        fpr, tpr = roc_curve(cases, others)
        # The area under the drawn curve is the area of the validation stage, so the curves are the ones its numbers describe.
        V.require(abs(float(np.trapz(tpr, fpr)) - row.auroc) < 1e-9, f"The ROC curve of {row.study_accession} {row.tissue} must have the stage's area")
        curves.append(pd.DataFrame(dict(stratum=f"{row.study_accession}|{row.platform}|{row.tissue}", study_accession=row.study_accession, platform=row.platform,
                                        tissue=row.tissue, fpr=fpr.round(8), tpr=tpr.round(8))))
    curves = pd.concat(curves, ignore_index=True)
    run.add("Fig5", "roc", curves)

    V.require(set(scores["tissue"]) <= set(TISSUE_COLORS), "A tissue of the score has no colour")
    order = policy[policy["test"].isin(["sepsis_score", "single_gene"])].sort_values("role", kind="stable")
    wide = genes.pivot_table(index="published_symbol", columns=keys, values="standardised_difference", aggfunc="first")
    wide = wide.reindex(index=order["published_symbol"], columns=pd.MultiIndex.from_frame(scores[keys]))
    wide.index = list(order["symbol"])
    tissues_of_study = scores.groupby(["study_accession", "platform"])["tissue"].transform("nunique")
    platforms_of_study = scores.groupby("study_accession")["platform"].transform("nunique")
    labels = [study + (f" {platform}" if platform_count > 1 else "") + (f" {TISSUE_LABELS[tissue]}" if tissue_count > 1 else "")
              for study, platform, tissue, platform_count, tissue_count in zip(scores["study_accession"], scores["platform"], scores["tissue"],
                                                                               platforms_of_study, tissues_of_study)]
    wide.columns = list(zip(labels, scores["tissue"]))
    column_tissue = dict(zip(wide.columns, scores["tissue"]))
    run.add("Fig5", "directions", genes[["study_accession", "platform", "tissue", "published_symbol", "role", "standardised_difference", "as_published"]].round(
        {"standardised_difference": 6}))

    # The marker sets of every sorted fraction, read again from the gene-level tables the validation stage read, on the same log2 values.
    markers_by_cell = {cell: list(part["hgnc_id"]) for cell, part in policy[policy["test"] == "cell_marker"].groupby("role")}
    fractions = set(markers_by_cell) | {"pbmc"}
    wanted = sorted({gene for ids in markers_by_cell.values() for gene in ids})
    marker_rows = []
    for study, platform in cells[["study_accession", "platform"]].drop_duplicates().itertuples(index=False):
        matrix = matrices[(matrices["accession"] == study) & (matrices["platform"] == platform)].iloc[0]
        V.require(matrix["value_scale"] in ("log2", "log2 counts per million"), f"The sorted-cell table of {study} must be on a log2 scale")
        # In double precision, as the validation stage reads them, so that the means agree to the last digits.
        values = pd.read_parquet(run.source(matrix["gene_matrix"].removeprefix("data/processed/")), filters=[("hgnc_id", "in", wanted)]).drop(
            columns="symbol").set_index("hgnc_id").astype("float64")
        own = tissues[(tissues["study_accession"] == study) & (tissues["platform"] == platform) & tissues["tissue"].isin(fractions)]
        sizes = own.groupby("tissue").size()
        kept = sizes[sizes >= 3].index
        for fraction in [cell for cell in CELL_LABELS if cell in kept]:
            samples = [sample for sample in own.loc[own["tissue"] == fraction, "sample_accession"] if sample in values.columns]
            for cell, ids in markers_by_cell.items():
                read = [gene for gene in ids if gene in values.index]
                if read:
                    marker_rows.append(dict(study_accession=study, platform=platform, fraction=fraction, samples=len(samples), marker_set=cell,
                                            mean=float(values.loc[read, samples].mean(axis=0, skipna=True).dropna().mean())))
    marker_table = pd.DataFrame(marker_rows)
    marker_table["standardised"] = marker_table.groupby(["study_accession", "marker_set"])["mean"].transform(
        lambda column: (column - column.mean()) / column.std(ddof=0))
    marker_table["row"] = [f"{CELL_LABELS[fraction]}, {study}" for fraction, study in zip(marker_table["fraction"], marker_table["study_accession"])]
    check = cells.merge(marker_table, left_on=["study_accession", "platform", "cell_type", "cell_type"],
                        right_on=["study_accession", "platform", "fraction", "marker_set"], how="left")
    V.require(check["mean"].notna().all() and np.allclose(check["mean"], check["mean_in_cell_type"], rtol=0, atol=1e-9),
              "The marker means of the sorted fractions must be those of the validation stage")
    run.add("Fig5", "markers", marker_table[["study_accession", "platform", "fraction", "samples", "marker_set", "mean", "standardised"]].round(6))

    agree, total = int(compared["concordant"].sum()), int(compared["compared"].sum())
    # Half-inch columns, so that the panels of the two rows can have widths of their own.
    run.figure("Fig5", (14, 12.6), dict(nrows=2, ncols=28, height_ratios=[5.6, 7.0]), [
        ("sex_plane", (0, slice(0, 10)), (5.0, 5.6), lambda canvas: sex_plane(canvas, plane, agree, total)),
        ("repeats", (0, slice(10, 19)), (4.5, 5.6), lambda canvas: repeats(canvas, pairs, kinds, config)),
        ("roc", (0, slice(19, 28)), (4.5, 5.6), lambda canvas: roc_curves(canvas, curves, scores)),
        ("directions", (1, slice(0, 17)), (8.5, 7.0), lambda canvas: directions(canvas, wide, list(order["role"]), column_tissue, display["difference_limit"])),
        ("markers", (1, slice(17, 28)), (5.5, 7.0), lambda canvas: markers(canvas, marker_table))])

    units["label"] = [name(study, platform) + (f"  {tissue}" if tissue else "") for study, platform, tissue in zip(units["study_accession"], units["platform"],
                                                                                                                units["tissue"])]
    run.add("SuppFig3", "splits", units[["study_accession", "platform", "tissue", "label", "samples", "female_gene_read", "lower_group", "higher_group", "gap",
                                         "separation", "midpoint", "marked", "concordant", "compared"]].round(6))
    run.figure("SuppFig3", (14, 12), dict(nrows=1, ncols=1), [("splits", (0, 0), (14, 12), lambda canvas: splits(canvas, placed, units,
                                                                                                                display["sex_split_columns"]))])
    by_pair = by_pair.sort_values(["cohort_group", "study_a", "study_b"]).reset_index(drop=True)
    by_pair["label"] = [f"{name(a, pa)}  and  {name(b, pb)}" for a, pa, b, pb in zip(by_pair["study_a"], by_pair["platform_a"], by_pair["study_b"], by_pair["platform_b"])]
    by_pair["similar_share"] = by_pair["most_similar"] / by_pair["pairs"]
    run.add("SuppFig4", "pairs_of_tables", by_pair[["cohort_group", "label", "pairs", "lowest_spearman", "median_spearman", "lowest_centred_correlation",
                                                    "median_centred_correlation", "most_similar"]].round(6))
    run.figure("SuppFig4", (11, 9), dict(nrows=1, ncols=1), [("pairs_of_tables", (0, 0), (11, 9), lambda canvas: pair_table(canvas, by_pair))])
    run.finish()


if __name__ == "__main__":
    main()
