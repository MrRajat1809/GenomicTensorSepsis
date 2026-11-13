"""Draw how the rows were mapped to genes and why the studies are not merged: mapping outcomes, the genes every table holds, the scale of every table
and a stacked view: Fig4."""
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import visualization_common as V

LOG_SCALES = ("log2", "log2_like", "log2 counts per million")
SCALE_COLORS = {"log2": "#457B9D", "log2_like": "#8ECAE6", "log2 counts per million": "#E76F51", "linear": "#6D597A", "other": "#777777"}
SCALE_LABELS = {"log2": "log2", "log2_like": "log2-like", "log2 counts per million": "log2 CPM", "linear": "Linear", "other": "Other scale"}
# The values of the log-scale tables are read on this range; the tails beyond it are left out of the drawn densities.
LOG_RANGE = (-8.0, 20.0)


def signed_log(values):
    """A logarithm that keeps zero and the sign, for the tables whose values are not on a log scale."""
    return np.sign(values) * np.log10(1 + np.abs(values))


def mapping(canvas, counts, config):
    """Probe-level rows of every technology, flowing to their mapping outcome."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_xlim(-0.47, 1.42)
    ax.set_ylim(-0.01, 1.01)
    flows = counts.stack().rename("value").reset_index()
    flows.columns = ["source", "target", "value"]
    flows["target"] = "map:" + flows["target"]
    families = list(counts.index)
    outcomes = [f"map:{name}" for name in config["mapping_order"] if counts[name].sum() > 0]
    colors = {**{name: config["family_colors"][name] for name in families}, **{f"map:{name}": config["mapping_colors"][name] for name in config["mapping_order"]}}
    labels = {**{name: f"{config['family_labels'][name]}\n{counts.loc[name].sum() / 1e6:.2f} M" for name in families},
              **{f"map:{name}": f"{config['mapping_labels'][name]}\n{counts[name].sum() / counts.to_numpy().sum():.0%}" for name in config["mapping_order"]}}
    V.alluvial(ax, [families, outcomes], [(flows, "target")], colors, labels, xs=[0.0, 1.0], gap=0.07)
    V.panel(canvas, ax, "A", "Rows mapped to genes")


def presence(canvas, held, order, genes, curve, config):
    """Which genes every gene-level table holds: one row for every table, grouped by technology, and one column for every gene, ordered by the number of
    tables that hold it; above, that number for every gene."""
    top, ax = canvas.subplots(2, 1, sharex=True, height_ratios=[1.0, 2.6])
    # Two blank rows between technologies keep the groups, and the names of the small ones, apart.
    rows, previous = [], None
    for index, family in order:
        if previous is not None and family != previous:
            rows += [(None, None)] * 2
        rows.append((index, family))
        previous = family
    image = np.full((len(rows), len(genes), 3), 255, dtype=np.uint8)
    position = {gene: index for index, gene in enumerate(genes)}
    for row, (index, family) in enumerate(rows):
        if index is not None:
            image[row, [position[gene] for gene in held[index]]] = (np.array(to_rgb(config["family_colors"][family])) * 255).astype(np.uint8)
    ax.imshow(image, aspect="auto", interpolation="antialiased", extent=(0, len(genes), len(rows), 0), rasterized=True)
    families = [family for _, family in rows]
    ticks, names = [], []
    for family in dict.fromkeys(name for name in families if name is not None):
        own = [row for row, name in enumerate(families) if name == family]
        ticks.append((own[0] + own[-1] + 1) / 2)
        names.append(f"{config['family_labels'][family]} ({len(own)})")
    ax.set_yticks(ticks, names)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Genes, ordered by the number of tables that hold them")
    ax.xaxis.set_major_formatter(lambda value, _: V.number(value))
    for side in ("left",):
        ax.spines[side].set_visible(False)
    counts = curve["tables_holding"].to_numpy()
    top.fill_between(np.arange(len(counts)) + 0.5, counts, step="mid", color=V.tint("#457B9D", 0.75), linewidth=0)
    top.step(np.arange(len(counts)) + 0.5, counts, where="mid", color="#457B9D", linewidth=1.3)
    tables = len(order)
    for share, offset in ((1.0, (8, 6)), (0.9, (14, 4)), (0.5, (14, 6))):
        need = int(np.ceil(tables * share))
        reach = int((counts >= need).sum())
        top.plot([reach], [need], marker="o", markersize=4.5, color="#457B9D", markeredgecolor="white", markeredgewidth=0.6)
        top.annotate(f"{V.number(reach)} genes in {'all' if share == 1 else 'at least ' + str(need)} tables", (reach, need), xytext=offset,
                     textcoords="offset points", ha="left", va="bottom", fontsize=8, color=V.INK)
    top.set_ylim(0, tables * 1.25)
    top.set_ylabel("Tables")
    V.clean_axis(top, grid="y")
    top.tick_params(axis="x", labelbottom=False)
    ax.set_xlim(0, len(genes))
    V.panel(canvas, top, "B", "The genes of every gene-level table")


def ridges(canvas, densities, others, config):
    """The values of every gene-level table as a density on its own scale, one ridge for every table ordered by its median: the log-scale tables on the
    left, the others on a signed logarithm on the right."""
    axes = canvas.subplots(1, 2, width_ratios=[4.0, 1.25])
    for ax, part, scale in ((axes[0], densities, "log"), (axes[1], others, "other")):
        step = 1.0
        height = 4.0 if scale == "log" else 2.2
        tables = list(dict.fromkeys(part["table"]))
        for rank, table in enumerate(tables):
            own = part[part["table"] == table]
            curve = own["density"].to_numpy()
            curve = curve / curve.max() * height if curve.max() > 0 else curve
            base = (len(tables) - 1 - rank) * step
            color = SCALE_COLORS[own["value_scale"].iloc[0]]
            ax.fill_between(own["value"], base, base + curve, color=color, alpha=0.9, linewidth=0, zorder=rank)
            ax.plot(own["value"], base + curve, color="white", linewidth=0.5, zorder=rank)
        ax.set_yticks([])
        ax.set_ylim(-0.5, (len(tables) - 1) * step + height + 0.5)
        ax.spines["left"].set_visible(False)
    axes[0].set_xlabel("Value as released (log scales)")
    axes[0].set_xlim(*LOG_RANGE)
    axes[1].set_xticks(signed_log(np.array([-100.0, 0.0, 100.0, 1e4])), ["-100", "0", "100", "10⁴"])
    axes[1].set_xlabel("Other scales", loc="right")
    V.panel(canvas, axes[0], "C", "Every table on its own scale")
    present = [name for name in config["scale_order"] if (densities["value_scale"] == name).any() or (others["value_scale"] == name).any()]
    canvas.legend([Patch(color=SCALE_COLORS[name]) for name in present], [SCALE_LABELS[name] for name in present], loc="outside lower center", ncols=3)


def stacked(canvas, points, column, explained, config, letter, title):
    """Illustration only: the samples of all log-scale tables stacked as they are, on the genes all of them hold."""
    ax = canvas.subplots()
    colors = config["family_colors"] if column == "family" else config["class_colors"]
    labels = config["family_labels"] if column == "family" else config["class_labels"]
    shown = points.iloc[np.random.default_rng(config["seed"]).permutation(len(points))]
    V.bubbles(ax, shown["pc1"], shown["pc2"], shown[column].map(colors), size=9, alpha=0.6, rasterized=True)
    ax.set_xlabel(f"Component 1 ({explained[0]:.0%} of variance)")
    ax.set_ylabel(f"Component 2 ({explained[1]:.0%} of variance)")
    names = [name for name in colors if (points[column] == name).any()]
    ax.legend([Line2D([], [], marker="o", linestyle="", markerfacecolor=colors[name], markeredgecolor="white", markersize=7) for name in names],
              [labels[name] for name in names], loc="best", handletextpad=0.2)
    V.clean_axis(ax)
    V.panel(canvas, ax, letter, title)


def main():
    run = V.Run(3)
    config = run.config
    view = config["stacked_view"]
    maps = run.table("annotation/feature_map_summary.parquet")
    summary = run.table("annotation/gene_matrix_summary.parquet").rename(columns={"accession": "study_accession"})
    samples = run.table("cohort/cohort_samples.parquet", ["study_accession", "sample_accession", "sample_class"])
    V.require(set(summary["value_scale"]) <= set(config["scale_order"]) and set(maps["family"]) <= set(config["family_order"]), "A scale or a family has no place")
    V.require((maps[config["mapping_order"]].sum(axis=1) == maps["features"]).all(), "Every probe-level row needs one mapping outcome of the configuration")

    counts = maps.groupby("family")[config["mapping_order"]].sum().reindex([name for name in config["family_order"] if (maps["family"] == name).any()])
    shares = counts.div(counts.sum(axis=1), axis=0)
    run.add("Fig4", "mapping", pd.concat([counts.add_suffix("_rows"), shares.round(6).add_suffix("_share")], axis=1).rename_axis("family").reset_index())
    ranked = summary.sort_values(["genes", "study_accession", "platform"], ascending=[False, True, True]).reset_index(drop=True)
    run.add("Fig4", "genes_per_table", ranked[["study_accession", "platform", "family", "genes", "protein_coding_genes"]])

    # First pass: which genes every table holds. Second pass: the values, once.
    summary["path"] = [run.source(name.removeprefix("data/processed/")) for name in summary["gene_matrix"]]
    held = {row.Index: set(pd.read_parquet(row.path, columns=["hgnc_id"])["hgnc_id"]) for row in summary.itertuples()}
    presence_count = pd.Series([gene for genes in held.values() for gene in genes]).value_counts()
    genes = sorted(presence_count.index, key=lambda gene: (-presence_count[gene], gene))
    curve = pd.DataFrame({"tables": range(1, len(summary) + 1)})
    curve["genes"] = [int((presence_count >= count).sum()) for count in curve["tables"]]
    run.add("Fig4", "shared_genes", curve)
    family_rank = {name: position for position, name in enumerate(config["family_order"])}
    order = [(row.Index, row.family) for row in summary.assign(rank=summary["family"].map(family_rank)).sort_values(
        ["rank", "genes", "study_accession", "platform"], ascending=[True, False, True, True]).itertuples()]
    holding = pd.DataFrame({"tables_holding": presence_count.reindex(genes).to_numpy()})

    chosen = summary[summary["value_scale"].isin(view["scales"]) & (summary["genes"] >= view["minimum_genes"])]
    common = sorted(set.intersection(*(held[index] for index in chosen.index)))
    V.require(len(common) >= view["genes"], "Too few genes are shared by the tables of the stacked view")
    bins = config["display"]["ridge_bins"]
    edges = np.linspace(*LOG_RANGE, bins + 1)
    quantiles, parts, origin, curves = [], [], [], []
    for row in summary.itertuples():
        frame = pd.read_parquet(row.path).drop(columns="symbol").set_index("hgnc_id")
        values = frame.to_numpy(dtype="float32")
        found = np.nanquantile(values, config["quantiles"])
        table = f"{row.study_accession}|{row.platform}"
        quantiles.append(dict(study_accession=row.study_accession, platform=row.platform, value_scale=row.value_scale, samples=values.shape[1],
                              **{f"q{int(round(level * 100)):02d}": float(value) for level, value in zip(config["quantiles"], found)}))
        flat = values[np.isfinite(values)]
        if row.value_scale in LOG_SCALES:
            grid, density, outside = (edges[:-1] + edges[1:]) / 2, V.smooth_density(flat, edges), float(((flat < edges[0]) | (flat > edges[-1])).mean())
        else:
            own = np.linspace(*np.quantile(signed_log(flat), [0.001, 0.999]), bins + 1)
            grid, density, outside = (own[:-1] + own[1:]) / 2, V.smooth_density(signed_log(flat), own), 0.0
        curves.append(pd.DataFrame(dict(table=table, value_scale=row.value_scale, median=float(found[2]), value=grid.round(6), density=density.round(8),
                                        share_outside=round(outside, 6))))
        if row.Index in chosen.index:
            parts.append(frame.loc[common].to_numpy(dtype="float32"))
            origin.extend((row.study_accession, row.platform, row.family, name) for name in frame.columns)
        print(f"[03]   {row.study_accession} {row.platform}: {values.shape[0]:,} genes, {values.shape[1]:,} samples", flush=True)
    scale_rank = {name: position for position, name in enumerate(config["scale_order"])}
    quantiles = pd.DataFrame(quantiles).assign(order=lambda frame: frame["value_scale"].map(scale_rank)).sort_values(["order", "q50"]).drop(columns="order")
    run.add("Fig4", "scales", quantiles)
    curves = pd.concat(curves, ignore_index=True)
    # Ridges from the highest median at the top to the lowest at the bottom.
    curves = curves.sort_values(["median", "table", "value"], ascending=[False, True, True], kind="stable")
    run.add("Fig4", "scale_densities", curves)
    log_curves = curves[curves["value_scale"].isin(LOG_SCALES)]
    other_curves = curves[~curves["value_scale"].isin(LOG_SCALES)]

    # The tables are stacked as they are: no centring within a study, no correction. That is the point of the panel.
    stack = np.concatenate(parts, axis=1).T
    complete = np.isfinite(stack).all(axis=0)
    stack = stack[:, complete]
    keep = np.sort(np.argsort(stack.var(axis=0))[::-1][:view["genes"]])
    stack = stack[:, keep].astype("float64")
    stack -= stack.mean(axis=0)
    values, vectors = np.linalg.eigh(stack.T @ stack / (len(stack) - 1))
    top = vectors[:, ::-1][:, :2]
    top = top * np.sign(top[np.abs(top).argmax(axis=0), range(2)])
    scores = stack @ top
    explained = values[::-1][:2] / values.sum()
    points = pd.DataFrame(origin, columns=["study_accession", "platform", "family", "sample_accession"])
    points = points.merge(samples, on=["study_accession", "sample_accession"], how="left", validate="one_to_one")
    points["table"] = points["study_accession"] + "|" + points["platform"]
    points["pc1"], points["pc2"] = scores[:, 0].round(4), scores[:, 1].round(4)
    V.require(points["sample_class"].notna().all() and len(points) == len(stack), "Every stacked sample must be a sample of the cohort table")
    run.add("Fig4", "stacked_view", points[["table", "sample_accession", "family", "sample_class", "pc1", "pc2"]].assign(
        explained_1=round(float(explained[0]), 6), explained_2=round(float(explained[1]), 6), genes_shared=int(complete.sum()), genes_used=len(keep)))

    run.figure("Fig4", (14, 11.6), dict(nrows=2, ncols=6, height_ratios=[5.4, 6.2]), [
        ("mapping", (0, slice(0, 2)), (4.6, 5.4), lambda canvas: mapping(canvas, counts, config)),
        ("presence", (0, slice(2, 6)), (9.4, 5.4), lambda canvas: presence(canvas, held, order, genes, holding, config)),
        ("scales", (1, slice(0, 2)), (4.6, 6.2), lambda canvas: ridges(canvas, log_curves, other_curves, config)),
        ("stacked_technology", (1, slice(2, 4)), (4.7, 6.2), lambda canvas: stacked(canvas, points, "family", explained, config, "D",
                                                                                     "Stacked without correction, by technology")),
        ("stacked_class", (1, slice(4, 6)), (4.7, 6.2), lambda canvas: stacked(canvas, points, "sample_class", explained, config, "E",
                                                                                "The same view, by class"))])
    run.finish()


if __name__ == "__main__":
    main()
