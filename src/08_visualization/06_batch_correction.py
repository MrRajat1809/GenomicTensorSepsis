"""Draw the removal of study differences: the corrected tables before and after, what their components follow, the held-out test, study mixing and
where every sample record went: Fig6, SuppFig6."""
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

import visualization_common as V

FACTORS = ("family", "study_accession", "sample_class", "sample_class_anchored")
METHOD_MARKERS = {"healthy_anchored": "o", "matched_groups": "^"}
GROUP_MARKERS = {"bulk_blood": "o", "pbmc": "^"}
FLOOR = "#9A9A9A"
FALLS, RISES = "#457B9D", "#DB4325"


def explained_by(scores, labels):
    """The share of the variance of a component that the means of the labels explain (the R squared of a one-way layout)."""
    table = pd.DataFrame({"score": scores, "label": np.asarray(labels)})
    means = table.groupby("label")["score"].transform("mean")
    return float(((means - table["score"].mean()) ** 2).sum() / ((table["score"] - table["score"].mean()) ** 2).sum())


def components(values, count):
    """The first two principal components of the samples (columns) on the given number of most variable genes (rows); the sign is fixed by the largest loading."""
    top = np.sort(np.argsort(values.var(axis=1))[::-1][:count])
    centred = values[top].T.astype("float64")
    centred -= centred.mean(axis=0)
    eigenvalues, vectors = np.linalg.eigh(centred.T @ centred / (len(centred) - 1))
    lead = vectors[:, ::-1][:, :2]
    lead = lead * np.sign(lead[np.abs(lead).argmax(axis=0), range(2)])
    return centred @ lead, eigenvalues[::-1][:2] / eigenvalues.sum()


def views(run, group, count):
    """Every sample of a corrected table on its first two components, before the correction (aligned values) and after it, on the genes both hold everywhere."""
    samples = run.table(f"correction/corrected/{group}_samples.parquet")
    frame = run.table(f"correction/corrected/{group}_expression.parquet")
    after = frame.drop(columns=["hgnc_id", "symbol"]).to_numpy(dtype="float32")
    V.require(list(frame.columns[2:]) == list(samples["sample_id"]), f"The corrected table of {group} must follow its sample table")
    keep = np.isfinite(after).all(axis=1)
    genes, after = list(frame.loc[keep, "hgnc_id"]), after[keep]
    del frame
    before = np.full_like(after, np.nan)
    position = {name: index for index, name in enumerate(samples["sample_id"])}
    for (study, platform), part in samples.groupby(["study_accession", "platform"], sort=False):
        aligned = run.table(f"correction/aligned/{group}/{study}__{platform}.parquet", ["hgnc_id", *part["sample_accession"]]).set_index("hgnc_id").reindex(genes)
        before[:, [position[name] for name in part["sample_id"]]] = aligned[list(part["sample_accession"])].to_numpy(dtype="float32")
    keep = np.isfinite(before).all(axis=1)
    points = samples[["sample_id", "study_accession", "platform", "family", "sample_class", "tissue", "unit_tier", "unit_method"]].copy()
    explained, shares = {}, []
    for stage, values in (("before", before[keep]), ("after", after[keep])):
        scores, explained[stage] = components(values, count)
        points[f"pc1_{stage}"], points[f"pc2_{stage}"] = scores[:, 0].round(4), scores[:, 1].round(4)
        # The class is also read over the studies corrected on their healthy controls alone: their sepsis samples took no part in their correction,
        # while the matched-group method lines the class means of the other studies up by construction.
        anchored = (points["unit_method"] == "healthy_anchored").to_numpy()
        for factor in FACTORS:
            own = anchored if factor == "sample_class_anchored" else np.ones(len(points), dtype=bool)
            labels = points["sample_class" if factor == "sample_class_anchored" else factor].to_numpy()[own]
            shares.append(dict(group=group, stage=stage, factor=factor, samples=int(own.sum()), component_1=round(explained_by(scores[own, 0], labels), 6),
                               component_2=round(explained_by(scores[own, 1], labels), 6), variance_1=round(float(explained[stage][0]), 6),
                               variance_2=round(float(explained[stage][1]), 6), genes_complete=int(keep.sum()), genes_used=min(count, int(keep.sum()))))
    return points, explained, pd.DataFrame(shares)


def scatter(ax, points, stage, column, config, explained, seed, size, shapes=False):
    """The samples as translucent bubbles in a random order, so that no study is drawn on top of all others; shapes tell the correction of a sample's
    study apart where colour alone would not."""
    colors = config["family_colors"] if column == "family" else config["class_colors"]
    labels = config["family_labels"] if column == "family" else config["class_labels"]
    shown = points.iloc[np.random.default_rng(seed).permutation(len(points))]
    for method, marker in (METHOD_MARKERS.items() if shapes else [(None, "o")]):
        part = shown if method is None else shown[shown["unit_method"] == method]
        V.bubbles(ax, part[f"pc1_{stage}"], part[f"pc2_{stage}"], part[column].map(colors), size=size if marker == "o" else size * 1.3, alpha=0.6,
                  marker=marker, rasterized=True)
    ax.set_xlabel(f"Component 1 ({explained[stage][0]:.0%} of variance)")
    ax.set_ylabel(f"Component 2 ({explained[stage][1]:.0%} of variance)")
    names = [name for name in colors if (points[column] == name).any()]
    handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor=colors[name], markeredgecolor="white", markersize=7) for name in names]
    texts = [labels[name] for name in names]
    if shapes:
        handles += [Line2D([], [], marker=marker, linestyle="", markerfacecolor=V.MUTED, markeredgecolor="white", markersize=7) for marker in METHOD_MARKERS.values()]
        texts += ["Study corrected on its healthy controls", "Study corrected by matched groups"]
    V.clean_axis(ax)
    return handles, texts


def view(canvas, points, explained, stage, column, letter, title, config):
    ax = canvas.subplots()
    handles, texts = scatter(ax, points, stage, column, config, explained, config["seed"], 9, shapes=column == "sample_class")
    V.panel(canvas, ax, letter, title)
    # The legend goes below the panel, where it covers no samples.
    canvas.legend(handles, texts, loc="outside lower center", ncols=2, handletextpad=0.2, columnspacing=1.0)


def factors(canvas, shares, config):
    """What the first two components follow: for technology, study and class, the share of the variance of each component that it explains, as an arrow
    from before the correction (open) to after it (filled), blue where the share falls and red where it rises."""
    axes = canvas.subplots(1, 2, sharey=True)
    labels = config["correction"]["factor_labels"]
    table = shares.set_index(["stage", "factor"])
    for ax, component, name in zip(axes, ("component_1", "component_2"), ("Component 1", "Component 2")):
        for row, factor in enumerate(FACTORS):
            before, after = float(table.loc[("before", factor), component]), float(table.loc[("after", factor), component])
            color = FALLS if after < before else RISES
            if abs(after - before) > 0.04:
                ax.annotate("", xy=(after, row), xytext=(before, row), arrowprops=dict(arrowstyle="-|>", color=color, linewidth=1.8, shrinkA=6, shrinkB=7,
                                                                                         mutation_scale=13), zorder=2)
            ax.scatter([before], [row], s=70, facecolor="white", edgecolor=V.MUTED, linewidth=1.3, zorder=3)
            V.bubbles(ax, [after], [row], color, size=90, alpha=0.9, zorder=4)
            # Each value stands on the outer side of its marker, away from the arrow.
            outward = 1 if after < before else -1
            ax.text(before + outward * 0.045, row, f"{before:.2f}", ha="left" if outward > 0 else "right", va="center", fontsize=8, color=V.INK)
            ax.text(after - outward * 0.045, row, f"{after:.2f}", ha="right" if outward > 0 else "left", va="center", fontsize=8, color=V.INK)
        ax.set_xlim(-0.22, 1.18)
        ax.set_ylim(len(FACTORS) - 0.5, -0.6)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1"])
        ax.set_title(name, fontsize=9, pad=4)
        ax.set_xlabel("Share of its variance explained")
        V.clean_axis(ax, grid="x", categorical="y")
    axes[0].set_yticks(range(len(FACTORS)), [labels[name] for name in FACTORS])
    axes[1].tick_params(axis="y", labelleft=False)
    V.panel(canvas, axes[0], "D", "Study and technology no longer drive the components", across=True)
    canvas.legend([Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=V.MUTED, markersize=8),
                   Line2D([], [], marker="o", linestyle="", markerfacecolor=FALLS, markeredgecolor="white", markersize=9),
                   Line2D([], [], marker="o", linestyle="", markerfacecolor=RISES, markeredgecolor="white", markersize=9)],
                  ["Before correction", "After, lower", "After, higher"], loc="outside lower center", ncols=3)


def held_out(canvas, table, config):
    """The hidden healthy controls of every held-out study, after each method: one line for every study across the methods, bold the median, grey the
    dispersion that each method leaves among the studies it learns from, dotted the error that sampling the hidden controls alone gives."""
    axes = canvas.subplots(1, 2)
    policy = config["correction"]
    order = policy["method_order"]
    V.require(set(table["group"]) <= set(GROUP_MARKERS), "A sample type of the held-out test has no marker")
    for ax, column, label in ((axes[0], "location_error", "Distance of the hidden controls' means\n(root mean square, within-study SDs)"),
                              (axes[1], "scale_error", "Error of their spread\n(root mean square of the log ratio)")):
        for (group, _, _), part in table.groupby(["group", "study_accession", "platform"]):
            values = part.set_index("method")[column].reindex(order).to_numpy(dtype=float)
            color = policy["group_colors"][group]
            ax.plot(range(len(order)), values, color=color, alpha=0.28, linewidth=0.9)
            V.bubbles(ax, range(len(order)), values, color, size=16, alpha=0.5, marker=GROUP_MARKERS[group])
        for group, marker in GROUP_MARKERS.items():
            own = table[table["group"] == group]
            if own.empty:
                continue
            middle = own.groupby("method")[column].median().reindex(order).to_numpy(dtype=float)
            ax.plot(range(len(order)), middle, color=policy["group_colors"][group], linewidth=2.6, marker=marker, markersize=7, markeredgecolor="white",
                    markeredgewidth=0.8, zorder=4)
            if column == "location_error":
                spread = own.groupby("method")["pool_dispersion"].median().reindex(order).to_numpy(dtype=float)
                ax.hlines(spread, np.arange(len(order)) - 0.18, np.arange(len(order)) + 0.18, color=FLOOR, linewidth=3.0, alpha=0.8, zorder=1)
        if column == "location_error":
            ax.axhline(table.drop_duplicates(["group", "study_accession", "platform"])["sampling_floor"].median(), color=V.MUTED, linestyle=":", linewidth=1.0)
        ax.set_xticks(range(len(order)), [policy["method_labels"][name] for name in order])
        ax.set_xlim(-0.35, len(order) - 0.65)
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
        V.clean_axis(ax, grid="y", categorical="x")
    V.panel(canvas, axes[0], "E", "Hidden healthy controls, put back by each method", across=True)
    counts = {group: table.loc[table["group"] == group, ["study_accession", "platform"]].drop_duplicates().shape[0] for group in GROUP_MARKERS}
    shown = [group for group in GROUP_MARKERS if counts[group]]
    handles = [Line2D([], [], color=policy["group_colors"][group], linewidth=2.6, marker=GROUP_MARKERS[group], markersize=7, markeredgecolor="white")
               for group in shown]
    canvas.legend(handles + [Line2D([], [], color=FLOOR, linewidth=3.0), Line2D([], [], color=V.MUTED, linestyle=":", linewidth=1.0)],
                  [f"{policy['groups'][group]}: median of {counts[group]} held-out studies (thin: each study)" for group in shown]
                  + ["Dispersion left among the studies learnt from", "Sampling floor of the hidden controls"], loc="outside lower center", ncols=2)


def mixing(canvas, table, config):
    """The study mixing of every corrected unit and class, before against after, in one space."""
    ax = canvas.subplots()
    policy = config["correction"]
    limit = max(1.0, float(table[["before", "after"]].max().max()) + 0.03)
    low = min(0.0, float(table[["before", "after"]].min().min()) - 0.03)
    ax.plot([low, limit], [low, limit], color=V.MUTED, linestyle="--", linewidth=0.8)
    for tier, color in policy["tier_colors"].items():
        part = table[table["tier"] == tier]
        V.bubbles(ax, part["before"], part["after"], color, size=18 + 2.2 * np.sqrt(part["samples"]), alpha=0.7)
    # A failed check leaves a study out only where the check decides for its tier; a failed check that does not decide is marked apart.
    failed = table[~table["passed"] & table["gate"]]
    ax.scatter(failed["before"], failed["after"], s=70, facecolor="none", edgecolor=V.INK, linewidth=0.9)
    noted = table[~table["passed"] & ~table["gate"]]
    ax.scatter(noted["before"], noted["after"], s=70, marker="s", facecolor="none", edgecolor=V.MUTED, linewidth=0.9)
    ax.set_xlim(low, limit)
    ax.set_ylim(low, limit)
    ax.set_xlabel("Own-study neighbours beyond chance, before the correction")
    ax.set_ylabel("Own-study neighbours beyond chance, after")
    V.clean_axis(ax, grid="both")
    V.panel(canvas, ax, "F", "Study mixing of every study and class")
    handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor=color, markeredgecolor="white", markersize=8) for color in policy["tier_colors"].values()]
    handles += [Line2D([], [], marker="o", linestyle="", markerfacecolor="none", markeredgecolor=V.INK, markersize=8),
                Line2D([], [], marker="s", linestyle="", markerfacecolor="none", markeredgecolor=V.MUTED, markersize=8)]
    canvas.legend(handles, [f"{policy['tier_labels'][tier]} ({int((table['tier'] == tier).sum())})" for tier in policy["tier_colors"]]
                  + [f"Left out: failed a deciding check ({len(failed)})", f"Kept: failed a check that does not decide ({len(noted)})"],
                  loc="outside lower center", ncols=2)


def ledger_flow(canvas, flows, config):
    """Every sample record of the resource, from the sample type it was corrected with to its outcome."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_xlim(-0.2, 1.45)
    ax.set_ylim(-0.01, 1.01)
    policy = config["correction"]
    first, second = flows
    groups = [name for name in policy["group_order"] if (first["target"] == name).any()]
    outcomes = list(dict.fromkeys(second.sort_values("rank")["target"]))
    totals = {**first.groupby("target")["value"].sum().to_dict(), **second.groupby("target")["value"].sum().to_dict()}
    colors = {**{name: policy["group_colors"][name] for name in groups}, **{name: policy["outcome_colors"].get(name, "#C8C8C8") for name in outcomes}}
    labels = {**{name: f"{policy['group_labels'][name]}\n{V.number(totals[name])}" for name in groups},
              **{name: f"{policy['outcome_labels'][name]}  {V.number(totals[name])}" for name in outcomes}}
    V.alluvial(ax, [groups, outcomes], [(second, "source")], colors, labels, xs=[0.0, 1.0], gap=0.035)
    V.panel(canvas, ax, "G", f"Where all {V.number(first['value'].sum())} sample records went")


def slopes(canvas, table, check, letter, title, label, config):
    """One check of every corrected unit: its value before and after the correction, joined, with the limit that a unit must reach after it."""
    ax = canvas.subplots()
    policy = config["correction"]
    part = table[table["check"] == check]
    limit = float(part["limit"].iloc[0])
    jitter = np.random.default_rng(config["seed"])
    for tier, color in policy["tier_colors"].items():
        own = part[part["tier"] == tier]
        offset = jitter.uniform(-0.06, 0.06, len(own))
        ax.plot(np.vstack([offset, 1 + offset]), own[["before", "after"]].to_numpy().T, color=color, alpha=0.35, linewidth=0.9)
        V.bubbles(ax, offset, own["before"], color, size=20, alpha=0.65)
        V.bubbles(ax, 1 + offset, own["after"], color, size=20, alpha=0.65)
    ax.axhline(limit, color="#DB4325", linestyle=":", linewidth=1.0)
    ax.text(1.48, limit, "limit", ha="right", va="bottom", fontsize=8, color="#DB4325")
    ax.set_xticks([0, 1], ["Before", "After"])
    ax.set_xlim(-0.3, 1.5)
    ax.set_ylim(min(float(part[["before", "after"]].min().min()), limit) - 0.02, 1.01)
    ax.set_ylabel(label)
    V.clean_axis(ax, grid="y", categorical="x")
    V.panel(canvas, ax, letter, title)
    handles = [Line2D([], [], marker="o", linestyle="-", color=color, markerfacecolor=color, markeredgecolor="white", markersize=7)
               for color in policy["tier_colors"].values()]
    canvas.legend(handles, [f"{policy['tier_labels'][tier]} ({int((part['tier'] == tier).sum())})" for tier in policy["tier_colors"]], loc="outside lower center",
                  ncols=1)


def pbmc_views(canvas, points, explained, config):
    axes = canvas.subplots(1, 3)
    titles = ("PBMC after alignment, before correction", "After correction, by technology", "After correction, by class")
    for ax, stage, column, title in zip(axes, ("before", "after", "after"), ("family", "family", "sample_class"), titles):
        handles, texts = scatter(ax, points, stage, column, config, explained, config["seed"], 14)
        ax.legend(handles, texts, loc="best", handletextpad=0.2)
        ax.set_title(title, pad=14)
    canvas.text(0.004, 0.996, "A", ha="left", va="top", fontsize=V.LETTER_SIZE, fontweight="bold", color="#111111", gid=V.LETTER)


def deposits(canvas, pairs):
    ax = canvas.subplots()
    top = float(pairs[["distance_before", "distance_after"]].max().max()) * 1.05
    ax.plot([0, top], [0, top], color=V.MUTED, linestyle="--", linewidth=0.8)
    nearest = pairs["rank_after"] == 1
    V.bubbles(ax, pairs.loc[nearest, "distance_before"], pairs.loc[nearest, "distance_after"], "#457B9D", size=12, alpha=0.5, rasterized=True)
    V.bubbles(ax, pairs.loc[~nearest, "distance_before"], pairs.loc[~nearest, "distance_after"], "#DB4325", size=36, alpha=0.95)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_xlabel("Distance between the two deposits before\n(root mean square over shared genes)")
    ax.set_ylabel("Distance after the correction")
    V.clean_axis(ax, grid="both")
    V.panel(canvas, ax, "D", "Two deposits of one sample")
    ax.legend([Line2D([], [], marker="o", linestyle="", markerfacecolor="#457B9D", markeredgecolor="white", markersize=7),
               Line2D([], [], marker="o", linestyle="", markerfacecolor="#DB4325", markeredgecolor="white", markersize=8)],
              [f"Still nearest its twin ({V.number(nearest.sum())})", f"Not nearest its twin ({V.number((~nearest).sum())})"], loc="upper left")


def main():
    run = V.Run(6)
    config = run.config
    policy = config["correction"]
    count = policy["view_genes"]
    bulk, bulk_explained, bulk_shares = views(run, "bulk_blood", count)
    pbmc, pbmc_explained, pbmc_shares = views(run, "pbmc", count)
    for group, points in (("bulk_blood", bulk), ("pbmc", pbmc)):
        V.require(set(points["unit_method"]) <= set(METHOD_MARKERS) and set(points["family"]) <= set(config["family_colors"])
                  and set(points["sample_class"]) <= set(config["class_colors"]), f"A method, technology or class of the {group} table has no marker or colour")
        run.add("Fig6" if group == "bulk_blood" else "SuppFig6", f"views_{group}", points.assign(group=group))
    shares = pd.concat([bulk_shares, pbmc_shares], ignore_index=True)
    run.add("Fig6", "factors", shares)
    held = run.table("correction/held_out_test.parquet")
    V.require(set(held["method"]) <= set(policy["method_order"]), "A held-out method has no place")
    run.add("Fig6", "held_out", held)
    checks = run.table("correction/unit_checks.parquet")
    applied = checks[checks["applicable"]]
    V.require(set(applied["tier"]) <= set(policy["tier_colors"]), "A tier of the unit checks has no colour")
    mixed = applied[applied["check"] == "study_mixing"]
    run.add("Fig6", "mixing", mixed)
    ledger = run.table("correction/correction_ledger.parquet")
    kinds = np.where(ledger["outcome"] == "in_corrected_table", "in_table_" + ledger["group"], np.where(ledger["outcome"] == "corrected_for_checks_only",
                                                                                                       "repeat_deposit", ledger["reason"]))
    counts = pd.Series(kinds).value_counts().reindex(policy["outcome_order"]).dropna().astype(int)
    V.require(int(counts.sum()) == len(ledger) and set(kinds) <= set(policy["outcome_order"]), "Every record of the ledger needs one place")
    run.add("Fig6", "outcomes", counts.rename_axis("outcome").rename("records").reset_index())
    V.require(set(ledger["group"]) <= set(policy["group_order"]), "A sample type of the ledger has no place")
    # The two corrected tables meet in one outcome, so that the ribbons of their sample types show how each table was filled.
    target = np.where(ledger["outcome"] == "in_corrected_table", "in_corrected_table", np.where(ledger["outcome"] == "corrected_for_checks_only",
                                                                                                "repeat_deposit", ledger["reason"]))
    rank = {name: position for position, name in enumerate(["in_corrected_table", *policy["outcome_order"]])}
    first = ledger.groupby("group").size().rename("value").reset_index().rename(columns={"group": "target"}).assign(source="records")
    second = pd.DataFrame({"source": ledger["group"].to_numpy(), "target": target}).value_counts().rename("value").reset_index()
    second["rank"] = second["target"].map(rank)
    run.add("Fig6", "ledger_flows", pd.concat([first.assign(step="records_to_group"), second.drop(columns="rank").assign(step="group_to_outcome")],
                                              ignore_index=True))
    run.add("SuppFig6", "checks", applied[applied["check"].isin(["sex_genes", "sepsis_score"])])
    pairs = run.table("correction/repeat_test.parquet")
    run.add("SuppFig6", "deposits", pairs)

    run.figure("Fig6", (14, 15), dict(nrows=3, ncols=6, height_ratios=[1.0, 0.9, 1.15]), [
        ("before_technology", (0, slice(0, 2)), (4.7, 4.8), lambda canvas: view(canvas, bulk, bulk_explained, "before", "family", "A",
                                                                                "Bulk blood, before correction", config)),
        ("after_technology", (0, slice(2, 4)), (4.7, 4.8), lambda canvas: view(canvas, bulk, bulk_explained, "after", "family", "B",
                                                                               "After correction, by technology", config)),
        ("after_class", (0, slice(4, 6)), (4.7, 4.8), lambda canvas: view(canvas, bulk, bulk_explained, "after", "sample_class", "C",
                                                                          "After correction, by class", config)),
        ("factors", (1, slice(0, 3)), (7, 4.3), lambda canvas: factors(canvas, bulk_shares, config)),
        ("held_out", (1, slice(3, 6)), (7, 4.3), lambda canvas: held_out(canvas, held, config)),
        ("mixing", (2, slice(0, 3)), (7, 5.5), lambda canvas: mixing(canvas, mixed, config)),
        ("outcomes", (2, slice(3, 6)), (7, 5.5), lambda canvas: ledger_flow(canvas, [first, second], config))])
    run.figure("SuppFig6", (14, 10), dict(nrows=2, ncols=6, height_ratios=[1.0, 1.1]), [
        ("pbmc_views", (0, slice(None)), (14, 4.7), lambda canvas: pbmc_views(canvas, pbmc, pbmc_explained, config)),
        ("sex_genes", (1, slice(0, 2)), (4.7, 5.2), lambda canvas: slopes(canvas, applied, "sex_genes", "B", "Sex genes, one split for all studies",
                                                                          "Agreement with the sex called inside the study", config)),
        ("sepsis_score", (1, slice(2, 4)), (4.7, 5.2), lambda canvas: slopes(canvas, applied, "sepsis_score", "C", "Score against other studies' controls",
                                                                             "AUROC against the healthy controls of the other studies", config)),
        ("deposits", (1, slice(4, 6)), (4.7, 5.2), lambda canvas: deposits(canvas, pairs))])
    run.finish()


if __name__ == "__main__":
    main()
