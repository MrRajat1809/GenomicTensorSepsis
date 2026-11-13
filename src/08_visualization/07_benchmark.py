"""Draw the benchmark of evaluation designs: twins of the test records in training, baseline classifiers, the sources of inflation, the study as a
shortcut and as a label that expression identifies, and sepsis signatures within single studies: Fig7, SuppFig7."""
import hashlib
import math

import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_hex, to_rgb
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Wedge
from matplotlib.textpath import TextPath
from matplotlib.ticker import MaxNLocator

import visualization_common as V

ROLE_MARKERS = {"": "o", "discovery": "s", "validation": "D"}
# The narrowest kernel of a density of AUROCs: five or so values close to one would otherwise make a needle.
NARROWEST = 0.012
# The rings of a donut of panel A, outside in, as outer radius and width in radii of the donut; the patient-safe ring holds no twin and is narrower.
RINGS = {"naive_record_split": (1.0, 0.31), "naive_leave_one_study_out": (0.67, 0.31), "patient_safe": (0.34, 0.15)}
# Inches per radius of a donut, slightly under the drawn scale, so that a share is placed on its arc only where it surely fits.
DONUT_INCHES = 1.08
# The width in inches of each plot of the sources of inflation, slightly under the drawn width, which sets the room left for the differences.
LABEL_INCHES = 1.6
SKIPPED = "#F3F3F3"
# A cloud below this share of its peak is drawn thinner than its outline, so the axis of a cloud reaches at least that far.
FADED = 0.03


def natural(accession):
    """A study accession as its prefix and number, so that GSE9692 comes before GSE13015."""
    prefix = accession.rstrip("0123456789")
    return prefix, int(accession[len(prefix):] or 0)


def darken(color, amount):
    red, green, blue = to_rgb(color)
    return to_hex((red * (1 - amount), green * (1 - amount), blue * (1 - amount)))


def ink_on(color):
    red, green, blue = to_rgb(color)
    return "white" if 0.299 * red + 0.587 * green + 0.114 * blue < 0.5 else V.INK


def bumps(values, grid, narrowest=NARROWEST):
    """A Gaussian kernel density of AUROCs on the grid, with the bandwidth of Scott's rule but no narrower than narrowest, reflected at one so that no
    weight falls beyond it."""
    values = np.asarray(values, dtype=float)
    width = max(narrowest, float(values.std(ddof=1)) * len(values) ** -0.2) if len(values) > 1 else narrowest
    points = np.r_[values, 2 - values]
    return np.exp(-0.5 * ((grid[:, None] - points[None, :]) / width) ** 2).sum(axis=1) / (len(values) * width * math.sqrt(2 * math.pi))


def fading_edge(groups, start, narrowest=NARROWEST):
    """The left end of an AUROC axis, from start down in steps of 0.05, at which the density of every group of values has fallen below FADED of its
    peak, so that every cloud fades out inside the plot rather than being cut by its edge."""
    full = np.linspace(0, 1, 1001)
    peaks = [float(bumps(values, full, narrowest).max()) for values in groups]
    edge = start
    while edge > 0 and any(bumps(values, np.array([edge]), narrowest)[0] >= FADED * peak for values, peak in zip(groups, peaks)):
        edge = round(edge - 0.05, 2)
    return max(edge, 0.0)


def roc_curve(cases, others):
    """The points of the ROC curve of a score that is higher in the cases, one point for every distinct value, so that ties move both rates at once."""
    values = np.concatenate([cases, others])
    labels = np.concatenate([np.ones(len(cases)), np.zeros(len(others))])
    order = np.argsort(-values, kind="mergesort")
    values, labels = values[order], labels[order]
    ends = np.r_[np.flatnonzero(np.diff(values)), len(values) - 1]
    hits = np.cumsum(labels)[ends]
    return np.r_[0.0, (ends + 1 - hits) / len(others)], np.r_[0.0, hits / len(cases)]


def area(fpr, tpr):
    fpr, tpr = np.asarray(fpr, dtype=float), np.asarray(tpr, dtype=float)
    return float(((fpr[1:] - fpr[:-1]) * (tpr[1:] + tpr[:-1]) / 2).sum())


def curve_rows(curves, keys, grid):
    """Every drawn ROC curve read at the false-positive rates of the grid: the rows of a median curve and its middle half."""
    return np.array([np.interp(grid, part["fpr"], part["tpr"]) for _, part in curves.groupby(keys, sort=False)])


def twins_in_training(frame, test, train):
    """The kind of twin that every test record has among the training records: a deposit of the same array, a sample of another array of the same
    identified patient, both, or none."""
    arrays, patients = frame["sample_id"].to_numpy(), frame["patient_key"].to_numpy()
    same_array = pd.Series(arrays).map(pd.Series(arrays[train]).value_counts()).fillna(0).to_numpy() > 0
    known = train & (patients != "")
    by_patient = pd.Series(patients).map(pd.Series(patients[known]).value_counts()).fillna(0).to_numpy()
    by_array = pd.Series(list(zip(patients, arrays))).map(pd.Series(list(zip(patients[known], arrays[known]))).value_counts()).fillna(0).to_numpy()
    same_patient = (patients != "") & (by_patient - by_array > 0)
    kind = np.select([same_array & same_patient, same_array, same_patient], ["both", "array_only", "patient_only"], "none")
    return kind[test], (same_array & test).sum(), (same_patient & test).sum()


def test_records(populations, splits, leakage, policy, settings):
    """Every test record of every repeat of the naive random split, of every fold of the naive leave-one-study-out split and of every fold of the
    patient-safe schemes of the cohort stage, with the kind of its twin in training; every fold must give the counts of the leakage table."""
    stored = leakage.set_index(["population", "task", "scheme", "fold"])
    rows = []
    for task in policy["task_order"]:
        naive = populations[(populations["population"] == "naive") & (populations["task"] == task)].reset_index(drop=True)
        folds = []
        for repeat in range(1, settings["repeats"] + 1):
            digests = np.array([hashlib.sha256(f"{settings['seed']}:{repeat}:{study}:{sample}".encode("utf-8")).hexdigest()
                                for study, sample in zip(naive["study_accession"], naive["sample_accession"])], dtype=object)
            test = np.zeros(len(naive), dtype=bool)
            test[np.argsort(digests, kind="stable")[:math.floor(settings["test_fraction"] * len(naive) + 0.5)]] = True
            folds.append(("naive", "naive_record_split", "naive_record_split", f"repeat_{repeat:02d}", naive, test, ~test))
        for study in sorted(set(naive["study_accession"])):
            test = (naive["study_accession"] == study).to_numpy()
            folds.append(("naive", "naive_leave_one_study_out", "naive_leave_one_study_out", study, naive, test, ~test))
        safe = populations[(populations["population"] == "gts") & (populations["task"] == task)].reset_index(drop=True)
        for (scheme, fold), part in splits[splits["task"] == task].groupby(["scheme", "fold"]):
            groups = safe["patient_safe_group"]
            folds.append(("gts", "patient_safe", scheme, fold, safe, groups.isin(part.loc[part["role"] != "train", "patient_safe_group"]).to_numpy(),
                          groups.isin(part.loc[part["role"] == "train", "patient_safe_group"]).to_numpy()))
        V.require(sum(fold[1] == "naive_record_split" for fold in folds)
                  == int(((leakage["population"] == "naive") & (leakage["task"] == task) & (leakage["scheme"] == "naive_record_split")).sum()),
                  f"Every repeat of the naive random split of the leakage table must be drawn ({task})")
        for population, split, scheme, fold, frame, test, train in folds:
            kind, same_array, same_patient = twins_in_training(frame, test, train)
            row = stored.loc[(population, task, scheme, fold)]
            V.require((int(test.sum()), int((kind != "none").sum()), int(same_array), int(same_patient))
                      == (row["test_records"], row["exposed"], row["exposed_same_array"], row["exposed_same_patient"]),
                      f"The twins rebuilt for fold {fold} of {scheme} ({task}) must be those of the leakage table")
            rows.append(frame.loc[test, ["study_accession", "sample_accession"]].assign(task=task, split=split, scheme=scheme, fold=fold, twin=kind))
    return pd.concat(rows, ignore_index=True)[["task", "split", "scheme", "fold", "study_accession", "sample_accession", "twin"]]


def text_points(text, config, size=8):
    """The width of a text in points."""
    return TextPath((0, 0), text, size=size, prop=FontProperties(family=config["font_family"])).get_extents().width


def text_size(text, config, size=8):
    """The width and height of a text, in radii of a donut."""
    return text_points(text, config, size) / 72 / DONUT_INCHES, size / 72 / DONUT_INCHES


def arc_order(policy):
    """The kinds of twin clockwise from twelve o'clock: the same patient, both, the same array, and last the records without a twin."""
    return [kind for kind in reversed(policy["twin_order"]) if kind != "none"] + ["none"]


def polar(centre, radius, angle):
    return centre[0] + radius * math.cos(math.radians(angle)), centre[1] + radius * math.sin(math.radians(angle))


def upright(angle):
    """A rotation in degrees turned so that the text never reads upside down."""
    return (angle + 90) % 180 - 90


def donut(ax, centre, shares, title, color, policy, config):
    """The rings of one task, outside in, every ring starting at twelve o'clock with the records that have a twin in training, clockwise; each coloured
    arc carries its share, along the arc where it fits, across the ring where only that fits, and outside the donut with a leader line otherwise, and
    the share of all records with a twin stands left of twelve o'clock in its ring."""
    colors, outside = policy["twin_colors"], []
    kinds = arc_order(policy)
    for split, (outer, width) in RINGS.items():
        own = shares[shares["split"] == split].set_index("twin")["share"].reindex(kinds, fill_value=0.0)
        start, middle = 0.0, outer - width / 2
        for kind, share in own.items():
            if share <= 0:
                continue
            ax.add_patch(Wedge(centre, outer, 90 - 360 * (start + share), 90 - 360 * start, width=width, facecolor=colors[kind], edgecolor="white",
                               linewidth=1.0))
            if kind != "none":
                text, angle = f"{share:.1%}", 90 - 360 * (start + share / 2)
                long, high = text_size(text, config)
                arc = 2 * math.pi * middle * share
                if arc >= long + 0.06 and width >= high + 0.06:
                    ax.text(*polar(centre, middle, angle), text, ha="center", va="center", rotation=upright(angle - 90), rotation_mode="anchor",
                            fontsize=8, color=ink_on(colors[kind]))
                elif arc >= high + 0.03 and width >= long + 0.04:
                    ax.text(*polar(centre, middle, angle), text, ha="center", va="center", rotation=upright(angle), rotation_mode="anchor", fontsize=8,
                            color=ink_on(colors[kind]))
                else:
                    outside.append((angle, middle, text, long, high))
            start += share
        total = 1.0 - float(own["none"])
        ax.text(centre[0] - 0.04, centre[1] + middle, f"{total:.1%}" if total else "0%", ha="right", va="center", fontsize=8, color=V.INK)
    ax.text(centre[0], centre[1] + 1.07, title, ha="center", va="bottom", fontsize=9, color=color)
    long, high = text_size(title, config, 9)
    # A share that fits on no arc stands outside the donut at its own angle, or further clockwise where the title or an earlier share stands.
    placed = [(centre[0] - long / 2, centre[0] + long / 2, centre[1] + 1.07, centre[1] + 1.07 + high)]
    for angle, middle, text, long, high in sorted(outside, reverse=True):
        for where in angle - 4 * np.arange(90):
            x, y = polar(centre, RINGS["naive_record_split"][0] + 0.07, where)
            box = (x, x + long, y - high / 2, y + high / 2) if math.cos(math.radians(where)) >= 0 else (x - long, x, y - high / 2, y + high / 2)
            if all(box[1] < other[0] - 0.03 or box[0] > other[1] + 0.03 or box[3] < other[2] - 0.03 or box[2] > other[3] + 0.03 for other in placed):
                break
        else:
            V.require(False, f"No room outside the donut of {title} for the share {text}")
        placed.append(box)
        ax.plot(*zip(polar(centre, middle, angle), polar(centre, RINGS["naive_record_split"][0] + 0.05, where)), color=V.INK, linewidth=0.6,
                solid_capstyle="butt", zorder=5)
        ax.text(x, y, text, ha="left" if math.cos(math.radians(where)) >= 0 else "right", va="center", fontsize=8, color=V.INK)


def twins_panel(canvas, shares, policy, config):
    """Three donuts, one for every task, in one plot: the rings are the naive random split, the naive leave-one-study-out split and the patient-safe
    folds; the fourth corner holds the key of the rings and the kinds of twin."""
    ax = canvas.subplots()
    ax.set_axis_off()
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0)
    ax.update_datalim([(-2.3, -2.38), (2.3, 2.42)])
    centres = [(-1.2, 1.2), (1.2, 1.2), (-1.2, -1.2)]
    for task, centre in zip(policy["task_order"], centres):
        donut(ax, centre, shares[shares["task"] == task], policy["task_labels"][task], policy["task_colors"][task], policy, config)
    # The key: a quarter of the rings at twelve o'clock, each named on its left.
    corner, scale = (1.72, -1.18), 0.66
    for split, (outer, width) in RINGS.items():
        ax.add_patch(Wedge(corner, outer * scale, 0, 90, width=width * scale, facecolor=policy["twin_colors"]["none"], edgecolor="white", linewidth=1.0))
        ax.text(corner[0] - 0.05, corner[1] + (outer - width / 2) * scale, policy["split_labels"][split], ha="right", va="center", fontsize=8, color=V.INK)
    kinds = arc_order(policy)
    ax.legend([Patch(facecolor=policy["twin_colors"][kind], edgecolor="white") for kind in kinds], [policy["twin_labels"][kind] for kind in kinds],
              loc="upper left", bbox_to_anchor=(0.3, -1.32), bbox_transform=ax.transData, ncols=1, handlelength=1.2, handletextpad=0.5, borderaxespad=0)
    V.panel(canvas, ax, "A", "Naive splits put twins on both sides")


def model_shades(color, policy):
    first, second = policy["model_order"]
    return {first: darken(color, 0.3), second: to_hex(V.tint(color, 0.25))}


def half_cloud(ax, values, row, side, color, grid, jitter):
    """One half of a mirrored raincloud, side -1 above the row and 1 below it: the density outermost, then the median and middle half, then every
    value nearest the row."""
    density = bumps(values, grid)
    density = density / density.max() * 0.29
    base = row + side * 0.17
    ax.fill_between(grid, base, base + side * density, color=color, alpha=0.32, linewidth=0)
    ax.plot(grid, base + side * density, color=color, linewidth=1.0)
    first, middle, third = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot([first, third], [row + side * 0.125] * 2, color=V.INK, linewidth=2.2, solid_capstyle="butt", zorder=4)
    ax.plot([middle], [row + side * 0.125], marker="o", markersize=4.2, markerfacecolor="white", markeredgecolor=V.INK, markeredgewidth=1.0, zorder=5)
    V.bubbles(ax, values, row + side * jitter.uniform(0.02, 0.085, len(values)), color, size=11, alpha=0.7, zorder=3)


def baselines(canvas, units, policy, seed):
    """A mirrored raincloud for every scheme of every task: logistic regression above the row, gradient boosting below, each with the density of the
    AUROCs of its units, their median and middle half, and every unit; the number of units at the right."""
    tasks, schemes, models = policy["task_order"], policy["scheme_order"], policy["model_order"]
    axes = canvas.subplots(1, len(tasks), sharey=True)
    jitter = np.random.default_rng(seed)
    for ax, task in zip(axes, tasks):
        color = policy["task_colors"][task]
        shades = model_shades(color, policy)
        own = units[units["task"] == task]
        # Every task on its own scale: the AUROCs against healthy controls lie near one, those of mortality near the middle.
        clouds = [part["auroc"].to_numpy(dtype=float) for _, part in own.groupby(["scheme", "model"])]
        low = fading_edge(clouds, float(np.floor((own["auroc"].min() - 0.02) * 20) / 20))
        grid = np.linspace(low, 1.0, 500)
        for row, scheme in enumerate(schemes):
            counts = set()
            for model, side in zip(models, (-1, 1)):
                values = own.loc[(own["scheme"] == scheme) & (own["model"] == model), "auroc"].to_numpy(dtype=float)
                half_cloud(ax, values, row, side, shades[model], grid, jitter)
                counts.add(len(values))
            ax.text(1.04, row, "n=" + "/".join(str(count) for count in sorted(counts)), transform=ax.get_yaxis_transform(), ha="left", va="center",
                    fontsize=8, color=V.MUTED)
        if low < 0.5:
            ax.axvline(0.5, color=V.MUTED, linestyle="--", linewidth=0.8, zorder=0)
        ax.set_title(policy["task_labels"][task], fontsize=9, color=color, pad=4)
        ax.set_xlim(low, 1.012)
        ticks = np.arange(np.ceil(low * 10) / 10, 1.0001, 0.1)
        ax.set_xticks(ticks, [f"{value:.1f}".rstrip("0").rstrip(".") for value in ticks])
        V.clean_axis(ax, grid="x", categorical="y")
    axes[0].set_yticks(range(len(schemes)), [policy["scheme_labels"][scheme] for scheme in schemes])
    axes[0].set_ylim(len(schemes) - 0.5, -0.5)
    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)
    axes[len(axes) // 2].set_xlabel("AUROC of the scored unit")
    V.panel(canvas, axes[0], "B", "Patient-safe splits lower the baseline AUROC", across=True)
    grey = model_shades("#7A7A7A", policy)
    handles = [(Patch(color=grey[model], alpha=0.45), Line2D([], [], marker="o", linestyle="", markerfacecolor=grey[model], markeredgecolor="white",
                                                                markersize=5)) for model in models]
    handles.append(Line2D([], [], color=V.INK, linewidth=2.2, marker="o", markersize=4.2, markerfacecolor="white", markeredgecolor=V.INK))
    canvas.legend(handles, [f"{policy['model_labels'][models[0]]}, above", f"{policy['model_labels'][models[1]]}, below", "Median and middle half"],
                  loc="outside lower center", ncols=3, columnspacing=1.6, handletextpad=0.5, handlelength=1.6)


def signed(value, sign="+"):
    """A difference in AUROC to three decimals, with its sign unless it rounds to zero."""
    text = f"{value:{sign}.3f}"
    return "0.000" if text in ("+0.000", "-0.000") else text


def label_room(anchors, low, high, inches):
    """The range of a plot the given number of inches wide that holds the values from low to high and every label beside them: an anchor in data
    units, the width of its label in points and the side it stands on."""
    lower, upper = low, high
    for _ in range(60):
        scale = (upper - lower) / inches / 72
        lower = min([low] + [anchor - width * scale for anchor, width, side in anchors if side < 0])
        upper = max([high] + [anchor + width * scale for anchor, width, side in anchors if side > 0])
    return lower, upper


def inflation(canvas, effects, policy, config):
    """Every source of inflation for every task and model as an arrow from the first AUROC of its comparison (open) to the second (filled), the 95%
    interval of the difference as a band around the second and the difference beside it; the number of studies at the right."""
    tasks, order, models = policy["task_order"], policy["effect_order"], policy["model_order"]
    axes = canvas.subplots(1, len(tasks), sharey=True)
    for ax, task in zip(axes, tasks):
        color = policy["task_colors"][task]
        shades = model_shades(color, policy)
        own = effects[effects["task"] == task].set_index(["effect", "model"])
        reach = np.r_[own["first"], own["second"], (own["first"] - own["high"]).dropna(), (own["first"] - own["low"]).dropna()].astype(float)
        low, high = float(reach.min()), float(reach.max())
        span = high - low
        anchors = []
        for row, effect in enumerate(order):
            studies = set()
            for model, offset in zip(models, (-0.17, 0.17)):
                entry = own.loc[(effect, model)]
                first, second, y = float(entry["first"]), float(entry["second"]), row + offset
                shade = shades[model]
                left, right = min(first, second), max(first, second)
                if pd.notna(entry["low"]):
                    band = (first - float(entry["high"]), first - float(entry["low"]))
                    ax.add_patch(Rectangle((band[0], y - 0.085), band[1] - band[0], 0.17, facecolor=shade, alpha=0.2, linewidth=0, zorder=1))
                    left, right = min(left, band[0]), max(right, band[1])
                if abs(second - first) > 0.025 * span:
                    ax.annotate("", xy=(second, y), xytext=(first, y), arrowprops=dict(arrowstyle="-|>", color=shade, linewidth=1.6, shrinkA=4, shrinkB=4.5,
                                                                                         mutation_scale=11), zorder=2)
                ax.scatter([first], [y], s=34, facecolor="white", edgecolor=shade, linewidth=1.2, zorder=3)
                ax.scatter([second], [y], s=40, facecolor=shade, edgecolor="white", linewidth=0.6, zorder=4)
                # The difference stands beyond the arrowhead and its band; a difference that rounds to zero stands right of its circles.
                towards = -1 if second < first and signed(entry["estimate"]) != "0.000" else 1
                label, anchor = signed(entry["estimate"]), left if towards < 0 else right
                ax.annotate(label, xy=(anchor, y), xytext=(towards * 6.5, 0), textcoords="offset points", ha="right" if towards < 0 else "left",
                            va="center", fontsize=8, color=darken(shade, 0.15))
                anchors.append((anchor, text_points(label, config) + 9, towards))
                studies.add(int(entry["studies"]))
            ax.text(1.03, row, "n=" + "/".join(str(count) for count in sorted(studies)), transform=ax.get_yaxis_transform(), ha="left", va="center",
                    fontsize=8, color=V.MUTED)
        # The plot is about LABEL_INCHES wide; the range leaves room for the labels and no tick lies above an AUROC of one.
        lower, upper = label_room(anchors, low - 0.04 * span, high + 0.04 * span, LABEL_INCHES)
        ax.set_xlim(lower, upper)
        ax.set_xticks([tick for tick in MaxNLocator(4, steps=[1, 2, 5, 10]).tick_values(lower, min(upper, 1.0)) if lower <= tick <= min(upper, 1.0)])
        ax.set_title(policy["task_labels"][task], fontsize=9, color=color, pad=4)
        ax.set_xlabel("AUROC")
        V.clean_axis(ax, grid="x", categorical="y")
    first = axes[0]
    first.set_yticks([])
    first.set_ylim(len(order) - 0.5, -0.5)
    # Every row names its effect and, below in grey, the two sides of its comparison.
    for row, effect in enumerate(order):
        first.text(-0.04, row - 0.11, policy["effect_labels"][effect], transform=first.get_yaxis_transform(), ha="right", va="center", fontsize=8)
        first.text(-0.04, row + 0.19, " → ".join(policy["effect_ends"][effect]), transform=first.get_yaxis_transform(), ha="right", va="center",
                   fontsize=8, color=V.MUTED)
    V.panel(canvas, first, "C", "Each source of inflation, measured on its own", across=True)
    grey = model_shades("#7A7A7A", policy)
    handles = [Line2D([], [], color=grey[model], linewidth=1.6, marker=">", markersize=6, markevery=[1], markerfacecolor=grey[model],
                      markeredgecolor=grey[model]) for model in models]
    handles += [Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="#7A7A7A", markersize=6),
                Line2D([], [], marker="o", linestyle="", markerfacecolor="#7A7A7A", markeredgecolor="white", markersize=6.5),
                Patch(color="#7A7A7A", alpha=0.25)]
    canvas.legend(handles, [policy["model_labels"][model] for model in models] + ["First AUROC", "Second AUROC", "95% interval"],
                  loc="outside lower center", ncols=5, columnspacing=1.4, handletextpad=0.5, handlelength=2.0)


def shortcut(canvas, curves, pooled, policy):
    """The ROC curve of the study-only predictor over the pooled test records of every repeat of the naive random split (thin) and the median curve of
    the repeats (bold) for every task; the diagonal is what it scores inside any one study."""
    ax = canvas.subplots()
    grid = np.linspace(0, 1, 201)
    handles, labels = [], []
    for task in policy["task_order"]:
        color = policy["task_colors"][task]
        own = curves[curves["task"] == task]
        for _, part in own.groupby("fold", sort=True):
            ax.plot(part["fpr"], part["tpr"], color=color, alpha=0.45, linewidth=0.9)
        middle = np.median(curve_rows(own, "fold", grid), axis=0)
        middle[0] = 0.0
        ax.plot(grid, middle, color=darken(color, 0.15), linewidth=2.4)
        handles.append(Line2D([], [], color=darken(color, 0.15), linewidth=2.4))
        labels.append(f"{policy['task_labels'][task]}  {pooled[task]:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1.1, zorder=0)
    handles += [Line2D([], [], color="#888888", alpha=0.7, linewidth=0.9), Line2D([], [], linestyle="--", color="#888888", linewidth=1.1)]
    labels += [f"One random split ({curves.groupby('task')['fold'].nunique().max()} per task)", "Within one study (0.5)"]
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_box_aspect(1)
    ax.set_xlabel("False-positive rate (1 - specificity)")
    ax.set_ylabel("True-positive rate (sensitivity)")
    V.clean_axis(ax, grid="both")
    V.panel(canvas, ax, "D", "The study alone predicts the pooled class")
    ax.legend(handles, labels, loc="lower right", title="AUROC of the pooled test records", title_fontsize=8)


def identification(canvas, recalls, accuracy, chance, config, seed):
    """Every study platform of the corrected bulk-blood table as a line from the share of its samples predicted as its own study platform from the
    aligned values to that from the corrected values, coloured by its technology, its ends as large as its samples; bold, the mean of the study
    platforms, which is the balanced accuracy."""
    ax = canvas.subplots()
    jitter = np.random.default_rng(seed)
    present = [family for family in config["family_order"] if (recalls["family"] == family).any()]
    for family in present:
        part = recalls[recalls["family"] == family]
        color = config["family_colors"][family]
        # The platforms spread sideways, so that the many shares of one stay apart.
        offset = jitter.uniform(-0.08, 0.08, len(part))
        size = 8 + 1.3 * np.sqrt(part["samples"].to_numpy(dtype=float))
        ax.plot(np.vstack([offset, 1 + offset]), part[["aligned", "corrected"]].to_numpy(dtype=float).T, color=color, alpha=0.35, linewidth=0.9)
        V.bubbles(ax, offset, part["aligned"], color, size=size, alpha=0.65)
        V.bubbles(ax, 1 + offset, part["corrected"], color, size=size, alpha=0.65)
    ax.plot([0, 1], [accuracy["aligned"], accuracy["corrected"]], color=V.INK, linewidth=2.6, marker="o", markersize=7, markerfacecolor=V.INK,
            markeredgecolor="white", markeredgewidth=0.8, zorder=5)
    ax.text(-0.15, accuracy["aligned"], f"{accuracy['aligned']:.3f}", ha="right", va="center", fontsize=8, color=V.INK)
    ax.text(1.15, accuracy["corrected"], f"{accuracy['corrected']:.3f}", ha="left", va="center", fontsize=8, color=V.INK)
    ax.axhline(chance, color=V.MUTED, linestyle=":", linewidth=1.0, zorder=0)
    ax.text(1.48, chance, "chance", ha="right", va="bottom", fontsize=8, color=V.MUTED)
    ax.set_xticks([0, 1], ["Aligned values", "Corrected values"])
    ax.set_xlim(-0.3, 1.5)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Share of its samples predicted as its own study platform")
    V.clean_axis(ax, grid="y", categorical="x")
    V.panel(canvas, ax, "E", "Expression identifies the study platform")
    handles = [Line2D([], [], marker="o", linestyle="-", color=config["family_colors"][family], markerfacecolor=config["family_colors"][family],
                      markeredgecolor="white", markersize=7) for family in present]
    handles.append(Line2D([], [], color=V.INK, linewidth=2.6, marker="o", markersize=7, markerfacecolor=V.INK, markeredgecolor="white"))
    canvas.legend(handles, [f"{config['family_labels'][family]} ({int((recalls['family'] == family).sum())})" for family in present]
                  + ["Mean of the study platforms"], loc="outside lower center", ncols=3, columnspacing=1.0, handletextpad=0.4)


def signatures(canvas, strata, medians, columns, policy, seed):
    """A raincloud for every signature and comparison: above, the density of all its tested strata; below, every stratum by its role and independence;
    on the row, an arrow from the median of all strata (open) to that of the independent strata (filled); right, the two medians and numbers."""
    ax = canvas.subplots()
    clouds = [part["auroc"].to_numpy(dtype=float) for _, part in strata.groupby(["signature", "variant", "comparison"])]
    low = fading_edge(clouds, min(0.4, float(np.floor(strata["auroc"].min() * 20) / 20)), 0.02)
    grid = np.linspace(low, 1.0, 500)
    jitter = np.random.default_rng(seed)
    labels = []
    for row, (signature, variant, comparison) in enumerate(columns):
        own = strata[(strata["signature"] == signature) & (strata["variant"] == variant) & (strata["comparison"] == comparison)]
        color = policy["task_colors"][comparison]
        dark = darken(color, 0.2)
        density = bumps(own["auroc"], grid, 0.02)
        density = density / density.max() * 0.36
        ax.fill_between(grid, row - 0.06, row - 0.06 - density, color=color, alpha=0.25, linewidth=0)
        ax.plot(grid, row - 0.06 - density, color=color, linewidth=1.0)
        spread = row + 0.11 + jitter.uniform(0, 0.25, len(own))
        for role, marker in ROLE_MARKERS.items():
            for independent in (True, False):
                chosen = ((own["derivation_role"] == role) & (own["independent"] == independent)).to_numpy()
                if chosen.any():
                    ax.scatter(own.loc[chosen, "auroc"], spread[chosen], s=26 if marker == "o" else 22, marker=marker,
                               facecolor=color if independent else "white", edgecolor=dark, linewidth=0.9, alpha=0.9, zorder=3)
        everything = float(medians.loc[(signature, variant, comparison, "all"), "auroc_median"])
        apart = float(medians.loc[(signature, variant, comparison, "independent"), "auroc_median"])
        top = row - 0.06 - float(np.interp(everything, grid, density))
        ax.plot([everything, everything], [row, top], color=V.INK, linewidth=1.3, alpha=0.6, zorder=4)
        ax.plot([apart, apart], [row, row + 0.38], color=V.INK, linewidth=1.3, zorder=4)
        if abs(apart - everything) > 0.012:
            ax.annotate("", xy=(apart, row), xytext=(everything, row), arrowprops=dict(arrowstyle="-|>", color=V.INK, linewidth=1.5, shrinkA=3.5,
                                                                                         shrinkB=4, mutation_scale=11), zorder=5)
        ax.scatter([everything], [row], s=34, facecolor="white", edgecolor=V.INK, linewidth=1.1, zorder=6)
        ax.scatter([apart], [row], s=38, facecolor=V.INK, edgecolor="white", linewidth=0.5, zorder=7)
        tested, independent = (int(medians.loc[(signature, variant, comparison, subset), "tested"]) for subset in ("all", "independent"))
        ax.text(1.015, row, f"{everything:.3f} → {apart:.3f}\n{tested} → {independent} strata", transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=8, color=V.INK, linespacing=1.35)
        labels.append(f"{policy['signature_labels'][signature]}\n{policy['comparison_labels'][comparison]}")
    ax.axvline(0.5, color=V.MUTED, linestyle="--", linewidth=0.8, zorder=0)
    ax.set_yticks(range(len(columns)), labels)
    ax.set_ylim(len(columns) - 0.5, -0.5)
    ax.set_xlim(low, 1.012)
    ax.set_xlabel("AUROC within one study platform and tissue")
    V.clean_axis(ax, grid="x", categorical="y")
    V.panel(canvas, ax, "F", "Sepsis signatures re-tested within single studies")
    grey = "#7A7A7A"
    handles = [Patch(color=grey, alpha=0.3), Line2D([], [], marker="o", linestyle="", markerfacecolor=grey, markeredgecolor=grey, markersize=6),
               Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=grey, markersize=6)]
    handles += [Line2D([], [], marker=ROLE_MARKERS[role], linestyle="", markerfacecolor="white", markeredgecolor=grey, markersize=5.5)
                for role in policy["role_labels"]]
    handles.append(Line2D([], [], color=V.INK, linewidth=1.5, marker=">", markersize=5, markevery=[1]))
    canvas.legend(handles, ["All tested strata", "Independent stratum", "Not independent", *policy["role_labels"].values(), "Median, all to independent"],
                  loc="outside lower center", ncols=4, columnspacing=1.2, handletextpad=0.4, handlelength=1.6)


def heatmap(canvas, tests, strata, columns, policy, config):
    """Every stratum (columns) of every signature, comparison and weighting (rows), the bulk-blood strata on the left and those of other material on the
    right: the AUROC of a tested stratum as its shade, the first reason of a skipped one as a hatch, a frame for a derivation study and a circle where
    the stratum is not independent; strips of the material and the technology on top, the key below."""
    ax, key = canvas.subplots(2, 1, height_ratios=[4.7, 1.0])
    cells = tests.set_index(["signature", "variant", "comparison", "study_accession", "platform", "tissue"])
    norm, shades = Normalize(0, 1), V.cmap(diverging=True)
    bulk = policy["material_order"][0]
    # The other material stands one empty column right of the bulk blood.
    place = np.where(strata["material"] == bulk, strata["column"], strata["column"] + int((strata["material"] == bulk).sum()) + 1).astype(int)
    # A row of another weighting names the weighting in place of the comparison, which the row of the primary formula above it names.
    rows = [policy["signature_labels"][signature] + "\n" + (policy["comparison_labels"][comparison] if variant == "primary" else policy["variant_labels"][variant])
            for signature, variant, comparison in columns]
    for x, stratum in zip(place, strata.itertuples(index=False)):
        ax.add_patch(Rectangle((x - 0.5, -2.25), 1, 0.7, facecolor=policy["material_colors"][stratum.material], linewidth=0))
        ax.add_patch(Rectangle((x - 0.5, -1.45), 1, 0.7, facecolor=config["family_colors"][stratum.family], linewidth=0))
        for y, column in enumerate(columns):
            index = (*column, stratum.study_accession, stratum.platform, stratum.tissue)
            if index not in cells.index:
                continue
            entry = cells.loc[index]
            if entry["status"] == "tested":
                ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=shades(norm(float(entry["auroc"]))), linewidth=0))
            else:
                ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=SKIPPED, edgecolor="#B0B0B0", hatch=policy["reason_hatches"][entry["reason"]],
                                       linewidth=0))
            if entry["derivation_role"]:
                ax.add_patch(Rectangle((x - 0.36, y - 0.36), 0.72, 0.72, facecolor="none", edgecolor="#111111", linewidth=1.1,
                                       linestyle="-" if entry["derivation_role"] == "discovery" else (0, (2, 1.2)), zorder=3))
            if not entry["independent"]:
                ax.plot([x], [y], marker="o", markersize=2.4, markerfacecolor="white", markeredgecolor="#111111", markeredgewidth=0.6, zorder=4)
    for low, high in ((0, int((strata["material"] == bulk).sum())), (int((strata["material"] == bulk).sum()) + 1, int(place.max()) + 1)):
        ax.vlines(np.arange(low, high + 1) - 0.5, -0.5, len(columns) - 0.5, color="white", linewidth=0.7, zorder=2)
        ax.hlines(np.arange(len(columns) + 1) - 0.5, low - 0.5, high - 0.5, color="white", linewidth=0.7, zorder=2)
    ax.set_xlim(-0.5, place.max() + 0.5)
    ax.set_ylim(len(columns) - 0.5, -2.3)
    ax.set_yticks([-1.9, -1.1, *range(len(columns))], ["Material", "Technology", *rows])
    ax.set_xticks(place, strata["label"], rotation=90)
    ax.tick_params(length=0)
    for side in ax.spines.values():
        side.set_visible(False)
    V.panel(canvas, ax, "A", "Every stratum of every signature")
    key.set_axis_off()
    materials = [material for material in policy["material_order"] if (strata["material"] == material).any()]
    families = [family for family in config["family_order"] if (strata["family"] == family).any()]
    reasons = [reason for reason in policy["reason_hatches"] if (tests["reason"] == reason).any()]
    groups = [[(Patch(facecolor=policy["material_colors"][material]), material) for material in materials],
              [(Patch(facecolor=config["family_colors"][family]), config["family_labels"][family]) for family in families],
              [(Patch(facecolor=SKIPPED, edgecolor="#B0B0B0", hatch=policy["reason_hatches"][reason], linewidth=0),
                policy["reason_labels"][reason].capitalize()) for reason in reasons],
              [(Patch(facecolor="none", edgecolor="#111111", linewidth=1.1, linestyle=style), label)
               for style, label in zip(("-", (0, (2, 1.2))), policy["role_labels"].values())]
              + [(Line2D([], [], marker="o", linestyle="", markersize=4, markerfacecolor="white", markeredgecolor="#111111", markeredgewidth=0.7),
                  "Not independent")]]
    # One group to a column of the key: the legend fills its columns top to bottom, so every group is padded to the longest.
    height = max(len(group) for group in groups)
    entries = [entry for group in groups for entry in group + [(Patch(facecolor="none", edgecolor="none"), "")] * (height - len(group))]
    key.legend(*zip(*entries), loc="center left", bbox_to_anchor=(0.08, 0.5), ncols=len(groups), columnspacing=2.2, handletextpad=0.5,
               handlelength=1.4, borderaxespad=0)
    bar = canvas.colorbar(ScalarMappable(norm, shades), cax=key.inset_axes([0.66, 0.5, 0.24, 0.13]), orientation="horizontal")
    bar.outline.set_visible(False)
    bar.set_ticks([0, 0.25, 0.5, 0.75, 1], labels=["0", "0.25", "0.5", "0.75", "1"])
    bar.set_label("AUROC within the stratum")


def signature_curves(canvas, curves, medians, primary, policy):
    """The ROC curve of every tested stratum of every primary signature and comparison, independent strata in the colour of the comparison and the others
    in grey, with the median curve of all strata and its middle half."""
    axes = canvas.subplots(1, len(primary), sharey=True)
    grid = np.linspace(0, 1, 201)
    keys = ["study_accession", "platform", "tissue"]
    for ax, (signature, variant, comparison) in zip(axes, primary):
        color = policy["task_colors"][comparison]
        own = curves[(curves["signature"] == signature) & (curves["variant"] == variant) & (curves["comparison"] == comparison)]
        for _, part in own.groupby(keys, sort=False):
            apart = part["independent"].iloc[0]
            ax.plot(part["fpr"], part["tpr"], color=color if apart else "#9A9A9A", alpha=0.5 if apart else 0.55, linewidth=0.8, zorder=2 if apart else 1)
        rows = curve_rows(own, keys, grid)
        middle = np.median(rows, axis=0)
        middle[0] = 0.0
        first, third = np.quantile(rows, [0.25, 0.75], axis=0)
        ax.fill_between(grid, first, third, color=color, alpha=0.16, linewidth=0, zorder=0)
        ax.plot(grid, middle, color=darken(color, 0.3), linewidth=2.2, zorder=3)
        ax.plot([0, 1], [0, 1], linestyle="--", color="#888888", linewidth=1.0, zorder=0)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_box_aspect(1)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
        ax.set_title(f"{policy['signature_labels'][signature]}\n{policy['comparison_labels'][comparison]}", fontsize=9, pad=4)
        V.clean_axis(ax, grid="both")
        counts = {subset: int(medians.loc[(signature, variant, comparison, subset), "tested"]) for subset in ("all", "independent")}
        values = {subset: float(medians.loc[(signature, variant, comparison, subset), "auroc_median"]) for subset in ("all", "independent")}
        ax.legend([Line2D([], [], color=darken(color, 0.3), linewidth=2.2), Line2D([], [], color=color, linewidth=0.9)],
                  [f"All {counts['all']}  {values['all']:.3f}", f"Independent {counts['independent']}  {values['independent']:.3f}"], loc="lower right",
                  handlelength=1.2, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3, frameon=True, facecolor="white", edgecolor="none",
                  framealpha=0.85, borderpad=0.3)
    axes[0].set_ylabel("True-positive rate")
    axes[len(axes) // 2].set_xlabel("False-positive rate")
    V.panel(canvas, axes[0], "B", "Every tested stratum as a ROC curve", across=True)
    canvas.legend([Line2D([], [], color="#4A4A4A", linewidth=0.9), Line2D([], [], color="#9A9A9A", linewidth=0.9),
                   Line2D([], [], color="#4A4A4A", linewidth=2.2), Patch(color="#7A7A7A", alpha=0.25),
                   Line2D([], [], linestyle="--", color="#888888", linewidth=1.0)],
                  ["Independent stratum", "Not independent", "Median curve of all strata", "Middle half of the curves", "Chance"],
                  loc="outside lower center", ncols=5, columnspacing=1.4, handletextpad=0.5, handlelength=1.6)


def main():
    run = V.Run(7)
    config = run.config
    policy = config["benchmark"]
    tasks, columns = policy["task_order"], [tuple(column) for column in policy["signature_columns"]]
    models = policy["model_order"]
    populations = run.table("benchmark/split_populations.parquet")
    leakage = run.table("benchmark/split_leakage.parquet")
    splits = run.table("cohort/benchmark_splits.parquet")
    units = run.table("benchmark/model_units.parquet")
    summary = run.table("benchmark/model_summary.parquet")
    effects = run.table("benchmark/model_inflation_effects.parquet")
    shortcut_folds = run.table("benchmark/shortcut_folds.parquet")
    shortcut_predictions = run.table("benchmark/shortcut_predictions.parquet")
    identification_table = run.table("benchmark/study_identification.parquet")
    predictions = run.table("benchmark/study_identification_predictions.parquet")
    families = run.table("annotation/matrix_index.parquet", ["accession", "platform", "family"])
    tests = run.table("benchmark/signature_tests.parquet")
    scores = run.table("benchmark/signature_scores.parquet")
    signature_summary = run.table("benchmark/signature_summary.parquet")
    leakage_receipt = run.receipt("outputs/benchmark/01.json")
    shortcut_receipt = run.receipt("outputs/benchmark/03.json")
    # The benchmark configuration is part of the code that the quality control of the benchmark stage binds, which the run has verified.
    settings = V.read_json(V.ROOT / "src/07_benchmark/benchmark_config.json")["splits"]

    # A: every test record of every fold, rebuilt from the populations and checked fold by fold against the leakage table, and the share of every kind
    # of twin in every ring, checked against the leakage receipt.
    V.require(set(policy["twin_order"]) == set(policy["twin_colors"]) == set(policy["twin_labels"]) and "none" in policy["twin_order"]
              and list(RINGS) == policy["split_order"], "Every kind of twin needs a colour and a label, and every ring a split (benchmark.twin_order, "
              "twin_colors, twin_labels, split_order)")
    records = test_records(populations, splits, leakage, policy, settings)
    counts = records.groupby(["task", "split", "twin"]).size().unstack(fill_value=0).reindex(columns=policy["twin_order"], fill_value=0).stack()
    counts = counts.rename("records").reset_index()
    counts["total"] = counts.groupby(["task", "split"])["records"].transform("sum")
    counts["share"] = counts["records"] / counts["total"]
    for (task, split), part in counts.groupby(["task", "split"]):
        exposed = int(part.loc[part["twin"] != "none", "records"].sum())
        if split == "patient_safe":
            stated = [entry for key, entry in leakage_receipt["leakage"][task].items() if key.startswith("gts:")]
            V.require(exposed == 0 and int(part["total"].iloc[0]) == sum(entry["test_records"] for entry in stated),
                      f"The patient-safe ring of {task} must hold every test record of the receipt, none with a twin")
        else:
            stated = leakage_receipt["leakage"][task][f"naive:{split}"]
            V.require(int(part["total"].iloc[0]) == stated["test_records"] and exposed == stated["exposed"]
                      and abs(exposed / int(part["total"].iloc[0]) - stated["share_exposed"]) < 1e-6,
                      f"The {split} ring of {task} must give the share of the leakage receipt")
    V.require(records.loc[records["split"] == "patient_safe", ["task", "scheme", "fold"]].drop_duplicates().shape[0]
              == leakage_receipt["patient_safe_folds"]["folds"], "Every patient-safe fold of the receipt must be drawn")
    run.add("Fig7", "twins", records)
    run.add("Fig7", "twin_shares", counts.round(6))

    # B: every scored unit of the two models, with the median and middle half that the clouds draw.
    scored = units[(units["status"] == "scored") & units["model"].isin(models)].copy()
    V.require(set(scored["task"]) == set(tasks) and set(scored["scheme"]) == set(policy["scheme_order"]),
              "The tasks and schemes of the baseline models must be those of the configuration (benchmark.task_order and scheme_order)")
    quantiles = scored.groupby(["task", "model", "scheme"])["auroc"].quantile([0.25, 0.5, 0.75]).unstack()
    quantiles.columns = ["auroc_q1", "auroc_median", "auroc_q3"]
    stated = summary[summary["model"].isin(models)].set_index(["task", "model", "scheme"])
    V.require(quantiles.index.sort_values().equals(stated.index.sort_values())
              and np.allclose(quantiles["auroc_median"], stated.loc[quantiles.index, "auroc_median"], rtol=0, atol=1e-6)
              and (scored.groupby(["task", "model", "scheme"]).size() == stated.loc[quantiles.index, "units_scored"]).all(),
              "The medians of the drawn units must be the medians of the model summary")
    run.add("Fig7", "baselines", scored[["task", "model", "scheme", "fold", "unit", "records", "auroc"]].round(6))
    run.add("Fig7", "baseline_quantiles", quantiles.join(stated["units_scored"]).reset_index().round(6))

    # C: the two ends of every comparison, the band of its interval around the second and its difference.
    effects = effects[effects["model"].isin(models)].copy()
    V.require(set(effects["effect"]) == set(policy["effect_order"])
              and len(effects) == len(tasks) * len(models) * len(policy["effect_order"])
              and effects.loc[effects["effect"] == "total_inflation", "low"].isna().all() and effects[["first", "second"]].notna().all().all()
              and effects.loc[effects["effect"] != "total_inflation", ["low", "high"]].notna().all().all()
              and np.allclose(effects["first"] - effects["second"], effects["estimate"], rtol=0, atol=2e-6),
              "Every task and model must have the total inflation, without an interval, and each effect of benchmark.effect_order once with its "
              "interval, its estimate the first AUROC minus the second")
    effects["band_low"], effects["band_high"] = effects["first"] - effects["high"], effects["first"] - effects["low"]
    run.add("Fig7", "inflation", effects[["task", "model", "effect", "matched", "status", "first", "second", "estimate", "low", "high", "band_low",
                                          "band_high", "studies"]].round(6))

    # D: the study-only score of every pooled test record of every repeat of the naive random split, as one ROC curve per repeat, each with the area of
    # its fold in the shortcut table; within every held-out study the score is one value, an AUROC of one half.
    chosen = shortcut_predictions[shortcut_predictions["scheme"] == "naive_record_split"]
    pooled_folds = shortcut_folds[(shortcut_folds["scheme"] == "naive_record_split") & (shortcut_folds["scoring"] == "pooled")].set_index(["task", "fold"])
    rows = []
    for (task, fold), part in chosen.groupby(["task", "fold"], sort=True):
        fpr, tpr = roc_curve(part.loc[part["label"] == 1, "score"].to_numpy(dtype=float), part.loc[part["label"] == 0, "score"].to_numpy(dtype=float))
        V.require(len(part) == int(pooled_folds.loc[(task, fold), "records"]) and abs(area(fpr, tpr) - float(pooled_folds.loc[(task, fold), "auroc"])) < 1e-6,
                  f"The ROC curve of the study-only score of {fold} ({task}) must enclose the AUROC of the shortcut table")
        rows.append(pd.DataFrame({"task": task, "fold": fold, "point": np.arange(len(fpr)), "fpr": fpr, "tpr": tpr}))
    curves = pd.concat(rows, ignore_index=True)
    within = shortcut_folds[(shortcut_folds["scheme"] == "naive_leave_one_study_out") & (shortcut_folds["scoring"] == "within_study")
                            & (shortcut_folds["status"] == "scored")]
    pooled = {}
    for task in tasks:
        stated = shortcut_receipt["study_only"][task]
        pooled[task] = round(float(pooled_folds.loc[task, "auroc"].mean()), 6)
        V.require(set(curves.loc[curves["task"] == task, "fold"]) == set(pooled_folds.loc[task].index)
                  and len(pooled_folds.loc[task]) == stated["random_split"]["repeats"] and abs(pooled[task] - stated["random_split"]["pooled_auroc_mean"]) < 1e-6
                  and (within.loc[within["task"] == task, "auroc"] == 0.5).all()
                  and int((within["task"] == task).sum()) == stated["held_out_study"]["studies_scored"],
                  f"The study-only scores of {task} must be those of the receipt of the shortcut test")
    run.add("Fig7", "shortcut_curves", curves.round(6))
    run.add("Fig7", "shortcut_scores", pd.concat([pooled_folds.reset_index().assign(scheme="naive_record_split", scoring="pooled"), within])[
        ["task", "scheme", "scoring", "fold", "unit", "records", "auroc"]].round(6))

    # E: the share of every study platform's samples predicted as its own, from the aligned and from the corrected values; their mean is the balanced
    # accuracy of the identification table.
    platforms = predictions[["unit", "study_accession", "platform"]].drop_duplicates().merge(
        families.rename(columns={"accession": "study_accession"}), on=["study_accession", "platform"], how="left", validate="one_to_one")
    V.require(platforms["family"].isin(config["family_order"]).all() and platforms["unit"].is_unique, "Every study platform identified needs its technology")
    V.require(set(predictions["predicted"]) <= set(platforms["unit"]) and set(predictions["table"]) == {"aligned", "corrected"},
              "A predicted study platform must be one of the study platforms identified, from the aligned and the corrected values")
    hits = predictions.assign(hit=predictions["unit"] == predictions["predicted"]).groupby(["unit", "table"])["hit"].agg(["mean", "size"]).unstack("table")
    V.require((hits[("size", "aligned")] == hits[("size", "corrected")]).all(), "Every sample must be predicted from both values")
    recalls = platforms.set_index("unit").assign(aligned=hits[("mean", "aligned")], corrected=hits[("mean", "corrected")],
                                                 samples=hits[("size", "aligned")]).reset_index()
    recalls = recalls.assign(order=[natural(study) for study in recalls["study_accession"]]).sort_values(["order", "platform"]).drop(columns="order")
    accuracy = {}
    for table in ("aligned", "corrected"):
        stored = identification_table[(identification_table["table"] == table) & (identification_table["fold"] == "all")]
        accuracy[table] = float(recalls[table].mean())
        V.require(len(stored) == 1 and abs(accuracy[table] - float(stored["balanced_accuracy"].iloc[0])) < 1e-6
                  and int(stored["units"].iloc[0]) == len(recalls)
                  and abs(accuracy[table] - shortcut_receipt["study_identification"]["balanced_accuracy"][table]) < 1e-6,
                  f"The mean share of own study platform predicted from the {table} values must be the balanced accuracy of the identification table")
    chance = float(identification_table["chance"].iloc[0])
    V.require(abs(chance - 1 / len(recalls)) < 1e-6 and abs(chance - shortcut_receipt["study_identification"]["chance"]) < 1e-6,
              "Chance must be one in the number of study platforms")
    run.add("Fig7", "identification", recalls[["unit", "study_accession", "platform", "family", "samples", "aligned", "corrected"]].round(6))
    run.add("Fig7", "identification_means", pd.DataFrame({"table": list(accuracy), "units": len(recalls), "mean": list(accuracy.values()),
                                                          "chance": chance}).round(6))

    # F: exactly the tested strata of the primary formulas, with the medians of the summary.
    keys = ["signature", "variant", "comparison"]
    V.require(set(map(tuple, tests[keys].drop_duplicates().to_numpy())) == set(columns),
              "Every signature, weighting and comparison of the signature tests needs a row in benchmark.signature_columns, and only those")
    V.require(set(tests.loc[tests["status"] != "tested", "reason"]) <= set(policy["reason_hatches"]) <= set(policy["reason_labels"])
              and set(tests["derivation_role"]) <= set(ROLE_MARKERS) and set(tests["tissue"]) <= set(policy["tissue_labels"]),
              "A reason, derivation role or tissue of the signature tests has no label or hatch")
    primary = [column for column in columns if column[1] == "primary"]
    strata = tests[(tests["status"] == "tested") & (tests["variant"] == "primary")]
    medians = signature_summary[(signature_summary["scope"] == "all_tissues") & (signature_summary["variant"] == "primary")].set_index(keys + ["strata_set"])
    for column in primary:
        own = strata[(strata[keys] == list(column)).all(axis=1)]
        for subset, part in (("all", own), ("independent", own[own["independent"]])):
            stated = medians.loc[(*column, subset)]
            V.require(int(stated["tested"]) == len(part) and abs(float(part["auroc"].median()) - float(stated["auroc_median"])) < 1e-6,
                      f"The median of the drawn strata of {column} ({subset}) must be that of the signature summary")
    run.add("Fig7", "signatures", strata[keys + ["study_accession", "platform", "tissue", "auroc", "independent", "derivation_role", "in_scope_group"]].round(6))
    run.add("Fig7", "signature_medians", medians.reset_index().loc[lambda frame: frame["strata_set"].isin(["all", "independent"]),
                                                                   keys + ["strata_set", "strata", "tested", "auroc_median"]].round(6))

    # SuppFig7 A: every stratum, ordered by material, technology and accession; a study names its platform where it has several, and a study platform
    # its tissue where it has several.
    stratum_keys = ["study_accession", "platform", "tissue"]
    layout = tests[stratum_keys + ["in_scope_group"]].drop_duplicates()
    V.require(not layout.duplicated(stratum_keys).any(), "A stratum must be inside or outside the bulk-blood scope for every signature")
    material = {tissue: group for group, tissues in config["tissue_groups"].items() for tissue in tissues}
    layout = layout.merge(families.rename(columns={"accession": "study_accession"}), on=["study_accession", "platform"], how="left", validate="many_to_one")
    layout["material"] = np.where(layout["in_scope_group"], policy["material_order"][0], layout["tissue"].map(material))
    V.require(layout["family"].isin(config["family_order"]).all() and layout["material"].isin(policy["material_order"]).all()
              and set(policy["material_order"]) == set(policy["material_colors"]), "Every stratum needs its technology and a material with a colour")
    several_platforms = layout.groupby("study_accession")["platform"].transform("nunique") > 1
    several_tissues = layout.groupby(["study_accession", "platform"])["tissue"].transform("nunique") > 1
    layout["label"] = [study + (f" {platform}" if platforms_ else "") + (f" {policy['tissue_labels'][tissue]}" if tissues_ else "")
                       for study, platform, tissue, platforms_, tissues_ in zip(layout["study_accession"], layout["platform"], layout["tissue"],
                                                                               several_platforms, several_tissues)]
    V.require(layout["label"].is_unique, "Every stratum needs a label of its own")
    rank = {family: position for position, family in enumerate(config["family_order"])}
    layout = layout.assign(block=layout["material"] != policy["material_order"][0], rank=layout["family"].map(rank),
                           order=[natural(study) for study in layout["study_accession"]])
    layout = layout.sort_values(["block", "rank", "order", "platform", "tissue"]).reset_index(drop=True)
    layout["column"] = layout.groupby("block").cumcount()
    run.add("SuppFig7", "heatmap", tests[keys + stratum_keys + ["value_scale", "positives", "negatives", "status", "reason", "auroc", "auroc_low",
                                                                "auroc_high", "derivation_role", "samples_shared_with_derivation", "repeat_deposits",
                                                                "independent", "in_scope_group"]].merge(
        layout[stratum_keys + ["family", "material", "label", "column"]], on=stratum_keys, validate="many_to_one").round(6))

    # SuppFig7 B: the ROC curve of every tested stratum of the primary formulas, from the stored scores, which carry the direction of the signature;
    # every curve must enclose the AUROC of its test.
    rows = []
    tested = tests[tests["status"] == "tested"].set_index(keys + stratum_keys)
    for key, part in scores[scores["variant"] == "primary"].groupby(keys + stratum_keys, sort=True):
        V.require(key in tested.index, f"Scores without a tested stratum: {key}")
        entry = tested.loc[key]
        fpr, tpr = roc_curve(part.loc[part["label"] == 1, "score"].to_numpy(dtype=float), part.loc[part["label"] == 0, "score"].to_numpy(dtype=float))
        V.require(int((part["label"] == 1).sum()) == entry["scored_positives"] and int((part["label"] == 0).sum()) == entry["scored_negatives"]
                  and abs(area(fpr, tpr) - float(entry["auroc"])) < 1e-6, f"The ROC curve of {key} must enclose the AUROC of the signature test")
        rows.append(pd.DataFrame({**dict(zip(keys + stratum_keys, key)), "point": np.arange(len(fpr)), "fpr": fpr, "tpr": tpr,
                                  "independent": bool(entry["independent"])}))
    signature_rocs = pd.concat(rows, ignore_index=True)
    V.require(signature_rocs[keys + stratum_keys].drop_duplicates().shape[0] == len(strata), "Every tested stratum of a primary formula needs its ROC curve")
    run.add("SuppFig7", "roc", signature_rocs.round(6))

    seed = config["seed"]
    # Twenty columns, so that the panels of every row can take two fifths and three fifths of the width.
    run.figure("Fig7", (14, 16.4), dict(nrows=3, ncols=20, height_ratios=[5.9, 5.1, 5.4]), [
        ("twins", (0, slice(0, 8)), (5.6, 5.9), lambda canvas: twins_panel(canvas, counts, policy, config)),
        ("baselines", (0, slice(8, 20)), (8.4, 5.9), lambda canvas: baselines(canvas, scored, policy, seed)),
        ("inflation", (1, slice(0, 12)), (8.4, 5.1), lambda canvas: inflation(canvas, effects, policy, config)),
        ("shortcut", (1, slice(12, 20)), (5.6, 5.1), lambda canvas: shortcut(canvas, curves, pooled, policy)),
        ("identification", (2, slice(0, 8)), (5.6, 5.4), lambda canvas: identification(canvas, recalls, accuracy, chance, config, seed)),
        ("signatures", (2, slice(8, 20)), (8.4, 5.4), lambda canvas: signatures(canvas, strata, medians, primary, policy, seed))])
    run.figure("SuppFig7", (14, 9.8), dict(nrows=2, ncols=1, height_ratios=[6.4, 3.4]), [
        ("heatmap", (0, 0), (14, 6.4), lambda canvas: heatmap(canvas, tests, layout, columns, policy, config)),
        ("roc", (1, 0), (14, 3.4), lambda canvas: signature_curves(canvas, signature_rocs, medians, primary, policy))])
    run.finish()


if __name__ == "__main__":
    main()
