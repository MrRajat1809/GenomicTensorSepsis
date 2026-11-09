"""Shared style, verified inputs, saving and provenance for the figures of the resource."""
import argparse
import hashlib
import importlib.metadata
import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as Curve
from matplotlib.text import Text
from scipy.cluster import hierarchy
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
UPSTREAM = ("data_acquisition", "expression_processing", "gene_annotation", "cohort_curation", "technical_validation", "batch_correction")
STAGES = {1: "01_selection_and_cohorts.py", 2: "02_recorded_metadata.py", 3: "03_genes_and_scales.py", 4: "04_technical_validation.py",
          5: "05_quality_flags.py", 6: "06_batch_correction.py"}
FIGURES = {1: ["Fig1", "Fig2", "SuppFig1"], 2: ["Fig3", "SuppFig2"], 3: ["Fig4"], 4: ["Fig5", "SuppFig3", "SuppFig4"], 5: ["SuppFig5"],
           6: ["Fig6", "SuppFig6"]}
INK, MUTED, GRID, MISSING, SPINE = "#222222", "#777777", "#DDDDDD", "#E0E0E0", "#333333"
WARM = ["#FDF4D6", "#FCE082", "#F3B032", "#E97233", "#D74241"]
# The panel letters are the only bold text: large, at the top left of their panel, as in the hand-made figures this style follows.
LETTER, LETTER_SIZE = "panel-letter", 16


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_config(config):
    layout, display = config["layout"], config["display"]
    require(config["version"] == "1.0.0" and config["dpi"] == 300 and 0 < config["minimum_font_size"] <= config["font_size"]
            and isinstance(config["seed"], int) and 0 < layout["minimum_blank_axes_width"] <= 1 and layout["tolerance_points"] >= 0
            and display["sex_split_columns"] >= 1 and display["difference_limit"] > 0 and display["calendar_days"] >= 1 and display["ridge_bins"] >= 20,
            "Unsupported visualization configuration")
    for section in (config, config["correction"]):
        for key, colors in section.items():
            if key.endswith("_colors"):
                require(len(set(colors.values())) == len(colors), f"Two categories share a colour in {key}")


def without_warnings(action, where):
    """Run one drawing step and stop on any warning: a warning of matplotlib means a figure that is not the one the code asks for."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        action()
    require(not caught, f"Warnings while drawing {where}: " + "; ".join(sorted({str(item.message) for item in caught})))


def drawn_texts(fig):
    """The texts that the last draw of the figure showed: visible and not empty, leaving out the axes drawn without axis lines and the tick labels
    outside the view."""
    hidden = set()
    for artist in fig.findobj(lambda artist: not artist.get_visible()):
        hidden.update(artist.findobj(Text))
    for ax in fig.get_axes():
        for axis in (ax.xaxis, ax.yaxis):
            if not ax.axison:
                hidden.update(axis.findobj(Text))
                continue
            scale = axis.get_transform()
            low, high = sorted(scale.transform(np.asarray(axis.get_view_interval(), dtype=float)))
            margin = (high - low) * 1e-9
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not low - margin <= scale.transform(np.array([tick.get_loc()], dtype=float))[0] <= high + margin:
                    hidden.update((tick.label1, tick.label2))
    return [text for text in fig.findobj(Text) if text not in hidden and text.get_visible() and text.get_text().strip()]


def check_layout(fig, where, config):
    """Every text that a reader sees is at least the smallest font size, lies inside the figure and covers no other text, and a diagram drawn
    without axis lines keeps the width of its panel."""
    renderer = fig.canvas.get_renderer()
    layout = config["layout"]
    tolerance = layout["tolerance_points"] * fig.dpi / 72
    frame = fig.bbox
    texts = drawn_texts(fig)
    upright, boxes = [], []
    for text in texts:
        label = text.get_text().replace("\n", " ")[:70]
        require(text.get_fontsize() >= config["minimum_font_size"], f"Text smaller than {config['minimum_font_size']} pt in {where}: {label!r}")
        # The extent of the text itself: that of an annotation also holds its arrow.
        box = Text.get_window_extent(text, renderer)
        require(box.x0 >= frame.x0 - tolerance and box.x1 <= frame.x1 + tolerance and box.y0 >= frame.y0 - tolerance and box.y1 <= frame.y1 + tolerance,
                f"Text outside the figure in {where}: {label!r}")
        # The box of a turned text is wider than its letters, so only upright and vertical texts are compared with each other.
        if text.get_rotation() % 90 == 0:
            upright.append(label)
            boxes.append((box.x0, box.x1, box.y0, box.y1))
    if boxes:
        x0, x1, y0, y1 = (np.array(column)[:, None] for column in zip(*boxes))
        overlap = ((np.minimum(x1, x1.T) - np.maximum(x0, x0.T)) > tolerance) & ((np.minimum(y1, y1.T) - np.maximum(y0, y0.T)) > tolerance)
        first, second = np.nonzero(np.triu(overlap, 1))
        require(not len(first), f"Texts overlap in {where}: " + "; ".join(f"{upright[a]!r} and {upright[b]!r}" for a, b in list(zip(first, second))[:3]))
    parents = {ax: sub for sub in fig.subfigs for ax in sub.axes}
    for ax in fig.get_axes():
        if not ax.axison and ax.get_visible() and (ax.patches or ax.texts or ax.lines or ax.collections):
            share = ax.get_window_extent(renderer).width / parents.get(ax, fig).bbox.width
            require(share >= layout["minimum_blank_axes_width"], f"A diagram of {where} keeps only {share:.0%} of the width of its panel: "
                    f"{ax.get_title(loc='left')!r}")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def record(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return {"sha256": digest.hexdigest(), "size_bytes": Path(path).stat().st_size}


def style(config):
    font_manager.findfont(font_manager.FontProperties(family=config["font_family"]), fallback_to_default=False)
    plt.rcParams.update({"font.family": config["font_family"], "font.size": config["font_size"], "font.weight": "normal", "axes.titleweight": "normal",
                         "axes.labelweight": "normal", "axes.titlesize": 11, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": SPINE, "xtick.color": INK, "ytick.color": INK,
                         "axes.labelpad": 6, "axes.linewidth": 0.9, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3.5,
                         "ytick.major.size": 3.5, "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True, "legend.frameon": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white", "figure.facecolor": "white"})


def clean_axis(ax, grid=None, categorical=None):
    """Draw a dotted grid along one axis only, and leave the categorical axis without spine or ticks."""
    if grid:
        ax.grid(axis=grid, color=GRID, linestyle=":", linewidth=0.8)
    if categorical:
        ax.spines["left" if categorical == "y" else "bottom"].set_visible(False)
        ax.tick_params(axis=categorical, length=0)


def panel(canvas, ax, letter, title, across=False):
    """A short title centred over the plot (across=True: over the whole panel, for a panel of several plots side by side) and, at the top left of the
    panel, its letter in large bold type. The title row is taller than the letter, so the letter stands beside the title and above the plot."""
    if across:
        # Left at its own place, so that the layout keeps a row for it above the plots.
        canvas.suptitle(title, fontsize=plt.rcParams["axes.titlesize"], color=INK)
    else:
        ax.set_title(title, loc="center", color=INK, pad=14)
    if letter:
        canvas.text(0.004, 0.996, letter, ha="left", va="top", fontsize=LETTER_SIZE, fontweight="bold", color="#111111", gid=LETTER)


def blank(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def cmap(diverging=False):
    return plt.get_cmap("RdBu_r") if diverging else LinearSegmentedColormap.from_list("warm", WARM)


def ramp(color):
    """From white to one colour, for a heatmap that keeps the colour of its category."""
    return LinearSegmentedColormap.from_list("ramp", ["#FFFFFF", tint(color, 0.55), color])


def bubbles(ax, x, y, color, size=30, alpha=0.75, **options):
    """Translucent bubbles with a thin white rim, the points of the hand-made PCA and SHAP plots; the rim grows with the bubble."""
    return ax.scatter(x, y, s=size, c=color, alpha=alpha, edgecolors="white", linewidths=0.08 * np.sqrt(size), **options)


def kde(values, grid, bandwidth=None):
    """A Gaussian kernel density on the grid, or None where the values cannot carry one."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 3 or np.ptp(values) == 0:
        return None
    return gaussian_kde(values, bw_method=bandwidth)(grid)


def smooth_density(values, edges, sigma=2.0):
    """The density of many values: a histogram on fixed bins smoothed by a Gaussian of sigma bins, with an area of one."""
    counts, _ = np.histogram(np.asarray(values)[np.isfinite(values)], bins=edges)
    smooth = gaussian_filter1d(counts.astype(float), sigma, mode="constant")
    return smooth / (smooth.sum() * (edges[1] - edges[0])) if smooth.sum() else smooth


def tree(ax, link, orientation="top", color="#555555", linewidth=0.8):
    """Draw a hierarchical clustering tree with its leaves at positions 0, 1, 2 ... and return the order of the leaves."""
    found = hierarchy.dendrogram(link, no_plot=True)
    for xs, ys in zip(found["icoord"], found["dcoord"]):
        xs = (np.asarray(xs) - 5) / 10
        if orientation == "top":
            ax.plot(xs, ys, color=color, linewidth=linewidth)
        else:
            ax.plot(ys, xs, color=color, linewidth=linewidth)
    return found["leaves"]


def ribbon(x0, top0, bottom0, x1, top1, bottom1):
    """A band between two vertical spans, with S-shaped edges."""
    middle = (x0 + x1) / 2
    vertices = [(x0, top0), (middle, top0), (middle, top1), (x1, top1), (x1, bottom1), (middle, bottom1), (middle, bottom0), (x0, bottom0), (x0, top0)]
    codes = [Curve.MOVETO, Curve.CURVE4, Curve.CURVE4, Curve.CURVE4, Curve.LINETO, Curve.CURVE4, Curve.CURVE4, Curve.CURVE4, Curve.CLOSEPOLY]
    return Curve(vertices, codes)


def alluvial(ax, columns, flows, colors, labels, gap=0.035, width=0.028, xs=None, alpha=0.42, inner="right"):
    """Nodes in columns, top to bottom, and ribbons between neighbouring columns, every height in proportion to its count.

    columns: one list of node keys for every column, keys unique over all columns; flows: for every two neighbouring columns, a frame with source, target
    and value and the end whose node colours the ribbons ("source" or "target"); colors and labels by key. A node of a middle column must send on what it
    receives. Returns the position of every node: x, bottom and height."""
    flows = [(frame[frame["value"] > 0], colour_by) for frame, colour_by in flows]
    totals = []
    for position, keys in enumerate(columns):
        if position == 0:
            values = flows[0][0].groupby("source")["value"].sum()
        else:
            values = flows[position - 1][0].groupby("target")["value"].sum()
            if position < len(flows):
                sent = flows[position][0].groupby("source")["value"].sum()
                require(values.reindex(keys, fill_value=0).equals(sent.reindex(keys, fill_value=0)), "A middle node must send on what it receives")
        require(set(values.index) <= set(keys), "Every node of a flow must be placed in its column")
        totals.append(values.reindex(keys, fill_value=0))
    largest = max(len(keys) for keys in columns)
    scale = (1 - gap * (largest - 1)) / max(float(total.sum()) for total in totals)
    xs = list(np.linspace(0, 1, len(columns))) if xs is None else xs
    nodes = {}
    for x, keys, total in zip(xs, columns, totals):
        height = float(total.sum()) * scale + gap * (len(keys) - 1)
        top = 0.5 + height / 2
        for key in keys:
            size = float(total[key]) * scale
            nodes[key] = (x, top - size, size)
            ax.add_patch(Rectangle((x - width / 2, top - size), width, size, facecolor=colors[key], edgecolor="none", zorder=3))
            top -= size + gap
    for position, (frame, colour_by) in enumerate(flows):
        frame = frame[frame["value"] > 0]
        order = {key: rank for rank, key in enumerate(columns[position] + columns[position + 1])}
        leaving = {key: nodes[key][1] + nodes[key][2] for key in columns[position]}
        arriving = {key: nodes[key][1] + nodes[key][2] for key in columns[position + 1]}
        # Sorted by source, then target: a source sends its ribbons in the order of the targets, and a target receives them in the order of the
        # sources, so that ribbons cross as little as the order of the nodes allows.
        for row in frame.assign(first=frame["source"].map(order), second=frame["target"].map(order)).sort_values(["first", "second"]).itertuples(index=False):
            size = float(row.value) * scale
            x0, x1 = nodes[row.source][0] + width / 2, nodes[row.target][0] - width / 2
            path = ribbon(x0, leaving[row.source], leaving[row.source] - size, x1, arriving[row.target], arriving[row.target] - size)
            leaving[row.source] -= size
            arriving[row.target] -= size
            ax.add_patch(PathPatch(path, facecolor=colors[getattr(row, colour_by)], edgecolor="none", alpha=alpha, zorder=2))
    for position, keys in enumerate(columns):
        last, first = position == len(columns) - 1, position == 0
        for key in keys:
            x, bottom, size = nodes[key]
            centre = bottom + size / 2
            text = labels[key]
            if first:
                ax.text(x - width / 2 - 0.012, centre, text, ha="right", va="center", fontsize=8, color=INK)
            elif last or inner == "right":
                ax.text(x + width / 2 + 0.012, centre, text, ha="left", va="center", fontsize=8, color=INK,
                        bbox=None if last else dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.8), zorder=4)
            else:
                ax.text(x - width / 2 - 0.012, centre, text, ha="right", va="center", fontsize=8, color=INK,
                        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="none", alpha=0.8), zorder=4)
    return nodes


def tint(color, amount):
    """Lighten a colour towards white by the given share."""
    red, green, blue = to_rgb(color)
    return (red + (1 - red) * amount, green + (1 - green) * amount, blue + (1 - blue) * amount)


def number(value):
    return f"{int(value):,}"


class Run:
    """One figure script: the verified tables it reads, the numbers it plots, the files it saves and the manifest entry of its stage."""
    def __init__(self, stage):
        parser = argparse.ArgumentParser(description=sys.modules["__main__"].__doc__)
        parser.add_argument("--config", default=str(HERE / "visualization_config.json"), help="Visualization configuration")
        parser.add_argument("--outdir", default="outputs/visualization", help="Figure directory inside outputs/")
        args = parser.parse_args()
        out = Path(args.outdir)
        self.stage, self.config_path = stage, Path(args.config).resolve()
        self.out = (out if out.is_absolute() else ROOT / out).resolve()
        require((ROOT / "outputs") in self.out.parents and self.out.parent == ROOT / "outputs" and self.out.name not in UPSTREAM,
                "Choose a figure directory directly inside outputs/ that is not the folder of another stage")
        self.config = read_json(self.config_path)
        validate_config(self.config)
        style(self.config)
        print(f"[{stage:02d}] Verifying the quality control of the stages v{self.config['version']}", flush=True)
        self.bound, self.upstream, self.receipts = {}, {}, {}
        for name in UPSTREAM:
            report = ROOT / "outputs" / name / "qc.json"
            qc = read_json(report)
            require(qc.get("status") == "PASS" and qc.get("failures") == 0, f"Stage {name} has not passed its quality control")
            for path, expected in {**qc["code"], **qc["dependencies"]}.items():
                require(record(ROOT / path) == expected, f"Changed since the quality control of {name}: {path}")
            for path in qc["dependencies"]:
                self.bound.update(read_json(ROOT / path).get("outputs", {}))
            self.receipts.update(qc["dependencies"])
            self.upstream[report.relative_to(ROOT).as_posix()] = record(report)
        self.inputs, self.rows, self.outputs = {}, [], {}
        self.reset()

    def reset(self):
        """Mark the figures of this script as being redrawn: the QC report stale, the manifest incomplete without the entry of this script, and the
        files of its figures removed, so that a renamed panel or a crash leaves nothing behind that looks complete."""
        (self.out / "panels").mkdir(parents=True, exist_ok=True)
        if (self.out / "qc.json").exists():
            write_json(self.out / "qc.json", dict(version=self.config["version"], status="STALE", reason=f"Stage {self.stage:02d} is being rerun"))
        manifest_path = self.out / "manifest.json"
        if manifest_path.exists():
            manifest = read_json(manifest_path)
            manifest["status"] = "INCOMPLETE"
            manifest.get("stages", {}).pop(str(self.stage), None)
            write_json(manifest_path, manifest)
        for figure in FIGURES[self.stage]:
            for path in [*self.out.glob(f"{figure}.*"), *(self.out / "panels").glob(f"{figure}_*")]:
                path.unlink()

    def source(self, name):
        """Return the path of a stage table after confirming that it is the file its stage receipt binds."""
        relative = f"data/processed/{name}"
        require(relative in self.bound, f"Not the output of a completed stage: {relative}")
        if relative not in self.inputs:
            require(record(ROOT / relative) == self.bound[relative], f"Changed since its stage completed: {relative}")
            self.inputs[relative] = self.bound[relative]
        return ROOT / relative

    def table(self, name, columns=None):
        return pd.read_parquet(self.source(name), columns=columns)

    def receipt(self, name):
        """Return a stage receipt, such as outputs/cohort_curation/01.json, after confirming that the quality control of its stage binds it unchanged."""
        require(name in self.receipts, f"Not a receipt that the quality control of a stage binds: {name}")
        if name not in self.inputs:
            require(record(ROOT / name) == self.receipts[name], f"Changed since its stage completed: {name}")
            self.inputs[name] = self.receipts[name]
        return read_json(ROOT / name)

    def add(self, figure, name, frame):
        """Keep the numbers behind one panel for figure_data.csv."""
        frame = frame.reset_index(drop=True).copy()
        frame.insert(0, "panel", name)
        frame.insert(0, "figure", figure)
        self.rows.append(frame)

    def save(self, fig, figure, name=None):
        """Save one figure as PDF and PNG, PNG last, and accept the files only if saving raised no warning and the layout of the PNG passes the check."""
        where = f"{figure} {name}" if name else figure
        for text in fig.findobj(Text):
            require(text.get_gid() == LETTER or text.get_fontweight() in ("normal", "regular", 400), f"Bold text other than a panel letter in {where}")
        width, height = fig.get_size_inches()
        saved = {}
        for extension in ("pdf", "png"):
            filename = f"panels/{figure}_{name}.{extension}" if name else f"{figure}.{extension}"
            path = self.out / filename
            saved[filename] = (path.with_name(path.name + ".tmp"), path)
            without_warnings(lambda: fig.savefig(saved[filename][0], format=extension, dpi=self.config["dpi"]), where)
        # The figure is drawn at the resolution of the files, so the last draw is the one of the PNG.
        check_layout(fig, where, self.config)
        for filename, (temporary, path) in saved.items():
            temporary.replace(path)
            self.outputs[filename] = dict(record(path), width_inches=float(width), height_inches=float(height), dpi=self.config["dpi"])
        plt.close(fig)

    def figure(self, figure, figsize, grid, panels):
        """Export every panel on its own, then draw the same panels into the assembled figure."""
        require(figure in FIGURES[self.stage] and len({name for name, *_ in panels}) == len(panels), f"Unexpected figure or repeated panel: {figure}")
        for name, _, size, draw in panels:
            fig = plt.figure(figsize=size, layout="constrained", dpi=self.config["dpi"])
            without_warnings(lambda: draw(fig), f"{figure} {name}")
            self.save(fig, figure, name)
        fig = plt.figure(figsize=figsize, layout="constrained", dpi=self.config["dpi"])
        spec = fig.add_gridspec(**grid)
        without_warnings(lambda: [draw(fig.add_subfigure(spec[slot])) for name, slot, _, draw in panels], figure)
        self.save(fig, figure)
        print(f"[{self.stage:02d}] {figure}: {len(panels)} panels, {figsize[0]:g} x {figsize[1]:g} in", flush=True)

    def finish(self):
        own = FIGURES[self.stage]
        data = pd.concat(self.rows, ignore_index=True, sort=False)
        require(set(data["figure"]) == set(own) and {name.split(".")[0] for name in self.outputs if "/" not in name} == set(own),
                f"Stage {self.stage:02d} must save and document exactly {own}")
        data = data.astype(object).where(data.notna(), "").astype(str)
        path = self.out / "figure_data.csv"
        manifest_path = self.out / "manifest.json"
        configuration = record(self.config_path)
        manifest = read_json(manifest_path) if manifest_path.exists() else {}
        if manifest.get("configuration") != configuration:
            # A changed configuration changes every figure, so the entries of all stages are dropped.
            manifest, kept = dict(version=self.config["version"], configuration=configuration, stages={}), None
        else:
            kept = pd.read_csv(path, dtype=str, keep_default_na=False) if path.exists() else None
        if kept is not None:
            data = pd.concat([kept[~kept["figure"].isin(own)], data], ignore_index=True, sort=False).fillna("")
        order = {name: position for position, name in enumerate(name for names in FIGURES.values() for name in names)}
        data = data.sort_values("figure", key=lambda column: column.map(order), kind="stable")
        temporary = path.with_name(path.name + ".tmp")
        data.to_csv(temporary, index=False, lineterminator="\n", encoding="utf-8")
        temporary.replace(path)
        mine = data[data["figure"].isin(own)]
        code = {name.relative_to(ROOT).as_posix(): record(name) for name in (HERE / "visualization_common.py", HERE / STAGES[self.stage])}
        runtime = {"python": sys.version.split()[0]}
        for name in ("matplotlib", "numpy", "pandas", "pyarrow", "scipy"):
            runtime[name] = importlib.metadata.version(name)
        for relative, expected in {**self.inputs, **{name: value for name, value in self.upstream.items()}}.items():
            require(record(ROOT / relative) == expected, f"Input changed during the run: {relative}")
        manifest["status"] = "INCOMPLETE"
        manifest["stages"][str(self.stage)] = dict(
            status="COMPLETE", generated_at_utc=datetime.now(timezone.utc).isoformat(), figures=own, inputs=self.inputs, upstream=self.upstream,
            outputs=dict(sorted(self.outputs.items())), code=code, data_rows=int(len(mine)),
            data_panels={figure: sorted(set(mine.loc[mine["figure"] == figure, "panel"])) for figure in own}, runtime=runtime)
        write_json(manifest_path, manifest)
        print(f"[{self.stage:02d}] Complete: {', '.join(own)}; {len(self.outputs)} files in {self.out.relative_to(ROOT).as_posix()}; "
              f"{len(mine):,} rows of figure data", flush=True)
