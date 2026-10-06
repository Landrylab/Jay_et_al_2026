#!/usr/bin/env python3
"""Render MLST, phylogeny/ancestry and metadata panels from numeric source tables."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_sources import read_source_table

MISSING_COLOUR = "#D6D6D6"
TABLES = {"samples": "ED2Data1", "cluster_labels": "ED2Data2", "palette": "ED2Data3", "tree_segments": "ED2Data4", "ancestry": "ED2Data5"}


def key(value):
    if pd.isna(value) or str(value).strip() == "":
        return None
    if isinstance(value, (float, np.floating)) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def cluster_colour_key(value, colours):
    return value if value in colours else "C" + value


def require(frame, columns, name):
    missing = set(columns).difference(frame.columns)
    if missing:
        raise ValueError(f"{name} lacks columns: {sorted(missing)}")


def prepare_data(tables):
    samples, labels, palette, segments, ancestry = (tables[name].copy() for name in TABLES)
    require(samples, ["Sample", "Tree tip order", "Cluster", "Published MLST clade", "Extended clade", "Fluconazole status", "Location", "Study"], "Samples")
    pc_columns = []
    for pc in ("PC1", "PC2"):
        candidates = [column for column in samples if column == pc or column.startswith(pc + " (")]
        if len(candidates) != 1:
            raise ValueError(f"Samples must contain exactly one {pc} coordinate column.")
        pc_columns.append(candidates[0])
    for column in ["Sample", "Cluster", "Published MLST clade", "Extended clade", "Fluconazole status", "Location", "Study"]:
        samples[column] = samples[column].map(key)
    if samples.empty or samples.Sample.isna().any() or samples.Sample.duplicated().any() or samples.Cluster.isna().any():
        raise ValueError("Sample identifiers must be unique and clusters nonmissing.")
    for column in ["Tree tip order"] + pc_columns:
        samples[column] = pd.to_numeric(samples[column], errors="raise")
        if not np.isfinite(samples[column].to_numpy(dtype=float)).all():
            raise ValueError(f"Non-finite {column} values.")
    samples = samples.sort_values("Tree tip order").reset_index(drop=True)
    if samples["Tree tip order"].tolist() != list(range(1, len(samples) + 1)):
        raise ValueError("Tree tip order must be a permutation of consecutive integers starting at one.")
    require(labels, ["Cluster", "Display label"], "Cluster labels")
    labels["Cluster"] = labels.Cluster.map(key)
    if labels.Cluster.duplicated().any() or set(labels.Cluster) != set(samples.Cluster) or labels["Display label"].isna().any():
        raise ValueError("Cluster display labels must contain each observed cluster exactly once.")
    require(palette, ["Track", "Category", "Colour (RGB)", "Marker"], "Palette")
    palette["Category"] = palette.Category.map(key)
    if palette.duplicated(["Track", "Category"]).any():
        raise ValueError("Palette track/category pairs must be unique.")
    colours = {}
    for track, group in palette.groupby("Track", sort=False):
        colours[track] = dict(zip(group.Category, group["Colour (RGB)"]))
        for colour in colours[track].values():
            to_rgba(colour)
    for track, values in [
        ("PCA cluster", [cluster_colour_key(c, colours.get("PCA cluster", {})) for c in samples.Cluster]),
        ("MLST / extended clade", samples["Published MLST clade"].tolist() + samples["Extended clade"].tolist()),
        *[(column, samples[column]) for column in ("Fluconazole status", "Location", "Study")],
    ]:
        missing = {value for value in values if pd.notna(value)}.difference(colours.get(track, {}))
        if missing:
            raise ValueError(f"Palette lacks {track} categories: {sorted(missing)}")
    require(segments, ["Start x", "Start y", "End x", "End y"], "Tree segments")
    segments = segments[["Start x", "Start y", "End x", "End y"]].apply(pd.to_numeric, errors="raise")
    if segments.empty or not np.isfinite(segments.to_numpy()).all():
        raise ValueError("Tree segments must contain finite coordinates.")
    if ((segments[["Start y", "End y"]] < 0) | (segments[["Start y", "End y"]] > len(samples))).any().any():
        raise ValueError("Tree segment positions fall outside the supplied tip rows.")
    require(ancestry, ["Sample", "Component", "Ancestry proportion"], "Ancestry")
    ancestry["Sample"] = ancestry.Sample.map(key)
    ancestry["Component"] = ancestry.Component.map(key)
    ancestry["Ancestry proportion"] = pd.to_numeric(ancestry["Ancestry proportion"], errors="raise")
    if ancestry.duplicated(["Sample", "Component"]).any() or set(ancestry.Sample) != set(samples.Sample):
        raise ValueError("Ancestry must have unique sample/component pairs for exactly the plotted samples.")
    components = [c for c in colours.get("ADMIXTURE component", {}) if c in set(ancestry.Component)]
    if set(components) != set(ancestry.Component) or not components:
        raise ValueError("Palette must identify every ancestry component.")
    q = ancestry.pivot(index="Sample", columns="Component", values="Ancestry proportion").reindex(index=samples.Sample, columns=components)
    values = q.to_numpy(dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any() or not np.allclose(values.sum(axis=1), 1, atol=1e-3, rtol=0):
        raise ValueError("Each sample requires all ancestry components in [0,1], summing to one within 0.001.")
    return samples, labels, palette, segments, q, colours, pc_columns


def overlap(samples, column, clusters, categories):
    """Count available metadata and retain NA fractions for unavailable clusters."""
    counts = pd.crosstab(samples.Cluster, samples[column]).reindex(index=clusters, columns=categories, fill_value=0)
    totals = counts.sum(axis=1)
    return counts, counts.div(totals.replace(0, np.nan), axis=0), totals


def render(tables):
    samples, labels, palette, segments, q, colours, pcs = prepare_data(tables)
    clusters = sorted(samples.Cluster.unique(), key=lambda c: (0, float(c)) if c.replace(".", "", 1).isdigit() else (1, c))
    metadata_order = labels.Cluster.tolist()
    mlst_categories = [c for c in colours["MLST / extended clade"] if c in set(samples["Published MLST clade"].dropna())]
    if not mlst_categories:
        mlst_categories = ["No available labels"]
    vertical_scale = max(1, len(clusters) / 24, len(q.columns) / 24)
    width, height = 180, 315 * vertical_scale
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))

    def axis(x, top, w, h):
        ax = fig.add_axes([x / width, 1 - (top + h) * vertical_scale / height, w / width, h * vertical_scale / height])
        ax.tick_params(labelsize=5, length=2, width=.5)
        for spine in ax.spines.values():
            spine.set_linewidth(.5)
        return ax

    cmap = plt.get_cmap("Blues").with_extremes(bad=MISSING_COLOUR)

    def heatmap(ax, column, levels, order, row_labels, title):
        _, fractions, _ = overlap(samples, column, order, levels)
        im = ax.imshow(np.ma.masked_invalid(fractions.to_numpy(dtype=float)), cmap=cmap, vmin=0, vmax=1, aspect="auto", interpolation="none")
        ax.set_xticks(range(len(levels)), levels, rotation=45 if len(levels) > 3 else 0, ha="right" if len(levels) > 3 else "center")
        ax.set_yticks(range(len(order)), row_labels)
        ax.tick_params(labelsize=5, length=0)
        ax.set_title(title, fontsize=6, loc="left", fontweight="bold")
        position = ax.get_position()
        cax = fig.add_axes([position.x1 + .025, position.y0 + position.height * .12, .011, position.height * .7])
        bar = fig.colorbar(im, cax=cax, ticks=[0, .5, 1])
        bar.set_label("Fraction of available labels", fontsize=5)
        bar.ax.tick_params(labelsize=5, length=2, width=.5)
        bar.outline.set_linewidth(.4)

    def scatter(ax, column, track, levels, title):
        handles = []
        for i, category in enumerate(levels):
            subset = samples[samples[column].eq(category)]
            if subset.empty:
                continue
            p = palette[palette.Track.eq(track) & palette.Category.eq(category)].iloc[0]
            marker = key(p.Marker) or ("o" if i < 8 else "^")
            ax.scatter(subset[pcs[0]], subset[pcs[1]], s=3.5, color=p["Colour (RGB)"], marker=marker, linewidth=.15, edgecolors="#333333", alpha=.9)
            handles.append(Line2D([], [], linestyle="none", marker=marker, color=p["Colour (RGB)"], markersize=2.5, label=category))
        ax.set_xlabel(pcs[0], fontsize=6)
        ax.set_ylabel(pcs[1], fontsize=6)
        ax.set_title(title, fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(.02, .98, f"n = {samples[column].notna().sum():,}", va="top", transform=ax.transAxes, fontsize=5.5)
        if handles:
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, .5), frameon=False, fontsize=5, title="MLST clade" if track == "MLST / extended clade" else track, title_fontsize=5, handlelength=.8, handletextpad=.3, labelspacing=.2)

    heatmap(axis(24, 13, 42, 40), "Published MLST clade", mlst_categories, clusters, ["PCA cluster " + c for c in clusters], "a  Published MLST clade")
    scatter(axis(24, 72, 42, 38), "Published MLST clade", "MLST / extended clade", mlst_categories, "")

    axes = [axis(x, 13, w, 97) for x, w in [(84, 24), (108, 28), (136, 7), (143, 18), (161, 7)]]
    n = len(samples)
    for ax in axes:
        ax.set_ylim(n, 0)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.tick_params(labelsize=5, length=2, width=.5)
    tree, mlst, cluster_track, admixture, clade_track = axes
    tree.add_collection(LineCollection(segments.to_numpy().reshape(-1, 2, 2), colors="#333333", linewidths=.25))
    tree.set_xlim(segments[["Start x", "End x"]].to_numpy().min() - .5, segments[["Start x", "End x"]].to_numpy().max() + .5)
    tree.axis("off")
    tree.set_title(f"b  Phylogeny\n(n = {n:,})", loc="left", fontsize=5.5)
    mlst.set_title("Published\nMLST", fontsize=5.5)
    mlst.set_xlim(-.5, len(mlst_categories) - .5)
    mlst.set_xticks(range(len(mlst_categories)), mlst_categories, rotation=90)
    lines, line_colours = [], []
    for i, value in enumerate(samples["Published MLST clade"]):
        if pd.notna(value):
            x = mlst_categories.index(value)
            lines.append([(x - .42, i + .5), (x + .42, i + .5)])
            line_colours.append(colours["MLST / extended clade"][value])
    mlst.add_collection(LineCollection(lines, colors=line_colours, linewidths=.6))
    for ax, column, track, title in [(cluster_track, "Cluster", "PCA cluster", "PCA\ncluster"), (clade_track, "Extended clade", "MLST / extended clade", "Extended\nclade")]:
        rgb = [to_rgba(colours[track].get(cluster_colour_key(value, colours[track]) if column == "Cluster" else value, MISSING_COLOUR)) for value in samples[column]]
        ax.imshow(np.array(rgb).reshape(n, 1, 4), origin="upper", extent=(0, 1, n, 0), aspect="auto", interpolation="none")
        ax.set_title(title, fontsize=5.5)
    admixture.set_title(f"ADMIXTURE\nK = {q.shape[1]}", fontsize=5.5)
    admixture.set_xlim(0, 1)
    admixture.set_xticks([0, .5, 1])
    # Saved proportions are normalized only for drawing unit-width stacked bars.
    values = q.to_numpy(dtype=float)
    values = values / values.sum(axis=1, keepdims=True)
    left = np.zeros(n)
    for j, component in enumerate(q.columns):
        right = left + values[:, j]
        polygons = [[(left[i], i), (right[i], i), (right[i], i + 1), (left[i], i + 1)] for i in range(n)]
        admixture.add_collection(PolyCollection(polygons, facecolors=colours["ADMIXTURE component"][component], edgecolors="none"))
        left = right
    # Labels follow contiguous runs, so separated branches are never combined.
    for ax, column, track in [(cluster_track, "Cluster", "PCA cluster"), (clade_track, "Extended clade", "MLST / extended clade")]:
        values_for_track = samples[column].fillna("Other").tolist()
        starts = [0] + [i for i in range(1, n) if values_for_track[i] != values_for_track[i - 1]] + [n]
        for start, stop in zip(starts[:-1], starts[1:]):
            value = values_for_track[start]
            if start:
                ax.axhline(start, color="white", linewidth=.2)
            if stop - start >= n * .018 and value != "Other":
                colour = to_rgba(colours[track][cluster_colour_key(value, colours[track]) if column == "Cluster" else value])
                ink = "white" if np.dot(colour[:3], [.2126, .7152, .0722]) < .45 else "black"
                ax.text(.5, (start + stop) / 2, value, ha="center", va="center", fontsize=6, fontweight="bold", color=ink)
    key_ax = axis(169, 13, 10, 97)
    key_ax.axis("off")
    cluster_key = [Patch(color=colours["PCA cluster"][cluster_colour_key(c, colours["PCA cluster"])], label=cluster_colour_key(c, colours["PCA cluster"])) for c in clusters]
    key_ax.legend(handles=cluster_key, loc="upper left", frameon=False, fontsize=5, title="PCA key", title_fontsize=5, handlelength=.7, handletextpad=.3, labelspacing=.7, borderaxespad=0)
    ancestry_ax = axis(84, 120, 94, 14)
    ancestry_ax.axis("off")
    ancestry_key = [Patch(color=colours["ADMIXTURE component"][c], label=c) for c in q.columns]
    ancestry_ax.legend(handles=ancestry_key, loc="upper left", ncol=min(12, len(q.columns)), frameon=False, fontsize=4.5, title="Ancestry component", title_fontsize=5, handlelength=.6, handletextpad=.2, columnspacing=.55, borderaxespad=0)

    for row, (letter, column, title) in enumerate([("c", "Fluconazole status", "Fluconazole phenotype"), ("d", "Location", "Geographic origin"), ("e", "Study", "Contributing study")], start=1):
        levels = list(colours.get(column, {})) or ["No available labels"]
        heatmap(axis(30, 145 + (row - 1) * 56, 37, 40), column, levels, metadata_order, labels["Display label"].tolist(), f"{letter}  {title}")
        scatter(axis(101, 145 + (row - 1) * 56, 46, 40), column, column, levels, "")
    fig.text(.5, .005, "Grey heatmap cells: no available metadata. PCA scatter plots omit missing labels.", ha="center", fontsize=5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Workbook or directory containing all five source tables.")
    for name, default in TABLES.items():
        flag = name.replace("_", "-")
        parser.add_argument("--" + flag, type=Path, help=f"Individual TSV/CSV or workbook overriding --source for {name}.")
        parser.add_argument("--" + flag + "-table", default=default)
    parser.add_argument("--output", required=True, type=Path, help="Output PDF, SVG, or PNG file.")
    args = parser.parse_args()
    tables = {}
    for name in TABLES:
        path = getattr(args, name) or args.source
        if path is None:
            parser.error(f"Provide --source or --{name.replace('_', '-')}.")
        tables[name] = read_source_table(path, getattr(args, name + "_table"))
    fig = render(tables)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
