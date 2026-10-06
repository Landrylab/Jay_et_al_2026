#!/usr/bin/env python3
"""Plot ADMIXTURE replicate errors without rerunning model fitting."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_sources import read_source_table


def prepare_runs(frame):
    frame = frame.rename(columns={"Replicate": "replicate", "Cross-validation error": "cv_error", "CV_error": "cv_error"}).copy()
    required = ["K", "replicate", "cv_error"]
    if not set(required).issubset(frame.columns):
        raise ValueError("CV data require K, Replicate, Cross-validation error (or K, replicate, cv_error).")
    frame = frame[required]
    for column in required:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if frame.empty or not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("CV data must contain finite values.")
    for column in ("K", "replicate"):
        if ((frame[column] < 1) | (frame[column] % 1 != 0)).any():
            raise ValueError(f"{column} must contain positive integers.")
    if (frame.cv_error < 0).any() or frame.duplicated(["K", "replicate"]).any():
        raise ValueError("Each K/replicate must occur once and CV errors must be nonnegative.")
    return frame.sort_values(["K", "replicate"])


def plot_cv(frame, highlight_k=None, highlight_label=None):
    frame = prepare_runs(frame)
    ks = sorted(frame.K.unique())
    if highlight_k is not None and highlight_k not in ks:
        raise ValueError("Highlighted K is absent from the supplied runs.")
    stats = [boxplot_stats(frame.loc[frame.K.eq(k), "cv_error"], whis=1.5)[0] for k in ks]
    fig, ax = plt.subplots(figsize=(max(7, len(ks) * .24), 3.2), layout="constrained")
    ax.bxp(stats, positions=ks, widths=.6, showfliers=False, patch_artist=True,
           boxprops={"facecolor": "#DCEAF2", "edgecolor": "#0072B2", "linewidth": .7},
           medianprops={"color": "#222222", "linewidth": .8},
           whiskerprops={"color": "#0072B2"}, capprops={"color": "#0072B2"})
    minima = []
    for k in ks:
        values = frame.loc[frame.K.eq(k), "cv_error"].to_numpy()
        offsets = np.linspace(-.17, .17, len(values)) if len(values) > 1 else np.zeros(1)
        ax.scatter(k + offsets, values, s=8, color="#0072B2", alpha=.6, linewidths=0, zorder=3)
        minima.append(values.min())
    ax.plot(ks, minima, "-o", color="#D55E00", markersize=3, linewidth=1, label="Lowest CV error per K")
    if highlight_k is not None:
        ax.axvline(highlight_k, color="#333333", linestyle="--", linewidth=.8)
        ax.scatter([highlight_k], [minima[ks.index(highlight_k)]], s=35, facecolors="none", edgecolors="#222222", zorder=5)
        ax.text(highlight_k - .2, .72, highlight_label or f"K = {highlight_k}", transform=ax.get_xaxis_transform(), va="top", ha="right", fontsize=9)
    ax.set(xlabel="Number of ancestral populations (K)", ylabel="Cross-validation error", xlim=(min(ks) - .7, max(ks) + .7))
    ax.set_xticks(ks)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="XLSX workbook, table directory, or a replicate TSV/CSV.")
    parser.add_argument("--table", default="ED1Data6", help="Named Excel Table or table-directory basename.")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF, SVG, or PNG file.")
    parser.add_argument("--highlight-k", type=int, help="Optional author-selected annotation; no knee is inferred.")
    parser.add_argument("--highlight-label", help="Text for the optional K annotation.")
    args = parser.parse_args()
    fig = plot_cv(read_source_table(args.source, args.table), args.highlight_k, args.highlight_label)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
