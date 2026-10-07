#!/usr/bin/env python3
"""Render grouped-validation and within-group RF figures from source tables."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_sources import read_source_table


STYLE = {"Variant based": ("#0072B2", "o"), "Gene-level variant": ("#D55E00", "s"),
         "Gene based": ("#D55E00", "s"),
         "Gene-level variant based": ("#D55E00", "s"), "Phylogeny based": ("#009E73", "^")}
GROUPED_METRICS = ["Balanced accuracy", "Precision", "Recall", "ROC AUC"]


def require(frame, columns, name):
    missing = set(columns).difference(frame.columns)
    if missing or frame.empty:
        raise ValueError(f"{name}: expected nonempty table with columns {sorted(columns)}; missing {sorted(missing)}")


def text(frame, columns):
    for column in columns:
        if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing identifier in {column}")
        frame[column] = frame[column].astype(str)


def numbers(frame, columns, nullable=()):
    for column in columns:
        raw = frame[column]
        absent = raw.isna() | raw.astype(str).str.strip().isin(["", "NA", "NaN"])
        value = pd.to_numeric(raw.where(~absent), errors="coerce")
        if ((~absent) & value.isna()).any() or (column not in nullable and value.isna().any()) or np.isinf(value.to_numpy(float)).any():
            raise ValueError(f"Invalid or missing numeric values in {column}")
        frame[column] = value


def positive_integers(frame, columns):
    numbers(frame, columns)
    if ((frame[columns] <= 0) | (frame[columns] % 1 != 0)).any().any():
        raise ValueError(f"Expected positive integers in {columns}")


def cohort_counts(cohorts, extra):
    require(cohorts, ["Background", "Isolates", "Resistant isolates", extra], "cohorts")
    text(cohorts, ["Background"])
    positive_integers(cohorts, ["Isolates", extra])
    numbers(cohorts, ["Resistant isolates"])
    if cohorts.Background.duplicated().any():
        raise ValueError("Backgrounds must be unique")
    resistant = cohorts["Resistant isolates"]
    if ((resistant < 0) | (resistant > cohorts.Isolates) | (resistant % 1 != 0)).any():
        raise ValueError("Resistant isolates must be an integer between zero and the group size")


def model_styles(models):
    fallback = plt.get_cmap("tab10")
    return {model: STYLE.get(model, (fallback(i % 10), ["o", "s", "^", "D", "v"][i % 5])) for i, model in enumerate(models)}


def configure_style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5,
                         "font.weight": "bold", "axes.labelweight": "bold",
                         "axes.titleweight": "bold", "axes.linewidth": .6,
                         "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.grid": False})


def save(fig, prefix, dpi):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in ["pdf", "png"]:
        fig.savefig(str(prefix) + "." + extension, dpi=dpi, facecolor="white")
    plt.close(fig)


def grouped_data(cohorts, metrics):
    cohorts, metrics = cohorts.copy(), metrics.copy()
    cohort_counts(cohorts, "Test fold")
    require(metrics, ["Background", "Feature set", "Repeat", *GROUPED_METRICS], "grouped metrics")
    text(metrics, ["Background", "Feature set"])
    positive_integers(metrics, ["Repeat"])
    numbers(metrics, GROUPED_METRICS, GROUPED_METRICS)
    if metrics.duplicated(["Background", "Feature set", "Repeat"]).any():
        raise ValueError("Duplicate background/model/repeat metric row")
    if set(metrics.Background) != set(cohorts.Background):
        raise ValueError("Metric and membership backgrounds differ")
    for column in GROUPED_METRICS:
        if not metrics[column].dropna().between(0, 1).all():
            raise ValueError(f"{column} must be in [0,1] or missing")
    models = metrics["Feature set"].drop_duplicates().tolist()
    repeats = set(metrics.Repeat)
    rows = []
    for background in cohorts.Background:
        for model in models:
            part = metrics.loc[metrics.Background.eq(background) & metrics["Feature set"].eq(model)]
            if set(part.Repeat) != repeats:
                raise ValueError("Each background/model combination must contain every repeat; use missing values for undefined metrics")
            for metric in GROUPED_METRICS:
                values = part[metric].dropna()
                rows.append({"Background": background, "Feature set": model, "Metric": metric,
                             "Mean": values.mean(), "Minimum": values.min(), "Maximum": values.max(),
                             "Defined repeats": len(values)})
    folds = sorted(cohorts["Test fold"].astype(int).unique())
    if len(folds) < 2:
        raise ValueError("Grouped validation requires at least two distinct test folds")
    totals = []
    for fold in folds:
        test = cohorts.loc[cohorts["Test fold"].eq(fold)]
        train = cohorts.loc[~cohorts["Test fold"].eq(fold)]
        totals.append({"Fold": fold, "Training isolates": int(train.Isolates.sum()),
                       "Training resistant": int(train["Resistant isolates"].sum()),
                       "Test isolates": int(test.Isolates.sum()), "Test resistant": int(test["Resistant isolates"].sum())})
    return cohorts, metrics, pd.DataFrame(rows), pd.DataFrame(totals)


def plot_grouped(cohorts, metrics, summary, totals, prefix, dpi=300):
    configure_style()
    groups = cohorts.Background.tolist()
    models = metrics["Feature set"].drop_duplicates().tolist()
    styles = model_styles(models)
    n, folds = len(groups), totals.Fold.tolist()
    width, matrix_height, score_height = 180., 6. * n, max(30., 7. * n + 1)
    score_top = 14 + matrix_height + 40
    height = score_top + score_height + 20
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))

    def label(x, y, value, **kwargs):
        return fig.text(x / width, 1 - y / height, str(value), va="center", fontsize=6, **kwargs)

    def axes(x, y, w, h):
        return fig.add_axes([x / width, 1 - (y + h) / height, w / width, h / height])

    label(1, 3, "a", weight="bold").set_fontsize(8)
    label(4, 10, "Genetic background")
    for x, name in [(42, "n"), (51, "R"), (60, "S")]:
        label(x, 10, name, ha="center")
    matrix = axes(69, 14, 106, matrix_height)
    matrix.set(xlim=(-.5, len(folds) - .5), ylim=(n - .5, -.5))
    matrix.set_axis_off()
    for j, fold in enumerate(folds):
        label(69 + 106 * (j + .5) / len(folds), 10, f"Fold {fold}", ha="center")
    for i, (_, row) in enumerate(cohorts.iterrows()):
        y = 14 + matrix_height * (i + .5) / n
        label(4, y, row.Background)
        for x, value in [(42, row.Isolates), (51, row["Resistant isolates"]), (60, row.Isolates - row["Resistant isolates"])]:
            label(x, y, int(value), ha="center")
        for j, fold in enumerate(folds):
            test = row["Test fold"] == fold
            matrix.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor="#0072B2" if test else "#E7E7E7", edgecolor="white", lw=.6))
            matrix.text(j, i, "Test" if test else "Train", ha="center", va="center", color="white" if test else "black", fontsize=6)
    for y, name, total_col, resistant_col in [(14 + matrix_height + 7, "Training, n (R/S)", "Training isolates", "Training resistant"),
                                              (14 + matrix_height + 14, "Test, n (R/S)", "Test isolates", "Test resistant")]:
        label(4, y, name)
        for j, (_, row) in enumerate(totals.iterrows()):
            number, resistant = int(row[total_col]), int(row[resistant_col])
            label(69 + 106 * (j + .5) / len(folds), y, f"{number} ({resistant}/{number - resistant})", ha="center").set_fontsize(5.7)
    label(4, score_top - 6, "Genetic background")
    offsets = np.linspace(-.23, .23, len(models)) if len(models) > 1 else [0]
    for k, metric in enumerate(GROUPED_METRICS):
        x = 29 + 37.5 * k
        label(1 if not k else x - 2, score_top - 12, chr(98 + k), weight="bold").set_fontsize(8)
        label(x + 16.75, score_top - 6, metric, ha="center").set_fontsize(6.1)
        ax = axes(x, score_top, 33.5, score_height)
        ax.set(xlim=(-.06, 1.06), ylim=(n - .35, -.65), xticks=[0, .5, 1], yticks=range(n), yticklabels=[])
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", length=2, width=.5, pad=2)
        ax.spines[["top", "right", "left"]].set_visible(False)
        if metric in ["Balanced accuracy", "ROC AUC"]:
            ax.axvline(.5, color="#666666", lw=.55, ls=(0, (2, 2)), zorder=0)
        for i, background in enumerate(groups):
            if not k:
                label(4, score_top + score_height * (i + .65) / (n + .3), background)
            defined = False
            for j, model in enumerate(models):
                row = summary.loc[summary.Background.eq(background) & summary["Feature set"].eq(model) & summary.Metric.eq(metric)].iloc[0]
                if pd.isna(row.Mean):
                    continue
                defined = True
                color, marker = styles[model]
                errors = np.maximum(0., [[row.Mean - row.Minimum], [row.Maximum - row.Mean]])
                ax.errorbar(row.Mean, i + offsets[j], xerr=errors,
                            color=color, marker=marker, markersize=3.2, markeredgecolor="white", markeredgewidth=.25,
                            elinewidth=.8, capsize=1.5, capthick=.65, linestyle="none")
            if not defined:
                ax.text(.5, i, "NA", ha="center", va="center", color="#555555", fontsize=6)
    handles = [Line2D([], [], marker=styles[m][1], color=styles[m][0], lw=.8, markersize=3.5, label=m) for m in models]
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(.55, 11 / height), ncol=min(3, len(models)), frameon=False, fontsize=6.1, handlelength=1.4)
    save(fig, prefix, dpi)
    summary.to_csv(str(prefix) + ".group_metrics.tsv", sep="\t", index=False, na_rep="NA", float_format="%.17g")
    totals.to_csv(str(prefix) + ".fold_composition.tsv", sep="\t", index=False)
    report = {"backgrounds": len(groups), "isolates": int(cohorts.Isolates.sum()), "folds": len(folds),
              "repeats": int(metrics.Repeat.nunique()), "feature_sets": models,
              "metric_rows": len(metrics), "dimensions_mm": [width, height],
              "whiskers": "Minimum–maximum of defined repeat values"}
    Path(str(prefix) + ".summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def within_data(cohorts, metrics, importance, chromosomes, highlights=None):
    cohorts, metrics, importance, chromosomes = [x.copy() for x in (cohorts, metrics, importance, chromosomes)]
    cohort_counts(cohorts, "Outer folds")
    require(metrics, ["Background", "Feature set", "Repeat", "Fold", "Metric", "Value"], "within-group metrics")
    require(importance, ["Background", "Variant ID", "Chromosome", "Position (bp)", "Mean Gini importance"], "variant importance")
    require(chromosomes, ["Chromosome", "Length (bp)"], "chromosomes")
    text(metrics, ["Background", "Feature set", "Metric"])
    text(importance, ["Background", "Variant ID", "Chromosome"])
    text(chromosomes, ["Chromosome"])
    positive_integers(metrics, ["Repeat", "Fold"])
    numbers(metrics, ["Value"], ["Value"])
    positive_integers(chromosomes, ["Length (bp)"])
    positive_integers(importance, ["Position (bp)"])
    numbers(importance, ["Mean Gini importance"])
    if set(metrics.Background) != set(cohorts.Background) or set(importance.Background) != set(cohorts.Background):
        raise ValueError("Cohort, metric, and importance backgrounds must agree")
    if metrics.duplicated(["Background", "Feature set", "Repeat", "Fold", "Metric"]).any():
        raise ValueError("Duplicated within-group fold metric")
    if not metrics.Value.dropna().between(0, 1).all():
        raise ValueError("Fold scores must be in [0,1] or missing")
    models = metrics["Feature set"].drop_duplicates().tolist()
    metric_names = metrics.Metric.drop_duplicates().tolist()
    for _, cohort in cohorts.iterrows():
        part = metrics.loc[metrics.Background.eq(cohort.Background)]
        all_keys = set(map(tuple, part[["Repeat", "Fold"]].to_numpy()))
        for repeat in part.Repeat.unique():
            if part.loc[part.Repeat.eq(repeat), "Fold"].nunique() != cohort["Outer folds"]:
                raise ValueError("Metric fold counts disagree with the cohort table")
        for model in models:
            for metric in metric_names:
                selected = part.loc[part["Feature set"].eq(model) & part.Metric.eq(metric)]
                if set(map(tuple, selected[["Repeat", "Fold"]].to_numpy())) != all_keys:
                    raise ValueError("Feature sets and metrics must cover identical folds within each background")
    if chromosomes.Chromosome.duplicated().any() or importance.duplicated(["Background", "Variant ID"]).any():
        raise ValueError("Duplicate chromosome or background/variant identifier")
    lengths = chromosomes.set_index("Chromosome")["Length (bp)"]
    if not set(importance.Chromosome).issubset(lengths.index):
        raise ValueError("An importance chromosome is absent from the chromosome table")
    if (importance["Position (bp)"] > importance.Chromosome.map(lengths)).any() or (importance["Mean Gini importance"] < 0).any():
        raise ValueError("Variant coordinate outside the reference or negative importance")
    if highlights is not None and not highlights.empty:
        highlights = highlights.copy()
        require(highlights, ["Background", "Variant ID", "Label"], "highlights")
        text(highlights, ["Background", "Variant ID", "Label"])
        if highlights.duplicated(["Background", "Variant ID"]).any():
            raise ValueError("Duplicated variant highlight")
        for optional in ["Window start (bp)", "Window end (bp)", "Label dx (pt)", "Label dy (pt)"]:
            if optional in highlights:
                numbers(highlights, [optional], [optional])
        if "Label alignment" in highlights:
            alignments = highlights["Label alignment"].dropna()
            if not alignments.isin(["left", "right", "center"]).all():
                raise ValueError("Label alignment must be left, right, or center")
        for _, peak in highlights.iterrows():
            point = importance.loc[importance.Background.eq(peak.Background) & importance["Variant ID"].eq(peak["Variant ID"])]
            if len(point) != 1:
                raise ValueError(f"Highlighted variant {peak['Variant ID']} is absent from its background")
            lo, hi = peak.get("Window start (bp)"), peak.get("Window end (bp)")
            if pd.notna(lo) != pd.notna(hi):
                raise ValueError("A highlight window needs both endpoints")
            if pd.notna(lo) and (lo < 0 or hi <= lo or hi > lengths.loc[point.Chromosome.iloc[0]]):
                raise ValueError("Invalid highlight display window")
    else:
        highlights = pd.DataFrame(columns=["Background", "Variant ID", "Label"])
    return cohorts, metrics, importance, chromosomes, highlights


def plot_within(cohorts, metrics, importance, chromosomes, highlights, prefix, dpi=300, chromosome_gap=150000, seed=42):
    configure_style()
    groups = cohorts.Background.tolist()
    models = metrics["Feature set"].drop_duplicates().tolist()
    metric_names = metrics.Metric.drop_duplicates().tolist()
    styles = model_styles(models)
    n = len(groups)
    height = 35 + 51.25 * n
    fig = plt.figure(figsize=(180 / 25.4, height / 25.4))
    grid = fig.add_gridspec(n, 2, width_ratios=[1, 1.45], left=.080, right=.976, top=.944, bottom=.086, hspace=.52, wspace=.37)
    lengths = chromosomes.set_index("Chromosome")["Length (bp)"]
    offsets = (lengths + chromosome_gap).cumsum().shift(fill_value=0)
    end = offsets.iloc[-1] + lengths.iloc[-1]
    rng = np.random.default_rng(seed)
    box_rows, panel_rows = [], []
    for row, background in enumerate(groups):
        ax = fig.add_subplot(grid[row, 0])
        data = metrics.loc[metrics.Background.eq(background)]
        positions = np.linspace(-.21, .21, len(models)) if len(models) > 1 else [0]
        for x, metric in enumerate(metric_names):
            for j, model in enumerate(models):
                values = data.loc[data.Metric.eq(metric) & data["Feature set"].eq(model), "Value"].dropna().to_numpy(float)
                centre = x + positions[j]
                if not len(values):
                    ax.text(centre, .5, "NA", ha="center", va="center", fontsize=5.5, rotation=90)
                    box_rows.append({"Background": background, "Feature set": model, "Metric": metric, "Defined folds": 0})
                    continue
                stats = boxplot_stats(values, whis=1.5)[0]
                color = styles[model][0]
                ax.bxp([stats], positions=[centre], widths=.68 / len(models), patch_artist=True, manage_ticks=False, showfliers=False,
                       boxprops={"facecolor": color, "edgecolor": "#404040", "linewidth": .6},
                       medianprops={"color": "black", "linewidth": .9}, whiskerprops={"color": "#555555", "linewidth": .6}, capprops={"color": "#555555", "linewidth": .6})
                ax.scatter(centre + rng.uniform(-.08, .08, len(values)), values, color="black", s=3, alpha=.35, linewidths=0, zorder=4)
                box_rows.append({"Background": background, "Feature set": model, "Metric": metric, "Defined folds": len(values),
                                 "Q1": stats["q1"], "Median": stats["med"], "Q3": stats["q3"],
                                 "Lower whisker": stats["whislo"], "Upper whisker": stats["whishi"], "Outliers": len(stats["fliers"])})
        labels = [m.replace("Balanced accuracy", "Balanced\naccuracy") for m in metric_names]
        ax.set(xticks=range(len(metric_names)), xticklabels=labels, xlim=(-.5, len(metric_names) - .5),
               ylim=(-.015, 1.025), ylabel="Held-out fold score", title=background)
        ax.set_yticks([0, .25, .5, .75, 1])
        for x, metric in enumerate(metric_names):
            if metric in ["Balanced accuracy", "ROC AUC"]:
                ax.plot([x - .43, x + .43], [.5, .5], color="#777777", ls=(0, (3, 2)), lw=.65, zorder=1)
        right = fig.add_subplot(grid[row, 1])
        values = importance.loc[importance.Background.eq(background)]
        for i, chrom in enumerate(lengths.index):
            part = values.loc[values.Chromosome.eq(chrom)]
            right.scatter(part["Position (bp)"] + offsets[chrom], part["Mean Gini importance"], s=2, color="#0072B2" if i % 2 == 0 else "#999999", linewidths=0, alpha=.85, rasterized=False)
        peaks = highlights.loc[highlights.Background.eq(background)]
        maximum = float(values["Mean Gini importance"].max())
        span = maximum if maximum > 0 else .001
        right.set(xlim=(-100000, end + 150000), ylim=(-span * .035, span * (1.4 if len(peaks) else 1.1)),
                  xlabel="Chromosome", ylabel="Mean RF Gini importance",
                  xticks=offsets.to_numpy() + lengths.to_numpy() / 2, xticklabels=lengths.index)
        right.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))
        formatter = ScalarFormatter(useMathText=False)
        formatter.set_scientific(False)
        right.yaxis.set_major_formatter(formatter)
        for j, (_, peak) in enumerate(peaks.iterrows()):
            point = values.loc[values["Variant ID"].eq(peak["Variant ID"])].iloc[0]
            offset = offsets[point.Chromosome]
            px, py = point["Position (bp)"] + offset, point["Mean Gini importance"]
            lo, hi = peak.get("Window start (bp)"), peak.get("Window end (bp)")
            if pd.notna(lo):
                region = values.loc[values.Chromosome.eq(point.Chromosome) & values["Position (bp)"].between(lo, hi)]
                top = max(py, region["Mean Gini importance"].max() if len(region) else py) + span * .045
                right.add_patch(Rectangle((offset + lo, 0), hi - lo, top, fill=False, edgecolor="#333333", lw=.8, ls=(0, (1.2, 1.8))))
            right.scatter([px], [py], s=19, color="#D55E00", edgecolor="black", linewidth=.4, zorder=5)
            dx, dy = peak.get("Label dx (pt)"), peak.get("Label dy (pt)")
            dx = float(dx) if pd.notna(dx) else -12
            dy = float(dy) if pd.notna(dy) else 20 + 18 * (j % 2)
            ha = peak.get("Label alignment", "right")
            if pd.isna(ha):
                ha = "right"
            if ha not in ["left", "right", "center"]:
                raise ValueError("Label alignment must be left, right, or center")
            label = str(peak.Label).replace(" ", "\n", 1)
            right.annotate(label, (px, py), xytext=(dx, dy), textcoords="offset points", ha=ha, va="bottom", fontsize=5.6,
                           arrowprops={"arrowstyle": "-", "color": "#333333", "linewidth": .55, "shrinkA": 1.5, "shrinkB": 3}, zorder=6)
        for target, letter, x in [(ax, chr(97 + row), .017), (right, chr(97 + n + row), .456)]:
            target.spines[["top", "right"]].set_visible(False)
            target.tick_params(length=2.5, width=.6)
            fig.text(x, target.get_position().y1 + .011, letter, fontsize=8, weight="bold", va="bottom")
        panel_rows.append({"Background": background, "variant_points": len(values), "highlighted_variants": len(peaks),
                           "metric_values": int(data.Value.notna().sum()), "total_mean_importance": float(values["Mean Gini importance"].sum())})
    handles = [Patch(facecolor=styles[m][0], edgecolor="#404040", linewidth=.6, label=m) for m in models]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.515, .020), frameon=False, ncol=min(3, len(models)), fontsize=7, handlelength=1.6)
    save(fig, prefix, dpi)
    pd.DataFrame(box_rows).to_csv(str(prefix) + ".box_summary.tsv", sep="\t", index=False, na_rep="NA", float_format="%.17g")
    report = {"backgrounds": len(groups), "feature_sets": models, "metric_names": metric_names,
              "dimensions_mm": [180, height], "chromosome_gap_bp": chromosome_gap, "panels": panel_rows,
              "importance": "Supplied mean Gini importance",
              "boxplot": "Linear-interpolation quartiles; whiskers to observations within 1.5 IQR; every defined fold shown"}
    Path(str(prefix) + ".summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    grouped = commands.add_parser("grouped", help="Held-out background membership and metrics")
    within = commands.add_parser("within", help="Within-background fold metrics and variant importance")
    for sub in (grouped, within):
        sub.add_argument("--source", type=Path, help="Source workbook or directory of named TSV/CSV tables")
        sub.add_argument("--cohorts", type=Path)
        sub.add_argument("--metrics", type=Path)
        sub.add_argument("--output-prefix", type=Path, required=True)
        sub.add_argument("--dpi", type=int, default=300)
    within.add_argument("--importance", type=Path)
    within.add_argument("--chromosomes", type=Path)
    within.add_argument("--highlights", type=Path)
    within.add_argument("--no-highlights", action="store_true")
    within.add_argument("--chromosome-gap", type=float, default=150000)
    within.add_argument("--jitter-seed", type=int, default=42)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("DPI must be positive")

    def read(argument, name):
        path = getattr(args, argument) or args.source
        if path is None:
            parser.error(f"Provide --{argument} or --source")
        return read_source_table(path, name)

    if args.command == "grouped":
        checked = grouped_data(read("cohorts", "ED8_PanelA"), read("metrics", "ED8_PanelsBE"))
        report = plot_grouped(*checked, args.output_prefix, args.dpi)
    else:
        if not np.isfinite(args.chromosome_gap) or args.chromosome_gap < 0:
            parser.error("Chromosome gap must be finite and nonnegative")
        peaks = None if args.no_highlights else read("highlights", "ED9_Highlights")
        checked = within_data(read("cohorts", "ED9_Cohorts"), read("metrics", "ED9_PanelsAD"),
                              read("importance", "ED9_PanelsEH"), read("chromosomes", "ED9_Chromosomes"), peaks)
        report = plot_within(*checked, args.output_prefix, args.dpi, args.chromosome_gap, args.jitter_seed)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
