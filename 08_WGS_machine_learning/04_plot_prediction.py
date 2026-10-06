#!/usr/bin/env python3
"""Summarize held-out random-forest predictions and plot validation performance."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import balanced_accuracy_score, precision_score, recall_score, roc_curve


INPUTS = {"variant": ("Variant based", "#0072B2"), "gene": ("Gene-level variant based", "#D55E00"), "pc": ("Phylogeny based", "#009E73")}
DESIGNS = {"stratified": "Stratified CV", "clade_grouped": "Clade-grouped CV"}
METRICS = {"balanced_accuracy": "Balanced accuracy", "precision": "Precision", "recall": "Recall"}


def load_results(root):
    metrics, predictions = [], []
    for design in DESIGNS:
        for input_key in INPUTS:
            folder = root / f"{design}_{input_key}"
            if not folder.exists() and input_key == "variant":
                folder = root / f"{design}_snp"
            metrics_path = folder / "cv_fold_metrics.tsv" if folder.exists() else root / f"{design}_{input_key}_cv_fold_metrics.tsv"
            predictions_path = folder / "cv_fold_predictions.tsv" if folder.exists() else root / f"{design}_{input_key}_cv_fold_predictions.tsv"
            m = pd.read_csv(metrics_path, sep="\t")
            p = pd.read_csv(predictions_path, sep="\t")
            if m.duplicated(["repeat", "fold"]).any() or p.duplicated(["repeat", "sample"]).any():
                raise ValueError(f"Duplicate held-out fold or sample in {folder}")
            if set(map(tuple, m[["repeat", "fold"]].values)) != set(map(tuple, p[["repeat", "fold"]].drop_duplicates().values)):
                raise ValueError(f"Metrics and predictions use different folds in {folder}")
            for data in [m, p]:
                data["input_key"], data["design"] = input_key, design
                data["analysis_id"] = f"{design}_{input_key}"
            metrics.append(m)
            predictions.append(p)
    return pd.concat(metrics, ignore_index=True), pd.concat(predictions, ignore_index=True)


def roc_summary(predictions):
    grid = np.linspace(0, 1, 501)
    rows = []
    for (design, input_key), data in predictions.groupby(["design", "input_key"], sort=False):
        curves = []
        for _, fold in data.groupby(["repeat", "fold"], sort=True):
            labels = fold.true_label.map({"S": 0, "R": 1})
            if labels.isna().any() or labels.nunique() != 2:
                raise ValueError("ROC curves require two known classes in every held-out fold")
            probabilities = fold.probability_R.to_numpy(float)
            if not np.isfinite(probabilities).all() or ((probabilities < 0) | (probabilities > 1)).any():
                raise ValueError("Resistance probabilities must be finite and between zero and one")
            fpr, tpr, _ = roc_curve(labels, probabilities)
            curve = np.interp(grid, fpr, tpr)
            curve[0], curve[-1] = 0, 1
            curves.append(curve)
        values = np.vstack(curves)
        for fpr, mean, sd in zip(grid, values.mean(axis=0), values.std(axis=0, ddof=1)):
            rows.append({"design": design, "input_key": input_key, "fpr": fpr, "mean_tpr": mean, "std_tpr": sd})
    return pd.DataFrame(rows)


def group_scores(predictions, metadata, sample_column="sample", group_column="cluster"):
    if metadata[sample_column].duplicated().any():
        raise ValueError("Metadata sample IDs must be unique")
    lookup = metadata.set_index(sample_column)[group_column]
    data = predictions.loc[predictions.design.eq("clade_grouped")].copy()
    data["cluster"] = data["sample"].map(lookup)
    if data.cluster.isna().any():
        raise ValueError("Every held-out sample requires a genetic-group assignment")
    rows = []
    for (input_key, repeat, cluster), part in data.groupby(["input_key", "repeat", "cluster"], sort=True):
        if part.fold.nunique() != 1:
            raise ValueError("Members of a grouped-validation cluster must share one outer fold")
        y = part.true_label.map({"S": 0, "R": 1}).to_numpy()
        pred = part.predicted_label.map({"S": 0, "R": 1}).to_numpy()
        if pd.isna(y).any() or pd.isna(pred).any():
            raise ValueError("Predicted and observed labels must be R or S")
        rows.append({"input_key": input_key, "repeat": repeat, "cluster": cluster,
                     "n_isolates": len(y), "n_resistant": int(y.sum()), "n_susceptible": int((y == 0).sum()),
                     "balanced_accuracy": balanced_accuracy_score(y, pred) if len(np.unique(y)) == 2 else np.nan,
                     "precision": precision_score(y, pred, zero_division=0) if pred.sum() else np.nan,
                     "recall": recall_score(y, pred, zero_division=0) if y.sum() else np.nan})
    return pd.DataFrame(rows)


def plot(metrics, curves, outdir, dpi):
    sns.set_theme(style="whitegrid", context="paper")
    matplotlib.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "axes.labelweight": "bold"})
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), layout="constrained")
    for design in DESIGNS:
        for input_key, (label, colour) in INPUTS.items():
            part = curves.loc[curves.design.eq(design) & curves.input_key.eq(input_key)]
            auc = metrics.loc[metrics.design.eq(design) & metrics.input_key.eq(input_key), "roc_auc"].mean()
            axes[0, 0].plot(part.fpr, part.mean_tpr, color=colour, ls="-" if design == "stratified" else "--",
                            label=f"{label}; {DESIGNS[design]} (AUC {auc:.3f})")
            axes[0, 0].fill_between(part.fpr, np.clip(part.mean_tpr - part.std_tpr, 0, 1),
                                   np.clip(part.mean_tpr + part.std_tpr, 0, 1), color=colour, alpha=.08)
    axes[0, 0].plot([0, 1], [0, 1], ":", color="grey")
    axes[0, 0].set(xlabel="False positive rate", ylabel="True positive rate", xlim=(0, 1), ylim=(0, 1))
    axes[0, 0].legend(fontsize=7, loc="lower right", frameon=False)
    order = [f"{design}_{key}" for design in DESIGNS for key in INPUTS]
    palette = {f"{design}_{key}": INPUTS[key][1] for design in DESIGNS for key in INPUTS}
    for ax, (metric, title) in zip(axes.flat[1:], METRICS.items()):
        sns.boxplot(data=metrics, x="analysis_id", y=metric, order=order, hue="analysis_id", hue_order=order,
                    palette=palette, dodge=False, legend=False, whis=1.5, showfliers=False, ax=ax)
        # A fixed jitter stream makes repeated exports reproducible.
        rng = np.random.default_rng(42)
        for i, name in enumerate(order):
            vals = metrics.loc[metrics.analysis_id.eq(name), metric]
            ax.scatter(i + rng.uniform(-.1, .1, len(vals)), vals, s=8, color="black", alpha=.4)
        ax.set(xlabel="", ylabel=title, ylim=(0, 1))
        ax.set_xticks(range(6), ["Variant", "Gene", "PC"] * 2)
        ax.axvline(2.5, color="grey", ls=":")
        if metric == "balanced_accuracy":
            ax.axhline(.5, color="grey", ls="--")
        ax.text(.25, 1.015, "Stratified CV", transform=ax.transAxes, ha="center", fontweight="bold")
        ax.text(.75, 1.015, "Clade-grouped CV", transform=ax.transAxes, ha="center", fontweight="bold")
    for label, ax in zip("BCDE", axes.flat):
        ax.set_title(label, loc="left", fontweight="bold", fontsize=18)
        ax.spines[["top", "right"]].set_visible(False)
    for suffix in ["png", "pdf"]:
        fig.savefig(outdir / f"random_forest_prediction.{suffix}", dpi=dpi)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--group-column", default="cluster")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    metrics, predictions = load_results(args.results_root)
    curves = roc_summary(predictions)
    args.outdir.mkdir(parents=True, exist_ok=False)
    metrics.to_csv(args.outdir / "fold_metrics.tsv", sep="\t", index=False)
    curves.to_csv(args.outdir / "mean_roc_coordinates.tsv", sep="\t", index=False)
    metrics.groupby(["design", "input_key"])[["roc_auc", *METRICS]].agg(["mean", "std"]).to_csv(args.outdir / "performance_summary.tsv", sep="\t")
    if args.metadata:
        scores = group_scores(predictions, pd.read_csv(args.metadata, sep="\t"), args.sample_column, args.group_column)
        scores.to_csv(args.outdir / "held_out_group_metrics_by_repeat.tsv", sep="\t", index=False)
        scores.groupby(["input_key", "cluster"])[list(METRICS)].agg(["mean", "std"]).to_csv(args.outdir / "held_out_group_summary.tsv", sep="\t")
    plot(metrics, curves, args.outdir, args.dpi)
    (args.outdir / "plot_info.json").write_text(json.dumps({"results_root": str(args.results_root.resolve()),
        "outer_fold_records": len(metrics), "roc_grid_points": 501, "roc_band": "Across-fold sample standard deviation; not a confidence interval"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
