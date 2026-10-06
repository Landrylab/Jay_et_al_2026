#!/usr/bin/env python3
"""Tune variant and gene-level random forests within genetic groups and summarize feature importance.

Example:
  python 03_run_within_group_importance.py --cache-dir /path/to/gene_cache \
    --metadata /path/to/metadata.tsv --groups /path/to/groups.tsv \
    --outdir /path/to/within_group_results

The group table requires group_key, group_label, cluster and n_splits columns.
Feature screening uses outer-training samples, and inner cross-validation
selects forest hyperparameters and the probability threshold.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold

from rf_common import RF_GRID, aligned_metadata, genotype_model, load_cache, metric_values, select_features, tune_genotypes


def read_groups(path):
    groups = pd.read_csv(path, sep="\t", dtype={"cluster": str})
    required = {"group_key", "group_label", "cluster", "n_splits"}
    if groups.empty or not required <= set(groups.columns) or groups[list(required)].isna().any().any():
        raise ValueError("Group table requires nonmissing group_key, group_label, cluster and n_splits")
    if groups.group_key.duplicated().any() or groups.cluster.duplicated().any():
        raise ValueError("Group keys and cluster values must be unique")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", str(key)) for key in groups.group_key):
        raise ValueError("Group keys must be safe directory names containing letters, digits, '_' or '-'")
    splits = pd.to_numeric(groups.n_splits, errors="raise")
    if not np.isfinite(splits).all() or (splits < 2).any() or (splits % 1 != 0).any():
        raise ValueError("n_splits must be an integer of at least two")
    groups["n_splits"] = splits.astype(int)
    return groups


def outer_splits(y, n_splits, repeats, seed):
    counts = np.bincount(np.asarray(y, dtype=int), minlength=2)
    if counts.min() < n_splits:
        raise ValueError("Both classes must contain at least n_splits samples")
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=repeats, random_state=seed)
    for number, (train, test) in enumerate(cv.split(np.zeros(len(y)), y), 1):
        if np.bincount(y[train], minlength=2).min() < 3:
            raise ValueError("Each outer-training set needs at least three samples of each class for inner CV")
        yield number, (number - 1) // n_splits + 1, (number - 1) % n_splits + 1, train, test


def aggregate_importance(by_fold, features, n_outer_models):
    """Average normalized Gini importance over every outer-training model.

    A feature contributes zero for models in which it was not retained. Only
    features retained in at least one model are included in the summary.
    """
    if n_outer_models < 1 or by_fold.empty:
        raise ValueError("Importance aggregation requires fitted outer models")
    if by_fold.duplicated(["fold_number", "variant_idx"]).any():
        raise ValueError("Duplicate feature importance within an outer model")
    if by_fold.fold_number.nunique() != n_outer_models:
        raise ValueError("Importance rows do not cover the expected outer models")
    if not np.isfinite(by_fold.gini_importance).all() or (by_fold.gini_importance < -1e-12).any():
        raise ValueError("Gini importances must be finite and nonnegative within numerical precision")
    result = by_fold.groupby("variant_idx", as_index=False, sort=True).agg(
        importance_sum=("gini_importance", "sum"), n_folds_present=("fold_number", "nunique"))
    result["mean_importance_across_folds"] = result.importance_sum / n_outer_models
    result["rank"] = result.mean_importance_across_folds.rank(method="first", ascending=False).astype(int)
    if features.variant_idx.duplicated().any():
        raise ValueError("Feature metadata must have unique variant_idx values")
    result = result.merge(features, on="variant_idx", how="left", validate="one_to_one")
    result.insert(0, "method_key", "tuned")
    return result


def run_group(X, features, samples, y, group, outdir, repeats=5, seed=42,
              n_estimators=500, n_jobs=-1, top_k=10000, max_missing=0.10, grid=None, input_type="gene"):
    outdir = Path(outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Group output directory is not empty: " + str(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    n_splits = int(group["n_splits"])
    folds = list(outer_splits(y, n_splits, repeats, seed))
    metrics, predictions, importances, tuning_rows, choices, selections = [], [], [], [], [], []
    for number, repeat, fold, train, test in folds:
        selected_scores = select_features(X, y, train, top_k=top_k, max_missing=max_missing)
        selected = selected_scores.variant_idx.to_numpy(int)
        Z = np.asarray(X[:, selected], dtype=np.float32)
        model_seed = seed + number
        best = tune_genotypes(Z[train], y[train], None, model_seed,
                              grid=RF_GRID if grid is None else grid,
                              n_estimators=n_estimators, n_jobs=n_jobs)
        classifier = genotype_model(best["params"], model_seed, n_estimators=n_estimators, n_jobs=n_jobs)
        classifier.fit(Z[train], y[train])
        probability = classifier.predict_proba(Z[test])[:, 1]
        threshold = best["threshold"]
        predicted = (probability >= threshold).astype(int)
        context = {"group": group["group_label"], "group_key": group["group_key"], "cluster": str(group["cluster"]),
                   "input_key": input_type, "method_key": "tuned", "fold_number": number, "repeat": repeat, "fold": fold}
        metrics.append({**context, "threshold": threshold, "inner_balanced_accuracy": best["inner_balanced_accuracy"],
                        "n_train": len(train), "n_test": len(test), "n_features": len(selected),
                        **metric_values(y[test], predicted, probability)})
        for index, sample_index in enumerate(test):
            predictions.append({**context, "sample": samples[sample_index], "true_label": "R" if y[sample_index] else "S",
                                "probability_R": probability[index], "threshold": threshold,
                                "predicted_label": "R" if predicted[index] else "S"})
        importance = pd.DataFrame({"variant_idx": selected, "gini_importance": classifier.named_steps["rf"].feature_importances_})
        for key, value in context.items():
            importance[key] = value
        importances.append(importance.merge(features, on="variant_idx", how="left", validate="one_to_one"))
        selection = selected_scores.copy()
        selection["selection_rank"] = np.arange(1, len(selection) + 1)
        for key, value in context.items():
            selection[key] = value
        selections.append(selection)
        tuning = best["tuning"]
        if isinstance(tuning, pd.DataFrame):
            tuning = tuning.to_dict("records")
        tuning_rows.extend({**context, **row} for row in tuning)
        choices.append({**context, "candidate_number": best["candidate_number"], **best["params"],
                        "mean_inner_roc_auc": best["mean_inner_roc_auc"], "selected_threshold": threshold,
                        "inner_balanced_accuracy_at_threshold": best["inner_balanced_accuracy"]})
        print(f"{group['group_label']}: outer model {number}/{len(folds)}", flush=True)
    by_fold = pd.concat(importances, ignore_index=True)
    summary = aggregate_importance(by_fold, features, len(folds))
    summary["group_key"] = group["group_key"]
    summary["group"] = group["group_label"]
    summary["cluster"] = str(group["cluster"])
    outputs = {
        "cv_fold_metrics.tsv": pd.DataFrame(metrics), "cv_fold_predictions.tsv": pd.DataFrame(predictions),
        "feature_importance_by_fold.tsv": by_fold, "feature_importance_summary.tsv": summary,
        "cv_selected_features_by_fold.tsv": pd.concat(selections, ignore_index=True),
        "inner_cv_tuning_results.tsv": pd.DataFrame(tuning_rows), "selected_hyperparameters_by_fold.tsv": pd.DataFrame(choices),
    }
    for name, frame in outputs.items():
        frame.to_csv(outdir / name, sep="\t", index=False)
    run_info = {"group": group["group_label"], "group_key": group["group_key"], "pca_cluster": str(group["cluster"]),
                "input": "Gene-level variant based" if input_type == "gene" else "Variant based", "samples": len(y), "R": int(y.sum()), "S": int((y == 0).sum()),
                "outer_folds": n_splits, "outer_repeats": repeats, "seed": seed, "n_estimators": n_estimators,
                "feature_selection": "outer-training samples only; variable features; missingness threshold; chi-square ranking followed by training-frequency tie break",
                "top_k": top_k, "max_missing": max_missing, "class_weight": "balanced_subsample",
                "tuning": "inner three-fold stratified cross-validation maximizing mean ROC AUC",
                "threshold": "inner out-of-fold balanced accuracy; ties nearest 0.5 then lower threshold",
                "importance": "normalized Gini importance averaged over all outer models; zero for models not retaining the feature",
                "synthetic_samples": False}
    (outdir / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--input-type", choices=["variant", "gene"], default="gene")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--groups", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--phenotype-column", default="rbin_flu")
    parser.add_argument("--group-column", default="cluster")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--trees", type=int, default=500)
    parser.add_argument("--jobs", type=int, default=-1)
    parser.add_argument("--top-k", type=int, default=10000)
    parser.add_argument("--max-missing", type=float, default=0.10)
    args = parser.parse_args()
    if min(args.repeats, args.trees, args.top_k) < 1 or not 0 <= args.max_missing <= 1:
        parser.error("Repeats, trees and top-k must be positive; missingness must lie between zero and one")
    try:
        groups = read_groups(args.groups)
        X, features, samples = load_cache(args.cache_dir)
        if args.input_type == "gene" and not {"gene_id", "gene_symbol", "chrom", "pos"} <= set(features.columns):
            raise ValueError("Gene cache feature metadata requires gene_id, gene_symbol, chrom and pos")
        y, assignments = aligned_metadata(args.metadata, samples, sample_column=args.sample_column,
                                           phenotype_column=args.phenotype_column, group_column=args.group_column)
        assignments = np.asarray(assignments, dtype=str)
        all_metrics, all_predictions, all_importances = [], [], []
        for group in groups.to_dict("records"):
            selected = np.flatnonzero(assignments == str(group["cluster"]))
            if selected.size == 0:
                raise ValueError("No cache samples match cluster " + str(group["cluster"]))
            output = run_group(np.asarray(X[selected], dtype=np.int8), features, [samples[i] for i in selected], y[selected],
                               group, args.outdir / group["group_key"], args.repeats, args.seed, args.trees,
                               args.jobs, args.top_k, args.max_missing, input_type=args.input_type)
            all_metrics.append(output["cv_fold_metrics.tsv"])
            all_predictions.append(output["cv_fold_predictions.tsv"])
            all_importances.append(output["feature_importance_summary.tsv"])
        for name, frames in [("combined_cv_fold_metrics.tsv", all_metrics), ("combined_cv_fold_predictions.tsv", all_predictions),
                             ("combined_feature_importance_summary.tsv", all_importances)]:
            pd.concat(frames, ignore_index=True).to_csv(args.outdir / name, sep="\t", index=False)
    except (OSError, ValueError) as error:
        parser.exit(2, str(error) + "\n")


if __name__ == "__main__":
    main()
