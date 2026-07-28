#!/usr/bin/env python3
"""Run the reported same-cluster paired random-forest analyses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold

from ml_common import (
    impute_hom_alt,
    load_cache,
    load_labelled_metadata,
    metric_values,
    screen_snps,
    write_json,
)


STRATEGY_LABELS = {
    "pair_stratified": "Pair-stratified CV",
    "group_stratified": "Genetic-cluster-grouped CV",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="rbin_flu")
    parser.add_argument("--group-col", default="cluster")
    parser.add_argument("--tree-distance-matrix", required=True, help="Square patristic-distance NPY matrix.")
    parser.add_argument(
        "--tree-distance-samples",
        required=True,
        help="One-column TSV listing samples in the distance-matrix order.",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--strategies",
        nargs="+",
        choices=sorted(STRATEGY_LABELS),
        default=list(STRATEGY_LABELS),
    )
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-repeats", type=int, default=5)
    parser.add_argument("--inner-splits", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-k", type=int, default=10000)
    parser.add_argument("--max-missing", type=float, default=0.10)
    parser.add_argument("--min-frequency", type=float, default=0.05)
    parser.add_argument("--chunk-size", type=int, default=50000)
    parser.add_argument("--n-jobs", type=int, default=1)
    return parser.parse_args()


def make_pairs(
    metadata: pd.DataFrame,
    sample_col: str,
    label_col: str,
    group_col: str,
    distances: np.ndarray,
    distance_samples: list[str],
) -> pd.DataFrame:
    distance_index = pd.Series(np.arange(len(distance_samples)), index=np.asarray(distance_samples, dtype=str))
    missing = sorted(set(metadata[sample_col].astype(str)).difference(distance_index.index))
    if missing:
        raise ValueError(f"{len(missing)} labelled samples are absent from the distance matrix.")
    rows: list[dict[str, object]] = []
    pair_number = 0
    for group, part in metadata.groupby(group_col, sort=True):
        resistant = part.loc[part[label_col].astype(str).eq("R")]
        susceptible = part.loc[part[label_col].astype(str).eq("S")]
        if resistant.empty or susceptible.empty:
            continue
        r_dist = distance_index.loc[resistant[sample_col].astype(str)].to_numpy(dtype=int)
        s_dist = distance_index.loc[susceptible[sample_col].astype(str)].to_numpy(dtype=int)
        cost = distances[np.ix_(r_dist, s_dist)]
        r_assignment, s_assignment = linear_sum_assignment(cost)
        for r_local, s_local in zip(r_assignment, s_assignment):
            pair_number += 1
            r_row = resistant.iloc[int(r_local)]
            s_row = susceptible.iloc[int(s_local)]
            rows.append(
                {
                    "pair_id": f"pair_{pair_number:04d}",
                    "resistant_sample": str(r_row[sample_col]),
                    "susceptible_sample": str(s_row[sample_col]),
                    "genetic_cluster": str(group),
                    "tree_patristic_distance": float(cost[int(r_local), int(s_local)]),
                }
            )
    pairs = pd.DataFrame(rows).sort_values(
        ["genetic_cluster", "tree_patristic_distance", "resistant_sample", "susceptible_sample"]
    )
    pairs["pair_rank_by_match_quality"] = pairs["tree_patristic_distance"].rank(method="first").astype(int)
    return pairs.reset_index(drop=True)


def paired_rows(
    pairs: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    cache_columns: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    metadata_index = pd.Series(np.arange(len(metadata)), index=metadata[sample_col].astype(str))
    meta_cache = pd.Series(cache_columns, index=metadata[sample_col].astype(str))
    row_cache_columns: list[int] = []
    labels: list[int] = []
    pair_ids: list[str] = []
    clusters: list[str] = []
    metadata_rows: list[int] = []
    for row in pairs.itertuples(index=False):
        for sample, label in [(row.resistant_sample, 1), (row.susceptible_sample, 0)]:
            row_cache_columns.append(int(meta_cache.loc[sample]))
            labels.append(label)
            pair_ids.append(str(row.pair_id))
            clusters.append(str(row.genetic_cluster))
            metadata_rows.append(int(metadata_index.loc[sample]))
    return (
        np.asarray(row_cache_columns, dtype=int),
        np.asarray(labels, dtype=np.int8),
        np.asarray(pair_ids, dtype=str),
        np.asarray(clusters, dtype=str),
        np.asarray(metadata_rows, dtype=int),
    )


def signed_design(
    imputed_rows: np.ndarray,
    row_indices: np.ndarray,
    labels: np.ndarray,
    pair_ids: np.ndarray,
    cluster_lookup: dict[str, str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    deltas: list[np.ndarray] = []
    used_pairs: list[str] = []
    for pair_id in pd.unique(pair_ids[row_indices]):
        members = row_indices[pair_ids[row_indices] == pair_id]
        resistant = members[labels[members] == 1]
        susceptible = members[labels[members] == 0]
        if len(resistant) == 1 and len(susceptible) == 1:
            deltas.append(imputed_rows[resistant[0]] - imputed_rows[susceptible[0]])
            used_pairs.append(str(pair_id))
    delta = np.asarray(deltas, dtype=np.float32)
    features = np.vstack([delta, -delta])
    y = np.r_[np.ones(len(delta), dtype=np.int8), np.zeros(len(delta), dtype=np.int8)]
    groups_pair = np.r_[used_pairs, used_pairs]
    used_clusters = [cluster_lookup[pair_id] for pair_id in used_pairs]
    groups_cluster = np.r_[used_clusters, used_clusters]
    return features, y, groups_pair, groups_cluster


def score_held_out_pairs(
    estimator: RandomForestClassifier,
    imputed_rows: np.ndarray,
    test_rows: np.ndarray,
    labels: np.ndarray,
    pair_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    probabilities: list[float] = []
    records: list[dict[str, object]] = []
    for pair_id in pd.unique(pair_ids[test_rows]):
        members = test_rows[pair_ids[test_rows] == pair_id]
        resistant = members[labels[members] == 1]
        susceptible = members[labels[members] == 0]
        if len(resistant) != 1 or len(susceptible) != 1:
            continue
        delta = (imputed_rows[resistant[0]] - imputed_rows[susceptible[0]]).reshape(1, -1)
        probability = float(estimator.predict_proba(delta)[0, 1])
        probabilities.append(probability)
        records.append(
            {
                "pair_id": str(pair_id),
                "predicted_probability_resistant_direction": probability,
                "pair_score": probability - 0.5,
                "pair_correct": probability > 0.5,
            }
        )
    probabilities_array = np.asarray(probabilities, dtype=float)
    y_true = np.r_[np.ones(len(probabilities_array)), np.zeros(len(probabilities_array))].astype(np.int8)
    scores = np.r_[probabilities_array, 1.0 - probabilities_array]
    correct = probabilities_array > 0.5
    predictions = np.r_[correct, ~correct].astype(np.int8)
    return y_true, scores, predictions, pd.DataFrame(records)


def summarize_importance(by_fold: pd.DataFrame, variants: pd.DataFrame, n_folds: int) -> pd.DataFrame:
    summary = (
        by_fold.groupby(["cv_strategy", "variant_idx"], as_index=False)
        .agg(
            n_folds_selected=("fold_number", "nunique"),
            n_folds_nonzero=("importance", lambda values: int((values > 0).sum())),
            mean_selection_rank=("selection_rank", "mean"),
            mean_importance_when_selected=("importance", "mean"),
            maximum_importance=("importance", "max"),
            sum_importance=("importance", "sum"),
            mean_chi2=("chi2", "mean"),
            mean_train_minor_encoded_frequency=("train_minor_encoded_frequency", "mean"),
            mean_train_missing_rate=("train_missing_rate", "mean"),
        )
    )
    summary["mean_importance_across_all_folds"] = summary["sum_importance"] / n_folds
    summary = summary.merge(variants, on="variant_idx", how="left", validate="many_to_one")
    summary["importance_rank"] = summary.groupby("cv_strategy")[
        "mean_importance_across_all_folds"
    ].rank(method="first", ascending=False).astype(int)
    summary["top_one_percent"] = False
    grouped = summary["cv_strategy"].eq("group_stratified")
    if grouped.any():
        cutoff = summary.loc[grouped, "mean_importance_across_all_folds"].quantile(0.99)
        summary.loc[grouped, "top_one_percent"] = (
            summary.loc[grouped, "mean_importance_across_all_folds"] >= cutoff
        )
    return summary.sort_values(["cv_strategy", "importance_rank"])


def main() -> None:
    args = parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    cache = load_cache(args.cache)
    metadata, y, cache_columns = load_labelled_metadata(
        args.metadata,
        cache.sample_ids,
        args.sample_col,
        args.label_col,
        args.group_col,
    )
    distance_samples_frame = pd.read_csv(args.tree_distance_samples, sep="\t", header=None)
    distance_samples = distance_samples_frame.iloc[:, 0].astype(str).tolist()
    distances = np.load(args.tree_distance_matrix)
    if distances.shape != (len(distance_samples), len(distance_samples)):
        raise ValueError("Distance matrix dimensions do not match the distance-sample list.")
    pairs = make_pairs(
        metadata,
        args.sample_col,
        args.label_col,
        args.group_col,
        distances,
        distance_samples,
    )
    pairs.to_csv(output / "matched_pairs.tsv", sep="\t", index=False)
    row_cache_columns, row_labels, pair_ids, row_clusters, metadata_rows = paired_rows(
        pairs, metadata, args.sample_col, cache_columns
    )
    cluster_lookup = dict(zip(pairs["pair_id"].astype(str), pairs["genetic_cluster"].astype(str)))

    fold_records: list[dict[str, object]] = []
    prediction_parts: list[pd.DataFrame] = []
    importance_parts: list[pd.DataFrame] = []
    grid_parts: list[pd.DataFrame] = []
    parameter_grid = {
        "n_estimators": [300],
        "max_features": ["sqrt", 0.1],
        "min_samples_leaf": [1, 2],
        "max_depth": [None],
    }

    for strategy in args.strategies:
        split_groups = pair_ids if strategy == "pair_stratified" else row_clusters
        for repeat in range(1, args.n_repeats + 1):
            splitter = StratifiedGroupKFold(
                n_splits=args.n_splits,
                shuffle=True,
                random_state=args.seed + repeat - 1,
            )
            for fold, (train_rows, test_rows) in enumerate(
                splitter.split(np.zeros(len(row_labels)), row_labels, split_groups),
                start=1,
            ):
                fold_number = (repeat - 1) * args.n_splits + fold
                fold_seed = args.seed + fold_number + (10000 if strategy == "group_stratified" else 0)
                print(
                    f"{strategy}: repeat {repeat}/{args.n_repeats}, fold {fold}/{args.n_splits}",
                    flush=True,
                )
                selected = screen_snps(
                    cache.snp_matrix,
                    row_cache_columns,
                    train_rows,
                    row_labels,
                    top_k=args.top_k,
                    max_missing=args.max_missing,
                    min_frequency=args.min_frequency,
                    chunk_size=args.chunk_size,
                )
                selected_idx = selected["variant_idx"].to_numpy(dtype=int)
                _, imputed_all = impute_hom_alt(
                    cache.snp_matrix,
                    selected_idx,
                    row_cache_columns,
                    train_rows,
                    np.arange(len(row_labels)),
                )
                train_features, train_y, train_pair_groups, train_cluster_groups = signed_design(
                    imputed_all,
                    train_rows,
                    row_labels,
                    pair_ids,
                    cluster_lookup,
                )
                inner_groups = train_pair_groups if strategy == "pair_stratified" else train_cluster_groups
                estimator = RandomForestClassifier(
                    class_weight="balanced_subsample",
                    random_state=fold_seed,
                    n_jobs=1,
                )
                inner = StratifiedGroupKFold(
                    n_splits=args.inner_splits,
                    shuffle=True,
                    random_state=fold_seed,
                )
                grid = GridSearchCV(
                    estimator,
                    parameter_grid,
                    scoring="roc_auc",
                    cv=inner,
                    n_jobs=args.n_jobs,
                    refit=True,
                    return_train_score=False,
                )
                grid.fit(train_features, train_y, groups=inner_groups)
                fitted = grid.best_estimator_
                y_true, scores, predictions, pair_predictions = score_held_out_pairs(
                    fitted,
                    imputed_all,
                    test_rows,
                    row_labels,
                    pair_ids,
                )
                values = metric_values(y_true, scores, predictions)
                fold_records.append(
                    {
                        "analysis_id": "same_cluster_paired_random_forest",
                        "cv_strategy": strategy,
                        "cv_strategy_display": STRATEGY_LABELS[strategy],
                        "model": "random_forest",
                        "model_display": "Random forest",
                        "repeat": repeat,
                        "fold": fold,
                        "fold_number": fold_number,
                        "n_train_pairs": len(pd.unique(pair_ids[train_rows])),
                        "n_test_pairs": len(pd.unique(pair_ids[test_rows])),
                        "n_train_clusters": len(pd.unique(row_clusters[train_rows])),
                        "n_test_clusters": len(pd.unique(row_clusters[test_rows])),
                        "test_clusters": ",".join(sorted(pd.unique(row_clusters[test_rows]))),
                        "n_features": len(selected_idx),
                        "best_inner_roc_auc": float(grid.best_score_),
                        "best_parameters": json.dumps(grid.best_params_, sort_keys=True),
                        "pairwise_accuracy": float(pair_predictions["pair_correct"].mean()),
                        **values,
                    }
                )
                pair_predictions.insert(0, "cv_strategy", strategy)
                pair_predictions.insert(1, "repeat", repeat)
                pair_predictions.insert(2, "fold", fold)
                pair_predictions.insert(3, "fold_number", fold_number)
                prediction_parts.append(pair_predictions)

                fold_importance = selected.copy()
                fold_importance.insert(0, "cv_strategy", strategy)
                fold_importance.insert(1, "repeat", repeat)
                fold_importance.insert(2, "fold", fold)
                fold_importance.insert(3, "fold_number", fold_number)
                fold_importance["importance"] = fitted.feature_importances_
                importance_parts.append(fold_importance)

                grid_frame = pd.DataFrame(grid.cv_results_)
                grid_frame.insert(0, "cv_strategy", strategy)
                grid_frame.insert(1, "repeat", repeat)
                grid_frame.insert(2, "fold", fold)
                grid_frame.insert(3, "fold_number", fold_number)
                grid_parts.append(grid_frame)

    folds = pd.DataFrame(fold_records)
    folds.to_csv(output / "fold_metrics.tsv", sep="\t", index=False)
    summary_rows: list[dict[str, object]] = []
    for strategy, part in folds.groupby("cv_strategy", sort=False):
        for metric in [
            "accuracy",
            "balanced_accuracy",
            "precision",
            "recall",
            "f1",
            "mcc",
            "roc_auc",
            "average_precision",
            "pairwise_accuracy",
        ]:
            values = pd.to_numeric(part[metric], errors="coerce")
            summary_rows.append(
                {
                    "analysis_id": "same_cluster_paired_random_forest",
                    "cv_strategy": strategy,
                    "cv_strategy_display": STRATEGY_LABELS[strategy],
                    "model": "random_forest",
                    "metric": metric,
                    "mean": values.mean(),
                    "sd": values.std(ddof=1),
                    "min": values.min(),
                    "max": values.max(),
                    "n_folds": values.notna().sum(),
                }
            )
    pd.DataFrame(summary_rows).to_csv(output / "metric_summary.tsv", sep="\t", index=False)
    pd.concat(prediction_parts, ignore_index=True).to_csv(
        output / "pair_predictions.tsv", sep="\t", index=False
    )
    by_fold = pd.concat(importance_parts, ignore_index=True)
    by_fold.to_csv(output / "feature_importance_by_fold.tsv", sep="\t", index=False)
    summarize_importance(
        by_fold,
        cache.snp_variants,
        args.n_splits * args.n_repeats,
    ).to_csv(output / "feature_importance_summary.tsv", sep="\t", index=False)
    pd.concat(grid_parts, ignore_index=True).to_csv(
        output / "inner_grid_results.tsv", sep="\t", index=False
    )
    write_json(
        output / "run_info.json",
        {
            "n_labelled_isolates_considered": len(metadata),
            "n_resistant_considered": int(y.sum()),
            "n_susceptible_considered": int((y == 0).sum()),
            "n_strict_pairs": len(pairs),
            "median_tree_patristic_distance": float(pairs["tree_patristic_distance"].median()),
            "maximum_tree_patristic_distance": float(pairs["tree_patristic_distance"].max()),
            "matching": "Hungarian minimum-cost one-to-one R/S matching within K=23 genetic cluster",
            "strategies": args.strategies,
            "outer_cv": f"{args.n_repeats} repeats of {args.n_splits}-fold StratifiedGroupKFold",
            "inner_cv_folds": args.inner_splits,
            "random_seed": args.seed,
            "snp_top_k": args.top_k,
            "snp_max_missing": args.max_missing,
            "snp_min_encoded_state_frequency": args.min_frequency,
            "random_forest_grid": parameter_grid,
            "feature_importance_for_marker_ranking": "group_stratified only",
            "importance_averaging": "unselected SNPs assigned zero in each of 25 outer folds",
        },
    )


if __name__ == "__main__":
    main()
