#!/usr/bin/env python3
"""Run the reported score-label and SNP-screen shuffled-label controls."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from ml_common import (
    load_cache,
    load_labelled_metadata,
    repeated_outer_splits,
    screen_snps,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="rbin_flu")
    parser.add_argument("--group-col", default="cluster")
    parser.add_argument("--all-sample-predictions", required=True)
    parser.add_argument("--grouped-selected-snps", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--score-permutations", type=int, default=1000)
    parser.add_argument("--feature-permutations", type=int, default=100)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-k", type=int, default=10000)
    parser.add_argument("--max-missing", type=float, default=0.10)
    parser.add_argument("--min-frequency", type=float, default=0.05)
    parser.add_argument("--chunk-size", type=int, default=50000)
    return parser.parse_args()


def shuffled_labels(
    labels: np.ndarray,
    groups: np.ndarray,
    mode: str,
    generator: np.random.Generator,
) -> np.ndarray:
    shuffled = labels.copy()
    if mode == "global":
        generator.shuffle(shuffled)
        return shuffled
    for group in np.unique(groups):
        indices = np.flatnonzero(groups == group)
        shuffled[indices] = generator.permutation(shuffled[indices])
    return shuffled


def score_label_nulls(
    predictions: pd.DataFrame,
    sample_col: str,
    group_col: str,
    n_permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {
        "analysis_id",
        "model",
        sample_col,
        group_col,
        "true_label",
        "prediction_score_R",
    }
    missing = required.difference(predictions.columns)
    if missing:
        raise ValueError(f"Prediction table is missing columns: {sorted(missing)}")
    predictions = predictions.loc[
        predictions["analysis_id"].isin(["snp_grouped", "pc_grouped", "snp_pc_grouped"])
    ].copy()
    summary_rows: list[dict[str, object]] = []
    replicate_rows: list[dict[str, object]] = []
    for analysis_id, model in predictions[["analysis_id", "model"]].drop_duplicates().itertuples(index=False):
        part = predictions.loc[
            predictions["analysis_id"].eq(analysis_id) & predictions["model"].eq(model)
        ].copy()
        sample_table = part[[sample_col, group_col, "true_label"]].drop_duplicates(sample_col)
        sample_table["y"] = sample_table["true_label"].map({"S": 0, "R": 1}).astype(int)
        sample_index = pd.Series(np.arange(len(sample_table)), index=sample_table[sample_col].astype(str))
        prediction_sample_rows = sample_index.loc[part[sample_col].astype(str)].to_numpy(dtype=int)
        observed_labels = sample_table["y"].to_numpy(dtype=int)
        groups = sample_table[group_col].astype(str).to_numpy()
        scores = pd.to_numeric(part["prediction_score_R"], errors="raise").to_numpy(dtype=float)
        observed_row_labels = observed_labels[prediction_sample_rows]
        observed_auc = float(roc_auc_score(observed_row_labels, scores))
        observed_ap = float(average_precision_score(observed_row_labels, scores))
        for mode_number, mode in enumerate(["global", "within_cluster"]):
            generator = np.random.default_rng(seed + 100000 * mode_number)
            null_auc: list[float] = []
            null_ap: list[float] = []
            for permutation in range(1, n_permutations + 1):
                permuted = shuffled_labels(observed_labels, groups, mode, generator)
                row_labels = permuted[prediction_sample_rows]
                auc = float(roc_auc_score(row_labels, scores))
                ap = float(average_precision_score(row_labels, scores))
                null_auc.append(auc)
                null_ap.append(ap)
                replicate_rows.append(
                    {
                        "analysis_id": analysis_id,
                        "model": model,
                        "permutation_mode": mode,
                        "permutation": permutation,
                        "roc_auc": auc,
                        "average_precision": ap,
                    }
                )
            summary_rows.append(
                {
                    "analysis_id": analysis_id,
                    "model": model,
                    "permutation_mode": mode,
                    "n_permutations": n_permutations,
                    "observed_roc_auc": observed_auc,
                    "null_roc_auc_mean": float(np.mean(null_auc)),
                    "null_roc_auc_sd": float(np.std(null_auc, ddof=1)),
                    "roc_auc_empirical_p_greater_equal": (
                        1 + int(np.sum(np.asarray(null_auc) >= observed_auc))
                    )
                    / (n_permutations + 1),
                    "observed_average_precision": observed_ap,
                    "null_average_precision_mean": float(np.mean(null_ap)),
                    "null_average_precision_sd": float(np.std(null_ap, ddof=1)),
                    "average_precision_empirical_p_greater_equal": (
                        1 + int(np.sum(np.asarray(null_ap) >= observed_ap))
                    )
                    / (n_permutations + 1),
                }
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(replicate_rows)


def recurrence_statistics(counts: pd.Series, n_folds: int) -> dict[str, float | int]:
    return {
        "max_n_folds_selected": int(counts.max()),
        "n_features_selected_all_folds": int((counts == n_folds).sum()),
        "n_features_selected_at_least_80pct_folds": int((counts >= 0.8 * n_folds).sum()),
        "n_features_selected_at_least_50pct_folds": int((counts >= 0.5 * n_folds).sum()),
        "mean_n_folds_selected_per_feature": float(counts.mean()),
    }


def feature_screen_nulls(
    cache,
    metadata: pd.DataFrame,
    y: np.ndarray,
    cache_columns: np.ndarray,
    group_col: str,
    selected_path: str | Path,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = pd.read_csv(selected_path, sep="\t")
    selected = selected.loc[selected["analysis_id"].eq("snp_grouped")].copy()
    observed_candidates = np.sort(selected["variant_idx"].astype(int).unique())
    observed_counts = selected.groupby("variant_idx")["fold_number"].nunique().reindex(
        observed_candidates, fill_value=0
    )
    n_folds = args.n_splits * args.n_repeats
    observed = recurrence_statistics(observed_counts, n_folds)
    candidate_matrix = np.asarray(cache.snp_matrix[observed_candidates, :], dtype=np.int8)
    groups = metadata[group_col].astype(str).to_numpy()
    fixed_splits = list(
        repeated_outer_splits(
            y,
            groups,
            "group_stratified",
            args.n_splits,
            args.n_repeats,
            args.seed,
        )
    )
    replicate_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for mode_number, mode in enumerate(["global", "within_cluster"]):
        generator = np.random.default_rng(args.seed + 500000 + 100000 * mode_number)
        mode_statistics: list[dict[str, float | int]] = []
        for permutation in range(1, args.feature_permutations + 1):
            permuted_y = shuffled_labels(y, groups, mode, generator)
            fold_selected: list[np.ndarray] = []
            for fold_number, _repeat, _fold, train_rows, _test_rows in fixed_splits:
                ranked = screen_snps(
                    candidate_matrix,
                    cache_columns,
                    train_rows,
                    permuted_y,
                    top_k=args.top_k,
                    max_missing=args.max_missing,
                    min_frequency=args.min_frequency,
                    chunk_size=args.chunk_size,
                )
                selected_global = observed_candidates[ranked["variant_idx"].to_numpy(dtype=int)]
                fold_selected.append(
                    pd.DataFrame(
                        {"variant_idx": selected_global, "fold_number": fold_number}
                    )
                )
            recurrence = (
                pd.concat(fold_selected, ignore_index=True)
                .groupby("variant_idx")["fold_number"]
                .nunique()
                .reindex(observed_candidates, fill_value=0)
            )
            statistics = recurrence_statistics(recurrence, n_folds)
            mode_statistics.append(statistics)
            replicate_rows.append(
                {
                    "permutation_mode": mode,
                    "permutation": permutation,
                    **statistics,
                }
            )
        mode_frame = pd.DataFrame(mode_statistics)
        row: dict[str, object] = {
            "analysis_id": "snp_grouped",
            "permutation_mode": mode,
            "n_permutations": args.feature_permutations,
            "observed_n_folds": n_folds,
            "observed_n_candidate_features": len(observed_candidates),
            "observed_top_k": args.top_k,
        }
        for statistic, observed_value in observed.items():
            null_values = pd.to_numeric(mode_frame[statistic], errors="raise").to_numpy()
            row[f"observed_{statistic}"] = observed_value
            row[f"null_{statistic}_mean"] = float(null_values.mean())
            row[f"null_{statistic}_sd"] = float(null_values.std(ddof=1))
            row[f"{statistic}_empirical_p_greater_equal"] = (
                1 + int(np.sum(null_values >= observed_value))
            ) / (args.feature_permutations + 1)
        summary_rows.append(row)
    return pd.DataFrame(summary_rows), pd.DataFrame(replicate_rows)


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
    predictions = pd.read_csv(args.all_sample_predictions, sep="\t", low_memory=False)
    score_summary, score_replicates = score_label_nulls(
        predictions,
        args.sample_col,
        args.group_col,
        args.score_permutations,
        args.seed,
    )
    score_summary.to_csv(output / "score_label_null_summary.tsv", sep="\t", index=False)
    score_replicates.to_csv(output / "score_label_null_replicates.tsv", sep="\t", index=False)

    feature_summary, feature_replicates = feature_screen_nulls(
        cache,
        metadata,
        y,
        cache_columns,
        args.group_col,
        args.grouped_selected_snps,
        args,
    )
    feature_summary.to_csv(output / "feature_screen_null_summary.tsv", sep="\t", index=False)
    feature_replicates.to_csv(output / "feature_screen_null_replicates.tsv", sep="\t", index=False)
    write_json(
        output / "run_info.json",
        {
            "random_seed": args.seed,
            "score_label_permutations": args.score_permutations,
            "score_label_modes": ["global", "within_cluster"],
            "score_label_models_retrained": False,
            "feature_screen_permutations": args.feature_permutations,
            "feature_screen_modes": ["global", "within_cluster"],
            "feature_screen_models_retrained": False,
            "feature_screen_candidate_universe": "union of observed grouped-SNP top-10000 features",
            "feature_screen_fold_structure": "fixed observed 5-repeat, 5-fold grouped splits",
        },
    )


if __name__ == "__main__":
    main()
