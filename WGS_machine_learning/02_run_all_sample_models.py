#!/usr/bin/env python3
"""Run the six reported all-sample SNP, PC and SNP+PC validation analyses."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from ml_common import (
    MODEL_LABELS,
    MODEL_ORDER,
    default_predictions,
    impute_hom_alt,
    inner_splits,
    load_cache,
    load_labelled_metadata,
    make_estimator,
    metric_values,
    model_parameter_grid,
    prediction_score,
    repeated_outer_splits,
    screen_snps,
    summarize_fold_metrics,
    write_json,
)


ANALYSES = {
    "snp_stratified": ("snp", "stratified", "SNP only, stratified CV"),
    "pc_stratified": ("pc", "stratified", "PC only, stratified CV"),
    "snp_pc_stratified": ("snp_pc", "stratified", "SNP+PC, stratified CV"),
    "snp_grouped": ("snp", "group_stratified", "SNP only, group-stratified CV"),
    "pc_grouped": ("pc", "group_stratified", "PC only, group-stratified CV"),
    "snp_pc_grouped": ("snp_pc", "group_stratified", "SNP+PC, group-stratified CV"),
}
PC_GRIDS = {
    "snp": [0],
    "pc": [2, 5, 10, 20, 40, 80],
    "snp_pc": [20, 40, 80],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="rbin_flu")
    parser.add_argument("--group-col", default="cluster")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--analyses",
        nargs="+",
        choices=sorted(ANALYSES),
        default=list(ANALYSES),
        help="Analyses to run; all six reported analyses are used by default.",
    )
    parser.add_argument("--models", nargs="+", choices=MODEL_ORDER, default=MODEL_ORDER)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-repeats", type=int, default=5)
    parser.add_argument("--inner-splits", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-k", type=int, default=10000)
    parser.add_argument("--max-missing", type=float, default=0.10)
    parser.add_argument("--min-frequency", type=float, default=0.05)
    parser.add_argument("--chunk-size", type=int, default=50000)
    return parser.parse_args()


def pc_projection(
    matrix: np.ndarray,
    cache_columns: np.ndarray,
    train_rows: np.ndarray,
    target_rows: np.ndarray,
    n_components: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    train = np.asarray(matrix[:, cache_columns[train_rows]], dtype=np.float32).T
    target = np.asarray(matrix[:, cache_columns[target_rows]], dtype=np.float32).T
    train[train < 0] = np.nan
    target[target < 0] = np.nan
    imputer = SimpleImputer(strategy="mean", keep_empty_features=True)
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(imputer.fit_transform(train))
    target_scaled = scaler.transform(imputer.transform(target))
    n_components = min(n_components, train_scaled.shape[0] - 1, train_scaled.shape[1])
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=seed)
    train_pcs = pca.fit_transform(train_scaled)
    target_pcs = pca.transform(target_scaled)
    return (
        train_pcs.astype(np.float32, copy=False),
        target_pcs.astype(np.float32, copy=False),
        float(pca.explained_variance_ratio_.sum()),
    )


def combine_features(
    feature_set: str,
    snp_train: np.ndarray | None,
    snp_target: np.ndarray | None,
    pc_train: np.ndarray | None,
    pc_target: np.ndarray | None,
    n_pcs: int,
) -> tuple[np.ndarray, np.ndarray]:
    if feature_set == "snp":
        return snp_train, snp_target
    if feature_set == "pc":
        return pc_train[:, :n_pcs], pc_target[:, :n_pcs]
    return (
        np.hstack([snp_train, pc_train[:, :n_pcs]]),
        np.hstack([snp_target, pc_target[:, :n_pcs]]),
    )


def tune_model(
    model: str,
    feature_set: str,
    pc_grid: list[int],
    inner_data: list[dict[str, object]],
    seed: int,
) -> tuple[dict[str, object], int, float, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for params, n_pcs in itertools.product(model_parameter_grid(model), pc_grid):
        scores: list[float] = []
        for split_number, split_data in enumerate(inner_data, start=1):
            train_features, valid_features = combine_features(
                feature_set,
                split_data["snp_train"],
                split_data["snp_valid"],
                split_data["pc_train"],
                split_data["pc_valid"],
                n_pcs,
            )
            estimator = make_estimator(model, params, seed + split_number)
            estimator.fit(train_features, split_data["y_train"])
            valid_scores = prediction_score(estimator, valid_features)
            scores.append(float(roc_auc_score(split_data["y_valid"], valid_scores)))
        rows.append(
            {
                "model": model,
                "n_pcs": n_pcs,
                "mean_inner_roc_auc": float(np.mean(scores)),
                "sd_inner_roc_auc": float(np.std(scores, ddof=1)),
                "parameters": json.dumps(params, sort_keys=True),
            }
        )
    grid = pd.DataFrame(rows).sort_values(
        ["mean_inner_roc_auc", "n_pcs", "parameters"],
        ascending=[False, True, True],
    )
    best = grid.iloc[0]
    return json.loads(best["parameters"]), int(best["n_pcs"]), float(best["mean_inner_roc_auc"]), grid


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
    groups = metadata[args.group_col].astype(str).to_numpy()

    fold_rows: list[dict[str, object]] = []
    prediction_parts: list[pd.DataFrame] = []
    selected_parts: list[pd.DataFrame] = []
    grid_parts: list[pd.DataFrame] = []

    for analysis_id in args.analyses:
        feature_set, cv_strategy, analysis_label = ANALYSES[analysis_id]
        for fold_number, repeat, fold, train_rows, test_rows in repeated_outer_splits(
            y, groups, cv_strategy, args.n_splits, args.n_repeats, args.seed
        ):
            fold_seed = args.seed + fold_number
            print(
                f"{analysis_id}: repeat {repeat}/{args.n_repeats}, fold {fold}/{args.n_splits}",
                flush=True,
            )
            selected = None
            selected_idx = np.asarray([], dtype=int)
            if feature_set in {"snp", "snp_pc"}:
                selected = screen_snps(
                    cache.snp_matrix,
                    cache_columns,
                    train_rows,
                    y,
                    top_k=args.top_k,
                    max_missing=args.max_missing,
                    min_frequency=args.min_frequency,
                    chunk_size=args.chunk_size,
                )
                selected_idx = selected["variant_idx"].to_numpy(dtype=int)
                selected_output = selected.merge(
                    cache.snp_variants,
                    on="variant_idx",
                    how="left",
                    validate="one_to_one",
                )
                selected_output.insert(0, "analysis_id", analysis_id)
                selected_output.insert(1, "repeat", repeat)
                selected_output.insert(2, "fold", fold)
                selected_output.insert(3, "fold_number", fold_number)
                selected_parts.append(selected_output)

            local_inner = inner_splits(
                y[train_rows],
                groups[train_rows],
                cv_strategy,
                args.inner_splits,
                fold_seed,
            )
            inner_data: list[dict[str, object]] = []
            max_pcs = max(PC_GRIDS[feature_set])
            for inner_train_local, inner_valid_local in local_inner:
                inner_train = train_rows[inner_train_local]
                inner_valid = train_rows[inner_valid_local]
                snp_train = snp_valid = None
                pc_train = pc_valid = None
                if feature_set in {"snp", "snp_pc"}:
                    snp_train, snp_valid = impute_hom_alt(
                        cache.snp_matrix,
                        selected_idx,
                        cache_columns,
                        inner_train,
                        inner_valid,
                    )
                if feature_set in {"pc", "snp_pc"}:
                    pc_train, pc_valid, _ = pc_projection(
                        cache.pc_matrix,
                        cache_columns,
                        inner_train,
                        inner_valid,
                        max_pcs,
                        fold_seed,
                    )
                inner_data.append(
                    {
                        "snp_train": snp_train,
                        "snp_valid": snp_valid,
                        "pc_train": pc_train,
                        "pc_valid": pc_valid,
                        "y_train": y[inner_train],
                        "y_valid": y[inner_valid],
                    }
                )

            final_snp_train = final_snp_test = None
            final_pc_train = final_pc_test = None
            explained_variance = np.nan
            if feature_set in {"snp", "snp_pc"}:
                final_snp_train, final_snp_test = impute_hom_alt(
                    cache.snp_matrix,
                    selected_idx,
                    cache_columns,
                    train_rows,
                    test_rows,
                )
            if feature_set in {"pc", "snp_pc"}:
                final_pc_train, final_pc_test, explained_variance = pc_projection(
                    cache.pc_matrix,
                    cache_columns,
                    train_rows,
                    test_rows,
                    max_pcs,
                    fold_seed,
                )

            for model in args.models:
                best_params, best_n_pcs, best_inner_auc, grid = tune_model(
                    model,
                    feature_set,
                    PC_GRIDS[feature_set],
                    inner_data,
                    fold_seed,
                )
                grid.insert(0, "analysis_id", analysis_id)
                grid.insert(1, "repeat", repeat)
                grid.insert(2, "fold", fold)
                grid.insert(3, "fold_number", fold_number)
                grid_parts.append(grid)

                train_features, test_features = combine_features(
                    feature_set,
                    final_snp_train,
                    final_snp_test,
                    final_pc_train,
                    final_pc_test,
                    best_n_pcs,
                )
                estimator = make_estimator(model, best_params, fold_seed)
                estimator.fit(train_features, y[train_rows])
                scores = prediction_score(estimator, test_features)
                predictions = default_predictions(model, scores)
                values = metric_values(y[test_rows], scores, predictions)
                fold_rows.append(
                    {
                        "analysis_id": analysis_id,
                        "analysis_label": analysis_label,
                        "feature_set": feature_set,
                        "cv_strategy": cv_strategy,
                        "model": model,
                        "model_display": MODEL_LABELS[model],
                        "repeat": repeat,
                        "fold": fold,
                        "fold_number": fold_number,
                        "n_train": len(train_rows),
                        "n_test": len(test_rows),
                        "train_resistant": int(y[train_rows].sum()),
                        "train_susceptible": int((y[train_rows] == 0).sum()),
                        "test_resistant": int(y[test_rows].sum()),
                        "test_susceptible": int((y[test_rows] == 0).sum()),
                        "n_train_groups": len(np.unique(groups[train_rows])),
                        "n_test_groups": len(np.unique(groups[test_rows])),
                        "test_groups": ",".join(sorted(np.unique(groups[test_rows]))),
                        "n_snp_features": len(selected_idx),
                        "n_pc_features": best_n_pcs,
                        "pca_explained_variance_sum": explained_variance,
                        "best_inner_roc_auc": best_inner_auc,
                        "best_parameters": json.dumps(best_params, sort_keys=True),
                        **values,
                    }
                )
                part = metadata.iloc[test_rows][[args.sample_col, args.group_col]].copy()
                part.insert(0, "analysis_id", analysis_id)
                part.insert(1, "feature_set", feature_set)
                part.insert(2, "cv_strategy", cv_strategy)
                part.insert(3, "model", model)
                part.insert(4, "model_display", MODEL_LABELS[model])
                part["repeat"] = repeat
                part["fold"] = fold
                part["fold_number"] = fold_number
                part["true_label"] = np.where(y[test_rows] == 1, "R", "S")
                part["predicted_label"] = np.where(predictions == 1, "R", "S")
                part["prediction_score_R"] = scores
                prediction_parts.append(part)

    folds = pd.DataFrame(fold_rows)
    folds.to_csv(output / "fold_metrics.tsv", sep="\t", index=False)
    summarize_fold_metrics(folds).to_csv(output / "metric_summary.tsv", sep="\t", index=False)
    pd.concat(prediction_parts, ignore_index=True).to_csv(
        output / "fold_predictions.tsv", sep="\t", index=False
    )
    if selected_parts:
        pd.concat(selected_parts, ignore_index=True).to_csv(
            output / "selected_snps_by_fold.tsv", sep="\t", index=False
        )
    pd.concat(grid_parts, ignore_index=True).to_csv(
        output / "inner_grid_results.tsv", sep="\t", index=False
    )
    write_json(
        output / "run_info.json",
        {
            "analyses": args.analyses,
            "models": args.models,
            "n_samples": len(y),
            "n_resistant": int(y.sum()),
            "n_susceptible": int((y == 0).sum()),
            "n_groups": int(len(np.unique(groups))),
            "outer_cv": f"{args.n_repeats} repeats of {args.n_splits}-fold validation",
            "inner_cv_folds": args.inner_splits,
            "random_seed": args.seed,
            "snp_top_k": args.top_k,
            "snp_max_missing": args.max_missing,
            "snp_min_encoded_state_frequency": args.min_frequency,
            "pc_marker_count": int(cache.manifest["n_pc_markers"]),
            "pc_only_component_grid": PC_GRIDS["pc"],
            "snp_pc_component_grid": PC_GRIDS["snp_pc"],
            "positive_class": "R",
        },
    )


if __name__ == "__main__":
    main()
