#!/usr/bin/env python3
"""Shared functions for the reported fluconazole-resistance ML analyses."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import (
    ParameterGrid,
    StratifiedGroupKFold,
    StratifiedKFold,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

try:
    from xgboost import XGBClassifier
except ImportError:  # pragma: no cover
    XGBClassifier = None


MODEL_ORDER = ["logistic_regression", "linear_svm", "random_forest", "xgboost"]
MODEL_LABELS = {
    "logistic_regression": "Logistic regression",
    "linear_svm": "Linear SVM",
    "random_forest": "Random forest",
    "xgboost": "XGBoost",
}
METRICS = [
    "accuracy",
    "balanced_accuracy",
    "precision",
    "recall",
    "f1",
    "mcc",
    "roc_auc",
    "average_precision",
]


@dataclass
class GenotypeCache:
    directory: Path
    sample_ids: np.ndarray
    snp_matrix: np.memmap
    snp_variants: pd.DataFrame
    pc_matrix: np.memmap
    pc_variants: pd.DataFrame
    manifest: dict[str, object]


def load_cache(cache_dir: str | Path) -> GenotypeCache:
    directory = Path(cache_dir)
    manifest = json.loads((directory / "cache_manifest.json").read_text())
    sample_ids = pd.read_csv(directory / "sample_ids.tsv", sep="\t")["sample"].astype(str).to_numpy()
    n_samples = int(manifest["n_samples"])
    n_snp = int(manifest["n_snp_records"])
    n_pc = int(manifest["n_pc_markers"])
    snp_matrix = np.memmap(
        directory / "snp_genotypes.int8.mmap",
        dtype=np.int8,
        mode="r",
        shape=(n_snp, n_samples),
    )
    pc_matrix = np.memmap(
        directory / "pc_genotypes.int8.mmap",
        dtype=np.int8,
        mode="r",
        shape=(n_pc, n_samples),
    )
    return GenotypeCache(
        directory=directory,
        sample_ids=sample_ids,
        snp_matrix=snp_matrix,
        snp_variants=pd.read_csv(directory / "snp_variants.tsv", sep="\t"),
        pc_matrix=pc_matrix,
        pc_variants=pd.read_csv(directory / "pc_variants.tsv", sep="\t"),
        manifest=manifest,
    )


def load_labelled_metadata(
    metadata_path: str | Path,
    sample_ids: np.ndarray,
    sample_col: str,
    label_col: str,
    group_col: str,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    meta = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    required = {sample_col, label_col, group_col}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")
    meta[sample_col] = meta[sample_col].astype(str)
    meta = meta.drop_duplicates(sample_col).set_index(sample_col).reindex(sample_ids)
    keep = meta[label_col].astype(str).isin(["R", "S"]) & meta[group_col].notna()
    meta = meta.loc[keep].reset_index().rename(columns={"index": sample_col})
    y = meta[label_col].astype(str).map({"S": 0, "R": 1}).to_numpy(dtype=np.int8)
    cache_index = pd.Series(np.arange(len(sample_ids)), index=sample_ids)
    columns = cache_index.loc[meta[sample_col].astype(str)].to_numpy(dtype=int)
    return meta, y, columns


def repeated_outer_splits(
    y: np.ndarray,
    groups: np.ndarray,
    strategy: str,
    n_splits: int,
    n_repeats: int,
    seed: int,
) -> Iterator[tuple[int, int, int, np.ndarray, np.ndarray]]:
    for repeat in range(1, n_repeats + 1):
        split_seed = seed + repeat - 1
        if strategy == "stratified":
            splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=split_seed)
            split_iter = splitter.split(np.zeros(len(y)), y)
        elif strategy == "group_stratified":
            splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=split_seed)
            split_iter = splitter.split(np.zeros(len(y)), y, groups)
        else:
            raise ValueError(f"Unsupported CV strategy: {strategy}")
        for fold, (train_idx, test_idx) in enumerate(split_iter, start=1):
            fold_number = (repeat - 1) * n_splits + fold
            yield fold_number, repeat, fold, train_idx, test_idx


def inner_splits(
    y: np.ndarray,
    groups: np.ndarray,
    strategy: str,
    n_splits: int,
    seed: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    if strategy == "stratified":
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        return list(splitter.split(np.zeros(len(y)), y))
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(splitter.split(np.zeros(len(y)), y, groups))


def screen_snps(
    matrix: np.ndarray,
    cache_columns: np.ndarray,
    train_rows: np.ndarray,
    y: np.ndarray,
    top_k: int = 10000,
    max_missing: float = 0.10,
    min_frequency: float = 0.05,
    chunk_size: int = 50000,
) -> pd.DataFrame:
    """Rank SNPs by a training-only 2x2 Pearson chi-square statistic."""
    train_columns = cache_columns[train_rows]
    y_train = y[train_rows]
    r_mask = y_train == 1
    s_mask = ~r_mask
    rows: list[pd.DataFrame] = []
    for start in range(0, matrix.shape[0], chunk_size):
        stop = min(start + chunk_size, matrix.shape[0])
        block = np.asarray(matrix[start:stop, train_columns], dtype=np.int8)
        observed = block >= 0
        n_observed = observed.sum(axis=1)
        missing_rate = 1.0 - n_observed / block.shape[1]
        count_one = ((block == 1) & observed).sum(axis=1)
        frequency = np.divide(
            count_one,
            n_observed,
            out=np.zeros_like(count_one, dtype=float),
            where=n_observed > 0,
        )
        minor_frequency = np.minimum(frequency, 1.0 - frequency)

        r = block[:, r_mask]
        s = block[:, s_mask]
        a = (r == 1).sum(axis=1).astype(float)
        b = (r == 0).sum(axis=1).astype(float)
        c = (s == 1).sum(axis=1).astype(float)
        d = (s == 0).sum(axis=1).astype(float)
        total = a + b + c + d
        denominator = (a + b) * (c + d) * (a + c) * (b + d)
        chi2 = np.divide(
            total * np.square(a * d - b * c),
            denominator,
            out=np.zeros_like(total),
            where=denominator > 0,
        )
        valid = (
            (missing_rate <= max_missing)
            & (minor_frequency >= min_frequency)
            & (minor_frequency > 0)
            & np.isfinite(chi2)
        )
        if valid.any():
            local = np.flatnonzero(valid)
            rows.append(
                pd.DataFrame(
                    {
                        "variant_idx": start + local,
                        "chi2": chi2[local],
                        "train_minor_encoded_frequency": minor_frequency[local],
                        "train_missing_rate": missing_rate[local],
                    }
                )
            )
    if not rows:
        raise ValueError("No SNPs passed the training-fold filters.")
    ranked = pd.concat(rows, ignore_index=True).sort_values(
        ["chi2", "train_minor_encoded_frequency", "variant_idx"],
        ascending=[False, False, True],
    )
    ranked = ranked.head(min(top_k, len(ranked))).reset_index(drop=True)
    ranked.insert(0, "selection_rank", np.arange(1, len(ranked) + 1))
    return ranked


def impute_hom_alt(
    matrix: np.ndarray,
    variant_indices: np.ndarray,
    cache_columns: np.ndarray,
    train_rows: np.ndarray,
    target_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    train = np.asarray(matrix[np.ix_(variant_indices, cache_columns[train_rows])], dtype=np.float32).T
    target = np.asarray(matrix[np.ix_(variant_indices, cache_columns[target_rows])], dtype=np.float32).T
    observed = train >= 0
    means = np.divide(
        np.where(observed, train, 0).sum(axis=0),
        observed.sum(axis=0),
        out=np.zeros(train.shape[1], dtype=float),
        where=observed.sum(axis=0) > 0,
    )
    modes = (means >= 0.5).astype(np.float32)
    train = np.where(train < 0, modes, train)
    target = np.where(target < 0, modes, target)
    return train.astype(np.float32, copy=False), target.astype(np.float32, copy=False)


def model_parameter_grid(model: str) -> list[dict[str, object]]:
    if model == "logistic_regression":
        return list(ParameterGrid({"C": [0.01, 0.1, 1.0, 10.0], "penalty": ["l1", "l2"]}))
    if model == "linear_svm":
        return list(ParameterGrid({"C": [0.01, 0.1, 1.0, 10.0]}))
    if model == "random_forest":
        return list(
            ParameterGrid(
                {
                    "n_estimators": [300],
                    "max_features": ["sqrt", 0.1],
                    "min_samples_leaf": [1, 2],
                    "max_depth": [None],
                }
            )
        )
    if model == "xgboost":
        return list(
            ParameterGrid(
                {
                    "n_estimators": [100, 200],
                    "max_depth": [2, 3],
                    "learning_rate": [0.05, 0.1],
                    "colsample_bytree": [0.5],
                    "subsample": [0.8],
                }
            )
        )
    raise ValueError(f"Unknown model: {model}")


def make_estimator(model: str, params: dict[str, object], seed: int):
    if model == "logistic_regression":
        classifier = LogisticRegression(
            C=float(params["C"]),
            penalty=str(params["penalty"]),
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        )
        return make_pipeline(StandardScaler(), classifier)
    if model == "linear_svm":
        classifier = LinearSVC(
            C=float(params["C"]),
            class_weight="balanced",
            dual="auto",
            max_iter=10000,
            random_state=seed,
        )
        return make_pipeline(StandardScaler(), classifier)
    if model == "random_forest":
        return RandomForestClassifier(
            **params,
            class_weight="balanced_subsample",
            n_jobs=1,
            random_state=seed,
        )
    if model == "xgboost":
        if XGBClassifier is None:
            raise ImportError("xgboost is required for the XGBoost analysis.")
        return XGBClassifier(
            **params,
            objective="binary:logistic",
            eval_metric="logloss",
            n_jobs=1,
            random_state=seed,
        )
    raise ValueError(f"Unknown model: {model}")


def prediction_score(estimator, features: np.ndarray) -> np.ndarray:
    if hasattr(estimator, "predict_proba"):
        return np.asarray(estimator.predict_proba(features)[:, 1], dtype=float)
    return np.asarray(estimator.decision_function(features), dtype=float)


def default_predictions(model: str, scores: np.ndarray) -> np.ndarray:
    threshold = 0.0 if model == "linear_svm" else 0.5
    return (scores >= threshold).astype(np.int8)


def metric_values(y_true: np.ndarray, scores: np.ndarray, predictions: np.ndarray) -> dict[str, float | int]:
    matrix = confusion_matrix(y_true, predictions, labels=[0, 1])
    result: dict[str, float | int] = {
        "accuracy": float(accuracy_score(y_true, predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "precision": float(precision_score(y_true, predictions, zero_division=0)),
        "recall": float(recall_score(y_true, predictions, zero_division=0)),
        "f1": float(f1_score(y_true, predictions, zero_division=0)),
        "mcc": float(matthews_corrcoef(y_true, predictions)),
        "roc_auc": float(roc_auc_score(y_true, scores)) if len(np.unique(y_true)) == 2 else math.nan,
        "average_precision": float(average_precision_score(y_true, scores)),
        "tn": int(matrix[0, 0]),
        "fp": int(matrix[0, 1]),
        "fn": int(matrix[1, 0]),
        "tp": int(matrix[1, 1]),
    }
    return result


def summarize_fold_metrics(folds: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["analysis_id", "analysis_label", "feature_set", "cv_strategy", "model", "model_display"]
    rows: list[dict[str, object]] = []
    for keys, part in folds.groupby(group_cols, sort=False):
        base = dict(zip(group_cols, keys))
        for metric in METRICS:
            values = pd.to_numeric(part[metric], errors="coerce")
            rows.append(
                {
                    **base,
                    "metric": metric,
                    "mean": float(values.mean()),
                    "sd": float(values.std(ddof=1)),
                    "min": float(values.min()),
                    "max": float(values.max()),
                    "n_folds": int(values.notna().sum()),
                }
            )
    return pd.DataFrame(rows)


def write_json(path: str | Path, payload: dict[str, object]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
