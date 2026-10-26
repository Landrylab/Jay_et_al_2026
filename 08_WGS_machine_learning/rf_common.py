"""Shared data handling, model selection and metrics for genomic random forests."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import chi2
from sklearn.impute import SimpleImputer
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
                             confusion_matrix, f1_score, matthews_corrcoef, precision_score,
                             recall_score, roc_auc_score)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


RF_GRID = [
    {"max_features": "sqrt", "min_samples_leaf": 2, "max_depth": None},
    {"max_features": "sqrt", "min_samples_leaf": 1, "max_depth": None},
    {"max_features": "sqrt", "min_samples_leaf": 5, "max_depth": None},
    {"max_features": "sqrt", "min_samples_leaf": 10, "max_depth": None},
    {"max_features": "sqrt", "min_samples_leaf": 2, "max_depth": 20},
    {"max_features": 0.05, "min_samples_leaf": 2, "max_depth": None},
    {"max_features": 0.10, "min_samples_leaf": 2, "max_depth": None},
    {"max_features": 0.20, "min_samples_leaf": 2, "max_depth": None},
    {"max_features": 0.10, "min_samples_leaf": 5, "max_depth": 20},
]
PC_GRID = [2, 5, 10, 20, 40, 80]
PC_RF_GRID = [RF_GRID[i] for i in [0, 2, 7]]


def load_cache(cache_dir):
    cache_dir = Path(cache_dir)
    shape = json.loads((cache_dir / "shape.json").read_text())
    sample_path = next((cache_dir / name for name in ["phenotyped_plink_samples.txt", "sample_ids.txt", "phenotyped_vcf_samples.txt"]
                        if (cache_dir / name).exists()), None)
    if sample_path is None:
        raise ValueError("Cache lacks a sample-order file")
    samples = sample_path.read_text().splitlines()
    features = pd.read_csv(cache_dir / "variants.tsv", sep="\t")
    if len(samples) != shape["n_samples"] or len(features) != shape["n_variants"] or len(set(samples)) != len(samples):
        raise ValueError("Cache sample/feature dimensions or sample identifiers are inconsistent")
    if "variant_idx" not in features:
        features["variant_idx"] = np.arange(len(features))
    if not np.array_equal(features.variant_idx, np.arange(len(features))):
        raise ValueError("Cache feature indices must match matrix column order")
    matrix = cache_dir / "genotypes.int8.memmap"
    if matrix.stat().st_size != len(samples) * len(features):
        raise ValueError("Cached genotype matrix has the wrong file size")
    x = np.memmap(matrix, dtype=np.int8, mode="r", shape=(len(samples), len(features)))
    return x, features, samples


def aligned_metadata(metadata, samples, sample_column="sample", phenotype_column="rbin_flu", group_column="cluster"):
    frame = pd.read_csv(metadata, sep="\t", dtype={sample_column: str}) if isinstance(metadata, (str, Path)) else metadata.copy()
    if not {sample_column, phenotype_column}.issubset(frame) or frame[sample_column].duplicated().any():
        raise ValueError("Metadata must contain unique sample IDs and the phenotype column")
    frame = frame.set_index(sample_column).loc[samples]
    labels = frame[phenotype_column].astype(str).str.upper().map({"R": 1, "S": 0, "1": 1, "0": 0})
    if labels.isna().any() or labels.nunique() != 2:
        raise ValueError("Every cache sample must have a binary phenotype, with both classes represented")
    groups = None
    if group_column is not None:
        if group_column not in frame or frame[group_column].isna().any():
            raise ValueError("Complete group annotations are required for grouped validation")
        groups = frame[group_column].astype(str).to_numpy()
    return labels.to_numpy(dtype=int), groups


def forest(params, seed, n_estimators=500, n_jobs=-1):
    return RandomForestClassifier(n_estimators=n_estimators, class_weight="balanced_subsample",
                                  random_state=seed, n_jobs=n_jobs, **params)


def genotype_model(params, seed, n_estimators=500, n_jobs=-1):
    return Pipeline([("imputer", SimpleImputer(missing_values=-1, strategy="most_frequent")),
                     ("rf", forest(params, seed, n_estimators, n_jobs))])


def fit_pca_projector(x_train, n_components, seed):
    pipeline = Pipeline([
        ("imputer", SimpleImputer(missing_values=-1, strategy="mean", keep_empty_features=True)),
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_components, svd_solver="randomized", random_state=seed)),
    ])
    return pipeline.fit(np.asarray(x_train, dtype=np.float32))


def repeated_outer_splits(y, n_splits=5, n_repeats=5, seed=42, design="stratified", groups=None):
    if design == "stratified":
        splitter = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
        for number, (train, test) in enumerate(splitter.split(np.zeros(len(y)), y), 1):
            yield number, (number - 1) // n_splits + 1, (number - 1) % n_splits + 1, train, test
    elif design == "clade_grouped":
        if groups is None:
            raise ValueError("Grouped cross-validation requires sample groups")
        for repeat in range(1, n_repeats + 1):
            splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed + repeat - 1)
            for fold, (train, test) in enumerate(splitter.split(np.zeros(len(y)), y, groups), 1):
                yield (repeat - 1) * n_splits + fold, repeat, fold, train, test
    else:
        raise ValueError("Unknown cross-validation design: " + design)


def inner_splits(y, groups, seed, n_splits=3):
    if groups is None:
        splits = list(StratifiedKFold(n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))
    else:
        splits = list(StratifiedGroupKFold(n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups))
    return validate_inner_splits(y, groups, splits, n_splits)


def validate_inner_splits(y, groups, splits, n_splits=3):
    """Validate explicit inner partitions in the outer-training sample order."""
    if len(splits) != n_splits:
        raise ValueError("The number of supplied inner folds differs from n_splits")
    checked, validation_counts = [], np.zeros(len(y), dtype=int)
    for train, valid in splits:
        train, valid = np.asarray(train, dtype=int), np.asarray(valid, dtype=int)
        combined = np.r_[train, valid]
        if (train.ndim != 1 or valid.ndim != 1 or len(train) == 0 or len(valid) == 0
                or len(combined) != len(y)
                or not np.array_equal(np.sort(combined), np.arange(len(y)))):
            raise ValueError("Inner training and validation indices must partition the outer training samples")
        if len(np.unique(y[train])) != 2 or len(np.unique(y[valid])) != 2:
            raise ValueError("Each inner training and validation fold must contain both phenotype classes")
        if groups is not None and set(groups[train]) & set(groups[valid]):
            raise ValueError("Group overlap in inner cross-validation")
        validation_counts[valid] += 1
        checked.append((train, valid))
    if not np.all(validation_counts == 1):
        raise ValueError("Every outer-training sample must occur in exactly one inner validation fold")
    return checked


def select_features(x, y, train_idx, top_k=10000, max_missing=0.10, chunk_size=8192):
    """Rank variable features using only an outer training partition.

    The chi-squared score uses modal imputation within this training partition.
    The selected columns remain fixed during its inner model-parameter search;
    imputation is fitted separately inside every subsequent model fit.
    """
    if top_k < 1 or chunk_size < 1 or not 0 <= max_missing <= 1:
        raise ValueError("Feature count/chunk size must be positive and missingness must lie in [0,1]")
    rows = []
    for start in range(0, x.shape[1], chunk_size):
        stop = min(start + chunk_size, x.shape[1])
        chunk = np.asarray(x[train_idx, start:stop], dtype=np.float32)
        called = chunk >= 0
        missing_rate = 1 - called.mean(axis=0)
        dosage_sum = np.where(called, chunk, 0).sum(axis=0)
        called_n = called.sum(axis=0)
        frequency = np.divide(dosage_sum, 2 * called_n,
                              out=np.zeros_like(dosage_sum, dtype=np.float32), where=called_n > 0)
        maf = np.minimum(frequency, 1 - frequency)
        maximum = np.where(called, chunk, -np.inf).max(axis=0)
        minimum = np.where(called, chunk, np.inf).min(axis=0)
        keep = (missing_rate <= max_missing) & (maf >= 0) & (called_n > 0) & (maximum > minimum)
        if not keep.any():
            continue
        subset = chunk[:, keep]
        for col in subset.T:
            observed = col[col >= 0]
            mode = np.bincount(observed.astype(int), minlength=3).argmax()
            col[col < 0] = mode
        scores, probabilities = chi2(subset, y[train_idx])
        rows.append(pd.DataFrame({"variant_idx": np.arange(start, stop)[keep],
                                  "chi2": scores, "chi2_pvalue": probabilities,
                                  "train_missing_rate": missing_rate[keep], "train_maf": maf[keep]}))
    if not rows:
        raise ValueError("No variable features pass outer-training missingness filtering")
    return pd.concat(rows, ignore_index=True).sort_values(["chi2", "train_maf"], ascending=[False, False]).head(top_k).reset_index(drop=True)


def choose_threshold(y_true, probabilities):
    probabilities = np.asarray(probabilities)
    if not np.isfinite(probabilities).all() or not ((probabilities >= 0) & (probabilities <= 1)).all():
        raise ValueError("Threshold selection requires complete out-of-fold probabilities in [0,1]")
    candidates = np.unique(np.r_[0.0, 0.5, 1.0, probabilities])
    scores = np.array([balanced_accuracy_score(y_true, probabilities >= threshold) for threshold in candidates])
    best = np.flatnonzero(np.isclose(scores, scores.max()))
    distances = np.abs(candidates[best] - 0.5)
    closest = best[np.flatnonzero(np.isclose(distances, distances.min()))]
    selected = int(closest[np.argmin(candidates[closest])])
    return float(candidates[selected]), float(scores[selected])


def tune_genotypes(z_train, y, groups, seed, grid=RF_GRID, n_estimators=500, n_jobs=-1, n_splits=3,
                   explicit_splits=None):
    splits = (inner_splits(y, groups, seed, n_splits) if explicit_splits is None
              else validate_inner_splits(y, groups, explicit_splits, n_splits))
    tuning, probabilities = [], []
    for candidate, params in enumerate(grid, 1):
        oof = np.full(len(y), np.nan)
        scores = []
        for train, valid in splits:
            clf = genotype_model(params, seed, n_estimators, n_jobs)
            clf.fit(z_train[train], y[train])
            p = clf.predict_proba(z_train[valid])[:, 1]
            oof[valid] = p
            scores.append(float(roc_auc_score(y[valid], p)))
        tuning.append({"candidate_number": candidate, **params,
                       "mean_inner_roc_auc": float(np.mean(scores)), "split_roc_auc": json.dumps(scores)})
        probabilities.append(oof)
    best = min(tuning, key=lambda row: (-row["mean_inner_roc_auc"], row["candidate_number"]))
    index = best["candidate_number"] - 1
    threshold, score = choose_threshold(y, probabilities[index])
    return {"params": dict(grid[index]), "threshold": threshold, "inner_balanced_accuracy": score,
            "mean_inner_roc_auc": best["mean_inner_roc_auc"], "candidate_number": best["candidate_number"],
            "tuning": tuning}


def tune_pcs(x_train, y, groups, seed, pc_grid=PC_GRID, rf_grid=PC_RF_GRID,
             n_estimators=500, n_jobs=-1, n_splits=3, explicit_splits=None):
    splits = (inner_splits(y, groups, seed, n_splits) if explicit_splits is None
              else validate_inner_splits(y, groups, explicit_splits, n_splits))
    limit = min(x_train.shape[1], min(len(train) for train, _ in splits))
    if not pc_grid or any(pc < 1 or pc > limit for pc in pc_grid):
        raise ValueError("Every candidate PC count must fit within each inner training matrix")
    candidates = [(n, params) for n in pc_grid for params in rf_grid]
    probabilities = [np.full(len(y), np.nan) for _ in candidates]
    scores = [[] for _ in candidates]
    for split_number, (train, valid) in enumerate(splits, 1):
        projector = fit_pca_projector(x_train[train], max(pc_grid), seed + split_number)
        z_train = projector.transform(np.asarray(x_train[train], dtype=np.float32))
        z_valid = projector.transform(np.asarray(x_train[valid], dtype=np.float32))
        for index, (n_pcs, params) in enumerate(candidates):
            clf = forest(params, seed, n_estimators, n_jobs)
            clf.fit(z_train[:, :n_pcs], y[train])
            p = clf.predict_proba(z_valid[:, :n_pcs])[:, 1]
            probabilities[index][valid] = p
            scores[index].append(float(roc_auc_score(y[valid], p)))
    tuning = [{"candidate_number": i + 1, "n_pcs": n, **params,
               "mean_inner_roc_auc": float(np.mean(scores[i])), "split_roc_auc": json.dumps(scores[i])}
              for i, (n, params) in enumerate(candidates)]
    best = min(tuning, key=lambda row: (-row["mean_inner_roc_auc"], row["candidate_number"]))
    index = best["candidate_number"] - 1
    threshold, score = choose_threshold(y, probabilities[index])
    return {"params": dict(candidates[index][1]), "n_pcs": candidates[index][0],
            "threshold": threshold, "inner_balanced_accuracy": score,
            "mean_inner_roc_auc": best["mean_inner_roc_auc"], "candidate_number": best["candidate_number"],
            "tuning": tuning}


def metric_values(y, pred, prob):
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {"accuracy": float(accuracy_score(y, pred)),
            "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
            "precision": float(precision_score(y, pred, zero_division=0)),
            "recall": float(recall_score(y, pred, zero_division=0)),
            "f1": float(f1_score(y, pred, zero_division=0)),
            "mcc": float(matthews_corrcoef(y, pred)),
            "roc_auc": float(roc_auc_score(y, prob)) if len(np.unique(y)) == 2 else float("nan"),
            "average_precision": float(average_precision_score(y, prob)),
            "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}
