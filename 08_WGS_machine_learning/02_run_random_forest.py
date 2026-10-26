#!/usr/bin/env python3
"""Evaluate variant-, gene- and PC-based random forests with repeated held-out folds.

Model parameters and decision thresholds are selected using inner validation.
Variant/gene feature selection uses only each outer training partition and is
held fixed during inner tuning. PC preprocessing is fitted inside every inner
training split and refitted on the outer training partition for final prediction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn

from rf_common import (PC_GRID, PC_RF_GRID, RF_GRID, aligned_metadata, fit_pca_projector,
                       forest, genotype_model, load_cache, metric_values, repeated_outer_splits,
                       select_features, tune_genotypes, tune_pcs)


DESIGN_LABELS = {"stratified": "Repeated stratified CV", "clade_grouped": "Clade-grouped CV"}
INPUT_LABELS = {"variant": "Variant based", "gene": "Gene-level variant based", "pc": "Phylogeny based"}


def run_analysis(x, features, samples, y, groups, input_type, design, destination,
                 n_splits=5, n_repeats=5, seed=42, inner_folds=3, n_estimators=500,
                 n_jobs=-1, top_k=10000, max_missing=0.10, chunk_size=8192,
                 pc_grid=PC_GRID, selected_fold_numbers=None, outer_splits=None,
                 inner_splits_by_fold=None):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    active_groups = groups if design == "clade_grouped" else None
    metrics, predictions, importances, tuning_rows, chosen_rows, selected_rows = [], [], [], [], [], []
    partitions = (repeated_outer_splits(y, n_splits, n_repeats, seed, design, active_groups)
                  if outer_splits is None else outer_splits)
    for number, repeat, fold, train, test in partitions:
        if selected_fold_numbers is not None and number not in selected_fold_numbers:
            continue
        train, test = np.asarray(train, dtype=int), np.asarray(test, dtype=int)
        if not np.array_equal(np.sort(np.r_[train, test]), np.arange(len(y))):
            raise ValueError("Outer training and test indices must partition all analysis samples")
        if len(np.unique(y[train])) != 2 or len(np.unique(y[test])) != 2:
            raise ValueError("Every outer fold must contain both classes for ROC-AUC evaluation")
        if active_groups is not None and set(active_groups[train]) & set(active_groups[test]):
            raise ValueError("Training and test groups overlap")
        fold_seed = seed + number
        train_groups = active_groups[train] if active_groups is not None else None
        fold_info = {"repeat": repeat, "fold": fold, "fold_number": number}
        explicit_inner = None
        if inner_splits_by_fold is not None:
            if number not in inner_splits_by_fold:
                raise ValueError("Missing explicit inner partitions for outer fold " + str(number))
            explicit_inner = inner_splits_by_fold[number]
        if input_type == "pc":
            x_train = np.asarray(x[train], dtype=np.float32)
            tuned = tune_pcs(x_train, y[train], train_groups, fold_seed, pc_grid=pc_grid,
                             n_estimators=n_estimators, n_jobs=n_jobs, n_splits=inner_folds,
                             explicit_splits=explicit_inner)
            n_features = tuned["n_pcs"]
            projector = fit_pca_projector(x_train, n_features, fold_seed)
            z_train = projector.transform(x_train)
            z_test = projector.transform(np.asarray(x[test], dtype=np.float32))
            clf = forest(tuned["params"], fold_seed, n_estimators, n_jobs)
            clf.fit(z_train, y[train])
            probability = clf.predict_proba(z_test)[:, 1]
            importance = pd.DataFrame({"feature_idx": np.arange(n_features),
                                       "feature": [f"PC{i+1}" for i in range(n_features)],
                                       "gini_importance": clf.feature_importances_})
        else:
            scores = select_features(x, y, train, top_k, max_missing, chunk_size)
            selected = scores.variant_idx.to_numpy(dtype=int)
            z = np.asarray(x[:, selected], dtype=np.float32)
            tuned = tune_genotypes(z[train], y[train], train_groups, fold_seed,
                                   n_estimators=n_estimators, n_jobs=n_jobs, n_splits=inner_folds,
                                   explicit_splits=explicit_inner)
            clf = genotype_model(tuned["params"], fold_seed, n_estimators, n_jobs)
            clf.fit(z[train], y[train])
            probability = clf.predict_proba(z[test])[:, 1]
            n_features = len(selected)
            selected_table = scores.merge(features, on="variant_idx", validate="one_to_one")
            selected_table["selection_rank"] = np.arange(1, n_features + 1)
            selected_rows.append(selected_table.assign(**fold_info))
            importance = features.set_index("variant_idx").loc[selected].reset_index()
            importance["gini_importance"] = clf.named_steps["rf"].feature_importances_
        predicted = (probability >= tuned["threshold"]).astype(int)
        common = {"cv_design": DESIGN_LABELS[design], "input_type": INPUT_LABELS[input_type],
                  "input_key": input_type, "design_key": design, **fold_info}
        metrics.append({**common, "threshold": tuned["threshold"], "n_train": len(train),
                        "n_test": len(test), "n_features": n_features,
                        **metric_values(y[test], predicted, probability)})
        for index, sample_index in enumerate(test):
            predictions.append({**common, "sample": samples[sample_index],
                                "true_label": "R" if y[sample_index] else "S",
                                "predicted_label": "R" if predicted[index] else "S",
                                "probability_R": float(probability[index]), "threshold": tuned["threshold"]})
        importances.append(importance.assign(**fold_info))
        tuning_rows.extend({**fold_info, **entry} for entry in tuned["tuning"])
        chosen_rows.append({**fold_info, "candidate_number": tuned["candidate_number"], **tuned["params"],
                            **({"n_pcs": tuned["n_pcs"]} if input_type == "pc" else {}),
                            "mean_inner_roc_auc": tuned["mean_inner_roc_auc"],
                            "selected_threshold": tuned["threshold"],
                            "inner_balanced_accuracy_at_threshold": tuned["inner_balanced_accuracy"]})
        print(f"{design}/{input_type}: repeat {repeat}, fold {fold}, AUC {metrics[-1]['roc_auc']:.3f}", flush=True)
    if not metrics:
        raise ValueError("No outer folds were selected")
    outputs = {"cv_fold_metrics.tsv": pd.DataFrame(metrics),
               "cv_fold_predictions.tsv": pd.DataFrame(predictions),
               "feature_importance_by_fold.tsv": pd.concat(importances, ignore_index=True),
               "inner_cv_tuning_results.tsv": pd.DataFrame(tuning_rows),
               "selected_hyperparameters_by_fold.tsv": pd.DataFrame(chosen_rows)}
    if selected_rows:
        outputs["cv_selected_features_by_fold.tsv"] = pd.concat(selected_rows, ignore_index=True)
    for name, table in outputs.items():
        table.to_csv(destination / name, sep="\t", index=False)
    info = {"input_type": input_type, "design": design, "n_samples": len(y),
            "n_resistant": int(y.sum()), "n_susceptible": int((1-y).sum()),
            "n_input_features": x.shape[1], "outer_folds": n_splits, "outer_repeats": n_repeats,
            "inner_folds": inner_folds, "seed": seed, "trees": n_estimators,
            "top_k": top_k, "maximum_training_missingness": max_missing,
            "feature_selection_scope": "outer training only, fixed across inner folds" if input_type != "pc" else "PCA preprocessing refitted within every training split",
            "rf_grid": RF_GRID if input_type != "pc" else PC_RF_GRID,
            "pc_grid": pc_grid if input_type == "pc" else None,
            "completed_fold_numbers": [row["fold_number"] for row in metrics],
            "outer_partition_source": "explicit supplied partitions" if outer_splits is not None else "seeded splitter",
            "inner_partition_source": "explicit supplied partitions" if inner_splits_by_fold is not None else "seeded splitter",
            "versions": {"numpy": np.__version__, "pandas": pd.__version__,
                         "scipy": scipy.__version__, "scikit_learn": sklearn.__version__}}
    (destination / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="model_inputs.tsv from input preparation")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input-types", nargs="+", choices=list(INPUT_LABELS), default=list(INPUT_LABELS))
    parser.add_argument("--designs", nargs="+", choices=list(DESIGN_LABELS), default=list(DESIGN_LABELS))
    parser.add_argument("--sample-column", default="sample")
    parser.add_argument("--phenotype-column", default="rbin_flu")
    parser.add_argument("--group-column", default="cluster")
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--n-repeats", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--trees", type=int, default=500)
    parser.add_argument("--jobs", type=int, default=-1)
    parser.add_argument("--top-k", type=int, default=10000)
    parser.add_argument("--max-missing", type=float, default=0.10)
    parser.add_argument("--chunk-size", type=int, default=8192)
    parser.add_argument("--pc-components", type=int, nargs="+", default=PC_GRID)
    parser.add_argument("--fold-numbers", type=int, nargs="+", help="Restrict execution to listed outer folds for a partial run")
    args = parser.parse_args()
    manifest = pd.read_csv(args.inputs, sep="\t")
    if not {"input_type", "cache_dir"}.issubset(manifest) or manifest.input_type.duplicated().any():
        parser.error("Input manifest must contain distinct input_type and cache_dir rows")
    args.outdir.mkdir(parents=True, exist_ok=False)
    for input_type in args.input_types:
        source = manifest.loc[manifest.input_type.eq(input_type)]
        if len(source) != 1:
            parser.error("Missing cache for " + input_type)
        path = Path(source.iloc[0].cache_dir).expanduser()
        if not path.is_absolute():
            path = args.inputs.resolve().parent / path
        x, features, samples = load_cache(path)
        y, groups = aligned_metadata(args.metadata, samples, args.sample_column, args.phenotype_column,
                                     args.group_column if "clade_grouped" in args.designs else None)
        for design in args.designs:
            run_analysis(x, features, samples, y, groups, input_type, design,
                         args.outdir / f"{design}_{input_type}", args.n_splits, args.n_repeats,
                         args.seed, args.inner_folds, args.trees, args.jobs, args.top_k,
                         args.max_missing, args.chunk_size, args.pc_components, args.fold_numbers)


if __name__ == "__main__":
    main()
