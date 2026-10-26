#!/usr/bin/env python3
"""Refit the six prediction analyses after fluconazole SDD exclusion.

Original outer and inner partitions are intersected with the retained samples.
Each outer fold is independently checkpointed for concurrent workers.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd

from rf_common import (PC_GRID, PC_RF_GRID, RF_GRID, aligned_metadata, load_cache,
                       metric_values, validate_inner_splits)

HERE = Path(__file__).resolve().parent
SUBMISSION = HERE.parents[1]
SPEC = importlib.util.spec_from_file_location("sdd_rf_driver", HERE / "02_run_random_forest.py")
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare(args):
    original_path = args.source_root / "run_configuration.json"
    original = json.loads(original_path.read_text())
    tasks = [t for t in original["tasks"] if t["scope"] == "prediction"]
    expected = {f"{design}_{kind}" for design in ["stratified", "clade_grouped"]
                for kind in ["variant", "gene", "pc"]}
    if {t["task"] for t in tasks} != expected:
        raise ValueError("The original configuration must contain all six prediction tasks")
    original_samples = tasks[0]["samples"]
    if any(t["samples"] != original_samples for t in tasks):
        raise ValueError("Original prediction tasks have different cohort sample orders")
    metadata = pd.read_csv(original["metadata"], sep="\t", dtype={"sample": str})
    if metadata["sample"].duplicated().any():
        raise ValueError("Original metadata contain duplicate sample identifiers")
    metadata = metadata.set_index("sample").loc[original_samples].reset_index()
    mic = pd.read_csv(args.mic_metadata, sep="\t", dtype={"sample": str})
    if mic["sample"].duplicated().any():
        raise ValueError("MIC metadata contain duplicate sample identifiers")
    values = pd.to_numeric(mic.set_index("sample").loc[original_samples, "Fluconazole"], errors="raise")
    if not np.isfinite(values).all() or not values.gt(0).all():
        raise ValueError("Every original cohort isolate must have a positive finite fluconazole MIC")
    excluded = values.gt(2) & values.le(4)
    keep_samples = [s for s in original_samples if not excluded.loc[s]]
    retained = metadata.loc[metadata["sample"].isin(keep_samples)].copy()
    y, groups = aligned_metadata(retained, keep_samples)
    removed = metadata.loc[metadata["sample"].isin(values.index[excluded])].copy()
    removed["Fluconazole_MIC"] = removed["sample"].map(values)
    if not removed.rbin_flu.eq("R").all():
        raise ValueError("Excluded SDD-range samples must carry the original R label")
    source = args.source_root / "source_data"
    outer_path, inner_path = source / "outer_sample_assignments.tsv.gz", source / "inner_sample_assignments.tsv.gz"
    outer = pd.read_csv(outer_path, sep="\t", dtype={"sample": str})
    inner = pd.read_csv(inner_path, sep="\t", dtype={"sample": str})
    out = args.outdir
    if (out / "run_configuration.json").exists():
        raise FileExistsError("A sensitivity configuration already exists: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    sources = out / "source_data"
    sources.mkdir(exist_ok=True)
    (out / "partitions").mkdir(exist_ok=True)
    retained.to_csv(sources / "retained_metadata.tsv", sep="\t", index=False)
    removed.to_csv(sources / "excluded_samples.tsv", sep="\t", index=False)
    lookup = {s: i for i, s in enumerate(keep_samples)}
    all_outer, all_inner, new_tasks = [], [], []
    for task in tasks:
        name = task["task"]
        saved = outer.loc[outer.task.eq(name)]
        parts, test_counts = [], np.zeros((task["n_repeats"], len(keep_samples)), dtype=int)
        for number in range(1, task["n_splits"] * task["n_repeats"] + 1):
            frame = saved.loc[saved.fold_number.eq(number)]
            if len(frame) != len(original_samples) or set(frame["sample"]) != set(original_samples):
                raise ValueError(f"Original outer assignments are incomplete: {name}/{number}")
            repeat, fold = int(frame.repeat.iloc[0]), int(frame.fold.iloc[0])
            indices = {}
            for label in ["Training", "Test"]:
                rows = frame.loc[frame.partition.eq(label) & frame["sample"].isin(lookup)].copy()
                idx = np.array([lookup[s] for s in rows["sample"]], dtype=int)
                indices[label] = idx
                all_outer.append(rows)
            train, test = indices["Training"], indices["Test"]
            if not np.array_equal(np.sort(np.r_[train, test]), np.arange(len(keep_samples))):
                raise ValueError(f"Retained outer partitions overlap or omit samples: {name}/{number}")
            if any(len(np.unique(y[idx])) != 2 for idx in [train, test]):
                raise ValueError(f"Outer partition lacks both classes: {name}/{number}")
            grouped = task["design"] == "clade_grouped"
            if grouped and set(groups[train]) & set(groups[test]):
                raise ValueError(f"Outer partition contains group leakage: {name}/{number}")
            test_counts[repeat - 1, test] += 1
            local = {keep_samples[idx]: j for j, idx in enumerate(train)}
            original_inner = inner.loc[inner.task.eq(name) & inner.fold_number.eq(number)]
            splits = []
            for inn in range(1, original["inner_folds"] + 1):
                saved_inner = original_inner.loc[original_inner.inner_fold.eq(inn)]
                before = set(frame.loc[frame.partition.eq("Training"), "sample"])
                if len(saved_inner) != len(before) or set(saved_inner["sample"]) != before:
                    raise ValueError(f"Original inner assignments are incomplete: {name}/{number}/{inn}")
                pair = []
                for label in ["Training", "Validation"]:
                    rows = saved_inner.loc[saved_inner.partition.eq(label) & saved_inner["sample"].isin(local)].copy()
                    pair.append(np.array([local[s] for s in rows["sample"]], dtype=int))
                    all_inner.append(rows)
                splits.append(tuple(pair))
            splits = validate_inner_splits(y[train], groups[train] if grouped else None,
                                           splits, original["inner_folds"])
            if task["input_type"] == "pc" and min(len(a) for a, _ in splits) < max(PC_GRID):
                raise ValueError("The retained inner-training samples cannot fit the original PC grid")
            parts.append(dict(number=number, repeat=repeat, fold=fold,
                              train=train.tolist(), test=test.tolist(),
                              inner=[dict(train=a.tolist(), valid=b.tolist()) for a, b in splits]))
        if not np.all(test_counts == 1):
            raise ValueError("Each retained sample must have one outer test prediction per repeat")
        saved_info = json.loads((args.source_root / "results" / name / "folds/01/run_info.json").read_text())
        if saved_info["rf_grid"] != (PC_RF_GRID if task["input_type"] == "pc" else RF_GRID):
            raise ValueError("Current candidate forest grid differs from the original run")
        if task["input_type"] == "pc" and saved_info["pc_grid"] != PC_GRID:
            raise ValueError("Current PC grid differs from the original run")
        save_json(out / "partitions" / (name + ".json"), parts)
        new_tasks.append({**task, "samples": keep_samples})
    pd.concat(all_outer, ignore_index=True).to_csv(sources / "outer_sample_assignments.tsv.gz", sep="\t", index=False)
    pd.concat(all_inner, ignore_index=True).to_csv(sources / "inner_sample_assignments.tsv.gz", sep="\t", index=False)
    files = [original_path, outer_path, inner_path, Path(original["metadata"]), args.mic_metadata]
    config = {k: original[k] for k in ["input_caches", "seed", "trees", "inner_folds", "top_k", "max_missing"]}
    config.update(source_root=str(args.source_root.resolve()), metadata=str((sources / "retained_metadata.tsv").resolve()),
                  mic_metadata=str(args.mic_metadata.resolve()), exclusion="Fluconazole MIC >2 and <=4 micrograms/mL",
                  tasks=new_tasks, original_n=len(original_samples), retained_n=len(keep_samples),
                  excluded_n=len(removed), retained_resistant=int(y.sum()), retained_susceptible=int((1-y).sum()),
                  partitions="Original outer and inner assignments intersected with retained samples",
                  gene_variant_frequency_scope="Original 871-isolate binary-phenotyped cohort",
                  pc_marker_scope="Original shared nuclear LD-pruned marker panel",
                  created_utc=datetime.now(timezone.utc).isoformat(),
                  sources=[dict(path=str(p.resolve()), sha256=digest(p)) for p in files])
    save_json(out / "run_configuration.json", config)
    print(json.dumps({"prepared": len(new_tasks), "outer_models": 150,
                      "original_n": len(original_samples), "excluded_n": len(removed),
                      "retained_n": len(keep_samples), "R": int(y.sum()), "S": int((1-y).sum())}), flush=True)


def work(args):
    root = args.outdir
    config = json.loads((root / "run_configuration.json").read_text())
    metadata = pd.read_csv(config["metadata"], sep="\t", dtype={"sample": str})
    loaded, active, finished = {}, {}, 0
    for task in config["tasks"]:
        if args.only_input and task["input_type"] != args.only_input:
            continue
        if args.only_design and task["design"] != args.only_design:
            continue
        target = root / "results" / task["task"] / "folds"
        target.mkdir(parents=True, exist_ok=True)
        partitions = json.loads((root / "partitions" / (task["task"] + ".json")).read_text())
        for part in partitions:
            number = part["number"]
            if args.fold_numbers and number not in args.fold_numbers:
                continue
            fold = target / f"{number:02d}"
            if (fold / "run_info.json").exists():
                continue
            lock = target / f".{number:02d}.lock"
            try:
                lock.mkdir()
            except FileExistsError:
                continue
            if (fold / "run_info.json").exists():
                lock.rmdir()
                continue
            save_json(lock / "owner.json", {"pid": os.getpid(), "started": time.time()})
            started = time.monotonic()
            try:
                if fold.exists():
                    raise RuntimeError("Inspect incomplete fold before restarting: " + str(fold))
                kind = task["input_type"]
                if kind not in loaded:
                    loaded[kind] = load_cache(config["input_caches"][kind])
                if kind not in active:
                    matrix, features, all_samples = loaded[kind]
                    lookup = {s: i for i, s in enumerate(all_samples)}
                    x = np.asarray(matrix[[lookup[s] for s in task["samples"]]], dtype=np.int8)
                    y, groups = aligned_metadata(metadata, task["samples"])
                    active[kind] = (x, features, y, groups)
                x, features, y, groups = active[kind]
                outer = [(number, part["repeat"], part["fold"], part["train"], part["test"])]
                inner = {number: [(s["train"], s["valid"]) for s in part["inner"]]}
                print(datetime.now().isoformat(timespec="seconds"), task["task"], number, "start", flush=True)
                DRIVER.run_analysis(x, features, task["samples"], y, groups, kind, task["design"], fold,
                                    n_splits=task["n_splits"], n_repeats=task["n_repeats"], seed=config["seed"],
                                    inner_folds=config["inner_folds"], n_estimators=config["trees"], n_jobs=args.jobs,
                                    top_k=config["top_k"], max_missing=config["max_missing"],
                                    selected_fold_numbers=[number], outer_splits=outer, inner_splits_by_fold=inner)
                save_json(fold / "sensitivity_timing.json", {"elapsed_seconds": time.monotonic() - started})
                print(datetime.now().isoformat(timespec="seconds"), task["task"], number, "complete", flush=True)
            except BaseException:
                (lock / "error.txt").write_text(traceback.format_exc())
                raise
            else:
                (lock / "owner.json").unlink()
                lock.rmdir()
                finished += 1
                if args.max_folds and finished >= args.max_folds:
                    return
    print("Worker has no unclaimed selected outer folds", flush=True)


def summarize(args):
    config = json.loads((args.outdir / "run_configuration.json").read_text())
    metrics, predictions = [], []
    for task in config["tasks"]:
        destination = args.outdir / "results" / task["task"]
        folders = [destination / "folds" / f"{n:02d}" for n in range(1, task["n_splits"] * task["n_repeats"] + 1)]
        if not all((p / "run_info.json").exists() for p in folders):
            raise ValueError("Incomplete sensitivity task: " + task["task"])
        for filename, collector in [("cv_fold_metrics.tsv", metrics), ("cv_fold_predictions.tsv", predictions)]:
            table = pd.concat([pd.read_csv(p / filename, sep="\t") for p in folders], ignore_index=True)
            table["task"] = task["task"]
            table.to_csv(destination / filename, sep="\t", index=False)
            collector.append(table)
    metrics = pd.concat(metrics, ignore_index=True)
    predictions = pd.concat(predictions, ignore_index=True)
    expected = len(config["tasks"]) * config["retained_n"] * config["tasks"][0]["n_repeats"]
    if len(predictions) != expected or predictions.duplicated(["task", "repeat", "sample"]).any():
        raise ValueError("Sensitivity predictions have missing or duplicated held-out sample/repeat pairs")
    metrics.to_csv(args.outdir / "results/combined_cv_fold_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(args.outdir / "results/combined_cv_fold_predictions.tsv", sep="\t", index=False)
    original = pd.read_csv(Path(config["source_root"]) / "results/combined_cv_fold_predictions.tsv", sep="\t")
    keys = ["task", "fold_number", "repeat", "fold", "sample"]
    paired = predictions.merge(original, on=keys, how="left", suffixes=("_sensitivity", "_original"), validate="one_to_one")
    if paired.probability_R_original.isna().any() or not paired.true_label_sensitivity.eq(paired.true_label_original).all():
        raise ValueError("Original and sensitivity predictions do not share held-out samples and labels")
    paired.to_csv(args.outdir / "results/paired_retained_predictions.tsv.gz", sep="\t", index=False)
    comparisons = []
    for (task, number, repeat, fold), frame in paired.groupby(["task", "fold_number", "repeat", "fold"], sort=False):
        y = frame.true_label_original.eq("R").to_numpy(dtype=int)
        scores = {}
        for side in ["original", "sensitivity"]:
            scores[side] = metric_values(y, frame[f"predicted_label_{side}"].eq("R").to_numpy(dtype=int),
                                        frame[f"probability_R_{side}"].to_numpy())
        row = dict(task=task, fold_number=int(number), repeat=int(repeat), fold=int(fold), n_test=len(frame))
        for metric in scores["original"]:
            row[metric + "_original_retained"] = scores["original"][metric]
            row[metric + "_sensitivity"] = scores["sensitivity"][metric]
            row[metric + "_difference"] = scores["sensitivity"][metric] - scores["original"][metric]
        comparisons.append(row)
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(args.outdir / "results/paired_fold_comparison.tsv", sep="\t", index=False)
    score_columns = [c for c in comparison if c.endswith(("_original_retained", "_sensitivity", "_difference"))]
    comparison.groupby("task")[score_columns].mean().reset_index().to_csv(args.outdir / "results/paired_mean_comparison.tsv", sep="\t", index=False)
    config.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat())
    save_json(args.outdir / "run_configuration.json", config)
    print("All six sensitivity tasks complete; paired retained-sample comparisons saved", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "work", "summarize"])
    parser.add_argument("--source-root", type=Path, default=SUBMISSION / "Analysis/RF")
    parser.add_argument("--outdir", type=Path, default=SUBMISSION / "Analysis/RF_SDD_sensitivity_20261007")
    parser.add_argument("--mic-metadata", type=Path, default=SUBMISSION / "Main_figures_archive/02_fungamr/source_data/panel_A_and_B/all_metadata.tsv")
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--only-input", choices=["variant", "gene", "pc"])
    parser.add_argument("--only-design", choices=["stratified", "clade_grouped"])
    parser.add_argument("--fold-numbers", type=int, nargs="+")
    parser.add_argument("--max-folds", type=int)
    args = parser.parse_args()
    {"prepare": prepare, "work": work, "summarize": summarize}[args.action](args)


if __name__ == "__main__":
    main()
