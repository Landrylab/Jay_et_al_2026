#!/usr/bin/env python3
"""Validate reviewed per-sample calls and assemble cohort tables and 10-kb tracks."""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


def load_aggregator():
    path = Path(__file__).with_name("05_aggregate_10kb_copy_number.py")
    spec = importlib.util.spec_from_file_location("aggregate_copy_number", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def path_relative(value, manifest):
    path = Path(value)
    return path if path.is_absolute() else manifest.parent / path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="TSV: sample_id, prefix, depth; paths relative to manifest")
    parser.add_argument("--review", required=True, help="TSV: sample_id, baseline_ploidy, reviewer, evidence")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--window-size", type=int, default=10000)
    parser.add_argument("--minimum-directional-support", type=float, default=0.7)
    args = parser.parse_args()
    if args.window_size <= 0 or not 0 < args.minimum_directional_support <= 1:
        parser.error("Invalid aggregation thresholds.")
    manifest_path = Path(args.manifest).resolve()
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    review = pd.read_csv(args.review, sep="\t", dtype=str, keep_default_na=False)
    for frame, columns in ((manifest, {"sample_id", "prefix", "depth"}), (review, {"sample_id", "baseline_ploidy", "reviewer", "evidence"})):
        if not columns.issubset(frame.columns) or frame.empty:
            raise ValueError(f"Input table requires non-empty rows and columns {sorted(columns)}")
        if frame["sample_id"].duplicated().any() or (frame[list(columns)] == "").any().any():
            raise ValueError("Duplicate samples or blank required cells.")
    if set(manifest.sample_id) != set(review.sample_id):
        raise ValueError("Each manifested sample requires one review record.")
    review = review.set_index("sample_id")
    aggregator = load_aggregator()
    baseline_rows, chromosome_columns, gene_columns, collected = [], {}, {}, []
    for row in manifest.itertuples(index=False):
        if Path(row.sample_id).name != row.sample_id or row.sample_id in (".", ".."):
            raise ValueError("Sample IDs must be safe single-component filenames.")
        prefix = path_relative(row.prefix, manifest_path)
        ploidy = pd.read_csv(str(prefix) + ".ploidy.tsv", sep="\t")
        chromosome = pd.read_csv(str(prefix) + ".chromosome_cn.tsv", sep="\t")
        gene = pd.read_csv(str(prefix) + ".gene_cn.tsv", sep="\t")
        segments = pd.read_csv(str(prefix) + ".segments.tsv", sep="\t")
        depth = pd.read_csv(path_relative(row.depth, manifest_path), sep="\t")
        baseline = float(review.at[row.sample_id, "baseline_ploidy"])
        if baseline not in (2, 3, 4, 6) or len(ploidy) != 1 or float(ploidy.baseline_ploidy.iloc[0]) != baseline:
            raise ValueError(f"{row.sample_id}: rerun 04_segment_and_assign_copy_number.R with --baseline-ploidy matching the review record; downstream states must be reassigned.")
        if chromosome.chr.duplicated().any() or gene.gene_id.duplicated().any():
            raise ValueError(f"{row.sample_id}: duplicate chromosome or gene records.")
        if not np.all(chromosome.baseline_ploidy == baseline):
            raise ValueError(f"{row.sample_id}: chromosome calls use a different baseline.")
        chromosome_values = chromosome.set_index("chr").chromosome_copy_number
        gene_values = gene.set_index("gene_id").integer_copy_number
        if chromosome_columns and set(chromosome_values.index) != set(next(iter(chromosome_columns.values())).index):
            raise ValueError("Chromosome sets differ between samples.")
        if gene_columns and set(gene_values.index) != set(next(iter(gene_columns.values())).index):
            raise ValueError("Gene sets differ between samples.")
        if not np.isfinite(chromosome_values).all():
            raise ValueError(f"{row.sample_id}: unresolved chromosome calls.")
        chromosome_columns[row.sample_id] = chromosome_values
        gene_columns[row.sample_id] = gene_values
        baseline_rows.append({"sample_id": row.sample_id, "baseline_ploidy": int(baseline), "candidate_baseline_ploidy": ploidy.get("candidate_baseline_ploidy", ploidy.baseline_ploidy).iloc[0], "assignment_method": "manual_review"})
        bins = aggregator.assign_bin_copy_numbers(depth, segments, baseline)
        windows = aggregator.aggregate_windows(bins, baseline, args.window_size, args.minimum_directional_support)
        collected.append((row.sample_id, windows))
    out = Path(args.output_dir)
    (out / "binCN_10kb").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(baseline_rows).to_csv(out / "baseline_ploidy.tsv", sep="\t", index=False)
    chromosomes = pd.DataFrame(chromosome_columns).T
    chromosomes.index.name = "sample_id"
    chromosomes.to_csv(out / "chromosome_cn.tsv", sep="\t", na_rep="NA")
    genes = pd.DataFrame(gene_columns)
    genes.index.name = "gene_id"
    genes.to_csv(out / "gene_cn.tsv", sep="\t", na_rep="NA")
    review.to_csv(out / "manual_review.tsv", sep="\t")
    for sample, windows in collected:
        windows.to_csv(out / "binCN_10kb" / f"{sample}.tsv", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()
