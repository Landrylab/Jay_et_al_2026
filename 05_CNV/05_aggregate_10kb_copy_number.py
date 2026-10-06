#!/usr/bin/env python3
"""Aggregate segment-assigned 1-kb bins into the 10-kb CNV figure track."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", required=True, help="Filtered 1-kb depth table")
    parser.add_argument("--segments", required=True, help="Segment copy-number table")
    parser.add_argument("--baseline-ploidy", required=True, type=float)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window-size", type=int, default=10_000)
    parser.add_argument("--minimum-directional-support", type=float, default=0.70)
    return parser.parse_args()


def assign_bin_copy_numbers(
    bins: pd.DataFrame,
    segments: pd.DataFrame,
    baseline: float,
) -> pd.DataFrame:
    result = bins.copy()
    result["copy_number_1kb"] = np.nan
    result["midpoint"] = (result["start"] + result["end"]) / 2

    for chromosome, chromosome_bins in result.groupby("chr", sort=False):
        chromosome_segments = segments.loc[
            segments["chr"].astype(str).eq(str(chromosome))
        ].sort_values(["start", "end"])
        if chromosome_segments.empty:
            continue
        starts = chromosome_segments["start"].to_numpy(dtype=float)
        ends = chromosome_segments["end"].to_numpy(dtype=float)
        states = chromosome_segments["integer_copy_number"].to_numpy(dtype=float)
        midpoints = chromosome_bins["midpoint"].to_numpy(dtype=float)
        indices = np.searchsorted(starts, midpoints, side="right") - 1
        clipped = np.clip(indices, 0, len(ends) - 1)
        valid = (
            (indices >= 0)
            & (indices < len(chromosome_segments))
            & (midpoints < ends[clipped])
        )
        if np.any(starts[1:] < ends[:-1]):
            raise ValueError(f"Overlapping copy-number segments on {chromosome}")
        assigned = np.full(len(chromosome_bins), np.nan, dtype=float)
        assigned[valid] = states[indices[valid]]
        result.loc[chromosome_bins.index, "copy_number_1kb"] = assigned
    return result


def aggregate_windows(
    bins: pd.DataFrame,
    baseline: float,
    window_size: int,
    support_fraction: float,
) -> pd.DataFrame:
    working = bins.copy()
    working["window_start0"] = (working["start"] // window_size) * window_size
    rows: list[dict[str, float | int | str]] = []

    for (chromosome, window_start0), group in working.groupby(
        ["chr", "window_start0"],
        sort=True,
    ):
        states = group["copy_number_1kb"].to_numpy(dtype=float)
        states = states[np.isfinite(states)]
        if not len(states):
            continue
        n_bins = len(states)
        n_above = int(np.sum(states > baseline))
        n_below = int(np.sum(states < baseline))
        threshold = int(math.ceil(support_fraction * n_bins))
        if n_above >= threshold:
            copy_number = int(np.rint(np.median(states[states > baseline])))
        elif n_below >= threshold:
            copy_number = int(np.rint(np.median(states[states < baseline])))
        else:
            copy_number = int(np.rint(baseline))

        rows.append(
            {
                "chr": str(chromosome),
                "start": int(window_start0) + 1,
                "end": int(min(window_start0 + window_size, group["end"].max())),
                "CN": copy_number,
                "baseline_ploidy": baseline,
                "n_bins": n_bins,
                "n_above": n_above,
                "n_below": n_below,
                "thresh": threshold,
            }
        )
    return pd.DataFrame(rows, columns=["chr", "start", "end", "CN", "baseline_ploidy", "n_bins", "n_above", "n_below", "thresh"])


def main() -> None:
    args = parse_args()
    if args.window_size < 1:
        raise ValueError("--window-size must be positive.")
    if not 0 < args.minimum_directional_support <= 1:
        raise ValueError("--minimum-directional-support must be in (0, 1].")
    if not np.isfinite(args.baseline_ploidy) or args.baseline_ploidy <= 0:
        raise ValueError("--baseline-ploidy must be finite and positive.")

    depth = pd.read_csv(args.depth, sep="\t")
    segments = pd.read_csv(args.segments, sep="\t")
    required_depth = {"chr", "start", "end"}
    required_segments = {"chr", "start", "end", "integer_copy_number"}
    if missing := required_depth.difference(depth.columns):
        raise ValueError(f"Depth table is missing: {sorted(missing)}")
    if missing := required_segments.difference(segments.columns):
        raise ValueError(f"Segment table is missing: {sorted(missing)}")

    for frame, columns in (
        (depth, ["start", "end"]),
        (segments, ["start", "end", "integer_copy_number"]),
    ):
        for column in columns:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
    depth = depth.loc[
        np.isfinite(depth["start"])
        & np.isfinite(depth["end"])
        & (depth["end"] > depth["start"])
    ].copy()
    segments = segments.loc[
        np.isfinite(segments["start"])
        & np.isfinite(segments["end"])
        & np.isfinite(segments["integer_copy_number"])
        & (segments["end"] > segments["start"])
    ].copy()

    assigned = assign_bin_copy_numbers(
        depth,
        segments,
        args.baseline_ploidy,
    )
    output = aggregate_windows(
        assigned,
        args.baseline_ploidy,
        args.window_size,
        args.minimum_directional_support,
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
