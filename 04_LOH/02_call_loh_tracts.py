#!/usr/bin/env python3
"""Call low-heterozygosity tracts from fixed-window SNP counts."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


WINDOW_BP = 10_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True, help="10 kb window TSV or TSV.GZ")
    parser.add_argument("--out", required=True, help="Output tract TSV")
    parser.add_argument(
        "--samples",
        help="Optional final sample TSV containing a sample column",
    )
    parser.add_argument(
        "--copy-number-dir",
        help="Optional directory containing SAMPLE.tsv 10 kb copy-number files",
    )
    parser.add_argument("--minimum-callable-snps", type=int, default=10)
    parser.add_argument("--maximum-heterozygous-snps", type=int, default=10)
    parser.add_argument("--minimum-low-windows", type=int, default=3)
    parser.add_argument("--maximum-intervening-windows", type=int, default=3)
    parser.add_argument("--minimum-purity", type=float, default=0.80)
    parser.add_argument("--expected-samples", type=int, help="Optional sample-count assertion")
    parser.add_argument(
        "--exclude-contig-regex",
        default=r"(?i)chrM(?:_|$)",
        help="Regular expression for contigs excluded from nuclear LOH calling",
    )
    return parser.parse_args()


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def load_windows(path: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    require_columns(
        frame,
        {
            "sample",
            "contig",
            "window_start0",
            "window_end0",
            "N_call",
            "N_het",
        },
        "Window table",
    )
    for column in ("window_start0", "window_end0", "N_call", "N_het"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame["sample"] = frame["sample"].astype(str).str.strip()
    frame["contig"] = frame["contig"].astype(str).str.strip()
    if frame.duplicated(["sample", "contig", "window_start0"]).any():
        raise ValueError("Window coordinates are duplicated within a sample and contig.")
    return frame


def load_sample_set(path: str) -> set[str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if "sample" not in frame.columns:
        raise ValueError("The sample table requires a sample column.")
    samples = frame["sample"].str.strip().dropna().tolist()
    if len(samples) != len(set(samples)):
        raise ValueError("The sample table contains duplicate identifiers.")
    return set(samples)


def load_copy_number(directory: str | None, sample: str) -> pd.DataFrame | None:
    if not directory:
        return None
    path = Path(directory) / f"{sample}.tsv"
    if not path.is_file() or path.stat().st_size == 0:
        return None
    raw = pd.read_csv(path, sep="\t")
    names = {column.lower(): column for column in raw.columns}
    contig_column = next(
        (
            names[name]
            for name in ("contig", "chr", "chrom", "chromosome")
            if name in names
        ),
        None,
    )
    start_column = names.get("start")
    end_column = names.get("end")
    cn_column = next(
        (
            names[name]
            for name in ("cn", "cn_bin", "copy_number", "copynumber")
            if name in names
        ),
        None,
    )
    ploidy_column = names.get("baseline_ploidy")
    if None in (contig_column, start_column, end_column, cn_column, ploidy_column):
        raise ValueError(f"{path} lacks required copy-number columns.")
    frame = pd.DataFrame(
        {
            "contig": raw[contig_column].astype(str).str.strip(),
            "start": pd.to_numeric(raw[start_column], errors="coerce"),
            "end": pd.to_numeric(raw[end_column], errors="coerce"),
            "CN": pd.to_numeric(raw[cn_column], errors="coerce"),
            "baseline_ploidy": pd.to_numeric(raw[ploidy_column], errors="coerce"),
        }
    ).dropna()
    if frame.empty:
        return None
    if frame["start"].min() >= 1:
        frame["start0"] = frame["start"].astype(int) - 1
        frame["end0"] = frame["end"].astype(int)
    else:
        frame["start0"] = frame["start"].astype(int)
        frame["end0"] = frame["end"].astype(int)
    return frame[frame["end0"] > frame["start0"]].copy()


def qualifying_intervals(
    low: np.ndarray,
    start_index: int,
    end_index: int,
    minimum_low_windows: int,
    minimum_purity: float,
) -> list[tuple[int, int]]:
    """Select non-overlapping qualifying intervals in genomic order.

    At the earliest qualifying low window, retain the longest interval whose
    endpoints are low and whose low-window count and purity pass the thresholds.
    Continue beyond that interval, or advance to the next low window when no
    interval qualifies from the current start.
    """
    positions = np.flatnonzero(low[start_index : end_index + 1]) + start_index
    intervals = []
    first = 0
    while first + minimum_low_windows <= len(positions):
        last_positions = positions[first + minimum_low_windows - 1 :]
        counts = np.arange(minimum_low_windows, len(positions) - first + 1)
        purity = counts / (last_positions - positions[first] + 1)
        valid = np.flatnonzero(purity >= minimum_purity)
        if not len(valid):
            first += 1
            continue
        last = first + minimum_low_windows - 1 + int(valid[-1])
        intervals.append((int(positions[first]), int(positions[last])))
        first = last + 1
    return intervals


def call_contig(
    frame: pd.DataFrame,
    minimum_callable: int,
    maximum_heterozygous: int,
    minimum_low_windows: int,
    maximum_intervening: int,
    minimum_purity: float,
) -> list[dict[str, float | int]]:
    if minimum_callable < 0 or maximum_heterozygous < 0:
        raise ValueError("SNP-count thresholds must be nonnegative.")
    if minimum_low_windows < 1 or maximum_intervening < 0:
        raise ValueError("Minimum low windows must be positive and intervening windows nonnegative.")
    if not 0 < minimum_purity <= 1:
        raise ValueError("Minimum purity must be in (0, 1].")
    frame = frame.sort_values("window_start0").reset_index(drop=True)
    starts = frame["window_start0"].to_numpy(dtype=int)
    ends = frame["window_end0"].to_numpy(dtype=int)
    if np.any(ends <= starts):
        raise ValueError("Window end coordinates must exceed start coordinates.")
    if len(frame) > 1 and np.any(starts[1:] != ends[:-1]):
        raise ValueError("Windows must be contiguous within each sample and contig.")

    callable_counts = frame["N_call"].to_numpy(dtype=float)
    heterozygous_counts = frame["N_het"].to_numpy(dtype=float)
    low = (callable_counts >= minimum_callable) & (
        heterozygous_counts <= maximum_heterozygous
    )
    tracts: list[dict[str, float | int]] = []
    index = 0
    while index < len(frame):
        if not low[index]:
            index += 1
            continue
        segment_start = index
        last_low = index
        consecutive_nonlow = 0
        cursor = index + 1
        while cursor < len(frame):
            if low[cursor]:
                last_low = cursor
                consecutive_nonlow = 0
            else:
                consecutive_nonlow += 1
                if consecutive_nonlow > maximum_intervening:
                    break
            cursor += 1
        for start_index, end_index in qualifying_intervals(
            low, segment_start, last_low, minimum_low_windows, minimum_purity
        ):
            count_low = int(low[start_index : end_index + 1].sum())
            tracts.append(
                {
                    "start0": int(starts[start_index]),
                    "end0": int(ends[end_index]),
                    "length_bp": int(ends[end_index] - starts[start_index]),
                    "n_windows_region": int(end_index - start_index + 1),
                    "n_low_windows": count_low,
                    "purity": count_low / (end_index - start_index + 1),
                    "median_N_het": float(
                        np.median(heterozygous_counts[start_index : end_index + 1])
                    ),
                }
            )
        index = max(cursor, last_low + 1)
    return tracts


def annotate_copy_number(
    copy_number: pd.DataFrame | None, contig: str, start0: int, end0: int
) -> tuple[float, float]:
    if copy_number is None:
        return float("nan"), float("nan")
    baseline = float(np.median(copy_number["baseline_ploidy"]))
    overlap = copy_number[
        (copy_number["contig"] == contig)
        & (copy_number["start0"] < end0)
        & (copy_number["end0"] > start0)
    ]
    median_cn = float(np.median(overlap["CN"])) if not overlap.empty else float("nan")
    return baseline, median_cn


def main() -> None:
    args = parse_args()
    if not 0 < args.minimum_purity <= 1:
        raise ValueError("--minimum-purity must be in (0, 1].")
    windows = load_windows(args.windows)
    windows = windows[
        ~windows["contig"].str.contains(
            re.compile(args.exclude_contig_regex), regex=True
        )
    ].copy()
    if args.samples:
        sample_set = load_sample_set(args.samples)
        windows = windows[windows["sample"].isin(sample_set)].copy()
        if set(windows["sample"]) != sample_set:
            raise ValueError("Window data do not cover every requested sample.")
    samples = windows["sample"].drop_duplicates().tolist()
    if args.expected_samples and len(samples) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} samples but found {len(samples)}."
        )

    output: list[dict[str, object]] = []
    for sample, sample_frame in windows.groupby("sample", sort=False):
        copy_number = load_copy_number(args.copy_number_dir, sample)
        for contig, contig_frame in sample_frame.groupby("contig", sort=False):
            for tract in call_contig(
                contig_frame,
                args.minimum_callable_snps,
                args.maximum_heterozygous_snps,
                args.minimum_low_windows,
                args.maximum_intervening_windows,
                args.minimum_purity,
            ):
                baseline, median_cn = annotate_copy_number(
                    copy_number, contig, int(tract["start0"]), int(tract["end0"])
                )
                output.append(
                    {
                        "sample": sample,
                        "baseline_ploidy": baseline,
                        "contig": contig,
                        **tract,
                        "median_CN": median_cn,
                        "minimum_callable_snps": args.minimum_callable_snps,
                        "maximum_heterozygous_snps": args.maximum_heterozygous_snps,
                        "maximum_intervening_windows": args.maximum_intervening_windows,
                        "minimum_purity": args.minimum_purity,
                        "minimum_low_windows": args.minimum_low_windows,
                        "window_bp": WINDOW_BP,
                    }
                )
    columns = [
        "sample", "baseline_ploidy", "contig", "start0", "end0", "length_bp",
        "n_windows_region", "n_low_windows", "purity", "median_N_het", "median_CN",
        "minimum_callable_snps", "maximum_heterozygous_snps",
        "maximum_intervening_windows", "minimum_purity", "minimum_low_windows", "window_bp",
    ]
    pd.DataFrame(output, columns=columns).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
