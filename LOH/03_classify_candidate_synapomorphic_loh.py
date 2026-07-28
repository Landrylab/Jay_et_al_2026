#!/usr/bin/env python3
"""Flag LOH tracts shared among neighboring tips in phylogenetic display order."""

from __future__ import annotations

import argparse
import bisect

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tracts", required=True, help="LOH tract TSV")
    parser.add_argument(
        "--tip-order",
        required=True,
        help="Phylogenetic display-order TSV containing a sample column",
    )
    parser.add_argument("--out", required=True, help="Output tract TSV")
    parser.add_argument("--neighbor-radius", type=int, default=2)
    parser.add_argument("--minimum-interval-jaccard", type=float, default=0.20)
    parser.add_argument("--minimum-support-fraction", type=float, default=0.80)
    parser.add_argument("--expected-samples", type=int, default=1_153)
    return parser.parse_args()


def interval_jaccard(
    first_start: int, first_end: int, second_start: int, second_end: int
) -> float:
    intersection = max(
        0, min(first_end, second_end) - max(first_start, second_start)
    )
    union = max(first_end, second_end) - min(first_start, second_start)
    return intersection / union if union > 0 else 0.0


def has_match(
    intervals: list[tuple[int, int]],
    start0: int,
    end0: int,
    minimum_jaccard: float,
) -> bool:
    if not intervals:
        return False
    starts = [interval[0] for interval in intervals]
    index = max(0, bisect.bisect_left(starts, start0) - 1)
    while index < len(intervals) and intervals[index][0] < end0:
        other_start, other_end = intervals[index]
        if (
            interval_jaccard(start0, end0, other_start, other_end)
            >= minimum_jaccard
        ):
            return True
        index += 1
    return False


def main() -> None:
    args = parse_args()
    if args.neighbor_radius < 1:
        raise ValueError("--neighbor-radius must be at least one.")
    for value, name in (
        (args.minimum_interval_jaccard, "--minimum-interval-jaccard"),
        (args.minimum_support_fraction, "--minimum-support-fraction"),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between zero and one.")

    order = pd.read_csv(args.tip_order, sep="\t", dtype=str)
    if "sample" not in order.columns:
        raise ValueError("The tip-order table requires a sample column.")
    tips = order["sample"].str.strip().dropna().tolist()
    if len(tips) != len(set(tips)):
        raise ValueError("The tip-order table contains duplicate samples.")
    if args.expected_samples and len(tips) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} tips but found {len(tips)}."
        )

    tracts = pd.read_csv(args.tracts, sep="\t")
    required = {"sample", "contig", "start0", "end0"}
    missing = sorted(required.difference(tracts.columns))
    if missing:
        raise ValueError(f"Tract table is missing columns: {', '.join(missing)}")
    tracts["sample"] = tracts["sample"].astype(str).str.strip()
    tracts["start0"] = pd.to_numeric(tracts["start0"], errors="raise").astype(int)
    tracts["end0"] = pd.to_numeric(tracts["end0"], errors="raise").astype(int)
    tracts = tracts[tracts["sample"].isin(tips)].copy()

    positions = {sample: index for index, sample in enumerate(tips)}
    neighbors = {
        sample: [
            tips[cursor]
            for cursor in range(
                max(0, index - args.neighbor_radius),
                min(len(tips), index + args.neighbor_radius + 1),
            )
            if cursor != index
        ]
        for sample, index in positions.items()
    }
    indexed: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for (sample, contig), group in tracts.groupby(["sample", "contig"], sort=False):
        indexed[(sample, contig)] = sorted(
            zip(group["start0"].astype(int), group["end0"].astype(int))
        )

    support_counts: list[int] = []
    neighbor_counts: list[int] = []
    support_fractions: list[float] = []
    flags: list[bool] = []
    for row in tracts.itertuples(index=False):
        nearby = neighbors[row.sample]
        support = sum(
            has_match(
                indexed.get((neighbor, row.contig), []),
                int(row.start0),
                int(row.end0),
                args.minimum_interval_jaccard,
            )
            for neighbor in nearby
        )
        fraction = support / len(nearby) if nearby else 0.0
        support_counts.append(support)
        neighbor_counts.append(len(nearby))
        support_fractions.append(fraction)
        flags.append(fraction >= args.minimum_support_fraction)

    tracts["phylogenetic_neighbor_support"] = support_counts
    tracts["phylogenetic_neighbors_tested"] = neighbor_counts
    tracts["phylogenetic_neighbor_support_fraction"] = support_fractions
    tracts["candidate_synapomorphic"] = flags
    tracts["neighbor_radius"] = args.neighbor_radius
    tracts["minimum_interval_jaccard"] = args.minimum_interval_jaccard
    tracts["minimum_support_fraction"] = args.minimum_support_fraction
    tracts.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
