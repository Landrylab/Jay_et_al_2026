#!/usr/bin/env python3
"""Summarize whole-chromosome and segmental copy-number events."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-cnv", required=True, help="Gene-by-sample CN matrix")
    parser.add_argument("--ploidy", required=True, help="Final sample ploidy table")
    parser.add_argument(
        "--chromosome-cn",
        required=True,
        help="Final sample-by-chromosome CN table",
    )
    parser.add_argument("--loci", required=True, help="Ordered gene-locus table")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--minimum-segment-bp", type=int, default=30_000)
    parser.add_argument("--expected-samples", type=int, default=1_153)
    return parser.parse_args()


def resolve_column(frame: pd.DataFrame, names: tuple[str, ...]) -> str:
    lowered = {column.lower(): column for column in frame.columns}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    raise ValueError(f"Missing required column; accepted names: {names}")


def chromosome_label(contig: str) -> str:
    token = contig.split("chr", 1)[-1].split("A_", 1)[0]
    return f"Chr{token}"


def format_copy_number(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def dominant_copy_number(counts: Counter[float]) -> float:
    highest = max(counts.values())
    tied = [value for value, count in counts.items() if count == highest]
    return sorted(tied, key=lambda value: (abs(value - 2), value), reverse=True)[0]


def load_inputs(args: argparse.Namespace):
    ploidy = pd.read_csv(args.ploidy, sep="\t")
    sample_column = resolve_column(ploidy, ("sample_id", "sample", "isolate"))
    ploidy_column = resolve_column(ploidy, ("ploidy_hat", "baseline_ploidy", "ploidy"))
    ploidy = ploidy.set_index(sample_column)
    baseline = pd.to_numeric(ploidy[ploidy_column], errors="raise")
    if len(baseline) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} ploidy calls but found {len(baseline)}."
        )

    chromosome_cn = pd.read_csv(args.chromosome_cn, sep="\t")
    chromosome_sample = resolve_column(
        chromosome_cn, ("sample_id", "sample", "isolate")
    )
    chromosome_cn = chromosome_cn.set_index(chromosome_sample)

    gene_cnv = pd.read_csv(args.gene_cnv, sep="\t", low_memory=False)
    gene_column = resolve_column(gene_cnv, ("gene_id", "gene", "locus_id"))
    gene_cnv = gene_cnv.set_index(gene_column)

    loci = pd.read_csv(args.loci, sep=None, engine="python")
    required_locus_columns = {
        "locus_id",
        "gene_ids",
        "chromosome",
        "start_1based_inclusive",
        "end_1based_inclusive",
    }
    missing = required_locus_columns.difference(loci.columns)
    if missing:
        raise ValueError(f"Locus table is missing: {sorted(missing)}")
    loci = loci.copy()
    loci["chromosome_label"] = loci.get(
        "chromosome_label", loci["chromosome"].map(chromosome_label)
    )
    loci = loci.sort_values(
        ["chromosome_label", "start_1based_inclusive", "end_1based_inclusive"]
    ).reset_index(drop=True)

    diploid_samples = [
        sample
        for sample, value in baseline.items()
        if value == 2
        and sample in chromosome_cn.index
        and sample in gene_cnv.columns
    ]
    if len(diploid_samples) != 1_127:
        raise ValueError(
            f"Expected 1,127 diploid samples but found {len(diploid_samples)}."
        )
    return baseline, chromosome_cn, gene_cnv, loci, diploid_samples


def locus_matrices(
    gene_cnv: pd.DataFrame,
    loci: pd.DataFrame,
    samples: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    state = np.zeros((len(loci), len(samples)), dtype=np.int8)
    altered_cn = np.full((len(loci), len(samples)), np.nan, dtype=float)
    baseline = np.full(len(samples), 2.0)

    for locus_index, row in loci.iterrows():
        genes = [
            gene
            for gene in str(row["gene_ids"]).split(";")
            if gene in gene_cnv.index
        ]
        if not genes:
            raise ValueError(f"No CN row found for locus {row['locus_id']}")
        values = gene_cnv.loc[genes, samples].to_numpy(dtype=float)
        maximum = np.nanmax(values, axis=0)
        minimum = np.nanmin(values, axis=0)
        gain = (maximum > baseline) & (minimum >= baseline)
        loss = (minimum < baseline) & (maximum <= baseline)
        state[locus_index, gain] = 1
        state[locus_index, loss] = -1
        altered_cn[locus_index, gain] = maximum[gain]
        altered_cn[locus_index, loss] = minimum[loss]
    return state, altered_cn


def call_segmental_events(
    chromosome_loci: pd.DataFrame,
    chromosome_state: np.ndarray,
    chromosome_cn: np.ndarray,
    samples: list[str],
    whole_samples: set[str],
    minimum_segment_bp: int,
) -> list[dict]:
    events: list[dict] = []
    for sample_index, sample in enumerate(samples):
        if sample in whole_samples:
            continue
        states = chromosome_state[:, sample_index]
        position = 0
        while position < len(states):
            direction_code = int(states[position])
            if direction_code == 0:
                position += 1
                continue
            run_end = position + 1
            while run_end < len(states) and int(states[run_end]) == direction_code:
                run_end += 1

            first = chromosome_loci.iloc[position]
            last = chromosome_loci.iloc[run_end - 1]
            start = int(first["start_1based_inclusive"])
            end = int(last["end_1based_inclusive"])
            length_bp = end - start + 1
            if length_bp >= minimum_segment_bp:
                values = chromosome_cn[position:run_end, sample_index]
                values = values[np.isfinite(values)]
                counts = Counter(float(value) for value in values)
                events.append(
                    {
                        "sample": sample,
                        "chromosome_label": str(first["chromosome_label"]),
                        "chromosome": str(first["chromosome"]),
                        "direction": "gain" if direction_code > 0 else "loss",
                        "start_1based_inclusive": start,
                        "end_1based_inclusive": end,
                        "length_bp": length_bp,
                        "n_consecutive_cnv_loci": run_end - position,
                        "dominant_copy_number": dominant_copy_number(counts),
                        "copy_numbers_observed": ";".join(
                            format_copy_number(value) for value in sorted(counts)
                        ),
                        "altered_locus_copy_number_counts": ";".join(
                            f"{format_copy_number(value)}:{counts[value]}"
                            for value in sorted(counts)
                        ),
                        "first_locus_id": str(first["locus_id"]),
                        "last_locus_id": str(last["locus_id"]),
                        "minimum_segment_bp": minimum_segment_bp,
                    }
                )
            position = run_end
    return events


def main() -> None:
    args = parse_args()
    if args.minimum_segment_bp < 1:
        raise ValueError("--minimum-segment-bp must be positive.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    baseline, chromosome_cn, gene_cnv, loci, samples = load_inputs(args)
    state, altered_cn = locus_matrices(gene_cnv, loci, samples)

    chromosome_order = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "ChrR"]
    event_rows: list[dict] = []
    sample_rows: list[dict] = []
    summary_rows: list[dict] = []

    for label in chromosome_order:
        locus_mask = loci["chromosome_label"].astype(str).eq(label).to_numpy()
        chromosome_loci = loci.loc[locus_mask].reset_index(drop=True)
        chromosome = str(chromosome_loci.iloc[0]["chromosome"])
        whole_calls: dict[str, str] = {}
        for sample in samples:
            copy_number = float(chromosome_cn.at[sample, chromosome])
            if copy_number != float(baseline.loc[sample]):
                whole_calls[sample] = "gain" if copy_number > 2 else "loss"
                sample_rows.append(
                    {
                        "sample": sample,
                        "chromosome_label": label,
                        "chromosome": chromosome,
                        "aneuploidy_type": "whole chromosome",
                        "direction": whole_calls[sample],
                        "baseline_ploidy": 2,
                        "aneuploid_copy_number": copy_number,
                        "minimum_segment_bp": args.minimum_segment_bp,
                    }
                )

        chromosome_events = call_segmental_events(
            chromosome_loci,
            state[locus_mask],
            altered_cn[locus_mask],
            samples,
            set(whole_calls),
            args.minimum_segment_bp,
        )
        event_rows.extend(chromosome_events)
        event_frame = pd.DataFrame(chromosome_events)
        segmental_samples: set[str] = set()
        if not event_frame.empty:
            for sample, group in event_frame.groupby("sample"):
                segmental_samples.add(sample)
                counts: Counter[float] = Counter()
                for encoded in group["altered_locus_copy_number_counts"]:
                    for item in str(encoded).split(";"):
                        value, count = item.split(":")
                        counts[float(value)] += int(count)
                directions = sorted(group["direction"].unique())
                sample_rows.append(
                    {
                        "sample": sample,
                        "chromosome_label": label,
                        "chromosome": chromosome,
                        "aneuploidy_type": "segmental",
                        "direction": directions[0] if len(directions) == 1 else "mixed",
                        "baseline_ploidy": 2,
                        "aneuploid_copy_number": dominant_copy_number(counts),
                        "minimum_segment_bp": args.minimum_segment_bp,
                    }
                )

        whole_samples = set(whole_calls)
        total = len(whole_samples) + len(segmental_samples)
        summary_rows.append(
            {
                "chromosome_label": label,
                "chromosome": chromosome,
                "diploid_samples_considered": len(samples),
                "whole_chromosome_isolates": len(whole_samples),
                "segmental_isolates_min_30kb": len(segmental_samples),
                "total_aneuploid_isolates": total,
                "whole_chromosome_gain_isolates": sum(
                    direction == "gain" for direction in whole_calls.values()
                ),
                "whole_chromosome_loss_isolates": sum(
                    direction == "loss" for direction in whole_calls.values()
                ),
                "minimum_segment_bp": args.minimum_segment_bp,
                "whole_chromosome_sample_names": ";".join(sorted(whole_samples)),
                "segmental_sample_names": ";".join(sorted(segmental_samples)),
            }
        )

    events = pd.DataFrame(event_rows)
    sample_calls = pd.DataFrame(sample_rows)
    summary = pd.DataFrame(summary_rows)
    state_summary = (
        sample_calls.groupby(
            ["chromosome_label", "aneuploidy_type", "aneuploid_copy_number"],
            as_index=False,
        )
        .size()
        .rename(columns={"size": "isolate_count"})
    )
    totals = summary.set_index("chromosome_label")["total_aneuploid_isolates"]
    state_summary["percent_of_aneuploid_isolates"] = state_summary.apply(
        lambda row: 100
        * row["isolate_count"]
        / totals.loc[row["chromosome_label"]],
        axis=1,
    )

    summary.to_csv(
        output_dir / "cnv_chromosome_aneuploidy_summary.tsv",
        sep="\t",
        index=False,
    )
    events.to_csv(
        output_dir / "cnv_segmental_events_min30kb.tsv",
        sep="\t",
        index=False,
    )
    sample_calls.to_csv(
        output_dir / "cnv_chromosome_event_calls.tsv",
        sep="\t",
        index=False,
    )
    state_summary.to_csv(
        output_dir / "cnv_chromosome_copy_number_summary.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
