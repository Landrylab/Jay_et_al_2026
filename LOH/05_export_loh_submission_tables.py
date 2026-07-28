#!/usr/bin/env python3
"""Export compact LOH supplementary tables and Figure 1F source data."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True, help="10 kb window TSV or TSV.GZ")
    parser.add_argument(
        "--tracts",
        required=True,
        help="LOH tracts with candidate_synapomorphic classifications",
    )
    parser.add_argument("--resolved", required=True, help="Clade-resolved LOH TSV")
    parser.add_argument("--samples", required=True, help="Final sample-order TSV")
    parser.add_argument("--fai", required=True, help="Reference FASTA index")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--source-data-out",
        help="Optional path for the compressed Figure 1F window table",
    )
    parser.add_argument("--expected-samples", type=int, default=1_153)
    return parser.parse_args()


def load_samples(path: str) -> list[str]:
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if "sample" not in frame.columns:
        raise ValueError("The final sample-order table requires a sample column.")
    samples = frame["sample"].str.strip().dropna().tolist()
    if len(samples) != len(set(samples)):
        raise ValueError("The final sample-order table contains duplicates.")
    return samples


def merged_length(frame: pd.DataFrame, contig_lengths: dict[str, int]) -> int:
    if frame.empty:
        return 0
    total = 0
    for _, part in frame.groupby("contig", sort=False):
        contig = str(part["contig"].iloc[0])
        intervals = sorted(
            (
                int(start0),
                min(int(end0), contig_lengths[contig]),
            )
            for start0, end0 in zip(part["start0"], part["end0"])
            if min(int(end0), contig_lengths[contig]) > int(start0)
        )
        if not intervals:
            continue
        current_start, current_end = intervals[0]
        for start0, end0 in intervals[1:]:
            if start0 <= current_end:
                current_end = max(current_end, end0)
            else:
                total += current_end - current_start
                current_start, current_end = start0, end0
        total += current_end - current_start
    return total


def main() -> None:
    args = parse_args()
    output_directory = Path(args.out_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    source_data_path = (
        Path(args.source_data_out)
        if args.source_data_out
        else output_directory / "Figure_1F_LOH_heterozygous_SNP_windows.tsv.gz"
    )
    source_data_path.parent.mkdir(parents=True, exist_ok=True)
    samples = load_samples(args.samples)
    if args.expected_samples and len(samples) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} samples but found {len(samples)}."
        )
    sample_set = set(samples)
    order = {sample: index for index, sample in enumerate(samples)}
    pd.DataFrame(
        {
            "phylogenetic_display_order": range(1, len(samples) + 1),
            "sample": samples,
        }
    ).to_csv(
        output_directory / "loh_phylogenetic_display_order.tsv",
        sep="\t",
        index=False,
    )

    fai = pd.read_csv(
        args.fai,
        sep="\t",
        header=None,
        names=["contig", "length", "offset", "line_bases", "line_width"],
    )
    nuclear = fai[
        ~fai["contig"].str.contains(r"(?i)chrM(?:_|$)", regex=True)
    ].copy()
    contig_set = set(nuclear["contig"])
    contig_lengths = dict(zip(nuclear["contig"], nuclear["length"].astype(int)))
    nuclear_total_bp = int(nuclear["length"].sum())

    windows = pd.read_csv(args.windows, sep="\t")
    windows = windows[
        windows["sample"].isin(sample_set) & windows["contig"].isin(contig_set)
    ].copy()
    if set(windows["sample"]) != sample_set:
        raise ValueError("Window data do not cover every final sample.")
    if "N_hom" not in windows.columns:
        windows["N_hom"] = windows["N_call"] - windows["N_het"]
    windows["HET_FRAC"] = windows["N_het"].div(
        windows["N_call"].where(windows["N_call"] > 0)
    )
    windows["_sample_order"] = windows["sample"].map(order)
    windows = windows.sort_values(
        ["_sample_order", "contig", "window_start0"]
    ).drop(columns="_sample_order")
    source_columns = [
        "sample",
        "contig",
        "window_start0",
        "window_end0",
        "N_call",
        "N_het",
        "N_hom",
        "HET_FRAC",
    ]
    windows[source_columns].to_csv(
        source_data_path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    tracts = pd.read_csv(args.tracts, sep="\t")
    tracts = tracts[
        tracts["sample"].isin(sample_set) & tracts["contig"].isin(contig_set)
    ].copy()
    tracts["_sample_order"] = tracts["sample"].map(order)
    tracts = tracts.sort_values(
        ["_sample_order", "contig", "start0", "end0"]
    ).drop(columns="_sample_order")
    if tracts.duplicated(["sample", "contig", "start0", "end0"]).any():
        raise ValueError("The tract table contains duplicate interval keys.")
    tracts.to_csv(
        output_directory / "loh_tracts_candidate_synapomorphic.tsv",
        sep="\t",
        index=False,
    )

    summaries: list[dict[str, object]] = []
    for sample in samples:
        sample_tracts = tracts[tracts["sample"] == sample]
        candidate = sample_tracts[sample_tracts["candidate_synapomorphic"].astype(bool)]
        total_loh_bp = merged_length(sample_tracts, contig_lengths)
        candidate_loh_bp = merged_length(candidate, contig_lengths)
        summaries.append(
            {
                "sample": sample,
                "n_loh_tracts": len(sample_tracts),
                "total_loh_bp_nonoverlapping": total_loh_bp,
                "fraction_nuclear_genome_loh": (
                    total_loh_bp / nuclear_total_bp
                ),
                "median_loh_length_bp": sample_tracts["length_bp"].median(),
                "n_candidate_synapomorphic_tracts": len(candidate),
                "candidate_synapomorphic_bp_nonoverlapping": candidate_loh_bp,
                "fraction_nuclear_genome_candidate_synapomorphic_loh": (
                    candidate_loh_bp / nuclear_total_bp
                ),
            }
        )
    sample_summary = pd.DataFrame(summaries)
    sample_summary.to_csv(
        output_directory / "loh_sample_summary.tsv", sep="\t", index=False
    )

    resolved = pd.read_csv(args.resolved, sep="\t")
    resolved.to_csv(
        output_directory / "loh_resolved_clade_regions.tsv",
        sep="\t",
        index=False,
    )

    parameters = pd.DataFrame(
        [
            ("cohort_samples", 1_153, "isolates"),
            ("window_size", 10_000, "bp"),
            ("minimum_callable_snps_per_window", 10, "SNP positions"),
            ("maximum_heterozygous_snps_in_low_window", 10, "SNP positions"),
            ("minimum_low_windows_per_tract", 3, "windows"),
            ("maximum_consecutive_intervening_windows", 3, "windows"),
            ("minimum_tract_purity", 0.80, "fraction"),
            ("phylogenetic_tip_neighbor_radius", 2, "tips on each side"),
            ("minimum_interval_jaccard", 0.20, "fraction"),
            ("minimum_phylogenetic_neighbor_support", 0.80, "fraction"),
            ("clade_resolution_maximum_gap", 10_000, "bp"),
            ("clade_resolution_maximum_midpoint_difference", 150_000, "bp"),
            ("clade_resolution_maximum_length_fold_difference", 3, "fold"),
            ("reporting_fraction_threshold_strictly_greater_than", 0.30, "fraction"),
            ("reporting_length_threshold_strictly_greater_than", 100_000, "bp"),
            ("nuclear_reference_length", nuclear_total_bp, "bp"),
        ],
        columns=["parameter", "value", "unit_or_definition"],
    )
    parameters.to_csv(
        output_directory / "loh_run_parameters.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
