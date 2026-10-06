#!/usr/bin/env python3
"""Count callable and heterozygous SNP positions in fixed genomic windows."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import numpy as np
import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count callable and heterozygous diploid SNP positions for a fixed "
            "sample set. Records at the same genomic position are counted once."
        )
    )
    parser.add_argument("--vcf", required=True, help="Indexed cohort VCF/BCF")
    parser.add_argument("--fai", required=True, help="Reference FASTA index")
    parser.add_argument(
        "--samples",
        required=True,
        help="One-column file containing the final sample identifiers",
    )
    parser.add_argument("--out", required=True, help="Output .tsv or .tsv.gz")
    parser.add_argument("--window-bp", type=int, default=10_000)
    parser.add_argument(
        "--exclude-contig-regex",
        default=r"(?i)chrM(?:_|$)",
        help="Regular expression for contigs excluded from the nuclear analysis",
    )
    parser.add_argument("--expected-samples", type=int, help="Optional sample-count assertion")
    return parser.parse_args()


def load_samples(path: str) -> list[str]:
    samples = [
        line.strip().split("\t")[0]
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if samples and samples[0].lower() in {"sample", "sample_id", "isolate"}:
        samples = samples[1:]
    if not samples:
        raise ValueError("The sample file is empty.")
    if len(samples) != len(set(samples)):
        raise ValueError("The sample file contains duplicate identifiers.")
    return samples


def load_contigs(
    fai_path: str, exclude_pattern: str, window_bp: int
) -> tuple[list[str], dict[str, int], dict[str, int]]:
    excluded = re.compile(exclude_pattern)
    ordered: list[str] = []
    lengths: dict[str, int] = {}
    windows: dict[str, int] = {}
    with open(fai_path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            contig, length_text = fields[:2]
            if excluded.search(contig):
                continue
            length = int(length_text)
            ordered.append(contig)
            lengths[contig] = length
            windows[contig] = math.ceil(length / window_bp)
    if not ordered:
        raise ValueError("No nuclear contigs remain after filtering the FASTA index.")
    return ordered, lengths, windows


def is_snp_record(record: pysam.VariantRecord) -> bool:
    if len(record.ref) != 1 or not record.alts:
        return False
    return all(len(alt) == 1 and alt not in {"*", "."} for alt in record.alts)


def called_and_heterozygous(sample_call: pysam.libcbcf.VariantRecordSample) -> tuple[bool, bool]:
    genotype = sample_call.get("GT")
    if genotype is None or len(genotype) != 2 or any(allele is None for allele in genotype):
        return False, False
    return True, genotype[0] != genotype[1]


def open_text_output(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def main() -> None:
    args = parse_args()
    if args.window_bp <= 0:
        raise ValueError("--window-bp must be positive.")

    samples = load_samples(args.samples)
    if args.expected_samples and len(samples) != args.expected_samples:
        raise ValueError(
            f"Expected {args.expected_samples} samples but found {len(samples)}."
        )

    contigs, contig_lengths, n_windows = load_contigs(
        args.fai, args.exclude_contig_regex, args.window_bp
    )
    vcf = pysam.VariantFile(args.vcf)
    missing_samples = [sample for sample in samples if sample not in vcf.header.samples]
    if missing_samples:
        preview = ", ".join(missing_samples[:10])
        raise ValueError(f"Samples absent from the VCF header: {preview}")

    n_call = {
        contig: np.zeros((len(samples), n_windows[contig]), dtype=np.int32)
        for contig in contigs
    }
    n_het = {
        contig: np.zeros((len(samples), n_windows[contig]), dtype=np.int32)
        for contig in contigs
    }

    for contig in contigs:
        current_position: int | None = None
        position_called = np.zeros(len(samples), dtype=bool)
        position_het = np.zeros(len(samples), dtype=bool)

        def flush_position() -> None:
            if current_position is None:
                return
            window_index = (current_position - 1) // args.window_bp
            n_call[contig][:, window_index] += position_called
            n_het[contig][:, window_index] += position_het

        for record in vcf.fetch(contig):
            if not is_snp_record(record):
                continue
            if current_position is None or record.pos != current_position:
                flush_position()
                current_position = record.pos
                position_called.fill(False)
                position_het.fill(False)
            for index, sample in enumerate(samples):
                called, heterozygous = called_and_heterozygous(record.samples[sample])
                position_called[index] |= called
                position_het[index] |= heterozygous
        flush_position()

    with open_text_output(args.out) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "contig",
                "window_start0",
                "window_end0",
                "N_call",
                "N_het",
                "N_hom",
                "HET_FRAC",
            ]
        )
        for sample_index, sample in enumerate(samples):
            for contig in contigs:
                for window_index in range(n_windows[contig]):
                    start0 = window_index * args.window_bp
                    end0 = start0 + args.window_bp
                    callable_count = int(n_call[contig][sample_index, window_index])
                    heterozygous_count = int(n_het[contig][sample_index, window_index])
                    homozygous_count = callable_count - heterozygous_count
                    fraction = (
                        heterozygous_count / callable_count
                        if callable_count
                        else float("nan")
                    )
                    writer.writerow(
                        [
                            sample,
                            contig,
                            start0,
                            end0,
                            callable_count,
                            heterozygous_count,
                            homozygous_count,
                            fraction,
                        ]
                    )


if __name__ == "__main__":
    main()
