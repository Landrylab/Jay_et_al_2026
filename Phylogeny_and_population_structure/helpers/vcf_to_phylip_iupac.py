#!/usr/bin/env python3
# =============================================================================
# Description
#   Convert diploid biallelic SNP genotypes to a sequential PHYLIP alignment.
#   Heterozygous genotypes use IUPAC ambiguity codes and missing calls use N.
#
# Inputs
#   Indexed VCF or BCF containing one biallelic SNP record per position.
#
# Outputs
#   Sequential PHYLIP alignment.
# =============================================================================

from __future__ import annotations

import argparse
from pathlib import Path

import pysam


IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def encode_call(record: pysam.VariantRecord, sample: str) -> str:
    genotype = record.samples[sample].get("GT")
    if genotype is None or len(genotype) != 2 or any(allele is None for allele in genotype):
        return "N"
    alleles = [record.alleles[index].upper() for index in genotype]
    if any(len(base) != 1 or base not in "ACGT" for base in alleles):
        return "N"
    if alleles[0] == alleles[1]:
        return alleles[0]
    return IUPAC.get(frozenset(alleles), "N")


def main() -> None:
    args = parse_args()
    alignment: dict[str, list[str]] = {}
    site_count = 0

    with pysam.VariantFile(args.vcf) as variants:
        samples = list(variants.header.samples)
        alignment = {sample: [] for sample in samples}
        for record in variants:
            if len(record.ref) != 1 or len(record.alts or ()) != 1 or len(record.alts[0]) != 1:
                raise ValueError(f"Non-biallelic SNP record encountered at {record.contig}:{record.pos}")
            for sample in samples:
                alignment[sample].append(encode_call(record, sample))
            site_count += 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(alignment)} {site_count}\n")
        for sample, sequence in alignment.items():
            handle.write(f"{sample} {''.join(sequence)}\n")


if __name__ == "__main__":
    main()
