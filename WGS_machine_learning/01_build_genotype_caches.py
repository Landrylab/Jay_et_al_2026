#!/usr/bin/env python3
"""Build compact genotype caches from the two VCF inputs used by the ML workflow."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="Cohort metadata TSV.")
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="rbin_flu")
    parser.add_argument("--group-col", default="cluster")
    parser.add_argument(
        "--snp-vcf",
        required=True,
        help="Indexed biallelic SNP VCF used for homozygous-alternate ML features.",
    )
    parser.add_argument(
        "--pc-vcf",
        required=True,
        help="Indexed LD-pruned biallelic SNP VCF used for fold-derived PCs.",
    )
    parser.add_argument("--output", required=True, help="Output cache directory.")
    parser.add_argument("--bcftools", default="bcftools")
    return parser.parse_args()


def run_text(command: list[str]) -> str:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return result.stdout


def vcf_samples(bcftools: str, vcf: Path) -> list[str]:
    return [line for line in run_text([bcftools, "query", "-l", str(vcf)]).splitlines() if line]


def vcf_record_count(bcftools: str, vcf: Path) -> int:
    value = run_text([bcftools, "index", "-n", str(vcf)]).strip()
    return int(value)


def parse_gt(value: str, encoding: str) -> int:
    genotype = value.split(":", 1)[0].replace("|", "/")
    alleles = genotype.split("/")
    if len(alleles) != 2 or any(allele == "." for allele in alleles):
        return -1
    try:
        numbers = [int(allele) for allele in alleles]
    except ValueError:
        return -1
    if any(number not in (0, 1) for number in numbers):
        raise ValueError("The cache builder requires decomposed biallelic records.")
    if encoding == "hom_alt":
        return int(numbers[0] == 1 and numbers[1] == 1)
    return int(numbers[0] + numbers[1])


def build_cache(
    bcftools: str,
    vcf: Path,
    selected_sample_indices: np.ndarray,
    selected_samples: list[str],
    output_matrix: Path,
    output_variants: Path,
    encoding: str,
) -> int:
    allocated_records = vcf_record_count(bcftools, vcf)
    matrix = np.memmap(
        output_matrix,
        dtype=np.int8,
        mode="w+",
        shape=(allocated_records, len(selected_samples)),
    )
    query_format = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    process = subprocess.Popen(
        [bcftools, "query", "-f", query_format, str(vcf)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    if process.stdout is None:
        raise RuntimeError("Unable to read bcftools output.")
    with output_variants.open("w") as variant_handle:
        variant_handle.write("variant_idx\tchrom\tpos\tref\talt\tvariant_id\n")
        record_index = 0
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, ref, alt = fields[:4]
            if "," in alt:
                raise ValueError(f"Multiallelic record encountered: {chrom}:{pos}:{ref}:{alt}")
            all_genotypes = fields[4:]
            selected = [all_genotypes[index] for index in selected_sample_indices]
            matrix[record_index, :] = [parse_gt(value, encoding) for value in selected]
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            variant_handle.write(
                f"{record_index}\t{chrom}\t{pos}\t{ref}\t{alt}\t{variant_id}\n"
            )
            record_index += 1
            if record_index % 100000 == 0:
                print(f"{vcf.name}: processed {record_index:,} records", flush=True)
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bcftools query failed for {vcf}:\n{stderr}")
    matrix.flush()
    if record_index != allocated_records:
        raise ValueError(
            f"Indexed record count ({allocated_records}) differs from parsed count ({record_index}) for {vcf}."
        )
    return record_index


def main() -> None:
    args = parse_args()
    metadata_path = Path(args.metadata)
    snp_vcf = Path(args.snp_vcf)
    pc_vcf = Path(args.pc_vcf)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    required = {args.sample_col, args.label_col, args.group_col}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")
    metadata[args.sample_col] = metadata[args.sample_col].astype(str)
    metadata = metadata.loc[
        metadata[args.label_col].astype(str).isin(["R", "S"]) & metadata[args.group_col].notna()
    ].drop_duplicates(args.sample_col)
    labelled = set(metadata[args.sample_col])

    snp_header_samples = vcf_samples(args.bcftools, snp_vcf)
    pc_header_samples = set(vcf_samples(args.bcftools, pc_vcf))
    selected_samples = [sample for sample in snp_header_samples if sample in labelled]
    if not selected_samples:
        raise ValueError("No labelled metadata samples were found in the SNP VCF.")
    missing_pc = sorted(set(selected_samples).difference(pc_header_samples))
    if missing_pc:
        raise ValueError(f"{len(missing_pc)} labelled samples are absent from the PC-marker VCF.")
    missing_metadata = sorted(labelled.difference(selected_samples))
    if missing_metadata:
        raise ValueError(f"{len(missing_metadata)} labelled metadata samples are absent from the SNP VCF.")

    snp_indices = np.asarray([snp_header_samples.index(sample) for sample in selected_samples], dtype=int)
    pc_header_list = vcf_samples(args.bcftools, pc_vcf)
    pc_indices = np.asarray([pc_header_list.index(sample) for sample in selected_samples], dtype=int)
    pd.DataFrame({"sample": selected_samples}).to_csv(output / "sample_ids.tsv", sep="\t", index=False)

    n_snp = build_cache(
        args.bcftools,
        snp_vcf,
        snp_indices,
        selected_samples,
        output / "snp_genotypes.int8.mmap",
        output / "snp_variants.tsv",
        "hom_alt",
    )
    n_pc = build_cache(
        args.bcftools,
        pc_vcf,
        pc_indices,
        selected_samples,
        output / "pc_genotypes.int8.mmap",
        output / "pc_variants.tsv",
        "additive",
    )
    manifest = {
        "metadata": str(metadata_path),
        "snp_vcf": str(snp_vcf),
        "pc_vcf": str(pc_vcf),
        "n_samples": len(selected_samples),
        "n_resistant": int((metadata.set_index(args.sample_col).loc[selected_samples, args.label_col] == "R").sum()),
        "n_susceptible": int((metadata.set_index(args.sample_col).loc[selected_samples, args.label_col] == "S").sum()),
        "n_snp_records": n_snp,
        "n_pc_markers": n_pc,
        "snp_encoding": "homozygous alternate=1; homozygous reference or heterozygous=0; missing=-1",
        "pc_encoding": "alternate-allele dosage 0/1/2; missing=-1",
        "sample_col": args.sample_col,
        "label_col": args.label_col,
        "group_col": args.group_col,
    }
    (output / "cache_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Cache written to {output}", flush=True)


if __name__ == "__main__":
    main()
