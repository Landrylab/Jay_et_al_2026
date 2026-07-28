#!/usr/bin/env python3
# =============================================================================
# Description
#   Helper logic for the GCTA MLMA-LOCO GWAS reproducibility bundle.
# Inputs
#   Local project metadata, PLINK files, MLMA outputs, SnpEff annotations, and GWAS tables.
# Outputs
#   Phenotype files, GWAS summaries, annotated hits, figures, or publication tables.
# =============================================================================

# =============================================================================
# Configuration
# =============================================================================

"""Build GCTA-ready MIC and R/S phenotype files for six antifungals."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DRUGS = [
    ("fluconazole", "log2_flu", "rbin_flu"),
    ("voriconazole", "log2_vor", "rbin_vor"),
    ("posaconazole", "log2_pos", "rbin_pos"),
    ("caspofungin", "log2_cas", "rbin_cas"),
    ("micafungin", "log2_mic", "rbin_mic"),
    ("anidulafungin", "log2_ani", "rbin_ani"),
]



# =============================================================================
# Functions
# =============================================================================

def is_missing(value: str | None) -> bool:
    if value is None:
        return True
    value = value.strip()
    return value == "" or value.upper() in {"NA", "NAN", "NULL"}


def clean_mic(value: str | None) -> str:
    if is_missing(value):
        return "-9"
    text = value.strip()
    try:
        float(text)
    except ValueError:
        return "-9"
    return text


def clean_binary(value: str | None) -> str:
    if is_missing(value):
        return "-9"
    text = value.strip().upper()
    if text in {"R", "1"}:
        return "1"
    if text in {"S", "0"}:
        return "0"
    return "-9"


def read_fam(fam_path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with fam_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                raise SystemExit(f"Malformed FAM line {line_number}: {line.rstrip()}")
            rows.append((parts[0], parts[1]))
    return rows


def read_metadata(metadata_path: Path) -> dict[str, dict[str, str]]:
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {"sample"} | {mic for _, mic, _ in DRUGS} | {binary for _, _, binary in DRUGS}
        missing = sorted(required - fields)
        if missing:
            raise SystemExit(f"Metadata is missing required columns: {', '.join(missing)}")

        metadata: dict[str, dict[str, str]] = {}
        for row in reader:
            sample = (row.get("sample") or "").strip()
            if sample:
                metadata[sample] = row
    return metadata


def write_pheno(
    out_path: Path,
    fam_rows: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
    column: str,
    cleaner,
    count_binary: bool = False,
) -> dict[str, int]:
    counts = {
        "n_total": 0,
        "n_nonmissing": 0,
        "n_missing": 0,
        "n_cases": 0,
        "n_controls": 0,
    }
    with out_path.open("w", encoding="utf-8", newline="") as out:
        for fid, iid in fam_rows:
            value = cleaner(metadata.get(iid, {}).get(column))
            counts["n_total"] += 1
            if value == "-9":
                counts["n_missing"] += 1
            else:
                counts["n_nonmissing"] += 1
                if count_binary and value == "1":
                    counts["n_cases"] += 1
                elif count_binary and value == "0":
                    counts["n_controls"] += 1
            out.write(f"{fid}\t{iid}\t{value}\n")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fam", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    fam_rows = read_fam(args.fam)
    metadata = read_metadata(args.metadata)

    fam_iids = {iid for _, iid in fam_rows}
    metadata_samples = set(metadata)
    missing_metadata = sorted(fam_iids - metadata_samples)
    extra_metadata = sorted(metadata_samples - fam_iids)

    summary_rows: list[dict[str, str | int]] = []
    for drug, mic_col, binary_col in DRUGS:
        mic_path = args.outdir / f"pheno.mic.{drug}.txt"
        mic_counts = write_pheno(mic_path, fam_rows, metadata, mic_col, clean_mic)
        summary_rows.append(
            {
                "mode": "mic",
                "drug": drug,
                "source_column": mic_col,
                "phenotype_file": mic_path.name,
                **mic_counts,
            }
        )

        binary_path = args.outdir / f"pheno.rs_binary.{drug}.txt"
        binary_counts = write_pheno(
            binary_path,
            fam_rows,
            metadata,
            binary_col,
            clean_binary,
            count_binary=True,
        )
        summary_rows.append(
            {
                "mode": "rs_binary",
                "drug": drug,
                "source_column": binary_col,
                "phenotype_file": binary_path.name,
                **binary_counts,
            }
        )

    summary_path = args.outdir / "phenotype_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as out:
        fieldnames = [
            "mode",
            "drug",
            "source_column",
            "phenotype_file",
            "n_total",
            "n_nonmissing",
            "n_missing",
            "n_cases",
            "n_controls",
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    drugs_path = args.outdir / "drugs.txt"
    with drugs_path.open("w", encoding="utf-8") as out:
        for drug, _, _ in DRUGS:
            out.write(f"{drug}\n")

    missing_path = args.outdir / "fam_samples_missing_metadata.txt"
    with missing_path.open("w", encoding="utf-8") as out:
        for iid in missing_metadata:
            out.write(f"{iid}\n")

    extra_path = args.outdir / "metadata_samples_not_in_fam.txt"
    with extra_path.open("w", encoding="utf-8") as out:
        for sample in extra_metadata:
            out.write(f"{sample}\n")

    print(f"Wrote phenotype files under {args.outdir}")
    print(f"FAM samples: {len(fam_rows)}")
    print(f"Metadata samples: {len(metadata)}")
    print(f"FAM samples missing metadata: {len(missing_metadata)}")
    print(f"Metadata samples not in FAM: {len(extra_metadata)}")
    print(f"Summary: {summary_path}")



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
