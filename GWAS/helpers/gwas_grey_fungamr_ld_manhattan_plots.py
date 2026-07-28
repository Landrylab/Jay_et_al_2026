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

"""Regenerate GCTA Manhattan plots with high-confidence FungAMR LD SNPs greyed."""

from __future__ import annotations

import csv
import importlib.util
import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
ANNOTATION_SCRIPT = Path(__file__).resolve().parent / "gwas_annotate_bonferroni_hits.py"
BASE_PLINK = (
    ROOT
    / "gcta_mlma_loco"
    / "plink_from_joint_vcf"
    / "joint_keep1153"
)
PLINK_BIN = os.environ.get("PLINK_BIN", "plink")
R2_THRESHOLD = 0.2
WINDOW_KB = 1000

RUNS = [
    ("maf0p05", "gcta_mlma_loco/gcta_mlma_loco_maf005_six_antifungals"),
    ("maf0p005", "gcta_mlma_loco/gcta_mlma_loco_maf0005_six_antifungals"),
]
DRUGS = [
    "fluconazole",
    "voriconazole",
    "posaconazole",
    "caspofungin",
    "micafungin",
    "anidulafungin",
]
MODES = ["mic", "rs_binary"]



# =============================================================================
# Functions
# =============================================================================

def load_annotation_module():
    spec = importlib.util.spec_from_file_location("gcta_hit_annotation", ANNOTATION_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {ANNOTATION_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_bim_ids(prefix: Path) -> set[str]:
    ids: set[str] = set()
    with Path(str(prefix) + ".bim").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                ids.add(fields[1])
    return ids


def read_bim_records(prefix: Path) -> dict[str, tuple[str, int]]:
    records: dict[str, tuple[str, int]] = {}
    with Path(str(prefix) + ".bim").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                bp = int(float(fields[3]))
            except ValueError:
                continue
            records[fields[1]] = (fields[0], bp)
    return records


def filtered_prefix(run_dir: Path) -> Path:
    bims = sorted(run_dir.glob("01_*_filtered_snps/*.bim"))
    if len(bims) != 1:
        raise RuntimeError(f"Expected one filtered BIM under {run_dir}, found {len(bims)}")
    return bims[0].with_suffix("")


def write_lines(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def update_best_detail(
    details: dict[tuple[str, str, str, str], dict[str, str]],
    snp_id: str,
    rec: dict[str, str],
    r2: float,
) -> None:
    key = (rec["drug"], snp_id, rec["fungamr_record"], rec["locus_id"])
    old = details.get(key)
    if old is not None and float(old["r2"]) >= r2:
        return
    details[key] = {
        "drug": rec["drug"],
        "snp": snp_id,
        "fungamr_record": rec["fungamr_record"],
        "fungamr_gene": rec["fungamr_gene"],
        "fungamr_substitution": rec["fungamr_substitution"],
        "fungamr_score": rec["fungamr_score"],
        "fungamr_locus": rec["locus_id"],
        "r2": f"{r2:.6g}",
    }


def compute_ld_details(
    filtered_records: dict[str, tuple[str, int]],
    base_records: dict[str, tuple[str, int]],
    records_by_locus: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str, str, str], dict[str, str]]:
    details: dict[tuple[str, str, str, str], dict[str, str]] = {}
    filtered_ids = set(filtered_records)
    fungamr_loci = set(records_by_locus)
    available_fungamr_loci = sorted(fungamr_loci & set(base_records))
    if not available_fungamr_loci:
        return details

    for locus in available_fungamr_loci:
        if locus in filtered_ids:
            for rec in records_by_locus[locus]:
                update_best_detail(details, locus, rec, 1.0)

    window_bp = WINDOW_KB * 1000
    target_by_chr: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for locus in available_fungamr_loci:
        chrom, bp = base_records[locus]
        target_by_chr[chrom].append((locus, bp))

    nearby_filtered_ids: set[str] = set()
    for snp_id, (chrom, bp) in filtered_records.items():
        for _locus, target_bp in target_by_chr.get(chrom, []):
            if abs(bp - target_bp) <= window_bp:
                nearby_filtered_ids.add(snp_id)
                break

    # Candidate extraction is limited to the LD search window around high-confidence FungAMR loci.
    available_extract_ids = sorted((nearby_filtered_ids | set(available_fungamr_loci)) & set(base_records))
    with tempfile.TemporaryDirectory(prefix="fungamr_ld_manhattan_", dir="/private/tmp") as td:
        tmpdir = Path(td)
        extract_path = tmpdir / "extract_ids.txt"
        fungamr_path = tmpdir / "fungamr_loci.txt"
        write_lines(extract_path, available_extract_ids)
        write_lines(fungamr_path, available_fungamr_loci)

        out_prefix = tmpdir / "fungamr_ld"
        subprocess.run(
            [
                str(PLINK_BIN),
                "--bfile",
                str(BASE_PLINK),
                "--allow-extra-chr",
                "--keep-allele-order",
                "--extract",
                str(extract_path),
                "--r2",
                "--ld-snp-list",
                str(fungamr_path),
                "--ld-window",
                "999999",
                "--ld-window-kb",
                str(WINDOW_KB),
                "--ld-window-r2",
                "0",
                "--out",
                str(out_prefix),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

        ld_path = Path(str(out_prefix) + ".ld")
        if not ld_path.exists():
            return details

        with ld_path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
            for row in reader:
                snp_a = row.get("SNP_A", "")
                snp_b = row.get("SNP_B", "")
                try:
                    r2 = float(row.get("R2", "nan"))
                except ValueError:
                    continue
                if not math.isfinite(r2) or r2 < R2_THRESHOLD:
                    continue
                for fung_locus, plotted_snp in ((snp_a, snp_b), (snp_b, snp_a)):
                    if fung_locus not in records_by_locus or plotted_snp not in filtered_ids:
                        continue
                    for rec in records_by_locus[fung_locus]:
                        update_best_detail(details, plotted_snp, rec, r2)

    return details


def rerender_trait(run_dir: Path, drug: str, mode: str, grey_file: Path) -> None:
    outdir = run_dir / "outputs" / drug / mode
    mlma = outdir / "gwas.loco.mlma"
    plotdir = outdir / "plots"
    label_file = plotdir / f"manuscript_candidate_labels.{drug}.{mode}.tsv"
    if not mlma.exists():
        raise FileNotFoundError(mlma)
    cmd = [
        "Rscript",
        "--vanilla",
        str(run_dir / "scripts" / "plot_mlma_results.R"),
        "--mlma",
        str(mlma),
        "--outdir",
        str(plotdir),
        "--trait",
        drug,
        "--mode",
        mode,
        "--grey-snps",
        str(grey_file),
        "--grey-label",
        "In LD with FungAMR mutation",
    ]
    if label_file.exists():
        cmd.extend(["--label-hits", str(label_file)])
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def main() -> None:
    if shutil.which(PLINK_BIN) is None and not Path(PLINK_BIN).exists():
        raise FileNotFoundError(PLINK_BIN)
    if not Path(str(BASE_PLINK) + ".bim").exists():
        raise FileNotFoundError(Path(str(BASE_PLINK) + ".bim"))

    annotation = load_annotation_module()
    fungamr_records, _by_locus_unused = annotation.load_fungamr_records()
    records_by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for rec in fungamr_records:
        records_by_locus[rec["locus_id"]].append(rec)

    base_records = read_bim_records(BASE_PLINK)
    summary_rows: list[dict[str, str]] = []
    detail_fields = [
        "drug",
        "snp",
        "fungamr_record",
        "fungamr_gene",
        "fungamr_substitution",
        "fungamr_score",
        "fungamr_locus",
        "r2",
    ]
    r2_label = f"r2_ge{str(R2_THRESHOLD).replace('.', 'p')}"

    for maf_label, folder in RUNS:
        run_dir = ROOT / folder
        prefix = filtered_prefix(run_dir)
        print(f"[{maf_label}] reading filtered SNPs from {prefix}.bim", flush=True)
        filtered_records = read_bim_records(prefix)
        details = compute_ld_details(filtered_records, base_records, records_by_locus)

        out_root = run_dir / "outputs" / "fungamr_ld_grey_snps"
        detail_by_drug: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in details.values():
            detail_by_drug[row["drug"]].append(row)

        for drug in DRUGS:
            drug_details = sorted(
                detail_by_drug.get(drug, []),
                key=lambda row: (row["snp"], row["fungamr_record"], row["fungamr_locus"]),
            )
            grey_snps = sorted({row["snp"] for row in drug_details})
            grey_file = out_root / f"{drug}.ld_{r2_label}.window{WINDOW_KB}kb.snps.txt"
            detail_file = out_root / f"{drug}.ld_{r2_label}.window{WINDOW_KB}kb.details.tsv"
            write_lines(grey_file, grey_snps)
            write_tsv(detail_file, drug_details, detail_fields)

            for mode in MODES:
                rerender_trait(run_dir, drug, mode, grey_file)

            summary_rows.append(
                {
                    "maf_label": maf_label,
                    "drug": drug,
                    "n_grey_snps": str(len(grey_snps)),
                    "n_ld_detail_rows": str(len(drug_details)),
                    "grey_snps_file": str(grey_file.relative_to(ROOT)),
                    "ld_details_file": str(detail_file.relative_to(ROOT)),
                }
            )
            print(f"  {drug}: {len(grey_snps)} grey SNPs", flush=True)

    write_tsv(
        ROOT / "gcta_mlma_loco/gcta_mlma_loco_fungamr_ld_grey_plot_manifest.tsv",
        summary_rows,
        [
            "maf_label",
            "drug",
            "n_grey_snps",
            "n_ld_detail_rows",
            "grey_snps_file",
            "ld_details_file",
        ],
    )
    print("Regenerated Manhattan plots with FungAMR LD SNPs greyed.", flush=True)



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
