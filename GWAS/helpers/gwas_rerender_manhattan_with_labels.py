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

"""Rerender all GCTA plots with grey FungAMR LD SNPs and gene labels."""

from __future__ import annotations

import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
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
MAX_WORKERS = 4



# =============================================================================
# Functions
# =============================================================================

def render_one(maf_label: str, folder: str, drug: str, mode: str) -> dict[str, str]:
    run_dir = ROOT / folder
    outdir = run_dir / "outputs" / drug / mode
    plotdir = outdir / "plots"
    mlma = outdir / "gwas.loco.mlma"
    grey = run_dir / "outputs" / "fungamr_ld_grey_snps" / f"{drug}.ld_r2_ge0p2.window1000kb.snps.txt"
    labels = plotdir / f"manuscript_candidate_labels.{drug}.{mode}.tsv"
    if not mlma.exists():
        raise FileNotFoundError(mlma)
    if not grey.exists():
        raise FileNotFoundError(grey)
    if not labels.exists():
        raise FileNotFoundError(labels)

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
        str(grey),
        "--grey-label",
        "In LD with FungAMR mutation",
        "--label-hits",
        str(labels),
    ]
    completed = subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return {
        "maf_label": maf_label,
        "drug": drug,
        "mode": mode,
        "manhattan_png": str((plotdir / f"manhattan.{drug}.{mode}.png").relative_to(ROOT)),
        "manhattan_pdf": str((plotdir / f"manhattan.{drug}.{mode}.pdf").relative_to(ROOT)),
        "label_file": str(labels.relative_to(ROOT)),
        "grey_snps_file": str(grey.relative_to(ROOT)),
        "status": "ok",
        "log_tail": completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else "",
    }


def main() -> None:
    tasks = [
        (maf_label, folder, drug, mode)
        for maf_label, folder in RUNS
        for drug in DRUGS
        for mode in MODES
    ]
    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(render_one, *task): task for task in tasks}
        for future in as_completed(futures):
            maf_label, _folder, drug, mode = futures[future]
            row = future.result()
            rows.append(row)
            print(f"{maf_label} {drug} {mode}: rendered")

    rows.sort(key=lambda row: (row["maf_label"], row["drug"], row["mode"]))
    manifest = ROOT / "gcta_mlma_loco/gcta_mlma_loco_labeled_manhattan_rerender_manifest.tsv"
    fields = [
        "maf_label",
        "drug",
        "mode",
        "manhattan_png",
        "manhattan_pdf",
        "label_file",
        "grey_snps_file",
        "status",
        "log_tail",
    ]
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {manifest}")



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
