#!/usr/bin/env python3
# =============================================================================
# Description
#   Annotate Bonferroni-significant GCTA MLMA-LOCO hits and render labelled
#   Manhattan/LD-context GWAS outputs.
#
# Inputs
#   gcta_mlma_loco/gcta_mlma_loco_maf005_six_antifungals/
#   gcta_mlma_loco/gcta_mlma_loco_maf0005_six_antifungals/
#   SnpEff-annotated VCF and PLINK resources referenced by helper defaults
#
# Outputs
#   gcta_mlma_loco/gcta_mlma_loco_significant_snp_annotation_csvs/
#   Label tables, annotated significant-hit CSVs, and rendered GWAS plots.
# =============================================================================

from __future__ import annotations

# =============================================================================
# Configuration
# =============================================================================

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HELPERS = Path(__file__).resolve().parent / "helpers"
PYTHON = sys.executable

HELPER_SEQUENCE = [
    "gwas_annotate_bonferroni_hits.py",
    "gwas_create_manhattan_label_tables.py",
    "gwas_grey_fungamr_ld_manhattan_plots.py",
    "gwas_rerender_manhattan_with_labels.py",
]


# =============================================================================
# Functions
# =============================================================================

def run_helper(script_name: str) -> None:
    cmd = [PYTHON, str(HELPERS / script_name)]
    os.chdir(ROOT)
    subprocess.check_call(cmd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate and render GCTA MLMA-LOCO significant-hit outputs.")
    parser.add_argument("--skip-ld-grey-plots", action="store_true")
    parser.add_argument("--skip-rerender", action="store_true")
    return parser.parse_args()


# =============================================================================
# Main block
# =============================================================================

def main() -> None:
    args = parse_args()
    for script_name in HELPER_SEQUENCE:
        if args.skip_ld_grey_plots and script_name == "gwas_grey_fungamr_ld_manhattan_plots.py":
            continue
        if args.skip_rerender and script_name == "gwas_rerender_manhattan_with_labels.py":
            continue
        run_helper(script_name)


if __name__ == "__main__":
    main()

