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

"""Summarize all per-trait GCTA MLMA-LOCO outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


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

def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                values[parts[0]] = parts[1]
    return values


def read_phenotype_summary(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(row["drug"], row["mode"]): row for row in reader}


def detect_columns(fieldnames: list[str]) -> dict[str, str]:
    candidates = {
        "chr": ["Chr", "CHR", "chr"],
        "snp": ["SNP", "snp"],
        "bp": ["bp", "BP", "pos", "POS"],
        "p": ["p", "P"],
        "a1": ["A1"],
        "a2": ["A2"],
        "freq": ["Freq", "freq"],
        "b": ["b", "BETA", "beta"],
        "se": ["se", "SE"],
    }
    out: dict[str, str] = {}
    for key, cols in candidates.items():
        for col in cols:
            if col in fieldnames:
                out[key] = col
                break
    missing = {"chr", "snp", "bp", "p"} - set(out)
    if missing:
        raise ValueError(f"Missing required MLMA columns: {sorted(missing)}")
    return out


def chr_label(value: str) -> str:
    if value == "8":
        return "ChrM"
    if value == "9":
        return "ChrR"
    return f"Chr{value}"


def summarize_mlma(
    path: Path,
    bonferroni_p: float,
    drug: str,
    mode: str,
) -> tuple[dict[str, str], list[dict[str, str]], list[dict[str, str]]]:
    top_row: dict[str, str] | None = None
    n_mlma_rows = 0
    n_variants = 0
    n_bonferroni = 0
    n_p_le_1e_5 = 0
    bonferroni_hits: list[dict[str, str]] = []
    top100: list[tuple[float, dict[str, str]]] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")
        cols = detect_columns(reader.fieldnames)

        for row in reader:
            n_mlma_rows += 1
            try:
                p_value = float(row[cols["p"]])
            except (TypeError, ValueError):
                continue
            if not (0 < p_value <= 1):
                continue

            n_variants += 1
            if p_value <= bonferroni_p:
                n_bonferroni += 1
                bonferroni_hits.append(format_hit(row, cols, drug, mode, p_value))
            if p_value <= 1e-5:
                n_p_le_1e_5 += 1

            if top_row is None or p_value < float(top_row[cols["p"]]):
                top_row = row

            top100.append((p_value, row))
            if len(top100) > 100:
                top100.sort(key=lambda item: item[0])
                top100 = top100[:100]

    if top_row is None:
        raise ValueError(f"No valid rows in {path}")

    top100.sort(key=lambda item: item[0])
    top100_rows = [format_hit(row, cols, drug, mode, p) for p, row in top100]

    summary = {
        "n_mlma_rows": str(n_mlma_rows),
        "n_valid_p_variants": str(n_variants),
        "n_invalid_p_variants": str(n_mlma_rows - n_variants),
        "n_bonferroni_hits": str(n_bonferroni),
        "n_p_le_1e_5": str(n_p_le_1e_5),
        "top_chr": top_row.get(cols["chr"], ""),
        "top_chr_label": chr_label(top_row.get(cols["chr"], "")),
        "top_snp": top_row.get(cols["snp"], ""),
        "top_bp": top_row.get(cols["bp"], ""),
        "top_p": top_row.get(cols["p"], ""),
        "top_a1": top_row.get(cols.get("a1", ""), ""),
        "top_a2": top_row.get(cols.get("a2", ""), ""),
        "top_freq": top_row.get(cols.get("freq", ""), ""),
        "top_b": top_row.get(cols.get("b", ""), ""),
        "top_se": top_row.get(cols.get("se", ""), ""),
    }
    return summary, bonferroni_hits, top100_rows


def format_hit(
    row: dict[str, str],
    cols: dict[str, str],
    drug: str,
    mode: str,
    p_value: float,
) -> dict[str, str]:
    return {
        "drug": drug,
        "mode": mode,
        "chr": row.get(cols["chr"], ""),
        "chr_label": chr_label(row.get(cols["chr"], "")),
        "snp": row.get(cols["snp"], ""),
        "bp": row.get(cols["bp"], ""),
        "a1": row.get(cols.get("a1", ""), ""),
        "a2": row.get(cols.get("a2", ""), ""),
        "freq": row.get(cols.get("freq", ""), ""),
        "b": row.get(cols.get("b", ""), ""),
        "se": row.get(cols.get("se", ""), ""),
        "p": f"{p_value:.15g}",
    }


def count_report_flags(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    text = path.read_text(encoding="utf-8", errors="replace")
    return (
        text.count("Log-likelihood ratio converged."),
        text.count("component(s) constrained"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args()
    pipeline_dir = args.run_dir.resolve()
    outputs = pipeline_dir / "outputs"
    phenotype_summary = read_phenotype_summary(pipeline_dir / "00_inputs" / "phenotype_summary.tsv")

    combined_rows: list[dict[str, str]] = []
    all_bonferroni_hits: list[dict[str, str]] = []
    all_top100: list[dict[str, str]] = []

    for drug in DRUGS:
        for mode in MODES:
            result_dir = outputs / drug / mode
            mlma = result_dir / "gwas.loco.mlma"
            plot_summary = result_dir / "plots" / f"plot_summary.{drug}.{mode}.txt"
            report = result_dir / "report_gcta_mlma_loco.txt"
            if not mlma.exists() or not plot_summary.exists():
                raise SystemExit(f"Missing required output for {drug} {mode}")

            plot_values = read_key_value_file(plot_summary)
            bonferroni_p = float(plot_values["Bonferroni_P"])
            mlma_summary, bonferroni_hits, top100_rows = summarize_mlma(
                mlma,
                bonferroni_p,
                drug,
                mode,
            )
            n_converged, n_constrained = count_report_flags(report)
            pheno = phenotype_summary[(drug, mode)]

            combined_rows.append(
                {
                    "drug": drug,
                    "mode": mode,
                    "n_total_samples": pheno["n_total"],
                    "n_nonmissing_phenotype": pheno["n_nonmissing"],
                    "n_missing_phenotype": pheno["n_missing"],
                    "n_cases_R": pheno["n_cases"],
                    "n_controls_S": pheno["n_controls"],
                    "lambda_gc": plot_values.get("Lambda_GC", ""),
                    "bonferroni_p": plot_values.get("Bonferroni_P", ""),
                    "bonferroni_log10p": plot_values.get("Bonferroni_log10P", ""),
                    "n_reml_converged_messages": str(n_converged),
                    "n_variance_constraint_messages": str(n_constrained),
                    **mlma_summary,
                }
            )
            all_bonferroni_hits.extend(bonferroni_hits)
            all_top100.extend(top100_rows)

    summary_path = outputs / "gcta_mlma_loco_combined_summary.tsv"
    summary_fields = [
        "drug",
        "mode",
        "n_total_samples",
        "n_nonmissing_phenotype",
        "n_missing_phenotype",
        "n_cases_R",
        "n_controls_S",
        "n_mlma_rows",
        "n_valid_p_variants",
        "n_invalid_p_variants",
        "lambda_gc",
        "bonferroni_p",
        "bonferroni_log10p",
        "n_bonferroni_hits",
        "n_p_le_1e_5",
        "n_reml_converged_messages",
        "n_variance_constraint_messages",
        "top_chr",
        "top_chr_label",
        "top_snp",
        "top_bp",
        "top_a1",
        "top_a2",
        "top_freq",
        "top_b",
        "top_se",
        "top_p",
    ]
    write_tsv(summary_path, combined_rows, summary_fields)

    hit_fields = ["drug", "mode", "chr", "chr_label", "snp", "bp", "a1", "a2", "freq", "b", "se", "p"]
    all_bonferroni_hits.sort(key=lambda row: float(row["p"]))
    write_tsv(outputs / "gcta_mlma_loco_bonferroni_hits.tsv", all_bonferroni_hits, hit_fields)

    all_top100.sort(key=lambda row: (row["drug"], row["mode"], float(row["p"])))
    write_tsv(outputs / "gcta_mlma_loco_top100_per_trait.tsv", all_top100, hit_fields)

    print(f"Wrote {summary_path}")
    print(f"Bonferroni hits: {len(all_bonferroni_hits)}")
    print(f"Top-100 rows: {len(all_top100)}")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
