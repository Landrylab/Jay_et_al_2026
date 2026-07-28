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

"""Build publication-ready tables for GCTA MLMA-LOCO GWAS outputs.

Only the manuscript-selected MAF thresholds are exported: 0.05 and 0.005.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
GWAS_ROOT = ROOT / "gcta_mlma_loco"
SUPP_DIR = GWAS_ROOT / "supplementary_data_tables"
EXT_DIR = GWAS_ROOT / "extended_tables"

RUNS = [
    {
        "maf_threshold": "0.05",
        "maf_label": "maf0p05",
        "folder": GWAS_ROOT / "gcta_mlma_loco_maf005_six_antifungals",
        "label": "MAF 0.05",
    },
    {
        "maf_threshold": "0.005",
        "maf_label": "maf0p005",
        "folder": GWAS_ROOT / "gcta_mlma_loco_maf0005_six_antifungals",
        "label": "MAF 0.005",
    },
]

DRUGS = ["fluconazole", "voriconazole", "posaconazole", "caspofungin", "micafungin", "anidulafungin"]
MODES = ["mic", "rs_binary"]
ALLOWED_MAFS = {run["maf_threshold"] for run in RUNS}
ALLOWED_LABELS = {run["maf_label"] for run in RUNS}



# =============================================================================
# Functions
# =============================================================================

def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", **kwargs)


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)


def add_neglog10p(df: pd.DataFrame, p_col: str = "p") -> pd.DataFrame:
    if p_col not in df.columns:
        return df
    out = df.copy()
    p = pd.to_numeric(out[p_col], errors="coerce")
    out["neg_log10_p"] = np.where(p > 0, -np.log10(p), np.nan)
    return out


def normalize_maf_threshold(value: object) -> str:
    value_str = str(value)
    if value_str in ALLOWED_MAFS:
        return value_str
    try:
        numeric = float(value_str)
    except ValueError:
        return value_str
    if np.isclose(numeric, 0.05):
        return "0.05"
    if np.isclose(numeric, 0.005):
        return "0.005"
    return value_str


def concat_with_run_metadata(filename: str, reader=read_table) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for run in RUNS:
        path = run["folder"] / filename
        df = reader(path)
        df.insert(0, "maf_label", run["maf_label"])
        df.insert(0, "maf_threshold", run["maf_threshold"])
        df.insert(2, "analysis_folder", rel(run["folder"]))
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def build_run_inventory(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in summary.to_dict("records"):
        run = next(item for item in RUNS if item["maf_threshold"] == str(row["maf_threshold"]))
        base = run["folder"] / "outputs" / row["drug"] / row["mode"]
        rows.append(
            {
                "maf_threshold": row["maf_threshold"],
                "maf_label": row["maf_label"],
                "drug": row["drug"],
                "mode": row["mode"],
                "model": "GCTA --mlma-loco",
                "analysis_folder": rel(run["folder"]),
                "phenotype_file": rel(run["folder"] / "00_inputs" / f"pheno.{row['mode']}.{row['drug']}.txt"),
                "plink_prefix_folder": rel(next(run["folder"].glob("01_*_filtered_snps"))),
                "ld_pruned_grm_folder": rel(run["folder"] / "02_ld_pruned_grm"),
                "qcovar_file": rel(run["folder"] / "02_ld_pruned_grm" / "covar_pca.pc1_pc10.txt"),
                "mlma_file": rel(base / "gwas.loco.mlma"),
                "report_file": rel(base / "report_gcta_mlma_loco.txt"),
                "manhattan_png": rel(base / "plots" / f"manhattan.{row['drug']}.{row['mode']}.png"),
                "qq_png": rel(base / "plots" / f"qq.{row['drug']}.{row['mode']}.png"),
                "n_total_samples": row.get("n_total_samples"),
                "n_nonmissing_phenotype": row.get("n_nonmissing_phenotype"),
                "n_cases_R": row.get("n_cases_R"),
                "n_controls_S": row.get("n_controls_S"),
                "n_valid_p_variants": row.get("n_valid_p_variants"),
                "lambda_gc": row.get("lambda_gc"),
                "bonferroni_p": row.get("bonferroni_p"),
                "n_bonferroni_hits": row.get("n_bonferroni_hits"),
                "top_snp": row.get("top_snp"),
                "top_chr_label": row.get("top_chr_label"),
                "top_bp": row.get("top_bp"),
                "top_p": row.get("top_p"),
            }
        )
    return pd.DataFrame(rows)


def build_key_signal_summary(summary: pd.DataFrame, annotated: pd.DataFrame) -> pd.DataFrame:
    ann = annotated.copy()
    ann["p_numeric"] = pd.to_numeric(ann.get("p"), errors="coerce")
    top_ann = (
        ann.sort_values(["maf_threshold", "drug", "mode", "p_numeric"], na_position="last")
        .groupby(["maf_threshold", "drug", "mode"], as_index=False)
        .first()
    )
    cols = [
        "maf_threshold",
        "drug",
        "mode",
        "snp",
        "gene_common_name",
        "gene_standard_name",
        "effect",
        "impact",
        "aa_substitution",
        "in_ld_with_fungamr_score1_or_2",
        "fungamr_ld_genes",
        "fungamr_ld_substitutions",
        "fungamr_ld_max_r2",
    ]
    cols = [c for c in cols if c in top_ann.columns]
    top_ann = top_ann[cols].rename(columns={c: f"top_bonferroni_{c}" for c in cols if c not in {"maf_threshold", "drug", "mode"}})
    return summary.merge(top_ann, on=["maf_threshold", "drug", "mode"], how="left")


def build_maf_comparison(summary: pd.DataFrame) -> pd.DataFrame:
    keep = [
        "drug",
        "mode",
        "n_valid_p_variants",
        "lambda_gc",
        "bonferroni_p",
        "n_bonferroni_hits",
        "n_p_le_1e_5",
        "top_snp",
        "top_chr_label",
        "top_bp",
        "top_p",
    ]
    wide = summary[["maf_threshold", *keep]].pivot(index=["drug", "mode"], columns="maf_threshold")
    wide.columns = [f"{metric}_maf{maf.replace('.', 'p')}" for metric, maf in wide.columns]
    wide = wide.reset_index()
    if {"n_valid_p_variants_maf0p005", "n_valid_p_variants_maf0p05"}.issubset(wide.columns):
        wide["delta_valid_variants_0p005_minus_0p05"] = (
            pd.to_numeric(wide["n_valid_p_variants_maf0p005"], errors="coerce")
            - pd.to_numeric(wide["n_valid_p_variants_maf0p05"], errors="coerce")
        )
    if {"n_bonferroni_hits_maf0p005", "n_bonferroni_hits_maf0p05"}.issubset(wide.columns):
        wide["delta_bonferroni_hits_0p005_minus_0p05"] = (
            pd.to_numeric(wide["n_bonferroni_hits_maf0p005"], errors="coerce")
            - pd.to_numeric(wide["n_bonferroni_hits_maf0p05"], errors="coerce")
        )
    return wide


def build_annotation_breakdown(annotated: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["maf_threshold", "maf_label", "drug", "mode", "effect", "impact", "is_in_coding_sequence"]
    for col in group_cols:
        if col not in annotated.columns:
            annotated[col] = ""
    return (
        annotated.groupby(group_cols, dropna=False)
        .agg(
            n_bonferroni_hits=("snp", "count"),
            n_unique_snps=("snp", "nunique"),
            n_unique_genes=("gene_common_name", lambda x: x.replace("", np.nan).dropna().nunique()),
            min_p=("p", lambda x: pd.to_numeric(x, errors="coerce").min()),
        )
        .reset_index()
        .sort_values(["maf_threshold", "drug", "mode", "n_bonferroni_hits"], ascending=[True, True, True, False])
    )


def build_gene_summary(annotated: pd.DataFrame) -> pd.DataFrame:
    ann = annotated.copy()
    ann["gene_label"] = ann["gene_common_name"].fillna("")
    ann.loc[ann["gene_label"].eq(""), "gene_label"] = ann["gene_standard_name"].fillna("")
    ann.loc[ann["gene_label"].eq(""), "gene_label"] = ann["snpeff_gene_name"].fillna("")
    ann.loc[ann["gene_label"].eq(""), "gene_label"] = ann["gene_context"].fillna("")
    ann.loc[ann["gene_label"].eq(""), "gene_label"] = "unassigned_or_intergenic"
    ann["p_numeric"] = pd.to_numeric(ann["p"], errors="coerce")
    out = (
        ann.groupby(["gene_label"], dropna=False)
        .agg(
            n_bonferroni_hit_rows=("snp", "count"),
            n_unique_snps=("snp", "nunique"),
            n_drug_mode_tests=("mode", lambda x: len(set(zip(ann.loc[x.index, "drug"], x)))),
            maf_thresholds=("maf_threshold", lambda x: ";".join(sorted(set(map(str, x))))),
            drugs=("drug", lambda x: ";".join(sorted(set(map(str, x))))),
            modes=("mode", lambda x: ";".join(sorted(set(map(str, x))))),
            effects=("effect", lambda x: ";".join(sorted(set(map(str, x.dropna()))))),
            impacts=("impact", lambda x: ";".join(sorted(set(map(str, x.dropna()))))),
            min_p=("p_numeric", "min"),
            example_snp=("snp", "first"),
            example_substitution=("aa_substitution", "first"),
        )
        .reset_index()
        .sort_values(["n_drug_mode_tests", "n_unique_snps", "min_p"], ascending=[False, False, True])
    )
    return out


def build_fungamr_ld_summary(annotated: pd.DataFrame) -> pd.DataFrame:
    ann = annotated.copy()
    if "in_ld_with_fungamr_score1_or_2" not in ann.columns:
        return pd.DataFrame()
    ann["ld_high_conf"] = ann["in_ld_with_fungamr_score1_or_2"].fillna("").str.lower().eq("yes")
    return (
        ann.groupby(["maf_threshold", "maf_label", "drug", "mode"], dropna=False)
        .agg(
            n_bonferroni_hits=("snp", "count"),
            n_hits_in_ld_with_high_conf_fungamr=("ld_high_conf", "sum"),
            fungamr_ld_genes=("fungamr_ld_genes", lambda x: ";".join(sorted(set(v for v in x.dropna().astype(str) if v)))),
            fungamr_ld_substitutions=("fungamr_ld_substitutions", lambda x: ";".join(sorted(set(v for v in x.dropna().astype(str) if v)))),
            max_fungamr_ld_r2=("fungamr_ld_max_r2", lambda x: pd.to_numeric(x, errors="coerce").max()),
        )
        .reset_index()
        .sort_values(["maf_threshold", "drug", "mode"])
    )


def build_qc(summary: pd.DataFrame, top100: pd.DataFrame, bonf_raw: pd.DataFrame, annotated: pd.DataFrame, inventory: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def add(status: str, severity: str, check: str, evidence: str, recommendation: str = "No action required.") -> None:
        rows.append(
            {
                "status": status,
                "severity": severity,
                "check": check,
                "evidence": evidence,
                "recommendation": recommendation,
            }
        )

    expected_rows = len(RUNS) * len(DRUGS) * len(MODES)
    add(
        "PASS" if len(summary) == expected_rows else "FAIL",
        "High",
        "Trait-level summary coverage",
        f"Observed {len(summary)} rows; expected {expected_rows} rows for 2 MAF thresholds x 6 drugs x 2 phenotypes.",
        "Investigate missing run summaries before manuscript use." if len(summary) != expected_rows else "No action required.",
    )
    bad_mafs = set(map(str, summary["maf_threshold"])) - ALLOWED_MAFS
    add(
        "PASS" if not bad_mafs else "FAIL",
        "High",
        "Requested MAF thresholds only",
        f"Unexpected MAF thresholds in summary: {sorted(bad_mafs)}.",
        "Remove non-requested MAF thresholds from publication table exports." if bad_mafs else "No action required.",
    )
    expected_top_rows = expected_rows * 100
    add(
        "PASS" if len(top100) == expected_top_rows else "WARN",
        "Medium",
        "Top-100 table coverage",
        f"Observed {len(top100)} rows; expected approximately {expected_top_rows} rows.",
        "Check traits with fewer than 100 valid variants if this warning appears.",
    )
    missing_mlma = inventory.loc[~inventory["mlma_file_exists"].astype(bool)]
    add(
        "PASS" if missing_mlma.empty else "FAIL",
        "High",
        "MLMA result files exist",
        f"Missing MLMA files: {len(missing_mlma)}.",
        "Regenerate or restore missing GCTA result files." if not missing_mlma.empty else "No action required.",
    )
    bad_lambda = pd.to_numeric(summary["lambda_gc"], errors="coerce").isna().sum()
    add(
        "PASS" if bad_lambda == 0 else "FAIL",
        "High",
        "Lambda GC numeric",
        f"Rows with missing/non-numeric lambda_gc: {bad_lambda}.",
    )
    add(
        "PASS" if len(annotated) >= len(bonf_raw) else "WARN",
        "Medium",
        "Annotated Bonferroni hit table available",
        f"Raw Bonferroni rows: {len(bonf_raw)}; annotated rows: {len(annotated)}.",
        "Use raw hits if annotation matching is intentionally filtered; otherwise inspect missing annotation rows.",
    )
    return pd.DataFrame(rows)


def main() -> None:
    SUPP_DIR.mkdir(parents=True, exist_ok=True)
    EXT_DIR.mkdir(parents=True, exist_ok=True)

    summary = concat_with_run_metadata("outputs/gcta_mlma_loco_combined_summary.tsv")
    top100 = concat_with_run_metadata("outputs/gcta_mlma_loco_top100_per_trait.tsv")
    bonf_raw = concat_with_run_metadata("outputs/gcta_mlma_loco_bonferroni_hits.tsv")
    pheno = concat_with_run_metadata("00_inputs/phenotype_summary.tsv")
    for df in (summary, top100, bonf_raw, pheno):
        df["maf_threshold"] = df["maf_threshold"].map(normalize_maf_threshold)

    top100["rank_within_maf_trait"] = (
        top100.groupby(["maf_threshold", "drug", "mode"]).cumcount() + 1
    )
    top100 = add_neglog10p(top100)
    bonf_raw = add_neglog10p(bonf_raw)
    summary = add_neglog10p(summary, "top_p")

    ann_frames = []
    ann_root = GWAS_ROOT / "gcta_mlma_loco_significant_snp_annotation_csvs"
    for run in RUNS:
        path = ann_root / run["maf_label"] / f"{run['maf_label']}.all_drugs_all_phenotypes.bonferroni_hits.annotated.csv"
        df = read_csv(path)
        ann_frames.append(df)
    annotated = pd.concat(ann_frames, ignore_index=True)
    annotated["maf_threshold"] = annotated["maf_threshold"].map(normalize_maf_threshold)
    annotated = annotated[annotated["maf_label"].isin(ALLOWED_LABELS)].copy()
    annotated = add_neglog10p(annotated)

    inventory = build_run_inventory(summary)
    for col in ["mlma_file", "report_file", "manhattan_png", "qq_png"]:
        inventory[f"{col}_exists"] = inventory[col].map(lambda x: (ROOT / x).exists())

    key_signal = build_key_signal_summary(summary, annotated)
    maf_comp = build_maf_comparison(summary)
    ann_breakdown = build_annotation_breakdown(annotated)
    gene_summary = build_gene_summary(annotated)
    fungamr_ld = build_fungamr_ld_summary(annotated)
    qc = build_qc(summary, top100, bonf_raw, annotated, inventory)

    supp_tables = [
        ("Supplementary_Data_Table_1_gcta_mlma_loco_run_inventory.tsv", "Run inventory and file provenance for selected GCTA MLMA-LOCO GWAS runs", inventory),
        ("Supplementary_Data_Table_2_trait_level_gwas_summary.tsv", "Trait-level GWAS summary metrics for MAF 0.05 and 0.005", summary),
        ("Supplementary_Data_Table_3_top100_variants_per_trait.tsv", "Top 100 variants per drug, phenotype and selected MAF threshold", top100),
        ("Supplementary_Data_Table_4_raw_bonferroni_hits.tsv", "Raw Bonferroni-significant variants for selected MAF thresholds", bonf_raw),
        ("Supplementary_Data_Table_5_annotated_bonferroni_hits.tsv", "SnpEff and FungAMR-LD annotation of Bonferroni-significant variants", annotated),
        ("Supplementary_Data_Table_6_phenotype_input_summary.tsv", "Phenotype input sample counts used by selected GWAS runs", pheno),
    ]
    ext_tables = [
        ("Extended_Table_1_gwas_qc_findings.tsv", "QC findings for selected GCTA MLMA-LOCO publication tables", qc),
        ("Extended_Table_2_key_signal_summary_by_trait.tsv", "Key lead signals and annotation context by drug, phenotype and MAF", key_signal),
        ("Extended_Table_3_maf_0p05_vs_0p005_comparison.tsv", "Trait-level comparison between MAF 0.05 and MAF 0.005 GWAS runs", maf_comp),
        ("Extended_Table_4_bonferroni_annotation_breakdown.tsv", "Bonferroni hit breakdown by annotation effect and impact", ann_breakdown),
        ("Extended_Table_5_gene_level_bonferroni_summary.tsv", "Gene-level recurrence summary for annotated Bonferroni-significant variants", gene_summary),
        ("Extended_Table_6_fungamr_ld_context_summary.tsv", "Summary of Bonferroni hits in LD with high-confidence FungAMR mutations", fungamr_ld),
    ]

    for filename, _title, df in supp_tables:
        write_tsv(df, SUPP_DIR / filename)
    for filename, _title, df in ext_tables:
        write_tsv(df, EXT_DIR / filename)

    supp_manifest = pd.DataFrame(
        [
            {
                "table": f"Supplementary Data Table {i}",
                "filename": filename,
                "title": title,
                "maf_thresholds_included": "0.05;0.005",
            }
            for i, (filename, title, _df) in enumerate(supp_tables, start=1)
        ]
    )
    ext_manifest = pd.DataFrame(
        [
            {
                "table": f"Extended Table {i}",
                "filename": filename,
                "title": title,
                "maf_thresholds_included": "0.05;0.005",
            }
            for i, (filename, title, _df) in enumerate(ext_tables, start=1)
        ]
    )
    write_tsv(supp_manifest, SUPP_DIR / "supplementary_data_table_manifest.tsv")
    write_tsv(ext_manifest, EXT_DIR / "extended_table_manifest.tsv")

    readme = f"""# GCTA MLMA-LOCO Publication Tables

Generated: {datetime.now().isoformat(timespec="seconds")}

Scope: local GCTA `--mlma-loco` GWAS outputs under `{rel(GWAS_ROOT)}`.

Included MAF thresholds: 0.05 and 0.005.

Model: GCTA mixed linear model association with leave-one-chromosome-out genetic relationship handling (`--mlma-loco`) and PC covariates from the run folders.

Folders:
- `supplementary_data_tables`: larger reproducibility/source tables.
- `extended_tables`: compact manuscript-facing summaries and QC tables.
"""
    (GWAS_ROOT / "README_publication_tables.md").write_text(readme, encoding="utf-8")

    report = f"""# GCTA MLMA-LOCO Table QC

Generated: {datetime.now().isoformat(timespec="seconds")}

- Trait-level summary rows: {len(summary)}
- Top-100 rows: {len(top100)}
- Raw Bonferroni-hit rows: {len(bonf_raw)}
- Annotated Bonferroni-hit rows: {len(annotated)}
- QC pass: {(qc['status'] == 'PASS').sum()}
- QC warn: {(qc['status'] == 'WARN').sum()}
- QC fail: {(qc['status'] == 'FAIL').sum()}

The exported publication tables include only MAF 0.05 and MAF 0.005 results.
"""
    (EXT_DIR / "gcta_mlma_loco_publication_table_qc_report.md").write_text(report, encoding="utf-8")

    print(
        {
            "supplementary_dir": rel(SUPP_DIR),
            "extended_dir": rel(EXT_DIR),
            "summary_rows": len(summary),
            "top100_rows": len(top100),
            "raw_bonferroni_rows": len(bonf_raw),
            "annotated_bonferroni_rows": len(annotated),
            "qc_pass": int((qc["status"] == "PASS").sum()),
            "qc_warn": int((qc["status"] == "WARN").sum()),
            "qc_fail": int((qc["status"] == "FAIL").sum()),
        }
    )



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
