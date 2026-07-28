#!/usr/bin/env python3
"""Export compact GWAS tables and figure source data for submission.

The input results directory must contain the two completed GCTA MLMA-LOCO
folders for MAF 0.05 and 0.005. The figure directory must contain the
FungAMR LD-detail files produced by step 06.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
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
RUNS = [
    ("0.05", "maf0p05", "gcta_mlma_loco_maf005_six_antifungals"),
    ("0.005", "maf0p005", "gcta_mlma_loco_maf0005_six_antifungals"),
]
SOURCE_FILES = {
    ("0.05", "rs_binary"): "Figure_2_GWAS_MAF0p05_binary.tsv.gz",
    ("0.05", "mic"): "Supplementary_Figure_GWAS_MAF0p05_MIC.tsv.gz",
    ("0.005", "rs_binary"): "Supplementary_Figure_GWAS_MAF0p005_binary.tsv.gz",
    ("0.005", "mic"): "Supplementary_Figure_GWAS_MAF0p005_MIC.tsv.gz",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            key, sep, value = line.rstrip("\n").partition("\t")
            if sep:
                values[key] = value
    return values


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def portable_figure_path(value: str) -> str:
    """Remove workstation-specific prefixes from a figure-manifest path."""
    path = Path(value)
    if "per_trait_plots" in path.parts:
        start = path.parts.index("per_trait_plots")
        return "/".join(path.parts[start:])
    return path.name


def export_figure_manifest(source: Path, destination: Path) -> None:
    rows = read_tsv(source)
    for row in rows:
        row["output"] = portable_figure_path(row.get("output", ""))
        for field in ("panel_A_sources", "panel_B_sources"):
            row[field] = ";".join(
                portable_figure_path(item)
                for item in row.get(field, "").split(";")
                if item
            )
    write_tsv(destination, rows)


def copy_with_run_metadata(
    results_root: Path,
    relative_path: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for maf, maf_label, folder in RUNS:
        source = results_root / folder / relative_path
        for row in read_tsv(source):
            rows.append({"maf_threshold": maf, "maf_label": maf_label, **row})
    return rows


def load_plot_labels(results_root: Path) -> tuple[list[dict[str, str]], dict[tuple[str, str, str, str], dict[str, str]]]:
    rows: list[dict[str, str]] = []
    lookup: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for maf, maf_label, folder in RUNS:
        run_dir = results_root / folder
        for drug in DRUGS:
            for mode in MODES:
                path = (
                    run_dir
                    / "outputs"
                    / drug
                    / mode
                    / "plots"
                    / f"manuscript_candidate_labels.{drug}.{mode}.tsv"
                )
                if not path.exists():
                    continue
                for row in read_tsv(path):
                    out = {
                        "maf_threshold": maf,
                        "maf_label": maf_label,
                        "drug": drug,
                        "phenotype": mode,
                        **row,
                    }
                    rows.append(out)
                    lookup[(maf_label, drug, mode, row["snp"])] = row
    return rows, lookup


def load_fungamr_ld(
    figure_root: Path,
) -> tuple[list[dict[str, str]], dict[tuple[str, str, str], bool]]:
    rows: list[dict[str, str]] = []
    lookup: dict[tuple[str, str, str], bool] = {}
    for maf, maf_label, _folder in RUNS:
        for drug in DRUGS:
            path = (
                figure_root
                / "fungamr_ld_grey_snps"
                / maf_label
                / f"{drug}.ld_r2_ge0p2.window1000kb.details.tsv"
            )
            if not path.exists():
                continue
            for row in read_tsv(path):
                rows.append(
                    {
                        "maf_threshold": maf,
                        "maf_label": maf_label,
                        **row,
                    }
                )
                lookup[(maf_label, drug, row["snp"])] = True
    return rows, lookup


def load_annotated_hits(results_root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    annotation_root = results_root / "gcta_mlma_loco_significant_snp_annotation_csvs"
    for maf, maf_label, _folder in RUNS:
        path = (
            annotation_root
            / maf_label
            / f"{maf_label}.all_drugs_all_phenotypes.bonferroni_hits.annotated.csv"
        )
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                rows.append({"maf_threshold": maf, "maf_label": maf_label, **row})
    return rows


def valid_mlma_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                p_value = float(row["p"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(p_value) and 0 < p_value <= 1:
                row["_p_numeric"] = repr(p_value)
                rows.append(row)
    rows.sort(key=lambda row: float(row["_p_numeric"]))
    return rows


def export_source_data(
    results_root: Path,
    source_dir: Path,
    label_lookup: dict[tuple[str, str, str, str], dict[str, str]],
    ld_lookup: dict[tuple[str, str, str], bool],
) -> list[dict[str, str]]:
    source_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, str]] = []
    fields = [
        "maf_threshold",
        "drug",
        "phenotype",
        "chr",
        "snp",
        "bp",
        "effect_allele",
        "other_allele",
        "effect_allele_frequency",
        "beta",
        "standard_error",
        "p_value",
        "minus_log10_p",
        "bonferroni_p",
        "lambda_gc",
        "qq_rank",
        "qq_expected_minus_log10_p",
        "in_fungamr_ld",
        "plot_label",
        "label_role",
    ]

    for maf, maf_label, folder in RUNS:
        run_dir = results_root / folder
        for mode in MODES:
            output_path = source_dir / SOURCE_FILES[(maf, mode)]
            n_rows = 0
            with gzip.open(output_path, "wt", encoding="utf-8", newline="") as out:
                writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                for drug in DRUGS:
                    result_dir = run_dir / "outputs" / drug / mode
                    mlma_rows = valid_mlma_rows(result_dir / "gwas.loco.mlma")
                    plot_values = read_key_values(
                        result_dir / "plots" / f"plot_summary.{drug}.{mode}.txt"
                    )
                    n_valid = len(mlma_rows)
                    for rank, row in enumerate(mlma_rows, start=1):
                        p_value = float(row["_p_numeric"])
                        label = label_lookup.get((maf_label, drug, mode, row["SNP"]), {})
                        writer.writerow(
                            {
                                "maf_threshold": maf,
                                "drug": drug,
                                "phenotype": mode,
                                "chr": row.get("Chr", ""),
                                "snp": row.get("SNP", ""),
                                "bp": row.get("bp", ""),
                                "effect_allele": row.get("A1", ""),
                                "other_allele": row.get("A2", ""),
                                "effect_allele_frequency": row.get("Freq", ""),
                                "beta": row.get("b", ""),
                                "standard_error": row.get("se", ""),
                                "p_value": row.get("p", ""),
                                "minus_log10_p": f"{-math.log10(p_value):.12g}",
                                "bonferroni_p": plot_values.get("Bonferroni_P", ""),
                                "lambda_gc": plot_values.get("Lambda_GC", ""),
                                "qq_rank": str(rank),
                                "qq_expected_minus_log10_p": (
                                    f"{-math.log10((rank - 0.5) / n_valid):.12g}"
                                ),
                                "in_fungamr_ld": (
                                    "yes" if ld_lookup.get((maf_label, drug, row["SNP"])) else "no"
                                ),
                                "plot_label": label.get("label", ""),
                                "label_role": label.get("gene_role", ""),
                            }
                        )
                        n_rows += 1
            manifest.append(
                {
                    "filename": output_path.name,
                    "maf_threshold": maf,
                    "phenotype": mode,
                    "drugs": ";".join(DRUGS),
                    "n_rows": str(n_rows),
                    "description": (
                        "Exact valid-P GWAS points used for Manhattan and quantile-quantile plots"
                    ),
                }
            )
    return manifest


def build_parameters(results_root: Path) -> list[dict[str, str]]:
    rows = [
        ("cohort_size", "1153", "isolates"),
        ("antifungals", "fluconazole;voriconazole;posaconazole;caspofungin;micafungin;anidulafungin", ""),
        ("phenotypes", "log2_MIC;binary_R_S", ""),
        ("binary_coding", "R=1;S=0;missing=-9", ""),
        ("association_model", "GCTA MLMA-LOCO linear mixed model", ""),
        ("genotype_model", "additive 0/1/2 effect-allele dosage", ""),
        ("gcta_version", "1.95.0", ""),
        ("plink_version", "1.9.0-b.8", ""),
        ("bcftools_version", "1.21", ""),
        ("snpeff_version", "5.2", ""),
        ("maf_thresholds", "0.05;0.005", ""),
        ("ld_pruning", "--indep-pairwise 50 5 0.2", ""),
        ("pca_covariates", "PC1-PC10", ""),
        ("multiple_testing", "Bonferroni 0.05 / number of valid SNP tests per analysis", ""),
        ("lambda_gc", "median chi-square(1 df) / median of chi-square(1 df)", ""),
        ("fungamr_confidence_scores", "1;2", ""),
        ("fungamr_ld_threshold", "r2>=0.2", ""),
        ("fungamr_ld_window", "1000", "kb"),
    ]
    for maf, maf_label, folder in RUNS:
        run_dir = results_root / folder
        filtered_bim = next(run_dir.glob("01_*_filtered_snps/*.bim"))
        pruned_bim = next(run_dir.glob("02_ld_pruned_grm/*.ldpruned.bim"))
        filtered_n = sum(1 for _ in filtered_bim.open("r", encoding="utf-8"))
        pruned_n = sum(1 for _ in pruned_bim.open("r", encoding="utf-8"))
        rows.extend(
            [
                (f"{maf_label}_association_variants", str(filtered_n), "variants"),
                (f"{maf_label}_ld_pruned_pca_variants", str(pruned_n), "variants"),
            ]
        )
    return [
        {"parameter": parameter, "value": value, "unit": unit}
        for parameter, value, unit in rows
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--figure-root", required=True, type=Path)
    parser.add_argument("--submission-root", required=True, type=Path)
    parser.add_argument(
        "--skip-source-data",
        action="store_true",
        help="Export compact tables without the four large compressed point tables.",
    )
    args = parser.parse_args()

    table_dir = args.submission_root / "Supp_data" / "tables"
    source_dir = args.submission_root / "Source_data"
    table_dir.mkdir(parents=True, exist_ok=True)

    summary = copy_with_run_metadata(
        args.results_root,
        "outputs/gcta_mlma_loco_combined_summary.tsv",
    )
    top100 = copy_with_run_metadata(
        args.results_root,
        "outputs/gcta_mlma_loco_top100_per_trait.tsv",
    )
    phenotype_all = copy_with_run_metadata(
        args.results_root,
        "00_inputs/phenotype_summary.tsv",
    )
    phenotype = [
        row for row in phenotype_all if row["maf_threshold"] == "0.05"
    ]
    labels, label_lookup = load_plot_labels(args.results_root)
    fungamr_ld, ld_lookup = load_fungamr_ld(args.figure_root)
    annotated_hits = load_annotated_hits(args.results_root)
    parameters = build_parameters(args.results_root)

    write_tsv(table_dir / "gwas_analysis_summary.tsv", summary)
    write_tsv(table_dir / "gwas_phenotype_summary.tsv", phenotype)
    write_tsv(table_dir / "gwas_top100_variants.tsv", top100)
    write_tsv(table_dir / "gwas_bonferroni_significant_variants.tsv", annotated_hits)
    write_tsv(table_dir / "gwas_plot_labels.tsv", labels)
    write_tsv(table_dir / "gwas_fungamr_ld_context.tsv", fungamr_ld)
    write_tsv(table_dir / "gwas_run_parameters.tsv", parameters)

    if not args.skip_source_data:
        manifest = export_source_data(
            args.results_root,
            source_dir,
            label_lookup,
            ld_lookup,
        )
        write_tsv(source_dir / "gwas_source_data_manifest.tsv", manifest)

    figure_manifest = args.figure_root / "final_gwas_paper_figure_manifest.tsv"
    if figure_manifest.exists():
        export_figure_manifest(
            figure_manifest,
            table_dir / "gwas_figure_manifest.tsv",
        )

    print(f"GWAS tables: {table_dir}")
    if not args.skip_source_data:
        print(f"GWAS source data: {source_dir}")


if __name__ == "__main__":
    main()
