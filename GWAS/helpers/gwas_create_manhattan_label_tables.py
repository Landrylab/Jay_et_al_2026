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

"""Create per-plot gene-label tables for GCTA Manhattan plots."""

from __future__ import annotations

import bisect
import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
GENE_TABLE = (
    ROOT
    / "fungamr_snpeff_mapping_outputs_20260423"
    / "nature_microbiology_plots"
    / "cnv_aneuploidy_supplement"
    / "genomewide_gene_locus_cnv_frequency_plot_table.csv"
)
ANNOTATION_ROOT = ROOT / "gcta_mlma_loco/gcta_mlma_loco_significant_snp_annotation_csvs"

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

AZOLE_DRUGS = {"fluconazole", "voriconazole", "posaconazole"}
ECHINOCANDIN_DRUGS = {"caspofungin", "micafungin", "anidulafungin"}

AZOLE_CANONICAL = ["ERG11", "TAC1", "UPC2", "CAP1", "MDR1", "MRR1", "MRR2", "CDR1", "CDR2"]
AZOLE_CANDIDATES = ["DNF2", "FAS1", "GIT2", "GIT3", "GIT4", "ROT2", "HUT1", "SNF6"]

ECHINOCANDIN_CANONICAL = ["FKS1"]
ECHINOCANDIN_CANDIDATES = ["CRZ1", "CCH1", "ECM17", "BCR1", "MNN15", "GUP1", "HSP90"]



# =============================================================================
# Functions
# =============================================================================

def chr_to_numeric(value: str) -> int | None:
    text = str(value).strip()
    text = text.removeprefix("Chr").removeprefix("chr")
    text = text.removeprefix("Ca22chr")
    text = text.split("A_C_albicans_SC5314", 1)[0]
    if text in {"M", "MT"}:
        return 8
    if text in {"R", "RA"}:
        return 9
    try:
        return int(float(text))
    except ValueError:
        return None


def split_names(value: str) -> set[str]:
    names: set[str] = set()
    for token in value.replace(";", ",").split(","):
        token = token.strip().upper()
        if token:
            names.add(token)
    return names


def read_gene_table() -> list[dict[str, str]]:
    with GENE_TABLE.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def row_chr(row: dict[str, str]) -> int | None:
    return chr_to_numeric(row.get("chromosome_label") or row.get("chromosome") or "")


def row_range(row: dict[str, str]) -> tuple[int, int, int] | None:
    chrom = row_chr(row)
    if chrom is None:
        return None
    try:
        start = int(float(row["start_1based_inclusive"]))
        end = int(float(row["end_1based_inclusive"]))
    except (KeyError, ValueError):
        return None
    return chrom, min(start, end), max(start, end)


def exact_symbol(row: dict[str, str]) -> set[str]:
    names = set()
    for col in ("GeneName_Common", "gff_gene_symbol"):
        value = row.get(col, "").strip().upper()
        if value:
            names.add(value)
    return names


def alias_symbol(row: dict[str, str]) -> set[str]:
    names = set()
    for col in ("GeneAlias", "GeneName_Standard", "gene_ids"):
        names.update(split_names(row.get(col, "")))
    return names


def find_gene_ranges(rows: list[dict[str, str]], gene: str) -> list[tuple[int, int, int]]:
    gene_upper = gene.upper()

    if gene_upper == "FKS1":
        # In this genome table, Candida albicans FKS1 is represented by the
        # common/gff symbol GSC1. GSL1 also carries FKS1 as an alias, so keep
        # this mapping explicit.
        matched = [row for row in rows if "GSC1" in exact_symbol(row)]
    else:
        matched = [row for row in rows if gene_upper in exact_symbol(row)]
        if not matched:
            matched = [row for row in rows if gene_upper in alias_symbol(row)]

    ranges = []
    for row in matched:
        interval = row_range(row)
        if interval is not None:
            ranges.append(interval)
    return ranges


def als_family_genes(rows: list[dict[str, str]]) -> list[str]:
    genes = sorted(
        {
            row.get("GeneName_Common", "").strip().upper()
            for row in rows
            if row.get("GeneName_Common", "").strip().upper().startswith("ALS")
            and row.get("GeneName_Common", "").strip()[3:].isdigit()
        },
        key=lambda x: (len(x), x),
    )
    return genes


def read_mlma(path: Path) -> tuple[dict[int, list[dict[str, object]]], dict[int, list[int]], int]:
    by_chr: dict[int, list[dict[str, object]]] = defaultdict(list)
    n_valid = 0
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        col = {name.lower(): i for i, name in enumerate(header)}
        chr_i = col.get("chr")
        snp_i = col.get("snp")
        bp_i = col.get("bp")
        p_i = col.get("p")
        if None in (chr_i, snp_i, bp_i, p_i):
            raise RuntimeError(f"Missing Chr/SNP/bp/p columns in {path}")
        for line in handle:
            fields = line.split()
            if len(fields) <= max(chr_i, snp_i, bp_i, p_i):
                continue
            chrom = chr_to_numeric(fields[chr_i])
            if chrom is None:
                continue
            try:
                bp = int(float(fields[bp_i]))
                p_value = float(fields[p_i])
            except ValueError:
                continue
            if not math.isfinite(p_value) or p_value <= 0 or p_value > 1:
                continue
            by_chr[chrom].append(
                {
                    "chr": chrom,
                    "bp": bp,
                    "snp": fields[snp_i],
                    "p": p_value,
                }
            )
            n_valid += 1
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda item: int(item["bp"]))
    positions_by_chr = {
        chrom: [int(item["bp"]) for item in variants] for chrom, variants in by_chr.items()
    }
    return by_chr, positions_by_chr, n_valid


def best_variant_in_ranges(
    by_chr: dict[int, list[dict[str, object]]],
    positions_by_chr: dict[int, list[int]],
    ranges: list[tuple[int, int, int]],
) -> dict[str, object] | None:
    best = None
    for chrom, start, end in ranges:
        variants = by_chr.get(chrom, [])
        positions = positions_by_chr.get(chrom, [])
        left = bisect.bisect_left(positions, start)
        right = bisect.bisect_right(positions, end)
        for item in variants[left:right]:
            if best is None or float(item["p"]) < float(best["p"]):
                best = item
    return best


def gene_sets_for_drug(drug: str, rows: list[dict[str, str]]) -> tuple[list[str], list[str]]:
    if drug in AZOLE_DRUGS:
        return AZOLE_CANONICAL, AZOLE_CANDIDATES
    if drug in ECHINOCANDIN_DRUGS:
        return ECHINOCANDIN_CANONICAL, ECHINOCANDIN_CANDIDATES + als_family_genes(rows)
    raise ValueError(f"Unknown drug: {drug}")


def write_label_table(
    path: Path,
    rows: list[dict[str, str]],
) -> None:
    fields = ["snp", "label", "gene_role", "p", "bonferroni_p", "significant", "chr", "bp"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def short_gene_id(value: str) -> str:
    for token in value.replace(",", ";").split(";"):
        token = token.strip()
        if token:
            return token
    return ""


def annotation_label(row: dict[str, str]) -> str:
    """Return a plot label only when a curated/common gene name exists."""
    common = row.get("gene_common_name", "").strip()
    if common:
        return common
    return ""


def echinocandin_significant_label_rows(maf_label: str, drug: str, mode: str) -> list[dict[str, str]]:
    annotated = ANNOTATION_ROOT / maf_label / drug / f"{drug}.{mode}.bonferroni_hits.annotated.csv"
    if not annotated.exists():
        return []

    best_by_label: dict[str, dict[str, str]] = {}
    with annotated.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                p_value = float(row["p"])
                bonf_p = float(row["bonferroni_p"])
            except (KeyError, ValueError):
                continue
            if p_value > bonf_p:
                continue

            label = annotation_label(row)
            if not label:
                continue
            role = "canonical" if label.upper() == "FKS1" else "candidate"
            chrom = chr_to_numeric(row.get("chr_label", "") or row.get("chr", ""))
            if chrom is None:
                continue
            try:
                bp = int(float(row.get("bp", "") or row.get("pos_vcf", "")))
            except ValueError:
                continue

            label_row = {
                "snp": row["snp"],
                "label": label,
                "gene_role": role,
                "p": f"{p_value:.10g}",
                "bonferroni_p": f"{bonf_p:.10g}",
                "significant": "yes",
                "chr": str(chrom),
                "bp": str(bp),
            }
            previous = best_by_label.get(label)
            if previous is None or p_value < float(previous["p"]):
                best_by_label[label] = label_row

    return sorted(best_by_label.values(), key=lambda row: (int(row["chr"]), int(row["bp"]), row["label"]))


def main() -> None:
    gene_rows = read_gene_table()
    gene_ranges = {}
    for gene in (
        AZOLE_CANONICAL
        + AZOLE_CANDIDATES
        + ECHINOCANDIN_CANONICAL
        + ECHINOCANDIN_CANDIDATES
        + als_family_genes(gene_rows)
    ):
        gene_ranges[gene] = find_gene_ranges(gene_rows, gene)

    manifest_rows: list[dict[str, str]] = []
    for maf_label, folder in RUNS:
        run_dir = ROOT / folder
        for drug in DRUGS:
            canonical, candidates = gene_sets_for_drug(drug, gene_rows)
            genes = list(dict.fromkeys(canonical + candidates))
            for mode in MODES:
                if drug in ECHINOCANDIN_DRUGS:
                    label_rows = echinocandin_significant_label_rows(maf_label, drug, mode)
                else:
                    mlma = run_dir / "outputs" / drug / mode / "gwas.loco.mlma"
                    if not mlma.exists():
                        continue
                    by_chr, positions_by_chr, n_variants = read_mlma(mlma)
                    bonf_p = 0.05 / n_variants

                    label_rows = []
                    for gene in genes:
                        ranges = gene_ranges.get(gene, [])
                        if not ranges:
                            continue
                        best = best_variant_in_ranges(by_chr, positions_by_chr, ranges)
                        if best is None:
                            continue
                        role = "canonical" if gene in canonical else "candidate"
                        significant = float(best["p"]) <= bonf_p
                        if role != "canonical" and not significant:
                            continue
                        label_rows.append(
                            {
                                "snp": str(best["snp"]),
                                "label": gene,
                                "gene_role": role,
                                "p": f"{float(best['p']):.10g}",
                                "bonferroni_p": f"{bonf_p:.10g}",
                                "significant": "yes" if significant else "no",
                                "chr": str(best["chr"]),
                                "bp": str(best["bp"]),
                            }
                        )

                    label_rows.sort(key=lambda row: (int(row["chr"]), int(row["bp"]), row["label"]))
                out = (
                    run_dir
                    / "outputs"
                    / drug
                    / mode
                    / "plots"
                    / f"manuscript_candidate_labels.{drug}.{mode}.tsv"
                )
                write_label_table(out, label_rows)
                manifest_rows.append(
                    {
                        "maf_label": maf_label,
                        "drug": drug,
                        "mode": mode,
                        "n_labels": str(len(label_rows)),
                        "n_canonical_labels": str(sum(row["gene_role"] == "canonical" for row in label_rows)),
                        "n_candidate_labels": str(sum(row["gene_role"] == "candidate" for row in label_rows)),
                        "label_file": str(out.relative_to(ROOT)),
                    }
                )
                print(f"{maf_label} {drug} {mode}: {len(label_rows)} labels")

    manifest = ROOT / "gcta_mlma_loco/gcta_mlma_loco_manhattan_label_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "maf_label",
            "drug",
            "mode",
            "n_labels",
            "n_canonical_labels",
            "n_candidate_labels",
            "label_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {manifest}")



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
