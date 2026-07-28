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

"""Annotate GCTA MLMA-LOCO Bonferroni hits with SnpEff and FungAMR LD context."""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SNPEFF_VCF = ROOT / "joint_all_snpEff.ann.vcf.gz"
DEFAULT_BASE_PLINK = (
    ROOT
    / "gcta_mlma_loco"
    / "plink_from_joint_vcf"
    / "joint_keep1153"
)
DEFAULT_PLINK_BIN = os.environ.get("PLINK_BIN", "plink")
DEFAULT_OUTDIR = ROOT / "gcta_mlma_loco/gcta_mlma_loco_significant_snp_annotation_csvs"

RUNS = [
    ("0.05", "maf0p05", "gcta_mlma_loco/gcta_mlma_loco_maf005_six_antifungals"),
    ("0.005", "maf0p005", "gcta_mlma_loco/gcta_mlma_loco_maf0005_six_antifungals"),
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

ANN_FIELDS = [
    "allele",
    "annotation",
    "impact",
    "gene_name",
    "gene_id",
    "feature_type",
    "feature_id",
    "transcript_biotype",
    "rank",
    "hgvs_c",
    "hgvs_p",
    "cdna_pos",
    "cds_pos",
    "aa_pos",
    "distance",
    "messages",
]

IMPACT_RANK = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
}

AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
    "*": "*",
}

NONCODING_ANNOTATIONS = {
    "intergenic_region",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
    "intragenic_variant",
    "conserved_intergenic_variant",
}

CODING_EFFECT_TERMS = {
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "initiator_codon_variant",
}



# =============================================================================
# Functions
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--snpeff-vcf", default=str(DEFAULT_SNPEFF_VCF))
    parser.add_argument("--plink-prefix", default=str(DEFAULT_BASE_PLINK))
    parser.add_argument("--plink-bin", default=str(DEFAULT_PLINK_BIN))
    parser.add_argument("--ld-r2-threshold", type=float, default=0.8)
    parser.add_argument("--ld-window-kb", type=int, default=1000)
    return parser.parse_args()


def resolve_plink_prefix(prefix: Path) -> Path:
    candidates = [
        prefix,
        ROOT / "gcta_mlma_loco" / "plink_from_joint_vcf" / "joint_keep1153",
    ]
    for candidate in candidates:
        if Path(str(candidate) + ".bim").exists():
            return candidate
    tried = "\n".join(str(path) for path in candidates)
    raise FileNotFoundError(f"No PLINK prefix found. Tried:\n{tried}")


def parse_snp_id(snp: str) -> tuple[str, str, str, str]:
    parts = snp.split(":")
    if len(parts) != 4:
        raise ValueError(f"Unexpected SNP ID: {snp}")
    return parts[0], parts[1], parts[2], parts[3]


def choose_effective_alt(ref: str, alt_field: str, a1: str, a2: str) -> str:
    alts = [alt.strip() for alt in alt_field.split(",") if alt.strip()]
    if len(alts) <= 1:
        return alts[0] if alts else alt_field

    for allele in (a1, a2):
        if allele != ref and allele in alts:
            return allele
    for allele in (a1, a2):
        if allele in alts:
            return allele
    return alts[0]


def variant_class(ref: str, alt: str) -> str:
    if "*" in {ref, alt}:
        return "spanning_deletion_or_complex"
    if len(ref) == 1 and len(alt) == 1:
        return "point_mutation"
    if len(ref) < len(alt):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    if len(ref) == len(alt):
        return "multi_base_substitution"
    return "complex"


def parse_ann_value(ann_value: str) -> list[dict[str, str]]:
    if not ann_value or ann_value == ".":
        return []

    entries = []
    for ann in ann_value.split(","):
        fields = ann.split("|")
        if len(fields) < len(ANN_FIELDS):
            fields.extend([""] * (len(ANN_FIELDS) - len(fields)))
        elif len(fields) > len(ANN_FIELDS):
            fields = fields[: len(ANN_FIELDS)]
        entries.append(dict(zip(ANN_FIELDS, fields)))
    return entries


def annotation_priority(entry: dict[str, str]) -> tuple[int, int, int, int, int, int, int]:
    terms = set(filter(None, entry["annotation"].split("&")))
    noncoding = terms.issubset(NONCODING_ANNOTATIONS) if terms else True
    return (
        1 if entry["hgvs_p"] else 0,
        IMPACT_RANK.get(entry["impact"], 0),
        1 if entry["transcript_biotype"] == "protein_coding" else 0,
        1 if entry["feature_type"] == "transcript" else 0,
        0 if noncoding else 1,
        1 if entry["gene_id"] or entry["gene_name"] else 0,
        len(entry["annotation"]),
    )


def choose_best_annotation(entries: list[dict[str, str]], alt: str) -> dict[str, str] | None:
    candidates = [entry for entry in entries if entry["allele"] == alt]
    if not candidates:
        return None
    return max(candidates, key=annotation_priority)


def aa_token_to_one_letter(token: str) -> str:
    return AA3_TO_1.get(token, token)


def hgvs_p_to_short(hgvs_p: str) -> str:
    if not hgvs_p:
        return ""

    import re

    text = hgvs_p.strip()
    if text.startswith("p."):
        text = text[2:]
    text = text.strip("()")
    if not text:
        return ""
    if text == "?":
        return "?"

    synonymous = re.fullmatch(r"([A-Z][a-z]{2})(\d+)=", text)
    if synonymous:
        aa, pos = synonymous.groups()
        aa1 = aa_token_to_one_letter(aa)
        return f"{aa1}{pos}{aa1}"

    protein_change = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)(.*)", text)
    if protein_change:
        ref_aa, pos, alt_aa, suffix = protein_change.groups()
        return f"{aa_token_to_one_letter(ref_aa)}{pos}{aa_token_to_one_letter(alt_aa)}{suffix}"

    return text


def blank_annotation(status: str) -> dict[str, str]:
    return {
        "annotation_status": status,
        "snpeff_gene_id": "",
        "snpeff_gene_name": "",
        "effect": "",
        "impact": "",
        "feature_type": "",
        "feature_id": "",
        "transcript_biotype": "",
        "hgvs_c": "",
        "hgvs_p": "",
        "aa_substitution": "",
        "cds_pos": "",
        "aa_pos": "",
        "distance_to_feature": "",
    }


def format_annotation(entry: dict[str, str] | None) -> dict[str, str]:
    if entry is None:
        return blank_annotation("alt_not_annotated")
    return {
        "annotation_status": "matched_ann",
        "snpeff_gene_id": entry["gene_id"],
        "snpeff_gene_name": entry["gene_name"],
        "effect": entry["annotation"],
        "impact": entry["impact"],
        "feature_type": entry["feature_type"],
        "feature_id": entry["feature_id"],
        "transcript_biotype": entry["transcript_biotype"],
        "hgvs_c": entry["hgvs_c"],
        "hgvs_p": entry["hgvs_p"],
        "aa_substitution": hgvs_p_to_short(entry["hgvs_p"]),
        "cds_pos": entry["cds_pos"],
        "aa_pos": entry["aa_pos"],
        "distance_to_feature": entry["distance"],
    }


def write_regions_file(targets: set[tuple[str, str, str, str]], tmpdir: Path) -> Path:
    regions = sorted({(chrom, pos) for chrom, pos, _ref, _alt in targets})
    out = tmpdir / "snpeff_targets.regions.tsv"
    with out.open("w", encoding="utf-8") as handle:
        for chrom, pos in regions:
            handle.write(f"{chrom}\t{pos}\t{pos}\n")
    return out


def load_snpeff_annotations(
    snpeff_vcf: Path,
    targets: set[tuple[str, str, str, str]],
) -> dict[tuple[str, str, str, str], dict[str, str]]:
    annotations = {}
    found = set()
    with tempfile.TemporaryDirectory(prefix="gcta_snpeff_regions_", dir="/private/tmp") as td:
        region_path = write_regions_file(targets, Path(td))
        proc = subprocess.run(
            [
                "bcftools",
                "query",
                "-R",
                str(region_path),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n",
                str(snpeff_vcf),
            ],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

    for line in proc.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, ref, alt_field, ann_value = fields[:5]
        ann_entries = parse_ann_value(ann_value)
        for alt in alt_field.split(","):
            key = (chrom, pos, ref, alt)
            if key not in targets:
                continue
            found.add(key)
            annotations[key] = format_annotation(choose_best_annotation(ann_entries, alt))

    for key in targets - found:
        annotations[key] = blank_annotation("variant_not_found_in_snpeff_vcf")
    return annotations


def normalize_gene_key(value: str) -> str:
    text = (value or "").strip()
    if not text:
        return ""
    if text.endswith("-T"):
        text = text[:-2]
    return text.upper()


def split_aliases(value: str) -> list[str]:
    if not value:
        return []
    tokens = []
    for chunk in value.replace(";", ",").split(","):
        token = chunk.strip()
        if token:
            tokens.append(token)
    return tokens


def load_gene_name_map() -> dict[str, dict[str, str]]:
    maps = {}

    def add(key: str, common: str, standard: str, aliases: str, source: str) -> None:
        norm = normalize_gene_key(key)
        if not norm or norm in maps:
            return
        maps[norm] = {
            "gene_common_name": common or "",
            "gene_standard_name": standard or "",
            "gene_aliases": aliases or "",
            "gene_name_source": source,
        }

    cnv_plot = (
        ROOT
        / "fungamr_snpeff_mapping_outputs_20260423"
        / "nature_microbiology_plots"
        / "cnv_aneuploidy_supplement"
        / "genomewide_gene_locus_cnv_frequency_plot_table.csv"
    )
    if cnv_plot.exists():
        with cnv_plot.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                common = row.get("display_gene") or row.get("GeneName_Common") or row.get("gff_gene_symbol") or ""
                standard = row.get("GeneName_Standard") or ""
                aliases = row.get("GeneAlias") or ""
                keys = [
                    row.get("locus_id", ""),
                    row.get("gene_id_hapA", ""),
                    row.get("gff_gene_symbol", ""),
                    row.get("GeneName_Common", ""),
                    row.get("display_gene", ""),
                ]
                keys.extend(split_aliases(row.get("gene_ids", "")))
                keys.extend(split_aliases(standard))
                keys.extend(split_aliases(aliases))
                for key in keys:
                    add(key, common, standard, aliases, str(cnv_plot))

    gene_cnv = ROOT / "fungamr_snpeff_mapping_outputs_20260423" / "geneCNV_table.annotated.tsv"
    if gene_cnv.exists():
        with gene_cnv.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                common = row.get("GeneName_Common") or ""
                standard = row.get("GeneName_Standard") or ""
                aliases = row.get("GeneAlias") or ""
                keys = [row.get("gene_id", ""), common]
                keys.extend(split_aliases(standard))
                keys.extend(split_aliases(aliases))
                for key in keys:
                    add(key, common, standard, aliases, str(gene_cnv))

    return maps


def resolve_gene_names(annotation: dict[str, str], gene_map: dict[str, dict[str, str]]) -> dict[str, str]:
    keys = [
        annotation.get("snpeff_gene_id", ""),
        annotation.get("snpeff_gene_name", ""),
        annotation.get("feature_id", ""),
    ]
    for key in keys:
        norm = normalize_gene_key(key)
        if norm in gene_map:
            return gene_map[norm]
    return {
        "gene_common_name": "",
        "gene_standard_name": "",
        "gene_aliases": "",
        "gene_name_source": "",
    }


def gene_context(effect: str, transcript_biotype: str, cds_pos: str) -> tuple[str, str]:
    terms = set(filter(None, effect.split("&")))
    in_cds = bool(cds_pos) or bool(terms & CODING_EFFECT_TERMS)
    if in_cds:
        return "yes", "coding_sequence"
    if "upstream_gene_variant" in terms:
        return "no", "upstream_gene_variant"
    if "downstream_gene_variant" in terms:
        return "no", "downstream_gene_variant"
    if "intron_variant" in terms:
        return "no", "intron_variant"
    if "intergenic_region" in terms:
        return "no", "intergenic"
    if transcript_biotype == "protein_coding":
        return "no", "protein_coding_gene_non_cds"
    return "no", "noncoding_or_unresolved"


def read_summary_thresholds(run_folder: Path) -> dict[tuple[str, str], dict[str, str]]:
    out = {}
    summary = run_folder / "outputs" / "gcta_mlma_loco_combined_summary.tsv"
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[(row["drug"], row["mode"])] = row
    return out


def collect_gwas_hits() -> tuple[list[dict[str, str]], set[tuple[str, str, str, str]]]:
    rows = []
    targets = set()
    for maf, maf_label, folder_name in RUNS:
        folder = ROOT / folder_name
        thresholds = read_summary_thresholds(folder)
        hit_path = folder / "outputs" / "gcta_mlma_loco_bonferroni_hits.tsv"
        with hit_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                chrom, pos, ref, alt_field = parse_snp_id(row["snp"])
                alt = choose_effective_alt(ref, alt_field, row.get("a1", ""), row.get("a2", ""))
                summary = thresholds.get((row["drug"], row["mode"]), {})
                row.update(
                    {
                        "maf_threshold": maf,
                        "maf_label": maf_label,
                        "analysis_folder": folder_name,
                        "chrom_vcf": chrom,
                        "pos_vcf": pos,
                        "ref": ref,
                        "alt": alt,
                        "alt_alleles_in_snp_id": alt_field,
                        "variant_class": variant_class(ref, alt),
                        "bonferroni_p": summary.get("bonferroni_p", ""),
                        "bonferroni_log10p": summary.get("bonferroni_log10p", ""),
                        "lambda_gc": summary.get("lambda_gc", ""),
                    }
                )
                rows.append(row)
                targets.add((chrom, pos, ref, alt))
    return rows, targets


def score_is_high_confidence(value: str) -> bool:
    try:
        return float(value) in {1.0, 2.0}
    except ValueError:
        return False


def load_fungamr_records() -> tuple[list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    records = []

    flu = ROOT / "results" / "known_flu_mutation_odds" / "flu_known_mutation_odds.tsv"
    if flu.exists():
        with flu.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                score = row.get("score_Fluconazole", "")
                if not score_is_high_confidence(score):
                    continue
                records.extend(make_fungamr_records(row, "fluconazole", score))

    for path in sorted((ROOT / "results" / "known_mutation_odds").glob("*/*_known_mutation_odds.tsv")):
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                score = row.get("drug_score", "")
                if not score_is_high_confidence(score):
                    continue
                records.extend(make_fungamr_records(row, (row.get("drug_name") or path.parent.name).lower(), score))

    by_locus = defaultdict(list)
    seen = set()
    deduped = []
    for rec in records:
        key = (rec["drug"], rec["locus_id"], rec["fungamr_record"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(rec)
        by_locus[rec["locus_id"]].append(rec)

    return deduped, by_locus


def make_fungamr_records(row: dict[str, str], drug: str, score: str) -> list[dict[str, str]]:
    gene = (row.get("ortho_homolog") or row.get("ANN_feature") or "").upper()
    mutation = row.get("mutation") or ""
    label = f"{gene} {mutation}".strip()
    loci = [item.strip() for item in (row.get("locus_id") or "").split(";") if item.strip()]
    out = []
    for locus in loci:
        out.append(
            {
                "drug": drug.lower(),
                "fungamr_record": label,
                "fungamr_gene": gene,
                "fungamr_substitution": mutation,
                "fungamr_score": score,
                "locus_id": locus,
            }
        )
    return out


def read_bim_ids(prefix: Path) -> set[str]:
    ids = set()
    with Path(str(prefix) + ".bim").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                ids.add(fields[1])
    return ids


def compute_ld(
    plink_bin: Path,
    plink_prefix: Path,
    hit_ids: set[str],
    fungamr_records: list[dict[str, str]],
    r2_threshold: float,
    window_kb: int,
) -> dict[tuple[str, str], list[dict[str, str]]]:
    bim_ids = read_bim_ids(plink_prefix)
    fungamr_loci = {rec["locus_id"] for rec in fungamr_records}
    available_fungamr_loci = fungamr_loci & bim_ids
    available_extract_ids = (hit_ids | fungamr_loci) & bim_ids

    records_by_locus = defaultdict(list)
    for rec in fungamr_records:
        records_by_locus[rec["locus_id"]].append(rec)

    ld_hits = defaultdict(list)
    for hit_id in hit_ids:
        for rec in records_by_locus.get(hit_id, []):
            ld_hits[(hit_id, rec["drug"])].append({**rec, "r2": "1"})

    if not available_fungamr_loci:
        return ld_hits

    with tempfile.TemporaryDirectory(prefix="gcta_fungamr_ld_", dir="/private/tmp") as td:
        tmpdir = Path(td)
        extract_path = tmpdir / "extract_ids.txt"
        fungamr_path = tmpdir / "fungamr_loci.txt"
        with extract_path.open("w", encoding="utf-8") as handle:
            for snp_id in sorted(available_extract_ids):
                handle.write(f"{snp_id}\n")
        with fungamr_path.open("w", encoding="utf-8") as handle:
            for snp_id in sorted(available_fungamr_loci):
                handle.write(f"{snp_id}\n")

        out_prefix = tmpdir / "fungamr_ld"
        subprocess.run(
            [
                str(plink_bin),
                "--bfile",
                str(plink_prefix),
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
                str(window_kb),
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
            return ld_hits

        with ld_path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
            for row in reader:
                snp_a = row.get("SNP_A", "")
                snp_b = row.get("SNP_B", "")
                try:
                    r2 = float(row.get("R2", "nan"))
                except ValueError:
                    continue
                if not math.isfinite(r2) or r2 < r2_threshold:
                    continue
                for fung_locus, hit_id in ((snp_a, snp_b), (snp_b, snp_a)):
                    if fung_locus not in records_by_locus or hit_id not in hit_ids:
                        continue
                    for rec in records_by_locus[fung_locus]:
                        ld_hits[(hit_id, rec["drug"])].append({**rec, "r2": f"{r2:.6g}"})

    return ld_hits


def collapse_ld_records(records: list[dict[str, str]]) -> dict[str, str]:
    if not records:
        return {
            "in_ld_with_fungamr_score1_or_2": "no",
            "fungamr_ld_records": "",
            "fungamr_ld_genes": "",
            "fungamr_ld_substitutions": "",
            "fungamr_ld_scores": "",
            "fungamr_ld_loci": "",
            "fungamr_ld_r2": "",
            "fungamr_ld_max_r2": "",
        }

    best = {}
    for rec in records:
        key = (rec["fungamr_record"], rec["fungamr_score"], rec["locus_id"])
        r2 = float(rec["r2"])
        if key not in best or r2 > float(best[key]["r2"]):
            best[key] = rec
    ordered = sorted(best.values(), key=lambda r: (-float(r["r2"]), r["fungamr_record"], r["locus_id"]))
    return {
        "in_ld_with_fungamr_score1_or_2": "yes",
        "fungamr_ld_records": ", ".join(dict.fromkeys(rec["fungamr_record"] for rec in ordered)),
        "fungamr_ld_genes": ", ".join(dict.fromkeys(rec["fungamr_gene"] for rec in ordered)),
        "fungamr_ld_substitutions": ", ".join(dict.fromkeys(rec["fungamr_substitution"] for rec in ordered)),
        "fungamr_ld_scores": ", ".join(dict.fromkeys(rec["fungamr_score"] for rec in ordered)),
        "fungamr_ld_loci": ", ".join(dict.fromkeys(rec["locus_id"] for rec in ordered)),
        "fungamr_ld_r2": ", ".join(f"{rec['fungamr_record']}={rec['r2']}" for rec in ordered),
        "fungamr_ld_max_r2": f"{max(float(rec['r2']) for rec in ordered):.6g}",
    }


def output_fieldnames() -> list[str]:
    return [
        "maf_threshold",
        "maf_label",
        "analysis_folder",
        "drug",
        "mode",
        "chr",
        "chr_label",
        "chrom_vcf",
        "bp",
        "pos_vcf",
        "snp",
        "ref",
        "alt",
        "alt_alleles_in_snp_id",
        "a1",
        "a2",
        "variant_class",
        "freq",
        "b",
        "se",
        "p",
        "bonferroni_p",
        "bonferroni_log10p",
        "lambda_gc",
        "annotation_status",
        "effect",
        "impact",
        "is_in_coding_sequence",
        "gene_context",
        "snpeff_gene_id",
        "snpeff_gene_name",
        "gene_common_name",
        "gene_standard_name",
        "gene_aliases",
        "feature_type",
        "feature_id",
        "transcript_biotype",
        "hgvs_c",
        "hgvs_p",
        "aa_substitution",
        "cds_pos",
        "aa_pos",
        "distance_to_feature",
        "ld_r2_threshold",
        "ld_window_kb",
        "in_ld_with_fungamr_score1_or_2",
        "fungamr_ld_records",
        "fungamr_ld_genes",
        "fungamr_ld_substitutions",
        "fungamr_ld_scores",
        "fungamr_ld_loci",
        "fungamr_ld_r2",
        "fungamr_ld_max_r2",
    ]


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    snpeff_vcf = Path(args.snpeff_vcf)
    plink_prefix = Path(args.plink_prefix)
    plink_prefix = resolve_plink_prefix(plink_prefix)
    plink_bin = Path(args.plink_bin)

    hits, targets = collect_gwas_hits()
    hit_ids = {row["snp"] for row in hits}
    print(f"Collected {len(hits)} Bonferroni-hit rows across {len(hit_ids)} unique SNPs.")

    annotations = load_snpeff_annotations(snpeff_vcf, targets)
    gene_map = load_gene_name_map()
    fungamr_records, _records_by_locus = load_fungamr_records()
    print(f"Loaded {len(fungamr_records)} FungAMR score 1/2 substitution loci.")

    ld_hits = compute_ld(
        plink_bin=plink_bin,
        plink_prefix=plink_prefix,
        hit_ids=hit_ids,
        fungamr_records=fungamr_records,
        r2_threshold=args.ld_r2_threshold,
        window_kb=args.ld_window_kb,
    )

    annotated_rows = []
    for row in hits:
        key = parse_snp_id(row["snp"])
        ann = annotations.get(key, blank_annotation("variant_not_found_in_snpeff_vcf"))
        resolved_names = resolve_gene_names(ann, gene_map)
        in_cds, context = gene_context(ann["effect"], ann["transcript_biotype"], ann["cds_pos"])
        ld = collapse_ld_records(ld_hits.get((row["snp"], row["drug"].lower()), []))
        annotated = dict(row)
        annotated.update(ann)
        annotated.update(resolved_names)
        annotated.update(
            {
                "is_in_coding_sequence": in_cds,
                "gene_context": context,
                "ld_r2_threshold": f"{args.ld_r2_threshold:g}",
                "ld_window_kb": str(args.ld_window_kb),
            }
        )
        annotated.update(ld)
        annotated_rows.append(annotated)

    annotated_rows.sort(
        key=lambda r: (
            float(r["maf_threshold"]),
            r["drug"],
            r["mode"],
            float(r["p"]),
            r["snp"],
        )
    )

    fieldnames = output_fieldnames()
    combined_path = outdir / "all_gcta_mlma_loco_bonferroni_hits.annotated.csv"
    write_csv(combined_path, annotated_rows, fieldnames)

    manifest_rows = []
    by_analysis = defaultdict(list)
    by_maf = defaultdict(list)
    for row in annotated_rows:
        by_analysis[(row["maf_label"], row["drug"], row["mode"])].append(row)
        by_maf[row["maf_label"]].append(row)

    for _maf, maf_label, _folder in RUNS:
        maf_rows = by_maf.get(maf_label, [])
        maf_path = outdir / maf_label / f"{maf_label}.all_drugs_all_phenotypes.bonferroni_hits.annotated.csv"
        write_csv(maf_path, maf_rows, fieldnames)
        manifest_rows.append(
            {
                "level": "maf_run",
                "maf_label": maf_label,
                "drug": "all",
                "mode": "all",
                "n_bonferroni_hits": str(len(maf_rows)),
                "n_hits_in_ld_with_fungamr_score1_or_2": str(
                    sum(row["in_ld_with_fungamr_score1_or_2"] == "yes" for row in maf_rows)
                ),
                "csv_path": str(maf_path.relative_to(ROOT)),
            }
        )
        for drug in DRUGS:
            for mode in MODES:
                rows = by_analysis.get((maf_label, drug, mode), [])
                out_path = outdir / maf_label / drug / f"{drug}.{mode}.bonferroni_hits.annotated.csv"
                write_csv(out_path, rows, fieldnames)
                manifest_rows.append(
                    {
                        "level": "maf_drug_phenotype",
                        "maf_label": maf_label,
                        "drug": drug,
                        "mode": mode,
                        "n_bonferroni_hits": str(len(rows)),
                        "n_hits_in_ld_with_fungamr_score1_or_2": str(
                            sum(1 for row in rows if row["in_ld_with_fungamr_score1_or_2"] == "yes")
                        ),
                        "csv_path": str(out_path.relative_to(ROOT)),
                    }
                )

    write_csv(
        outdir / "annotation_csv_manifest.csv",
        manifest_rows,
        [
            "level",
            "maf_label",
            "drug",
            "mode",
            "n_bonferroni_hits",
            "n_hits_in_ld_with_fungamr_score1_or_2",
            "csv_path",
        ],
    )

    readme = outdir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# GCTA MLMA-LOCO Bonferroni Hit Annotation CSVs",
                "",
                "These CSVs annotate Bonferroni-significant SNPs from the completed GCTA MLMA-LOCO Manhattan plot runs.",
                "",
                f"- Combined CSV: `{combined_path.name}`",
                "- Per-analysis CSVs are under `<maf_label>/<drug>/`.",
                "- Per-MAF-run aggregate CSVs are under `<maf_label>/`.",
                "- SnpEff annotations are from `joint_all_snpEff.ann.vcf.gz`.",
                "- Gene common names are resolved from the existing CNV/gene annotation tables when possible.",
                "- FungAMR LD records are restricted to score 1 or 2 substitutions for the same drug.",
                f"- LD flag uses PLINK r2 >= {args.ld_r2_threshold:g} within {args.ld_window_kb} kb.",
                "- If a hit is linked to multiple FungAMR records, records are comma-collapsed in the same row.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote combined CSV: {combined_path}")
    print(f"Wrote per-analysis CSVs and manifest under: {outdir}")



# =============================================================================
# Main block
# =============================================================================

if __name__ == "__main__":
    main()
