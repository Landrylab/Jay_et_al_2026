#!/usr/bin/env python3
"""Build grouped paired-RF feature rankings and regional carrier support tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact, mannwhitneyu


NONSYNONYMOUS_EFFECTS = {
    "missense_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "conservative_inframe_insertion",
    "conservative_inframe_deletion",
    "disruptive_inframe_insertion",
    "disruptive_inframe_deletion",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paired-feature-summary", required=True)
    parser.add_argument("--binary-gwas", required=True, help="MAF >=0.05 binary fluconazole GWAS TSV.")
    parser.add_argument("--candidate-windows", required=True)
    parser.add_argument("--annotated-vcf", required=True, help="SnpEff-annotated cohort VCF.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fungamr-carriers", required=True, help="TSV with a sample column.")
    parser.add_argument("--sample-col", default="sample")
    parser.add_argument("--label-col", default="rbin_flu")
    parser.add_argument("--log2-mic-col", default="log2_flu")
    parser.add_argument("--output", required=True)
    parser.add_argument("--minimum-carriers", type=int, default=5)
    parser.add_argument("--minimum-log2-mic-difference", type=float, default=2.0)
    parser.add_argument("--minimum-resistant-fraction-difference", type=float, default=0.2)
    parser.add_argument(
        "--exclude-regional-region",
        default="chrom 5 peak",
        help="Candidate region retained in the landscape but omitted from the regional support panel.",
    )
    return parser.parse_args()


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return pd.Series(result, index=p_values.index)
    valid_values = values[valid]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[np.flatnonzero(valid)] = restored
    return pd.Series(result, index=p_values.index)


def normalize_variant_id(chrom: object, pos: object, ref: object, alt: object) -> str:
    return f"{chrom}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def build_landscape(
    feature_path: str | Path,
    gwas_path: str | Path,
    windows: pd.DataFrame,
) -> pd.DataFrame:
    features = pd.read_csv(feature_path, sep="\t", low_memory=False)
    features = features.loc[features["cv_strategy"].eq("group_stratified")].copy()
    if "variant_id" not in features:
        features["variant_id"] = [
            normalize_variant_id(chrom, pos, ref, alt)
            for chrom, pos, ref, alt in features[["chrom", "pos", "ref", "alt"]].itertuples(index=False)
        ]
    features["rf_importance_rank"] = features["mean_importance_across_all_folds"].rank(
        method="first", ascending=False
    ).astype(int)
    rf_cutoff = features["mean_importance_across_all_folds"].quantile(0.99)
    features["rf_top_one_percent"] = features["mean_importance_across_all_folds"] >= rf_cutoff

    gwas = pd.read_csv(gwas_path, sep="\t", low_memory=False)
    id_column = next((column for column in ["variant_id", "SNP", "snp"] if column in gwas), None)
    p_column = next((column for column in ["P", "p_value", "p"] if column in gwas), None)
    if id_column is None or p_column is None:
        raise ValueError("GWAS table requires a variant identifier and P-value column.")
    gwas = gwas[[id_column, p_column]].rename(columns={id_column: "variant_id", p_column: "gwas_p"})
    gwas["gwas_p"] = pd.to_numeric(gwas["gwas_p"], errors="coerce")
    gwas["gwas_minus_log10_p"] = -np.log10(gwas["gwas_p"].clip(lower=np.finfo(float).tiny))
    gwas_cutoff = gwas["gwas_minus_log10_p"].quantile(0.99)
    gwas["gwas_top_one_percent"] = gwas["gwas_minus_log10_p"] >= gwas_cutoff
    landscape = features.merge(gwas, on="variant_id", how="left")
    landscape["gwas_top_one_percent"] = landscape["gwas_top_one_percent"].fillna(False)
    landscape["top_one_percent_overlap"] = (
        landscape["rf_top_one_percent"] & landscape["gwas_top_one_percent"]
    )
    landscape["candidate_region"] = ""
    for row in windows.itertuples(index=False):
        mask = (
            landscape["chrom"].astype(str).eq(str(row.chrom))
            & landscape["pos"].between(int(row.start), int(row.end))
        )
        landscape.loc[mask, "candidate_region"] = str(row.region_id)
    return landscape.sort_values(["chrom", "pos", "variant_id"])


def ann_fields(record: pysam.VariantRecord, alt: str) -> list[dict[str, str]]:
    annotations = record.info.get("ANN", ())
    if isinstance(annotations, str):
        annotations = [annotations]
    parsed: list[dict[str, str]] = []
    for annotation in annotations:
        fields = str(annotation).split("|")
        if len(fields) < 11 or fields[0] != alt:
            continue
        effects = set(fields[1].split("&"))
        if not effects.intersection(NONSYNONYMOUS_EFFECTS):
            continue
        parsed.append(
            {
                "effect": fields[1],
                "impact": fields[2],
                "gene": fields[3] or fields[4],
                "gene_id": fields[4],
                "transcript": fields[6],
                "hgvs_c": fields[9],
                "hgvs_p": fields[10],
            }
        )
    return parsed


def regional_variant_rows(
    vcf_path: str | Path,
    windows: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    label_col: str,
    mic_col: str,
) -> pd.DataFrame:
    vcf = pysam.VariantFile(str(vcf_path))
    metadata = metadata.set_index(sample_col)
    available_samples = [sample for sample in vcf.header.samples if sample in metadata.index]
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, int, str, str, str, str]] = set()
    for window in windows.itertuples(index=False):
        for record in vcf.fetch(str(window.chrom), int(window.start) - 1, int(window.end)):
            if len(record.alts or ()) != 1:
                continue
            alt = str(record.alts[0])
            for annotation in ann_fields(record, alt):
                key = (
                    str(window.region_id),
                    int(record.pos),
                    str(record.ref),
                    alt,
                    annotation["gene"],
                    annotation["effect"],
                )
                if key in seen:
                    continue
                seen.add(key)
                carriers: list[str] = []
                noncarriers: list[str] = []
                for sample in available_samples:
                    genotype = record.samples[sample].get("GT")
                    if genotype is None or any(allele is None for allele in genotype):
                        continue
                    if 1 in genotype:
                        carriers.append(sample)
                    else:
                        noncarriers.append(sample)
                if not carriers or not noncarriers:
                    continue
                carrier_meta = metadata.loc[carriers]
                noncarrier_meta = metadata.loc[noncarriers]
                carrier_r = int(carrier_meta[label_col].astype(str).eq("R").sum())
                carrier_s = int(carrier_meta[label_col].astype(str).eq("S").sum())
                noncarrier_r = int(noncarrier_meta[label_col].astype(str).eq("R").sum())
                noncarrier_s = int(noncarrier_meta[label_col].astype(str).eq("S").sum())
                carrier_fraction = carrier_r / (carrier_r + carrier_s)
                noncarrier_fraction = noncarrier_r / (noncarrier_r + noncarrier_s)
                carrier_mic = pd.to_numeric(carrier_meta[mic_col], errors="coerce").dropna()
                noncarrier_mic = pd.to_numeric(noncarrier_meta[mic_col], errors="coerce").dropna()
                fisher = fisher_exact(
                    [[carrier_r, carrier_s], [noncarrier_r, noncarrier_s]],
                    alternative="two-sided",
                )
                mann_p = (
                    float(mannwhitneyu(carrier_mic, noncarrier_mic, alternative="two-sided").pvalue)
                    if len(carrier_mic) and len(noncarrier_mic)
                    else np.nan
                )
                rows.append(
                    {
                        "region_id": str(window.region_id),
                        "chrom": str(record.chrom),
                        "pos": int(record.pos),
                        "ref": str(record.ref),
                        "alt": alt,
                        "variant_id": normalize_variant_id(record.chrom, record.pos, record.ref, alt),
                        **annotation,
                        "called_isolates": len(carriers) + len(noncarriers),
                        "carrier_total": len(carriers),
                        "carrier_resistant": carrier_r,
                        "carrier_susceptible": carrier_s,
                        "carrier_resistant_fraction": carrier_fraction,
                        "noncarrier_total": len(noncarriers),
                        "noncarrier_resistant": noncarrier_r,
                        "noncarrier_susceptible": noncarrier_s,
                        "noncarrier_resistant_fraction": noncarrier_fraction,
                        "delta_resistant_fraction": carrier_fraction - noncarrier_fraction,
                        "fisher_odds_ratio": float(fisher.statistic),
                        "fisher_exact_p": float(fisher.pvalue),
                        "carrier_log2_mic_n": len(carrier_mic),
                        "carrier_median_log2_mic": float(carrier_mic.median()) if len(carrier_mic) else np.nan,
                        "noncarrier_log2_mic_n": len(noncarrier_mic),
                        "noncarrier_median_log2_mic": (
                            float(noncarrier_mic.median()) if len(noncarrier_mic) else np.nan
                        ),
                        "delta_median_log2_mic": (
                            float(carrier_mic.median() - noncarrier_mic.median())
                            if len(carrier_mic) and len(noncarrier_mic)
                            else np.nan
                        ),
                        "mann_whitney_p": mann_p,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    windows = pd.read_csv(args.candidate_windows, sep="\t")
    required_windows = {"region_id", "chrom", "start", "end"}
    missing_windows = required_windows.difference(windows.columns)
    if missing_windows:
        raise ValueError(f"Candidate-window table is missing columns: {sorted(missing_windows)}")
    landscape = build_landscape(args.paired_feature_summary, args.binary_gwas, windows)
    landscape.to_csv(output / "grouped_paired_rf_gwas_landscape.tsv", sep="\t", index=False)

    metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    carriers = pd.read_csv(args.fungamr_carriers, sep="\t")
    if args.sample_col not in carriers:
        raise ValueError(f"FungAMR carrier table requires column: {args.sample_col}")
    metadata[args.sample_col] = metadata[args.sample_col].astype(str)
    metadata = metadata.loc[
        metadata[args.label_col].astype(str).isin(["R", "S"])
        & ~metadata[args.sample_col].isin(carriers[args.sample_col].astype(str))
    ].copy()
    tested = regional_variant_rows(
        args.annotated_vcf,
        windows,
        metadata,
        args.sample_col,
        args.label_col,
        args.log2_mic_col,
    )
    tested["fisher_exact_q_bh"] = bh_adjust(tested["fisher_exact_p"])
    tested["mann_whitney_q_bh"] = bh_adjust(tested["mann_whitney_p"])
    tested.to_csv(output / "regional_variant_support_all_tested.tsv", sep="\t", index=False)
    displayed = tested.loc[
        (tested["carrier_total"] >= args.minimum_carriers)
        & (tested["delta_median_log2_mic"] >= args.minimum_log2_mic_difference)
        & (tested["delta_resistant_fraction"] >= args.minimum_resistant_fraction_difference)
    ].copy()
    displayed = displayed.loc[
        ~displayed["region_id"].astype(str).eq(args.exclude_regional_region)
    ]
    displayed.to_csv(output / "regional_variant_support_displayed.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"parameter": "labelled_isolates_after_FungAMR_exclusion", "value": len(metadata)},
            {
                "parameter": "resistant_isolates_after_FungAMR_exclusion",
                "value": int(metadata[args.label_col].astype(str).eq("R").sum()),
            },
            {
                "parameter": "susceptible_isolates_after_FungAMR_exclusion",
                "value": int(metadata[args.label_col].astype(str).eq("S").sum()),
            },
            {"parameter": "minimum_carriers_for_display", "value": args.minimum_carriers},
            {
                "parameter": "minimum_delta_median_log2_mic_for_display",
                "value": args.minimum_log2_mic_difference,
            },
            {
                "parameter": "minimum_delta_resistant_fraction_for_display",
                "value": args.minimum_resistant_fraction_difference,
            },
            {"parameter": "tested_nonsynonymous_variants", "value": len(tested)},
            {"parameter": "displayed_variants", "value": len(displayed)},
        ]
    ).to_csv(output / "marker_support_run_parameters.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
