#!/usr/bin/env python3
"""Combine chromosome-excluded gene tests and report carrier effects and thresholds."""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_ndtr


NULL_CHISQ_MEDIAN = 0.45493642311957275


def read_table(path):
    separator = "\t"
    if str(path).endswith(".mlma"):
        with Path(path).open() as handle:
            separator = "\t" if "\t" in handle.readline() else r"\s+"
    return pd.read_csv(path, sep=separator)


def resolved(value, manifest):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(manifest).resolve().parent / path


def summarize(raw, qc, drug, n, expected_count=None):
    required_raw = {"Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p", "loco_excluded_chromosome"}
    required_qc = {"gene_marker_id", "gene_id", "gene_symbol", "chromosome", "bp", "n_carriers", "n_noncarriers", "n_missing", "n_called"}
    if not required_raw.issubset(raw) or not required_qc.issubset(qc):
        raise ValueError("Missing association/QC columns")
    raw = raw.copy()
    qc = qc.copy()
    eligible = qc.n_carriers.gt(0) & qc.n_noncarriers.gt(0)
    qc = qc.loc[eligible].copy()
    if raw.empty or raw.SNP.duplicated().any() or qc.gene_marker_id.duplicated().any() or qc.gene_id.duplicated().any():
        raise ValueError("Each gene marker must have exactly one association and QC row")
    if set(raw.SNP) != set(qc.gene_marker_id):
        raise ValueError("Association results do not match all eligible gene markers")
    if expected_count is not None and len(raw) != int(expected_count):
        raise ValueError("Test count differs from the prepared manifest")
    raw = raw.rename(columns={"bp": "gcta_bp"})
    result = raw.merge(qc, left_on="SNP", right_on="gene_marker_id", validate="one_to_one")
    if not (result.Chr.eq(result.chromosome) & result.gcta_bp.eq(result.bp)
            & result.loco_excluded_chromosome.eq(result.chromosome)).all():
        raise ValueError("Gene coordinate or excluded chromosome mismatch")
    if not (result.A1.astype(str).eq("2") & result.A2.astype(str).eq("1")).all():
        raise ValueError("Expected carrier effect allele '2' and noncarrier allele '1'")
    if not (result.n_carriers.add(result.n_noncarriers).eq(result.n_called)
            & result.n_called.add(result.n_missing).eq(int(n))).all():
        raise ValueError("Carrier/missing counts disagree with the cohort size")
    frequency = result.n_carriers / result.n_called
    if (result.Freq - frequency).abs().max() >= 1e-6:
        raise ValueError("Reported allele frequency differs from carrier frequency")
    for field in ("b", "se", "p"):
        result[field] = pd.to_numeric(result[field], errors="raise")
        if not np.isfinite(result[field]).all():
            raise ValueError("Nonfinite association statistic: " + field)
    if not (result.se.gt(0) & result.p.between(0, 1)).all():
        raise ValueError("Invalid standard error or P value")
    z = result.b.to_numpy(float) / result.se.to_numpy(float)
    p = result.p.to_numpy(float)
    nonzero = p > 0
    neglog = np.empty_like(p)
    neglog[nonzero] = -np.log10(p[nonzero])
    neglog[~nonzero] = -(math.log(2) + log_ndtr(-np.abs(z[~nonzero]))) / math.log(10)
    m = len(result)
    result["neglog10P"] = neglog
    result["drug"] = drug
    result["phenotype"] = "log2_MIC"
    result["cohort_n"] = int(n)
    result["carrier_frequency_called"] = frequency
    result["per_drug_bonferroni_p"] = 0.05 / m
    result["per_drug_suggestive_p"] = 1.0 / m
    result["per_drug_bonferroni_significant"] = result.p.lt(0.05 / m)
    result["per_drug_suggestive_only"] = result.p.lt(1.0 / m) & ~result.per_drug_bonferroni_significant
    result["gcta_effect_units"] = "per encoded allele; carrier dosage 2, noncarrier dosage 0"
    result["carrier_vs_noncarrier_log2mic"] = 2 * result.b
    result["carrier_vs_noncarrier_se"] = 2 * result.se
    result["carrier_vs_noncarrier_fold_mic"] = np.exp2(2 * result.b)
    summary = {"drug": drug, "cohort_n": int(n), "tested_gene_markers": m,
               "per_drug_bonferroni_p": 0.05/m, "per_drug_suggestive_p": 1.0/m,
               "lambda_gc": float(np.median(z*z) / NULL_CHISQ_MEDIAN),
               "significant_genes": int(result.per_drug_bonferroni_significant.sum()),
               "suggestive_only_genes": int(result.per_drug_suggestive_only.sum())}
    return result.sort_values(["p", "chromosome", "bp", "gene_id"]), summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True,
                        help="association_manifest.tsv from the chromosome-specific GCTA runner")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    required = {"drug", "n", "tested_gene_markers", "result_path", "qc_path"}
    if manifest.empty or not required.issubset(manifest) or manifest.drug.duplicated().any():
        parser.error("Manifest must contain one row per drug and all required columns")
    if not manifest.drug.str.fullmatch(r"[A-Za-z0-9_-]+").all():
        parser.error("Drug names must contain letters, digits, underscores or hyphens")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    frames, summaries = [], []
    for row in manifest.to_dict("records"):
        result, summary = summarize(read_table(resolved(row["result_path"], args.manifest)),
                                    read_table(resolved(row["qc_path"], args.manifest)),
                                    row["drug"], row["n"], row["tested_gene_markers"])
        destination = outdir / (row["drug"] + ".tsv.gz")
        result.to_csv(destination, sep="\t", index=False, compression="gzip")
        summary["result_path"] = str(destination)
        summaries.append(summary)
        frames.append(result)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(outdir / "all_gene_level_variant_results.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(summaries).to_csv(outdir / "scan_manifest.tsv", sep="\t", index=False)
    combined.loc[combined.per_drug_bonferroni_significant | combined.per_drug_suggestive_only].to_csv(
        outdir / "significant_and_suggestive_genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
