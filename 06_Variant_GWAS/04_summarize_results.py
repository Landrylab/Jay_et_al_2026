#!/usr/bin/env python3
"""Validate variant association results and calculate drug-specific thresholds."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


def read_table(path):
    return pd.read_csv(path, sep=r"\s+" if str(path).endswith(".mlma") else "\t")


def resolve_input(value, manifest):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(manifest).resolve().parent / path


def summarize_frame(frame, drug, n, expected_count=None):
    frame = frame.copy()
    needed = {"SNP", "Chr", "bp", "A1", "A2", "b", "se", "p", "marker_type"}
    if not needed.issubset(frame):
        raise ValueError("Association table lacks columns: " + ", ".join(sorted(needed - set(frame))))
    if frame.empty or frame.SNP.duplicated().any():
        raise ValueError("Association table must contain distinct tested variants")
    if "allele_key" in frame and frame.allele_key.duplicated().any():
        raise ValueError("Association table contains duplicate biological alleles")
    if expected_count is not None and len(frame) != int(expected_count):
        raise ValueError("Result count differs from prepared eligible-variant count")
    for name in ("b", "se", "p"):
        frame[name] = pd.to_numeric(frame[name], errors="raise")
        if not np.isfinite(frame[name]).all():
            raise ValueError("Association output contains non-finite " + name)
    if not (frame.se.gt(0) & frame.p.gt(0) & frame.p.le(1)).all():
        raise ValueError("Association output contains invalid standard errors or P values")
    if not frame.marker_type.isin(["SNP", "INDEL"]).all():
        raise ValueError("marker_type must be SNP or INDEL")
    count = len(frame)
    bonferroni, suggestive = 0.05 / count, 1.0 / count
    frame["drug"] = drug
    frame["n"] = int(n)
    frame["bonferroni_p_per_scan"] = bonferroni
    frame["suggestive_p_per_scan"] = suggestive
    frame["per_scan_significant"] = frame.p.lt(bonferroni)
    frame["suggestive_only"] = frame.p.lt(suggestive) & ~frame.per_scan_significant
    frame["minus_log10_p"] = -np.log10(frame.p)
    lambda_gc = np.median(chi2.isf(frame.p.to_numpy(), df=1)) / chi2.ppf(0.5, df=1)
    summary = {"drug": drug, "n": int(n), "tested_variants": count,
               "tested_SNP": int(frame.marker_type.eq("SNP").sum()),
               "tested_INDEL": int(frame.marker_type.eq("INDEL").sum()),
               "per_scan_bonferroni_p": bonferroni, "per_scan_suggestive_p": suggestive,
               "per_scan_hits": int(frame.per_scan_significant.sum()),
               "suggestive_only_hits": int(frame.suggestive_only.sum()),
               "lambda_gc": float(lambda_gc)}
    return frame, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True,
                        help="TSV containing drug,n,result_path and optionally tested_variants,qc_path")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    if not {"drug", "n", "result_path"}.issubset(manifest) or manifest.empty or manifest.drug.duplicated().any():
        parser.error("Manifest must contain one distinct row per drug with n and result_path")
    if not manifest.drug.str.fullmatch(r"[A-Za-z0-9_-]+").all():
        parser.error("Drug names must contain letters, digits, underscores or hyphens")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    summaries, highlights = [], []
    for row in manifest.to_dict("records"):
        frame = read_table(resolve_input(row["result_path"], args.manifest))
        if pd.notna(row.get("qc_path")):
            qc = read_table(resolve_input(row["qc_path"], args.manifest))
            if qc.SNP.duplicated().any() or frame.SNP.duplicated().any() or set(frame.SNP) != set(qc.SNP):
                raise ValueError("QC and association tables must describe exactly the same distinct variants")
            indexed = qc.set_index("SNP").loc[frame.SNP]
            for column in ("Chr", "bp"):
                if not np.array_equal(frame[column].astype(str).to_numpy(), indexed[column].astype(str).to_numpy()):
                    raise ValueError("Association/QC allele or coordinate mismatch: " + column)
            same = (frame.A1.to_numpy() == indexed.A1.to_numpy()) & (frame.A2.to_numpy() == indexed.A2.to_numpy())
            swapped = (frame.A1.to_numpy() == indexed.A2.to_numpy()) & (frame.A2.to_numpy() == indexed.A1.to_numpy())
            if not (same | swapped).all():
                raise ValueError("Association/QC alleles do not match")
            # The association effect is expressed per output A1; source counts retain their own allele labels.
            qc = qc.rename(columns={name: "input_" + name for name in
                                   ("A1", "A2", "n_A1_hom", "n_A2_hom", "A1_allele_count") if name in qc})
            extra = [name for name in qc if name == "SNP" or name not in frame]
            frame = frame.merge(qc[extra], on="SNP", validate="one_to_one")
        frame, summary = summarize_frame(frame, row["drug"], row["n"], row.get("tested_variants"))
        dest = outdir / (row["drug"] + ".tsv.gz")
        frame.to_csv(dest, sep="\t", index=False, compression="gzip")
        summary["result_path"] = str(dest)
        summaries.append(summary)
        highlights.append(frame.loc[frame.per_scan_significant | frame.suggestive_only])
    pd.DataFrame(summaries).to_csv(outdir / "primary_scan_manifest.tsv", sep="\t", index=False)
    pd.concat(highlights, ignore_index=True).sort_values(["drug", "p"]).to_csv(
        outdir / "significant_and_suggestive_variants.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
