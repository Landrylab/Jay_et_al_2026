#!/usr/bin/env python3
"""Plot Manhattan and quantile–quantile summaries of additive log2-MIC GWAS.

The manifest is the primary_scan_manifest.tsv produced by
04_summarize_results.py. All eligible variants are plotted without rescaling
their P values. Relative result paths are resolved against the manifest folder.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import chi2


CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "9", "8"]
DRUG_ORDER = ["fluconazole", "caspofungin", "voriconazole", "micafungin",
              "posaconazole", "anidulafungin"]
STYLE = {"font.family": "DejaVu Sans", "font.size": 10,
         "font.weight": "bold", "axes.labelweight": "bold",
         "axes.titleweight": "bold", "pdf.fonttype": 42,
         "axes.spines.top": False, "axes.spines.right": False}


def valid_probabilities(values):
    p = np.asarray(values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
        raise ValueError("P values must be a nonempty vector of finite values in (0, 1]")
    return p


def qq_statistics(values):
    p = np.sort(valid_probabilities(values))
    ranks = np.arange(1, len(p) + 1, dtype=float)
    expected = -np.log10((ranks - 0.5) / len(p))
    observed = -np.log10(p)
    lambda_gc = float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1))
    return expected, observed, lambda_gc


def chromosome_label(value):
    text = str(value)
    if text in CHROMOSOMES:
        return text
    if text.lower() in {"r", "chrr"}:
        return "9"
    if text.lower() in {"m", "mt", "chrm", "chrmt", "ca22chrm_c_albicans_sc5314"}:
        return "8"
    match = re.fullmatch(r"(?:Ca22)?chr([1-7R])(?:A_C_albicans_SC5314)?", text, re.IGNORECASE)
    if match:
        return "9" if match.group(1).upper() == "R" else match.group(1)
    return None


def read_chromosome_lengths(fai):
    lengths = {}
    with Path(fai).open() as handle:
        for number, line in enumerate(handle, 1):
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed FASTA index at line {number}")
            chromosome = chromosome_label(parts[0])
            if chromosome is None:
                continue
            if chromosome in lengths:
                raise ValueError("FASTA index contains duplicate chromosome labels")
            length = int(parts[1])
            if length < 1:
                raise ValueError("Chromosome lengths must be positive")
            lengths[chromosome] = length
    if set(lengths) != set(CHROMOSOMES):
        raise ValueError("FASTA index must include chromosomes 1–7, R (9) and MT (8)")
    return lengths


def genomic_positions(frame, lengths):
    chromosomes = frame.Chr.map(chromosome_label)
    if chromosomes.isna().any():
        raise ValueError("Association table contains an unrecognized chromosome")
    bp = pd.to_numeric(frame.bp, errors="raise").to_numpy(dtype=float)
    upper = chromosomes.map(lengths).to_numpy(dtype=float)
    if not np.isfinite(bp).all() or not ((bp >= 1) & (bp <= upper) & (bp == np.floor(bp))).all():
        raise ValueError("Variant coordinates must be integer positions within the reference chromosomes")
    offsets, ticks = {}, []
    cumulative = 0
    for chromosome in CHROMOSOMES:
        offsets[chromosome] = cumulative
        ticks.append(cumulative + lengths[chromosome] / 2)
        cumulative += lengths[chromosome]
    positions = bp + chromosomes.map(offsets).to_numpy(dtype=float)
    return chromosomes, positions, ticks, cumulative


def validate_frame(frame, row):
    if not {"SNP", "Chr", "bp", "p"}.issubset(frame):
        raise ValueError("Association table requires SNP, Chr, bp and p columns")
    if frame.SNP.duplicated().any() or frame.SNP.isna().any():
        raise ValueError("Association table must contain distinct variant identifiers")
    p = valid_probabilities(frame.p)
    if len(p) != int(row["tested_variants"]):
        raise ValueError("Association count differs from the manifest")
    expected, observed, lambda_gc = qq_statistics(p)
    for name, value in [("per_scan_bonferroni_p", 0.05 / len(p)),
                        ("per_scan_suggestive_p", 1 / len(p)),
                        ("lambda_gc", lambda_gc)]:
        if not math.isclose(float(row[name]), value, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(f"Manifest {name} disagrees with the association table")
    return expected, observed, lambda_gc


def draw_qq(ax, expected, observed):
    xmax = max(1, math.ceil(float(expected.max())))
    ymax = max(xmax, math.ceil(float(observed.max()) * 1.05))
    ax.scatter(expected, observed, s=7, color="black", linewidths=0,
               rasterized=True, zorder=3)
    ax.plot([0, xmax], [0, xmax], color="black", linestyle="--", linewidth=1)
    ax.set(xlim=(0, xmax), ylim=(0, ymax), xlabel="Expected −log₁₀(P)",
           ylabel="Observed −log₁₀(P)")


def draw_manhattan(ax, frame, lengths):
    chromosomes, positions, ticks, total_length = genomic_positions(frame, lengths)
    p = valid_probabilities(frame.p)
    m = len(p)
    bonferroni, suggestive = 0.05 / m, 1 / m
    y = -np.log10(p)
    significant = p < bonferroni
    suggestive_only = (p < suggestive) & ~significant
    neutral = ~(significant | suggestive_only)
    colors = ["#527c97", "#a6a6a6"]
    for index, chromosome in enumerate(CHROMOSOMES):
        mask = (chromosomes == chromosome).to_numpy() & neutral
        ax.scatter(positions[mask], y[mask], s=5, color=colors[index % 2],
                   linewidths=0, rasterized=True, zorder=2)
    for mask, color in [(suggestive_only, "#e69f00"), (significant, "#cc397b")]:
        ax.scatter(positions[mask], y[mask], s=9, color=color, linewidths=0,
                   rasterized=True, zorder=3)
    ax.axhline(-math.log10(bonferroni), color="black", linestyle="--", linewidth=1)
    ax.axhline(-math.log10(suggestive), color="black", linestyle=":", linewidth=1)
    ymax = max(float(y.max()), -math.log10(bonferroni))
    ax.set(xlim=(0, total_length), ylim=(0, max(1, ymax * 1.10)),
           xlabel="Chromosome", ylabel="−log₁₀(P)")
    ax.set_xticks(ticks, ["1", "2", "3", "4", "5", "6", "7", "R", "MT"])


def load_scans(manifest):
    path = Path(manifest).expanduser().resolve()
    table = pd.read_csv(path, sep="\t")
    required = {"drug", "n", "tested_variants", "per_scan_bonferroni_p",
                "per_scan_suggestive_p", "lambda_gc", "result_path"}
    if not required.issubset(table) or table.empty or table.drug.duplicated().any():
        raise ValueError("Manifest requires one distinct row per drug and all summary columns")
    order = {name: i for i, name in enumerate(DRUG_ORDER)}
    rows = sorted(table.to_dict("records"), key=lambda row: (order.get(row["drug"], len(order)), row["drug"]))
    scans = []
    for row in rows:
        result = Path(row["result_path"]).expanduser()
        if not result.is_absolute():
            result = path.parent / result
        frame = pd.read_csv(result, sep="\t")
        expected, observed, lambda_gc = validate_frame(frame, row)
        scans.append((row, frame, expected, observed, lambda_gc))
    return scans


def render(manifest, outdir, fai=None, plots="both", dpi=300):
    scans = load_scans(manifest)
    include_manhattan = plots in {"both", "manhattan"}
    if include_manhattan and fai is None:
        raise ValueError("A FASTA index (--fai) is required for Manhattan plots")
    lengths = read_chromosome_lengths(fai) if include_manhattan else None
    if include_manhattan:
        for _, frame, *_ in scans:
            genomic_positions(frame, lengths)
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    ncols = min(2, len(scans))
    nrows = math.ceil(len(scans) / ncols)
    kinds = ["manhattan", "qq"] if plots == "both" else [plots]
    with plt.rc_context(STYLE):
        for kind in kinds:
            fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows), squeeze=False)
            for index, (row, frame, expected, observed, lambda_gc) in enumerate(scans):
                ax = axes.flat[index]
                if kind == "qq":
                    draw_qq(ax, expected, observed)
                    title = f"{row['drug'].capitalize()} · λ = {lambda_gc:.3f}"
                else:
                    draw_manhattan(ax, frame, lengths)
                    title = row["drug"].capitalize()
                ax.set_title(title + f"\nn = {int(row['n']):,}; variants = {len(frame):,}", fontsize=13, pad=13)
                ax.text(-0.11, 1.14, chr(ord("A") + index), transform=ax.transAxes,
                        fontsize=18, fontweight="bold", ha="left", va="bottom")
            for index in range(len(scans), nrows * ncols):
                axes.flat[index].set_visible(False)
            bottom = 0.05 if kind == "manhattan" else 0.01
            fig.tight_layout(rect=(0.01, bottom, 0.99, 0.99), h_pad=3.0, w_pad=2.5)
            if kind == "manhattan":
                handles = [Line2D([0], [0], color="black", linestyle="--", label="Bonferroni: P = 0.05/M"),
                           Line2D([0], [0], color="black", linestyle=":", label="Suggestive: P = 1/M")]
                fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False)
            for extension in ("png", "pdf"):
                fig.savefig(outdir / f"variant_gwas_{kind}.{extension}", dpi=dpi, facecolor="white")
            plt.close(fig)
    audit = {"expected_quantile_formula": "-log10((rank - 0.5) / M), ascending P ranks 1 through M",
             "lambda_gc_formula": "median(chi2.isf(P, 1)) / chi2.ppf(0.5, 1)",
             "lambda_scope": "all eligible variants per drug", "p_values_rescaled": False,
             "qq_confidence_bands": False, "all_points_plotted": True,
             "scans": [{"drug": row["drug"], "n": int(row["n"]),
                        "plotted_variants": len(frame), "lambda_gc": float(lam)}
                       for row, frame, _, _, lam in scans]}
    (outdir / "plot_summary.json").write_text(json.dumps(audit, indent=2) + "\n")
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fai", type=Path, help="Reference FASTA index for Manhattan coordinates")
    parser.add_argument("--plots", choices=["both", "manhattan", "qq"], default="both")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.dpi < 1:
        parser.error("DPI must be positive")
    try:
        audit = render(args.manifest, args.outdir, args.fai, args.plots, args.dpi)
    except (ValueError, OSError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
