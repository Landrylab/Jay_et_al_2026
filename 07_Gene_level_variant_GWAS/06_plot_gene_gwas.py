#!/usr/bin/env python3
"""Plot gene-level Manhattan and QQ panels from aggregate association results.

Example:
  python 06_plot_gene_gwas.py --results /path/to/all_gene_level_variant_results.tsv.gz \
    --outdir /path/to/plots --gene-evidence /path/to/reviewed_gene_evidence.tsv

Evidence annotations require drug, gene_id, display_symbol and evidence_level.
Only direct/response annotations within each drug's current top 20 are labelled.
Without evidence annotations, highlights denote statistical thresholds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DRUG_ORDER = ["fluconazole", "caspofungin", "voriconazole", "micafungin", "posaconazole", "anidulafungin"]
CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "R", "MT"]
GAP = 360_000
CHI2_MEDIAN = 0.45493642311957275
BLUE, GREY = "#0072B2", "#999999"
STAT_COLORS = {"bonferroni": "#CC79A7", "suggestive": "#E69F00"}
EVIDENCE_COLORS = {"direct": "#8E145F", "response": "#D55E00"}


def chromosome_label(value):
    text = str(value)
    if text in CHROMOSOMES:
        return text
    if text in ("9", "9.0"):
        return "R"
    if text.lower() in {"8", "8.0", "m", "mt", "chrm", "chrmt", "ca22chrm_c_albicans_sc5314"}:
        return "MT"
    match = re.fullmatch(r"(?:Ca22)?chr([1-7R])(?:A_C_albicans_SC5314)?", text, flags=re.I)
    return match.group(1).upper() if match else None


def validate_results(data):
    required = {"gene_id", "gene_symbol", "chromosome", "bp", "p", "b", "se", "drug"}
    missing = required - set(data.columns)
    if missing or data.empty:
        raise ValueError("Results must be nonempty and include: " + ", ".join(sorted(required)))
    out = data.copy()
    if out[["drug", "gene_id"]].isna().any().any() or out.duplicated(["drug", "gene_id"]).any():
        raise ValueError("Each drug/gene_id must identify one nonmissing association")
    for column in ("bp", "p", "b", "se"):
        out[column] = pd.to_numeric(out[column], errors="raise")
        if not np.isfinite(out[column]).all():
            raise ValueError("Nonfinite result values in " + column)
    if not out.p.between(0, 1, inclusive="right").all():
        raise ValueError("P values must satisfy 0 < P <= 1")
    if (out.se <= 0).any() or (out.bp < 1).any() or (out.bp % 1 != 0).any():
        raise ValueError("SE must be positive and genomic coordinates positive integers")
    out["plot_chromosome"] = out.chromosome.map(chromosome_label)
    if out.plot_chromosome.isna().any():
        raise ValueError("Unrecognized chromosome in results")
    out["gene_symbol"] = out.gene_symbol.fillna(out.gene_id)
    calculated = -np.log10(out.p.to_numpy(float))
    if "neglog10P" in out and not np.allclose(out.neglog10P, calculated, rtol=1e-8, atol=1e-10):
        raise ValueError("neglog10P disagrees with raw P values")
    out["neglog10P"] = calculated
    count = out.groupby("drug").gene_id.transform("size")
    expected_bonf = 0.05 / count
    if "per_drug_bonferroni_p" in out and not np.allclose(out.per_drug_bonferroni_p, expected_bonf, rtol=1e-8, atol=1e-15):
        raise ValueError("Bonferroni thresholds disagree with tested gene counts")
    out["bonferroni_p"] = expected_bonf
    out["suggestive_p"] = 1.0 / count
    out["bonferroni_significant"] = out.p < out.bonferroni_p
    out["suggestive_only"] = (out.p < out.suggestive_p) & ~out.bonferroni_significant
    out = out.sort_values(["drug", "p", "gene_id"]).reset_index(drop=True)
    out["rank"] = out.groupby("drug").cumcount() + 1
    return out


def qq_statistics(p, effect, standard_error):
    p, effect, standard_error = [np.asarray(x, dtype=float) for x in (p, effect, standard_error)]
    if p.ndim != 1 or not p.size or p.shape != effect.shape or p.shape != standard_error.shape:
        raise ValueError("QQ inputs must be nonempty vectors of equal length")
    if not np.isfinite(p).all() or np.any((p <= 0) | (p > 1)):
        raise ValueError("QQ probabilities must satisfy 0 < P <= 1")
    if not np.isfinite(effect).all() or not np.isfinite(standard_error).all() or np.any(standard_error <= 0):
        raise ValueError("QQ effects must be finite and standard errors positive")
    expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    observed = -np.log10(np.sort(p))
    lambda_gc = float(np.median(np.square(effect / standard_error)) / CHI2_MEDIAN)
    return expected, observed, lambda_gc


def genomic_positions(data, fai=None):
    lengths = {}
    if fai is not None:
        with Path(fai).open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError("FASTA-index rows require a name and length")
                chrom = chromosome_label(fields[0])
                if chrom is None:
                    continue
                if chrom in lengths or int(fields[1]) <= 0:
                    raise ValueError("Reference chromosomes must be unique with positive lengths")
                lengths[chrom] = int(fields[1])
    else:
        lengths = data.groupby("plot_chromosome").bp.max().astype(int).to_dict()
    observed = set(data.plot_chromosome)
    if not observed <= set(lengths):
        raise ValueError("Reference lengths are missing a result chromosome")
    if (data.bp > data.plot_chromosome.map(lengths)).any():
        raise ValueError("Gene coordinates extend beyond reference chromosome lengths")
    order = [chrom for chrom in CHROMOSOMES if chrom in lengths]
    offsets, centers, cursor = {}, {}, 0
    for chrom in order:
        offsets[chrom] = cursor
        centers[chrom] = cursor + lengths[chrom] / 2
        cursor += lengths[chrom] + GAP
    out = data.copy()
    out["plot_bp"] = out.bp + out.plot_chromosome.map(offsets)
    return out, centers, cursor - GAP


def apply_highlights(data, evidence=None):
    out = data.copy()
    out["display_symbol"] = out.gene_symbol
    out["highlight"] = ""
    if evidence is None:
        out.loc[out.suggestive_only, "highlight"] = "suggestive"
        out.loc[out.bonferroni_significant, "highlight"] = "bonferroni"
        return out
    required = {"drug", "gene_id", "display_symbol", "evidence_level"}
    if not required <= set(evidence.columns) or evidence.duplicated(["drug", "gene_id"]).any():
        raise ValueError("Evidence requires unique drug/gene_id rows with display_symbol and evidence_level")
    if evidence[list(required)].isna().any().any() or not evidence.evidence_level.isin(["direct", "response", "pathway", "none"]).all():
        raise ValueError("Evidence fields must be nonmissing and use direct/response/pathway/none levels")
    lookup = evidence.set_index(["drug", "gene_id"])
    for index, row in out.loc[out["rank"] <= 20].iterrows():
        key = (row.drug, row.gene_id)
        if key in lookup.index:
            annotation = lookup.loc[key]
            out.at[index, "display_symbol"] = annotation.display_symbol
            if annotation.evidence_level in EVIDENCE_COLORS:
                out.at[index, "highlight"] = annotation.evidence_level
    return out


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)


def place_labels(ax, hits, total_length, ymax, thresholds, colors):
    occupied = []
    xmin, xmax = ax.get_xlim()
    span = xmax - xmin
    renderer = ax.figure.canvas.get_renderer()
    font = FontProperties(size=9.5, weight="bold", style="italic")
    placements = []
    for row in hits.sort_values(["p", "gene_id"]).itertuples(index=False):
        width, height, descent = renderer.get_text_width_height_descent(str(row.display_symbol), font, ismath=False)
        half_width = min(0.45 * span, 0.5 * (width + 10) * span / ax.bbox.width)
        half_height = max(0.018 * ymax, 0.5 * (height + descent + 6) * ymax / ax.bbox.height)
        label_x = min(max(float(row.plot_bp), xmin + half_width), xmax - half_width)
        label_y = float(row.neglog10P) + 0.045 * ymax
        for _ in range(200):
            rectangle = (label_x - half_width, label_x + half_width, label_y - half_height, label_y + half_height)
            collision = any(rectangle[0] < old[1] and rectangle[1] > old[0] and rectangle[2] < old[3] and rectangle[3] > old[2] for old in occupied)
            collision |= any(rectangle[2] - 0.01 * ymax < y < rectangle[3] + 0.01 * ymax for y in thresholds)
            if not collision:
                break
            label_y += 0.045 * ymax
        occupied.append(rectangle)
        placements.append((row, label_x, label_y))
    for row, x, y in placements:
        ax.annotate(str(row.display_symbol), xy=(row.plot_bp, row.neglog10P), xytext=(x, y),
                    ha="center", va="center", fontsize=9.5, fontstyle="italic", fontweight="bold",
                    color=colors[row.highlight] if row.highlight in EVIDENCE_COLORS else "#222222",
                    arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.5},
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.1})
    return max([y + 0.04 * ymax for _, _, y in placements] + [ymax])


def render(data, outdir, evidence=None, fai=None, dpi=300):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = apply_highlights(validate_results(data), evidence)
    data, centers, total_length = genomic_positions(data, fai)
    drugs = [drug for drug in DRUG_ORDER if drug in set(data.drug)]
    drugs += sorted(set(data.drug) - set(drugs))
    ncols, nrows = min(2, len(drugs)), math.ceil(len(drugs) / 2)
    mpl.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "bold", "axes.labelweight": "bold", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(nrows, ncols, figsize=(10 * ncols, 5.5 * nrows + 1), squeeze=False)
    qq_fig, qq_axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5 * nrows + 0.5), squeeze=False)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.955, bottom=0.105, hspace=0.38, wspace=0.2)
    qq_fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.08, hspace=0.3, wspace=0.28)
    summaries, qq_rows = [], []
    colors = EVIDENCE_COLORS if evidence is not None else STAT_COLORS
    for drug, ax, qax in zip(drugs, axes.flat, qq_axes.flat):
        part = data.loc[data.drug == drug]
        count = len(part)
        for index, chrom in enumerate(centers):
            chrom_data = part.loc[part.plot_chromosome == chrom]
            ax.scatter(chrom_data.plot_bp, chrom_data.neglog10P, s=15, c=BLUE if index % 2 == 0 else GREY, alpha=0.8, linewidths=0, rasterized=True)
        bonf, suggestive = -np.log10(0.05 / count), np.log10(count)
        ax.axhline(bonf, color="#333333", lw=1, linestyle=(0, (5, 3)))
        ax.axhline(suggestive, color="#666666", lw=1, linestyle=(0, (1, 3)))
        for level, color in colors.items():
            selected = part.loc[part.highlight == level]
            ax.scatter(selected.plot_bp, selected.neglog10P, s=50, c=color, edgecolor="#263238", linewidth=0.6, zorder=5)
        ax.set_xlim(-0.025 * total_length, 1.025 * total_length)
        ymax = max(bonf + 1, float(part.neglog10P.max()) * 1.18)
        hits = part.loc[part.highlight != ""]
        ymax = place_labels(ax, hits, total_length, ymax, [bonf, suggestive], colors)
        ax.set_ylim(0, ymax)
        ax.set_xticks(list(centers.values()), list(centers))
        ax.set_title(drug.capitalize(), fontsize=19, pad=12)
        ax.set_xlabel("Chromosome", fontsize=13)
        ax.set_ylabel("−log₁₀(P)", fontsize=13)
        clean(ax)

        ordered = part.sort_values(["p", "gene_id"])
        expected, observed, inflation = qq_statistics(ordered.p, ordered.b, ordered.se)
        qax.scatter(expected, observed, s=12, color=BLUE, alpha=0.8, linewidths=0, rasterized=True)
        significant = ordered.bonferroni_significant.to_numpy()
        qax.scatter(expected[significant], observed[significant], s=28, color=STAT_COLORS["bonferroni"], edgecolor="#263238", linewidth=0.6, zorder=5)
        limit = max(1, float(expected.max()) * 1.04)
        qax.plot([0, limit], [0, limit], linestyle="--", color="#333333", lw=1)
        qax.set_xlim(0, limit)
        qax.set_ylim(0, max(limit, float(observed.max()) * 1.07))
        qax.set_title(f"{drug.capitalize()} [λ = {inflation:.3f}]", fontsize=17, pad=12)
        qax.set_xlabel("Expected −log₁₀(P)", fontsize=12)
        qax.set_ylabel("Observed −log₁₀(P)", fontsize=12)
        clean(qax)
        qq_rows.append(pd.DataFrame({"drug": drug, "gene_id": ordered.gene_id.to_numpy(), "p": ordered.p.to_numpy(), "rank": np.arange(1, count + 1), "expected_neglog10p": expected, "observed_neglog10p": observed}))
        summaries.append({"drug": drug, "tested_genes": count, "bonferroni_p": 0.05 / count, "suggestive_p": 1 / count,
                          "bonferroni_significant_genes": int(significant.sum()), "suggestive_only_genes": int(part.suggestive_only.sum()),
                          "labelled_genes": len(hits), "lambda_gc": inflation,
                          "highlight_mode": "reviewed_top20_evidence" if evidence is not None else "statistical_thresholds"})
    for axis in list(axes.flat)[len(drugs):] + list(qq_axes.flat)[len(drugs):]:
        axis.set_visible(False)
    labels = {"direct": "Published azole/echinocandin phenotype or target evidence", "response": "Other drug, cell-wall or indirect response evidence"} if evidence is not None else {"bonferroni": "Bonferroni-significant gene", "suggestive": "Suggestive gene (not Bonferroni-significant)"}
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=8, label=labels[level]) for level, color in colors.items()]
    handles += [Line2D([0], [0], color="#333333", linestyle="--", label="Bonferroni threshold (0.05/M)"), Line2D([0], [0], color="#666666", linestyle=":", label="Suggestive threshold (1/M)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=11, bbox_to_anchor=(0.5, 0.015))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.955, bottom=0.105, hspace=0.38, wspace=0.2)
    qq_fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=STAT_COLORS["bonferroni"], markersize=7, label="Bonferroni-significant gene")], loc="lower center", frameon=False, fontsize=11)
    qq_fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.08, hspace=0.3, wspace=0.28)
    for figure, stem in [(fig, "gene_association_manhattan"), (qq_fig, "gene_association_qq")]:
        for suffix in ("pdf", "png"):
            figure.savefig(outdir / f"{stem}.{suffix}", dpi=dpi, facecolor="white")
        plt.close(figure)
    pd.DataFrame(summaries).to_csv(outdir / "plot_summary.tsv", sep="\t", index=False)
    data.loc[data.highlight != "", ["drug", "gene_id", "display_symbol", "p", "rank", "highlight"]].to_csv(outdir / "labelled_genes.tsv", sep="\t", index=False)
    pd.concat(qq_rows, ignore_index=True).to_csv(outdir / "qq_data.tsv.gz", sep="\t", index=False)
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--gene-evidence", type=Path)
    parser.add_argument("--fai", type=Path, help="Optional chromosome lengths; otherwise maximum observed gene coordinates are used")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.dpi < 1:
        parser.error("--dpi must be positive")
    try:
        evidence = pd.read_csv(args.gene_evidence, sep="\t") if args.gene_evidence else None
        summaries = render(pd.read_csv(args.results, sep="\t"), args.outdir, evidence, args.fai, args.dpi)
        source_paths = [path for path in (args.results, args.gene_evidence, args.fai) if path is not None]
        provenance = {"inputs": [{"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
                      "lambda_definition": "median((b/se)^2)/0.45493642311957275",
                      "qq_expected_probabilities": "(rank-0.5)/M", "chromosome_gap_bp": GAP, "panels": summaries}
        (args.outdir / "plot_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    except (ValueError, OSError) as error:
        parser.exit(2, str(error) + "\n")


if __name__ == "__main__":
    main()
