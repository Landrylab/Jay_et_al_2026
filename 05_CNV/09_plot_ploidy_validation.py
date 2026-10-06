#!/usr/bin/env python3
"""Plot allele-fraction density, genomic copy number, and FACS source tables."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from figure_sources import read_source_table


BLUE, GREY, RED = "#0072B2", "#7A7A7A", "#B22222"
TABLES = {"samples": "ED3Samples", "densities": "ED3Densities",
          "genome": "ED3GenomeBins", "facs": "ED3FACS", "composition": "ED3Composition"}


def require(frame, columns, name):
    missing = set(columns).difference(frame.columns)
    if missing or frame.empty:
        raise ValueError(f"{name}: expected nonempty table with columns {sorted(columns)}; missing {sorted(missing)}")


def numbers(frame, columns, nullable=()):
    for column in columns:
        raw = frame[column]
        absent = raw.isna() | raw.astype(str).str.strip().isin(["", "NA", "NaN"])
        result = pd.to_numeric(raw.where(~absent), errors="coerce")
        if ((~absent) & result.isna()).any() or (column not in nullable and result.isna().any()):
            raise ValueError(f"Invalid or missing numeric values in {column}")
        if np.isinf(result.to_numpy(float)).any():
            raise ValueError(f"Nonfinite values in {column}")
        frame[column] = result


def text_columns(frame, columns):
    for column in columns:
        if frame[column].isna().any() or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing identifier in {column}")
        frame[column] = frame[column].astype(str)


def validate_data(samples, density, genome, facs, composition=None, reference=None):
    samples, density, genome, facs = [x.copy() for x in (samples, density, genome, facs)]
    require(samples, ["Panel", "Isolate", "Assigned ploidy (N)"], "samples")
    require(density, ["Panel", "Allele fraction", "REF density", "ALT density"], "densities")
    require(genome, ["Panel", "Chromosome", "Bin start (bp)", "Bin end (bp)", "Absolute copy number", "Median minor AF"], "genome")
    require(facs, ["Sample", "Bin start (a.u.)", "Bin end (a.u.)", "Events", "Events (%)"], "FACS")
    text_columns(samples, ["Panel", "Isolate"])
    text_columns(density, ["Panel"])
    text_columns(genome, ["Panel", "Chromosome"])
    text_columns(facs, ["Sample"])
    numbers(samples, ["Assigned ploidy (N)"])
    if samples.Panel.duplicated().any() or samples.Isolate.duplicated().any():
        raise ValueError("Sample panels and isolates must be unique")
    ploidies = samples["Assigned ploidy (N)"]
    if ((ploidies < 1) | (ploidies % 1 != 0)).any():
        raise ValueError("Assigned ploidy must be a positive integer")
    panels = set(samples.Panel)
    if set(density.Panel) != panels or set(genome.Panel) != panels:
        raise ValueError("Density and genome tables must cover exactly the sample panels")
    numbers(density, ["Allele fraction", "REF density", "ALT density"])
    if not density["Allele fraction"].between(0, 1).all() or (density[["REF density", "ALT density"]] < 0).any().any():
        raise ValueError("Density coordinates must have AF in [0,1] and nonnegative densities")
    if density.duplicated(["Panel", "Allele fraction"]).any() or density.groupby("Panel").size().min() < 2:
        raise ValueError("Each density curve needs at least two unique AF coordinates")
    numbers(genome, ["Bin start (bp)", "Bin end (bp)", "Absolute copy number", "Median minor AF"], ["Median minor AF"])
    if ((genome["Bin start (bp)"] < 0) | (genome["Bin end (bp)"] <= genome["Bin start (bp)"]) | (genome["Absolute copy number"] < 0)).any():
        raise ValueError("Invalid genomic interval or copy number")
    if (genome[["Bin start (bp)", "Bin end (bp)"]] % 1 != 0).any().any():
        raise ValueError("Genomic coordinates must be integer base positions")
    if not genome["Median minor AF"].dropna().between(0, .5).all():
        raise ValueError("Minor AF must be in [0,0.5] or missing")
    if genome.duplicated(["Panel", "Chromosome", "Bin start (bp)"]).any():
        raise ValueError("Duplicated genomic bin")
    for _, group in genome.groupby(["Panel", "Chromosome"], sort=False):
        ordered = group.sort_values("Bin start (bp)")
        if np.any(ordered["Bin start (bp)"].to_numpy()[1:] < ordered["Bin end (bp)"].to_numpy()[:-1]):
            raise ValueError("Genomic bins overlap")
    numbers(facs, ["Bin start (a.u.)", "Bin end (a.u.)", "Events", "Events (%)"])
    if reference not in set(facs.Sample):
        raise ValueError(f"FACS reference sample {reference!r} is absent")
    if ((facs.Events < 0) | (facs.Events % 1 != 0)).any():
        raise ValueError("FACS events must be nonnegative integer counts")
    facs = facs.sort_values(["Sample", "Bin start (a.u.)"]).reset_index(drop=True)
    reference_bins = None
    event_counts = {}
    for sample, part in facs.groupby("Sample", sort=False):
        lo, hi = part["Bin start (a.u.)"].to_numpy(), part["Bin end (a.u.)"].to_numpy()
        if np.any(hi <= lo) or np.any(lo[1:] != hi[:-1]):
            raise ValueError(f"FACS bins for {sample} must be contiguous, ordered intervals")
        bins = np.column_stack([lo, hi])
        if reference_bins is not None and not np.array_equal(bins, reference_bins):
            raise ValueError("All FACS recordings must use the same bins")
        reference_bins = bins
        total = int(part.Events.sum())
        if not total or not np.allclose(part["Events (%)"], 100 * part.Events / total, rtol=1e-12, atol=1e-12):
            raise ValueError(f"FACS percentages for {sample} must use every recorded event")
        event_counts[sample] = total
    mapped = []
    for _, sample in samples.iterrows():
        if "FACS sample" in samples:
            value = sample["FACS sample"]
            match = None if pd.isna(value) or str(value).strip() == "" else str(value)
            if match is not None and match not in event_counts:
                raise ValueError(f"Mapped FACS sample {match!r} is absent")
        else:
            candidates = [sample.Isolate, sample.Isolate.split("_", 1)[0]]
            matches = list(dict.fromkeys(x for x in candidates if x in event_counts and x != reference))
            if len(matches) > 1:
                raise ValueError("Ambiguous FACS sample: supply the 'FACS sample' column")
            match = matches[0] if matches else None
        mapped.append(match)
    present = [x for x in mapped if x is not None]
    if len(set(present)) != len(present) or set(event_counts) - {reference} != set(present):
        raise ValueError("Each non-reference FACS recording must map to exactly one sample row")
    samples["FACS sample"] = mapped
    if composition is not None:
        composition = composition.copy()
        require(composition, ["Chromosome", "Aneuploidy type", "Copy number", "Isolates", "Chromosome denominator", "Aneuploid isolates (%)"], "composition")
        text_columns(composition, ["Chromosome", "Aneuploidy type"])
        numbers(composition, ["Copy number", "Isolates", "Chromosome denominator", "Aneuploid isolates (%)"])
        if not set(composition["Aneuploidy type"]).issubset({"whole chromosome", "segmental"}):
            raise ValueError("Aneuploidy type must be whole chromosome or segmental")
        integer_columns = ["Copy number", "Isolates", "Chromosome denominator"]
        if (composition[integer_columns] < 0).any().any() or (composition[integer_columns] % 1 != 0).any().any():
            raise ValueError("Composition counts and copy numbers must be nonnegative integers")
        if composition.duplicated(["Chromosome", "Aneuploidy type", "Copy number"]).any():
            raise ValueError("Duplicated composition category")
        for chrom, part in composition.groupby("Chromosome", sort=False):
            denominator = part["Chromosome denominator"].iloc[0]
            if denominator <= 0 or part["Chromosome denominator"].nunique() != 1 or part.Isolates.sum() != denominator:
                raise ValueError(f"Composition categories do not reconcile for chromosome {chrom}")
            if not np.allclose(part["Aneuploid isolates (%)"], 100 * part.Isolates / denominator, rtol=1e-12, atol=1e-12):
                raise ValueError(f"Incorrect composition percentages for chromosome {chrom}")
    return samples, density, genome, facs, composition, event_counts


def chromosome_coordinates(genome):
    """Use source chromosome order and maximum supplied bin ends."""
    lengths = genome.groupby("Chromosome", sort=False)["Bin end (bp)"].max()
    offsets = lengths.cumsum().shift(fill_value=0)
    return lengths, offsets


def plot(samples, density, genome, facs, composition, event_counts, output_prefix,
         reference=None, unavailable="Sample not available for testing", dpi=300,
         facs_min=None, facs_max=None, facs_ymax=None, cn_max=None, chromosome_layout="shared"):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6,
                         "axes.labelweight": "bold", "axes.linewidth": .55,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    if facs_min is None:
        facs_min = float(facs["Bin start (a.u.)"].min())
    if facs_max is None:
        facs_max = float(facs["Bin end (a.u.)"].max())
    if facs_ymax is None:
        facs_ymax = max(1., float(facs["Events (%)"].max()) * 1.05)
    if cn_max is None:
        cn_max = max(float(samples["Assigned ploidy (N)"].max()), float(genome["Absolute copy number"].max())) * 1.05
    if facs_min >= facs_max:
        raise ValueError("The fluorescence minimum must be below its maximum")
    common_lengths, common_offsets = chromosome_coordinates(genome)
    facs_legend_drawn = False
    n = len(samples)
    height = 22 + 58 * n + (55 if composition is not None else 0)
    fig = plt.figure(figsize=(180 / 25.4, height / 25.4))
    for x, heading in [(.17, "Allele-fraction distribution"), (.53, "Copy number and minor AF"), (.90, "FACS fluorescence")]:
        fig.text(x, .986, heading, ha="center", va="top", fontsize=6, weight="bold")
    grid = fig.add_gridspec(n + (composition is not None), 3, width_ratios=[1, 2.4, 1],
                          height_ratios=[1] * n + ([.95] if composition is not None else []),
                          left=.085, right=.98, top=.94, bottom=.07, hspace=.48, wspace=.55)
    summaries = []
    for row, (_, sample) in enumerate(samples.iterrows()):
        panel, ploidy = sample.Panel, int(sample["Assigned ploidy (N)"])
        d = density.loc[density.Panel.eq(panel)].sort_values("Allele fraction")
        ax = fig.add_subplot(grid[row, 0])
        ax.plot(d["Allele fraction"], d["ALT density"], color=RED, lw=.8, label="ALT")
        ax.plot(d["Allele fraction"], d["REF density"], color="#1E4E9A", lw=.8, label="REF")
        for x in np.arange(1, ploidy) / ploidy:
            ax.axvline(x, color="#777777", lw=.4, ls="--", alpha=.75)
        ax.set(xlim=(0, 1), xticks=[0, .5, 1], xlabel="Allele fraction", ylabel="Density")
        display = sample.get("Display label", sample.Isolate)
        if pd.isna(display) or not str(display).strip():
            display = sample.Isolate
        ax.set_title(f"{display} ({ploidy}N)", fontsize=6, pad=10, loc="left")
        ax.text(-.48, 1.1, panel, transform=ax.transAxes, fontsize=8, weight="bold")
        if row == 0:
            ax.legend(frameon=False, ncol=2, fontsize=5.5, loc="upper center", bbox_to_anchor=(.5, 1.24), handlelength=1.4)
        nested = grid[row, 1].subgridspec(2, 1, height_ratios=[1.8, 1], hspace=.2)
        cn_ax, af_ax = fig.add_subplot(nested[0]), fig.add_subplot(nested[1])
        g = genome.loc[genome.Panel.eq(panel)]
        lengths, offsets = (common_lengths, common_offsets) if chromosome_layout == "shared" else chromosome_coordinates(g)
        for i, (chrom, length) in enumerate(lengths.items()):
            values = g.loc[g.Chromosome.eq(chrom)]
            x = values["Bin start (bp)"] + offsets[chrom]
            cn_ax.scatter(x, values["Absolute copy number"], s=.45, color=BLUE if i % 2 == 0 else GREY, alpha=.72, linewidths=0)
            valid = values["Median minor AF"].notna()
            af_ax.scatter(x[valid], values.loc[valid, "Median minor AF"], s=.4, color="#1E4E9A", alpha=.72, linewidths=0)
            if i:
                for target in (cn_ax, af_ax):
                    target.axvline(offsets[chrom], ls=":", color="#333333", lw=.4)
        for target in (cn_ax, af_ax):
            target.set_xlim(-.003 * lengths.sum(), 1.003 * lengths.sum())
        cn_ax.set(ylim=(0, cn_max), ylabel="Copy number", xticks=[])
        cn_ax.axhline(ploidy, color="#333333", lw=.6)
        af_ax.set(ylim=(0, .5), ylabel="Minor AF", xlabel="Chromosome",
                  xticks=offsets.to_numpy() + lengths.to_numpy() / 2, xticklabels=lengths.index)
        af_ax.set_yticks([0, .25, .5])
        for y in np.arange(1, ploidy // 2 + 1) / ploidy:
            af_ax.axhline(y, color="#777777", ls="--", lw=.35, zorder=0)
        f_ax = fig.add_subplot(grid[row, 2])
        match = None if pd.isna(sample["FACS sample"]) else str(sample["FACS sample"])
        if match is None:
            f_ax.set_axis_off()
            words = unavailable.rsplit(" for ", 1)
            f_ax.text(.5, .5, "\nfor ".join(words), ha="center", va="center", transform=f_ax.transAxes, color="#666666", wrap=True)
        else:
            for name, color, linestyle in [(reference, "#555555", "--"), (match, BLUE, "-")]:
                part = facs.loc[facs.Sample.eq(name)]
                edges = np.r_[part["Bin start (a.u.)"].to_numpy(), part["Bin end (a.u.)"].iloc[-1]]
                f_ax.stairs(part["Events (%)"], edges, baseline=None, color=color, lw=.7, linestyle=linestyle,
                            label=reference if name == reference else "Isolate")
            f_ax.set(xlim=(facs_min, facs_max), ylim=(0, facs_ymax), xlabel="Fluorescence (a.u.)", ylabel="Events (%)")
            f_ax.xaxis.set_major_formatter(FuncFormatter(
                lambda value, position: f"{value / 1000:g}k" if abs(value) >= 1000 else f"{value:g}"))
            f_ax.text(.96, .95, f"n={event_counts[match]:,}", ha="right", va="top", transform=f_ax.transAxes, fontsize=5.5)
            if not facs_legend_drawn:
                f_ax.legend(frameon=False, fontsize=5.3, loc="upper center", bbox_to_anchor=(.5, 1.18), ncol=2, handlelength=1.5)
                facs_legend_drawn = True
        for target in (ax, cn_ax, af_ax, f_ax):
            target.spines[["top", "right"]].set_visible(False)
            target.tick_params(labelsize=5.5, length=2, width=.5)
        summaries.append({"panel": panel, "isolate": sample.Isolate, "ploidy": ploidy,
                          "density_points": len(d), "CN_bins": len(g),
                          "minor_AF_bins": int(g["Median minor AF"].notna().sum()),
                          "CN_bins_outside_axis": int((g["Absolute copy number"] > cn_max).sum()),
                          "FACS_sample": match, "FACS_events": event_counts.get(match)})
    if composition is not None:
        ax = fig.add_subplot(grid[n, :])
        chroms = composition.Chromosome.drop_duplicates().tolist()
        states = sorted(composition["Copy number"].unique())
        palette = {0: "#053061", 1: "#2166AC", 3: "#B2182B", 4: "#67001F"}
        for i, state in enumerate(states):
            palette.setdefault(state, plt.get_cmap("tab10")(i % 10))
        for j, chrom in enumerate(chroms):
            bottom = 0
            for kind in ["whole chromosome", "segmental"]:
                for state in states:
                    part = composition.loc[composition.Chromosome.eq(chrom) & composition["Aneuploidy type"].eq(kind) & composition["Copy number"].eq(state)]
                    if part.empty:
                        continue
                    item = part.iloc[0]
                    value = item["Aneuploid isolates (%)"]
                    ax.bar(j, value, bottom=bottom, width=.6, color=palette[state], edgecolor="white", lw=.5,
                           hatch="///" if kind == "segmental" else None)
                    if value >= 4:
                        ax.text(j, bottom + value / 2, str(int(item.Isolates)), color="white", ha="center", va="center", fontsize=5.5, weight="bold",
                                bbox={"facecolor": palette[state], "edgecolor": "none", "pad": .15})
                    bottom += value
            denominator = int(composition.loc[composition.Chromosome.eq(chrom), "Chromosome denominator"].iloc[0])
            ax.text(j, 102, f"n={denominator}", ha="center", fontsize=5.5)
        ax.set(ylim=(0, 106), xticks=range(len(chroms)), xticklabels=chroms,
               ylabel="Aneuploid isolates (%)", xlabel="Chromosome")
        ax.set_yticks(range(0, 101, 20))
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-.08, 1.12, chr(97 + n), transform=ax.transAxes, fontsize=8, weight="bold")
        handles = [Patch(facecolor=palette[x], label=f"CN{int(x)}") for x in states]
        handles += [Patch(facecolor="#737373", label="Whole chromosome"), Patch(facecolor="#737373", hatch="///", edgecolor="white", label="Segmental")]
        ax.legend(handles=handles, frameon=False, ncol=len(handles), loc="lower center", bbox_to_anchor=(.5, 1.15), fontsize=5.5)
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        fig.savefig(str(output_prefix) + "." + suffix, dpi=dpi, facecolor="white")
    plt.close(fig)
    report = {"samples": summaries, "FACS_reference": reference, "FACS_event_counts": event_counts,
              "FACS_normalization": "All histogram bins; no renormalization to displayed range",
              "FACS_display_limits": [facs_min, facs_max, 0, facs_ymax], "CN_display_limits": [0, cn_max],
              "chromosome_layout": chromosome_layout,
              "composition_rows": 0 if composition is None else len(composition),
              "dimensions_mm": [180, height]}
    Path(str(output_prefix) + ".summary.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Source workbook or directory of named TSV/CSV tables")
    for key, name in TABLES.items():
        parser.add_argument("--" + key, type=Path, help=f"Separate TSV/CSV for {name}; overrides --source")
    parser.add_argument("--no-composition", action="store_true")
    parser.add_argument("--reference-sample", required=True, help="Reference recording identifier in the FACS table")
    parser.add_argument("--unavailable-text", default="Sample not available for testing")
    parser.add_argument("--facs-map", type=Path, help="Optional TSV/CSV: Panel, FACS sample; blank means unavailable")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--facs-min", type=float, help="Minimum displayed fluorescence; default includes all supplied bins")
    parser.add_argument("--facs-max", type=float, help="Maximum displayed fluorescence; default includes all supplied bins")
    parser.add_argument("--facs-ymax", type=float, help="Maximum displayed events (%); default covers every bin")
    parser.add_argument("--cn-max", type=float, help="Maximum displayed CN; default covers every source value")
    parser.add_argument("--chromosome-layout", choices=["shared", "per-panel"], default="shared",
                        help="Use common chromosome offsets, or each panel's maximum bin ends")
    args = parser.parse_args()
    if args.dpi <= 0 or any(x is not None and (not np.isfinite(x) or x <= 0) for x in [args.facs_max, args.facs_ymax, args.cn_max]):
        parser.error("DPI and axis maxima must be positive")
    if args.facs_min is not None and (not np.isfinite(args.facs_min) or (args.facs_max is not None and args.facs_min >= args.facs_max)):
        parser.error("The fluorescence minimum must be finite and below its maximum")
    data = {}
    for key, name in TABLES.items():
        if key == "composition" and args.no_composition:
            data[key] = None
            continue
        source = getattr(args, key) or args.source
        if source is None:
            parser.error(f"Provide --{key} or --source")
        data[key] = read_source_table(source, name)
    if args.facs_map:
        mapping = read_source_table(args.facs_map, "ED3FACSMap")
        require(mapping, ["Panel", "FACS sample"], "FACS map")
        mapping["Panel"] = mapping.Panel.astype(str)
        if mapping.Panel.duplicated().any() or set(mapping.Panel) != set(data["samples"].Panel.astype(str)):
            raise ValueError("FACS map must have exactly one row per sample panel")
        data["samples"]["FACS sample"] = data["samples"].Panel.astype(str).map(mapping.set_index("Panel")["FACS sample"])
    checked = validate_data(data["samples"], data["densities"], data["genome"], data["facs"], data["composition"], args.reference_sample)
    report = plot(*checked, args.output_prefix, reference=args.reference_sample,
                  unavailable=args.unavailable_text, dpi=args.dpi, facs_min=args.facs_min, facs_max=args.facs_max,
                  facs_ymax=args.facs_ymax, cn_max=args.cn_max, chromosome_layout=args.chromosome_layout)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
