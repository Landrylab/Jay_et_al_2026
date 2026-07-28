#!/usr/bin/env python3
# =============================================================================
# Description
#   Render a separate set of final GWAS paper figures with ordinary Manhattan
#   points in black and grey points marking GWAS SNPs in LD with high-confidence
#   FungAMR mutations, including FungAMR mutations below the GWAS MAF threshold.
#
# Inputs
#   GCTA MLMA-LOCO outputs, filtered GWAS BIM files, unfiltered joint PLINK
#   files, and FungAMR score 1/2 nonsynonymous mapping records.
#
# Outputs
#   final_gwas_paper_figures_lowmaf_fungamr_ld_black_points/
# =============================================================================

from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[2]
GWAS_RESULTS_ROOT = Path(
    os.environ.get("GWAS_RESULTS_ROOT", ROOT / "gcta_mlma_loco")
)
BASE_PLINK = Path(
    os.environ.get(
        "GWAS_BASE_PLINK",
        ROOT / "gcta_mlma_loco" / "plink_from_joint_vcf" / "joint_keep1153",
    )
)
FUNGAMR_MAPPING = Path(
    os.environ.get(
        "FUNGAMR_MAPPING",
        ROOT / "reference_data" / "fungamr_score_1_or_2_nonsynonymous_matches.csv",
    )
)
PLOT_SCRIPT = ROOT / "scripts_for_git" / "GWAS" / "helpers" / "gwas_plot_mlma_results.R"
PLINK_BIN = Path(os.environ.get("PLINK_BIN", "plink"))

OUTDIR = Path(
    os.environ.get(
        "GWAS_FIGURE_OUTDIR",
        ROOT / "final_gwas_paper_figures_lowmaf_fungamr_ld_black_points",
    )
)
PER_TRAIT_ROOT = OUTDIR / "per_trait_plots"
LD_ROOT = OUTDIR / "fungamr_ld_grey_snps"

R2_THRESHOLD = 0.2
WINDOW_KB = 1000
DPI = 300

RUNS = [
    ("maf0p05", "MAF 0.05", GWAS_RESULTS_ROOT / "gcta_mlma_loco_maf005_six_antifungals"),
    ("maf0p005", "MAF 0.005", GWAS_RESULTS_ROOT / "gcta_mlma_loco_maf0005_six_antifungals"),
]

DRUGS = [
    "fluconazole",
    "voriconazole",
    "posaconazole",
    "caspofungin",
    "micafungin",
    "anidulafungin",
]
AZOLE_DRUGS = ["fluconazole", "posaconazole", "voriconazole"]
ECHINOCANDIN_DRUGS = ["caspofungin", "micafungin", "anidulafungin"]
MODES = {"mic": "MIC", "rs_binary": "R/S binary"}
PLOT_TYPES = ["manhattan", "qq"]


def score_is_high_confidence(value: str) -> bool:
    try:
        return float(value) in {1.0, 2.0}
    except ValueError:
        return False


def normalize_drug(value: str) -> str:
    return value.strip().lower()


def normalize_gene(value: str) -> str:
    return value.strip().upper()


def locus_from_mapping_row(row: dict[str, str]) -> str:
    chrom = (row.get("chrom") or "").strip()
    if chrom.startswith("Ca22chr"):
        chrom_label = chrom
    elif chrom.upper() in {"M", "MT"}:
        chrom_label = "Ca22chrMA_C_albicans_SC5314"
    else:
        chrom_label = f"Ca22chr{chrom}A_C_albicans_SC5314"
    return ":".join(
        [
            chrom_label,
            (row.get("pos") or "").strip(),
            (row.get("ref") or "").strip(),
            (row.get("alt") or "").strip(),
        ]
    )


def load_fungamr_records() -> tuple[list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    records: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    with FUNGAMR_MAPPING.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            drug = normalize_drug(row.get("antifungal", ""))
            score = (row.get("confidence_score") or "").strip()
            if drug not in DRUGS or not score_is_high_confidence(score):
                continue
            gene = normalize_gene(row.get("gene") or row.get("fungamr_gene") or "")
            mutation = (row.get("mutation_found") or "").strip()
            locus_id = locus_from_mapping_row(row)
            key = (drug, locus_id, gene, mutation, score)
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "drug": drug,
                    "fungamr_record": f"{gene} {mutation}".strip(),
                    "fungamr_gene": gene,
                    "fungamr_substitution": mutation,
                    "fungamr_score": score,
                    "locus_id": locus_id,
                }
            )

    by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for rec in records:
        by_locus[rec["locus_id"]].append(rec)
    return records, by_locus


def read_bim_records(prefix: Path) -> dict[str, tuple[str, int]]:
    records: dict[str, tuple[str, int]] = {}
    with Path(str(prefix) + ".bim").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                bp = int(float(fields[3]))
            except ValueError:
                continue
            records[fields[1]] = (fields[0], bp)
    return records


def filtered_prefix(run_dir: Path) -> Path:
    bims = sorted(run_dir.glob("01_*_filtered_snps/*.bim"))
    if len(bims) != 1:
        raise RuntimeError(f"Expected one filtered BIM under {run_dir}, found {len(bims)}")
    return bims[0].with_suffix("")


def write_lines(path: Path, values: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def update_best_detail(
    details: dict[tuple[str, str, str, str], dict[str, str]],
    snp_id: str,
    rec: dict[str, str],
    r2: float,
) -> None:
    key = (rec["drug"], snp_id, rec["fungamr_record"], rec["locus_id"])
    old = details.get(key)
    if old is not None and float(old["r2"]) >= r2:
        return
    details[key] = {
        "drug": rec["drug"],
        "snp": snp_id,
        "fungamr_record": rec["fungamr_record"],
        "fungamr_gene": rec["fungamr_gene"],
        "fungamr_substitution": rec["fungamr_substitution"],
        "fungamr_score": rec["fungamr_score"],
        "fungamr_locus": rec["locus_id"],
        "r2": f"{r2:.6g}",
    }


def compute_ld_details(
    filtered_records: dict[str, tuple[str, int]],
    base_records: dict[str, tuple[str, int]],
    records_by_locus: dict[str, list[dict[str, str]]],
) -> dict[tuple[str, str, str, str], dict[str, str]]:
    details: dict[tuple[str, str, str, str], dict[str, str]] = {}
    filtered_ids = set(filtered_records)
    available_fungamr_loci = sorted(set(records_by_locus) & set(base_records))
    if not available_fungamr_loci:
        return details

    for locus in available_fungamr_loci:
        if locus in filtered_ids:
            for rec in records_by_locus[locus]:
                update_best_detail(details, locus, rec, 1.0)

    window_bp = WINDOW_KB * 1000
    target_by_chr: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for locus in available_fungamr_loci:
        chrom, bp = base_records[locus]
        target_by_chr[chrom].append((locus, bp))

    nearby_filtered_ids: set[str] = set()
    for snp_id, (chrom, bp) in filtered_records.items():
        for _locus, target_bp in target_by_chr.get(chrom, []):
            if abs(bp - target_bp) <= window_bp:
                nearby_filtered_ids.add(snp_id)
                break

    available_extract_ids = sorted((nearby_filtered_ids | set(available_fungamr_loci)) & set(base_records))
    if not available_extract_ids:
        return details

    with tempfile.TemporaryDirectory(prefix="fungamr_lowmaf_ld_", dir="/private/tmp") as tmp:
        tmpdir = Path(tmp)
        extract_path = tmpdir / "extract_ids.txt"
        fungamr_path = tmpdir / "fungamr_loci.txt"
        write_lines(extract_path, available_extract_ids)
        write_lines(fungamr_path, available_fungamr_loci)

        out_prefix = tmpdir / "fungamr_ld"
        subprocess.run(
            [
                str(PLINK_BIN),
                "--bfile",
                str(BASE_PLINK),
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
                str(WINDOW_KB),
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
            return details

        with ld_path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
            for row in reader:
                snp_a = row.get("SNP_A", "")
                snp_b = row.get("SNP_B", "")
                try:
                    r2 = float(row.get("R2", "nan"))
                except ValueError:
                    continue
                if not math.isfinite(r2) or r2 < R2_THRESHOLD:
                    continue
                for fung_locus, plotted_snp in ((snp_a, snp_b), (snp_b, snp_a)):
                    if fung_locus not in records_by_locus or plotted_snp not in filtered_ids:
                        continue
                    for rec in records_by_locus[fung_locus]:
                        update_best_detail(details, plotted_snp, rec, r2)

    return details


def render_trait(run_key: str, run_dir: Path, drug: str, mode: str, grey_file: Path) -> None:
    source_outdir = run_dir / "outputs" / drug / mode
    mlma = source_outdir / "gwas.loco.mlma"
    source_plotdir = source_outdir / "plots"
    out_plotdir = PER_TRAIT_ROOT / run_key / drug / mode / "plots"
    label_file = source_plotdir / f"manuscript_candidate_labels.{drug}.{mode}.tsv"
    if not mlma.exists():
        raise FileNotFoundError(mlma)

    cmd = [
        "Rscript",
        "--vanilla",
        str(PLOT_SCRIPT),
        "--mlma",
        str(mlma),
        "--outdir",
        str(out_plotdir),
        "--trait",
        drug,
        "--mode",
        mode,
        "--grey-snps",
        str(grey_file),
        "--grey-label",
        "In LD with high-confidence FungAMR mutation",
        "--monochrome-points",
    ]
    if run_key == "maf0p05":
        cmd.extend(["--manhattan-ymax", "22"])
    if label_file.exists():
        cmd.extend(["--label-hits", str(label_file)])
    subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def plot_path(run_key: str, drug: str, mode: str, plot_type: str) -> Path:
    return (
        PER_TRAIT_ROOT
        / run_key
        / drug
        / mode
        / "plots"
        / f"{plot_type}.{drug}.{mode}.png"
    )


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def white_canvas(width: int, height: int) -> Image.Image:
    return Image.new("RGB", (width, height), "white")


def read_rgb(path: Path) -> Image.Image:
    img = Image.open(path).convert("RGBA")
    canvas = white_canvas(img.width, img.height)
    canvas.paste(img, mask=img.getchannel("A"))
    return canvas


def add_panel_letter(img: Image.Image, letter: str) -> None:
    draw = ImageDraw.Draw(img)
    draw.text((45, 35), letter, fill="black", font=load_font(96))


def extent_center_width(img: Image.Image, target_width: int) -> Image.Image:
    if img.width == target_width:
        return img
    canvas = white_canvas(target_width, img.height)
    canvas.paste(img, ((target_width - img.width) // 2, 0))
    return canvas


def extent_top_height(img: Image.Image, target_height: int) -> Image.Image:
    if img.height == target_height:
        return img
    canvas = white_canvas(img.width, target_height)
    canvas.paste(img, (0, 0))
    return canvas


def stack_plots(paths: list[Path], panel_letter: str) -> Image.Image:
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("\n".join(str(path) for path in missing))
    imgs = [read_rgb(path) for path in paths]
    add_panel_letter(imgs[0], panel_letter)
    target_width = max(img.width for img in imgs)
    imgs = [extent_center_width(img, target_width) for img in imgs]
    total_height = sum(img.height for img in imgs)
    canvas = white_canvas(target_width, total_height)
    y = 0
    for img in imgs:
        canvas.paste(img, (0, y))
        y += img.height
    return canvas


def compose_figure(run_key: str, maf_label: str, mode: str, mode_label: str, plot_type: str) -> dict[str, str]:
    left_paths = [plot_path(run_key, drug, mode, plot_type) for drug in AZOLE_DRUGS]
    right_paths = [plot_path(run_key, drug, mode, plot_type) for drug in ECHINOCANDIN_DRUGS]
    left = stack_plots(left_paths, "A")
    right = stack_plots(right_paths, "B")
    target_height = max(left.height, right.height)
    left = extent_top_height(left, target_height)
    right = extent_top_height(right, target_height)
    spacer_width = 90
    composite = white_canvas(left.width + spacer_width + right.width, target_height)
    composite.paste(left, (0, 0))
    composite.paste(right, (left.width + spacer_width, 0))

    out_name = f"gwas_{plot_type}_{run_key}_{mode}_panel_A_azoles_panel_B_echinocandins.png"
    out_path = OUTDIR / out_name
    composite.save(out_path, dpi=(DPI, DPI))
    return {
        "figure": out_name,
        "plot_type": plot_type,
        "maf": maf_label,
        "mode": mode_label,
        "output": str(out_path),
        "panel_A_sources": ";".join(str(path) for path in left_paths),
        "panel_B_sources": ";".join(str(path) for path in right_paths),
    }


def main() -> None:
    if not PLINK_BIN.exists():
        raise FileNotFoundError(PLINK_BIN)
    if not Path(str(BASE_PLINK) + ".bim").exists():
        raise FileNotFoundError(Path(str(BASE_PLINK) + ".bim"))
    if not FUNGAMR_MAPPING.exists():
        raise FileNotFoundError(FUNGAMR_MAPPING)
    if not PLOT_SCRIPT.exists():
        raise FileNotFoundError(PLOT_SCRIPT)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fungamr_records, records_by_locus = load_fungamr_records()
    base_records = read_bim_records(BASE_PLINK)
    available_loci = sorted(set(records_by_locus) & set(base_records))
    print(
        f"Loaded {len(fungamr_records)} high-confidence FungAMR mutation records; "
        f"{len(available_loci)} loci are present in the unfiltered PLINK set.",
        flush=True,
    )

    detail_fields = [
        "drug",
        "snp",
        "fungamr_record",
        "fungamr_gene",
        "fungamr_substitution",
        "fungamr_score",
        "fungamr_locus",
        "r2",
    ]
    r2_label = f"r2_ge{str(R2_THRESHOLD).replace('.', 'p')}"
    ld_manifest: list[dict[str, str]] = []

    for run_key, _maf_label, run_dir in RUNS:
        prefix = filtered_prefix(run_dir)
        print(f"[{run_key}] computing LD using filtered GWAS SNPs from {prefix}.bim", flush=True)
        filtered_records = read_bim_records(prefix)
        details = compute_ld_details(filtered_records, base_records, records_by_locus)

        detail_by_drug: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in details.values():
            detail_by_drug[row["drug"]].append(row)

        for drug in DRUGS:
            drug_details = sorted(
                detail_by_drug.get(drug, []),
                key=lambda row: (row["snp"], row["fungamr_record"], row["fungamr_locus"]),
            )
            grey_snps = sorted({row["snp"] for row in drug_details})
            grey_file = LD_ROOT / run_key / f"{drug}.ld_{r2_label}.window{WINDOW_KB}kb.snps.txt"
            detail_file = LD_ROOT / run_key / f"{drug}.ld_{r2_label}.window{WINDOW_KB}kb.details.tsv"
            write_lines(grey_file, grey_snps)
            write_tsv(detail_file, drug_details, detail_fields)
            print(f"  {drug}: {len(grey_snps)} grey SNPs", flush=True)

            for mode in MODES:
                render_trait(run_key, run_dir, drug, mode, grey_file)

            ld_manifest.append(
                {
                    "maf_label": run_key,
                    "drug": drug,
                    "n_grey_snps": str(len(grey_snps)),
                    "n_ld_detail_rows": str(len(drug_details)),
                    "grey_snps_file": str(grey_file),
                    "ld_details_file": str(detail_file),
                    "fungamr_source": str(FUNGAMR_MAPPING),
                    "ld_reference_plink": str(BASE_PLINK),
                }
            )

    figure_manifest: list[dict[str, str]] = []
    for run_key, maf_label, _run_dir in RUNS:
        for mode, mode_label in MODES.items():
            for plot_type in PLOT_TYPES:
                print(f"Composing {plot_type} {run_key} {mode}", flush=True)
                figure_manifest.append(compose_figure(run_key, maf_label, mode, mode_label, plot_type))

    write_tsv(
        OUTDIR / "fungamr_ld_grey_plot_manifest.tsv",
        ld_manifest,
        [
            "maf_label",
            "drug",
            "n_grey_snps",
            "n_ld_detail_rows",
            "grey_snps_file",
            "ld_details_file",
            "fungamr_source",
            "ld_reference_plink",
        ],
    )
    write_tsv(
        OUTDIR / "final_gwas_paper_figure_manifest.tsv",
        figure_manifest,
        ["figure", "plot_type", "maf", "mode", "output", "panel_A_sources", "panel_B_sources"],
    )
    print(f"Wrote edited final GWAS figures to {OUTDIR}", flush=True)


if __name__ == "__main__":
    main()
