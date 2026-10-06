#!/usr/bin/env python3
"""Test gene-level carrier markers with chromosome-excluded common-SNP GRMs.

GCTA fits an intercept-only null mixed model separately for each drug and tested
chromosome. Carrier markers use A1='2' and A2='1'; raw GCTA estimates refer to one
encoded A1 copy, with carriers represented by two copies. Output P values and
effect estimates are retained as reported by GCTA for subsequent summarization.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pandas as pd


def resolve_path(value, manifest):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (Path(manifest).resolve().parent / path).resolve()


def parse_variance_components(log_path):
    text = Path(log_path).read_text(errors="replace")
    if "Log-likelihood ratio converged." not in text or "Analysis finished" not in text:
        raise ValueError(f"GCTA did not complete with REML convergence: {log_path}")
    rows = [line.split() for line in text.splitlines()
            if re.match(r"^\s*\d+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*$", line)]
    if not rows:
        raise ValueError(f"Variance components missing from GCTA log: {log_path}")
    vg, ve = float(rows[-1][2]), float(rows[-1][3])
    if not np.isfinite([vg, ve]).all() or vg < 0 or ve <= 0:
        raise ValueError("GCTA returned invalid genetic/residual variances")
    return vg, ve


def read_burden_inputs(row, manifest):
    bfile = resolve_path(row["genotype_prefix"], manifest)
    phenotype = resolve_path(row["phenotype_path"], manifest)
    qc_path = resolve_path(row["qc_path"], manifest)
    bim = pd.read_csv(str(bfile) + ".bim", sep=r"\s+", header=None,
                      names=["Chr", "SNP", "cm", "bp", "A1", "A2"], dtype=str)
    fam = pd.read_csv(str(bfile) + ".fam", sep=r"\s+", header=None, dtype=str)
    pheno = pd.read_csv(phenotype, sep=r"\s+", header=None, dtype={0: str, 1: str})
    if fam.shape[1] != 6 or pheno.shape[1] != 3 or fam[[0, 1]].duplicated().any():
        raise ValueError("Burden FAM and phenotype files have invalid format")
    if len(fam) != int(row["n"]) or not fam[[0, 1]].equals(pheno[[0, 1]]):
        raise ValueError("Burden genotypes and phenotypes must contain the same sample order and count")
    if not np.isfinite(pd.to_numeric(pheno[2], errors="raise")).all() or pheno[2].eq(-9).any():
        raise ValueError("Burden phenotype files must contain finite observed log2 MIC values")
    if len(bim) != int(row["tested_gene_markers"]) or bim.empty or bim.SNP.duplicated().any():
        raise ValueError("BIM gene-marker count or distinctness differs from the manifest")
    if not bim.A1.eq("2").all() or not bim.A2.eq("1").all() or not bim.Chr.str.fullmatch(r"[1-9][0-9]*").all():
        raise ValueError("Gene markers require numeric chromosomes, carrier A1='2' and noncarrier A2='1'")
    bed = Path(str(bfile) + ".bed")
    if bed.stat().st_size != 3 + len(bim) * ((len(fam) + 3) // 4):
        raise ValueError("Burden BED dimensions differ from the BIM/FAM files")
    with bed.open("rb") as handle:
        if handle.read(3) != bytes([108, 27, 1]):
            raise ValueError("Burden BED must be SNP-major")
    if not qc_path.exists():
        raise ValueError("Gene marker QC table is missing")
    return bfile, phenotype, qc_path, bim, fam


def check_result(result, expected):
    if result.SNP.duplicated().any() or set(result.SNP) != set(expected.SNP):
        raise ValueError("GCTA results do not match the chromosome's eligible gene markers")
    aligned = expected.set_index("SNP").loc[result.SNP]
    for name in ("Chr", "bp", "A1", "A2"):
        if result[name].astype(str).tolist() != aligned[name].astype(str).tolist():
            raise ValueError("GCTA changed gene-marker coordinates or effect alleles: " + name)
    for column in ("b", "se", "p", "Freq"):
        if not np.isfinite(result[column]).all():
            raise ValueError("GCTA returned nonfinite association statistics")
    if not (result.se.gt(0) & result.p.between(0, 1) & result.Freq.gt(0) & result.Freq.lt(1)).all():
        raise ValueError("GCTA returned invalid or nonpolymorphic gene-marker statistics")


def run(manifest, loco_manifest, outdir, gcta="gcta64", threads=2):
    executable = shutil.which(str(gcta))
    if executable is None:
        raise ValueError("GCTA executable was not found: " + str(gcta))
    if threads < 1:
        raise ValueError("Thread count must be positive")
    prepared = pd.read_csv(manifest, sep="\t")
    required = {"drug", "n", "tested_gene_markers", "genotype_prefix", "phenotype_path", "qc_path"}
    if not required.issubset(prepared) or prepared.empty or prepared.drug.duplicated().any():
        raise ValueError("Prepared manifest must contain one distinct row per drug and required input paths/counts")
    if not prepared.drug.str.fullmatch(r"[A-Za-z0-9_-]+").all():
        raise ValueError("Drug names must contain letters, digits, underscores or hyphens")
    grms = pd.read_csv(loco_manifest, sep="\t", dtype={"excluded_chr": str})
    if not {"excluded_chr", "grm_prefix"}.issubset(grms) or grms.excluded_chr.duplicated().any():
        raise ValueError("LOCO manifest requires distinct excluded_chr and grm_prefix columns")
    grm_paths = {row.excluded_chr: resolve_path(row.grm_prefix, loco_manifest)
                 for row in grms.itertuples()}
    inputs = [(row, read_burden_inputs(row, manifest)) for row in prepared.to_dict("records")]
    for row, (_, _, _, bim, fam) in inputs:
        for chromosome in bim.Chr.unique():
            if chromosome not in grm_paths:
                raise ValueError("Missing chromosome-excluded GRM for chromosome " + chromosome)
            prefix = grm_paths[chromosome]
            ids = pd.read_csv(str(prefix) + ".grm.id", sep=r"\s+", header=None, dtype=str)
            pairs = set(map(tuple, ids.to_numpy()))
            if ids.shape[1] != 2 or ids.duplicated().any() or not set(map(tuple, fam[[0, 1]].to_numpy())).issubset(pairs):
                raise ValueError("GRM IDs are not distinct or do not contain the complete phenotype cohort")
            expected_size = len(ids) * (len(ids) + 1) // 2 * 4
            for suffix in (".grm.bin", ".grm.N.bin"):
                if Path(str(prefix) + suffix).stat().st_size != expected_size:
                    raise ValueError("GRM binary dimensions disagree with the sample IDs")
    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    scans, fits = [], []
    for row, (bfile, phenotype, qc_path, bim, fam) in inputs:
        destination = outdir / row["drug"]
        destination.mkdir(exist_ok=True)
        pieces = []
        for chromosome in sorted(bim.Chr.unique(), key=int):
            selected = bim.loc[bim.Chr.eq(chromosome)]
            work = destination / f"chr{chromosome}"
            work.mkdir(exist_ok=True)
            extract = work / "gene_markers.extract.txt"
            extract.write_text("\n".join(selected.SNP) + "\n")
            prefix = work / "gwas"
            command = [executable, "--bfile", str(bfile), "--extract", str(extract),
                       "--grm", str(grm_paths[chromosome]), "--mlma", "--pheno", str(phenotype),
                       "--autosome-num", str(max(map(int, grm_paths))), "--maf", "0",
                       "--thread-num", str(threads), "--out", str(prefix)]
            (work / "gcta_command.json").write_text(json.dumps(command, indent=2) + "\n")
            with (work / "gcta_console.log").open("w") as handle:
                subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
            vg, ve = parse_variance_components(str(prefix) + ".log")
            result = pd.read_csv(str(prefix) + ".mlma", sep=r"\s+", dtype={"Chr": str, "SNP": str, "A1": str, "A2": str})
            check_result(result, selected)
            result["loco_excluded_chromosome"] = chromosome
            result["relatedness_adjustment"] = "full nuclear" if chromosome == "8" else "LOCO"
            pieces.append(result)
            fits.append({"drug": row["drug"], "excluded_chromosome": chromosome,
                         "tested_genes": len(result), "Vg": vg, "Ve": ve,
                         "result_path": str(prefix) + ".mlma", "grm_prefix": str(grm_paths[chromosome])})
        merged = pd.concat(pieces, ignore_index=True)
        result_path = destination / "gene_burden_loco.mlma"
        merged.to_csv(result_path, sep="\t", index=False)
        scans.append({"drug": row["drug"], "n": len(fam), "tested_gene_markers": len(merged),
                      "result_path": str(result_path), "qc_path": str(qc_path)})
    pd.DataFrame(scans).to_csv(outdir / "association_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(fits).to_csv(outdir / "chromosome_fit_manifest.tsv", sep="\t", index=False)
    return scans


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--loco-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--gcta", default="gcta64")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    try:
        result = run(args.manifest, args.loco_manifest, args.outdir, args.gcta, args.threads)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
