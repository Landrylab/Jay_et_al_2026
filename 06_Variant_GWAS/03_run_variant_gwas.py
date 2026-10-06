#!/usr/bin/env python3
"""Fit additive quantitative-trait associations with an externally supplied GRM."""

import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess


def read_manifest(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"drug", "n", "tested_variants", "genotype_prefix", "phenotype_path", "qc_path"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Prepared manifest is empty or lacks required columns: " + ", ".join(sorted(required)))
    names = [row["drug"] for row in rows]
    if len(set(names)) != len(names):
        raise ValueError("Each drug must occur once in the prepared manifest")
    for name in names:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
            raise ValueError("Drug names must contain letters, digits, underscores or hyphens")
    return rows


def resolve_input(value, manifest):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (Path(manifest).resolve().parent / path).resolve()


def build_command(gcta, bfile, grm, phenotype, threads, out_prefix):
    return [str(gcta), "--mlma", "--bfile", str(bfile), "--grm", str(grm),
            "--pheno", str(phenotype), "--autosome-num", "9", "--maf", "0", "--thread-num", str(threads), "--out", str(out_prefix)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-manifest", type=Path, required=True)
    parser.add_argument("--grm", type=Path, required=True, help="GCTA GRM prefix for the full genotype cohort")
    parser.add_argument("--gcta", default="gcta64")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    executable = shutil.which(args.gcta)
    if executable is None:
        parser.error("GCTA executable not found; supply --gcta path/to/gcta64")
    rows = read_manifest(args.prepared_manifest)
    grm = args.grm.expanduser().resolve()
    for suffix in (".grm.bin", ".grm.N.bin", ".grm.id"):
        if not Path(str(grm) + suffix).is_file():
            parser.error("Missing GRM file: " + str(grm) + suffix)
    grm_ids = set()
    with Path(str(grm) + ".grm.id").open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) != 2 or tuple(fields) in grm_ids:
                raise ValueError("GRM IDs must be unique FID/IID pairs")
            grm_ids.add(tuple(fields))
    for row in rows:
        prefix = resolve_input(row["genotype_prefix"], args.prepared_manifest)
        phenotype = resolve_input(row["phenotype_path"], args.prepared_manifest)
        for suffix in (".bed", ".bim", ".fam"):
            if not Path(str(prefix) + suffix).is_file():
                raise FileNotFoundError(str(prefix) + suffix)
        pairs = [tuple(line.split()[:2]) for line in Path(str(prefix) + ".fam").read_text().splitlines()]
        phenotypes = [tuple(line.split()[:2]) for line in phenotype.read_text().splitlines()]
        if len(pairs) != int(row["n"]) or len(set(pairs)) != len(pairs) or pairs != phenotypes:
            raise ValueError("Sample count or phenotype order differs from the prepared genotype file")
        if not set(pairs).issubset(grm_ids):
            raise ValueError("Some association samples are absent from the supplied GRM")
        row["genotype_prefix"] = str(prefix)
        row["phenotype_path"] = str(phenotype)
        row["qc_path"] = str(resolve_input(row["qc_path"], args.prepared_manifest))
        if not Path(row["qc_path"]).is_file():
            raise FileNotFoundError(row["qc_path"])
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    completed = []
    for row in rows:
        dest = outdir / row["drug"]
        dest.mkdir(exist_ok=True)
        command = build_command(executable, row["genotype_prefix"], grm,
                                row["phenotype_path"], args.threads, dest / "gwas")
        (dest / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        with (dest / "console.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        result = dest / "gwas.mlma"
        logfile = dest / "gwas.log"
        if not result.is_file() or not logfile.is_file():
            raise RuntimeError("GCTA did not write its result and log files")
        if "Log-likelihood ratio converged." not in logfile.read_text():
            raise RuntimeError("GCTA did not report REML convergence; inspect " + str(logfile))
        completed.append({**row, "result_path": str(result), "grm": str(grm),
                          "model": "GCTA MLMA, additive dosage, intercept and supplied GRM"})
        print("Completed " + row["drug"], flush=True)
    with (outdir / "scan_manifest.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(completed[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(completed)


if __name__ == "__main__":
    main()
