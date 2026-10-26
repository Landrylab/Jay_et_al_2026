#!/usr/bin/env python3
"""Compare rare-variant gene GWAS after excluding intermediate/SDD MIC ranges.

The fixed analysis retains original allele membership and subsets gene indicators.
The rebuilt analysis reapplies the 5% carrier-frequency criterion in each retained
cohort. Both use the original, full-cohort-standardized relationship matrices.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SUBMISSION = HERE.parents[1]
WORKSPACE = SUBMISSION.parents[1]
sys.path.insert(0, str(HERE.parent / "06_Variant_GWAS/helpers"))
from master_inputs import write_bed

RANGES = {"fluconazole": (2.0, 4.0), "voriconazole": (0.125, 0.5),
          "caspofungin": (0.25, 0.5), "micafungin": (0.25, 0.5),
          "anidulafungin": (0.25, 0.5)}


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


BUILD = module("sensitivity_gene_build", "02_build_gene_markers.py")
RUN = module("sensitivity_gene_run", "04_run_gene_gwas.py")


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 ** 2), b""):
            value.update(chunk)
    return value.hexdigest()


def json_write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def prepare(args):
    out = args.outdir
    if (out / "run_configuration.json").exists():
        raise FileExistsError("Sensitivity inputs already prepared: " + str(out))
    out.mkdir(parents=True, exist_ok=True)
    original_manifest_path = args.original_root / "inputs/prepared_scan_manifest.tsv"
    original = pd.read_csv(original_manifest_path, sep="\t").set_index("drug")
    metadata = pd.read_csv(args.mic_metadata, sep="\t", dtype={"sample": str})
    if metadata["sample"].duplicated().any():
        raise ValueError("MIC metadata must have unique sample identifiers")
    metadata = metadata.set_index("sample")
    phenotypes = out / "retained_phenotypes"
    phenotypes.mkdir(exist_ok=True)
    summaries, manifests, excluded_rows, retained_rows = [], [], [], []
    paths = [original_manifest_path, args.mic_metadata, args.loco_manifest,
             HERE / "02_build_gene_markers.py", HERE / "04_run_gene_gwas.py",
             HERE / "05_summarize_gene_results.py", Path(__file__).resolve(), args.gcta]
    for drug in args.drugs:
        row = original.loc[drug]
        source = RUN.resolve_path(row.genotype_prefix, original_manifest_path).parent
        samples = pd.read_csv(source / "samples.tsv", sep="\t")["sample"].astype(str).tolist()
        phenotype_path = RUN.resolve_path(row.phenotype_path, original_manifest_path)
        pheno = pd.read_csv(phenotype_path, sep=r"\s+", header=None,
                            names=["FID", "sample", "log2_mic"], dtype={"FID": str, "sample": str})
        if len(samples) != len(set(samples)) or samples != pheno["sample"].tolist() or not pheno.FID.equals(pheno["sample"]):
            raise ValueError("Original sample and phenotype identifiers/order disagree: " + drug)
        mic = pd.to_numeric(metadata.loc[samples, drug.capitalize()], errors="raise").to_numpy(float)
        if not np.isfinite(mic).all() or not (mic > 0).all():
            raise ValueError("Original cohort MICs must be positive and finite: " + drug)
        if not np.allclose(np.log2(mic), pheno.log2_mic.to_numpy(float), rtol=0, atol=1e-10):
            raise ValueError("Archived MICs disagree with GWAS log2 MIC phenotypes: " + drug)
        lower, upper = RANGES[drug]
        exclude = (mic > lower) & (mic <= upper)
        keep = ~exclude
        retained = [s for s, take in zip(samples, keep) if take]
        pheno.loc[keep].to_csv(phenotypes / f"pheno.mic.{drug}.txt", sep="\t", index=False, header=False)
        for i, sample in enumerate(samples):
            record = {"drug": drug, "sample": sample, "MIC": float(mic[i]), "log2_MIC": float(pheno.log2_mic.iloc[i])}
            (excluded_rows if exclude[i] else retained_rows).append(record)
        codes_path = source / "codes.npy"
        codes = np.load(codes_path)
        genes = pd.read_csv(source / "gene_metadata.tsv", sep="\t")
        if codes.shape != (len(genes), len(samples)) or not np.isin(codes, [-1, 0, 1]).all():
            raise ValueError("Original gene-code dimensions or values are invalid: " + drug)
        codes = codes[:, keep]
        genes["n_carriers"] = (codes == 1).sum(axis=1)
        genes["n_noncarriers"] = (codes == 0).sum(axis=1)
        genes["n_missing"] = (codes < 0).sum(axis=1)
        genes["n_called"] = genes.n_carriers + genes.n_noncarriers
        genes["carrier_frequency_called"] = genes.n_carriers / genes.n_called.replace(0, np.nan)
        polymorphic = genes.n_carriers.gt(0) & genes.n_noncarriers.gt(0)
        genes["is_polymorphic_gene_marker"] = polymorphic
        genes["passes_gene_marker_maf_filter"] = polymorphic
        tested, selected_codes = genes.loc[polymorphic].copy(), codes[polymorphic]
        dest = out / "fixed/inputs" / drug
        dest.mkdir(parents=True, exist_ok=True)
        genes.to_csv(dest / "gene_qc.tsv", sep="\t", index=False)
        tested.to_csv(dest / "gene_metadata.tsv", sep="\t", index=False)
        pd.DataFrame({"sample": retained}).to_csv(dest / "samples.tsv", sep="\t", index=False)
        np.save(dest / "codes.npy", selected_codes)
        pheno.loc[keep].to_csv(dest / "phenotype.log2_mic.txt", sep="\t", index=False, header=False)
        tested = tested.assign(cm=0, alt="2", ref="1")
        dosage = np.where(selected_codes < 0, -1, selected_codes * 2).astype(np.int8)
        write_bed(dest / "genotypes", tested, retained, [dosage])
        manifests.append({"drug": drug, "n": len(retained), "tested_gene_markers": len(tested),
                          "mitochondrial_gene_markers": int(tested.chromosome.eq(8).sum()),
                          "qualifying_variants": int(row.qualifying_variants),
                          "genotype_prefix": str((dest / "genotypes").resolve()),
                          "phenotype_path": str((dest / "phenotype.log2_mic.txt").resolve()),
                          "qc_path": str((dest / "gene_qc.tsv").resolve()),
                          "original_variant_membership": str((source / "variant_membership.tsv.gz").resolve()),
                          "original_frequency_cohort_n": len(samples), "maximum_carrier_frequency": 0.05})
        summaries.append({"drug": drug, "original_n": len(samples), "excluded_n": int(exclude.sum()),
                          "retained_n": len(retained), "exclusion_lower_exclusive": lower,
                          "exclusion_upper_inclusive": upper, "original_tested_genes": len(genes),
                          "fixed_tested_genes": len(tested), "fixed_constant_genes_removed": int((~polymorphic).sum())})
        paths.extend([codes_path, source / "samples.tsv", source / "gene_metadata.tsv",
                      source / "variant_membership.tsv.gz", phenotype_path])
    pd.DataFrame(manifests).to_csv(out / "fixed/inputs/prepared_scan_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(out / "cohort_summary.tsv", sep="\t", index=False)
    pd.DataFrame(excluded_rows).to_csv(out / "excluded_samples.tsv", sep="\t", index=False)
    pd.DataFrame(retained_rows).to_csv(out / "retained_samples.tsv", sep="\t", index=False)
    BUILD.build(args.master_prefix, args.catalog, args.gene_map, phenotypes,
                out / "rebuilt/inputs", drugs=args.drugs, cap=0.05)
    rebuilt = pd.read_csv(out / "rebuilt/inputs/prepared_scan_manifest.tsv", sep="\t").set_index("drug")
    for row in summaries:
        drug = row["drug"]
        a = pd.read_csv(out / "fixed/inputs" / drug / "samples.tsv", sep="\t")
        b = pd.read_csv(out / "rebuilt/inputs" / drug / "samples.tsv", sep="\t")
        if not a.equals(b):
            raise ValueError("Fixed and rebuilt sample order differs: " + drug)
        row["rebuilt_tested_genes"] = int(rebuilt.loc[drug, "tested_gene_markers"])
        row["rebuilt_qualifying_variants"] = int(rebuilt.loc[drug, "qualifying_variants"])
    pd.DataFrame(summaries).to_csv(out / "cohort_summary.tsv", sep="\t", index=False)
    paths.extend([args.catalog, args.gene_map,
                  *[Path(str(args.master_prefix) + suffix)
                    for suffix in [".dosage.npy", ".variants.tsv.gz", ".samples.tsv"]]])
    grms = pd.read_csv(args.loco_manifest, sep="\t")
    for value in grms.grm_prefix:
        prefix = RUN.resolve_path(value, args.loco_manifest)
        for suffix in [".grm.id", ".grm.bin", ".grm.N.bin"]:
            paths.append(Path(str(prefix) + suffix))
    config = dict(drugs=args.drugs, original_root=str(args.original_root.resolve()),
                  gcta=str(args.gcta.resolve()), threads=args.threads, loco_manifest=str(args.loco_manifest.resolve()),
                  exclusion_ranges=RANGES, master_prefix=str(args.master_prefix.resolve()),
                  fixed_definitions="Original allele memberships and gene indicators; retained sample columns",
                  rebuilt_definitions="Original selected functional annotations; carrier eligibility recomputed at 5% of retained cohort",
                  relatedness="Original full-cohort-standardized GRMs subset internally by GCTA",
                  created_utc=datetime.now(timezone.utc).isoformat(),
                  sources=[dict(path=str(p.resolve()), sha256=sha256(p)) for p in dict.fromkeys(paths)])
    json_write(out / "run_configuration.json", config)
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def scan_and_summarize(manifest, grms, destination, gcta, threads):
    RUN.run(manifest, grms, destination / "scans", gcta=gcta, threads=threads)
    subprocess.run([sys.executable, str(HERE / "05_summarize_gene_results.py"),
                    "--manifest", str(destination / "scans/association_manifest.tsv"),
                    "--outdir", str(destination / "results")], check=True)


def control(args, config):
    dest = args.outdir / "control"
    report_path = dest / "comparison.json"
    if report_path.exists():
        report = json.loads(report_path.read_text())
        if report["status"] != "PASS":
            raise ValueError("The original-data reproduction control failed")
        return report
    (dest / "inputs").mkdir(parents=True, exist_ok=True)
    source = Path(config["original_root"])
    original_manifest = source / "inputs/prepared_scan_manifest.tsv"
    original = pd.read_csv(original_manifest, sep="\t")
    for field in ["genotype_prefix", "phenotype_path", "qc_path"]:
        original[field] = original[field].map(lambda value: str(RUN.resolve_path(value, original_manifest)))
    manifest = dest / "inputs/prepared_scan_manifest.tsv"
    original.loc[original.drug.eq("fluconazole")].to_csv(manifest, sep="\t", index=False)
    started = time.monotonic()
    scan_and_summarize(manifest, Path(config["loco_manifest"]), dest, config["gcta"], config["threads"])
    before = pd.read_csv(source / "results/fluconazole.tsv.gz", sep="\t").set_index("gene_id")
    after = pd.read_csv(dest / "results/fluconazole.tsv.gz", sep="\t").set_index("gene_id")
    if set(before.index) != set(after.index):
        raise ValueError("Reproduction control gene sets differ")
    after = after.loc[before.index]
    differences = {field: float(np.abs(before[field] - after[field]).max()) for field in ["b", "se", "p", "Freq"]}
    exact = all(value == 0 for value in differences.values())
    close = all(np.allclose(before[f], after[f], rtol=1e-10, atol=1e-15) for f in differences)
    report = dict(status="PASS" if close else "FAIL", exact=exact, tested_genes=len(before),
                  maximum_absolute_differences=differences, elapsed_seconds=time.monotonic() - started,
                  gcta=str(config["gcta"]), gcta_sha256=sha256(config["gcta"]))
    json_write(report_path, report)
    if not close:
        raise ValueError("Original-data reproduction differs beyond tolerance; inspect " + str(report_path))
    print("Original fluconazole reproduction control: " + json.dumps(report), flush=True)
    return report


def run(args):
    config = json.loads((args.outdir / "run_configuration.json").read_text())
    control(args, config)
    for mode in ["fixed", "rebuilt"]:
        dest = args.outdir / mode
        if (dest / "completion.json").exists():
            continue
        started = time.monotonic()
        scan_and_summarize(dest / "inputs/prepared_scan_manifest.tsv", Path(config["loco_manifest"]),
                           dest, config["gcta"], config["threads"])
        json_write(dest / "completion.json", dict(status="PASS", elapsed_seconds=time.monotonic() - started,
                                                  drugs=config["drugs"]))
        print("Completed " + mode, flush=True)
    config.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat())
    json_write(args.outdir / "run_configuration.json", config)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "run", "all"])
    parser.add_argument("--drugs", nargs="+", choices=list(RANGES), default=["fluconazole"])
    parser.add_argument("--original-root", type=Path, default=SUBMISSION / "Analysis/Gene_GWAS")
    parser.add_argument("--outdir", type=Path, default=SUBMISSION / "Analysis/Gene_GWAS_category_sensitivity_20261007")
    parser.add_argument("--mic-metadata", type=Path, default=SUBMISSION / "Main_figures_archive/02_fungamr/source_data/panel_A_and_B/all_metadata.tsv")
    parser.add_argument("--master-prefix", type=Path, default=SUBMISSION / "Variant_data/cohort1153.snps_indels.master")
    parser.add_argument("--catalog", type=Path, default=SUBMISSION / "Variant_data/selected_gene_variants.tsv.gz")
    parser.add_argument("--gene-map", type=Path, default=SUBMISSION / "Variant_data/gene_map.tsv")
    parser.add_argument("--loco-manifest", type=Path, default=SUBMISSION / "Analysis/GRM/loco_manifest.tsv")
    parser.add_argument("--gcta", type=Path, default=WORKSPACE / "f/GWAS_attempts/All_six_antifungals_MLMA_100kb_sliding_window_no_MAF_review/tools/gcta/gcta-1.95.0-macOS-arm64/bin/gcta64")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if len(args.drugs) != len(set(args.drugs)):
        parser.error("Each drug may be selected once")
    if args.action in ["prepare", "all"]:
        prepare(args)
    if args.action in ["run", "all"]:
        run(args)


if __name__ == "__main__":
    main()
