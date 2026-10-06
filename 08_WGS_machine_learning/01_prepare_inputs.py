#!/usr/bin/env python3
"""Derive aligned prediction inputs from one canonical ALT-dosage catalogue."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def initialize_cache(path, samples, features):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    if (path / "shape.json").exists():
        raise ValueError("A completed cache already exists: " + str(path))
    (path / "sample_ids.txt").write_text("\n".join(samples) + "\n")
    features.to_csv(path / "variants.tsv", sep="\t", index=False)
    return np.memmap(path / "genotypes.int8.memmap", mode="w+", dtype=np.int8,
                     shape=(len(samples), len(features)))


def finish_cache(path, matrix, encoding):
    matrix.flush()
    (Path(path) / "genotype_encoding.txt").write_text(encoding + "\n")
    (Path(path) / "shape.json").write_text(json.dumps(
        {"n_samples": matrix.shape[0], "n_variants": matrix.shape[1]}, indent=2) + "\n")


def labelled_samples(metadata, all_samples, sample_order=None):
    metadata = pd.read_csv(metadata, sep="\t", dtype={"sample": str})
    if metadata['sample'].duplicated().any():
        raise ValueError("Metadata sample identifiers must be unique")
    eligible = metadata.loc[metadata.rbin_flu.isin(["R", "S"])].set_index("sample")
    samples = [s for s in all_samples if s in eligible.index]
    if sample_order is not None:
        supplied = Path(sample_order).read_text().splitlines()
        if len(supplied) != len(set(supplied)) or set(supplied) != set(samples):
            raise ValueError("Sample order must contain every labelled genotype sample exactly once")
        samples = supplied
    if set(eligible.loc[samples, "rbin_flu"]) != {"R", "S"}:
        raise ValueError("Both binary phenotype classes must be represented")
    return samples, eligible.loc[samples].reset_index()


def feature_table(variants):
    table = variants[["index", "id", "chrom", "pos", "ref", "alt", "type", "is_mito"]].copy()
    table = table.rename(columns={"index": "master_index"})
    table.insert(0, "variant_idx", np.arange(len(table)))
    table["feature"] = table['id']
    return table


def cache_dosages(dosage, variants, sample_indices, samples, path, marker_indices=None, chunk_size=4096):
    indices = np.arange(len(variants)) if marker_indices is None else np.asarray(marker_indices, dtype=int)
    features = feature_table(variants.iloc[indices].reset_index(drop=True))
    matrix = initialize_cache(path, samples, features)
    for start in range(0, len(indices), chunk_size):
        selected = indices[start:start + chunk_size]
        calls = np.asarray(dosage[np.ix_(selected, sample_indices)], dtype=np.int8)
        if not np.isin(calls, [-1, 0, 1, 2]).all():
            raise ValueError("Master dosages must be -1, 0, 1 or 2")
        matrix[:, start:start + len(selected)] = calls.T
    finish_cache(path, matrix, "target ALT copies 0/1/2; any missing genotype=-1")
    return features


def validate_pc_ids(path, variants):
    ids = Path(path).read_text().splitlines()
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Shared LD-pruned marker identifiers must be nonempty and unique")
    indexed = variants.set_index("id", drop=False)
    selected = indexed.loc[ids]
    is_snp = selected.ref.str.len().eq(1) & selected.alt.str.len().eq(1)
    if not is_snp.all() or selected.is_mito.astype(bool).any():
        raise ValueError("The shared GRM/PC marker panel must contain nuclear SNPs only")
    return selected['index'].to_numpy(int)


def collapse_gene_calls(calls):
    """Presence overrides missingness; complete absence is needed for zero."""
    return np.where((calls > 0).any(axis=0), 1,
                    np.where((calls < 0).any(axis=0), -1, 0)).astype(np.int8)


def cache_genes(dosage, variants, annotations, genes, sample_indices, samples, path, cap=0.05):
    required = {"index", "id", "gene_id", "impact"}
    if not required.issubset(annotations):
        raise ValueError("Shared selected annotations lack canonical index/id/gene_id/impact")
    selected = annotations.loc[annotations.impact.isin(["HIGH", "MODERATE"])].copy()
    if selected['index'].duplicated().any():
        raise ValueError("Use the shared catalogue with one selected HIGH/MODERATE-impact annotation per allele")
    selected = selected.merge(variants[["index", "id", "is_mito"]], on=["index", "id"], validate="one_to_one")
    gene_lookup = genes.set_index("gene_id")
    if not set(selected.gene_id).issubset(gene_lookup.index):
        raise ValueError("A selected annotation refers to an unmapped gene")
    indices = selected['index'].to_numpy(int)
    counts = np.zeros(len(selected), dtype=np.int64)
    missing = np.zeros(len(selected), dtype=np.int64)
    for start in range(0, len(indices), 4096):
        block = np.asarray(dosage[np.ix_(indices[start:start+4096], sample_indices)])
        counts[start:start+len(block)] = (block > 0).sum(axis=1)
        missing[start:start+len(block)] = (block < 0).sum(axis=1)
    selected["n_carriers"] = counts
    selected["n_missing"] = missing
    selected["carrier_frequency_total_cohort"] = counts / len(samples)
    selected["is_rare_for_collapse"] = (counts > 0) & (counts / len(samples) <= cap)
    included = selected.loc[selected.is_rare_for_collapse].copy()
    matrices, feature_rows, qc_rows = [], [], []
    for gene_id, members in included.groupby("gene_id", sort=True):
        row = gene_lookup.loc[gene_id]
        calls = np.asarray(dosage[np.ix_(members['index'].to_numpy(int), sample_indices)])
        codes = collapse_gene_calls(calls)
        matrices.append(codes)
        feature_rows.append({"gene_id": gene_id, "gene_symbol": "FKS1" if gene_id.removesuffix("_A") == "C1_02420C" else row.get("gene_symbol", ""),
                             "chrom": row["chrom"], "pos": int(row["pos"]),
                             "id": "GENE:" + str(gene_id), "feature": "GENE:" + str(gene_id),
                             "n_collapsed_variants": len(members)})
        n_called = int((codes >= 0).sum())
        n_carriers = int((codes == 1).sum())
        qc_rows.append({"gene_id": gene_id, "n_called": n_called, "n_carriers": n_carriers,
                        "n_noncarriers": int((codes == 0).sum()), "n_missing": int((codes < 0).sum()),
                        "carrier_frequency_called": n_carriers/n_called if n_called else float('nan'),
                        "n_collapsed_variants": len(members)})
    if not matrices:
        raise ValueError("No qualifying gene-carrier features")
    features = pd.DataFrame(feature_rows)
    features.insert(0, "variant_idx", np.arange(len(features)))
    matrix = initialize_cache(path, samples, features)
    matrix[:] = np.vstack(matrices).T
    finish_cache(path, matrix, "gene carrier=1; complete absence=0; missing=-1; presence overrides missing")
    pd.DataFrame(qc_rows).to_csv(Path(path)/"gene_marker_qc.tsv.gz", sep="\t", index=False)
    included.to_csv(Path(path)/"variant_membership.tsv.gz", sep="\t", index=False)
    info = {"n_genes": len(features), "qualifying_variant_memberships": len(included),
            "variant_carrier_frequency_cap": cap, "gene_frequency_filter": None,
            "frequency_scope": "entire binary-phenotyped cohort before cross-validation or subgroup selection",
            "variant_carrier_frequency_denominator": len(samples), "chromosome_scope": "nuclear and mitochondrial",
            "annotation_scope": "shared selected allele-specific HIGH/MODERATE annotations",
            "encoding": "carrier=1; complete absence=0; missing=-1; presence overrides missing"}
    (Path(path)/"gene_input_definition.json").write_text(json.dumps(info, indent=2)+"\n")
    return features, info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-prefix", type=Path, required=True)
    parser.add_argument("--selected-gene-annotations", type=Path, required=True)
    parser.add_argument("--gene-map", type=Path, required=True,
                        help="Shared gene map with gene_id, gene_symbol, chromosome, start and end")
    parser.add_argument("--pc-marker-ids", type=Path, required=True,
                        help="The exact nuclear LD-pruned marker list used for the GWAS GRM")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--sample-order", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--input-types", nargs="+", choices=["variant", "gene", "pc"], default=["variant", "gene", "pc"])
    args = parser.parse_args()
    prefix = str(args.master_prefix)
    dosage = np.load(prefix + ".dosage.npy", mmap_mode="r")
    variants = pd.read_csv(prefix + ".variants.tsv.gz", sep="\t", usecols=["index", "id", "chrom", "pos", "ref", "alt", "type", "is_mito"],
                           dtype={"chrom": str, "ref": str, "alt": str})
    all_samples = pd.read_csv(prefix + ".samples.tsv", sep="\t")['sample'].astype(str).tolist()
    if dosage.shape != (len(variants), len(all_samples)) or not np.array_equal(variants['index'], np.arange(len(variants))):
        raise ValueError("Master matrix dimensions and variant indices must agree")
    samples, metadata = labelled_samples(args.metadata, all_samples, args.sample_order)
    lookup = {sample:i for i,sample in enumerate(all_samples)}
    sample_indices = np.array([lookup[sample] for sample in samples])
    args.outdir.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.outdir/"labelled_metadata.tsv", sep="\t", index=False)
    provenance = {"master_prefix": prefix, "master_sample_order_sha256": sha256(prefix+".samples.tsv"),
                  "master_variant_catalogue_sha256": sha256(prefix+".variants.tsv.gz"),
                  "n_labelled_samples": len(samples), "R": int(metadata.rbin_flu.eq("R").sum()),
                  "S": int(metadata.rbin_flu.eq("S").sum()), "master_encoding": "specified ALT copies; missing=-1"}
    for kind in args.input_types:
        path = args.outdir/kind
        if kind == "gene":
            annotations = pd.read_csv(args.selected_gene_annotations, sep="\t", low_memory=False)
            genes = pd.read_csv(args.gene_map, sep="\t").fillna("")
            if "chromosome" in genes: genes["chrom"] = genes["chromosome"].astype(str)
            if "pos" not in genes: genes["pos"] = ((genes.start + genes.end) // 2).astype(int)
            features, info = cache_genes(dosage, variants, annotations, genes, sample_indices, samples, path)
            info.update(selected_annotation_catalogue=str(args.selected_gene_annotations.resolve()),
                        selected_annotation_sha256=sha256(args.selected_gene_annotations),
                        gene_map=str(args.gene_map.resolve()), gene_map_sha256=sha256(args.gene_map))
        else:
            indices = validate_pc_ids(args.pc_marker_ids, variants) if kind == "pc" else None
            features = cache_dosages(dosage, variants, sample_indices, samples, path, indices)
            info = {"chromosome_scope": "nuclear" if kind == "pc" else "nuclear and mitochondrial",
                    "encoding": "specified ALT copies 0/1/2; missing=-1"}
        info.update(provenance)
        if kind == "pc":
            info.update({"shared_pc_grm_marker_ids": str(args.pc_marker_ids),
                         "shared_pc_grm_marker_ids_sha256": sha256(args.pc_marker_ids),
                         "panel_filters": "nuclear SNPs; MAF>=0.05; missingness<=0.10; LD50/5/0.2 in full cohort"})
        (path/"input_provenance.json").write_text(json.dumps(info, indent=2)+"\n")
        print(f"{kind}: {len(samples)} samples, {len(features)} features", flush=True)
    rows=[]
    for kind in ["variant", "gene", "pc"]:
        if (args.outdir/kind/"shape.json").exists():
            shape=json.loads((args.outdir/kind/"shape.json").read_text())
            rows.append({"input_type":kind,"cache_dir":str((args.outdir/kind).resolve()),
                         "n_samples":shape["n_samples"],"n_features":shape["n_variants"]})
    pd.DataFrame(rows).to_csv(args.outdir/"model_inputs.tsv",sep="\t",index=False)


if __name__ == "__main__":
    main()
