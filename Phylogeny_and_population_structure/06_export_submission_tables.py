#!/usr/bin/env python3
# =============================================================================
# Description
#   Export tabular data underlying phylogeny, PCA, and clustering figures.
#
# Inputs
#   Newick tree, sample metadata, phenotype table, mutation-match table,
#   ADMIXTURE FAM and Q files, and PCA/clustering output directory.
#
# Outputs
#   Journal-ready TSV files with one row per plotted sample or statistic.
# =============================================================================

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
from Bio import Phylo


ANTIFUNGALS = [
    ("Fluconazole", "rbin_flu"),
    ("Voriconazole", "rbin_vor"),
    ("Posaconazole", "rbin_pos"),
    ("Caspofungin", "rbin_cas"),
    ("Micafungin", "rbin_mic"),
    ("Anidulafungin", "rbin_ani"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--phenotypes", required=True, type=Path)
    parser.add_argument("--mutation-matches", required=True, type=Path)
    parser.add_argument("--admixture-fam", required=True, type=Path)
    parser.add_argument("--admixture-q", required=True, type=Path)
    parser.add_argument("--pca-clustering-dir", required=True, type=Path)
    parser.add_argument("--pca-panel-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def require_unique(table: pd.DataFrame, column: str, label: str) -> None:
    duplicated = table[column].duplicated(keep=False)
    if duplicated.any():
        examples = ", ".join(table.loc[duplicated, column].astype(str).head(5))
        raise ValueError(f"Duplicate {label} identifiers: {examples}")


def load_tree_order(path: Path) -> list[str]:
    tree = Phylo.read(path, "newick")
    tips = [terminal.name for terminal in tree.get_terminals()]
    if any(tip is None for tip in tips):
        raise ValueError("Unnamed tree tips encountered.")
    if len(tips) != len(set(tips)):
        raise ValueError("Duplicate tree-tip identifiers encountered.")
    return tips


def load_admixture(fam_path: Path, q_path: Path) -> pd.DataFrame:
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str)
    q = pd.read_csv(q_path, sep=r"\s+", header=None)
    if len(fam) != len(q):
        raise ValueError("ADMIXTURE FAM and Q row counts differ.")
    if q.shape[1] != 28:
        raise ValueError(f"Expected 28 ADMIXTURE components; found {q.shape[1]}.")
    row_sums = q.sum(axis=1)
    if ((row_sums - 1.0).abs() > 1e-3).any():
        raise ValueError("One or more ADMIXTURE Q rows do not sum to one within tolerance.")
    q.columns = [f"ADMIXTURE_K28_component_{index:02d}" for index in range(1, q.shape[1] + 1)]
    q.insert(0, "sample", fam.iloc[:, 1].astype(str).values)
    require_unique(q, "sample", "ADMIXTURE sample")
    return q


def select_mutations(matches: pd.DataFrame) -> pd.DataFrame:
    required = {"sample", "gene", "mutation_found", "confidence_score"}
    missing = required.difference(matches.columns)
    if missing:
        raise ValueError(f"Mutation table lacks required columns: {sorted(missing)}")
    matches = matches[matches["confidence_score"].isin([1, 2])].copy()
    selected: list[pd.DataFrame] = []
    for gene, limit in (("ERG11", 4), ("FKS1", 2)):
        subset = matches[matches["gene"].eq(gene)].drop_duplicates(["sample", "mutation_found"])
        counts = (
            subset.groupby("mutation_found", as_index=False)["sample"]
            .nunique()
            .rename(columns={"sample": "carrier_count"})
            .sort_values(["carrier_count", "mutation_found"], ascending=[False, True])
            .head(limit)
        )
        counts.insert(0, "gene", gene)
        selected.append(counts)
    return pd.concat(selected, ignore_index=True)


def export_phylogeny_tables(args: argparse.Namespace, output_dir: Path) -> None:
    tips = load_tree_order(args.tree)
    order = pd.DataFrame({"tree_tip_order": range(1, len(tips) + 1), "sample": tips})

    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    phenotypes = pd.read_csv(args.phenotypes, sep="\t", dtype=str)
    matches = pd.read_csv(args.mutation_matches, dtype={"sample": str})
    admixture = load_admixture(args.admixture_fam, args.admixture_q)
    require_unique(metadata, "sample", "metadata sample")
    require_unique(phenotypes, "sample", "phenotype sample")
    if "Clade" in metadata.columns:
        metadata["Clade"] = metadata["Clade"].fillna("Other").replace({"": "Other", "NA": "Other", "nan": "Other"})

    metadata_columns = [
        column
        for column in ["sample", "sample_key", "Study_clean", "Loc", "MLST_clean", "Clade"]
        if column in metadata.columns
    ]
    phenotype_columns = ["sample"] + [column for _, column in ANTIFUNGALS if column in phenotypes.columns]
    table = order.merge(metadata[metadata_columns], on="sample", how="left", validate="one_to_one")
    table = table.merge(phenotypes[phenotype_columns], on="sample", how="left", validate="one_to_one")
    table = table.merge(admixture, on="sample", how="left", validate="one_to_one")

    selected = select_mutations(matches)
    mutation_columns: list[str] = []
    for row in selected.itertuples(index=False):
        column = f"{row.gene}_{row.mutation_found}"
        carriers = set(
            matches.loc[
                matches["gene"].eq(row.gene)
                & matches["mutation_found"].eq(row.mutation_found)
                & matches["confidence_score"].isin([1, 2]),
                "sample",
            ].astype(str)
        )
        table[column] = table["sample"].isin(carriers).astype(int)
        mutation_columns.append(column)

    required_complete = ["sample", "Clade"]
    missing_samples = set(tips).difference(metadata["sample"].astype(str))
    if missing_samples:
        raise ValueError(f"Tree samples absent from metadata: {sorted(missing_samples)[:5]}")
    missing_admixture = table.filter(like="ADMIXTURE_K28_component_").isna().all(axis=1)
    if missing_admixture.any():
        examples = ", ".join(table.loc[missing_admixture, "sample"].head(5))
        raise ValueError(f"Tree samples absent from ADMIXTURE Q data: {examples}")
    if table[required_complete].isna().any().any():
        raise ValueError("Missing required phylogeny annotation values encountered.")

    table.to_csv(output_dir / "phylogeny_figure_sample_metadata.tsv", sep="\t", index=False)
    selected.to_csv(output_dir / "phylogeny_selected_mutation_counts.tsv", sep="\t", index=False)
    order.to_csv(output_dir / "phylogeny_tree_tip_order.tsv", sep="\t", index=False)


def export_pca_tables(args: argparse.Namespace, output_dir: Path) -> None:
    pca_dir = args.pca_clustering_dir
    scores = pd.read_csv(pca_dir / "pca_scores_pc1_pc50.tsv", sep="\t", dtype={"sample": str})
    assignments = pd.read_csv(pca_dir / "final_cluster_assignments.tsv", sep="\t", dtype={"sample": str})
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    require_unique(scores, "sample", "PCA sample")
    require_unique(assignments, "sample", "cluster-assignment sample")
    require_unique(metadata, "sample", "metadata sample")
    pc_columns = [column for column in scores.columns if column.startswith("PC")]
    if pc_columns != [f"PC{index}" for index in range(1, 51)]:
        raise ValueError("PCA score table does not contain PC1-PC50 in order.")

    metadata_columns = [
        column
        for column in ["sample", "sample_key", "Study_clean", "Loc", "MLST_clean", "Clade"]
        if column in metadata.columns
    ]
    combined = scores.merge(assignments, on="sample", how="left", validate="one_to_one")
    combined = combined.merge(metadata[metadata_columns], on="sample", how="left", validate="one_to_one")
    if len(combined) != 1153:
        raise ValueError(f"Expected 1,153 PCA samples; found {len(combined):,}.")
    if combined["cluster"].isna().any():
        raise ValueError("Missing final cluster assignments encountered.")
    combined.to_csv(output_dir / "pca_scores_and_cluster_assignments.tsv", sep="\t", index=False)

    names = [
        "pca_variance_explained.tsv",
        "clustering_model_selection_statistics.tsv",
        "clustering_candidate_clusterwise_jaccard.tsv",
        "clustering_candidate_feature_stability.tsv",
        "final_clusterwise_jaccard.tsv",
        "cluster_centroids_scaled_pc1_pc20.tsv",
        "cluster_centroid_hierarchy_ward_d2.tsv",
        "clustering_selection_summary.tsv",
        "clustering_run_parameters.tsv",
    ]
    for name in names:
        shutil.copy2(pca_dir / name, output_dir / name)
    shutil.copy2(args.pca_panel_summary, output_dir / "variant_and_pca_panel_summary.tsv")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    export_phylogeny_tables(args, args.output_dir)
    export_pca_tables(args, args.output_dir)


if __name__ == "__main__":
    main()
