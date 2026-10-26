# Prioritizing resistance gene candidates

These workflows reproduce the gene-enrichment, joint co-carriage adjustment and individual-variant analyses underlying Extended Data Figure 10. They use prepared, fixed gene and variant carriage inputs, fit the statistical tests again and retain the full testing families. They do not write to the submission figure folders.

## Dependencies and execution

Use Python 3.10 or later, the packages in `requirements.txt`, and R with its base `stats` package. `Rscript` must be on `PATH`, or supply its executable path with `--rscript`. No additional R packages are needed. Keep `helpers/` beside the numbered scripts.

From this folder:

```sh
python3 -m pip install -r requirements.txt
unzip input_data.zip
python3 01_test_gene_enrichment.py --input-dir inputs --outdir results/gene_screen
python3 02_adjust_cocarriage.py \
  --gene-results results/gene_screen/gene_results.csv \
  --gene-calls results/gene_screen/sample_gene_calls.csv.gz \
  --outdir results/cocarriage
python3 03_test_variants.py \
  --input-dir inputs \
  --gene-results-dir results/gene_screen \
  --conditional-results results/cocarriage/conditional_results.csv \
  --outdir results/variants
python3 04_validate_results.py --input-dir inputs --results-dir results
```

Each output directory is supplied explicitly. The scripts retain their input files and accept paths outside this folder. Extracting the archive creates `inputs/`, including reference tables and a checksum manifest. Outputs from a different dataset need not match these study-specific references; the fourth entry point is a regression check for the supplied study inputs.

The co-carriage command accepts `--thresholds` followed by distinct finite values between 0 and 1, including the primary 0.20 threshold. Defaults are `0.20 0.15 0.25 0.30`. More precise sensitivity values keep distinct model identifiers and separate testing families.

## Workflows

1. **Gene enrichment.** The inventory contains the union of genes from the GWAS and within-background random forest networks. The GWAS union is assessed across six antifungals; genes contributed only by the random forest networks are assessed for the three azoles. Each gene–drug hypothesis is counted once. The script uses the prepared drug-specific classifications to select isolates without high-confidence FungAMR markers, following the marker definitions in Methods. It computes two-sided Fisher tests and exact conditional comparisons stratified by PCA cluster, with PCA-cluster-by-study sensitivity analyses. BH adjustment covers the complete gene–drug family separately for each method. The separate pooled-gene comparisons have their own testing family.
2. **Co-carriage adjustment.** The primary PCA-cluster screen selects associations with q < 0.05 and odds ratio > 1. For each target, partners are other qualifying genes for the same antifungal with an absolute within-cluster residual carriage correlation of at least 0.20. Sensitivity thresholds are 0.15, 0.25 and 0.30. Exact tests condition jointly on PCA cluster and the full selected-partner carriage pattern; PCA-only tests use the identical complete-call samples. BH adjustment retains the complete set of target tests separately for every threshold and model type. Genes with no selected partner or no informative comparison remain in the outputs.
3. **Individual variants.** The prepared variant calls must reconstruct each original qualifying gene burden exactly. All variants in the original screened gene–drug families are tested using exact PCA-cluster-stratified comparisons. BH adjustment is within each gene–drug family. The script then reports variants in associations that remain significant after primary joint co-carriage adjustment and have at least one selected partner. Individual-variant tests adjust for PCA clusters; they do not adjust for co-carried genes. The co-carriage adjustment determines which genes are followed up.
4. **Validation.** The checker verifies every packaged input checksum and compares all rows, columns, estimates, confidence intervals, counts and adjusted P values in the four principal result tables against the deposited references. It also checks the study's selection counts.

The supplied inputs reproduce 759 gene–drug tests, 36 qualifying gene–drug associations, eight retained associations and nine enriched variants within the retained genes. The full variant output retains 1,692 tests in 36 gene–drug families, including 34 enriched variants before the conditional gene selection. This full output preserves the original testing families.

## Inputs and their preparation

The archive contains prepared analysis inputs rather than sequencing reads or the full genotype master:

- `high_confidence_classifications.csv`: one row per eligible isolate–antifungal observation, including `sample`, `antifungal`, `phenotype`, `high_confidence_status`, `cluster_label` and `Study_clean`. The study table has 1,901 observations for 494 isolates. Marker classifications follow the manuscript Methods and the reusable `09_FungAMR` workflow.
- `network_gene_inventory.csv`: one row per candidate gene, with stable `gene_id`, display symbol, GWAS/RF membership and original network edges. The supplied inventory has 162 genes. Original network-selection outputs are produced by the association and within-background prediction workflows.
- `carriers/<drug>.npz`: `codes` is a gene-by-sample int8 matrix, accompanied by unique `gene_ids` and `sample_ids`. Codes are 1 for carrier, 0 for noncarrier and −1 for unknown. These are fixed, original drug-cohort rare HIGH/MODERATE burdens; variants are not reselected after restricting the test cohort. The original gene-burden workflow is in `07_Gene_level_variant_GWAS`. Presence takes precedence over missingness; without a present allele, any unknown qualifying allele makes the gene call unknown.
- `variant_metadata.csv`: one row per variant–gene–drug test, including unique `key`, `family`, `variant_id`, `variant_idx`, `gene_id`, `drug`, protein consequence and the original membership annotations.
- `sample_variant_calls.csv.gz`: one row per `key` and isolate, with `sample`, `phenotype`, `cluster_label` and `carrier`. These binary calls were derived from target-ALT dosage in the standardized master for the exact fixed gene-burden membership. The variant workflow verifies sample sets, phenotype/cluster joins, missingness codes and gene-burden reconstruction before testing.
- `reference_results/`: saved statistical reference tables used only by the validation command. No analysis stage reads its results from this directory.
- `manifest.json`: provenance, expected study counts and SHA-256 hashes for every packaged input and reference table.

For new data, supply equivalent prepared tables and matrices with unique identifiers and the same column contracts. Keep the original burden-selection cohort and the subsequent eligible test cohort distinct. Recalculate network inventories, marker classifications, fixed burdens and variant metadata from that dataset's upstream analyses; do not substitute the packaged study references as new observations.

## Outputs and interpretation

`gene_screen/gene_results.csv` and `sample_gene_calls.csv.gz` retain the complete gene screen and sample calls. `cocarriage/conditional_results.csv`, `model_specifications.csv` and `pairwise_pca_carriage_correlations.csv` record every target model and partner selection. The legacy field `q_36_tests` denotes the follow-up BH family; the supplied study has 36 tests, while the code calculates its family size from the supplied screen.

`variants/variant_results.csv` retains all tested variant families. `retained_gene_associations.csv` contains the conditional gene shortlist, `selected_gene_variant_results.csv` contains every variant tested in those genes, `prioritized_variants.csv` contains the nine qualifying variants, and `retained_gene_summary.csv` includes retained genes with no individually qualifying variant.

Raw `R` and `S` encodings and legacy machine fields are retained. Their display meanings are **R phenotype** and **S phenotype**. Unknown carriage is never recoded as absence. Untestable hypotheses enter BH adjustment as P = 1 and retain missing reported P/q values. Boundary odds ratios of zero or infinity are retained when exact inference is available. These are exploratory associations from the discovery cohort; neither conditional retention nor individual-variant enrichment establishes a causal mechanism or independent validation.
