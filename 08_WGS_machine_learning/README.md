# Genomic random forests

This workflow compares variant, gene-level variant and SNP-derived principal-component inputs for binary resistance prediction. All inputs derive from the same standardized SNP-and-indel master catalogue. Shared model code in `rf_common.py` keeps screening, tuning, threshold selection and metrics consistent across prediction and within-group importance analyses.

## Inputs

Install the Python packages in `requirements.txt`. Required resources are:

- The master prefix produced by the standardized variant workflow: `.dosage.npy` (variants × samples, int8), `.samples.tsv`, and `.variants.tsv.gz`. Dosages count copies of the specified ALT allele (0, 1, 2), with missing calls −1.
- The shared `selected_gene_variants.tsv.gz` annotation catalogue, containing canonical `index`, `id`, `gene_id` and `impact` columns and one selected HIGH/MODERATE-impact annotation per allele.
- The shared `gene_map.tsv`, with `gene_id`, `gene_symbol`, `chromosome`, `start`, `end` and `description` columns.
- The exact `common_nuclear_ld_pruned.ids` file used to build GWAS relatedness matrices. It is selected from the full cohort using nuclear SNPs, MAF ≥0.05, genotype missingness ≤0.10 and LD pruning with windows/steps/r² = 50/5/0.2.
- Metadata with unique `sample` identifiers, `rbin_flu` labels (`R` or `S`) and previously assigned genetic `cluster` labels.
- An optional explicit sample-order file, listing every binary-phenotyped genotype sample once.
- A cohort-metadata table for prediction and a group table for within-group analysis. `groups.example.tsv` illustrates group labels and outer-fold counts.

All paths are arguments or manifest entries. Original master genotype calls and variant representations remain traceable through `master_index` and the canonical allele identifier.

## Input definitions

**Variant input:** every eligible master SNP and indel, including mitochondrial variants, is initially retained. There is no GWAS MAF or minimum-carrier filter. At a decomposed multiallelic site, zero means absence of the specified ALT allele, not necessarily a reference-homozygous genotype.

**Gene input:** nuclear and mitochondrial HIGH/MODERATE-impact alleles are eligible when `0 < carrier_count / labelled_cohort_size <= 0.05`. The denominator includes all labelled isolates, including missing calls. Any observed eligible allele makes the gene indicator one; absence requires all constituent variants to be called and absent; otherwise it is missing. No 1% gene-carrier or noncarrier-frequency threshold is imposed. Eligibility uses the entire original labelled cohort before cross-validation or subgroup selection. Within-group analyses subset this same matrix without recomputing gene definitions.

**PC-marker input:** the exact shared GWAS relatedness marker list selects rows of the master matrix. The preparation step checks that all listed markers are nuclear SNPs and records a checksum of the list. It does not independently rerun pruning. The cache contains marker dosages; mean imputation, scaling and PCA are fitted within training partitions.

## Preparation and execution

```bash
python3 01_prepare_inputs.py \
  --master-prefix path/to/cohort.snps_indels.master \
  --selected-gene-annotations path/to/selected_gene_variants.tsv.gz \
  --gene-map path/to/gene_map.tsv \
  --pc-marker-ids path/to/common_nuclear_ld_pruned.ids \
  --metadata path/to/metadata.tsv \
  --sample-order path/to/phenotyped_sample_order.txt \
  --outdir path/to/rf_inputs

python3 06_run_analysis_suite.py --prepare \
  --inputs path/to/rf_inputs/model_inputs.tsv \
  --metadata path/to/metadata.tsv \
  --cohort-metadata path/to/subset_metadata.tsv \
  --groups groups.example.tsv \
  --outdir path/to/rf_analysis

python3 06_run_analysis_suite.py --outdir path/to/rf_analysis --jobs 2
python3 06_run_analysis_suite.py --outdir path/to/rf_analysis --summarize
```

The suite uses 5 repeats of five-fold stratified and clade-grouped validation for all three prediction inputs, and 5 repeats of the group-specific fold counts for both variant and gene inputs. It records outer and inner sample assignments before fitting. Each outer model is checkpointed in its own directory. Multiple worker processes may run the third command concurrently; atomic fold claims prevent duplication. A completed checkpoint is skipped on subsequent runs. An interrupted incomplete checkpoint must be inspected before restarting that fold. `--only-input` and `--only-scope` restrict workers to particular tasks.

For independent custom runs, script `02_run_random_forest.py` accepts the same input manifest, while `03_run_within_group_importance.py` accepts one cache and `--input-type variant` or `gene`. Script 01 can prepare selected `--input-types` while other shared inputs are still being generated.

## Common prediction procedure

The prediction cohort is supplied explicitly. Background selection can use a minimum of 15 phenotyped isolates and five resistant isolates when constructing its metadata. Gene-variant eligibility remains defined in the original binary-phenotyped cohort used for input preparation. Within-group sample order is inherited from that same original cohort.

- Outer stratified folds use `RepeatedStratifiedKFold`, seed 42. Grouped folds use `StratifiedGroupKFold`, seeds 42–46, keeping each genetic cluster entirely in training or testing. All three representations use identical samples and outer partitions.
- Variant/gene features must vary among called values and have no more than 10% missingness in the outer training partition. There is no positive lower frequency or minimum-carrier threshold. Training-modal imputation supplies scikit-learn's `chi2` score; up to 10,000 features are retained, ranked by decreasing score then training frequency.
- The frequency tie-breaker is `min(sum_called_values/(2*n_called), 1-sum_called_values/(2*n_called))`, including for binary gene indicators. It is a ranking criterion, not a gene-frequency filter.
- The selected feature set stays fixed during inner tuning. Modal imputation is fitted independently within each inner training split and refitted on the full outer training set. PC inputs undergo training-mean imputation, unit-variance scaling and randomized PCA separately in each training partition.
- Forests use Gini splits, 500 bootstrap-sampled trees and `balanced_subsample` class weights. Three-fold inner validation follows the outer stratified or grouped design. The nine ordered forest configurations are defined in `rf_common.py`. PC inputs evaluate six component counts (2, 5, 10, 20, 40, 80) crossed with three forest configurations. Highest mean inner ROC AUC wins; ties follow candidate order.
- The selected candidate's inner out-of-fold probabilities determine the threshold maximizing balanced accuracy. Candidate thresholds include the observed probabilities and 0, 0.5 and 1. Ties favour the threshold nearest 0.5, then the lower value. The final forest is refitted on all outer training samples and evaluated once on the held-out samples.

Inner splits and forests use seed `42 + outer_fold_number`. Inner PCA projections add the inner-fold number to this seed; outer PCA uses the outer-fold seed. Run records include the library versions and complete tuning results. Screening remains at the outer-training level during inner tuning, matching the stated analysis design.

## Outputs and importance

Each task exports held-out predictions, metrics, selected features, hyperparameters, thresholds, tuning scores and feature importance. The suite aggregates per-fold records and produces complete top-20 lists. Gini importance measures the sample-weighted decrease in class impurity attributable to a feature, normalized within trees and across the fitted forest. Mean importance divides the summed contribution by every outer model, assigning zero when a feature was not retained. Rankings have no effect direction and do not establish causality.

Script 07 recomputes prediction metrics from held-out probabilities and summarizes grouped splits. Script 04 provides prediction plots. Script 08 prepares reference-mapped gene importance and the current top-20 network from completed within-group results, using a reviewed evidence table with `gene_id`, `evidence_level` and `evidence_term` columns. Unreviewed loci retain their reference descriptions. Only current top-20 genes with supported evidence receive literature highlights. Script 05 draws the four importance panels and complete two-bank network at publication size; gene and description text use the matching marker colour. The first two configured groups occupy the left bank and the other two the right bank, with shared genes placed once in the first matching bank and all their links retained. Genes with more linked groups appear higher within each bank. Single-link genes follow configured group order, with importance rank ordering genes within each group. Evidence categories are `direct`, `response`, `pathway` and `none`, as defined in the annotation procedure.

```bash
python3 07_summarize_diagnostics.py --root path/to/rf_analysis
python3 08_prepare_figure_data.py \
  --rf-root path/to/rf_analysis --groups groups.example.tsv \
  --gene-map path/to/gene_map.tsv --reference-fai path/to/reference.fasta.fai \
  --evidence path/to/reviewed_gene_evidence.tsv --outdir path/to/within_figure
python3 05_plot_feature_importance.py \
  --importance-data path/to/within_figure/source_data/plot_data/ML_within_clade_gene_feature_importance_multipanel_plot_data.tsv \
  --network-data path/to/within_figure/source_data/plot_data/Within_cluster_gene_burden_top20_evidence_network_associated_genes.tsv \
  --outdir path/to/within_figure
```

Optional `--label-offsets` supplies a JSON mapping each group and labelled gene to a three-element array of x/y offsets in points and horizontal alignment. It changes label placement without modifying ranks, membership or values.

## Supplementary prediction figures

`09_plot_supplementary.py` renders held-out-background validation and within-background validation/variant-importance figures from supplied source tables. It does not refit models, select variants or rescale importance values. Inputs may be an Excel workbook containing the named tables below, a directory containing one `TABLE.tsv`, `TABLE.tsv.gz`, `TABLE.csv` or `TABLE.csv.gz` file per table, or individual files supplied with the corresponding argument. Excel input requires `openpyxl`. The shared `figure_sources.py` module remains alongside the numbered workflow directories.

```bash
python3 09_plot_supplementary.py grouped \
  --source path/to/grouped_source.xlsx \
  --output-prefix path/to/figures/grouped_validation

python3 09_plot_supplementary.py within \
  --source path/to/within_source_tables \
  --output-prefix path/to/figures/within_background
```

Table names are independent of worksheet names or positions. Background, feature-set and chromosome identifiers are matched exactly; text identifiers retain leading zeroes in delimited input. Table row order determines background and chromosome display order; the first occurrence of each feature set or metric determines its display order. Metric values use the interval [0, 1]. Blank, `NA` and `NaN` numeric values represent undefined metrics; undefined observations are omitted from summaries and displayed as `NA` when none are defined.

| Command | Named table | Required columns |
|---|---|---|
| `grouped` | `ED8_PanelA` | `Background`, `Isolates`, `Resistant isolates`, `Test fold` |
| `grouped` | `ED8_PanelsBE` | `Background`, `Feature set`, `Repeat`, `Balanced accuracy`, `Precision`, `Recall`, `ROC AUC` |
| `within` | `ED9_Cohorts` | `Background`, `Isolates`, `Resistant isolates`, `Outer folds` |
| `within` | `ED9_PanelsAD` | `Background`, `Feature set`, `Repeat`, `Fold`, `Metric`, `Value` |
| `within` | `ED9_PanelsEH` | `Background`, `Variant ID`, `Chromosome`, `Position (bp)`, `Mean Gini importance` |
| `within` | `ED9_Chromosomes` | `Chromosome`, `Length (bp)` |
| `within` | `ED9_Highlights` | `Background`, `Variant ID`, `Label` |

The grouped membership table contains one row per background and one positive test-fold identifier. This representation applies when the same background-to-fold assignment is used across repeats and feature sets. There must be at least two test folds, and every background/feature-set combination must contain the same repeat identifiers. Training membership is the complement of each test fold. Points show arithmetic means over defined repeats and horizontal whiskers span their minimum and maximum; these are not confidence intervals. Resistant and susceptible totals are derived from the supplied membership counts.

Within-background metrics contain one row per background, feature set, repeat, fold and metric. Every feature set and metric must cover the same fold keys within a background, and each repeat must contain its specified number of outer folds. Fold and repeat identifiers are positive integers. Boxes show the median and linear-interpolation quartiles; whiskers extend to observations within 1.5 interquartile ranges. Every defined fold value is overlaid, including outliers. Repeated folds can share isolates and are not independent biological replicates.

Variant positions and chromosome lengths are positive integers in base pairs, with one-based variant positions. Importance values must be finite and nonnegative. Each background/variant pair and each chromosome identifier must be unique. Chromosome order and lengths come from `ED9_Chromosomes`; `--chromosome-gap` sets the display gap in base pairs. The highlight table refers to existing background/variant pairs. Optional `Window start (bp)` and `Window end (bp)` columns define a display rectangle on that variant's chromosome, while `Label dx (pt)`, `Label dy (pt)` and `Label alignment` control label placement. Alignment is `left`, `right` or `center`. Use `--no-highlights` to omit highlights without supplying this table.

`--cohorts` and `--metrics` override their respective tables; the within command also accepts `--importance`, `--chromosomes` and `--highlights`. `--dpi` sets PNG resolution, and `--jitter-seed` controls within-panel point jitter. Each command writes PDF and PNG figures plus a JSON summary. The grouped command also writes `.group_metrics.tsv` and `.fold_composition.tsv`; the within command writes `.box_summary.tsv`. These exports record plotted summaries without altering the source files.

References: [Breiman (2001)](https://doi.org/10.1023/A:1010933404324), [scikit-learn](https://jmlr.org/papers/v12/pedregosa11a.html), and [impurity importance](https://proceedings.neurips.cc/paper/2013/hash/e3796ae838835da0b6f6ea37bcf8bcb7-Abstract.html).

## Included input data

`input_data.zip` contains prediction results, gene importance and network annotations, and the supplementary figure source workbooks. Run from this directory:

```bash
unzip input_data.zip
python3 04_plot_prediction.py \
  --results-root . --outdir plots/prediction
python3 05_plot_feature_importance.py \
  --importance-data gene_importance.tsv \
  --network-data gene_network.tsv \
  --label-offsets label_offsets.json \
  --outdir plots/gene_importance
python3 09_plot_supplementary.py grouped \
  --source Source_Data_Extended_Data_Fig_8.xlsx \
  --output-prefix plots/grouped_validation
python3 09_plot_supplementary.py within \
  --source Source_Data_Extended_Data_Fig_9.xlsx \
  --output-prefix plots/within_background
```

Use new output directories. The prediction command draws the numerical panels. Keep `../figure_sources.py` alongside the workflow directories for the supplementary commands.
