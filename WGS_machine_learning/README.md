# Whole-genome fluconazole-resistance machine learning

This directory reproduces only the machine-learning analyses reported in Figures 3 and 4.

## Reported analyses

The all-sample workflow compares three feature sets in 871 labelled isolates (264 resistant and 607 susceptible):

- **SNP only:** each feature records whether an isolate is homozygous for the alternate allele.
- **PC only:** features are principal components calculated inside each training split from an LD-pruned genotype panel.
- **SNP+PC:** the fold-selected SNP indicators and training-derived PCs are combined.

Each feature set is evaluated with stratified and genetic-cluster-grouped cross-validation. The grouped analysis keeps every member of a K=23 PCA cluster in a single fold. Logistic regression, linear support vector machine, random forest and XGBoost are compared in each analysis.

The paired workflow forms unique 1:1 resistant-susceptible matches within the same K=23 genetic cluster by minimizing phylogenetic patristic distance. A resistant-minus-susceptible SNP vector is labelled 1 and the reversed vector is labelled 0. Pair-stratified validation keeps both directions from a pair together; genetic-cluster-grouped validation also keeps all pairs from a cluster together. Only feature importance from the genetic-cluster-grouped paired random forest is used for marker prioritization.

## Software

Python 3.10 or later is recommended. Required Python packages are listed in `requirements.txt`. BCFtools is required only for cache construction.

## Input files

- A metadata TSV containing `sample`, binary fluconazole label `rbin_flu` (`R` or `S`) and K=23 cluster `cluster`.
- An indexed, decomposed biallelic SNP VCF for the homozygous-alternate feature matrix.
- An indexed, LD-pruned biallelic SNP VCF for fold-derived PCA.
- A square NumPy patristic-distance matrix and a one-column TSV giving its sample order.
- The MAF ≥0.05 binary fluconazole GWAS table.
- A SnpEff-annotated cohort VCF.
- A candidate-window TSV with `region_id`, `chrom`, `start` and `end`.
- A TSV listing samples carrying a high-confidence FungAMR fluconazole mutation.

The two VCFs must use the final 1,153-isolate cohort. Only the 871 isolates with binary fluconazole labels enter model fitting.

## Reproduction sequence

```bash
python3 01_build_genotype_caches.py \
  --metadata path/to/metadata.tsv \
  --snp-vcf path/to/final_cohort_biallelic_snps.vcf.gz \
  --pc-vcf path/to/ld_pruned_pc_markers.vcf.gz \
  --output path/to/ml_cache

python3 02_run_all_sample_models.py \
  --cache path/to/ml_cache \
  --metadata path/to/metadata.tsv \
  --output path/to/all_sample_models

python3 03_run_paired_random_forest.py \
  --cache path/to/ml_cache \
  --metadata path/to/metadata.tsv \
  --tree-distance-matrix path/to/patristic_distances.npy \
  --tree-distance-samples path/to/patristic_distance_samples.tsv \
  --output path/to/paired_random_forest

python3 04_run_permutation_controls.py \
  --cache path/to/ml_cache \
  --metadata path/to/metadata.tsv \
  --all-sample-predictions path/to/all_sample_models/fold_predictions.tsv \
  --grouped-selected-snps path/to/all_sample_models/selected_snps_by_fold.tsv \
  --output path/to/permutation_controls

python3 05_build_marker_support_tables.py \
  --paired-feature-summary path/to/paired_random_forest/feature_importance_summary.tsv \
  --binary-gwas path/to/fluconazole_maf0.05_binary_gwas.tsv \
  --candidate-windows path/to/candidate_windows.tsv \
  --annotated-vcf path/to/snpeff_annotated_cohort.vcf.gz \
  --metadata path/to/metadata.tsv \
  --fungamr-carriers path/to/high_confidence_fungamr_fluconazole_carriers.tsv \
  --output path/to/marker_support
```

The full analysis is computationally intensive. The `--analyses` and `--models` arguments in script 02 allow independent jobs to be run and combined after completion without changing the reported analysis definitions.

## Fold safety and output interpretation

SNP filtering, chi-square ranking, imputation, scaling and PCA are estimated from training data only. Hyperparameters and PC count are chosen by three-fold inner cross-validation, and final performance is reported over five repeats of five-fold outer validation.

The score-label permutation control keeps held-out prediction scores fixed and shuffles labels; it is not a model-retraining null. The feature-screen permutation control repeats the grouped SNP screen after label shuffling but does not refit the classifiers.

Random-forest importance is the mean decrease in node impurity. For the reported group-stratified paired ranking, a SNP not selected in an outer fold is assigned zero importance for that fold before averaging across all 25 folds. Importance is a predictive ranking and does not establish causality.
