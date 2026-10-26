# Gene-level variant GWAS

This analysis uses the same annotated SNP/indel master as variant GWAS and RF. Annotation is performed once in `02_Variant_calling_and_filtering`; do not independently annotate an alternative genotype panel. Required inputs are the master dosage/metadata/sample files, `selected_gene_variants.tsv.gz`, `gene_map.tsv`, the variant-GWAS phenotype directory and the shared GRM directory. Scripts use GCTA 1.95.0 and Python with NumPy, pandas and SciPy. Keep the numbered workflow directories together.

## Gene indicators

The shared annotation catalogue assigns each allele one mapped HIGH/MODERATE consequence, prioritizing HIGH, then a protein-level consequence, and finally deterministic gene/effect/HGVS ordering. Both nuclear and mitochondrial genes are eligible.

Within each drug cohort, retain alleles carried by at least one isolate and no more than 5% of the total cohort. This denominator includes isolates with missing genotype calls. No lower MAF, five-carrier minimum or variant-missingness threshold is imposed. A gene indicator is present when any retained allele is carried; absent when all retained alleles are called and absent; and missing when none is observed but one or more calls are missing. Presence takes precedence over missingness. A gene is tested when at least one called carrier and noncarrier are present. There is no 1% gene-frequency filter.

`codes.npy` stores genes × samples as 0 (noncarrier), 1 (carrier) or -1 (missing); accompanying metadata and sample tables retain ordering. BED files encode carriers as two A1 copies and noncarriers as zero, with A1='2' and A2='1'. To express the carrier–noncarrier effect, multiply GCTA's raw coefficient and SE by two. This leaves P values unchanged.

```bash
python3 02_build_gene_markers.py \
  --master-prefix path/to/cohort.snps_indels.master \
  --catalog path/to/selected_gene_variants.tsv.gz \
  --gene-map path/to/gene_map.tsv \
  --phenotype-dir path/to/phenotypes --outdir path/to/gene_inputs

python3 04_run_gene_gwas.py \
  --manifest path/to/gene_inputs/prepared_scan_manifest.tsv \
  --loco-manifest path/to/GRM/loco_manifest.tsv \
  --gcta path/to/gcta64 --threads 2 --outdir path/to/gene_associations

python3 05_summarize_gene_results.py \
  --manifest path/to/gene_associations/association_manifest.tsv \
  --outdir path/to/gene_results
```

## Relatedness and testing

The full and LOCO matrices are produced once by `06_Variant_GWAS/02_build_relationship_matrix.py` from the shared LD-pruned nuclear SNP file. Nuclear genes use the matrix excluding their own chromosome. Mitochondrial genes use the full nuclear matrix. GCTA fits null REML variance components separately per drug/chromosome. Intercept, quantitative phenotype handling, mean imputation and variant-effect testing follow variant GWAS.

Bonferroni and suggestive thresholds use 0.05/M and 1/M, with M the number of tested genes across all chromosomes for each drug. Lambda GC is `median((b/se)^2)/0.45493642311957275`.

```bash
python3 06_plot_gene_gwas.py \
  --results path/to/gene_results/all_gene_level_variant_results.tsv.gz \
  --gene-evidence path/to/reviewed_gene_evidence.tsv \
  --outdir path/to/gene_plots
```

The optional evidence table contains `drug`, `gene_id`, `display_symbol` and `evidence_level`. Categories `direct` and `response` highlight previously supported genes among the current top 20 for each drug, independently of statistical significance. Every top-20 membership is retained in exported source data.

References: [SnpEff](https://doi.org/10.4161/fly.19695), [GCTA](https://doi.org/10.1016/j.ajhg.2010.11.011), [LOCO mixed models](https://doi.org/10.1038/ng.2876), [genomic control](https://doi.org/10.1111/j.0006-341X.1999.00997.x).

## Intermediate/SDD MIC-range sensitivity analysis

`07_run_category_sensitivity.py` repeats the quantitative log₂ MIC gene GWAS after excluding measurements in the selected CLSI intermediate or susceptible-dose-dependent ranges: fluconazole >2 to ≤4 µg/ml, voriconazole >0.125 to ≤0.5 µg/ml, and caspofungin, micafungin and anidulafungin >0.25 to ≤0.5 µg/ml. These exclusions are inferred from recorded MICs. They use the next-higher-category rule for values between dilution categories in [CLSI M27M44S, third edition (2022), Table 1](https://clsi.org/media/osthhxax/m27m44sed3e_sample.pdf).

Two analyses separate sample exclusion from gene-marker composition. `fixed` retains the original variants assigned to each gene and subsets the original carrier indicators. `rebuilt` reapplies the ≤5% variant-carrier criterion in the reduced drug-specific cohort. Each analysis tests genes with both carriers and called noncarriers and updates the multiple-testing thresholds using its tested gene count. Both retain the original log₂ MIC values and full-cohort-standardized nuclear GRMs, with LOCO adjustment for nuclear genes and the full nuclear GRM for mitochondrial genes. Variance components are refitted for each reduced cohort.

```bash
python3 07_run_category_sensitivity.py prepare \
  --drugs fluconazole voriconazole caspofungin micafungin anidulafungin \
  --original-root path/to/original_gene_GWAS \
  --mic-metadata path/to/all_metadata.tsv \
  --master-prefix path/to/cohort1153.snps_indels.master \
  --catalog path/to/selected_gene_variants.tsv.gz \
  --gene-map path/to/gene_map.tsv \
  --loco-manifest path/to/loco_manifest.tsv \
  --gcta path/to/gcta64 --threads 2 \
  --outdir path/to/gene_GWAS_category_sensitivity
python3 07_run_category_sensitivity.py run \
  --outdir path/to/gene_GWAS_category_sensitivity
python3 08_compare_category_sensitivity.py \
  --source-root path/to/original_gene_GWAS \
  --analysis-root path/to/gene_GWAS_category_sensitivity
```

The runner first reproduces the original fluconazole scan as a control, then writes both sensitivity analyses to separate directories. The comparison exports gene ranks, P values, carrier effects and counts, top-20 overlap, and changes in Bonferroni and suggestive classifications. This assesses robustness to outcome-based sample exclusion in the quantitative GWAS. Posaconazole retains its original cohort as the reference outside these selected exclusion ranges.

## Included input data

`input_data.zip` contains gene-association values, highlight annotations and reference chromosome lengths for the Manhattan and QQ panels. Run from this directory:

```bash
unzip input_data.zip
python3 06_plot_gene_gwas.py \
  --results gene_results.tsv.gz \
  --gene-evidence gene_evidence.tsv \
  --fai reference.fasta.fai --outdir plots
```

Outputs are individual association panels.
