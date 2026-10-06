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

Bonferroni and suggestive thresholds use 0.05/M and 1/M, with M the number of tested genes across all chromosomes for each drug. Lambda GC is `median((b/se)^2)/0.45493642311957275`. No P-value rescaling is performed.

```bash
python3 06_plot_gene_gwas.py \
  --results path/to/gene_results/all_gene_level_variant_results.tsv.gz \
  --gene-evidence path/to/reviewed_gene_evidence.tsv \
  --outdir path/to/gene_plots
```

The optional evidence table contains `drug`, `gene_id`, `display_symbol` and `evidence_level`. Categories `direct` and `response` highlight previously supported genes among the current top 20 for each drug, independently of statistical significance. Every top-20 membership is retained in exported source data.

References: [SnpEff](https://doi.org/10.4161/fly.19695), [GCTA](https://doi.org/10.1016/j.ajhg.2010.11.011), [LOCO mixed models](https://doi.org/10.1038/ng.2876), [genomic control](https://doi.org/10.1111/j.0006-341X.1999.00997.x).

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
