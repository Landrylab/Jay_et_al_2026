# Variant GWAS of quantitative antifungal MIC

All inputs derive from the normalized, annotated SNP/indel master produced in `02_Variant_calling_and_filtering`. Keep the numbered workflow directories together. Scripts use Python, NumPy, pandas, SciPy, PLINK 1.9 and GCTA 1.95.0. All paths below are placeholders.

## Genotypes and phenotypes

The master supplies `.dosage.npy` (alleles × samples; signed int8), `.samples.tsv`, and `.variants.tsv.gz`. Dosages count copies of the specified ALT: 0, 1 or 2, with -1 denoting missing. At multiallelic sites, zero means no copies of the target ALT and can include another ALT. BED A1 is always the target ALT; its A2 column records REF for identification, while the genotype encoding represents any called non-target allele.

Provide one headerless, three-column phenotype file per drug, named `pheno.mic.DRUG.txt`, containing FID, IID and numeric **log2 MIC**. FID and IID must match. Missing values are absent, nonfinite or -9. Values are already transformed and are not logged again. Supported drugs are fluconazole, voriconazole, posaconazole, caspofungin, micafungin and anidulafungin.

## Shared common-SNP panel and relatedness

```bash
python3 00_build_common_snp_panel.py \
  --master-prefix path/to/cohort.snps_indels.master \
  --out-prefix path/to/common_nuclear_ld_pruned --plink path/to/plink

python3 02_build_relationship_matrix.py \
  --master-prefix path/to/cohort.snps_indels.master \
  --shared-panel-prefix path/to/common_nuclear_ld_pruned \
  --outdir path/to/GRM
```

The shared panel contains nuclear SNPs at MAF ≥0.05 and missingness ≤10% across all master samples, pruned with `--indep-pairwise 50 5 0.2`. Its BED/BIM/FAM files, retained IDs and master row indices define the exact same SNPs used for GWAS relatedness and PC-based RF. Chromosomes 1–7 and 9 (ChrR) are nuclear; chromosome 8 is mitochondrial. No mitochondrial allele contributes to the GRM or PCA.

For each retained SNP, `Z=(dosage−2p)/sqrt(2p(1−p))`, with p calculated in the complete master cohort and missing Z set to zero. The full matrix is `K=Z.T@Z/M`. Chromosome-excluded matrices subtract each chromosome's cross-product and divide by the remaining SNP count. No additional markers are removed during GRM construction. `loco_manifest.tsv` maps nuclear chromosomes to their excluded matrices and mitochondrial markers to the full nuclear matrix for the gene-level analysis.

## Association inputs and fitting

```bash
python3 01_prepare_inputs.py \
  --master-prefix path/to/cohort.snps_indels.master \
  --phenotype-dir path/to/phenotypes --outdir path/to/variant_inputs

python3 03_run_variant_gwas.py \
  --prepared-manifest path/to/variant_inputs/prepared_scan_manifest.tsv \
  --grm path/to/GRM/common_nuclear --gcta path/to/gcta64 \
  --threads 2 --outdir path/to/variant_associations

python3 04_summarize_results.py \
  --manifest path/to/variant_associations/scan_manifest.tsv \
  --outdir path/to/variant_results
```

Eligibility is calculated within each drug cohort: MAF ≥0.005, missing genotypes ≤10%, and at least five minor-allele carriers. Frequency uses called diploid alleles; minor-allele carriers have at least one copy of whichever target-versus-nontarget allele is less frequent. There is no minimum count outside the modal genotype class. A dosage constant across all called samples has no estimable additive effect and is retained in a separate audit table. Both nuclear and mitochondrial sequence SNPs and indels are eligible. Association markers are not LD pruned.

GCTA MLMA fits an intercept, an additive target-ALT effect and a random genetic effect using the supplied full nuclear GRM. No PC covariates are included. Null-model variance components are estimated by REML; missing association genotypes are mean-imputed within the phenotype cohort. Scripts check REML convergence, sample ordering, allele identity, test completeness and result validity.

## Thresholds and plots

A two-sided association P value is obtained from `(effect/SE)^2` under a one-degree-of-freedom chi-squared null. Within each drug, Bonferroni and suggestive cutoffs are 0.05/M and 1/M over all tested nuclear and mitochondrial SNP/indel alleles. QQ expected probabilities use `(rank−0.5)/M`. Lambda GC is the median P-derived chi-squared statistic divided by its theoretical null median (0.454936).

```bash
python3 05_plot_variant_associations.py \
  --manifest path/to/variant_results/primary_scan_manifest.tsv \
  --fai path/to/reference.fasta.fai --outdir path/to/variant_plots
```

Relatedness and mixed-model methods follow [Yang et al. 2011](https://doi.org/10.1016/j.ajhg.2010.11.011) and [Yang et al. 2014](https://doi.org/10.1038/ng.2876); genomic control follows [Devlin and Roeder 1999](https://doi.org/10.1111/j.0006-341X.1999.00997.x).

## Included input data

`input_data.zip` contains variant-association values, the scan manifest and reference chromosome lengths for the Manhattan and QQ panels. Run from this directory:

```bash
unzip input_data.zip
python3 05_plot_variant_associations.py \
  --manifest primary_scan_manifest.tsv \
  --fai reference.fasta.fai --outdir plots
```

Use a new output directory. Outputs are individual association panels.
