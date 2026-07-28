# GWAS Scripts

This folder contains a compact reproducibility bundle for the GCTA MLMA-LOCO GWAS analyses.

## Numbered Entry Points

1. `01_build_plink_from_joint_vcf.sh`
   - Converts `joint.keep1153.snps.setid.renamed.vcf.gz` into a full PLINK bed/bim/fam dataset.
   - Variant IDs are normalized to `CHR:POS:REF:ALT`.

2. `02_run_gcta_mlma_loco.sh`
   - Runs GCTA `--mlma-loco` for six antifungals and both MIC and binary R/S phenotypes.
   - Uses MAF filtering, an LD-pruned GRM for PC1-PC10, and per-trait plotting.
   - GCTA constructs the chromosome-specific LOCO relationship matrices from the full MAF-filtered association panel.
   - Run separately with `MAF=0.05` and `MAF=0.005` for the publication tables.

3. `03_annotate_and_render_gwas_hits.py`
   - Annotates Bonferroni-significant hits and renders labelled Manhattan/LD-context outputs.

4. `04_build_gwas_publication_tables.py`
   - Builds the compact supplementary GWAS tables using the MAF 0.05 and MAF 0.005 results.

5. `05_compile_final_gwas_figures.R`
   - Rebuilds the final GWAS figure assembly when all component plots are present.

6. `06_render_final_gwas_figures_lowmaf_ld_black_points.py`
   - Renders the final monochrome GWAS figure set and marks SNPs in linkage disequilibrium with high-confidence FungAMR mutations, including mutations below the GWAS MAF threshold.

7. `07_export_submission_data.py`
   - Exports compact GWAS tables and four compressed point tables containing the values used for the Manhattan and quantile-quantile plots.

## Typical Sequence

```bash
bash 01_build_plink_from_joint_vcf.sh
MAF=0.05 bash 02_run_gcta_mlma_loco.sh
MAF=0.005 bash 02_run_gcta_mlma_loco.sh
python3 03_annotate_and_render_gwas_hits.py
python3 04_build_gwas_publication_tables.py \
  --results-root path/to/gcta_mlma_loco \
  --figure-root path/to/final_gwas_paper_figures_lowmaf_fungamr_ld_black_points \
  --submission-root path/to/Nature_submission
Rscript 05_compile_final_gwas_figures.R
python3 06_render_final_gwas_figures_lowmaf_ld_black_points.py
python3 07_export_submission_data.py \
  --results-root path/to/gcta_mlma_loco \
  --figure-root path/to/final_gwas_paper_figures_lowmaf_fungamr_ld_black_points \
  --submission-root path/to/Nature_submission
```

The helper scripts under `helpers/` preserve the core plotting, phenotype, annotation, and table-building logic used by the local analysis.

For step 06, local paths can be supplied with `GWAS_RESULTS_ROOT`, `GWAS_BASE_PLINK`, `FUNGAMR_MAPPING` and `GWAS_FIGURE_OUTDIR`.

## Methods implemented

- Variants are represented in PLINK by A/C/G/T allele pairs and tested under an additive 0/1/2 dosage model.
- Analyses are restricted to MAF 0.05 and 0.005.
- The LD-pruned set is generated with `--indep-pairwise 50 5 0.2` and used for PC1-PC10.
- GCTA `--mlma-loco` fits a linear mixed model and excludes the tested chromosome from the relationship term.
- Binary resistance is coded R=1, S=0 and missing=-9 and is analysed as a binary-coded outcome in the same linear mixed model.
- Bonferroni thresholds are `0.05 / N`, where `N` is the number of valid SNP tests in each analysis.
- Lambda GC is the median one-degree-of-freedom chi-square statistic divided by its null median.
- FungAMR highlighting uses confidence scores 1 or 2 and PLINK LD of r2 at least 0.2 within 1 Mb.

## Primary references

- GCTA: Yang et al. (2011), https://doi.org/10.1016/j.ajhg.2010.11.011
- Mixed-model association and LOCO: Yang et al. (2014), https://doi.org/10.1038/ng.2876
- PLINK: Chang et al. (2015), https://doi.org/10.1186/s13742-015-0047-8
- BCFtools: Danecek et al. (2021), https://doi.org/10.1093/gigascience/giab008
- Genomic control: Devlin and Roeder (1999), https://doi.org/10.1111/j.0006-341X.1999.00997.x
- Linkage disequilibrium r2: Hill and Robertson (1968), https://doi.org/10.1007/BF01245622
- SnpEff: Cingolani et al. (2012), https://doi.org/10.4161/fly.19695
- FungAMR: https://doi.org/10.1038/s41564-025-02084-7
