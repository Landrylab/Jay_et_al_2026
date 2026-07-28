# Reproducibility Scripts

This directory contains compact script bundles for variant processing, phylogenetic inference, population-structure analysis, LOH detection, ploidy and copy-number analysis, WGS machine learning, and GCTA MLMA-LOCO GWAS. The numbered scripts are the public entry points. The `helpers/` folders contain supporting functions.

Input paths are supplied through command-line arguments or environment variables. Example placeholders use the following generic file names:

- `path/to/cohort.vcf.gz`
- `path/to/reference.fasta`
- `path/to/metadata.tsv`
- `path/to/alignment.phy`

No external database or connector access is required. Command-line tools such as `plink`, `bcftools`, `gcta64` or `gcta`, `Rscript`, and Python scientific packages must be available in the active environment. The LOH workflow uses pysam/HTSlib, NumPy and pandas. The CNV workflow uses deepTools, SAMtools, BEDTools, R/DNAcopy, NumPy and pandas.

## Folder Layout

- `WGS_machine_learning/`: the reported fluconazole binary R/S analyses only—SNP, PC and SNP+PC models under stratified or K=23 group-stratified cross-validation; same-cluster paired random forests; shuffled-label controls; group-stratified paired-RF feature ranking; and regional carrier support.
- `GWAS/`: PLINK conversion, GCTA MLMA-LOCO runs, significant-hit annotation, GWAS plots, and publication tables for the MAF 0.05 and MAF 0.005 analyses.
- `Variant_calling_and_filtering/`: per-chromosome GATK gVCF calling, joint genotyping, and cohort-level SNP filtering.
- `Phylogeny_and_population_structure/`: VCF-to-PHYLIP conversion, IQ-TREE inference, MAF and LD filtering, PCA, clustering, and submission-table exports.
- `LOH/`: fixed-window heterozygous-SNP counting, nuclear LOH tract detection, phylogenetic-neighbor sharing, clade-level region resolution and submission-table exports.
- `CNV/`: GC- and chromosome-end-bias correction, joint read-depth and allele-fraction segmentation, baseline-ploidy and integer-copy-number assignment, final manual-call reconciliation, diploid chromosome-event summaries and Source Data export.

No absolute personal paths are embedded in the shareable scripts. Repository-specific defaults may be overridden with the documented environment variables.

Run `./validate_release.sh` before publication to check shell, Python, and R syntax.
