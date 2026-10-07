# Reproducibility Scripts

This directory contains script bundles for paired-end read processing, variant processing, phylogenetic inference, population-structure analysis, LOH detection, ploidy and copy-number analysis, variant and gene-level mixed-model association, and WGS machine learning. The numbered scripts are the public entry points. The `helpers/` folders contain supporting functions.

Input paths are supplied through command-line arguments or environment variables. Example placeholders use the following generic file names:

- `path/to/cohort.vcf.gz`
- `path/to/reference.fasta`
- `path/to/metadata.tsv`
- `path/to/alignment.phy`

No external database or connector access is required. Command-line tools such as `plink`, `bcftools`, `gcta64` or `gcta`, `Rscript`, and Python scientific packages must be available in the active environment. The LOH workflow uses pysam/HTSlib, NumPy and pandas. The CNV workflow uses deepTools, SAMtools, BEDTools, R/DNAcopy, NumPy and pandas.

## Folder Layout

- `01_Genome_processing/`: paired-end trimming with Trimmomatic or fastp, read-quality reports, reference preparation, alignment with read groups, duplicate removal and read-depth measurement.
- `02_Variant_calling_and_filtering/`: per-chromosome GATK gVCF calling, joint genotyping, cohort-level SNP filtering, and sequence-indel extraction with explicit allele and genotype representation.
- `03_Phylogeny_and_population_structure/`: VCF-to-PHYLIP conversion, IQ-TREE inference, MAF and LD filtering, PCA, clustering, and submission-table exports.
- `04_LOH/`: fixed-window heterozygous-SNP counting, nuclear LOH tract detection, phylogenetic-neighbor sharing, clade-level region resolution and submission-table exports.
- `05_CNV/`: GC- and chromosome-end-bias correction, joint read-depth and allele-fraction segmentation, baseline-ploidy and integer-copy-number assignment, final manual-call reconciliation, diploid chromosome-event summaries and Source Data export.
- `06_Variant_GWAS/`: drug-specific filtering of the shared SNP/indel master, a common nuclear SNP panel shared with population PCs, full and chromosome-excluded relatedness matrices, additive log2-MIC association with GCTA MLMA, multiple-testing thresholds, Manhattan plots and Q–Q diagnostics.
- `07_Gene_level_variant_GWAS/`: rare-variant carrier indicators from the shared annotated master, quantitative MIC association using the shared relatedness matrices, carrier-effect summaries, and nuclear and mitochondrial gene-level Manhattan and Q–Q plots.
- `08_WGS_machine_learning/`: variant, gene-level variant and training-derived PC inputs; random forests with inner tuning and threshold selection under stratified or genetic-cluster-grouped validation; within-group gene importance; prediction and genomic-importance plots.
- `09_FungAMR/`: allele-specific matching of the standardized master to curated resistance substitutions, drug-matched carrier classifications, and Figure 2b source tables.
- `10_VEP_FungAMR_ML/`: the Figure 6 notebook and accompanying VEP/FungAMR input files.

No absolute personal paths are embedded in the shareable scripts. Repository-specific defaults may be overridden with the documented environment variables.

## Included input data

Folders 01–09 provide `input_data.zip` archives containing reference files or prepared numeric inputs. Extract an archive beside its scripts and run the commands in that folder's README from the same folder. Input filenames are relative to that working folder. Folder 10 retains its original input archives. Sequencing reads, BAMs and full genotype matrices are supplied separately for upstream analyses.

## Figures from numeric source data

The plotting entry points in the population-structure, copy-number and machine-learning folders accept named Excel tables, individual TSV/CSV files, or a directory of delimited tables. Each folder documents the table names, required columns and plotting options. Sample identifiers, cohort sizes, group labels and measurements are read from the supplied data. No analysis results are required to match a predetermined study size.

Install the plotting dependencies with `python3 -m pip install -r plotting_requirements.txt`. The shared `figure_sources.py` reader selects Excel tables by name independently of worksheet layout. In a table directory, each file uses its table name as the filename followed by `.tsv`, `.csv`, `.tsv.gz` or `.csv.gz`. Delimited identifiers remain text, including leading zeroes; plotting scripts explicitly convert measurement fields. Output paths are supplied on the command line. The figure scripts consume numeric data.
