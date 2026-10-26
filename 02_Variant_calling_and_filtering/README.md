# Variant calling and cohort filtering

The scripts perform GATK variant discovery, joint genotyping, cohort-level SNP filtering and sequence-indel extraction. Paths are supplied through environment variables or command-line arguments. Prepared BAMs can be produced with the sibling `01_Genome_processing/` workflow.

Required software: GATK 4.1.4.1, Picard, BCFtools/HTSlib 1.21, SAMtools 1.21, and Slurm for the array-job wrappers.

## Included input data

`input_data.zip` contains the SC5314 Assembly 22 haplotype-A reference FASTA, its matching release s07-m01-r183 GFF annotation, the cohort sample list in master-column order, and the reference `gene_map.tsv` used by downstream analyses. From this directory:

```bash
unzip input_data.zip
samtools faidx reference.fasta
samtools dict -o reference.dict reference.fasta
```

The files extract beside the scripts. Use `reference.gff3` for annotation and `sample_keep.txt` for cohort selection. BAMs, gVCFs, cohort VCFs and their local file manifests are supplied separately.

Expected configuration:

```bash
export REFERENCE_FASTA=reference.fasta
export BAM_MANIFEST=path/to/bam_manifest.tsv
export GVCF_DIR=path/to/gvcfs
export SAMPLE_MAP=path/to/sample_map.tsv
export JOINT_DIR=path/to/joint_genotyping
export JOINT_VCF=path/to/cohort.vcf.gz
export SAMPLE_KEEP=sample_keep.txt
export FILTERED_DIR=path/to/filtered_variants
```

`bam_manifest.tsv` contains a header followed by `sample` and `bam` columns. `sample_map.tsv` contains tab-delimited `sample` and `gvcf` columns without a header, as required by `GenomicsDBImport`.

The calling, gathering, joint-genotyping and concatenation scripts use the exact ordered contig names from `${REFERENCE_FASTA}.fai`. `REFERENCE_FAI` can select another FASTA index. When no index is available, the scripts read `@SQ` records from the reference sequence dictionary (`reference.dict` beside `reference.fasta`), or from `REFERENCE_DICT`. All four steps must use the same reference configuration. An explicit whitespace-separated `CHROMOSOMES` list overrides the index or dictionary and defines the order of the selected shards. Contig names must be unique and usable as filename components. Chromosome-array indices are one-based; the HaplotypeCaller array contains one task per sample and selected contig, while the joint-genotyping array contains one task per selected contig.

## Numbered Entry Points

1. `01_haplotypecaller_array.sh`: call chromosome-level gVCFs for diploid isolates.
2. `02_gather_sample_gvcfs.sh`: gather chromosome shards into one gVCF per isolate.
3. `03_joint_genotype_array.sh`: import the cohort into GenomicsDB and jointly genotype each chromosome.
4. `04_concatenate_joint_vcfs.sh`: concatenate chromosome-level cohort VCFs in reference order.
5. `05_filter_cohort_vcf.sh`: subset the final cohort, trim absent alleles, decompose multiallelic records, recalculate cohort tags, and retain high-confidence biallelic SNPs.
6. `06_extract_sequence_indels.py`: split sequence-indel ALT alleles without realignment, preserve parent-record identifiers, and export allele-specific PLINK genotypes, an allele manifest and an initial count summary. Sample order can follow the VCF, a sample list or a supplied FAM.
7. `07_build_master_variants.py`: normalize the shared SNP/indel catalogue and retain original records, sample order and target-ALT dosages.
8. `08_annotate_master_variants.py`: annotate the normalized catalogue and prepare shared allele-specific gene annotations.
9. `09_audit_master_dosages.py`: independently compare deterministic sampled dosage rows with their original genotypes and source-record mappings.

## Variant representations

The SNP filtering workflow decomposes multiallelic records and recalculates AC, AN, AF, NS, F_MISSING and MAF after sample selection. It retains segregating biallelic SNPs with at most 10% genotype missingness and mapping quality of at least 20 or an absent MQ value. It does not apply genotype-quality masking. Its output is a SNP panel; use the joint variant source, before SNP-only selection, for sequence-indel extraction.

Indels are sequence alleles whose REF and ALT lengths differ. The extraction command uses `bcftools norm -m -any -N --multi-overlaps 0`; coordinates and REF/ALT representations are retained without reference-driven left alignment. Symbolic alleles, spanning-deletion `*` alleles and equal-length substitutions are excluded from this output. Distinct representations in repetitive sequence are not collapsed. Repeated exact chromosome/position/REF/ALT keys are accepted only when their selected-sample genotypes agree.

For each extracted allele, complete diploid calls encode the number of copies of that target ALT. Other called ALT alleles contribute zero copies of the target. Partial missing and non-diploid calls are missing. The original parent record and ALT index remain available in the allele manifest. Optional `--exclude-bim` skips exact allele keys listed in another panel, without recoding that panel's genotypes.

The indel extractor inventories the available sequence alleles without imposing missingness, allele-frequency or carrier-count thresholds. Selecting samples does not remove alleles that have no carriers in that subset. The SNP filtering summary and indel summary distinguish allele records from genomic positions.

```bash
python3 06_extract_sequence_indels.py \
  --vcf path/to/joint_calls.vcf.gz \
  --out-prefix path/to/indels/sequence_indels
```

By default, all VCF samples are retained in VCF order. Use `--samples path/to/sample_keep.txt` for one sample ID per line in the desired order, or `--sample-fam path/to/cohort.fam` to retain the exact identifiers and metadata of an existing six-column FAM. These options are mutually exclusive. A headerless, two-column `--chromosome-map` can translate contig names for the BIM while the allele manifest preserves original names. Output files are not overwritten.

The extractor requires Python 3 and BCFtools with `norm --multi-overlaps` support.

See the [BCFtools normalization documentation](https://samtools.github.io/bcftools/bcftools.html#norm) for multiallelic splitting and overlapping-allele handling. Run the indel script with `--help` for its input and output options.

## Shared SNP and indel catalogue

`07_build_master_variants.py` creates an analysis catalogue directly from the complete joint VCF and a cohort sample list. It retains sequence SNPs and unequal-length sequence indels, including eligible ALT alleles at mixed and multiallelic sites. Symbolic alleles, spanning-deletion `*` alleles, equal-length multinucleotide substitutions and reference-identical alleles are excluded before normalization. No minor-allele-frequency, carrier-count, missingness, mapping-quality or VCF FILTER threshold is applied to the master catalogue. Nuclear and mitochondrial alleles are retained. Complex multinucleotide substitutions are not atomized into SNPs.

BCFtools checks REF alleles against the supplied FASTA and left-aligns and minimally represents indels. Stable identifiers contain the reference contig, position, REF and ALT. Canonical duplicates are collapsed only when all selected-isolate target-ALT dosages agree; discordant duplicates terminate the workflow. A separate table maps every retained source allele to its canonical allele, and the original multiallelic genotypes are retained in a compact BCF. The original-record table records the source-record ordinal, original coordinate/alleles and target ALT index.

Complete diploid genotypes are encoded as 0, 1 or 2 copies of the specified ALT. Any partial missing, fully missing or non-diploid call is represented as missing (-1 in the matrix). Other called ALT alleles contribute zero target copies. For example, original genotype `1/2` yields dosage 1 for ALT1 and dosage 1 for ALT2. Original `1/.` remains missing. Dosage zero therefore means that the target ALT is absent, not necessarily that the original genotype was homozygous reference. Phasing remains available in the original-genotype BCF; analysis dosages are unphased.

```bash
python3 07_build_master_variants.py \
  --vcf path/to/joint_calls.vcf.gz \
  --samples sample_keep.txt \
  --reference reference.fasta \
  --out-prefix path/to/variants/cohort.snps_indels.master
```

The sample list contains one unique identifier per nonempty line and defines the exact output column order. The selected source BCF and normalized intermediates must preserve this order. Each output prefix has its own `<prefix>.master_work/` directory. Its manifest records the source VCF, reference and sample-list hashes, ordered sample identifiers, BCFtools version and artifact checksums. A fresh run requires an unused output prefix.

Use the same command with `--resume` to validate and reuse a completed preparation stage. Every cached artifact and its sample order are checked before downstream output is produced. Preparation must have finished successfully and written its manifest; an incomplete preparation cannot be resumed. A completed catalogue is checked and returned without rewriting its files or annotation audit. Changed inputs, altered cached artifacts, missing outputs or incompatible manifests require a new output prefix. Each prefix is locked while in use.

`08_annotate_master_variants.py` builds a SnpEff database from the supplied reference and GFF, restricted to contigs in the reference, and annotates the normalized allele catalogue. Candida albicans nuclear chromosomes use translation table 12, while its mitochondrial chromosome uses table 4. Gene and transcript biotypes are explicitly encoded from the reference feature types; nuclear tRNAs use SnpEff’s non_coding category and mitochondrial tRNAs use Mt_tRNA. Pseudogene identifiers and their noncoding biotype are retained. Independent CDS and protein references can be supplied for database validation; their presence or absence is recorded in the annotation audit. Full allele-specific ANN records are retained. For explicitly noncoding reference features, protein-only consequences are represented as noncoding-transcript exon effects with MODIFIER impact and protein HGVS fields are empty; legitimate HIGH transcript-ablation consequences remain eligible. Raw SnpEff output and exact record-level standardization are retained for traceability. One HIGH/MODERATE annotation per allele is selected for shared gene-level inputs using this priority: HIGH before MODERATE, membership in the reference gene map, nonempty protein HGVS before empty HGVS, then gene ID, effect, protein HGVS and transcript ID in lexical order. Shared gene assignment precedes any analysis-specific carrier-frequency filtering.

```bash
python3 08_annotate_master_variants.py \
  --out-prefix path/to/variants/cohort.snps_indels.master \
  --reference reference.fasta \
  --gff reference.gff3 \
  --snpeff-jar path/to/snpEff.jar \
  --base-config path/to/snpEff.config \
  --java path/to/java
```

The final annotated VCF includes recalculated AC, AN, AF, NS, F_MISSING and MAF. Its accompanying dosage NPY has **variants as rows and isolates as columns** and int8 values -1/0/1/2. Row order follows `*.variants.tsv.gz`; column order follows `*.samples.tsv`. `selected_gene_variants.tsv.gz` and `gene_map.tsv` provide a common annotation interface for gene-level association and random-forest inputs. Downstream analyses calculate their own cohort-specific eligibility rules from these shared files. The SNP panel used for phylogenetic and population-structure figures remains a separate output of step 05.

Dependencies for the shared catalogue are Python with NumPy and cyvcf2, BCFtools/HTSlib, Java and SnpEff.

## Audit master dosages

```bash
python3 09_audit_master_dosages.py \
  --prefix path/to/variants/cohort.snps_indels.master \
  --out results/master_dosage_validation.json
```

The audit requires Python with NumPy and BCFtools. It uses the prefix's dosage matrix, sample table, variant metadata, normalization map, original-record table and indexed original-genotype BCF. `--original-vcf` can select an indexed original VCF/BCF with the same record and sample order; `--bcftools` selects the executable.

By default, 4,000 rows are sampled without replacement using seed 20260923, or all rows when the catalogue is smaller. `--rows` and `--seed` configure this check. The audit verifies matrix dimensions and sample order, sampled allele counts and frequencies, original record ordinals, every original mapping for selected canonical alleles, and target-ALT dosages including partial missing and non-diploid calls. Duplicate source representations are checked separately. It records source-file checksums and exits with failure for a disagreement. A successful sampled audit does not assert that every unsampled genotype was checked. Source data are read without modification.
