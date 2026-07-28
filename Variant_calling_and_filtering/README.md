# Variant calling and cohort filtering

The scripts reproduce the GATK and BCFtools stages described in the Methods. Paths are passed through environment variables so the workflow can run on local or cluster storage without embedded workstation paths.

Required software: GATK 4.1.4.1, Picard, BCFtools/HTSlib 1.21, SAMtools 1.21, and Slurm for the array-job wrappers.

Expected configuration:

```bash
export REFERENCE_FASTA=path/to/reference.fasta
export BAM_MANIFEST=path/to/bam_manifest.tsv
export GVCF_DIR=path/to/gvcfs
export SAMPLE_MAP=path/to/sample_map.tsv
export JOINT_DIR=path/to/joint_genotyping
export JOINT_VCF=path/to/cohort.vcf.gz
export SAMPLE_KEEP=path/to/sample_keep.txt
export FILTERED_DIR=path/to/filtered_variants
```

`bam_manifest.tsv` contains a header followed by `sample` and `bam` columns. `sample_map.tsv` contains tab-delimited `sample` and `gvcf` columns without a header, as required by `GenomicsDBImport`.

## Numbered Entry Points

1. `01_haplotypecaller_array.sh`: call chromosome-level gVCFs for diploid isolates.
2. `02_gather_sample_gvcfs.sh`: gather chromosome shards into one gVCF per isolate.
3. `03_joint_genotype_array.sh`: import the cohort into GenomicsDB and jointly genotype each chromosome.
4. `04_concatenate_joint_vcfs.sh`: concatenate chromosome-level cohort VCFs in reference order.
5. `05_filter_cohort_vcf.sh`: subset the final cohort, trim absent alleles, decompose multiallelic records, recalculate cohort tags, and retain high-confidence biallelic SNPs.

The public filtering workflow decomposes multiallelic records and does not apply genotype-quality masking. The manuscript-consistent output contains 1,093,947 biallelic records at 1,053,142 distinct genomic positions.
