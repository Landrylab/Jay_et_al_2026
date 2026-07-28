# Ploidy and copy-number analysis

This workflow corrects GC and chromosome-end coverage trends, integrates read depth with heterozygous-SNP allele fractions, assigns baseline ploidy and integer copy number, applies recorded manual-review corrections, summarizes chromosome and gene copy number, and exports the 10 kb figure track and compact submission tables.

## Inputs

- `path/to/alignment.bam`: coordinate-sorted isolate alignment.
- `path/to/reference.2bit`: SC5314 Assembly 22 reference in two-bit format.
- `path/to/reference.fasta.fai`: reference index.
- `path/to/mask.bed`: rDNA, centromere, telomere and repeat mask.
- `path/to/heterozygous_snps.tsv`: biallelic heterozygous SNP allele counts.
- `path/to/genes.bed`: gene intervals.
- `path/to/final_ploidy.tsv`, `path/to/final_chromosome_cn.tsv` and `path/to/final_gene_cn.tsv`: final manually reviewed calls.

## Run order

```bash
bash 01_correct_gc_bias.sh \
  path/to/alignment.bam \
  path/to/reference.2bit \
  14500000 \
  path/to/gc_corrected.bam

bash 02_summarize_1kb_coverage.sh \
  path/to/gc_corrected.bam \
  path/to/reference.fasta.fai \
  path/to/coverage_1kb.tsv

Rscript 03_correct_chromosome_end_bias.R \
  --input path/to/coverage_1kb.tsv \
  --mask path/to/mask.bed \
  --output path/to/coverage_1kb.end_corrected.tsv

Rscript 04_segment_and_assign_copy_number.R \
  --depth path/to/coverage_1kb.end_corrected.tsv \
  --baf path/to/heterozygous_snps.tsv \
  --genes path/to/genes.bed \
  --output-prefix path/to/sample

python3 05_aggregate_10kb_copy_number.py \
  --depth path/to/coverage_1kb.end_corrected.tsv \
  --segments path/to/sample.segments.tsv \
  --baseline-ploidy 2 \
  --output path/to/copy_number_10kb/sample.tsv

python3 06_summarize_aneuploidy.py \
  --gene-cnv path/to/final_gene_cn.tsv \
  --ploidy path/to/final_ploidy.tsv \
  --chromosome-cn path/to/final_chromosome_cn.tsv \
  --loci path/to/gene_loci.tsv \
  --output-dir path/to/cnv_summary

python3 07_export_submission_data.py \
  --bin-cn-dir path/to/copy_number_10kb \
  --ploidy path/to/final_ploidy.tsv \
  --chromosome-cn path/to/final_chromosome_cn.tsv \
  --gene-cnv path/to/final_gene_cn.tsv \
  --loci path/to/gene_loci.tsv \
  --manual-review path/to/manual_reviewed_samples.tsv \
  --output-dir path/to/submission_tables \
  --source-data path/to/Figure_1G_CNV_10kb_windows.tsv.gz
```

The end-bias correction uses 10 folds and LOESS span 0.3. Bins require at least 75% covered bases and median matched-base depth of at least 80. The 10 kb aggregation requires directional support from at least 70% of contributing 1 kb bins. Segmental events in the diploid-isolate summary span at least 30 kb.

## References

- Ramírez F et al. *Nucleic Acids Research* (2016), deepTools2. DOI: [10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257).
- Benjamini Y and Speed TP. *Nucleic Acids Research* (2012), GC-bias correction. DOI: [10.1093/nar/gks001](https://doi.org/10.1093/nar/gks001).
- Li H et al. *Bioinformatics* (2009), SAMtools. DOI: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352).
- Quinlan AR and Hall IM. *Bioinformatics* (2010), BEDTools. DOI: [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033).
- Cleveland WS. *Journal of the American Statistical Association* (1979), locally weighted regression. DOI: [10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038).
- Olshen AB et al. *Biostatistics* (2004), circular binary segmentation. DOI: [10.1093/biostatistics/kxh008](https://doi.org/10.1093/biostatistics/kxh008).
- Abbey DA et al. *Genome Medicine* (2014), joint copy-number and allele-fraction visualization in *Candida*. DOI: [10.1186/s13073-014-0100-8](https://doi.org/10.1186/s13073-014-0100-8).
- Skrzypek MS et al. *Nucleic Acids Research* (2017), Assembly 22 and the Candida Genome Database. DOI: [10.1093/nar/gkw924](https://doi.org/10.1093/nar/gkw924).
