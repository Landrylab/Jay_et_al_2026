# Ploidy and copy-number analysis

This workflow summarizes nuclear read depth, corrects GC and chromosome-end trends, segments depth and allele-fraction tracks, and produces reviewed segment, chromosome, gene and 10-kb copy-number tables. All example paths are placeholders. Commands run from this directory; the two `04` scripts perform calling and diagnostic plotting, respectively.

## Dependencies and coordinates

Use deepTools (`computeGCBias`, `correctGCBias`), SAMtools, BEDTools, BCFtools, Python 3 with NumPy/pandas, and R with data.table and [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html). `faToTwoBit` converts the reference FASTA for deepTools. Record installed software versions with each analysis. Review plots use base R graphics.

Reference BED intervals, gene intervals, depth bins and segments use **zero-based, half-open** coordinates. Input SNP `pos` values use **one-based** VCF positions; the caller converts them internally before interval overlap. Gene-locus exports and 10-kb tracks use **one-based, inclusive** coordinates. Contig names must agree across BAM, VCF, FASTA, annotation and mask files.

## Reference and input preparation

```bash
samtools faidx path/to/reference.fasta
faToTwoBit path/to/reference.fasta path/to/reference.2bit
python3 00_prepare_reference.py \
  --fai path/to/reference.fasta.fai \
  --gff path/to/reference.gff3 \
  --mask-bed path/to/rdna_centromere_telomere_repeat_mask.bed \
  --output-dir path/to/reference_inputs

python3 02_extract_heterozygous_snp_counts.py \
  --vcf path/to/undecomposed_joint_calls.vcf.gz \
  --sample SAMPLE_ID \
  --output path/to/SAMPLE_ID.heterozygous_snps.tsv
```

Reference preparation selects nuclear contigs present in the FASTA index, converts matching GFF features to gene intervals, groups genes with identical coordinates into loci, and merges exclusion intervals. Configure `--gene-types` and `--mask-types` for the annotation vocabulary. The defaults recognize `gene,ORF,pseudogene` and `centromere,telomere,repeat_region,rRNA_gene`, respectively. Supply explicit BED intervals for regions not represented by these annotation types; the script does not infer missing rDNA, repeat or chromosome-end boundaries. The default mitochondrial exclusion is `chrM|mitochond` (case insensitive).

Allele-count preparation uses FORMAT/AD from the alignment-based joint variant calls. It requires a **biallelic SNP in the undecomposed record** and a called diploid heterozygote (`0/1`, `1/0`, or phased equivalent), with at least one read supporting each allele. It excludes multiallelic records, indels, partial calls, missing AD and repeated qualifying positions. Counts reflect the diploid genotyping callset even when the underlying isolate is polyploid; no SNP-frequency or LD pruning is applied. A companion JSON records retention and exclusion counts. Optional depth thresholds are configurable. This step extracts existing allele counts rather than recounting reads from the GC-resampled BAM.

## Coverage correction and candidate calls

```bash
bash 01_correct_gc_bias.sh \
  path/to/alignment.bam path/to/reference.2bit \
  EFFECTIVE_GENOME_SIZE path/to/SAMPLE_ID.gc_corrected.bam

bash 02_summarize_1kb_coverage.sh \
  path/to/SAMPLE_ID.gc_corrected.bam path/to/reference.fasta.fai \
  path/to/SAMPLE_ID.coverage_1kb.tsv

Rscript 03_correct_chromosome_end_bias.R \
  --input path/to/SAMPLE_ID.coverage_1kb.tsv \
  --mask path/to/reference_inputs/mask.bed \
  --output path/to/SAMPLE_ID.end_corrected.tsv

Rscript 04_segment_and_assign_copy_number.R \
  --depth path/to/SAMPLE_ID.end_corrected.tsv \
  --baf path/to/SAMPLE_ID.heterozygous_snps.tsv \
  --genes path/to/reference_inputs/genes.tsv \
  --output-prefix path/to/candidate/SAMPLE_ID

Rscript 04_plot_copy_number_review.R \
  --depth path/to/SAMPLE_ID.end_corrected.tsv \
  --baf path/to/SAMPLE_ID.heterozygous_snps.tsv \
  --prefix path/to/candidate/SAMPLE_ID \
  --output path/to/review/SAMPLE_ID.pdf
```

`EFFECTIVE_GENOME_SIZE` is a numeric value appropriate to the reference and mappability policy. Coverage tables retain read-overlap counts, median matched-base depth and covered-base fraction. Bins require at least 75% covered bases and median matched-base depth of at least 80, and must not overlap the mask. The chromosome-end correction models log(read-overlap count + 1) by distance to the chromosome ends, using quadratic LOESS (span 0.3) and ten interleaved held-out folds. Chromosome lengths come from the full coverage table before coverage filtering and masking.

The caller uses CBS separately for normalized depth and absolute allele-fraction deviation from 0.5 (alpha 0.01; minimum width three observations; seed 42). Internal change points are placed between adjacent segment boundary observations and combined across the two tracks. Constant or short tracks have no internal breakpoints. Consensus intervals summarize depth and SNP evidence; intervals without depth evidence remain uncalled.

Candidate baseline ploidy is scored for 2N, 3N, 4N and 6N by mean squared distance between observed minor-allele fractions and the nearest expected state. This is an **allele-fraction heuristic**, not a calibrated likelihood or an independent measurement of absolute DNA content. Higher-ploidy state sets overlap lower-ploidy state sets, and some profiles cannot distinguish these alternatives. All candidate scores and equally scoring ploidies are exported; ties select the first, lower-ploidy candidate for inspection. Baseline ploidy requires review of allele-fraction peaks, depth patterns, segmentation, and independent evidence where available. Relative read depth alone cannot distinguish genome-wide dosage multiples. Without informative SNPs, an explicit reviewed baseline is required.

For each baseline, candidate segment states 0–8 minimize squared error to expected log2 depth plus the nearest expected minor-allele-fraction error when at least ten SNPs support the interval. The zero-copy depth expectation uses a finite log2 floor of −4. The two error terms have equal numerical weights; no probabilistic calibration is implied. Chromosome CN is the length-weighted modal segment state; gene CN is the state of the segment with greatest overlap when it covers at least 50% of the gene. Sparse or uncalled intervals require inspection.

## Review and cohort tables

Record one review row per sample in a TSV with columns `sample_id`, `baseline_ploidy`, `reviewer`, and `evidence`. Evidence can identify the diagnostic plot and an independent ploidy measurement. The reviewer supplies the accepted baseline rather than relying on the candidate heuristic alone.

Run the caller with the accepted baseline to obtain internally consistent segment, chromosome and gene calls:

```bash
Rscript 04_segment_and_assign_copy_number.R \
  --depth path/to/SAMPLE_ID.end_corrected.tsv \
  --baf path/to/SAMPLE_ID.heterozygous_snps.tsv \
  --genes path/to/reference_inputs/genes.tsv \
  --baseline-ploidy 4 \
  --output-prefix path/to/reviewed/SAMPLE_ID
```

The number 4 illustrates a reviewed tetraploid baseline; use the reviewed value for the sample. Optional `--segment-overrides path/to/segment_review.tsv` accepts `segment_id`, `integer_copy_number`, `reviewer`, and `evidence`; use segment IDs from the same input data, parameters and seed. These explicit states are applied before deriving gene and chromosome summaries. A baseline override recalculates every segment's expected depth and copy-number score; it does not simply replace a summary value. Replot accepted calls for inspection.

Create a sample manifest with columns `sample_id`, `prefix`, and `depth`. `prefix` identifies the reviewed caller output prefix; `depth` identifies the corrected 1-kb table. Relative paths are resolved against the manifest directory.

```bash
python3 08_collate_reviewed_calls.py \
  --manifest path/to/sample_manifest.tsv \
  --review path/to/sample_review.tsv \
  --output-dir path/to/cohort

python3 06_summarize_aneuploidy.py \
  --gene-cnv path/to/cohort/gene_cn.tsv \
  --ploidy path/to/cohort/baseline_ploidy.tsv \
  --chromosome-cn path/to/cohort/chromosome_cn.tsv \
  --loci path/to/reference_inputs/gene_loci.tsv \
  --output-dir path/to/cnv_summary

python3 07_export_submission_data.py \
  --bin-cn-dir path/to/cohort/binCN_10kb \
  --ploidy path/to/cohort/baseline_ploidy.tsv \
  --chromosome-cn path/to/cohort/chromosome_cn.tsv \
  --gene-cnv path/to/cohort/gene_cn.tsv \
  --loci path/to/reference_inputs/gene_loci.tsv \
  --manual-review path/to/cohort/manual_review.tsv \
  --output-dir path/to/source_tables \
  --source-data path/to/CNV_10kb_windows.tsv.gz
```

Collation checks sample, gene and chromosome consistency and rejects a mismatch between reviewed baseline and caller output. It writes a sample-by-chromosome matrix, gene-by-sample matrix, baseline and review tables, and one 10-kb file per sample. These files contain `chr`, `start`, `end`, `CN` and `baseline_ploidy`, as consumed by `../04_LOH/02_call_loh_tracts.py`.

Ten-kb aggregation requires at least 70% of available, called 1-kb bins to support a gain or loss relative to the accepted baseline; other windows take the baseline state. Windows with no called bins are omitted rather than interpreted as copy neutral. Summary tables describe whole-chromosome events and consecutive altered gene-locus runs of at least 30 kb among diploid samples, excluding runs on a chromosome with a whole-chromosome event in the same sample. Source export preserves the reviewed, segment-derived 10-kb calls and records their support counts. `--expected-samples` is optional and can enforce a cohort size in summary/export steps.

## Ploidy-validation figure

`09_plot_ploidy_validation.py` draws one row per supplied isolate, with allele-fraction density, absolute copy number above minor-allele fraction, and a FACS histogram against a supplied reference recording. An optional final panel shows chromosome-specific aneuploidy composition. The command reads prepared numerical values. Install the plotting dependencies listed in `../plotting_requirements.txt`.

```bash
python3 09_plot_ploidy_validation.py \
  --source path/to/ploidy_source_data.xlsx \
  --reference-sample REFERENCE_SAMPLE \
  --output-prefix path/to/plots/ploidy_validation
```

`--source` accepts a workbook containing the named Excel Tables below, or a directory containing files named after those tables with `.tsv`, `.csv`, `.tsv.gz` or `.csv.gz` extensions. Each table can instead be supplied separately through `--samples`, `--densities`, `--genome`, `--facs` and `--composition`. Column names are case sensitive. Tables can be stacked on a worksheet; their positions and worksheet names do not affect loading.

| Named table | Required columns | Row definition |
| --- | --- | --- |
| `ED3Samples` | `Panel`, `Isolate`, `Assigned ploidy (N)` | One isolate and positive integer ploidy per unique panel. Optional `Display label` and `FACS sample` columns control labels and recording matching. |
| `ED3Densities` | `Panel`, `Allele fraction`, `REF density`, `ALT density` | One supplied curve coordinate per panel and allele fraction. Allele fractions span 0–1; densities are nonnegative and retain their original normalization and smoothing. |
| `ED3GenomeBins` | `Panel`, `Chromosome`, `Bin start (bp)`, `Bin end (bp)`, `Absolute copy number`, `Median minor AF` | One nonoverlapping genomic bin. Starts are zero-based and ends exclusive. CN is already on the assigned absolute-ploidy scale. Minor AF may be blank; otherwise it must lie in 0–0.5. |
| `ED3FACS` | `Sample`, `Bin start (a.u.)`, `Bin end (a.u.)`, `Events`, `Events (%)` | One fluorescence bin per recording. All recordings share contiguous bin edges. Counts are integers; percentages must equal 100 times the bin count divided by the sum over every supplied bin for that recording. |
| `ED3Composition` | `Chromosome`, `Aneuploidy type`, `Copy number`, `Isolates`, `Chromosome denominator`, `Aneuploid isolates (%)` | One chromosome/type/CN category. Type is `whole chromosome` or `segmental`; category counts sum to the chromosome denominator, and percentages use that denominator. |

Use `--no-composition` when the final panel is not needed. The FACS table must contain the identifier passed to `--reference-sample`. An explicit `FACS sample` value in the sample table takes precedence; a blank value marks an unavailable sample. Alternatively, `--facs-map path/to/mapping.tsv` supplies columns `Panel` and `FACS sample` for every panel. Without an explicit mapping, exact isolate identifiers are matched first, followed by the identifier prefix before the first underscore. Ambiguous matches and unused non-reference recordings are rejected. Unmatched isolates display “Sample not available for testing”; `--unavailable-text` changes that text.

Default axes include every supplied CN value and fluorescence bin. Optional `--cn-max`, `--facs-min`, `--facs-max` and `--facs-ymax` set display limits without changing values or normalization. Record the relevant fluorescence channel and units with the source table; the command plots the supplied histogram without assuming a stain or instrument.

Chromosomes follow their first appearance in the genomic table. The default `--chromosome-layout shared` gives every isolate the same offsets, based on the maximum supplied bin end for each chromosome across all panels. `--chromosome-layout per-panel` instead computes those maxima independently for each isolate, for source profiles originally drawn with isolate-specific chromosome extents. Neither option substitutes reference lengths for absent terminal bins.

The output prefix produces a vector PDF, a PNG and a JSON summary with panel identities, genomic point counts, FACS event totals, display limits and omitted CN points above the display maximum. The supplied density coordinates and all finite genomic observations are used directly. No source workbook or analysis result is overwritten. Use an output prefix outside the input-data directories.

## References

- Ramírez F et al. deepTools2. [10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257).
- Benjamini Y and Speed TP. GC-bias correction. [10.1093/nar/gks001](https://doi.org/10.1093/nar/gks001).
- Li H et al. SAMtools. [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352).
- Quinlan AR and Hall IM. BEDTools. [10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033).
- Cleveland WS. Locally weighted regression. [10.1080/01621459.1979.10481038](https://doi.org/10.1080/01621459.1979.10481038).
- Olshen AB et al. Circular binary segmentation. [10.1093/biostatistics/kxh008](https://doi.org/10.1093/biostatistics/kxh008).
- Abbey DA et al. Candida copy-number and allele-fraction visualization. [10.1186/s13073-014-0100-8](https://doi.org/10.1186/s13073-014-0100-8).
- Skrzypek MS et al. Candida Genome Database and Assembly 22. [10.1093/nar/gkw924](https://doi.org/10.1093/nar/gkw924).

## Included input data

`input_data.zip` contains the unchanged ploidy-validation source workbook: allele-fraction densities, genomic copy number and minor-allele fractions, fluorescence histograms, sample labels and chromosome-level composition. Extract beside the scripts and run commands from this folder:

```sh
unzip input_data.zip
python 09_plot_ploidy_validation.py --source Source_Data_Extended_Data_Fig_3.xlsx --reference-sample fRS1 --facs-min 0 --facs-max 35000 --facs-ymax 35 --cn-max 8 --chromosome-layout per-panel --output-prefix outputs/ploidy_validation
```

These are prepared plotting inputs for the validation figure.
