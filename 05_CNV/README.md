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

`09_plot_ploidy_validation.py` draws the genomic examples, a flow-cytometry heatmap of every supplied recording, and an optional chromosome-specific aneuploidy composition panel. The current four genomic examples occupy panels a–d, the heatmap is panel e, and composition is panel f. Install the dependencies in `../plotting_requirements.txt`.

```bash
python3 09_plot_ploidy_validation.py \
  --source path/to/figure_inputs \
  --output-prefix path/to/plots/ploidy_validation
```

`--source` accepts a directory of the named CSV/TSV tables below, a workbook containing those named tables, or the current submission workbook. The current workbook adapter reads the genomic named tables and the three headed sections on `Flow cytometry and composition`; it does not require named tables on that sheet. Individual CSV/TSV inputs may be supplied with `--samples`, `--densities`, `--genome`, `--flow-wells` (alias `--wells`), `--flow-bins` and `--composition`.

| Table | Required columns | Row definition |
| --- | --- | --- |
| `ED3Samples` | `Panel`, `Isolate`, `Assigned ploidy (N)` | One genomic example per panel. Optional `Display label` and `Flow well` columns override labels and exact cohort-sample matching. An empty `Flow well` means no flow recording. |
| `ED3Densities` | `Panel`, `Allele fraction`, `REF density`, `ALT density` | One supplied density coordinate. Original smoothing and normalization are retained. |
| `ED3GenomeBins` | `Panel`, `Chromosome`, `Bin start (bp)`, `Bin end (bp)`, `Absolute copy number`, `Median minor AF` | One nonoverlapping genomic bin. Starts are zero-based and ends exclusive. CN is already on the assigned absolute-ploidy scale. Minor AF may be missing. |
| `ED3FlowWells` | `sample`, `well`, `recorded_ploidy`, `group`, `cohort_sample`, `control_strain`, `events` | One recording. `group` is `Cohort isolates` or `Reference controls`. Cohort rows have unique `cohort_sample` values; reference rows have `control_strain` instead. `recorded_ploidy` uses labels such as `2N`. Optional `display_label` controls the complete isolate-name display. |
| `ED3FlowBins` | `well`, `bin_start`, `bin_end`, `events`, `percent_of_all_events` | One fluorescence bin per recording. All recordings must have matching, contiguous bin edges. Counts sum to the well's event total; percentages use all recorded events. |
| `ED3Composition` | `Chromosome`, `Aneuploidy type`, `Copy number`, `Isolates`, `Chromosome denominator`, `Aneuploid isolates (%)` | One chromosome/type/CN category. Types are `whole chromosome` and `segmental`. Category counts sum to each chromosome's denominator. |

Heatmap rows are ordered by increasing recorded ploidy, with cohort isolates followed by controls of the same ploidy. Supplied order within those groups is retained. Genomic examples map to the heatmap by exact `Isolate`/`cohort_sample` equality, unless `Flow well` is supplied. Additional tested samples are retained even when they are not genomic examples. The plate-recorded heatmap ploidy is distinct from the genomic assignment; the plotter does not infer a new assignment from fluorescence.

Default display limits are 0–35,000 fluorescence units, 0–45% per bin on a log(1 + percentage) colour scale, and 0–8 absolute copies. `--flow-min`, `--flow-max`, `--colour-max` and `--cn-max` change display limits without changing normalization or source values. The plotter rejects a colour maximum that would clip an observed bin. The supplied full 0–100,000 range remains in the denominator when the fluorescence display is zoomed. No gating is added.

Chromosomes follow their source order. `--chromosome-layout per-panel` (default) uses each example's maximum supplied chromosome-bin ends; `shared` uses the maximum across examples. Neither mode substitutes reference lengths for absent terminal bins. Use `--no-composition` to omit the final panel.

Outputs are a vector PDF, a 600-dpi PNG (adjustable with `--dpi`) and a JSON summary of event totals, well order, group counts, genomic point counts, arrow positions and display settings. The four-example layout is 215 × 309 mm. `--font-family` accepts an installed font; the portable default is Matplotlib's DejaVu Sans. The figure contains no column titles. Source workbooks and analysis results are read-only. Keep output prefixes outside the input-data directories.

## Confirmed cohort ploidy assignments

`10_apply_reviewed_ploidy.py` applies explicit reviewed assignments to a cohort baseline table and recomputes ploidy totals. The supplied review records BL042 as 3N. The packaged cohort has 1,153 isolates: 1,127 2N, eight 3N, eleven 4N and seven 6N.

```bash
python3 10_apply_reviewed_ploidy.py \
  --ploidy cohort_baseline_ploidy.tsv \
  --review ploidy_reviews.tsv \
  --output-dir results/reviewed_ploidy
```

This updates baseline metadata and produces reviewed-baseline, review and count tables. It does not rescale or relabel existing integer copy-number tracks. When an accepted baseline changes, rerun `04_segment_and_assign_copy_number.R` with `--baseline-ploidy`, then collate with `08_collate_reviewed_calls.py`. Submission export rejects a disagreement between an accepted baseline and an existing track's baseline. The genomic BL042 example already uses 3N; diploid-only chromosome-composition counts are unaffected by this correction.

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

`input_data.zip` contains the current submission source workbook, six plain CSV plotting tables, the written legend, the confirmed cohort baseline/review/count tables, and manifests with source hashes. The CSV plotting route reproduces all displayed row labels. The official workbook can also be read directly; its control strain names are retained when it has no optional `display_label` field.

```sh
unzip input_data.zip
python3 09_plot_ploidy_validation.py --source . --output-prefix outputs/ploidy_validation
python3 10_apply_reviewed_ploidy.py --ploidy cohort_baseline_ploidy.tsv --review ploidy_reviews.tsv --output-dir outputs/reviewed_ploidy
```

To use the workbook's numeric sections with the supplied species-qualified control labels:

```sh
python3 09_plot_ploidy_validation.py --source Source_Data_Extended_Data_Fig_3.xlsx --flow-wells ED3FlowWells.csv --output-prefix outputs/ploidy_validation
```

The flow-cytometry inputs contain all 23 tested members of the study cohort and four selected control recordings: fRS1 and fRS757 (2N), A48P1 (3N) and F18_P35 (4N). No 1N control is included. There are 27 recordings, 5,400 fluorescence-bin records and 116,351 recorded events. The bins are 500 units wide and cover 0–100,000 units. The acquisition channel, control species and experimental details are documented in `legend.md` and the source workbook. Later altered-concentration acquisitions are not pooled with these inputs.
