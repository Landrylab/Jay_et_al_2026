#!/usr/bin/env bash
# Align each run with a unique read group and a shared sample name per isolate.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
need_command "$BWA_MEM2"
need_command "$SAMTOOLS"
require_reference_indexes
mkdir -p -- "$OUTDIR/aligned"
while IFS=$'\t' read -r sample run library flowcell lane read1 read2 trim_method; do
    need_file "$OUTDIR/trimmed/$run/R1.fastq.gz"
    need_file "$OUTDIR/trimmed/$run/R2.fastq.gz"
    no_output "$OUTDIR/aligned/$run"
    make_stage
    # ID and PU retain flowcell/lane provenance; SM connects runs of an isolate.
    read_group="@RG\\tID:${flowcell}.${lane}.${run}\\tSM:${sample}\\tLB:${library}\\tPL:ILLUMINA\\tPU:${flowcell}.${lane}"
    "$BWA_MEM2" mem -t "$THREADS" -R "$read_group" "$REFERENCE" \
        "$OUTDIR/trimmed/$run/R1.fastq.gz" "$OUTDIR/trimmed/$run/R2.fastq.gz" \
        2> "$STAGE/alignment.log" | \
        "$SAMTOOLS" sort -@ "$THREADS" -m "$SORT_MEMORY" -T "$STAGE/sort" \
            -o "$STAGE/sorted.bam" -
    "$SAMTOOLS" calmd -bAQr "$STAGE/sorted.bam" "$REFERENCE" \
        > "$STAGE/calmd.bam" 2> "$STAGE/calmd.log"
    "$SAMTOOLS" quickcheck "$STAGE/calmd.bam"
    rm -- "$STAGE/sorted.bam"
    mv -- "$STAGE" "$OUTDIR/aligned/$run"
done < <(manifest_rows)
