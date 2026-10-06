#!/usr/bin/env bash
# Remove duplicate reads after combining the runs belonging to each isolate.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
need_command "$SAMTOOLS"
need_command "${PICARD_CMD[0]}"
mkdir -p -- "$OUTDIR/bam"
while IFS= read -r sample; do
    no_output "$OUTDIR/bam/$sample"
    make_stage
    inputs=()
    while IFS=$'\t' read -r this_sample run library flowcell lane read1 read2 trim_method; do
        [[ "$this_sample" == "$sample" ]] || continue
        need_file "$OUTDIR/aligned/$run/calmd.bam"
        inputs+=("$OUTDIR/aligned/$run/calmd.bam")
    done < <(manifest_rows)
    if [[ ${#inputs[@]} -eq 1 ]]; then
        merged="${inputs[0]}"
    else
        merged="$STAGE/merged.bam"
        "$SAMTOOLS" merge -@ "$THREADS" -o "$merged" "${inputs[@]}"
    fi
    "${PICARD_CMD[@]}" MarkDuplicates "INPUT=$merged" "OUTPUT=$STAGE/$sample.bam" \
        "METRICS_FILE=$STAGE/duplicate_metrics.txt" "TMP_DIR=$STAGE" \
        REMOVE_DUPLICATES=true CREATE_INDEX=false ASSUME_SORT_ORDER=coordinate \
        2> "$STAGE/markduplicates.log"
    "$SAMTOOLS" quickcheck "$STAGE/$sample.bam"
    "$SAMTOOLS" index -@ "$THREADS" "$STAGE/$sample.bam"
    if [[ ${#inputs[@]} -gt 1 ]]; then rm -- "$merged"; fi
    mv -- "$STAGE" "$OUTDIR/bam/$sample"
done < <(awk -F '\t' 'NR > 1 && !seen[$1]++ {print $1}' "$RUNS_TSV")
