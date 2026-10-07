#!/usr/bin/env bash
# Summarize per-base depth from the indexed, duplicate-filtered isolate BAMs.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
need_command "$MOSDEPTH"
mkdir -p -- "$OUTDIR/depth"
while IFS= read -r sample; do
    bam="$OUTDIR/bam/$sample/$sample.bam"
    need_file "$bam"; need_file "$bam.bai"
    no_output "$OUTDIR/depth/$sample"
    make_stage
    # Default fragment-overlap handling is retained.
    "$MOSDEPTH" --threads "$THREADS" --mapq "$MOSDEPTH_MAPQ" \
        "$STAGE/$sample" "$bam" 2> "$STAGE/mosdepth.log"
    need_file "$STAGE/$sample.mosdepth.summary.txt"
    mv -- "$STAGE" "$OUTDIR/depth/$sample"
done < <(awk -F '\t' 'NR > 1 && !seen[$1]++ {print $1}' "$RUNS_TSV")
