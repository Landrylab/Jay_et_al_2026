#!/usr/bin/env bash
# Assess raw and paired trimmed reads; summarize reports across runs.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
need_command "$FASTQC"
need_command "$MULTIQC"
no_output "$OUTDIR/read_qc"
make_stage
while IFS=$'\t' read -r sample run library flowcell lane read1 read2 trim_method; do
    read1=$(read_path "$read1"); read2=$(read_path "$read2")
    need_file "$read1"; need_file "$read2"
    need_file "$OUTDIR/trimmed/$run/R1.fastq.gz"
    need_file "$OUTDIR/trimmed/$run/R2.fastq.gz"
    mkdir -p -- "$STAGE/raw/$run" "$STAGE/trimmed/$run"
    # Stable unique filenames distinguish runs and processing stages in MultiQC.
    mkdir -p -- "$STAGE/input_links/raw/$run" "$STAGE/input_links/trimmed/$run"
    for mate in 1 2; do
        if [[ "$mate" == 1 ]]; then raw="$read1"; else raw="$read2"; fi
        case "$raw" in *.gz) extension="fastq.gz" ;; *) extension="fastq" ;; esac
        ln -s -- "$raw" "$STAGE/input_links/raw/$run/${run}.raw.R${mate}.$extension"
        ln -s -- "$OUTDIR/trimmed/$run/R${mate}.fastq.gz" \
            "$STAGE/input_links/trimmed/$run/${run}.trimmed.R${mate}.fastq.gz"
    done
    "$FASTQC" -t "$THREADS" -o "$STAGE/raw/$run" "$STAGE/input_links/raw/$run/"*
    "$FASTQC" -t "$THREADS" -o "$STAGE/trimmed/$run" "$STAGE/input_links/trimmed/$run/"*
done < <(manifest_rows)
"$MULTIQC" "$STAGE/raw" "$STAGE/trimmed" "$OUTDIR/trimmed" \
    --outdir "$STAGE/summary" --filename multiqc_report.html
rm -rf -- "$STAGE/input_links"
mv -- "$STAGE" "$OUTDIR/read_qc"
