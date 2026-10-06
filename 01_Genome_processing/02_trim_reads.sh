#!/usr/bin/env bash
# Trim paired-end reads with the method specified for each sequencing run.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
mkdir -p -- "$OUTDIR/trimmed"
while IFS=$'\t' read -r sample run library flowcell lane read1 read2 trim_method; do
    read1=$(read_path "$read1"); read2=$(read_path "$read2")
    need_file "$read1"; need_file "$read2"
    no_output "$OUTDIR/trimmed/$run"
    make_stage
    if [[ "$trim_method" == trimmomatic ]]; then
        need_command "${TRIMMOMATIC_CMD[0]}"
        need_file "$ADAPTERS"
        "${TRIMMOMATIC_CMD[@]}" PE -threads "$THREADS" -phred33 \
            "$read1" "$read2" "$STAGE/R1.fastq.gz" "$STAGE/R1.unpaired.fastq.gz" \
            "$STAGE/R2.fastq.gz" "$STAGE/R2.unpaired.fastq.gz" \
            "ILLUMINACLIP:$ADAPTERS:2:30:10" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31 \
            2> "$STAGE/trimming.log"
    else
        need_command "$FASTP"
        "$FASTP" -w "$THREADS" -i "$read1" -I "$read2" \
            -o "$STAGE/R1.fastq.gz" -O "$STAGE/R2.fastq.gz" \
            -n 5 -q 15 -u 40 -l 31 --detect_adapter_for_pe -5 \
            --cut_front_window_size 1 --cut_front_mean_quality 3 -r \
            --cut_right_window_size 4 --cut_right_mean_quality 15 -g -x \
            --html "$STAGE/fastp.html" --json "$STAGE/fastp.json" \
            2> "$STAGE/trimming.log"
    fi
    need_file "$STAGE/R1.fastq.gz"; need_file "$STAGE/R2.fastq.gz"
    mv -- "$STAGE" "$OUTDIR/trimmed/$run"
done < <(manifest_rows)
