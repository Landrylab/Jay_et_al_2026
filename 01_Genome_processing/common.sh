#!/usr/bin/env bash
# Shared configuration and validation for paired-end read processing.
set -euo pipefail

die() { printf '%s\n' "ERROR: $*" >&2; exit 1; }
need_file() { [[ -s "$1" ]] || die "Missing or empty input: $1"; }
need_command() { command -v "$1" >/dev/null 2>&1 || die "Command unavailable: $1"; }
no_output() { [[ ! -e "$1" && ! -L "$1" ]] || die "Output already exists: $1"; }
absolute_path() {
    case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$2" "$1" ;; esac
}

load_config() {
    [[ $# -eq 1 ]] || die "Provide one configuration file."
    need_file "$1"
    CONFIG_DIR=$(cd -- "$(dirname -- "$1")" && pwd)
    # The configuration is a trusted shell file maintained by the analyst.
    source "$CONFIG_DIR/$(basename -- "$1")"
    : "${RUNS_TSV:?Set RUNS_TSV}" "${REFERENCE:?Set REFERENCE}" "${OUTDIR:?Set OUTDIR}"
    THREADS=${THREADS:-8}
    SORT_MEMORY=${SORT_MEMORY:-2G}
    MOSDEPTH_MAPQ=${MOSDEPTH_MAPQ:-0}
    BWA_MEM2=${BWA_MEM2:-bwa-mem2}
    SAMTOOLS=${SAMTOOLS:-samtools}
    FASTP=${FASTP:-fastp}
    FASTQC=${FASTQC:-fastqc}
    MULTIQC=${MULTIQC:-multiqc}
    MOSDEPTH=${MOSDEPTH:-mosdepth}
    if ! declare -p PICARD_CMD >/dev/null 2>&1; then PICARD_CMD=(picard); fi
    if ! declare -p TRIMMOMATIC_CMD >/dev/null 2>&1; then TRIMMOMATIC_CMD=(trimmomatic); fi
    [[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer."
    [[ "$MOSDEPTH_MAPQ" =~ ^[0-9]+$ ]] && [[ "$MOSDEPTH_MAPQ" -le 255 ]] || die "MOSDEPTH_MAPQ must be between 0 and 255."
    RUNS_TSV=$(absolute_path "$RUNS_TSV" "$CONFIG_DIR")
    REFERENCE=$(absolute_path "$REFERENCE" "$CONFIG_DIR")
    if [[ -n "${REFERENCE_DICT:-}" ]]; then
        REFERENCE_DICT=$(absolute_path "$REFERENCE_DICT" "$CONFIG_DIR")
    else
        reference_name=$(basename -- "$REFERENCE")
        reference_name=${reference_name%.gz}
        REFERENCE_DICT="$(dirname -- "$REFERENCE")/${reference_name%.*}.dict"
    fi
    OUTDIR=$(absolute_path "$OUTDIR" "$CONFIG_DIR")
    ADAPTERS=$(absolute_path "${ADAPTERS:-path/to/adapters.fa}" "$CONFIG_DIR")
    need_file "$RUNS_TSV"
    RUNS_DIR=$(cd -- "$(dirname -- "$RUNS_TSV")" && pwd)
    awk -F '\t' '
        NR == 1 {
            if ($0 != "sample_id\trun_id\tlibrary_id\tflowcell\tlane\tread1\tread2\ttrim_method") exit 1
            next
        }
        NF != 8 {exit 1}
        {
            for (i = 1; i <= 8; i++) if ($i == "" || $i ~ /\r/) exit 1
            for (i = 1; i <= 5; i++) if ($i !~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/) exit 1
            if (seen[$2]++) exit 1
            if ($8 != "trimmomatic" && $8 != "fastp") exit 1
            rows++
        }
        END {if (rows < 1) exit 1}
    ' "$RUNS_TSV" || die "Invalid run manifest: use the exact header, eight non-empty tab-separated fields, unique run IDs and a valid trim_method."
    mkdir -p -- "$OUTDIR"
}

manifest_rows() { awk 'NR > 1 {print}' "$RUNS_TSV"; }
read_path() { absolute_path "$1" "$RUNS_DIR"; }

make_stage() {
    STAGE=$(mktemp -d "$OUTDIR/.stage.XXXXXX")
    trap 'rm -rf -- "$STAGE"' EXIT
}

require_reference_indexes() {
    need_file "$REFERENCE"
    need_file "$REFERENCE.fai"
    need_file "$REFERENCE.bwt.2bit.64"
    need_file "$REFERENCE.0123"
    need_file "$REFERENCE.amb"
    need_file "$REFERENCE.ann"
    need_file "$REFERENCE.pac"
}
