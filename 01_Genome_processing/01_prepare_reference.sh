#!/usr/bin/env bash
# Prepare alignment indexes and the SAM reference sequence dictionary.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/common.sh"
load_config "$@"
need_file "$REFERENCE"
need_command "$BWA_MEM2"
need_command "$SAMTOOLS"

present=0
for suffix in 0123 amb ann bwt.2bit.64 pac; do
    [[ ! -e "$REFERENCE.$suffix" ]] || present=$((present + 1))
done
if [[ "$present" -eq 0 ]]; then
    "$BWA_MEM2" index "$REFERENCE"
elif [[ "$present" -ne 5 ]]; then
    die "Partial BWA-MEM2 index found; use a clean reference copy or complete the index separately."
fi
if [[ ! -e "$REFERENCE.fai" ]]; then "$SAMTOOLS" faidx "$REFERENCE"; fi
require_reference_indexes
if [[ ! -e "$REFERENCE_DICT" && ! -L "$REFERENCE_DICT" ]]; then
    mkdir -p -- "$(dirname -- "$REFERENCE_DICT")"
    "$SAMTOOLS" dict -o "$REFERENCE_DICT" "$REFERENCE"
fi
need_file "$REFERENCE_DICT"
