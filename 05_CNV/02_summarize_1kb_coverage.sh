#!/usr/bin/env bash
# =============================================================================
# Description
#   Summarize read overlap, matched-base depth, and covered-base fraction in
#   non-overlapping 1 kb bins.
# =============================================================================

set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 INPUT_BAM REFERENCE_FAI OUTPUT_TSV" >&2
  exit 2
fi

input_bam=$1
reference_fai=$2
output_tsv=$3
mkdir -p "$(dirname "${output_tsv}")"

for executable in samtools bedtools awk paste cut mktemp; do
  command -v "${executable}" >/dev/null 2>&1 || {
    echo "Required executable not found: ${executable}" >&2
    exit 1
  }
done

work_dir=$(mktemp -d)
trap 'rm -rf "${work_dir}"' EXIT

genome_file="${work_dir}/genome.tsv"
bins_file="${work_dir}/bins.bed"
overlap_file="${work_dir}/read_overlap.tsv"
depth_bed="${work_dir}/depth.bed"
depth_summary="${work_dir}/depth_summary.tsv"
body_file="${work_dir}/coverage_body.tsv"

awk 'BEGIN {OFS="\t"} $1 !~ /chrM/ {print $1, $2}' "${reference_fai}" > "${genome_file}"
bedtools makewindows -g "${genome_file}" -w 1000 > "${bins_file}"

bedtools coverage \
  -a "${bins_file}" \
  -b "${input_bam}" \
  -counts \
  -sorted \
  -g "${genome_file}" > "${overlap_file}"

samtools depth -aa -b "${bins_file}" "${input_bam}" |
  awk 'BEGIN {OFS="\t"} {print $1, $2 - 1, $2, $3, ($3 > 0 ? 1 : 0)}' \
  > "${depth_bed}"

bedtools map \
  -a "${bins_file}" \
  -b "${depth_bed}" \
  -c 4,5 \
  -o median,mean \
  -null 0 \
  -sorted \
  -g "${genome_file}" > "${depth_summary}"

paste "${overlap_file}" <(cut -f4,5 "${depth_summary}") > "${body_file}"

{
  printf 'chr\tstart\tend\treads_overlapping\tmedian_matched\tfraction_covered\n'
  awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}' "${body_file}"
} > "${output_tsv}"
