#!/usr/bin/env bash
# =============================================================================
# Description
#   Gather chromosome-level gVCF shards into one gVCF per isolate.
#
# Inputs
#   BAM_MANIFEST, REFERENCE_FASTA, GVCF_DIR, CHROMOSOMES, SLURM_ARRAY_TASK_ID.
#
# Outputs
#   One gathered gVCF and index per isolate.
# =============================================================================

set -euo pipefail

: "${BAM_MANIFEST:=path/to/bam_manifest.tsv}"
: "${REFERENCE_FASTA:=path/to/reference.fasta}"
: "${GVCF_DIR:=path/to/gvcfs}"
: "${PICARD_BIN:=picard}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/helpers/reference_contigs.sh"
load_chromosomes
validate_array_index

sample=$(awk -F '\t' -v row="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == row {print $1}' "${BAM_MANIFEST}")
if [[ -z "${sample}" ]]; then
  printf 'No sample found for array index %s.\n' "${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

sample_dir="${GVCF_DIR}/${sample}"
inputs=()
for chromosome in "${chromosome_list[@]}"; do
  shard="${sample_dir}/${sample}.${chromosome}.g.vcf.gz"
  [[ -s "${shard}" ]] || { printf 'Missing gVCF shard: %s\n' "${shard}" >&2; exit 3; }
  inputs+=("I=${shard}")
done

"${PICARD_BIN}" GatherVcfs \
  "${inputs[@]}" \
  O="${GVCF_DIR}/${sample}.g.vcf.gz" \
  CREATE_INDEX=true
