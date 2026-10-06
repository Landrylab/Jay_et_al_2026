#!/usr/bin/env bash
# =============================================================================
# Description
#   Call one chromosome of one diploid isolate in GVCF mode with GATK.
#
# Inputs
#   BAM_MANIFEST, REFERENCE_FASTA, GVCF_DIR, CHROMOSOMES, SLURM_ARRAY_TASK_ID.
#
# Outputs
#   One chromosome-level gVCF and index per array task.
# =============================================================================

set -euo pipefail

: "${BAM_MANIFEST:=path/to/bam_manifest.tsv}"
: "${REFERENCE_FASTA:=path/to/reference.fasta}"
: "${GVCF_DIR:=path/to/gvcfs}"
: "${GATK_BIN:=gatk}"
: "${PAIRHMM_THREADS:=2}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/helpers/reference_contigs.sh"
load_chromosomes
validate_array_index
sample_count=$(( $(wc -l < "${BAM_MANIFEST}") - 1 ))
chromosome_count=${#chromosome_list[@]}
task_index=$((SLURM_ARRAY_TASK_ID - 1))
sample_index=$((task_index / chromosome_count + 2))
chromosome_index=$((task_index % chromosome_count))

if (( sample_index > sample_count + 1 )); then
  printf 'Array index exceeds manifest dimensions.\n' >&2
  exit 2
fi

IFS=$'\t' read -r sample bam < <(sed -n "${sample_index}p" "${BAM_MANIFEST}")
chromosome=${chromosome_list[$chromosome_index]}
sample_dir="${GVCF_DIR}/${sample}"
mkdir -p "${sample_dir}"

"${GATK_BIN}" HaplotypeCaller \
  -R "${REFERENCE_FASTA}" \
  -I "${bam}" \
  -L "${chromosome}" \
  -O "${sample_dir}/${sample}.${chromosome}.g.vcf.gz" \
  -ERC GVCF \
  --sample-ploidy 2 \
  --native-pair-hmm-threads "${PAIRHMM_THREADS}"
