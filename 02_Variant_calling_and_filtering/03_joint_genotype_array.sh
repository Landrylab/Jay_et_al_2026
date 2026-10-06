#!/usr/bin/env bash
# =============================================================================
# Description
#   Import gathered gVCFs into GenomicsDB and jointly genotype one chromosome.
#
# Inputs
#   REFERENCE_FASTA, SAMPLE_MAP, JOINT_DIR, CHROMOSOMES, SLURM_ARRAY_TASK_ID.
#
# Outputs
#   One chromosome-level multisample VCF and index per array task.
# =============================================================================

set -euo pipefail

: "${REFERENCE_FASTA:=path/to/reference.fasta}"
: "${SAMPLE_MAP:=path/to/sample_map.tsv}"
: "${JOINT_DIR:=path/to/joint_genotyping}"
: "${GATK_BIN:=gatk}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/helpers/reference_contigs.sh"
load_chromosomes
validate_array_index
if (( SLURM_ARRAY_TASK_ID > ${#chromosome_list[@]} )); then
  printf 'No chromosome found for array index %s.\n' "${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
chromosome=${chromosome_list[$((SLURM_ARRAY_TASK_ID - 1))]}

workspace="${JOINT_DIR}/genomicsdb.${chromosome}"
mkdir -p "${JOINT_DIR}"

"${GATK_BIN}" GenomicsDBImport \
  --sample-name-map "${SAMPLE_MAP}" \
  --genomicsdb-workspace-path "${workspace}" \
  -L "${chromosome}" \
  --reader-threads "${SLURM_CPUS_PER_TASK:-1}"

"${GATK_BIN}" GenotypeGVCFs \
  -R "${REFERENCE_FASTA}" \
  -V "gendb://${workspace}" \
  -L "${chromosome}" \
  -O "${JOINT_DIR}/cohort.${chromosome}.vcf.gz" \
  --include-non-variant-sites true
