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
: "${CHROMOSOMES:=1 2 3 4 5 6 7 8 9}"
: "${GATK_BIN:=gatk}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

chromosome=$(tr ' ' '\n' <<<"${CHROMOSOMES}" | sed '/^$/d' | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n "${chromosome}" ]] || { printf 'No chromosome found for array index %s.\n' "${SLURM_ARRAY_TASK_ID}" >&2; exit 2; }

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
