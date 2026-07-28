#!/usr/bin/env bash
# =============================================================================
# Description
#   Concatenate chromosome-level joint-genotyping VCFs in reference order.
#
# Inputs
#   JOINT_DIR, CHROMOSOMES.
#
# Outputs
#   cohort.joint.vcf.gz and CSI index.
# =============================================================================

set -euo pipefail

: "${JOINT_DIR:=path/to/joint_genotyping}"
: "${CHROMOSOMES:=1 2 3 4 5 6 7 8 9}"
: "${BCFTOOLS_BIN:=bcftools}"

input_list="${JOINT_DIR}/chromosome_vcfs.list"
: > "${input_list}"
for chromosome in ${CHROMOSOMES}; do
  chromosome_vcf="${JOINT_DIR}/cohort.${chromosome}.vcf.gz"
  [[ -s "${chromosome_vcf}" ]] || { printf 'Missing chromosome VCF: %s\n' "${chromosome_vcf}" >&2; exit 2; }
  printf '%s\n' "${chromosome_vcf}" >> "${input_list}"
done

"${BCFTOOLS_BIN}" concat \
  --file-list "${input_list}" \
  --output-type z \
  --output "${JOINT_DIR}/cohort.joint.vcf.gz"
"${BCFTOOLS_BIN}" index --csi --force "${JOINT_DIR}/cohort.joint.vcf.gz"
