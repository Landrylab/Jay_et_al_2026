#!/usr/bin/env bash
# =============================================================================
# Description
#   Concatenate chromosome-level joint-genotyping VCFs in reference order.
#
# Inputs
#   REFERENCE_FASTA, JOINT_DIR, CHROMOSOMES.
#
# Outputs
#   cohort.joint.vcf.gz and CSI index.
# =============================================================================

set -euo pipefail

: "${JOINT_DIR:=path/to/joint_genotyping}"
: "${REFERENCE_FASTA:=path/to/reference.fasta}"
: "${BCFTOOLS_BIN:=bcftools}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/helpers/reference_contigs.sh"
load_chromosomes

input_list="${JOINT_DIR}/chromosome_vcfs.list"
inputs=()
for chromosome in "${chromosome_list[@]}"; do
  chromosome_vcf="${JOINT_DIR}/cohort.${chromosome}.vcf.gz"
  [[ -s "${chromosome_vcf}" ]] || { printf 'Missing chromosome VCF: %s\n' "${chromosome_vcf}" >&2; exit 2; }
  inputs+=("${chromosome_vcf}")
done
printf '%s\n' "${inputs[@]}" > "${input_list}"

"${BCFTOOLS_BIN}" concat \
  --file-list "${input_list}" \
  --output-type z \
  --output "${JOINT_DIR}/cohort.joint.vcf.gz"
"${BCFTOOLS_BIN}" index --csi --force "${JOINT_DIR}/cohort.joint.vcf.gz"
