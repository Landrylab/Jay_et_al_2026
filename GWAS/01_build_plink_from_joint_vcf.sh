#!/usr/bin/env bash
# =============================================================================
# Description
#   Convert the joint cohort VCF into a full unpruned PLINK bed/bim/fam dataset
#   for downstream GCTA MLMA-LOCO and ML genotype-cache construction.
#
# Inputs
#   joint.keep1153.snps.setid.renamed.vcf.gz
#
# Outputs
#   gcta_mlma_loco/plink_from_joint_vcf/joint_keep1153.setid.vcf.gz
#   gcta_mlma_loco/plink_from_joint_vcf/joint_keep1153.{bed,bim,fam}
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

VCF_IN="${VCF_IN:-${ROOT}/joint.keep1153.snps.setid.renamed.vcf.gz}"
OUTDIR="${OUTDIR:-${ROOT}/gcta_mlma_loco/plink_from_joint_vcf}"
PLINK_PREFIX="${PLINK_PREFIX:-${OUTDIR}/joint_keep1153}"
VCF_SETID="${VCF_SETID:-${OUTDIR}/joint_keep1153.setid.vcf.gz}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
PLINK_BIN="${PLINK_BIN:-plink}"
THREADS="${THREADS:-2}"


# =============================================================================
# Functions
# =============================================================================

require_cmd() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: required command not found: ${cmd}" >&2
    exit 127
  fi
}

require_file() {
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: missing required file: ${path}" >&2
    exit 2
  fi
}

run_cmd() {
  printf '[cmd]'
  printf ' %q' "$@"
  echo
  "$@"
}


# =============================================================================
# Main block
# =============================================================================

mkdir -p "${OUTDIR}"
require_cmd "${BCFTOOLS_BIN}"
require_cmd "${PLINK_BIN}"
require_file "${VCF_IN}"

run_cmd "${BCFTOOLS_BIN}" annotate \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  -Oz \
  -o "${VCF_SETID}" \
  "${VCF_IN}"

run_cmd "${BCFTOOLS_BIN}" index -t "${VCF_SETID}"

run_cmd "${PLINK_BIN}" \
  --vcf "${VCF_SETID}" \
  --double-id \
  --allow-extra-chr \
  --keep-allele-order \
  --make-bed \
  --threads "${THREADS}" \
  --out "${PLINK_PREFIX}"

echo "PLINK prefix: ${PLINK_PREFIX}"

