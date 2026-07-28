#!/usr/bin/env bash
# =============================================================================
# Description
#   Apply MAF filtering and LD pruning to a decomposed biallelic SNP VCF.
#
# Inputs
#   FILTERED_VCF, PCA_DIR.
#
# Outputs
#   PLINK files for common and LD-pruned SNP panels plus a count summary.
# =============================================================================

set -euo pipefail

: "${FILTERED_VCF:=path/to/cohort.filtered.decomposed.biallelic_snps.vcf.gz}"
: "${PCA_DIR:=path/to/pca}"
: "${PLINK_BIN:=plink}"
: "${BCFTOOLS_BIN:=bcftools}"
: "${MAF_THRESHOLD:=0.05}"
: "${LD_WINDOW_VARIANTS:=50}"
: "${LD_STEP_VARIANTS:=5}"
: "${LD_R2_THRESHOLD:=0.2}"

mkdir -p "${PCA_DIR}"
common_prefix="${PCA_DIR}/cohort.common_snps"
prune_prefix="${PCA_DIR}/cohort.common_snps.prune"
ld_prefix="${PCA_DIR}/cohort.ld_pruned_snps"
id_stream="${PCA_DIR}/cohort.unique_ids.vcf"

[[ -e "${id_stream}" ]] && { printf 'Temporary stream path already exists: %s\n' "${id_stream}" >&2; exit 2; }
mkfifo "${id_stream}"
producer_pid=""
cleanup_stream() {
  if [[ -n "${producer_pid}" ]]; then
    kill "${producer_pid}" 2>/dev/null || true
  fi
  unlink "${id_stream}" 2>/dev/null || true
}
trap cleanup_stream EXIT

"${BCFTOOLS_BIN}" annotate \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  --output-type v \
  --output "${id_stream}" \
  "${FILTERED_VCF}" &
producer_pid=$!

"${PLINK_BIN}" \
  --vcf "${id_stream}" \
  --double-id \
  --allow-extra-chr \
  --keep-allele-order \
  --maf "${MAF_THRESHOLD}" \
  --make-bed \
  --out "${common_prefix}"

wait "${producer_pid}"
producer_pid=""
unlink "${id_stream}"
trap - EXIT

"${PLINK_BIN}" \
  --bfile "${common_prefix}" \
  --allow-extra-chr \
  --indep-pairwise "${LD_WINDOW_VARIANTS}" "${LD_STEP_VARIANTS}" "${LD_R2_THRESHOLD}" \
  --out "${prune_prefix}"

"${PLINK_BIN}" \
  --bfile "${common_prefix}" \
  --allow-extra-chr \
  --keep-allele-order \
  --extract "${prune_prefix}.prune.in" \
  --make-bed \
  --out "${ld_prefix}"

common_count=$(wc -l < "${common_prefix}.bim" | tr -d ' ')
retained_count=$(wc -l < "${ld_prefix}.bim" | tr -d ' ')
removed_count=$((common_count - retained_count))
sample_count=$(wc -l < "${ld_prefix}.fam" | tr -d ' ')

printf 'metric\tvalue\ncohort_samples\t%s\ncommon_snps_maf_ge_%s\t%s\nld_pruned_snps_retained\t%s\nld_pruned_snps_removed\t%s\n' \
  "${sample_count}" "${MAF_THRESHOLD}" "${common_count}" "${retained_count}" "${removed_count}" \
  > "${PCA_DIR}/pca_panel_summary.tsv"
