#!/usr/bin/env bash
# =============================================================================
# Description
#   Calculate the first 50 principal components from the LD-pruned SNP panel.
#
# Inputs
#   PCA_DIR.
#
# Outputs
#   PLINK eigenvectors and eigenvalues.
# =============================================================================

set -euo pipefail

: "${PCA_DIR:=path/to/pca}"
: "${PLINK_BIN:=plink}"
: "${PCA_COMPONENTS:=50}"

input_prefix="${PCA_DIR}/cohort.ld_pruned_snps"
output_prefix="${PCA_DIR}/cohort.pca"

"${PLINK_BIN}" \
  --bfile "${input_prefix}" \
  --allow-extra-chr \
  --pca "${PCA_COMPONENTS}" \
  --out "${output_prefix}"
