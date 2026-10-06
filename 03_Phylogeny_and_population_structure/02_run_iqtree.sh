#!/usr/bin/env bash
# =============================================================================
# Description
#   Infer a maximum-likelihood tree and calculate SH-aLRT and UFBoot2 support.
#
# Inputs
#   ALIGNMENT, PHYLOGENY_DIR.
#
# Outputs
#   IQ-TREE result files under PHYLOGENY_DIR.
# =============================================================================

set -euo pipefail

: "${ALIGNMENT:=path/to/alignment.phy}"
: "${PHYLOGENY_DIR:=path/to/phylogeny}"
: "${IQTREE_BIN:=iqtree2}"
: "${THREADS:=8}"
: "${RANDOM_SEED:=20250812}"
: "${MODEL:=GTR+F+I+R6}"

mkdir -p "${PHYLOGENY_DIR}"
ml_prefix="${PHYLOGENY_DIR}/cohort.ml_shalrt"
ufboot_prefix="${PHYLOGENY_DIR}/cohort.ufboot"

"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -m "${MODEL}" \
  -nt "${THREADS}" \
  -seed "${RANDOM_SEED}" \
  -alrt 1000 \
  -pre "${ml_prefix}"

"${IQTREE_BIN}" \
  -s "${ALIGNMENT}" \
  -m "${MODEL}" \
  -nt "${THREADS}" \
  -seed "${RANDOM_SEED}" \
  -te "${ml_prefix}.treefile" \
  -bb 1000 \
  -bnni \
  -pre "${ufboot_prefix}"
