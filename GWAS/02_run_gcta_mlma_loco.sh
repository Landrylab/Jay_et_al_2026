#!/usr/bin/env bash
# =============================================================================
# Description
#   Run GCTA MLMA-LOCO for six antifungals using the full PLINK dataset from
#   step 01. The workflow writes phenotype files, applies the requested MAF
#   filter, constructs an LD-pruned GRM for PCA, runs MIC and binary R/S
#   association models, and renders per-trait plots. GCTA constructs the
#   chromosome-specific LOCO relationship matrices from the association panel.
#
# Inputs
#   gcta_mlma_loco/plink_from_joint_vcf/joint_keep1153.{bed,bim,fam}
#   all_metadata.tsv
#
# Outputs
#   gcta_mlma_loco/gcta_mlma_loco_maf005_six_antifungals/
#   gcta_mlma_loco/gcta_mlma_loco_maf0005_six_antifungals/
#   Per-trait MLMA outputs, plots, logs, GRM files, and PCA covariates.
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
HELPER_DIR="${SCRIPT_DIR}/helpers"

PLINK_PREFIX="${PLINK_PREFIX:-${ROOT}/gcta_mlma_loco/plink_from_joint_vcf/joint_keep1153}"
METADATA_TSV="${METADATA_TSV:-${ROOT}/all_metadata.tsv}"
PLINK_BIN="${PLINK_BIN:-plink}"
GCTA_BIN="${GCTA_BIN:-gcta64}"
THREADS="${THREADS:-2}"
MAF="${MAF:-0.05}"
PRUNE_WINDOW="${PRUNE_WINDOW:-50}"
PRUNE_STEP="${PRUNE_STEP:-5}"
PRUNE_R2="${PRUNE_R2:-0.2}"
N_PCS="${N_PCS:-10}"

if ! command -v "${GCTA_BIN}" >/dev/null 2>&1 && command -v gcta >/dev/null 2>&1; then
  GCTA_BIN="gcta"
fi

case "${MAF}" in
  0.05) RUN_LABEL="gcta_mlma_loco_maf005_six_antifungals"; MAF_LABEL="maf0p05" ;;
  0.005) RUN_LABEL="gcta_mlma_loco_maf0005_six_antifungals"; MAF_LABEL="maf0p005" ;;
  *)
    echo "ERROR: MAF must be 0.05 or 0.005 for the reported analyses." >&2
    exit 2
    ;;
esac

PIPELINE_DIR="${ROOT}/gcta_mlma_loco/${RUN_LABEL}"
INPUT_DIR="${PIPELINE_DIR}/00_inputs"
FILTER_DIR="${PIPELINE_DIR}/01_${MAF_LABEL}_filtered_snps"
GRM_DIR="${PIPELINE_DIR}/02_ld_pruned_grm"
OUTPUT_ROOT="${PIPELINE_DIR}/outputs"
LOG_DIR="${PIPELINE_DIR}/logs"

FILTERED_PREFIX="${FILTER_DIR}/joint_keep1153.${MAF_LABEL}.snps"
PRUNE_PREFIX="${GRM_DIR}/joint_keep1153.${MAF_LABEL}.snps.prune"
PRUNED_PREFIX="${GRM_DIR}/joint_keep1153.${MAF_LABEL}.snps.ldpruned"
GRM_PREFIX="${GRM_DIR}/joint_keep1153.${MAF_LABEL}.snps.ldpruned.grm"
PCA_PREFIX="${GRM_DIR}/joint_keep1153.${MAF_LABEL}.snps.ldpruned.pca"
PCA_COVAR="${GRM_DIR}/covar_pca.pc1_pc${N_PCS}.txt"

DRUGS=(fluconazole voriconazole posaconazole caspofungin micafungin anidulafungin)


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

phenotype_file_for() {
  local drug="$1"
  local mode="$2"
  if [[ "${mode}" == "mic" ]]; then
    echo "${INPUT_DIR}/pheno.mic.${drug}.txt"
  else
    echo "${INPUT_DIR}/pheno.rs_binary.${drug}.txt"
  fi
}

run_gcta_trait() {
  local drug="$1"
  local mode="$2"
  local pheno
  pheno="$(phenotype_file_for "${drug}" "${mode}")"
  local outdir="${OUTPUT_ROOT}/${drug}/${mode}"
  local plotdir="${outdir}/plots"
  local report="${outdir}/report_gcta_mlma_loco.txt"
  local outprefix="${outdir}/gwas"

  mkdir -p "${outdir}" "${plotdir}"
  require_file "${pheno}"

  {
    echo "trait=${drug}"
    echo "mode=${mode}"
    echo "phenotype=${pheno}"
    echo "bfile=${FILTERED_PREFIX}"
    echo "loco_relationship_source=${FILTERED_PREFIX}"
    echo "pca_grm=${GRM_PREFIX}"
    echo "qcovar=${PCA_COVAR}"
    echo "covariates=PC1-PC${N_PCS}"
    echo "genotype_coding=additive_0_1_2"
    echo "maf_filter=${MAF}"
    echo "model=GCTA --mlma-loco"
    if [[ "${mode}" == "rs_binary" ]]; then
      echo "phenotype_coding=1_R_0_S_minus9_missing"
    fi
    echo
  } > "${report}"

  run_cmd "${GCTA_BIN}" \
    --mlma-loco \
    --bfile "${FILTERED_PREFIX}" \
    --pheno "${pheno}" \
    --qcovar "${PCA_COVAR}" \
    --thread-num "${THREADS}" \
    --out "${outprefix}" 2>&1 | tee -a "${report}"

  local mlma=""
  for candidate in "${outprefix}.loco.mlma" "${outprefix}.mlma" "${outdir}/"*.loco.mlma "${outdir}/"*.mlma; do
    if [[ -s "${candidate}" ]]; then
      mlma="${candidate}"
      break
    fi
  done
  if [[ -z "${mlma}" ]]; then
    echo "ERROR: MLMA output not found for ${drug} ${mode}" >&2
    exit 4
  fi

  run_cmd Rscript --vanilla "${HELPER_DIR}/gwas_plot_mlma_results.R" \
    --mlma "${mlma}" \
    --outdir "${plotdir}" \
    --trait "${drug}" \
    --mode "${mode}" 2>&1 | tee -a "${report}"
}


# =============================================================================
# Main block
# =============================================================================

mkdir -p "${INPUT_DIR}" "${FILTER_DIR}" "${GRM_DIR}" "${OUTPUT_ROOT}" "${LOG_DIR}"
require_cmd python3
require_cmd "${PLINK_BIN}"
require_cmd "${GCTA_BIN}"
require_cmd Rscript
require_file "${METADATA_TSV}"

for suffix in bed bim fam; do
  require_file "${PLINK_PREFIX}.${suffix}"
done

run_cmd python3 "${HELPER_DIR}/gwas_make_phenotype_files.py" \
  --fam "${PLINK_PREFIX}.fam" \
  --metadata "${METADATA_TSV}" \
  --outdir "${INPUT_DIR}"

run_cmd "${PLINK_BIN}" \
  --bfile "${PLINK_PREFIX}" \
  --allow-extra-chr \
  --keep-allele-order \
  --maf "${MAF}" \
  --snps-only just-acgt \
  --make-bed \
  --threads "${THREADS}" \
  --out "${FILTERED_PREFIX}"

run_cmd "${PLINK_BIN}" \
  --bfile "${FILTERED_PREFIX}" \
  --allow-extra-chr \
  --indep-pairwise "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" \
  --threads "${THREADS}" \
  --out "${PRUNE_PREFIX}"

run_cmd "${PLINK_BIN}" \
  --bfile "${FILTERED_PREFIX}" \
  --allow-extra-chr \
  --extract "${PRUNE_PREFIX}.prune.in" \
  --make-bed \
  --threads "${THREADS}" \
  --out "${PRUNED_PREFIX}"

run_cmd "${GCTA_BIN}" \
  --bfile "${PRUNED_PREFIX}" \
  --make-grm-bin \
  --thread-num "${THREADS}" \
  --out "${GRM_PREFIX}"

run_cmd "${GCTA_BIN}" \
  --grm "${GRM_PREFIX}" \
  --pca "${N_PCS}" \
  --out "${PCA_PREFIX}"

awk -v OFS="\t" -v n_pcs="${N_PCS}" '
  {
    out = $1 OFS $2
    for (i = 3; i <= 2 + n_pcs; i++) {
      out = out OFS $i
    }
    print out
  }
' "${PCA_PREFIX}.eigenvec" > "${PCA_COVAR}"

for drug in "${DRUGS[@]}"; do
  run_gcta_trait "${drug}" mic
done

for drug in "${DRUGS[@]}"; do
  run_gcta_trait "${drug}" rs_binary
done

run_cmd python3 "${HELPER_DIR}/gwas_summarize_gcta_outputs.py" \
  --run-dir "${PIPELINE_DIR}"

echo "GCTA MLMA-LOCO outputs: ${PIPELINE_DIR}"
