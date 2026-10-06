#!/usr/bin/env bash
# =============================================================================
# Description
#   Prepare an ADMIXTURE-compatible copy of the final LD-pruned PLINK panel and
#   run 10-fold cross-validation for K=2-30 with 10 recorded seeds per K.
#
# Inputs
#   PCA_PREFIX: PLINK prefix containing 1,153 samples and 23,356 SNPs.
#   OUTPUT_DIR: directory for the ADMIXTURE-ready panel and replicate outputs.
#
# Outputs
#   One run.log and the corresponding Q/P files for every K and replicate.
#   Set TASK_ID=1..290 to run one K/replicate combination in an array job.
# =============================================================================

set -euo pipefail

: "${PCA_PREFIX:=path/to/pca/cohort.ld_pruned_snps}"
: "${OUTPUT_DIR:=path/to/admixture_cv}"
: "${PLINK_BIN:=plink}"
: "${ADMIXTURE_BIN:=admixture}"
: "${THREADS:=8}"
: "${CV_FOLDS:=10}"
: "${EXPECTED_SAMPLES:=1153}"
: "${EXPECTED_SNPS:=23356}"
: "${K_MIN:=2}"
: "${K_MAX:=30}"
: "${REPLICATES:=10}"
: "${TASK_ID:=}"

for extension in bed bim fam; do
  if [[ ! -s "${PCA_PREFIX}.${extension}" ]]; then
    printf 'Missing PLINK input: %s\n' "${PCA_PREFIX}.${extension}" >&2
    exit 2
  fi
done

sample_count=$(wc -l < "${PCA_PREFIX}.fam" | tr -d ' ')
snp_count=$(wc -l < "${PCA_PREFIX}.bim" | tr -d ' ')
if [[ "${sample_count}" -ne "${EXPECTED_SAMPLES}" ]]; then
  printf 'Expected %s samples; found %s.\n' "${EXPECTED_SAMPLES}" "${sample_count}" >&2
  exit 2
fi
if [[ "${snp_count}" -ne "${EXPECTED_SNPS}" ]]; then
  printf 'Expected %s SNPs; found %s.\n' "${EXPECTED_SNPS}" "${snp_count}" >&2
  exit 2
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(cd "${OUTPUT_DIR}" && pwd)
admixture_prefix="${OUTPUT_DIR}/admixture_input"
chromosome_updates="${OUTPUT_DIR}/variant_numeric_chromosomes.tsv"

if [[ ! -s "${admixture_prefix}.bed" || ! -s "${admixture_prefix}.bim" || ! -s "${admixture_prefix}.fam" ]]; then
  awk '
    NR == FNR {
      if (!($1 in chromosome_number)) chromosome_number[$1] = ++chromosome_count
      next
    }
    { print $2 "\t" chromosome_number[$1] }
  ' "${PCA_PREFIX}.bim" "${PCA_PREFIX}.bim" > "${chromosome_updates}"

  "${PLINK_BIN}" \
    --bfile "${PCA_PREFIX}" \
    --allow-extra-chr \
    --keep-allele-order \
    --update-chr "${chromosome_updates}" 2 1 \
    --make-bed \
    --out "${admixture_prefix}"
fi

prepared_samples=$(wc -l < "${admixture_prefix}.fam" | tr -d ' ')
prepared_snps=$(wc -l < "${admixture_prefix}.bim" | tr -d ' ')
if [[ "${prepared_samples}" -ne "${EXPECTED_SAMPLES}" || "${prepared_snps}" -ne "${EXPECTED_SNPS}" ]]; then
  printf 'ADMIXTURE input dimensions changed during chromosome recoding.\n' >&2
  exit 2
fi

run_one() {
  local k="$1"
  local replicate="$2"
  local seed=$((200000 + 100 * k + replicate))
  local run_dir="${OUTPUT_DIR}/K${k}/rep${replicate}"
  mkdir -p "${run_dir}"
  (
    cd "${run_dir}"
    "${ADMIXTURE_BIN}" \
      -s "${seed}" \
      --cv="${CV_FOLDS}" \
      -j"${THREADS}" \
      "${admixture_prefix}.bed" \
      "${k}" \
      > run.log 2>&1
  )
}

task_count=$(((K_MAX - K_MIN + 1) * REPLICATES))
if [[ -n "${TASK_ID}" ]]; then
  if [[ "${TASK_ID}" -lt 1 || "${TASK_ID}" -gt "${task_count}" ]]; then
    printf 'TASK_ID must be between 1 and %s.\n' "${task_count}" >&2
    exit 2
  fi
  task_index=$((TASK_ID - 1))
  k=$((K_MIN + task_index / REPLICATES))
  replicate=$((task_index % REPLICATES + 1))
  run_one "${k}" "${replicate}"
else
  for k in $(seq "${K_MIN}" "${K_MAX}"); do
    for replicate in $(seq 1 "${REPLICATES}"); do
      run_one "${k}" "${replicate}"
    done
  done
fi
