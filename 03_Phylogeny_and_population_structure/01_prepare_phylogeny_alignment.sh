#!/usr/bin/env bash
# =============================================================================
# Description
#   Create a one-record-per-position biallelic SNP panel and PHYLIP alignment.
#
# Inputs
#   COHORT_VCF, PHYLOGENY_DIR.
#
# Outputs
#   cohort.phylogeny.biallelic_snps.vcf.gz and cohort.phylogeny.phy.
# =============================================================================

set -euo pipefail

: "${COHORT_VCF:=path/to/cohort.vcf.gz}"
: "${PHYLOGENY_DIR:=path/to/phylogeny}"
: "${BCFTOOLS_BIN:=bcftools}"
: "${PYTHON_BIN:=python3}"

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "${PHYLOGENY_DIR}"
panel="${PHYLOGENY_DIR}/cohort.phylogeny.biallelic_snps.vcf.gz"

"${BCFTOOLS_BIN}" view -a -Ou "${COHORT_VCF}" \
  | "${BCFTOOLS_BIN}" +fill-tags -Ou -- -t AC,AN,AF,NS,F_MISSING,MAF \
  | "${BCFTOOLS_BIN}" view \
      -v snps -m2 -M2 \
      -i 'AC>0 && AC<AN && F_MISSING<=0.10 && (INFO/MQ>=20 || INFO/MQ=".")' \
      -Oz -o "${panel}"
"${BCFTOOLS_BIN}" index --csi --force "${panel}"

"${PYTHON_BIN}" "${script_dir}/helpers/vcf_to_phylip_iupac.py" \
  --vcf "${panel}" \
  --output "${PHYLOGENY_DIR}/cohort.phylogeny.phy"
