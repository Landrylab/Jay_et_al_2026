#!/usr/bin/env bash
# =============================================================================
# Description
#   Subset the joint cohort, trim absent alleles, decompose multiallelic records,
#   recalculate cohort-dependent tags, and retain high-confidence biallelic SNPs.
#
# Inputs
#   JOINT_VCF, SAMPLE_KEEP, FILTERED_DIR.
#
# Outputs
#   cohort.filtered.decomposed.biallelic_snps.vcf.gz and index.
# =============================================================================

set -euo pipefail

: "${JOINT_VCF:=path/to/cohort.vcf.gz}"
: "${SAMPLE_KEEP:=path/to/sample_keep.txt}"
: "${FILTERED_DIR:=path/to/filtered_variants}"
: "${BCFTOOLS_BIN:=bcftools}"

mkdir -p "${FILTERED_DIR}"
output_vcf="${FILTERED_DIR}/cohort.filtered.decomposed.biallelic_snps.vcf.gz"

"${BCFTOOLS_BIN}" view -S "${SAMPLE_KEEP}" -a -Ou "${JOINT_VCF}" \
  | "${BCFTOOLS_BIN}" norm -m -both -Ou \
  | "${BCFTOOLS_BIN}" +fill-tags -Ou -- -t AC,AN,AF,NS,F_MISSING,MAF \
  | "${BCFTOOLS_BIN}" view \
      -v snps -m2 -M2 \
      -i 'AC>0 && AC<AN && F_MISSING<=0.10 && (INFO/MQ>=20 || INFO/MQ=".")' \
      -Oz -o "${output_vcf}"

"${BCFTOOLS_BIN}" index --csi --force "${output_vcf}"

record_count=$("${BCFTOOLS_BIN}" index --nrecords "${output_vcf}")
position_count=$("${BCFTOOLS_BIN}" query -f '%CHROM\t%POS\n' "${output_vcf}" | uniq | wc -l | tr -d ' ')
printf 'biallelic_snp_records\t%s\ndistinct_positions\t%s\n' "${record_count}" "${position_count}" \
  > "${FILTERED_DIR}/variant_filtering_summary.tsv"
