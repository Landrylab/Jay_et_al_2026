#!/usr/bin/env bash
# Load exact contig names in reference order or the configured explicit order.

load_chromosomes() {
  local contig_text reference_fai reference_dict chromosome
  reference_fai=${REFERENCE_FAI:-${REFERENCE_FASTA}.fai}
  reference_dict=${REFERENCE_DICT:-${REFERENCE_FASTA%.*}.dict}
  if [[ -n "${CHROMOSOMES:-}" ]]; then
    contig_text=$(printf '%s\n' "${CHROMOSOMES}" | awk '{for (i=1;i<=NF;i++) print $i}')
  elif [[ -s "${reference_fai}" ]]; then
    contig_text=$(awk -F '\t' 'NF < 2 || $1 == "" || $2 !~ /^[0-9]+$/ || $2 < 1 {exit 2} {print $1}' "${reference_fai}") || {
      printf 'Invalid reference FASTA index: %s\n' "${reference_fai}" >&2; return 2;
    }
  elif [[ -s "${reference_dict}" ]]; then
    contig_text=$(awk -F '\t' '$1 == "@SQ" {name=""; for (i=2;i<=NF;i++) if ($i ~ /^SN:/) name=substr($i,4); if (name == "") exit 2; print name}' "${reference_dict}") || {
      printf 'Invalid reference sequence dictionary: %s\n' "${reference_dict}" >&2; return 2;
    }
  else
    printf 'Supply CHROMOSOMES, a reference FASTA index, or REFERENCE_DICT.\n' >&2
    return 2
  fi
  if ! printf '%s\n' "${contig_text}" | awk 'NF != 1 || seen[$0]++ || $0 ~ /[[:space:]\/]/ || $0 == "." || $0 == ".." {exit 2} END {if (NR == 0) exit 2}'; then
    printf 'Contig names must be nonempty, unique and valid shard filename components.\n' >&2
    return 2
  fi
  chromosome_list=()
  while IFS= read -r chromosome; do
    chromosome_list+=("${chromosome}")
  done <<<"${contig_text}"
}

validate_array_index() {
  if [[ ! "${SLURM_ARRAY_TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'SLURM_ARRAY_TASK_ID must be a positive integer.\n' >&2
    return 2
  fi
}
