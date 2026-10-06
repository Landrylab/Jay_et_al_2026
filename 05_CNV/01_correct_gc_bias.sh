#!/usr/bin/env bash
# =============================================================================
# Description
#   Estimate and correct fragment-level GC bias for one coordinate-sorted BAM.
# =============================================================================

set -euo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 INPUT_BAM REFERENCE_2BIT EFFECTIVE_GENOME_SIZE OUTPUT_BAM" >&2
  exit 2
fi

input_bam=$1
reference_2bit=$2
effective_genome_size=$3
output_bam=$4

for executable in computeGCBias correctGCBias samtools; do
  command -v "${executable}" >/dev/null 2>&1 || {
    echo "Required executable not found: ${executable}" >&2
    exit 1
  }
done

output_dir=$(dirname "${output_bam}")
output_name=$(basename "${output_bam}" .bam)
mkdir -p "${output_dir}"

frequency_table="${output_dir}/${output_name}.gc_bias_frequencies.tsv"
diagnostic_plot="${output_dir}/${output_name}.gc_bias.png"

computeGCBias \
  --bamfile "${input_bam}" \
  --genome "${reference_2bit}" \
  --effectiveGenomeSize "${effective_genome_size}" \
  --GCbiasFrequenciesFile "${frequency_table}" \
  --biasPlot "${diagnostic_plot}"

correctGCBias \
  --bamfile "${input_bam}" \
  --genome "${reference_2bit}" \
  --effectiveGenomeSize "${effective_genome_size}" \
  --GCbiasFrequenciesFile "${frequency_table}" \
  --correctedFile "${output_bam}"

samtools index "${output_bam}"
