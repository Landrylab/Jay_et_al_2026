# Paths may contain spaces. Relative paths are resolved from this file.
RUNS_TSV="path/to/runs.tsv"
REFERENCE="reference.fasta"
# Optional override; the default replaces the reference extension with .dict.
# REFERENCE_DICT="path/to/reference.dict"
ADAPTERS="path/to/TruSeq3-PE.fa"
OUTDIR="path/to/genome_processing_results"

THREADS=8
SORT_MEMORY="2G" # Per sorting thread.
MOSDEPTH_MAPQ=0

BWA_MEM2="bwa-mem2"
SAMTOOLS="samtools"
FASTP="fastp"
FASTQC="fastqc"
MULTIQC="multiqc"
MOSDEPTH="mosdepth"
# Arrays also allow java -jar invocation, for example:
# PICARD_CMD=(java -Xmx16g -jar "path/to/picard.jar")
PICARD_CMD=(picard)
TRIMMOMATIC_CMD=(trimmomatic)
