# LOH detection and quantification

This workflow counts callable and heterozygous SNP positions in non-overlapping 10 kb windows, calls low-heterozygosity tracts, identifies tracts shared among neighboring samples in the displayed phylogenetic order, resolves recurrent candidate regions within annotated clades, and exports submission tables.

## Inputs

- `path/to/cohort.vcf.gz`: indexed, filtered cohort SNP VCF/BCF.
- `path/to/reference.fasta.fai`: reference FASTA index.
- `path/to/final_samples.tsv`: table with a `sample` column containing the 1,153 final isolates.
- `path/to/phylogenetic_display_order.tsv`: ladderized Main Figure 1 tree-tip display order with a `sample` column.
- `path/to/metadata.tsv`: metadata with `sample` and `Clade` columns.
- `path/to/copy_number_10kb/`: optional per-sample copy-number annotations.

Records sharing a genomic position are counted once, which prevents decomposed multiallelic records from contributing multiple alignment positions. Mitochondrial windows are excluded.

## Run order

```bash
python3 01_count_heterozygous_snps_by_window.py \
  --vcf path/to/cohort.vcf.gz \
  --fai path/to/reference.fasta.fai \
  --samples path/to/final_samples.tsv \
  --out path/to/loh_windows.10kb.tsv.gz

python3 02_call_loh_tracts.py \
  --windows path/to/loh_windows.10kb.tsv.gz \
  --samples path/to/final_samples.tsv \
  --copy-number-dir path/to/copy_number_10kb \
  --out path/to/loh_tracts.10kb.tsv

python3 03_classify_candidate_synapomorphic_loh.py \
  --tracts path/to/loh_tracts.10kb.tsv \
  --tip-order path/to/phylogenetic_display_order.tsv \
  --out path/to/loh_tracts.candidate_synapomorphic.tsv

Rscript 04_resolve_candidate_loh_by_clade.R \
  --metadata path/to/metadata.tsv \
  --tracts path/to/loh_tracts.candidate_synapomorphic.tsv \
  --fai path/to/reference.fasta.fai \
  --out path/to/loh_resolved_clade_regions.tsv

python3 05_export_loh_submission_tables.py \
  --windows path/to/loh_windows.10kb.tsv.gz \
  --tracts path/to/loh_tracts.candidate_synapomorphic.tsv \
  --resolved path/to/loh_resolved_clade_regions.tsv \
  --samples path/to/final_samples.tsv \
  --fai path/to/reference.fasta.fai \
  --out-dir path/to/supplementary_tables \
  --source-data-out path/to/source_data/Figure_1F_LOH_heterozygous_SNP_windows.tsv.gz
```

The copy-number directory is optional and affects only tract annotation. LOH detection does not depend on copy-number availability.

## Recorded thresholds

A low-heterozygosity window contains at least 10 callable SNP positions and no more than 10 heterozygous SNP positions. A tract contains at least three low-heterozygosity windows, allows at most three consecutive intervening windows, and has at least 0.80 low-window purity. A tract is classified as candidate synapomorphic when at least 0.80 of the two neighboring displayed tips on either side carry an overlapping tract on the same chromosome with interval Jaccard similarity of at least 0.20.

For annotated clades, intervals are grouped when they overlap or lie within 10 kb, differ in midpoint by no more than 150 kb, and differ in length by no more than threefold. Consensus starts and ends are the 25th and 75th percentiles, respectively. The reporting flag requires presence in more than 30% of a clade and a resolved length greater than 100 kb.

## References

- Skrzypek MS et al. *Nucleic Acids Research* (2017), Assembly 22 and the Candida Genome Database. DOI: [10.1093/nar/gkw924](https://doi.org/10.1093/nar/gkw924).
- Bonfield JK et al. *GigaScience* (2021), HTSlib. DOI: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007).
- Abbey D et al. *G3* (2011), genome-wide LOH detection in *Candida albicans*. DOI: [10.1534/g3.111.000885](https://doi.org/10.1534/g3.111.000885).
- Ford CB et al. *eLife* (2015), LOH in clinical *Candida albicans* isolates. DOI: [10.7554/eLife.00662](https://doi.org/10.7554/eLife.00662).
- Jaccard P. *Bulletin de la Société Vaudoise des Sciences Naturelles* (1901), Jaccard similarity. DOI: [10.5169/seals-266450](https://doi.org/10.5169/seals-266450).
- Harris CR et al. *Nature* (2020), NumPy. DOI: [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2).
- McKinney W. *Proceedings of the 9th Python in Science Conference* (2010), pandas. DOI: [10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a).
- R Core Team. *R: A Language and Environment for Statistical Computing*. [R Project](https://www.R-project.org/). The base R software citation does not assign a DOI.
