# FungAMR substitution matching and Figure 2 assessment

`01_match_master_fungamr.py` matches exact drug-specific FungAMR substitutions to the shared standardized SNP/indel catalogue. Its default Figure 2 marker definition uses resolved positive confidence scores **1–4**, following the marker rules documented in Methods. `02_assess_fungamr.py` calculates the selected cohort's classifications, heatmap, gene totals and Fisher tests. Keep both scripts and `fungamr_evidence.py` together.

## Matching from the full catalogue

```bash
python 01_match_master_fungamr.py \
  --master-prefix /path/to/cohort.snps_indels.master \
  --gene-map /path/to/gene_map.tsv \
  --fungamr /path/to/FungAMR.csv \
  --phenotypes /path/to/isolate_phenotypes.tsv \
  --metadata /path/to/isolate_metadata.tsv \
  --confidence-scores 1 2 3 4 \
  --exclude-studies Chew2021 Gong2023 \
  --outdir /path/to/fungamr_source_data
```

Dependencies are NumPy, pandas and SciPy. The master prefix supplies `.gene_annotations.tsv.gz`, `.dosage.npy` and `.samples.tsv`. The reference gene map supplies systematic identifiers and primary gene symbols. The FungAMR database retains its original column names. The phenotype TSV requires unique `sample` identifiers and `rbin_flu`, `rbin_vor`, `rbin_pos`, `rbin_cas`, `rbin_ani` and `rbin_mic`, encoded as R/S or 1/0. The assessment also accepts `R phenotype` and `S phenotype`. Other values are omitted separately for each drug. Metadata requires `sample` and `Study_clean`; studies must agree when also present in the phenotype table. The matching script defaults metadata to the phenotype file when `--metadata` is omitted.

`--confidence-scores` accepts selected positive integers from 1 to 8; the default is `1 2 3 4`. Selecting `1 2` permits a narrower sensitivity analysis with the current Methods-defined marker rules. It does not restore other historical assumptions automatically.

## Marker and cohort rules

1. Keep standalone exact single-amino-acid *Candida albicans* substitution evidence. Harmonize gene aliases before resolving evidence.
2. Resolve all nonzero signed evidence within each canonical gene, exact substitution and drug by the smallest absolute score. Positive resistance evidence wins equal-magnitude ties: +8/−2 resolves to −2, while +2/−2 resolves to +2. Apply the selected positive score set only after this resolution.
3. Apply the marker rules documented in Methods. FKS1 S645P uses its original drug-specific resolved scores; there is no echinocandin-wide override. Multiple-substitution contexts never qualify individual component substitutions. CNV and expression mechanisms are outside this figure's exact-substitution marker definition.
4. Match exact reference gene and protein consequence to allele-specific annotations. A called alternate-allele dosage greater than zero defines a carrier. Missing dosage does not count as a detected substitution.
5. Produce the full 1,153-isolate carrier catalogue before applying cohort exclusions. Exclude Chew2021 and Gong2023 only for the Figure 2b assessment. The remaining 769 isolates include 494 with at least one corresponding phenotype, contributing 1,901 isolate–drug observations.
6. Use **R phenotype** and **S phenotype** in all displayed phenotype labels. Raw `R`/`S` values remain in computational tables. Machine field names containing `resistant` or `susceptible` are retained for compatibility and map to R phenotype and S phenotype, respectively; these are not alternative display labels.

## Assessment and display calculations

The heatmap uses database-mutation-present or database-mutation-absent group totals as denominators for each drug. Each sample appears once per drug regardless of multiple qualifying substitutions. Absence means no detected qualifying drug-specific substitution.

Each gene bar is the percentage of all isolates with R phenotype for that drug carrying at least one qualifying substitution in that gene. Gene-carrier sets overlap. Bar totals and source counts remain unique sample counts. For the stacked Figure 2b display, each carrier contributes a weight of 1/k to each of its k qualifying substitutions within that gene. The five largest weighted substitution contributions are shown separately, with remaining contributions grouped as Others. Ties use protein-position order. These display weights partition the gene total; they are not estimates of causal effects. Stack preparation and rendering are maintained with the figure provenance.

Two-sided Fisher exact tests use `[[present_R, present_S], [absent_R, absent_S]]`. Sample odds ratios use those four counts. Benjamini–Hochberg correction includes drugs with all row and column margins nonzero. Under the current Figure 2 definition, all six drugs are evaluable. Display percentages use decimal half-up rounding to one decimal; exact counts and unrounded statistics are retained.

## Outputs

The matching stage writes:

- `fungamr_signed_evidence_resolution.csv`: all eligible exact gene–substitution–drug evidence groups and their resolved signed scores under the Methods-defined marker rules.
- `fungamr_high_confidence_reference_records.csv`: eligible database evidence records supporting selected scores.
- `master_fungamr_allele_matches.csv`: exact allele matches with canonical identifiers, transcript annotations and protein consequences.
- `master_fungamr_allele_callability.csv`: called, missing and carrier counts.
- `fungamr_high_confidence_carriers.csv`: full-cohort sample/allele/gene/substitution/drug matches.
- `master_fungamr_audit.json`: definitions, input hashes and catalogue totals.

The assessment stage writes `panel_B_cohort_selection.csv`, `panel_B_isolate_classifications.csv`, `antifungal_resistance_explained_confusion_metrics.csv`, `antifungal_resistance_explained_gene_contributions.csv`, `Figure_2_phylogeny_and_FungAMR_panel_B_verified_cells.csv`, `Figure_2_phylogeny_and_FungAMR_panel_B_verified_genes.csv`, `panel_B_fisher_enrichment.csv` and `panel_B_assessment_audit.json`. The legacy filenames and machine fields are preserved for downstream compatibility; user-facing phenotype values in verified heatmap cells are R phenotype/S phenotype.

## Included reproducible assessment inputs

`input_data.zip` contains the current full-cohort `matches.csv`, the original cohort `phenotypes.tsv` with study labels, and a manifest with checksums and rule definitions. Extract it beside the scripts:

```bash
unzip input_data.zip
python 02_assess_fungamr.py \
  --matches matches.csv \
  --phenotypes phenotypes.tsv \
  --metadata phenotypes.tsv \
  --confidence-scores 1 2 3 4 \
  --exclude-studies Chew2021 Gong2023 \
  --outdir results
```

This reproduces the Figure 2b assessment from the supplied signed-resolved carrier catalogue. Matching new genotypes requires the full master and reference inputs described above. Other manuscript analyses may retain their separately documented FungAMR thresholds; the Figure 2 threshold should not be silently substituted into those workflows.

## Validate the portable assessment

```bash
python validate_current_figure2.py --report results/validation.json
```

The validator reads `input_data.zip` beside the script, verifies its member checksums, replays the assessment in a temporary directory, reconciles isolate counts, and independently enumerates Fisher probabilities and the Benjamini–Hochberg adjustment. It preserves the input archive and source files. `--inputs` selects another archive or an extracted directory containing the same manifest and two tables.

To compare the replay with deposited Figure 2b source tables, supply their directory:

```bash
python validate_current_figure2.py \
  --reference-dir path/to/figure2_source_tables \
  --report results/validation.json
```

The report identifies every compared table. Without `--reference-dir`, the check verifies the portable assessment internally and does not claim agreement with an external reference. Validation does not rerun the upstream genotype-to-FungAMR matching stage.
