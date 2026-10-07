# FungAMR substitution matching and cohort assessment

`01_match_master_fungamr.py` identifies drug-matched high-confidence substitutions from the standardized, annotated SNP-and-indel catalogue. It uses the complete allele-specific annotation table, including every annotated transcript, and counts carriers from the same target-alternate-allele dosage matrix used by the association and prediction workflows.

## Inputs and execution

```bash
python 01_match_master_fungamr.py \
  --master-prefix /path/to/cohort.snps_indels.master \
  --gene-map /path/to/gene_map.tsv \
  --fungamr /path/to/FungAMR.csv \
  --phenotypes /path/to/isolate_phenotypes.tsv \
  --metadata /path/to/isolate_metadata.tsv \
  --exclude-studies Chew2021 Gong2023 \
  --outdir /path/to/fungamr_source_data
```

Python dependencies are NumPy, pandas and SciPy. Keep both scripts and `fungamr_evidence.py` together. The master prefix supplies `.gene_annotations.tsv.gz`, `.dosage.npy` and `.samples.tsv`. The reference gene map supplies systematic identifiers and primary gene symbols. The FungAMR CSV or TSV retains its published column names. The phenotype TSV requires a unique `sample` column and `rbin_flu`, `rbin_vor`, `rbin_pos`, `rbin_cas`, `rbin_ani` and `rbin_mic`, encoded as R/S or 1/0. Other phenotype values are excluded separately for each drug. The metadata TSV requires unique `sample` identifiers and exact study labels in `Study_clean`; study assignments must agree if also present in the phenotype table. The matching script can omit `--metadata` when the phenotype table already supplies these assignments.

`01_match_master_fungamr.py` always produces the full matching and carrier catalogue for the master genotype cohort, then calls `02_assess_fungamr.py` for the selected phenotype assessment. Study exclusions never filter reference records, allele matches, callability or the full carrier catalogue. Omit `--exclude-studies` to assess all studies.

To reassess an existing catalogue without repeating matching or changing its files:

```bash
python 02_assess_fungamr.py \
  --matches /path/to/fungamr_source_data/fungamr_score_1_or_2_nonsynonymous_matches.csv \
  --phenotypes /path/to/isolate_phenotypes.tsv \
  --metadata /path/to/isolate_metadata.tsv \
  --exclude-studies Chew2021 Gong2023 \
  --outdir /path/to/fungamr_source_data
```

The standalone assessment uses the phenotype table as its cohort; supply the phenotype table corresponding to the matched genotype cohort. The integrated matching workflow explicitly restricts assessment to its master sample list. An unknown exclusion label, duplicate sample identifier, missing study assignment or conflicting study assignment stops assessment with an error.

## Matching and counting

- Retain exact single-amino-acid *Candida albicans* records and resolve gene aliases before comparing evidence. Within each canonical reference gene, exact substitution and drug, retain the nonzero confidence score with the smallest absolute value; positive resistance evidence wins equal-magnitude ties. For example, +8 and −2 resolve to −2, whereas +2 and −2 resolve to +2. Only resolved scores +1 or +2 define high-confidence resistance markers.
- Match each database substitution to its reference gene and exact protein consequence in a protein-changing allele-specific annotation. Gene aliases are resolved to reference symbols.
- Retain the V488I analysis exclusion. No gene–drug pair is excluded: FKS1–fluconazole evidence is assessed by the same signed-score and exact-single-substitution rules as other pairs. Evidence for a component within a multiple-substitution record remains a separate combination context and does not qualify the component as an exact single-substitution marker.
- Count an isolate as a carrier when its complete genotype contains at least one copy of the matched alternate allele. Missing dosage remains unavailable and does not count as a detected substitution.
- Consolidate duplicate evidence records for the same isolate, allele, gene, substitution and drug, retaining the resolved positive confidence score.
- Exclude specified studies from the assessment cohort before calculating any phenotype counts, heatmap percentages, gene contributions or statistical tests. Each eligible isolate is counted once per drug, regardless of the number of matched mutations or evidence records.
- Calculate each drug's heatmap using mutation-present or mutation-absent isolates as denominators. Mutation absence denotes no detected qualifying substitution, rather than proof of susceptibility.
- Calculate gene bars against all resistant isolates with the corresponding phenotype. Carrier sets for different genes may overlap.
- Test each drug's association between mutation presence and resistance with a two-sided Fisher exact test on `[[present_R, present_S], [absent_R, absent_S]]`. The reported odds ratio is the sample odds ratio returned by SciPy. Apply Benjamini–Hochberg correction across drugs with nonzero row and column margins; an empty margin makes the odds ratio, p value and adjusted p value unavailable.

## Outputs

The matching stage writes `fungamr_signed_evidence_resolution.csv`, recording all eligible exact mutation–drug evidence groups, raw supporting scores/IDs, the resolved signed score and high-confidence resistance status. It also writes eligible database reference records, allele-level matches with canonical IDs, transcript IDs and protein consequences, allele callability and full isolate-level carrier records. `master_fungamr_audit.json` records input hashes and matching definitions.

The assessment stage writes only these files:

| File | Contents |
| --- | --- |
| `panel_B_cohort_selection.csv` | One row per phenotype-cohort isolate: study, exclusion decision and availability of drug phenotypes. `included_in_assessment` marks retained isolates; `included_in_panel_B` additionally requires at least one usable phenotype. |
| `panel_B_isolate_classifications.csv` | Retained isolate–drug observations with R/S phenotype and mutation-presence classification. |
| `antifungal_resistance_explained_confusion_metrics.csv` | Counts, denominators and unrounded percentages for each drug. |
| `antifungal_resistance_explained_gene_contributions.csv` | Resistant carriers per gene and drug, with the selected cohort's resistant denominator. |
| `Figure_2_phylogeny_and_FungAMR_panel_B_verified_cells.csv` | Heatmap counts, denominators and displayed percentages. |
| `Figure_2_phylogeny_and_FungAMR_panel_B_verified_genes.csv` | Gene-bar counts, denominators and displayed percentages. |
| `panel_B_fisher_enrichment.csv` | The four contingency counts, group totals, resistance percentages, odds ratio, two-sided Fisher p, BH-adjusted q, test status and correction-family size. |
| `panel_B_assessment_audit.json` | Exclusions, cohort totals, methods, package versions, input hashes and assessment-output hashes. |

Percentages for figure labels use decimal half-up rounding to one decimal place; an empty denominator is reported as `n/a`. Unavailable Fisher-test results are written as `NA`. This stage produces assessment tables.

## Included input data

`input_data.zip` contains the signed-evidence-resolved carrier matches (1,266 rows; updated 6 October 2026) and a phenotype table with study labels for the Figure 2b assessment. Extract the files beside the scripts and run from this folder:

```bash
unzip input_data.zip
python3 02_assess_fungamr.py \
  --matches matches.csv \
  --phenotypes phenotypes.tsv \
  --metadata phenotypes.tsv \
  --exclude-studies Chew2021 Gong2023 \
  --outdir results
```

Matching a new cohort with script 01 requires the annotated master catalogue and FungAMR database described above.

The pre-update input archive and source outputs are retained under `outputs/fungamr_signed_score_update_20261006/backups/` in the project root. The input archive is an assessment snapshot; matching new variants requires the full database and master catalogue.

On 6 October 2026, the FKS1–fluconazole pair exclusion was removed. Four exact database pairs are now represented in the signed-resolution table, but none resolves to +1 or +2; the high-confidence carrier snapshot and assessment results are unchanged. The immediately preceding files are retained under `outputs/fungamr_fks1_rule_removal_20261006/backups/`.
