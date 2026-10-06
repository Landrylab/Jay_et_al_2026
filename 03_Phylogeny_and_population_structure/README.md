# Phylogeny and population structure

These scripts prepare phylogenetic and PCA inputs, compare population clusters with sample annotations, and render figures from supplied numeric tables. Set paths and expected input dimensions for the dataset being analysed. Plotting saved sources does not require rerunning sequence processing, clustering or ADMIXTURE.

## Dependencies

Analysis: BCFtools/HTSlib, PLINK, IQ-TREE, ADMIXTURE, Python with `pysam`, `pandas` and `biopython`, and R with `cluster`, `clusterCrit`, `data.table`, `fpc` and `ggplot2`. Requirements depend on the entrypoint: the cluster-label export needs only R and `data.table`; the two numeric-source Python plotters need the packages in `../plotting_requirements.txt`. Reading a workbook uses `openpyxl`. Neither numeric-source plotter requires Biopython or a Newick file.

## Entry points

1. `01_prepare_phylogeny_alignment.sh`: prepare the undecomposed biallelic SNP panel and convert genotypes to IUPAC-coded PHYLIP.
2. `02_run_iqtree.sh`: infer the maximum-likelihood tree and branch support.
3. `03_prepare_pca_panel.sh`: apply MAF and LD filters to the decomposed biallelic SNP panel.
4. `04_run_pca.sh`: calculate principal components in PLINK.
5. `05_cluster_pcs.R`: evaluate clustering candidates, export assignments and diagnostics, and record the centroid hierarchy and stability estimates. Configure expected sample/SNP dimensions and clustering parameters for the input panel.
6. `06_export_submission_tables.py`: reconcile tree annotations, PCA and analysis tables. Its analysis-specific input contract is separate from the numeric-source plotting routes below.
7. `07_plot_pca_and_clustering.R`: render PCA and clustering diagnostics from the analysis TSVs.
8. `08_compare_clusters_to_labels.R`: compare assignments with annotations; export cluster-level counts, available-case fractions, sample coordinates and cluster/clade mappings. The full mode also calculates label agreement and PC-label associations.
9. `09_run_admixture_cv.sh`: run replicated ADMIXTURE cross-validation. Configure the input panel, expected dimensions, K range and number of replicates for the analysis.
10. `10_summarize_admixture_cv.R`: parse ADMIXTURE logs into replicate and K-level tables. Its dimension and run-count checks must match the supplied analysis.
11. `11_plot_admixture_cv.py`: render the cross-validation panel corresponding to Extended Data Figure 1g from saved replicate errors.
12. `12_plot_clade_overlap.py`: render all panels of Extended Data Figure 2 from sample coordinates/annotations, display labels, colours, tree segments and ancestry proportions.

## Cluster-label exports

Run the descriptive export without rerunning clustering or association tests:

```sh
Rscript 08_compare_clusters_to_labels.R \
  --scores inputs/pca_scores_and_cluster_assignments.tsv \
  --metadata inputs/sample_metadata.tsv \
  --output-dir outputs/cluster_labels \
  --heatmaps-only
```

The score table requires `sample`, `cluster`, `PC1` and `PC2`. Metadata require `sample`, `Loc`, `Study_clean`, `rbin_flu`, `MLST_clean` and `Clade`; empty annotation fields are permitted. Sample identifiers must be unique and the two tables must contain the same sample set. Identifiers are read as text, preserving leading zeroes. `--expected-samples` and `--final-k` are optional consistency checks; otherwise sizes are inferred from the input. Cluster identifiers need not be consecutive integers.

`figure2_clade_metadata_heatmaps.tsv` contains `label`, `cluster`, `cluster_label`, `extrapolated_clade`, `display_label`, `external_category`, `count`, `cluster_total_complete`, `within_cluster_fraction` and `K`. Each actual cluster remains a distinct aggregation unit. A clade is a display annotation and does not combine clusters sharing that clade. Clusters without a named clade use their cluster identifier; multiple different clade annotations within a cluster also use that identifier. All observed clusters are retained for every metadata variable.

Fractions use the number of samples in that cluster with an available value for the relevant label. A zero denominator yields `NA`, displayed as grey by the plotter; an observed category absent from a cluster with available metadata has fraction zero. If no samples have any value for a metadata variable, a row with an unavailable category is retained for each cluster. `cluster_label_confusion_k<K>.tsv` applies the same rules to all labels, including MLST. The sample export retains the original assignments and coordinates.

Omit `--heatmaps-only` for the full comparisons. `--clustering-pcs` selects the PCs to scale, defaulting to the available leading PCs up to a maximum of 20. `--k-min` and `--k-max` control alternative K-means fits; the default upper bound is limited by the number of distinct PC vectors. The supplied final assignments are used for the final-K comparison. `--seed` and `--nstart` control reruns. Label comparisons use available cases. Normalized mutual information uses the arithmetic mean of the two entropies; Benjamini–Hochberg correction is applied separately within each label across the tested PCs. Undefined comparisons are reported as unavailable.

## Plotting saved figure sources

Run these commands from this folder, using paths to the supplied source workbooks:

```sh
python 11_plot_admixture_cv.py \
  --source inputs/Source_Data_Extended_Data_Fig_1.xlsx \
  --output outputs/admixture_cv.pdf

python 12_plot_clade_overlap.py \
  --source inputs/Source_Data_Extended_Data_Fig_2.xlsx \
  --output outputs/clade_overlap.pdf
```

Use `.pdf`, `.svg` or `.png` for `--output`. Output directories are created as needed. The inputs are read without modification. These entrypoints reproduce the plotted numeric content and annotations; typography can differ across systems. They do not select K, infer clades, refit a tree, or calculate new ancestry estimates.

The CV plot reads named Excel Table `ED1Data6`: `K`, `Replicate`, `Cross-validation error`. It also accepts a standalone TSV/CSV with those headers, or the analysis headers `K`, `replicate`, `cv_error`. Each K/replicate pair must occur once. Boxes show the interquartile range and median; whiskers extend to observed values within 1.5 times the interquartile range. All runs are drawn as points, including outliers, and the orange line connects the lowest error per K. An optional `--highlight-k` integer and `--highlight-label` string add an author-selected annotation; no knee-selection procedure is performed. No K or input cohort size is assumed by the plotter.

The clade plot uses these named Excel Tables:

| Table | Required columns and meaning |
| --- | --- |
| `ED2Data1` | `Sample`, `Tree tip order`, `Cluster`, `PC1` and `PC2` (optional variance percentages in their headers), `Published MLST clade`, `Extended clade`, `Fluconazole status`, `Location`, `Study`. One row per sample. |
| `ED2Data2` | `Cluster`, `Display label`. One row per observed cluster; row order defines the metadata heatmap order. |
| `ED2Data3` | `Track`, `Category`, `Colour (RGB)`, `Marker`. Track names are `PCA cluster`, `MLST / extended clade`, `ADMIXTURE component`, `Fluconazole status`, `Location`, `Study`. Category row order defines the displayed category/component order. Markers may be empty. |
| `ED2Data4` | `Start x`, `Start y`, `End x`, `End y`. Recorded tree line segments; tip centres must correspond to tree order minus one half. These coordinates encode the displayed tree topology and do not imply substitution distances. |
| `ED2Data5` | `Sample`, `Component`, `Ancestry proportion`. One row per sample and ancestry component. Sample IDs join proportions to the tree order. |

Blank metadata cells represent missing information; nonempty categories, including `Other`, remain labels. Palette entries must cover the observed categories. PCA cluster colours match the exact cluster identifier where present, otherwise its `C`-prefixed form. Sample tree orders must be consecutive and start at one. Ancestry vectors must be complete, bounded between zero and one, and sum to one within rounding tolerance; only drawn bar widths are normalized to unit sum. The number of clusters, components, samples and annotation counts are derived from the supplied data. PCA coordinates and their axis descriptions are read directly from the table.

For standalone tabular sources, export the five tables to `ED2Data1.tsv` through `ED2Data5.tsv` and pass their directory as `--source`. Alternatively supply `--samples`, `--cluster-labels`, `--palette`, `--tree-segments` and `--ancestry` paths individually. Each accepts TSV, CSV or XLSX. `--samples-table`, `--cluster-labels-table`, `--palette-table`, `--tree-segments-table` and `--ancestry-table` override the corresponding named table. The shared `../figure_sources.py` loader also accepts compressed TSV/CSV files. Preserve identifier strings and numeric precision when exporting tables.

The MLST heatmap and scatter show available published labels. The phylogeny, published-label marks, PCA-cluster strip, ancestry bars and extended-clade strip share the supplied sample order. Each metadata panel contains a cluster heatmap and a PCA scatter; missing metadata are omitted only from the corresponding scatter and denominator. Grey heatmap cells denote unavailable fractions.

## Clustering method references

- Hartigan, J. A. & Wong, M. A. A k-means clustering algorithm. *Applied Statistics* **28**, 100-108 (1979). https://doi.org/10.2307/2346830
- Tibshirani, R., Walther, G. & Hastie, T. Estimating the number of clusters in a data set via the gap statistic. *J. R. Stat. Soc. B* **63**, 411-423 (2001). https://doi.org/10.1111/1467-9868.00293
- Rousseeuw, P. J. Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. *J. Comput. Appl. Math.* **20**, 53-65 (1987). https://doi.org/10.1016/0377-0427(87)90125-7
- Hubert, L. & Arabie, P. Comparing partitions. *J. Classification* **2**, 193-218 (1985). https://doi.org/10.1007/BF01908075
- Caliński, T. & Harabasz, J. A dendrite method for cluster analysis. *Commun. Stat.* **3**, 1-27 (1974). https://doi.org/10.1080/03610927408827101
- Davies, D. L. & Bouldin, D. W. A cluster separation measure. *IEEE Trans. Pattern Anal. Mach. Intell.* **1**, 224-227 (1979). https://doi.org/10.1109/TPAMI.1979.4766909
- Hennig, C. Cluster-wise assessment of cluster stability. *Comput. Stat. Data Anal.* **52**, 258-271 (2007). https://doi.org/10.1016/j.csda.2006.11.025
- Vinh, N. X., Epps, J. & Bailey, J. Information theoretic measures for clusterings comparison: variants, properties, normalization and correction for chance. *J. Mach. Learn. Res.* **11**, 2837-2854 (2010). https://jmlr.org/papers/v11/vinh10a.html
- Benjamini, Y. & Hochberg, Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* **57**, 289-300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Alexander, D. H., Novembre, J. & Lange, K. Fast model-based estimation of ancestry in unrelated individuals. *Genome Res.* **19**, 1655-1664 (2009). https://doi.org/10.1101/gr.094052.109
- Alexander, D. H. & Lange, K. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. *BMC Bioinformatics* **12**, 246 (2011). https://doi.org/10.1186/1471-2105-12-246

## Included input data

`input_data.zip` contains the numeric inputs for the PCA/clustering panels, CV replicates, and the clade-overlap source workbook. Extract beside the scripts and run commands from this folder:

```sh
unzip input_data.zip
Rscript 07_plot_pca_and_clustering.R --tables-dir . --output-dir outputs/pca
python 11_plot_admixture_cv.py --source admixture_cv.tsv --output outputs/admixture_cv.pdf
python 12_plot_clade_overlap.py --source Source_Data_Extended_Data_Fig_2.xlsx --output outputs/clade_overlap.pdf
```

The PCA command writes separate panels from the saved analysis.
