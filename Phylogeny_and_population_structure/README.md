# Phylogeny and population structure

This folder reproduces the phylogenetic and PCA stages described in the Methods. All entry scripts use configurable paths and generic defaults.

Required software: BCFtools/HTSlib 1.21, PLINK v1.9.0-b.8, IQ-TREE v2.1.3, ADMIXTURE v1.3.0, Python 3 with `pysam`, and R 4.2.1 with `cluster`, `clusterCrit`, `data.table`, and `fpc`.

## Numbered Entry Points

1. `01_prepare_phylogeny_alignment.sh`: create an undecomposed, one-record-per-position biallelic SNP panel and convert it to PHYLIP with IUPAC heterozygote codes.
2. `02_run_iqtree.sh`: infer the maximum-likelihood tree and calculate SH-aLRT and UFBoot2 support.
3. `03_prepare_pca_panel.sh`: apply MAF and LD filters to the decomposed biallelic SNP panel.
4. `04_run_pca.sh`: calculate 50 PCs in PLINK.
5. `05_cluster_pcs.R`: verify the 1,153-sample, 23,356-SNP PCA input; evaluate clustering across `k=2..30`; select a consensus `k`; and export assignments, diagnostics, scaled PC1-PC20 centroids, the Ward.D2 hierarchy, and cluster-wise stability.
6. `06_export_submission_tables.py`: reconcile the tree, annotations, PC scores, cluster assignments, model-selection statistics, and centroid hierarchy into journal TSV files.
7. `07_plot_pca_and_clustering.R`: render Supplementary Figure 1 panels A-F from the packaged TSV tables.
8. `08_compare_clusters_to_labels.R`: compare K=2-30 and final K=23 assignments with published MLST clades and isolate metadata; calculate complete-case ARI, purity and arithmetic-mean normalized mutual information; test PC1-PC20 against each label; and export the source data for Supplementary Figure 2.
9. `09_run_admixture_cv.sh`: verify an ADMIXTURE-ready panel of 1,153 isolates and 23,356 SNPs, then run 10-fold cross-validation for K=2-30 with 10 seeds per K. The script can run all 290 fits or one array task selected by `TASK_ID=1..290`.
10. `10_summarize_admixture_cv.R`: verify cohort and SNP dimensions in every completed ADMIXTURE log and export replicate-level and K-level source data for Supplementary Figure 3.

Large VCF and PHYLIP files belong in a permanent sequence-data repository rather than the journal Source Data workbook. The final Newick tree and tabular figure annotations are small enough for the supplementary package.

The undecomposed phylogeny panel contains 957,775 one-record-per-position sites. For PCA, the MAF >= 0.05 filter retained 248,293 SNPs; PLINK `--indep-pairwise 50 5 0.2` removed 224,937 and retained 23,356 SNPs.

Default clustering settings use centered and scaled PC1-PC20, 150 k-means starts, 150 gap-statistic reference replicates, 20 row-subsample replicates per K, 150 candidate Jaccard bootstrap replicates, 200 feature-perturbation replicates, and 200 final cluster-bootstrap replicates. Row and feature perturbations each retain 80% of the relevant dimension; row subsampling is performed without replacement. These values and the PCA input dimensions are written to `clustering_run_parameters.tsv` with every run.

The centroid dendrogram uses Ward.D2 hierarchical clustering on Euclidean distances among the 23 centroids calculated in centered and scaled PC1-PC20 space. `cluster_centroids_scaled_pc1_pc20.tsv` contains the plotted centroids, and `cluster_centroid_hierarchy_ward_d2.tsv` records all 22 merges and Ward fusion heights.

Cluster-label comparisons use available cases for each label. The packaged metadata contain 139 isolates with published MLST clades, 871 with binary fluconazole phenotype, 1,147 with geographic origin and 1,153 with contributing-study labels. Normalized mutual information is mutual information divided by the arithmetic mean of the two label entropies. Benjamini-Hochberg adjustment is applied separately within each external label across PC1-PC20. Named clades in the Supplementary Figure 2 heat-map table exclude samples recorded as `Other`; the complete sample-level table retains all 1,153 isolates.

The gap statistic and mean silhouette width define the two candidates for direct adjudication. If they disagree, the higher K is selected only when it improves the between-to-total sum-of-squares ratio by at least 0.02, mean silhouette width by at least 0.02, or the Calinski-Harabasz index by at least 10%, while its mean clusterwise Jaccard similarity, row-subsample ARI, and feature-perturbation ARI are each no more than 0.02 below the lower-K solution. In the archived publication analysis, the gap and silhouette candidates were K=2 and K=23, respectively, and this rule selected K=23.

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
