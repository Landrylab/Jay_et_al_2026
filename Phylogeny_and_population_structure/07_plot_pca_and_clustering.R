#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Generate the six PCA and clustering panels described in Supplementary
#   Figure 1.
#
# Inputs
#   Submission-table directory produced by 06_export_submission_tables.py.
#
# Outputs
#   PDF and PNG files for panels A-F.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

tables_dir <- get_arg("--tables-dir", "path/to/supplementary_tables")
output_dir <- get_arg("--output-dir", "path/to/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(file.path(tables_dir, "pca_scores_and_cluster_assignments.tsv"))
variance <- fread(file.path(tables_dir, "pca_variance_explained.tsv"))
diagnostics <- fread(file.path(tables_dir, "clustering_model_selection_statistics.tsv"))
selection <- fread(file.path(tables_dir, "clustering_selection_summary.tsv"))
final_jaccard <- fread(file.path(tables_dir, "final_clusterwise_jaccard.tsv"))
centroids <- fread(file.path(tables_dir, "cluster_centroids_scaled_pc1_pc20.tsv"))
parameters <- fread(file.path(tables_dir, "clustering_run_parameters.tsv"))

consensus_k <- as.integer(selection[criterion == "consensus_k", value])
silhouette_k <- as.integer(selection[criterion == "maximum_mean_silhouette", value])
row_subsample_k <- as.integer(selection[criterion == "maximum_row_subsample_ARI", value])
jaccard_threshold <- as.numeric(parameters[parameter == "clusterwise_jaccard_stability_threshold", value])
clustering_pc_count <- as.integer(parameters[parameter == "clustering_PCs", value])
if (length(row_subsample_k) != 1L || is.na(row_subsample_k)) {
  row_subsample_k <- diagnostics[which.max(row_subsample_mean_ARI), K]
}
if (!identical(sort(unique(as.integer(scores$cluster))), seq_len(consensus_k))) {
  stop("Cluster assignments do not span the expected final labels.", call. = FALSE)
}
if (nrow(centroids) != consensus_k) {
  stop("Centroid table does not contain one row per final cluster.", call. = FALSE)
}
scores[, cluster_label := factor(cluster_label, levels = paste0("C", sort(unique(cluster))))]

theme_publication <- function() {
  theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot, width = width, height = height, units = "in")
  ggsave(file.path(output_dir, paste0(name, ".png")), plot, width = width, height = height, units = "in", dpi = 300)
}

pca_plot <- ggplot(scores, aes(PC1, PC2, colour = cluster_label)) +
  geom_point(size = 1.1, alpha = 0.82) +
  scale_colour_viridis_d(option = "turbo", end = 0.95) +
  labs(x = "PC1", y = "PC2", colour = "Cluster") +
  theme_publication()
save_figure(pca_plot, "figure1A_pca_pc1_pc2_by_cluster", 7.2, 5.6)

scree_plot <- ggplot(variance, aes(x = factor(PC, levels = PC), y = variance_explained_percent_of_reported_PCs)) +
  geom_col(fill = "#457B9D", width = 0.8) +
  labs(x = "Principal component", y = "Variance among the 50 reported PCs (%)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save_figure(scree_plot, "figure1B_pca_scree_pc1_pc50", 9.0, 4.8)

metric_plot <- function(column, y_label, name, line_colour, marker_k, marker_colour) {
  plot <- ggplot(diagnostics, aes(x = K, y = .data[[column]])) +
    geom_line(linewidth = 0.5, colour = line_colour) +
    geom_point(size = 1.4, colour = line_colour) +
    geom_vline(xintercept = marker_k, linetype = 2, linewidth = 0.6, colour = marker_colour) +
    scale_x_continuous(breaks = seq(min(diagnostics$K), max(diagnostics$K), by = 2)) +
    labs(x = "Number of clusters (K)", y = y_label) +
    theme_publication()
  save_figure(plot, name, 7.2, 4.5)
}

metric_plot(
  "mean_silhouette_width",
  "Mean silhouette width",
  "figure1C_clustering_mean_silhouette",
  "#4D4D4D",
  silhouette_k,
  "#009E73"
)
metric_plot(
  "row_subsample_mean_ARI",
  "Mean row-subsample ARI",
  "figure1D_clustering_row_subsample_ari",
  "#4D4D4D",
  row_subsample_k,
  "#7B2CBF"
)

jaccard_plot <- ggplot(
  final_jaccard,
  aes(x = factor(cluster_label, levels = paste0("C", seq_len(consensus_k))), y = jaccard)
) +
  geom_col(fill = "#BDBDBD", colour = "#000000", linewidth = 0.25, width = 0.8) +
  geom_hline(yintercept = jaccard_threshold, linetype = 2, linewidth = 0.6) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "K=23 cluster", y = "Cluster-wise Jaccard stability") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save_figure(jaccard_plot, "figure1E_final_clusterwise_jaccard", 7.2, 4.8)

centroid_pc_names <- paste0("PC", seq_len(clustering_pc_count))
if (!all(centroid_pc_names %in% names(centroids))) {
  stop("Centroid table is missing one or more clustering PCs.", call. = FALSE)
}
centroid_matrix <- as.matrix(centroids[, ..centroid_pc_names])
rownames(centroid_matrix) <- centroids$cluster_label
centroid_hclust <- hclust(dist(centroid_matrix, method = "euclidean"), method = "ward.D2")

save_dendrogram <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(path, width = 7.2, height = 5.2, useDingbats = FALSE)
  } else {
    png(path, width = 7.2, height = 5.2, units = "in", res = 300, bg = "white")
  }
  old_par <- par(mar = c(5.5, 4.5, 1.0, 1.0))
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  plot(
    centroid_hclust,
    main = "",
    xlab = "K=23 cluster",
    sub = "",
    ylab = "Ward.D2 fusion height",
    hang = -1,
    cex = 0.75
  )
}

save_dendrogram(file.path(output_dir, "figure1F_cluster_centroid_hierarchy_ward_d2.pdf"), "pdf")
save_dendrogram(file.path(output_dir, "figure1F_cluster_centroid_hierarchy_ward_d2.png"), "png")
