#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Evaluate unsupervised k-means clustering on centered and scaled PC1-PC20,
#   select a consensus number of clusters, and export figure-ready diagnostics.
#
# Inputs
#   PLINK eigenvector, eigenvalue, and BIM files supplied through command-line
#   options. The defaults assert the final 1,153-isolate, 23,356-SNP panel.
#
# Outputs
#   PC scores, variance summaries, model-selection diagnostics, stability
#   summaries, and final cluster assignments.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(clusterCrit)
  library(fpc)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

eigenvec_path <- get_arg("--eigenvec", "path/to/cohort.pca.eigenvec")
eigenval_path <- get_arg("--eigenval", "path/to/cohort.pca.eigenval")
bim_path <- get_arg("--bim", "path/to/cohort.ld_pruned_snps.bim")
output_dir <- get_arg("--output-dir", "path/to/clustering")
seed <- as.integer(get_arg("--seed", "42"))
pc_count <- as.integer(get_arg("--clustering-pcs", "20"))
expected_samples <- as.integer(get_arg("--expected-samples", "1153"))
expected_input_snps <- as.integer(get_arg("--expected-input-snps", "23356"))
expected_reported_pcs <- as.integer(get_arg("--expected-reported-pcs", "50"))
expected_consensus_k <- as.integer(get_arg("--expected-consensus-k", "23"))
k_min <- as.integer(get_arg("--k-min", "2"))
k_max <- as.integer(get_arg("--k-max", "30"))
gap_replicates <- as.integer(get_arg("--gap-replicates", "150"))
row_replicates <- as.integer(get_arg("--row-replicates", "20"))
jaccard_replicates <- as.integer(get_arg("--jaccard-replicates", "150"))
feature_replicates <- as.integer(get_arg("--feature-replicates", "200"))
final_jaccard_replicates <- as.integer(get_arg("--final-jaccard-replicates", "200"))
nstart <- as.integer(get_arg("--nstart", "150"))
row_subsample_fraction <- as.numeric(get_arg("--row-subsample-fraction", "0.8"))
feature_subsample_fraction <- as.numeric(get_arg("--feature-subsample-fraction", "0.8"))
jaccard_stability_threshold <- as.numeric(get_arg("--jaccard-stability-threshold", "0.75"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(eigenvec_path, header = FALSE)
eigenvalues <- scan(eigenval_path, quiet = TRUE)
input_snp_count <- length(readLines(bim_path, warn = FALSE))
reported_pc_count <- ncol(scores) - 2L
if (reported_pc_count < 2L) stop("At least two PCs are required.", call. = FALSE)
if (pc_count > reported_pc_count) stop("Requested clustering PC count exceeds available PCs.", call. = FALSE)
if (nrow(scores) != expected_samples) {
  stop(sprintf("Expected %d PCA samples; found %d.", expected_samples, nrow(scores)), call. = FALSE)
}
if (input_snp_count != expected_input_snps) {
  stop(sprintf("Expected %d LD-pruned SNPs; found %d.", expected_input_snps, input_snp_count), call. = FALSE)
}
if (reported_pc_count != expected_reported_pcs) {
  stop(sprintf("Expected %d reported PCs; found %d.", expected_reported_pcs, reported_pc_count), call. = FALSE)
}
if (row_subsample_fraction <= 0 || row_subsample_fraction >= 1) {
  stop("Row-subsample fraction must be between zero and one.", call. = FALSE)
}
if (feature_subsample_fraction <= 0 || feature_subsample_fraction >= 1) {
  stop("Feature-subsample fraction must be between zero and one.", call. = FALSE)
}

setnames(scores, c("FID", "sample", paste0("PC", seq_len(reported_pc_count))))
if (anyDuplicated(scores$sample)) stop("Duplicate sample identifiers found in eigenvector file.", call. = FALSE)
if (length(eigenvalues) < reported_pc_count) stop("Fewer eigenvalues than eigenvectors.", call. = FALSE)

pc_names <- paste0("PC", seq_len(reported_pc_count))
pc_scores <- scores[, c("sample", pc_names), with = FALSE]
fwrite(pc_scores, file.path(output_dir, "pca_scores_pc1_pc50.tsv"), sep = "\t")

variance_fraction <- eigenvalues[seq_len(reported_pc_count)] / sum(eigenvalues[seq_len(reported_pc_count)])
variance_table <- data.table(
  PC = pc_names,
  eigenvalue = eigenvalues[seq_len(reported_pc_count)],
  variance_explained_fraction_of_reported_PCs = variance_fraction,
  variance_explained_percent_of_reported_PCs = 100 * variance_fraction,
  cumulative_variance_percent_of_reported_PCs = 100 * cumsum(variance_fraction)
)
fwrite(variance_table, file.path(output_dir, "pca_variance_explained.tsv"), sep = "\t")

x_raw <- as.matrix(scores[, paste0("PC", seq_len(pc_count)), with = FALSE])
rownames(x_raw) <- scores$sample
x <- scale(x_raw, center = TRUE, scale = TRUE)
if (any(!is.finite(x))) stop("Non-finite values encountered after PC scaling.", call. = FALSE)

adjusted_rand_index <- function(partition_a, partition_b) {
  contingency <- table(partition_a, partition_b)
  choose_two <- function(value) value * (value - 1) / 2
  observed <- sum(choose_two(contingency))
  row_pairs <- sum(choose_two(rowSums(contingency)))
  column_pairs <- sum(choose_two(colSums(contingency)))
  total_pairs <- choose_two(sum(contingency))
  expected <- row_pairs * column_pairs / total_pairs
  maximum <- (row_pairs + column_pairs) / 2
  denominator <- maximum - expected
  if (denominator == 0) return(if (observed == maximum) 1 else 0)
  (observed - expected) / denominator
}

set.seed(seed)
k_values <- k_min:k_max
distance_matrix <- dist(x)
kmeans_fits <- setNames(vector("list", length(k_values)), as.character(k_values))
silhouette_values <- setNames(numeric(length(k_values)), as.character(k_values))
ch_values <- setNames(numeric(length(k_values)), as.character(k_values))
db_values <- setNames(numeric(length(k_values)), as.character(k_values))
within_ss <- setNames(numeric(length(k_values)), as.character(k_values))
between_ratio <- setNames(numeric(length(k_values)), as.character(k_values))

for (k in k_values) {
  fit <- kmeans(x, centers = k, nstart = nstart, iter.max = 300)
  kmeans_fits[[as.character(k)]] <- fit
  silhouette_values[[as.character(k)]] <- mean(silhouette(fit$cluster, distance_matrix)[, "sil_width"], na.rm = TRUE)
  internal <- intCriteria(as.matrix(x), as.integer(fit$cluster), c("Calinski_Harabasz", "Davies_Bouldin"))
  ch_values[[as.character(k)]] <- internal$calinski_harabasz
  db_values[[as.character(k)]] <- internal$davies_bouldin
  within_ss[[as.character(k)]] <- sum(fit$withinss)
  between_ratio[[as.character(k)]] <- fit$betweenss / fit$totss
}

gap_fun <- function(data, k) kmeans(data, centers = k, nstart = nstart, iter.max = 300)
set.seed(seed)
gap_result <- clusGap(x, FUN = gap_fun, K.max = k_max, B = gap_replicates)
gap_values <- gap_result$Tab[, "gap"]
gap_se <- gap_result$Tab[, "SE.sim"]
k_gap <- max(k_min, min(k_max, as.integer(maxSE(gap_values, gap_se, method = "firstSEmax"))))
k_silhouette <- as.integer(names(which.max(silhouette_values)))

row_subsample_ari <- setNames(numeric(length(k_values)), as.character(k_values))
for (k in k_values) {
  base_clusters <- kmeans_fits[[as.character(k)]]$cluster
  set.seed(seed + k)
  ari <- replicate(row_replicates, {
    rows <- sample.int(nrow(x), size = floor(row_subsample_fraction * nrow(x)), replace = FALSE)
    fit <- kmeans(x[rows, , drop = FALSE], centers = k, nstart = nstart, iter.max = 300)
    adjusted_rand_index(base_clusters[rows], fit$cluster)
  })
  row_subsample_ari[[as.character(k)]] <- mean(ari, na.rm = TRUE)
}
k_row_stability <- as.integer(names(which.max(row_subsample_ari)))

diagnostics <- data.table(
  K = k_values,
  mean_silhouette_width = unname(silhouette_values[as.character(k_values)]),
  gap_statistic = gap_values[k_values],
  gap_standard_error = gap_se[k_values],
  row_subsample_mean_ARI = unname(row_subsample_ari[as.character(k_values)]),
  calinski_harabasz = unname(ch_values[as.character(k_values)]),
  davies_bouldin = unname(db_values[as.character(k_values)]),
  total_within_cluster_ss = unname(within_ss[as.character(k_values)]),
  between_to_total_ss_ratio = unname(between_ratio[as.character(k_values)])
)
fwrite(diagnostics, file.path(output_dir, "clustering_model_selection_statistics.tsv"), sep = "\t")

match_cluster_labels <- function(partition, reference) {
  overlap <- table(partition, reference)
  candidates <- as.data.table(as.table(overlap))
  setnames(candidates, c("partition_cluster", "reference_cluster", "overlap"))
  setorder(candidates, -overlap, partition_cluster, reference_cluster)
  mapping <- integer(0)
  used_reference <- character(0)
  for (row_index in seq_len(nrow(candidates))) {
    row <- candidates[row_index]
    source <- as.character(row$partition_cluster)
    target <- as.character(row$reference_cluster)
    if (!(source %in% names(mapping)) && !(target %in% used_reference)) {
      mapping[source] <- as.integer(target)
      used_reference <- c(used_reference, target)
    }
  }
  mapping
}

cluster_jaccard <- function(k, runs, reference_clusters = NULL) {
  set.seed(seed)
  result <- clusterboot(
    x,
    B = runs,
    clustermethod = kmeansCBI,
    k = k,
    runs = runs,
    count = FALSE,
    seed = seed
  )
  table <- data.table(K = k, cluster = seq_along(result$bootmean), jaccard = as.numeric(result$bootmean))
  if (!is.null(reference_clusters)) {
    mapping <- match_cluster_labels(result$partition, reference_clusters)
    table[, reference_cluster := unname(mapping[as.character(cluster)])]
  }
  table
}

feature_stability <- function(k, runs) {
  set.seed(seed + 1000L + k)
  base_clusters <- kmeans_fits[[as.character(k)]]$cluster
  feature_count <- max(2L, floor(feature_subsample_fraction * ncol(x)))
  ari <- replicate(runs, {
    columns <- sort(sample.int(ncol(x), feature_count, replace = FALSE))
    fit <- kmeans(x[, columns, drop = FALSE], centers = k, nstart = nstart, iter.max = 300)
    adjusted_rand_index(base_clusters, fit$cluster)
  })
  data.table(K = k, feature_subset_mean_ARI = mean(ari), feature_subset_sd_ARI = sd(ari))
}

candidate_k <- sort(unique(c(k_gap, k_silhouette)))
candidate_jaccard_detail <- rbindlist(lapply(candidate_k, cluster_jaccard, runs = jaccard_replicates))
candidate_feature <- rbindlist(lapply(candidate_k, feature_stability, runs = feature_replicates))
candidate_comparison <- data.table(
  K = candidate_k,
  candidate_mean_silhouette_width = unname(silhouette_values[as.character(candidate_k)]),
  calinski_harabasz = unname(ch_values[as.character(candidate_k)]),
  davies_bouldin = unname(db_values[as.character(candidate_k)]),
  mean_clusterwise_jaccard = vapply(
    candidate_k,
    function(k) mean(candidate_jaccard_detail[K == k]$jaccard, na.rm = TRUE),
    numeric(1)
  ),
  median_clusterwise_jaccard = vapply(
    candidate_k,
    function(k) median(candidate_jaccard_detail[K == k]$jaccard, na.rm = TRUE),
    numeric(1)
  ),
  proportion_clusters_jaccard_ge_0.75 = vapply(
    candidate_k,
    function(k) mean(candidate_jaccard_detail[K == k]$jaccard >= jaccard_stability_threshold, na.rm = TRUE),
    numeric(1)
  )
)
fwrite(candidate_comparison, file.path(output_dir, "clustering_candidate_clusterwise_jaccard.tsv"), sep = "\t")
fwrite(candidate_feature, file.path(output_dir, "clustering_candidate_feature_stability.tsv"), sep = "\t")

if (k_gap == k_silhouette) {
  consensus_k <- k_gap
  decision_reason <- sprintf("Gap firstSEmax and maximum silhouette agreed at K=%d.", consensus_k)
} else {
  lower_k <- min(k_gap, k_silhouette)
  higher_k <- max(k_gap, k_silhouette)
  lower <- diagnostics[K == lower_k]
  higher <- diagnostics[K == higher_k]
  lower_jaccard <- mean(candidate_jaccard_detail[K == lower_k]$jaccard, na.rm = TRUE)
  higher_jaccard <- mean(candidate_jaccard_detail[K == higher_k]$jaccard, na.rm = TRUE)
  lower_feature <- candidate_feature[K == lower_k]$feature_subset_mean_ARI
  higher_feature <- candidate_feature[K == higher_k]$feature_subset_mean_ARI
  separation_gain <-
    (higher$between_to_total_ss_ratio - lower$between_to_total_ss_ratio >= 0.02) ||
    (higher$mean_silhouette_width - lower$mean_silhouette_width >= 0.02) ||
    ((higher$calinski_harabasz - lower$calinski_harabasz) / lower$calinski_harabasz >= 0.10)
  stability_retained <-
    higher_jaccard >= lower_jaccard - 0.02 &&
    higher$row_subsample_mean_ARI >= lower$row_subsample_mean_ARI - 0.02 &&
    higher_feature >= lower_feature - 0.02
  if (separation_gain && stability_retained) {
    consensus_k <- higher_k
    decision_reason <- sprintf("Higher candidate K=%d improved separation while retaining resampling stability.", consensus_k)
  } else {
    consensus_k <- lower_k
    decision_reason <- sprintf("Lower candidate K=%d was retained because the higher-K separation gain or stability criterion was not met.", consensus_k)
  }
}

if (consensus_k != expected_consensus_k) {
  stop(
    sprintf("Expected consensus K=%d; selected K=%d.", expected_consensus_k, consensus_k),
    call. = FALSE
  )
}

set.seed(seed)
final_fit <- kmeans(x, centers = consensus_k, nstart = nstart, iter.max = 300)
centroid_order <- order(apply(final_fit$centers, 1L, function(row) paste(sprintf("%+.8f", row), collapse = "|")))
cluster_map <- setNames(seq_along(centroid_order), centroid_order)
stable_clusters <- unname(cluster_map[as.character(final_fit$cluster)])
assignments <- data.table(sample = scores$sample, cluster = stable_clusters, cluster_label = paste0("C", stable_clusters))
setorder(assignments, sample)
fwrite(assignments, file.path(output_dir, "final_cluster_assignments.tsv"), sep = "\t")

centroids_scaled <- as.data.table(x)
centroids_scaled[, cluster := stable_clusters]
centroids_scaled <- centroids_scaled[, lapply(.SD, mean), by = cluster, .SDcols = paste0("PC", seq_len(pc_count))]
centroids_scaled[, cluster_label := paste0("C", cluster)]
setcolorder(centroids_scaled, c("cluster", "cluster_label", paste0("PC", seq_len(pc_count))))
setorder(centroids_scaled, cluster)
fwrite(centroids_scaled, file.path(output_dir, "cluster_centroids_scaled_pc1_pc20.tsv"), sep = "\t")

centroid_matrix <- as.matrix(centroids_scaled[, paste0("PC", seq_len(pc_count)), with = FALSE])
rownames(centroid_matrix) <- centroids_scaled$cluster_label
centroid_hclust <- hclust(dist(centroid_matrix, method = "euclidean"), method = "ward.D2")
member_counts <- integer(nrow(centroid_hclust$merge))
node_label <- function(index) {
  if (index < 0L) rownames(centroid_matrix)[-index] else paste0("M", index)
}
node_members <- function(index) {
  if (index < 0L) 1L else member_counts[index]
}
hierarchy_rows <- vector("list", nrow(centroid_hclust$merge))
for (index in seq_len(nrow(centroid_hclust$merge))) {
  left <- centroid_hclust$merge[index, 1]
  right <- centroid_hclust$merge[index, 2]
  member_counts[index] <- node_members(left) + node_members(right)
  hierarchy_rows[[index]] <- data.table(
    merge_step = index,
    node = paste0("M", index),
    left_child = node_label(left),
    right_child = node_label(right),
    ward_height = centroid_hclust$height[index],
    n_clusters = member_counts[index]
  )
}
centroid_hierarchy <- rbindlist(hierarchy_rows)
fwrite(centroid_hierarchy, file.path(output_dir, "cluster_centroid_hierarchy_ward_d2.tsv"), sep = "\t")

final_jaccard_raw <- cluster_jaccard(consensus_k, final_jaccard_replicates, final_fit$cluster)
final_jaccard_raw[, stable_cluster := unname(cluster_map[as.character(reference_cluster)])]
cluster_sizes <- assignments[, .(size = .N), by = cluster]
final_jaccard <- final_jaccard_raw[, .(K, cluster = stable_cluster, jaccard)]
final_jaccard <- merge(final_jaccard, cluster_sizes, by = "cluster", all.x = TRUE, sort = FALSE)
final_jaccard[, cluster_label := paste0("C", cluster)]
final_jaccard[, stable_ge_0.75 := jaccard >= jaccard_stability_threshold]
setcolorder(final_jaccard, c("K", "cluster", "cluster_label", "size", "jaccard", "stable_ge_0.75"))
setorder(final_jaccard, cluster)
fwrite(final_jaccard, file.path(output_dir, "final_clusterwise_jaccard.tsv"), sep = "\t")

selection_summary <- data.table(
  criterion = c(
    "gap_firstSEmax", "maximum_mean_silhouette", "maximum_row_subsample_ARI",
    "consensus_k", "decision_reason"
  ),
  value = c(k_gap, k_silhouette, k_row_stability, consensus_k, decision_reason)
)
fwrite(selection_summary, file.path(output_dir, "clustering_selection_summary.tsv"), sep = "\t")

run_parameters <- data.table(
  parameter = c(
    "input_cohort_samples", "input_LD_pruned_SNPs", "reported_PCs", "clustering_PCs",
    "PCs_centered", "PCs_scaled", "k_min", "k_max", "seed",
    "kmeans_nstart", "gap_reference_replicates", "row_subsample_replicates_per_K",
    "row_subsample_fraction", "row_subsample_with_replacement",
    "candidate_jaccard_bootstrap_replicates", "candidate_clusterboot_kmeans_starts",
    "feature_subsample_replicates", "feature_subsample_fraction",
    "final_jaccard_bootstrap_replicates", "final_clusterboot_kmeans_starts",
    "clusterwise_jaccard_stability_threshold"
  ),
  value = c(
    nrow(scores), input_snp_count, reported_pc_count, pc_count, TRUE, TRUE,
    k_min, k_max, seed, nstart, gap_replicates, row_replicates,
    row_subsample_fraction, FALSE, jaccard_replicates, jaccard_replicates,
    feature_replicates, feature_subsample_fraction,
    final_jaccard_replicates, final_jaccard_replicates,
    jaccard_stability_threshold
  )
)
fwrite(run_parameters, file.path(output_dir, "clustering_run_parameters.tsv"), sep = "\t")
