#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Compare K-means solutions with sampling labels, test PC-label associations,
#   and summarize the relationship between final clusters and curated clades.
#
# Inputs
#   Final PCA/cluster table and sample annotation table supplied through
#   command-line options.
#
# Outputs
#   Long-form cluster-label agreement, selected-K summaries and confusion counts,
#   PC-label ANOVA results, cluster-clade counts, and run parameters.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_path)), "helpers", "cluster_label_tables.R"))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

scores_path <- get_arg("--scores", "path/to/pca_scores_and_cluster_assignments.tsv")
metadata_path <- get_arg("--metadata", "path/to/phylogeny_figure_sample_metadata.tsv")
output_dir <- get_arg("--output-dir", "path/to/population_structure_comparison")
expected_samples <- as.integer(get_arg("--expected-samples", NA_character_))
final_k <- as.integer(get_arg("--final-k", NA_character_))
k_min <- as.integer(get_arg("--k-min", "2"))
k_max_requested <- get_arg("--k-max")
pc_count_requested <- get_arg("--clustering-pcs")
heatmaps_only <- "--heatmaps-only" %in% args
seed <- as.integer(get_arg("--seed", "42"))
nstart <- as.integer(get_arg("--nstart", "150"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(scores_path, na.strings = c("", "NA", "NaN"), colClasses = c(sample = "character", cluster = "character"))
metadata <- fread(metadata_path, na.strings = c("", "NA", "NaN"), colClasses = "character")
available_pcs <- grep("^PC[0-9]+$", names(scores), value = TRUE)
pc_count <- if (is.null(pc_count_requested)) min(20L, length(available_pcs)) else as.integer(pc_count_requested)
if (is.na(pc_count) || pc_count < 2L) stop("At least PC1 and PC2 are required.", call. = FALSE)
required_score_columns <- c("sample", "cluster", paste0("PC", seq_len(if (heatmaps_only) 2L else pc_count)))
required_metadata_columns <- c("sample", "Loc", "Study_clean", "rbin_flu", "MLST_clean", "Clade")
missing_scores <- setdiff(required_score_columns, names(scores))
missing_metadata <- setdiff(required_metadata_columns, names(metadata))
if (length(missing_scores)) stop("PCA table lacks required columns: ", paste(missing_scores, collapse = ", "), call. = FALSE)
if (length(missing_metadata)) stop("Metadata table lacks required columns: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
if (!is.na(expected_samples) && nrow(scores) != expected_samples) stop("Unexpected PCA sample count.", call. = FALSE)
if (!is.na(expected_samples) && nrow(metadata) != expected_samples) stop("Unexpected metadata sample count.", call. = FALSE)
if (anyDuplicated(scores$sample)) stop("Duplicate PCA sample identifiers found.", call. = FALSE)
if (anyDuplicated(metadata$sample)) stop("Duplicate metadata sample identifiers found.", call. = FALSE)
if (!setequal(scores$sample, metadata$sample)) stop("PCA and metadata sample sets differ.", call. = FALSE)
if (!nrow(scores) || anyNA(scores$cluster) || any(!nzchar(as.character(scores$cluster)))) stop("Cluster identifiers must be nonmissing.", call. = FALSE)
observed_k <- uniqueN(scores$cluster)
if (is.na(final_k)) final_k <- observed_k
if (final_k != observed_k) stop("Final K differs from the number of observed clusters.", call. = FALSE)
final_clusters <- as.character(scores$cluster)

metadata <- metadata[match(scores$sample, sample)]
if (!identical(scores$sample, metadata$sample)) stop("Sample-order reconciliation failed.", call. = FALSE)

fluconazole_label <- normalize_label(metadata$rbin_flu)
fluconazole_label[fluconazole_label %in% c("0", "0.0", "S", "Susceptible")] <- "Susceptible"
fluconazole_label[fluconazole_label %in% c("1", "1.0", "R", "Resistant")] <- "Resistant"
fluconazole_label[!fluconazole_label %in% c("Resistant", "Susceptible")] <- NA_character_

label_values <- list(
  published_MLST_clade = normalize_label(metadata$MLST_clean),
  geographic_origin = normalize_label(metadata$Loc),
  contributing_study = normalize_label(metadata$Study_clean),
  fluconazole_resistance = fluconazole_label
)

adjusted_rand_index <- function(partition_a, partition_b) {
  if (length(partition_a) < 2L) return(NA_real_)
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

cluster_purity <- function(clusters, labels) {
  if (!length(labels)) return(NA_real_)
  contingency <- table(clusters, labels)
  sum(apply(contingency, 1L, max)) / sum(contingency)
}

normalized_mutual_information <- function(clusters, labels) {
  if (!length(labels)) return(NA_real_)
  contingency <- table(clusters, labels)
  joint <- contingency / sum(contingency)
  row_prob <- rowSums(joint)
  column_prob <- colSums(joint)
  expected <- outer(row_prob, column_prob)
  nonzero <- joint > 0
  mutual_information <- sum(joint[nonzero] * log(joint[nonzero] / expected[nonzero]))
  cluster_entropy <- -sum(row_prob[row_prob > 0] * log(row_prob[row_prob > 0]))
  label_entropy <- -sum(column_prob[column_prob > 0] * log(column_prob[column_prob > 0]))
  denominator <- (cluster_entropy + label_entropy) / 2
  if (denominator == 0) return(1)
  mutual_information / denominator
}

pc_names <- paste0("PC", seq_len(pc_count))
if (!heatmaps_only) {
  x <- scale(as.matrix(scores[, ..pc_names]), center = TRUE, scale = TRUE)
  if (any(!is.finite(x))) stop("Non-finite values encountered after PC scaling.", call. = FALSE)

  set.seed(seed)
  k_max <- if (is.null(k_max_requested)) min(30L, nrow(unique(as.data.frame(x))) - 1L) else as.integer(k_max_requested)
  if (is.na(k_min) || is.na(k_max) || k_min < 2L || k_max < k_min || k_max >= nrow(unique(as.data.frame(x)))) stop("K range must be between 2 and one less than the number of distinct PC vectors; use --heatmaps-only for descriptive exports.", call. = FALSE)
  k_values <- k_min:k_max
  kmeans_fits <- setNames(vector("list", length(k_values)), as.character(k_values))
  for (k in k_values) {
    kmeans_fits[[as.character(k)]] <- kmeans(x, centers = k, nstart = nstart, iter.max = 300)
  }
  if (final_k %in% k_values) kmeans_fits[[as.character(final_k)]]$cluster <- final_clusters

  agreement_by_k <- rbindlist(lapply(names(label_values), function(label_name) {
    labels <- label_values[[label_name]]
    keep <- !is.na(labels)
    rbindlist(lapply(k_values, function(k) {
      clusters <- kmeans_fits[[as.character(k)]]$cluster[keep]
      data.table(
        label = label_name,
        K = k,
        n_complete = sum(keep),
        n_categories = uniqueN(labels[keep]),
        adjusted_rand_index = adjusted_rand_index(clusters, labels[keep])
      )
    }))
  }))
  fwrite(agreement_by_k, file.path(output_dir, "cluster_label_agreement_by_k.tsv"), sep = "\t")

  agreement_final <- rbindlist(lapply(names(label_values), function(label_name) {
    labels <- label_values[[label_name]]
    keep <- !is.na(labels)
    data.table(
      label = label_name,
      K = final_k,
      n_complete = sum(keep),
      n_categories = uniqueN(labels[keep]),
      adjusted_rand_index = adjusted_rand_index(final_clusters[keep], labels[keep]),
      cluster_purity = cluster_purity(final_clusters[keep], labels[keep]),
      normalized_mutual_information_arithmetic = normalized_mutual_information(final_clusters[keep], labels[keep])
    )
  }))
  fwrite(agreement_final, file.path(output_dir, paste0("cluster_label_agreement_k", final_k, ".tsv")), sep = "\t")

}
confusion_final <- rbindlist(lapply(names(label_values), function(label_name) {
  result <- as.data.table(cluster_metadata_table(final_clusters, label_values[[label_name]], metadata$Clade, label_name))
  result[, K := final_k]
  result
}))
fwrite(confusion_final, file.path(output_dir, paste0("cluster_label_confusion_k", final_k, ".tsv")), sep = "\t", na = "NA")

figure2_samples <- data.table(
  sample = scores$sample,
  PC1 = scores$PC1,
  PC2 = scores$PC2,
  cluster = final_clusters,
  cluster_label = paste0("C", final_clusters),
  published_MLST_clade = label_values$published_MLST_clade,
  extrapolated_clade = normalize_label(metadata$Clade),
  fluconazole_resistance = fluconazole_label,
  geographic_origin = label_values$geographic_origin,
  contributing_study = label_values$contributing_study
)
figure2_samples[is.na(extrapolated_clade), extrapolated_clade := "Other"]
fwrite(figure2_samples, file.path(output_dir, "figure2_sample_data.tsv"), sep = "\t", na = "")

clade_metadata <- confusion_final[label != "published_MLST_clade"]
fwrite(clade_metadata, file.path(output_dir, "figure2_clade_metadata_heatmaps.tsv"), sep = "\t", na = "NA")

if (!heatmaps_only) {
  anova_results <- rbindlist(lapply(names(label_values), function(label_name) {
    labels <- label_values[[label_name]]
    keep <- !is.na(labels)
    result <- rbindlist(lapply(pc_names, function(pc_name) {
      response <- scores[[pc_name]][keep]
      group <- factor(labels[keep])
      if (nlevels(group) < 2L || length(response) <= nlevels(group)) return(data.table(
        label = label_name, PC = pc_name, n_complete = length(response), n_categories = nlevels(group),
        df_between = NA_real_, df_within = NA_real_, F_statistic = NA_real_, p_value = NA_real_, eta_squared = NA_real_))
      fit <- aov(response ~ group)
      table <- summary(fit)[[1L]]
      between_ss <- table[1L, "Sum Sq"]
      total_ss <- sum(table[, "Sum Sq"])
      data.table(
        label = label_name,
        PC = pc_name,
        n_complete = length(response),
        n_categories = nlevels(group),
        df_between = table[1L, "Df"],
        df_within = table[2L, "Df"],
        F_statistic = table[1L, "F value"],
        p_value = table[1L, "Pr(>F)"],
        eta_squared = between_ss / total_ss
      )
    }))
    result[, Benjamini_Hochberg_FDR := p.adjust(p_value, method = "BH")]
    result[, significant_FDR_0.05 := Benjamini_Hochberg_FDR < 0.05]
    result
  }))
  fwrite(anova_results, file.path(output_dir, "pc_label_anova.tsv"), sep = "\t")
}

clade <- normalize_label(metadata$Clade)
clade[is.na(clade)] <- "Other"
cluster_clade <- data.table(cluster = final_clusters, clade = clade)[, .(count = .N), by = .(cluster, clade)]
cluster_clade[, cluster_total := sum(count), by = cluster]
cluster_clade[, within_cluster_fraction := count / cluster_total]
cluster_clade[, is_majority_clade := count == max(count), by = cluster]
cluster_clade[, `:=`(K = final_k, cluster_label = paste0("C", cluster))]
setorder(cluster_clade, cluster, -count, clade)
setcolorder(cluster_clade, c("K", "cluster", "cluster_label", "clade", "count", "cluster_total", "within_cluster_fraction", "is_majority_clade"))
fwrite(cluster_clade, file.path(output_dir, "cluster_clade_mapping.tsv"), sep = "\t")

if (!heatmaps_only) {
  run_parameters <- data.table(
    parameter = c(
      "input_cohort_samples", "clustering_PCs", "k_min", "k_max", "final_K",
      "seed", "kmeans_nstart", "NMI_normalization", "ANOVA_PCs",
      "multiple_testing_correction", "multiple_testing_scope"
    ),
    value = c(
      nrow(scores), pc_count, k_min, k_max, final_k, seed, nstart,
      "mutual information divided by arithmetic mean of the two entropies",
      paste0("PC1-PC", pc_count), "Benjamini-Hochberg", "within each external label across PCs"
    )
  )
  fwrite(run_parameters, file.path(output_dir, "cluster_label_run_parameters.tsv"), sep = "\t")

}
