#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Parse ADMIXTURE cross-validation logs and export replicate-level and K-level
#   source-data tables after verifying the final cohort and SNP dimensions.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

input_dir <- get_arg("--input-dir", "path/to/admixture_cv")
output_dir <- get_arg("--output-dir", "path/to/admixture_cv/summary")
expected_samples <- as.integer(get_arg("--expected-samples", "1153"))
expected_snps <- as.integer(get_arg("--expected-snps", "23356"))
k_min <- as.integer(get_arg("--k-min", "2"))
k_max <- as.integer(get_arg("--k-max", "30"))
expected_replicates <- as.integer(get_arg("--replicates", "10"))

extract_one <- function(log_path) {
  lines <- readLines(log_path, warn = FALSE)
  k_match <- regexec("/K([0-9]+)/rep([0-9]+)/run\\.log$", log_path)
  k_parts <- regmatches(log_path, k_match)[[1L]]
  if (length(k_parts) != 3L) stop("Unexpected log path: ", log_path, call. = FALSE)

  seed_line <- grep("^Random seed:", lines, value = TRUE)
  size_line <- grep("^Size of G:", lines, value = TRUE)
  cv_line <- grep("^CV error \\(K=[0-9]+\\):", lines, value = TRUE)
  if (length(seed_line) != 1L || length(size_line) != 1L || length(cv_line) != 1L) {
    stop("Incomplete ADMIXTURE log: ", log_path, call. = FALSE)
  }

  size_match <- regexec("Size of G: ([0-9]+)x([0-9]+)", size_line)
  size_parts <- regmatches(size_line, size_match)[[1L]]
  cv_value <- as.numeric(sub(".*: *", "", cv_line))
  data.table(
    K = as.integer(k_parts[2L]),
    replicate = as.integer(k_parts[3L]),
    seed = as.integer(sub("Random seed: *", "", seed_line)),
    n_samples = as.integer(size_parts[2L]),
    n_snps = as.integer(size_parts[3L]),
    cv_error = cv_value
  )
}

logs <- list.files(
  input_dir,
  pattern = "run\\.log$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(logs)) stop("No ADMIXTURE logs found.", call. = FALSE)

replicates <- rbindlist(lapply(logs, extract_one))
setorder(replicates, K, replicate)
expected_k <- k_min:k_max
if (!identical(sort(unique(replicates$K)), expected_k)) stop("K range is incomplete.", call. = FALSE)
if (replicates[, any(.N != expected_replicates), by = K]$V1) stop("Replicate count is incomplete.", call. = FALSE)
if (anyDuplicated(replicates[, .(K, replicate)])) stop("Duplicate K/replicate logs found.", call. = FALSE)
if (any(replicates$n_samples != expected_samples)) stop("Unexpected ADMIXTURE sample count.", call. = FALSE)
if (any(replicates$n_snps != expected_snps)) stop("Unexpected ADMIXTURE SNP count.", call. = FALSE)

summary_by_k <- replicates[, {
  best_index <- which.min(cv_error)
  .(
    n_replicates = .N,
    mean_cv_error = mean(cv_error),
    sd_cv_error = sd(cv_error),
    median_cv_error = median(cv_error),
    minimum_cv_error = min(cv_error),
    maximum_cv_error = max(cv_error),
    best_replicate = replicate[best_index],
    best_seed = seed[best_index]
  )
}, by = K]

cv_at_k2 <- summary_by_k[K == k_min, mean_cv_error]
cv_at_kmax <- summary_by_k[K == k_max, mean_cv_error]
summary_by_k[, fraction_total_mean_CV_reduction :=
  (cv_at_k2 - mean_cv_error) / (cv_at_k2 - cv_at_kmax)]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(replicates, file.path(output_dir, "admixture_cv_replicates.tsv"), sep = "\t")
fwrite(summary_by_k, file.path(output_dir, "admixture_cv_summary.tsv"), sep = "\t")

