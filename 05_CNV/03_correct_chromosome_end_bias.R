#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Filter masked or weakly covered 1 kb bins and correct the remaining depth
#   for distance-to-chromosome-end effects by held-out LOESS prediction.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

input_path <- get_arg("--input")
mask_path <- get_arg("--mask")
output_path <- get_arg("--output")
k_folds <- as.integer(get_arg("--folds", "10"))
loess_span <- as.numeric(get_arg("--span", "0.3"))
minimum_fraction_covered <- as.numeric(get_arg("--minimum-fraction-covered", "0.75"))
minimum_median_matched <- as.numeric(get_arg("--minimum-median-matched", "80"))

if (is.null(input_path) || is.null(mask_path) || is.null(output_path)) {
  stop("--input, --mask and --output are required.", call. = FALSE)
}
if (k_folds < 2L) stop("--folds must be at least 2.", call. = FALSE)
if (!is.finite(loess_span) || loess_span <= 0 || loess_span > 1) {
  stop("--span must be in (0, 1].", call. = FALSE)
}

depth <- fread(input_path)
required <- c(
  "chr", "start", "end", "reads_overlapping",
  "median_matched", "fraction_covered"
)
missing_columns <- setdiff(required, names(depth))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

numeric_columns <- setdiff(required, "chr")
depth[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
depth[, chromosome_length := max(end), by = chr]
depth <- depth[
  !grepl("chrM", chr, ignore.case = TRUE) &
    is.finite(start) &
    is.finite(end) &
    is.finite(reads_overlapping) &
    is.finite(median_matched) &
    is.finite(fraction_covered) &
    end > start &
    fraction_covered >= minimum_fraction_covered &
    median_matched >= minimum_median_matched
]

mask <- if (file.info(mask_path)$size == 0) data.table(chr = character(), start = integer(), end = integer()) else fread(
  mask_path,
  header = FALSE,
  select = 1:3,
  col.names = c("chr", "start", "end")
)
mask[, `:=`(start = as.integer(start), end = as.integer(end))]
mask <- mask[is.finite(start) & is.finite(end) & end > start]

depth[, row_id := .I]
depth_intervals <- depth[, .(
  chr,
  interval_start = as.integer(start),
  interval_end = as.integer(end) - 1L,
  row_id
)]
mask_intervals <- mask[, .(
  chr,
  interval_start = as.integer(start),
  interval_end = as.integer(end) - 1L
)]
setkey(depth_intervals, chr, interval_start, interval_end)
setkey(mask_intervals, chr, interval_start, interval_end)
overlaps <- foverlaps(depth_intervals, mask_intervals, nomatch = 0L)
if (nrow(overlaps)) depth <- depth[!row_id %in% unique(overlaps$row_id)]
depth[, row_id := NULL]

setorder(depth, chr, start, end)
if (!nrow(depth)) stop("No depth bins pass the coverage and mask filters.")
depth[, distance_to_end := pmax(0, pmin(start, chromosome_length - end))]
depth[, fold := ((seq_len(.N) - 1L) %% k_folds) + 1L]
depth[, expected_depth := NA_real_]

predict_expected <- function(training, held_out) {
  good <- is.finite(training$reads_overlapping) &
    is.finite(training$distance_to_end)
  fallback <- median(training$reads_overlapping[good], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) fallback <- 1
  if (sum(good) < 10L || uniqueN(training$distance_to_end[good]) < 2L) {
    return(rep(fallback, nrow(held_out)))
  }

  model <- tryCatch(
    loess(
      log(reads_overlapping + 1) ~ distance_to_end,
      data = training[good],
      span = loess_span,
      degree = 2,
      family = "gaussian",
      control = loess.control(surface = "direct")
    ),
    error = function(condition) NULL
  )
  if (is.null(model)) return(rep(fallback, nrow(held_out)))

  prediction <- tryCatch(
    predict(model, newdata = data.frame(distance_to_end = held_out$distance_to_end)),
    error = function(condition) rep(NA_real_, nrow(held_out))
  )
  predicted_depth <- exp(prediction) - 1
  predicted_depth[!is.finite(predicted_depth) | predicted_depth <= 0] <- fallback
  predicted_depth
}

for (fold_number in seq_len(k_folds)) {
  held_out <- depth[fold == fold_number]
  depth[
    fold == fold_number,
    expected_depth := predict_expected(depth[fold != fold_number], held_out)
  ]
}

original_median <- median(depth$reads_overlapping, na.rm = TRUE)
depth[, reads_overlapping_end_corrected :=
  reads_overlapping / pmax(expected_depth, 1e-6)]
corrected_median <- median(depth$reads_overlapping_end_corrected, na.rm = TRUE)
if (is.finite(corrected_median) && corrected_median > 0) {
  depth[, reads_overlapping_end_corrected :=
    reads_overlapping_end_corrected * original_median / corrected_median]
}

depth[, `:=`(
  end_correction_folds = k_folds,
  end_correction_span = loess_span,
  minimum_fraction_covered = minimum_fraction_covered,
  minimum_median_matched = minimum_median_matched
)]
depth[, c("chromosome_length", "fold", "expected_depth") := NULL]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(depth, output_path, sep = "\t", quote = FALSE, na = "NA")
