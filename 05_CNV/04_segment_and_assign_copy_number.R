#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Segment corrected read depth and heterozygous-SNP allele fractions, infer a
#   candidate baseline ploidy, and assign integer segment, chromosome, and gene
#   copy numbers for manual review.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(DNAcopy)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[[index + 1L]]
}

depth_path <- get_arg("--depth")
baf_path <- get_arg("--baf")
gene_path <- get_arg("--genes")
output_prefix <- get_arg("--output-prefix")
segment_alpha <- as.numeric(get_arg("--segment-alpha", "0.01"))
minimum_segment_bins <- as.integer(get_arg("--minimum-segment-bins", "3"))
minimum_baf_sites <- as.integer(get_arg("--minimum-baf-sites", "10"))
reviewed_baseline <- as.integer(get_arg("--baseline-ploidy", NA_character_))
segment_override_path <- get_arg("--segment-overrides")
seed <- as.integer(get_arg("--seed", "42"))
set.seed(seed)
if (!is.na(reviewed_baseline) && !(reviewed_baseline %in% c(2L, 3L, 4L, 6L))) {
  stop("--baseline-ploidy must be 2, 3, 4 or 6.", call. = FALSE)
}
if (!is.finite(segment_alpha) || segment_alpha <= 0 || segment_alpha >= 1 ||
    minimum_segment_bins < 2 || minimum_baf_sites < 1) stop("Invalid segmentation thresholds.")

if (any(vapply(
  list(depth_path, baf_path, gene_path, output_prefix),
  is.null,
  logical(1)
))) {
  stop(
    "--depth, --baf, --genes and --output-prefix are required.",
    call. = FALSE
  )
}

depth <- fread(depth_path)
baf <- fread(baf_path)
genes <- fread(gene_path)

required_depth <- c("chr", "start", "end", "reads_overlapping_end_corrected")
required_baf <- c("chr", "pos", "ref_count", "alt_count")
required_genes <- c("chr", "start", "end", "gene_id")
for (entry in list(
  list(table = depth, columns = required_depth, name = "depth"),
  list(table = baf, columns = required_baf, name = "BAF"),
  list(table = genes, columns = required_genes, name = "genes")
)) {
  missing_columns <- setdiff(entry$columns, names(entry$table))
  if (length(missing_columns)) {
    stop(
      entry$name,
      " input is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

depth[, `:=`(
  start = as.integer(start),
  end = as.integer(end),
  corrected_depth = as.numeric(reads_overlapping_end_corrected)
)]
depth <- depth[
  !grepl("chrM", chr, ignore.case = TRUE) &
    is.finite(start) &
    is.finite(end) &
    is.finite(corrected_depth) &
    corrected_depth > 0 &
    end > start
]
depth[, midpoint := floor((start + end) / 2)]
if (!nrow(depth)) stop("No finite, positive depth bins remain.")
setorder(depth, chr, start, end)
if (depth[, any(start < shift(end, fill = -1L)), by = chr][, any(V1)]) {
  stop("Depth intervals must not overlap.")
}
genome_median <- median(depth$corrected_depth, na.rm = TRUE)
depth[, log2_ratio := log2(corrected_depth / genome_median)]

baf[, `:=`(
  pos = as.integer(pos),
  ref_count = as.numeric(ref_count),
  alt_count = as.numeric(alt_count)
)]
baf <- baf[
  !grepl("chrM", chr, ignore.case = TRUE) &
    is.finite(pos) &
    is.finite(ref_count) &
    is.finite(alt_count) &
    pos >= 1 & ref_count > 0 & alt_count > 0 &
    ref_count + alt_count > 0
]
if (anyDuplicated(baf[, .(chr, pos)])) stop("Repeated SNP coordinates in BAF input.")
setorder(baf, chr, pos)
baf[, pos0 := pos - 1L]
baf[, allele_fraction := alt_count / (ref_count + alt_count)]
baf[, minor_allele_fraction := pmin(allele_fraction, 1 - allele_fraction)]
baf[, baf_deviation := abs(allele_fraction - 0.5)]

segment_track <- function(values, chromosomes, positions) {
  output <- lapply(unique(chromosomes), function(chromosome) {
    use <- chromosomes == chromosome
    vals <- values[use]
    locs <- positions[use]
    if (length(vals) < 2L * minimum_segment_bins || uniqueN(vals) == 1L) {
      return(data.table(chrom = chromosome, loc.start = min(locs), loc.end = max(locs)))
    }
    cna <- CNA(genomdat = vals, chrom = rep(chromosome, length(vals)), maploc = locs, data.type = "logratio")
    as.data.table(segment(smooth.CNA(cna), alpha = segment_alpha,
      min.width = minimum_segment_bins, verbose = 0)$output)[, .(chrom, loc.start, loc.end)]
  })
  if (!length(output)) return(data.table(chrom = character(), loc.start = numeric(), loc.end = numeric()))
  rbindlist(output)
}

depth_segments <- segment_track(depth$log2_ratio, depth$chr, depth$midpoint)
baf_segments <- segment_track(baf$baf_deviation, baf$chr, baf$pos0)

internal_breaks <- function(table, chromosome) {
  track <- table[as.character(chrom) == chromosome][order(loc.start)]
  if (nrow(track) < 2L) return(numeric())
  floor((head(track$loc.end, -1L) + tail(track$loc.start, -1L)) / 2)
}

build_consensus_intervals <- function(chromosome) {
  chromosome_depth <- depth[chr == chromosome]
  if (!nrow(chromosome_depth)) return(data.table())
  chromosome_start <- min(chromosome_depth$start)
  chromosome_end <- max(chromosome_depth$end)

  depth_breaks <- internal_breaks(depth_segments, chromosome)
  baf_breaks <- internal_breaks(baf_segments, chromosome)
  breaks <- sort(unique(c(
    chromosome_start,
    depth_breaks[depth_breaks > chromosome_start & depth_breaks < chromosome_end],
    baf_breaks[baf_breaks > chromosome_start & baf_breaks < chromosome_end],
    chromosome_end
  )))
  if (length(breaks) < 2L) return(data.table())
  data.table(
    chr = chromosome,
    start = head(breaks, -1L),
    end = tail(breaks, -1L)
  )[end > start]
}

consensus <- rbindlist(
  lapply(unique(depth$chr), build_consensus_intervals),
  use.names = TRUE
)
consensus[, segment_id := sprintf("segment_%06d", .I)]

summarize_segment <- function(chromosome, start_value, end_value) {
  depth_subset <- depth[
    chr == chromosome & midpoint >= start_value & midpoint < end_value
  ]
  baf_subset <- baf[
    chr == chromosome & pos0 >= start_value & pos0 < end_value
  ]
  list(
    n_depth_bins = nrow(depth_subset),
    median_log2_ratio = median(depth_subset$log2_ratio, na.rm = TRUE),
    n_informative_snps = nrow(baf_subset),
    median_allele_fraction = median(baf_subset$allele_fraction, na.rm = TRUE),
    median_minor_allele_fraction =
      median(baf_subset$minor_allele_fraction, na.rm = TRUE)
  )
}

segment_summaries <- rbindlist(lapply(seq_len(nrow(consensus)), function(index) {
  row <- consensus[index]
  as.data.table(summarize_segment(row$chr, row$start, row$end))
}))
consensus <- cbind(consensus, segment_summaries)

candidate_ploidies <- c(2L, 3L, 4L, 6L)
informative_minor_af <- baf[
  is.finite(minor_allele_fraction) &
    minor_allele_fraction >= 0.05 &
    minor_allele_fraction <= 0.5,
  minor_allele_fraction
]

ploidy_score <- function(ploidy) {
  expected <- unique(pmin(
    seq_len(ploidy - 1L) / ploidy,
    1 - seq_len(ploidy - 1L) / ploidy
  ))
  if (!length(informative_minor_af)) return(Inf)
  mean(vapply(
    informative_minor_af,
    function(value) min((value - expected)^2),
    numeric(1)
  ))
}

scores <- vapply(candidate_ploidies, ploidy_score, numeric(1))
candidate_baseline <- if (any(is.finite(scores))) candidate_ploidies[which.min(scores)] else NA_integer_
if (is.na(candidate_baseline) && is.na(reviewed_baseline)) {
  stop("Baseline ploidy is unidentifiable without informative SNPs; supply a reviewed --baseline-ploidy.")
}
baseline_ploidy <- if (is.na(reviewed_baseline)) candidate_baseline else reviewed_baseline
matching_ploidies <- if (any(is.finite(scores))) candidate_ploidies[abs(scores - min(scores)) <= 1e-8] else candidate_ploidies

expected_minor_af <- function(copy_number) {
  if (copy_number < 2L) return(numeric())
  unique(pmin(
    seq_len(copy_number - 1L) / copy_number,
    1 - seq_len(copy_number - 1L) / copy_number
  ))
}

copy_number_score <- function(
  copy_number,
  median_log2_ratio,
  median_minor_af,
  n_snp
) {
  expected_depth <- if (copy_number == 0L) {
    -4
  } else {
    log2(copy_number / baseline_ploidy)
  }
  score <- (median_log2_ratio - expected_depth)^2
  expected_baf <- expected_minor_af(copy_number)
  if (
    n_snp >= minimum_baf_sites &&
      is.finite(median_minor_af) &&
      length(expected_baf)
  ) {
    score <- score + min((median_minor_af - expected_baf)^2)
  }
  score
}

consensus[, integer_copy_number := vapply(seq_len(.N), function(index) {
  if (!is.finite(median_log2_ratio[index])) return(NA_integer_)
  candidate_copy_numbers <- 0:8
  candidate_scores <- vapply(
    candidate_copy_numbers,
    copy_number_score,
    numeric(1),
    median_log2_ratio = median_log2_ratio[index],
    median_minor_af = median_minor_allele_fraction[index],
    n_snp = n_informative_snps[index]
  )
  candidate_copy_numbers[which.min(candidate_scores)]
}, integer(1))]

consensus[, assignment_method := "depth_and_allele_fraction"]
if (!is.null(segment_override_path)) {
  overrides <- fread(segment_override_path)
  required <- c("segment_id", "integer_copy_number", "reviewer", "evidence")
  if (!all(required %in% names(overrides)) || anyDuplicated(overrides$segment_id)) {
    stop("Segment overrides require unique segment_id, integer_copy_number, reviewer and evidence columns.")
  }
  if (any(!overrides$segment_id %in% consensus$segment_id) ||
      any(!is.finite(overrides$integer_copy_number) | overrides$integer_copy_number < 0 |
        overrides$integer_copy_number != floor(overrides$integer_copy_number)) ||
      any(is.na(overrides$reviewer) | !nzchar(overrides$reviewer) | is.na(overrides$evidence) | !nzchar(overrides$evidence))) {
    stop("Invalid segment override or missing review evidence.")
  }
  consensus[overrides, on = "segment_id", `:=`(integer_copy_number = i.integer_copy_number, assignment_method = "manual_review")]
}

consensus[, nearest_expected_minor_af := vapply(
  seq_len(.N),
  function(index) {
    if (is.na(integer_copy_number[index])) return(NA_real_)
    states <- expected_minor_af(integer_copy_number[index])
    observed <- median_minor_allele_fraction[index]
    if (!length(states) || !is.finite(observed)) return(NA_real_)
    states[which.min(abs(observed - states))]
  },
  numeric(1)
)]

consensus[, segment_length_bp := end - start]
chromosome_calls <- consensus[, {
  weights <- segment_length_bp[is.finite(integer_copy_number)]
  states <- integer_copy_number[is.finite(integer_copy_number)]
  totals <- tapply(weights, states, sum)
  .(
    baseline_ploidy = baseline_ploidy,
    chromosome_copy_number = if (length(totals)) as.integer(names(totals)[which.max(totals)]) else NA_integer_,
    chromosome_length_bp = sum(weights)
  )
}, by = chr]

genes[, `:=`(
  start = as.integer(start),
  end = as.integer(end)
)]
gene_calls <- genes[, {
  overlaps <- consensus[
    chr == .BY$chr & start < .BY$end & end > .BY$start
  ]
  if (!nrow(overlaps)) {
    list(
      segment_id = NA_character_,
      overlap_fraction = 0,
      integer_copy_number = NA_integer_
    )
  } else {
    overlap_bp <- pmax(
      0,
      pmin(overlaps$end, .BY$end) - pmax(overlaps$start, .BY$start)
    )
    best <- which.max(overlap_bp)
    list(
      segment_id = overlaps$segment_id[best],
      overlap_fraction = overlap_bp[best] / (.BY$end - .BY$start),
      integer_copy_number = if (
        overlap_bp[best] / (.BY$end - .BY$start) >= 0.5
      ) overlaps$integer_copy_number[best] else NA_integer_
    )
  }
}, by = .(gene_id, chr, start, end)]

ploidy_table <- data.table(
  baseline_ploidy = baseline_ploidy,
  candidate_baseline_ploidy = candidate_baseline,
  score_2N = scores[candidate_ploidies == 2L],
  score_3N = scores[candidate_ploidies == 3L],
  score_4N = scores[candidate_ploidies == 4L],
  score_6N = scores[candidate_ploidies == 6L],
  n_informative_snps = length(informative_minor_af),
  equally_scoring_ploidies = paste(matching_ploidies, collapse = ";"),
  baseline_scoring_method = "nearest_expected_allele_fraction_squared_error",
  status = if (!is.na(reviewed_baseline)) "reviewed_baseline" else if (length(matching_ploidies) > 1L) "ambiguous_candidate_for_manual_review" else "candidate_for_manual_review",
  seed = seed
)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(consensus, paste0(output_prefix, ".segments.tsv"), sep = "\t", na = "NA")
fwrite(chromosome_calls, paste0(output_prefix, ".chromosome_cn.tsv"), sep = "\t")
fwrite(gene_calls, paste0(output_prefix, ".gene_cn.tsv"), sep = "\t", na = "NA")
fwrite(ploidy_table, paste0(output_prefix, ".ploidy.tsv"), sep = "\t")
