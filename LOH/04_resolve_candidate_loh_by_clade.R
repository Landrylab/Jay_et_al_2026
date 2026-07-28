#!/usr/bin/env Rscript

# Resolve recurrent candidate synapomorphic LOH tracts within annotated clades.

options(stringsAsFactors = FALSE)

parse_args <- function(values) {
  if (length(values) %% 2 != 0) {
    stop("Arguments must be supplied as --name value pairs.")
  }
  output <- list()
  for (index in seq(1, length(values), by = 2)) {
    key <- sub("^--", "", values[index])
    output[[key]] <- values[index + 1]
  }
  output
}

required_value <- function(arguments, name) {
  value <- arguments[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

as_flag <- function(values) {
  toupper(trimws(as.character(values))) %in% c("TRUE", "T", "1")
}

merge_covered_bp <- function(start0, end0) {
  if (!length(start0)) return(0)
  ordering <- order(start0, end0)
  starts <- start0[ordering]
  ends <- end0[ordering]
  current_start <- starts[1]
  current_end <- ends[1]
  total <- 0
  if (length(starts) > 1) {
    for (index in 2:length(starts)) {
      if (starts[index] <= current_end) {
        current_end <- max(current_end, ends[index])
      } else {
        total <- total + current_end - current_start
        current_start <- starts[index]
        current_end <- ends[index]
      }
    }
  }
  total + current_end - current_start
}

covered_bp_by_sample <- function(frame) {
  if (!nrow(frame)) return(0)
  sum(vapply(
    split(frame, frame$contig),
    function(part) merge_covered_bp(part$start0, part$end0),
    numeric(1)
  ))
}

same_locus <- function(
    first_start, first_end, second_start, second_end,
    maximum_gap_bp, maximum_midpoint_difference_bp,
    maximum_length_fold_difference
) {
  first_length <- first_end - first_start
  second_length <- second_end - second_start
  if (first_length <= 0 || second_length <= 0) return(FALSE)
  separation <- max(second_start - first_end, first_start - second_end, 0)
  if (separation > maximum_gap_bp) return(FALSE)
  midpoint_difference <- abs(
    (first_start + first_end) / 2 - (second_start + second_end) / 2
  )
  if (midpoint_difference > maximum_midpoint_difference_bp) return(FALSE)
  max(first_length, second_length) / min(first_length, second_length) <=
    maximum_length_fold_difference
}

assign_groups <- function(
    start0, end0, maximum_gap_bp,
    maximum_midpoint_difference_bp, maximum_length_fold_difference
) {
  count <- length(start0)
  if (count < 2) return(seq_len(count))
  parent <- seq_len(count)
  find_root <- function(index) {
    while (parent[index] != index) {
      parent[index] <<- parent[parent[index]]
      index <- parent[index]
    }
    index
  }
  union_sets <- function(first, second) {
    first_root <- find_root(first)
    second_root <- find_root(second)
    if (first_root != second_root) parent[second_root] <<- first_root
  }
  ordering <- order(start0, end0)
  for (first_position in seq_len(count - 1)) {
    first <- ordering[first_position]
    for (second_position in (first_position + 1):count) {
      second <- ordering[second_position]
      if (start0[second] > end0[first] + maximum_gap_bp) break
      if (same_locus(
        start0[first], end0[first], start0[second], end0[second],
        maximum_gap_bp, maximum_midpoint_difference_bp,
        maximum_length_fold_difference
      )) {
        union_sets(first, second)
      }
    }
  }
  roots <- vapply(seq_len(count), find_root, integer(1))
  match(roots, unique(roots))
}

keep_best_per_key <- function(frame, key_columns) {
  if (!nrow(frame)) return(frame)
  keys <- apply(frame[, key_columns, drop = FALSE], 1, paste, collapse = "||")
  groups <- split(seq_len(nrow(frame)), keys)
  keep <- vapply(groups, function(indices) {
    part <- frame[indices, , drop = FALSE]
    ranking <- order(
      -part$fraction_samples_in_clade,
      -part$n_supporting_samples,
      -part$resolved_length_bp,
      part$resolved_start0,
      part$resolved_end0
    )
    indices[ranking[1]]
  }, integer(1))
  frame[sort(keep), , drop = FALSE]
}

prune_shared_boundaries <- function(frame) {
  repeat {
    previous_count <- nrow(frame)
    frame <- keep_best_per_key(
      frame, c("clade", "contig", "resolved_start0")
    )
    frame <- keep_best_per_key(
      frame, c("clade", "contig", "resolved_end0")
    )
    if (nrow(frame) == previous_count) break
  }
  rownames(frame) <- NULL
  frame
}

arguments <- parse_args(commandArgs(trailingOnly = TRUE))
metadata_path <- required_value(arguments, "metadata")
tract_path <- required_value(arguments, "tracts")
fai_path <- required_value(arguments, "fai")
output_path <- required_value(arguments, "out")

maximum_gap_bp <- as.integer(ifelse(
  is.null(arguments$maximum_gap_bp), "10000", arguments$maximum_gap_bp
))
maximum_midpoint_difference_bp <- as.integer(ifelse(
  is.null(arguments$maximum_midpoint_difference_bp),
  "150000",
  arguments$maximum_midpoint_difference_bp
))
maximum_length_fold_difference <- as.numeric(ifelse(
  is.null(arguments$maximum_length_fold_difference),
  "3",
  arguments$maximum_length_fold_difference
))
reporting_minimum_fraction <- as.numeric(ifelse(
  is.null(arguments$reporting_minimum_fraction),
  "0.30",
  arguments$reporting_minimum_fraction
))
reporting_minimum_length_bp <- as.integer(ifelse(
  is.null(arguments$reporting_minimum_length_bp),
  "100000",
  arguments$reporting_minimum_length_bp
))

metadata <- read.delim(metadata_path, check.names = FALSE)
tracts <- read.delim(tract_path, check.names = FALSE)
required_metadata <- c("sample", "Clade")
required_tracts <- c(
  "sample", "contig", "start0", "end0", "length_bp",
  "candidate_synapomorphic"
)
if (!all(required_metadata %in% names(metadata))) {
  stop("Metadata must contain sample and Clade columns.")
}
if (!all(required_tracts %in% names(tracts))) {
  stop("Tract table is missing required columns.")
}

metadata$sample <- trimws(as.character(metadata$sample))
metadata$Clade <- trimws(as.character(metadata$Clade))
metadata <- metadata[
  !is.na(metadata$Clade) & nzchar(metadata$Clade),
  c("sample", "Clade"),
  drop = FALSE
]
names(metadata)[2] <- "clade"
metadata <- unique(metadata)

fai <- read.delim(
  fai_path, header = FALSE, sep = "\t", quote = "",
  col.names = c("contig", "length", "offset", "line_bases", "line_width")
)
nuclear_fai <- fai[!grepl("(?i)chrM(?:_|$)", fai$contig, perl = TRUE), ]
nuclear_total_bp <- sum(as.numeric(nuclear_fai$length))
tracts <- tracts[tracts$contig %in% nuclear_fai$contig, , drop = FALSE]
tracts$candidate_synapomorphic <- as_flag(tracts$candidate_synapomorphic)
contig_lengths <- setNames(
  as.numeric(nuclear_fai$length),
  nuclear_fai$contig
)
burden_tracts <- tracts
burden_tracts$end0 <- pmin(
  burden_tracts$end0,
  contig_lengths[burden_tracts$contig]
)
burden_tracts <- burden_tracts[
  burden_tracts$end0 > burden_tracts$start0,
  ,
  drop = FALSE
]

sample_burden <- lapply(seq_len(nrow(metadata)), function(index) {
  sample_id <- metadata$sample[index]
  sample_tracts <- burden_tracts[
    burden_tracts$sample == sample_id,
    ,
    drop = FALSE
  ]
  data.frame(
    sample = sample_id,
    clade = metadata$clade[index],
    fraction_nuclear_genome_any_loh =
      covered_bp_by_sample(sample_tracts) / nuclear_total_bp
  )
})
sample_burden <- do.call(rbind, sample_burden)
clade_burden <- aggregate(
  fraction_nuclear_genome_any_loh ~ clade,
  data = sample_burden,
  FUN = mean
)
names(clade_burden)[2] <- "mean_fraction_nuclear_genome_any_loh"

samples_per_clade <- aggregate(
  sample ~ clade, data = metadata,
  FUN = function(values) length(unique(values))
)
names(samples_per_clade)[2] <- "n_samples_in_clade"

annotated <- merge(tracts, metadata, by = "sample", sort = FALSE)
candidate <- annotated[
  annotated$candidate_synapomorphic,
  ,
  drop = FALSE
]
if (!nrow(candidate)) stop("No candidate synapomorphic LOH tracts remain.")

resolved_parts <- lapply(
  split(candidate, list(candidate$clade, candidate$contig), drop = TRUE),
  function(part) {
    part$group <- assign_groups(
      part$start0,
      part$end0,
      maximum_gap_bp,
      maximum_midpoint_difference_bp,
      maximum_length_fold_difference
    )
    do.call(rbind, lapply(split(part, part$group), function(group) {
      resolved_start <- round(quantile(
        group$start0, 0.25, names = FALSE, type = 7
      ))
      resolved_end <- round(quantile(
        group$end0, 0.75, names = FALSE, type = 7
      ))
      data.frame(
        clade = group$clade[1],
        contig = group$contig[1],
        resolved_start0 = resolved_start,
        resolved_end0 = resolved_end,
        resolved_length_bp = resolved_end - resolved_start,
        observed_start0_min = min(group$start0),
        observed_start0_max = max(group$start0),
        observed_end0_min = min(group$end0),
        observed_end0_max = max(group$end0),
        n_input_tracts = nrow(group),
        n_supporting_samples = length(unique(group$sample)),
        median_tract_purity = median(group$purity, na.rm = TRUE),
        median_heterozygous_snps_per_window =
          median(group$median_N_het, na.rm = TRUE)
      )
    }))
  }
)
resolved <- do.call(rbind, resolved_parts)
rownames(resolved) <- NULL
resolved <- merge(resolved, samples_per_clade, by = "clade", sort = FALSE)
resolved$fraction_samples_in_clade <-
  resolved$n_supporting_samples / resolved$n_samples_in_clade
resolved <- merge(resolved, clade_burden, by = "clade", sort = FALSE)
resolved$fraction_nuclear_genome_resolved_region <-
  (
    pmin(
      resolved$resolved_end0,
      contig_lengths[resolved$contig]
    ) - resolved$resolved_start0
  ) / nuclear_total_bp
resolved$resolved_effective_length_bp <- pmax(
  0,
  pmin(
    resolved$resolved_end0,
    contig_lengths[resolved$contig]
  ) - resolved$resolved_start0
)
resolved <- prune_shared_boundaries(resolved)
resolved$passes_reporting_filter <-
  resolved$fraction_samples_in_clade > reporting_minimum_fraction &
  resolved$resolved_length_bp > reporting_minimum_length_bp
resolved$maximum_gap_bp <- maximum_gap_bp
resolved$maximum_midpoint_difference_bp <- maximum_midpoint_difference_bp
resolved$maximum_length_fold_difference <- maximum_length_fold_difference
resolved$reporting_minimum_fraction_strict <- reporting_minimum_fraction
resolved$reporting_minimum_length_bp_strict <- reporting_minimum_length_bp
resolved <- resolved[
  order(
    resolved$clade,
    resolved$contig,
    resolved$resolved_start0,
    resolved$resolved_end0
  ),
]
write.table(
  resolved,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
