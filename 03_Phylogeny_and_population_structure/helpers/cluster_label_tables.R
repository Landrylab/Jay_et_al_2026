# Cluster identity defines the aggregation unit; clades supply display labels.
normalize_label <- function(values) {
  values <- trimws(as.character(values))
  values[is.na(values) | values %in% c("", "NA", "NaN", "nan", ".")] <- NA_character_
  values
}

cluster_display_labels <- function(clusters, clades) {
  clusters <- as.character(clusters)
  levels <- unique(clusters)
  numeric_levels <- suppressWarnings(as.numeric(levels))
  levels <- if (all(is.finite(numeric_levels))) levels[order(numeric_levels)] else sort(levels)
  clades <- normalize_label(clades)
  clades[clades == "Other" & !is.na(clades)] <- NA_character_
  assigned <- vapply(levels, function(cluster) {
    values <- unique(clades[clusters == cluster & !is.na(clades)])
    if (length(values) == 1L) values else NA_character_
  }, character(1))
  repeated <- duplicated(assigned) | duplicated(assigned, fromLast = TRUE)
  display <- ifelse(is.na(assigned), paste("PCA cluster", levels), paste("Clade", assigned))
  display[repeated & !is.na(assigned)] <- paste0(display[repeated & !is.na(assigned)], " (C", levels[repeated & !is.na(assigned)], ")")
  data.frame(cluster = levels, extrapolated_clade = unname(assigned), display_label = unname(display), stringsAsFactors = FALSE)
}

cluster_metadata_table <- function(clusters, labels, clades, label_name) {
  clusters <- as.character(clusters)
  labels <- normalize_label(labels)
  mapping <- cluster_display_labels(clusters, clades)
  categories <- sort(unique(labels[!is.na(labels)]))
  # A missing-category row keeps every cluster visible when no label is available.
  if (!length(categories)) categories <- NA_character_
  result <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(index) {
    selected <- clusters == mapping$cluster[index]
    total <- sum(selected & !is.na(labels))
    counts <- vapply(categories, function(category) {
      if (is.na(category)) 0L else sum(selected & !is.na(labels) & labels == category, na.rm = TRUE)
    }, integer(1))
    data.frame(
      label = label_name, cluster = mapping$cluster[index],
      cluster_label = paste0("C", mapping$cluster[index]),
      extrapolated_clade = mapping$extrapolated_clade[index],
      display_label = mapping$display_label[index], external_category = categories,
      count = unname(counts), cluster_total_complete = total,
      within_cluster_fraction = if (total > 0L) unname(counts) / total else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}
