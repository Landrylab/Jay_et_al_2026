#!/usr/bin/env Rscript
# =============================================================================
# Description
#   Helper logic for the GCTA MLMA-LOCO GWAS reproducibility bundle.
# Inputs
#   Local project metadata, PLINK files, MLMA outputs, SnpEff annotations, and GWAS tables.
# Outputs
#   Phenotype files, GWAS summaries, annotated hits, figures, or publication tables.
# =============================================================================

# =============================================================================
# Configuration
# =============================================================================


args <- commandArgs(trailingOnly = TRUE)


# =============================================================================
# Functions
# =============================================================================

arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

mlma <- arg_value("--mlma")
outdir <- arg_value("--outdir")
trait <- arg_value("--trait", "trait")
mode <- arg_value("--mode", "mode")
grey_snps_path <- arg_value("--grey-snps")
grey_label <- arg_value("--grey-label", "In LD with FungAMR mutation")
label_hits_path <- arg_value("--label-hits")
monochrome_points <- "--monochrome-points" %in% args
manhattan_ymax_arg <- arg_value("--manhattan-ymax")
manhattan_ymax <- if (is.null(manhattan_ymax_arg)) NA_real_ else suppressWarnings(as.numeric(manhattan_ymax_arg))
if (!is.null(manhattan_ymax_arg) && (!is.finite(manhattan_ymax) || manhattan_ymax <= 0)) {
  stop("--manhattan-ymax must be a positive numeric value")
}

if (is.null(mlma) || is.null(outdir)) {
  stop("Usage: plot_mlma_results.R --mlma <gwas.mlma> --outdir <plotdir> --trait <drug> --mode <mic|rs_binary>")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_mlma <- function(path) {
  read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

pick_col <- function(d, candidates, label) {
  hit <- candidates[candidates %in% names(d)]
  if (!length(hit)) {
    stop(sprintf("No %s column found. Tried: %s", label, paste(candidates, collapse = ", ")))
  }
  hit[[1]]
}

chr_to_numeric <- function(x) {
  text <- as.character(x)
  text <- sub("^Chr", "", text, ignore.case = TRUE)
  text <- sub("^Ca22chr", "", text, ignore.case = TRUE)
  text <- sub("A_C_albicans_SC5314.*$", "", text)
  text[text %in% c("M", "MT")] <- "8"
  text[text %in% c("R", "RA")] <- "9"
  suppressWarnings(as.numeric(text))
}

d <- read_mlma(mlma)
chr_col <- pick_col(d, c("Chr", "CHR", "chr", "Chromosome", "chromosome"), "chromosome")
bp_col <- pick_col(d, c("bp", "BP", "pos", "POS", "Position", "position"), "position")
p_col <- pick_col(d, c("p", "P", "pval", "PVAL", "p.value", "P.value"), "P value")
snp_col <- if ("SNP" %in% names(d)) "SNP" else if ("snp" %in% names(d)) "snp" else NA

d[[p_col]] <- suppressWarnings(as.numeric(d[[p_col]]))
d[[bp_col]] <- suppressWarnings(as.numeric(d[[bp_col]]))
d$chr_numeric <- chr_to_numeric(d[[chr_col]])

d <- d[
  is.finite(d[[p_col]]) &
    d[[p_col]] > 0 &
    d[[p_col]] <= 1 &
    is.finite(d[[bp_col]]) &
    is.finite(d$chr_numeric),
  ,
  drop = FALSE
]

if (!nrow(d)) {
  stop("No valid association rows after filtering MLMA results.")
}

grey_snps <- character()
if (!is.null(grey_snps_path) && nzchar(grey_snps_path) && file.exists(grey_snps_path)) {
  grey_snps <- unique(readLines(grey_snps_path, warn = FALSE))
  grey_snps <- grey_snps[nzchar(grey_snps)]
}

if (!is.na(snp_col) && length(grey_snps)) {
  d$fungamr_ld_grey <- d[[snp_col]] %in% grey_snps
} else {
  d$fungamr_ld_grey <- FALSE
}
n_fungamr_ld_grey <- sum(d$fungamr_ld_grey, na.rm = TRUE)

label_hits <- data.frame()
if (!is.null(label_hits_path) && nzchar(label_hits_path) && file.exists(label_hits_path) && !is.na(snp_col)) {
  label_input <- read.table(label_hits_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("snp", "label") %in% names(label_input))) {
    stop("Label hits file must contain columns: snp, label")
  }
  if (!"gene_role" %in% names(label_input)) {
    label_input$gene_role <- "candidate"
  }
  label_hits <- label_input[, c("snp", "label", "gene_role"), drop = FALSE]
  label_hits <- label_hits[!duplicated(label_hits$label), , drop = FALSE]
}

desired_chr_order <- c(1:7, 9, 8)
present_chr <- unique(d$chr_numeric)
chr_levels <- desired_chr_order[desired_chr_order %in% present_chr]
extra_chr <- sort(setdiff(present_chr, desired_chr_order))
chr_levels <- c(chr_levels, extra_chr)
d <- d[order(match(d$chr_numeric, chr_levels), d[[bp_col]]), , drop = FALSE]

label_chr <- function(chr) {
  if (chr == 8) {
    return("ChrM")
  }
  if (chr == 9) {
    return("ChrR")
  }
  paste0("Chr", chr)
}
chr_labels <- vapply(chr_levels, label_chr, character(1))

gap_size <- 5e5
offset <- 0
centers <- numeric(length(chr_levels))
plot_x <- numeric(nrow(d))

for (i in seq_along(chr_levels)) {
  chr <- chr_levels[[i]]
  idx <- which(d$chr_numeric == chr)
  this_bp <- d[[bp_col]][idx]
  plot_x[idx] <- this_bp + offset
  centers[[i]] <- (min(plot_x[idx]) + max(plot_x[idx])) / 2
  offset <- max(plot_x[idx]) + gap_size
}

d$plot_x <- plot_x
d$raw_logp <- -log10(d[[p_col]])
d$logp <- d$raw_logp
if (is.finite(manhattan_ymax)) {
  d$logp <- pmin(d$raw_logp, manhattan_ymax * 0.99)
}

if (nrow(label_hits)) {
  label_hits <- merge(
    label_hits,
    d[, c(snp_col, "plot_x", "logp", "fungamr_ld_grey"), drop = FALSE],
    by.x = "snp",
    by.y = snp_col,
    all = FALSE,
    sort = FALSE
  )
  label_hits <- label_hits[order(label_hits$plot_x, -label_hits$logp), , drop = FALSE]
}

bonf_p <- 0.05 / nrow(d)
bonf_logp <- -log10(bonf_p)
lambda_gc <- median(qchisq(1 - d[[p_col]], df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

plot_prefix <- paste0(trait, ".", mode)
out_manhattan_pdf <- file.path(outdir, paste0("manhattan.", plot_prefix, ".pdf"))
out_manhattan_png <- file.path(outdir, paste0("manhattan.", plot_prefix, ".png"))
out_qq_pdf <- file.path(outdir, paste0("qq.", plot_prefix, ".pdf"))
out_qq_png <- file.path(outdir, paste0("qq.", plot_prefix, ".png"))
out_summary <- file.path(outdir, paste0("plot_summary.", plot_prefix, ".txt"))
out_top <- file.path(outdir, paste0("top_hits.", plot_prefix, ".tsv"))
plot_title <- paste0(toupper(substr(trait, 1, 1)), substr(trait, 2, nchar(trait)))

palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#00A6D6", "#1f78b4", "#8dd3c7"
)
point_col <- palette[((match(d$chr_numeric, chr_levels) - 1) %% length(palette)) + 1]
if (monochrome_points) {
  point_col <- rep("black", nrow(d))
}
point_col[d$fungamr_ld_grey] <- "#BDBDBD"

label_cex <- if (nrow(label_hits) >= 30) 0.85 else if (nrow(label_hits) >= 20) 0.95 else 1.20
canonical_label_color <- "#C0392B"
candidate_label_color <- "#E69F00"

draw_candidate_labels <- function() {
  if (!nrow(label_hits)) {
    return(invisible(NULL))
  }

  labs <- label_hits
  many_candidate_labels <- nrow(labs[labs$gene_role != "canonical", , drop = FALSE]) >= 20
  labs$gene_role <- tolower(labs$gene_role)
  labs$gene_role[!labs$gene_role %in% c("canonical", "candidate")] <- "candidate"
  usr <- par("usr")
  y_span <- usr[4] - usr[3]
  x_span <- usr[2] - usr[1]
  pad_x <- 0.008 * x_span
  pad_y <- 0.012 * y_span
  label_floor <- max(
    bonf_logp + 0.035 * y_span,
    unname(quantile(d$logp, probs = 0.995, na.rm = TRUE)) + 0.025 * y_span
  )
  label_floor <- min(label_floor, usr[4] - 0.18 * y_span)
  local_label_y <- function(x, w, h, anchor_y) {
    span <- w / 2 + 1.5 * pad_x
    nearby <- abs(d$plot_x - x) <= span
    local_max <- if (any(nearby)) max(d$logp[nearby], na.rm = TRUE) else 0
    max(
      anchor_y + 0.035 * y_span,
      label_floor,
      local_max + h / 2 + 0.020 * y_span
    )
  }

  labs$width <- strwidth(labs$label, units = "user", cex = label_cex, font = 2)
  labs$height <- strheight(labs$label, units = "user", cex = label_cex, font = 2)
  clamp_label_x <- function(x, w) {
    max(usr[1] + w / 2 + pad_x, min(usr[2] - w / 2 - pad_x, x))
  }

  placed <- list()
  overlaps_any <- function(x, y, w, h) {
    if (!length(placed)) {
      return(FALSE)
    }
    for (box in placed) {
      if (
        x - w / 2 - pad_x < box$x + box$w / 2 &&
          x + w / 2 + pad_x > box$x - box$w / 2 &&
          y - h / 2 - pad_y < box$y + box$h / 2 &&
          y + h / 2 + pad_y > box$y - box$h / 2
      ) {
        return(TRUE)
      }
    }
    FALSE
  }

  remember_boxes <- function(x) {
    if (!nrow(x)) {
      return(invisible(NULL))
    }
    for (i in seq_len(nrow(x))) {
      placed[[length(placed) + 1]] <<- list(
        x = x$label_x[i],
        y = x$label_y[i],
        w = x$width[i],
        h = x$height[i]
      )
    }
    invisible(NULL)
  }

  resolve_candidate_boxes <- function(x) {
    if (!nrow(x)) {
      return(x)
    }
    x <- x[order(x$label_y, x$label_x, x$label), , drop = FALSE]
    for (i in seq_len(nrow(x))) {
      y <- x$label_y[i]
      tries <- 0
      while (overlaps_any(x$label_x[i], y, x$width[i], x$height[i]) && tries < 100) {
        y <- y + 0.050 * y_span
        if (y + x$height[i] / 2 > usr[4]) {
          y <- usr[4] - x$height[i] / 2 - pad_y
          break
        }
        tries <- tries + 1
      }
      x$label_y[i] <- y
      placed[[length(placed) + 1]] <<- list(
        x = x$label_x[i],
        y = x$label_y[i],
        w = x$width[i],
        h = x$height[i]
      )
    }
    x
  }

  spread_ordered_labels <- function(x, gap_multiplier = 1.8) {
    if (!nrow(x)) {
      return(x)
    }
    x <- x[order(x$plot_x, x$label), , drop = FALSE]
    if (nrow(x) == 1) {
      x$label_x <- clamp_label_x(x$plot_x, x$width)
      x$label_y <- mapply(local_label_y, x$label_x, x$width, x$height, x$logp)
      x$label_y <- pmin(x$label_y, usr[4] - 0.08 * y_span)
      return(x)
    }

    label_x <- x$plot_x
    for (k in 2:nrow(x)) {
      min_gap <- (x$width[[k - 1]] + x$width[[k]]) / 2 + gap_multiplier * pad_x
      label_x[[k]] <- max(label_x[[k]], label_x[[k - 1]] + min_gap)
    }

    label_x <- label_x + mean(x$plot_x - label_x)
    for (iter in seq_len(4)) {
      left_shift <- max((usr[1] + x$width / 2 + pad_x) - label_x)
      if (is.finite(left_shift) && left_shift > 0) {
        label_x <- label_x + left_shift
      }
      right_shift <- max(label_x - (usr[2] - x$width / 2 - pad_x))
      if (is.finite(right_shift) && right_shift > 0) {
        label_x <- label_x - right_shift
      }
      for (k in 2:nrow(x)) {
        min_gap <- (x$width[[k - 1]] + x$width[[k]]) / 2 + gap_multiplier * pad_x
        label_x[[k]] <- max(label_x[[k]], label_x[[k - 1]] + min_gap)
      }
    }

    x$label_x <- label_x
    cluster_floor <- max(mapply(local_label_y, x$label_x, x$width, x$height, x$logp))
    if (many_candidate_labels && nrow(x) >= 3) {
      n_levels <- min(6, max(3, ceiling(nrow(x) / 5)))
      level <- (seq_len(nrow(x)) - 1) %% n_levels
      x$label_y <- cluster_floor + level * (x$height + 2.5 * pad_y)
      x$label_y <- pmin(x$label_y, usr[4] - 0.06 * y_span)
    } else {
      x$label_y <- min(cluster_floor, usr[4] - 0.08 * y_span)
    }
    x
  }

  layout_candidates <- function(x) {
    if (!nrow(x)) {
      return(x)
    }
    x <- x[order(x$plot_x, x$label), , drop = FALSE]
    cluster_gap <- 0.045 * x_span
    cluster_id <- cumsum(c(TRUE, diff(x$plot_x) > cluster_gap))
    out <- list()
    for (cluster in unique(cluster_id)) {
      out[[length(out) + 1]] <- spread_ordered_labels(x[cluster_id == cluster, , drop = FALSE])
    }
    x <- do.call(rbind, out)
    resolve_candidate_boxes(x)
  }

  layout_canonical <- function(x) {
    if (!nrow(x)) {
      return(x)
    }
    x <- x[order(x$plot_x, -x$logp, x$label), , drop = FALSE]
    x$label_x <- mapply(clamp_label_x, x$plot_x, x$width)
    x$label_y <- mapply(local_label_y, x$label_x, x$width, x$height, x$logp)
    x$label_y <- pmin(x$label_y, usr[4] - 0.08 * y_span)
    for (i in seq_len(nrow(x))) {
      y <- x$label_y[i]
      tries <- 0
      while (overlaps_any(x$label_x[i], y, x$width[i], x$height[i]) && tries < 80) {
        y <- y + 0.040 * y_span
        if (y + x$height[i] / 2 > usr[4]) {
          y <- usr[4] - x$height[i] / 2 - pad_y
          break
        }
        tries <- tries + 1
      }
      x$label_y[i] <- y
      placed[[length(placed) + 1]] <<- list(
        x = x$label_x[i],
        y = x$label_y[i],
        w = x$width[i],
        h = x$height[i]
      )
    }
    x
  }

  candidate_labs <- layout_candidates(labs[labs$gene_role != "canonical", , drop = FALSE])
  canonical_labs <- layout_canonical(labs[labs$gene_role == "canonical", , drop = FALSE])
  labs <- rbind(candidate_labs, canonical_labs)
  labs <- labs[order(labs$plot_x, labs$label_y), , drop = FALSE]
  labs$label_color <- ifelse(labs$gene_role == "canonical", canonical_label_color, candidate_label_color)

  idx <- match(labs$snp, d[[snp_col]])
  idx <- idx[!is.na(idx)]
  if (length(idx)) {
    points(
      d$plot_x[idx],
      d$logp[idx],
      pch = 21,
      bg = point_col[idx],
      col = labs$label_color[match(d[[snp_col]][idx], labs$snp)],
      lwd = 1.8,
      cex = 1.25
    )
  }

  dx <- labs$label_x - labs$plot_x
  dy <- labs$label_y - labs$logp
  half_w <- labs$width / 2 + 0.004 * x_span
  half_h <- labs$height / 2 + 0.004 * y_span
  scale_x <- ifelse(abs(dx) > 0, half_w / abs(dx), Inf)
  scale_y <- ifelse(abs(dy) > 0, half_h / abs(dy), Inf)
  scale <- pmin(scale_x, scale_y)
  scale[!is.finite(scale)] <- 0
  line_end_x <- labs$label_x - dx * scale
  line_end_y <- labs$label_y - dy * scale

  segments(
    x0 = labs$plot_x,
    y0 = labs$logp + 0.08,
    x1 = line_end_x,
    y1 = line_end_y,
    col = labs$label_color,
    lwd = 1
  )
  text(
    labs$label_x,
    labs$label_y,
    labels = labs$label,
    col = "white",
    cex = label_cex + 0.10,
    font = 2
  )
  halo_dx <- 0.0012 * x_span
  halo_dy <- 0.0012 * y_span
  halo_offsets <- expand.grid(dx = c(-halo_dx, 0, halo_dx), dy = c(-halo_dy, 0, halo_dy))
  halo_offsets <- halo_offsets[halo_offsets$dx != 0 | halo_offsets$dy != 0, , drop = FALSE]
  for (j in seq_len(nrow(halo_offsets))) {
    text(
      labs$label_x + halo_offsets$dx[j],
      labs$label_y + halo_offsets$dy[j],
      labels = labs$label,
      col = "white",
      cex = label_cex,
      font = 2
    )
  }
  text(
    labs$label_x,
    labs$label_y,
    labels = labs$label,
    col = labs$label_color,
    cex = label_cex,
    font = 2
  )
}

draw_manhattan <- function() {
  par(mar = c(5.8, 6, 3.2, 1), mgp = c(3.2, 1, 0))
  ylim_extra <- if (nrow(label_hits)) max(2, 0.55 * max(d$logp, na.rm = TRUE)) else 0
  dynamic_ymax <- max(d$logp, na.rm = TRUE) + ylim_extra
  plot_ymax <- if (is.finite(manhattan_ymax)) manhattan_ymax else dynamic_ymax
  plot(
    d$plot_x,
    d$logp,
    pch = 16,
    cex = 0.90,
    col = point_col,
    type = "n",
    xaxt = "n",
    xlab = "",
    ylab = "",
    main = plot_title,
    bty = "l",
    ylim = c(0, plot_ymax),
    cex.axis = 1.15,
    cex.lab = 1.35,
    cex.main = 1.65,
    font.axis = 2,
    font.lab = 2,
    font.main = 2
  )
  background_idx <- which(!d$fungamr_ld_grey)
  if (length(background_idx)) {
    points(
      d$plot_x[background_idx],
      d$logp[background_idx],
      pch = 16,
      cex = 0.90,
      col = point_col[background_idx]
    )
  }
  abline(h = bonf_logp, lty = 3, lwd = 1.5, col = "#C0392B")
  axis(1, at = centers, labels = chr_labels, tick = FALSE, line = 0.5, cex.axis = 1.15, font = 2)
  mtext("Chromosome", side = 1, line = 3.1, cex = 1.35, font = 2)
  mtext(expression(bold(-log[10](italic(P)))), side = 2, line = 3.7, cex = 1.55)
  if (n_fungamr_ld_grey > 0) {
    points(
      d$plot_x[d$fungamr_ld_grey],
      d$logp[d$fungamr_ld_grey],
      pch = 16,
      cex = 0.90,
      col = "#BDBDBD"
    )
  }
  draw_candidate_labels()
  legend_labels <- sprintf("Bonferroni threshold=%.3g", bonf_p)
  legend_lty <- 3
  legend_lwd <- 1.5
  legend_pch <- NA
  legend_col <- "#C0392B"
  if (n_fungamr_ld_grey > 0) {
    legend_labels <- c(legend_labels, grey_label)
    legend_lty <- c(legend_lty, NA)
    legend_lwd <- c(legend_lwd, NA)
    legend_pch <- c(legend_pch, 16)
    legend_col <- c(legend_col, "#BDBDBD")
  }
  legend(
    "topright",
    legend = legend_labels,
    lty = legend_lty,
    lwd = legend_lwd,
    pch = legend_pch,
    col = legend_col,
    bty = "n",
    cex = 1.05,
    text.font = 2
  )
}

draw_qq <- function() {
  obs_p <- sort(d[[p_col]])
  exp_p <- ppoints(length(obs_p))
  obs <- -log10(obs_p)
  exp <- -log10(exp_p)
  max_xy <- max(c(obs, exp), na.rm = TRUE)

  par(mar = c(5.8, 6, 3.2, 1), mgp = c(3.2, 1, 0))
  plot(
    exp,
    obs,
    pch = 16,
    cex = 1.10,
    xlim = c(0, max_xy),
    ylim = c(0, max_xy),
    xlab = "",
    ylab = "",
    main = plot_title,
    bty = "l",
    cex.axis = 1.15,
    cex.lab = 1.35,
    cex.main = 1.65,
    font.axis = 2,
    font.lab = 2,
    font.main = 2
  )
  abline(0, 1, lwd = 1.2)
  mtext(expression(bold(Expected~~-log[10](italic(P)))), side = 1, line = 3.1, cex = 1.35)
  mtext(expression(bold(Observed~~-log[10](italic(P)))), side = 2, line = 3.7, cex = 1.35)
  usr <- par("usr")
  text(
    x = usr[1] + 0.08 * (usr[2] - usr[1]),
    y = usr[4] - 0.08 * (usr[4] - usr[3]),
    labels = sprintf("lambda = %.3f", lambda_gc),
    adj = c(0, 1),
    cex = 1.2,
    font = 2
  )
}


# =============================================================================
# Main block
# =============================================================================

pdf(out_manhattan_pdf, width = 13, height = 6.5, useDingbats = FALSE)
draw_manhattan()
dev.off()

png(out_manhattan_png, width = 2600, height = 1300, res = 200)
draw_manhattan()
dev.off()

pdf(out_qq_pdf, width = 6.5, height = 6.5, useDingbats = FALSE)
draw_qq()
dev.off()

png(out_qq_png, width = 1400, height = 1400, res = 200)
draw_qq()
dev.off()

summary_lines <- c(
  sprintf("Trait\t%s", trait),
  sprintf("Mode\t%s", mode),
  sprintf("Input\t%s", mlma),
  sprintf("N_variants\t%d", nrow(d)),
  sprintf("Bonferroni_P\t%.15g", bonf_p),
  sprintf("Bonferroni_log10P\t%.6f", bonf_logp),
  sprintf("Lambda_GC\t%.6f", lambda_gc),
  sprintf("Manhattan_Ymax\t%s", ifelse(is.finite(manhattan_ymax), as.character(manhattan_ymax), "")),
  sprintf("FungAMR_LD_grey_snps_file\t%s", ifelse(is.null(grey_snps_path), "", grey_snps_path)),
  sprintf("N_FungAMR_LD_grey_points\t%d", n_fungamr_ld_grey),
  sprintf("Candidate_label_hits_file\t%s", ifelse(is.null(label_hits_path), "", label_hits_path)),
  sprintf("N_candidate_label_hits\t%d", nrow(label_hits))
)
writeLines(summary_lines, out_summary)

top_cols <- c(chr_col, snp_col, bp_col, p_col)
top_cols <- top_cols[!is.na(top_cols)]
top <- d[order(d[[p_col]]), top_cols, drop = FALSE]
write.table(head(top, 100), out_top, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("Wrote %s\n", out_manhattan_pdf))
cat(sprintf("Wrote %s\n", out_manhattan_png))
cat(sprintf("Wrote %s\n", out_qq_pdf))
cat(sprintf("Wrote %s\n", out_qq_png))
cat(sprintf("Wrote %s\n", out_summary))
