#!/usr/bin/env Rscript
# Plot the evidence used to review candidate baseline ploidy and segment calls.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing argument: ", flag)
  args[index + 1L]
}
depth <- read.delim(get_arg("--depth"), check.names = FALSE)
baf <- read.delim(get_arg("--baf"), check.names = FALSE)
prefix <- get_arg("--prefix")
output <- get_arg("--output")
segments <- read.delim(paste0(prefix, ".segments.tsv"))
ploidy <- read.delim(paste0(prefix, ".ploidy.tsv"))
baseline <- ploidy$baseline_ploidy[1]
depth$log2_ratio <- log2(depth$reads_overlapping_end_corrected / median(depth$reads_overlapping_end_corrected))
baf$allele_fraction <- baf$alt_count / (baf$ref_count + baf$alt_count)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
pdf(output, width = 11, height = 9, onefile = TRUE)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (candidate in c(2, 3, 4, 6)) {
  if (nrow(baf)) {
    hist(baf$allele_fraction, breaks = seq(0, 1, by = 0.01), main = paste(candidate, "N expected allele fractions"),
      xlab = "Alternate-allele fraction", col = "grey80", border = NA, xlim = c(0, 1))
  } else {
    plot.new(); title(main = paste(candidate, "N: no informative SNPs"))
    plot.window(xlim = c(0, 1), ylim = c(0, 1)); axis(1)
  }
  abline(v = seq_len(candidate - 1) / candidate, col = "#9B1B30", lwd = 2)
}
for (chromosome in unique(depth$chr)) {
  d <- depth[depth$chr == chromosome, ]
  b <- baf[baf$chr == chromosome, ]
  s <- segments[segments$chr == chromosome, ]
  limits <- range(c(d$start, d$end)) / 1e6
  par(mfrow = c(3, 1), mar = c(4.5, 4, 2.5, 1))
  plot((d$start + d$end) / 2e6, d$log2_ratio, pch = 16, cex = 0.35,
    xlim = limits, xlab = "Position (Mb)", ylab = "Log2 depth ratio",
    main = paste(chromosome, "| baseline", baseline, "N |", ploidy$status[1]))
  abline(h = 0, col = "grey60")
  segments(s$start / 1e6, s$median_log2_ratio, s$end / 1e6, s$median_log2_ratio, col = "#225EA8", lwd = 2)
  plot(b$pos / 1e6, b$allele_fraction, pch = 16, cex = 0.35, xlim = limits,
    ylim = c(0, 1), xlab = "Position (Mb)", ylab = "Alternate-allele fraction")
  abline(h = seq_len(baseline - 1) / baseline, col = "grey65", lty = 2)
  abline(v = s$start / 1e6, col = "#225EA850", lty = 3)
  plot(NA, xlim = limits, ylim = c(0, max(8, s$integer_copy_number, na.rm = TRUE)),
    xlab = "Position (Mb)", ylab = "Segment copy number")
  segments(s$start / 1e6, s$integer_copy_number, s$end / 1e6, s$integer_copy_number, lwd = 3, col = "#225EA8")
  abline(h = baseline, col = "grey60", lty = 2)
  abline(v = s$start / 1e6, col = "#225EA850", lty = 3)
}
dev.off()
