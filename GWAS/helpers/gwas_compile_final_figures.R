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


suppressPackageStartupMessages(library(magick))

root <- normalizePath(".", mustWork = TRUE)
outdir <- file.path(root, "final_gwas_paper_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

runs <- list(
  maf0p05 = list(folder = "gcta_mlma_loco/gcta_mlma_loco_maf005_six_antifungals", label = "MAF 0.05"),
  maf0p005 = list(folder = "gcta_mlma_loco/gcta_mlma_loco_maf0005_six_antifungals", label = "MAF 0.005")
)

modes <- list(
  mic = "MIC",
  rs_binary = "R/S binary"
)

plot_types <- c("manhattan", "qq")
azole_drugs <- c("fluconazole", "posaconazole", "voriconazole")
echinocandin_drugs <- c("caspofungin", "micafungin", "anidulafungin")


# =============================================================================
# Functions
# =============================================================================

plot_path <- function(run_folder, drug, mode, plot_type) {
  file.path(
    root,
    run_folder,
    "outputs",
    drug,
    mode,
    "plots",
    sprintf("%s.%s.%s.png", plot_type, drug, mode)
  )
}

annotate_panel_letter <- function(img, letter) {
  image_annotate(
    img,
    letter,
    gravity = "northwest",
    location = "+45+35",
    size = 96,
    weight = 700,
    color = "black"
  )
}

read_plot_stack <- function(paths, panel_letter = NULL) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing source plot(s):\n", paste(missing, collapse = "\n"))
  }
  imgs <- lapply(paths, image_read)
  if (!is.null(panel_letter) && length(imgs)) {
    imgs[[1]] <- annotate_panel_letter(imgs[[1]], panel_letter)
  }
  widths <- vapply(imgs, function(img) image_info(img)$width[[1]], numeric(1))
  target_width <- max(widths)
  imgs <- Map(
    function(img, width) {
      if (width == target_width) {
        img
      } else {
        image_extent(img, geometry = sprintf("%dx", target_width), gravity = "center", color = "white")
      }
    },
    imgs,
    widths
  )
  image_append(do.call(c, imgs), stack = TRUE)
}

compose_figure <- function(run_key, run_info, mode, mode_label, plot_type) {
  left_paths <- vapply(azole_drugs, plot_path, character(1), run_folder = run_info$folder, mode = mode, plot_type = plot_type)
  right_paths <- vapply(echinocandin_drugs, plot_path, character(1), run_folder = run_info$folder, mode = mode, plot_type = plot_type)

  left_stack <- read_plot_stack(left_paths, panel_letter = "A")
  right_stack <- read_plot_stack(right_paths, panel_letter = "B")
  left_info <- image_info(left_stack)
  right_info <- image_info(right_stack)
  target_height <- max(left_info$height[[1]], right_info$height[[1]])

  left_stack <- image_extent(
    left_stack,
    geometry = sprintf("%dx%d", left_info$width[[1]], target_height),
    gravity = "north",
    color = "white"
  )
  right_stack <- image_extent(
    right_stack,
    geometry = sprintf("%dx%d", right_info$width[[1]], target_height),
    gravity = "north",
    color = "white"
  )

  left_col <- left_stack
  right_col <- right_stack

  spacer <- image_blank(width = 90, height = image_info(left_col)$height[[1]], color = "white")
  composite <- image_append(c(left_col, spacer, right_col), stack = FALSE)

  out_name <- sprintf("gwas_%s_%s_%s_panel_A_azoles_panel_B_echinocandins.png", plot_type, run_key, mode)
  out_path <- file.path(outdir, out_name)
  image_write(composite, path = out_path, format = "png")

  data.frame(
    figure = out_name,
    plot_type = plot_type,
    maf = run_info$label,
    mode = mode_label,
    output = out_path,
    panel_A_sources = paste(left_paths, collapse = ";"),
    panel_B_sources = paste(right_paths, collapse = ";"),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Main block
# =============================================================================

manifest <- list()
for (run_key in names(runs)) {
  for (mode in names(modes)) {
    for (plot_type in plot_types) {
      message(sprintf("Composing %s %s %s", plot_type, run_key, mode))
      manifest[[length(manifest) + 1]] <- compose_figure(run_key, runs[[run_key]], mode, modes[[mode]], plot_type)
    }
  }
}

manifest <- do.call(rbind, manifest)
manifest_path <- file.path(outdir, "final_gwas_paper_figure_manifest.tsv")
write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", manifest_path)
