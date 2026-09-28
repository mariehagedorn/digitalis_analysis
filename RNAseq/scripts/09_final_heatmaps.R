#!/usr/bin/env Rscript

# ==============================================================================
# Flavonoid-related expression heatmap for Digitalis purpurea RNA-seq data
# ==============================================================================
#
# Purpose
# -------
# Create a two-panel figure summarizing expression of a predefined flavonoid-
# related candidate-gene panel before (t0) and after 24 h (t24) of continuous-
# light treatment.
#
# Panel A: sample-wise expression heatmap
#   - colours: row-wise z-scores of VST-transformed gene expression
#   - numbers: gene-level TPM values
#   - stars in t24 cells: significance of the within-genotype t24-vs-t0 contrast
#
# Panel B: genotype-specific t24-vs-t0 log2 fold changes
#   - *  panel-level FDR < 0.05
#   - ** genome-wide FDR < 0.05
#
# Required R packages
# -------------------
# DESeq2, pheatmap, grid, gridExtra
#
# Expected directory structure
# ----------------------------
# <base_dir>/
#   results/
#     dds_paired_model.rds
#     full_flavonoid_time_effects_v3.csv
#     tximport_gene_level.rds
#   metadata/
#     full_flavonoid_candidate_panel_v3.csv
#   plots/
#
# The candidate-panel CSV is also searched for in R_scripts_used_final/, next to
# this script, and in the current working directory.
#
# Usage
# -----
# Rscript plot_flavonoid_expression_heatmap.R /path/to/deseq2
#
# Alternatively, define DIGITALIS_DESEQ2_DIR and run the script without an
# argument. If neither is supplied, the current working directory is used.
#
# Output
# ------
# PNG and PDF are written to <base_dir>/plots/.
# SVG is written to <base_dir>/plots/final/.
# ============================================================================== 

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(grid)
  library(gridExtra)
})

# ==============================================================================
# 1. Paths and input data
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

base_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  Sys.getenv("DIGITALIS_DESEQ2_DIR", unset = ".")
}

base_dir <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(
    normalizePath(
      sub("^--file=", "", script_arg[1]),
      winslash = "/",
      mustWork = FALSE
    )
  )
} else {
  getwd()
}

results_dir <- file.path(base_dir, "results")
metadata_dir <- file.path(base_dir, "metadata")
scripts_dir <- file.path(base_dir, "R_scripts_used_final")
plots_dir <- file.path(base_dir, "plots")
final_plots_dir <- file.path(plots_dir, "final")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_plots_dir, recursive = TRUE, showWarnings = FALSE)

find_existing_file <- function(candidate_paths, description) {
  existing <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing) == 0) {
    stop(
      description,
      " was not found. Checked paths: ",
      paste(candidate_paths, collapse = "; ")
    )
  }

  existing[1]
}

panel_file <- find_existing_file(
  c(
    file.path(metadata_dir, "full_flavonoid_candidate_panel_v3.csv"),
    file.path(scripts_dir, "full_flavonoid_candidate_panel_v3.csv"),
    file.path(script_dir, "full_flavonoid_candidate_panel_v3.csv"),
    file.path(getwd(), "full_flavonoid_candidate_panel_v3.csv")
  ),
  "The flavonoid candidate panel"
)

dds_file <- file.path(results_dir, "dds_paired_model.rds")
time_results_file <- file.path(results_dir, "full_flavonoid_time_effects_v3.csv")
txi_file <- file.path(results_dir, "tximport_gene_level.rds")

for (required_file in c(dds_file, time_results_file, txi_file)) {
  if (!file.exists(required_file)) {
    stop("Required input file not found: ", required_file)
  }
}

dds <- readRDS(dds_file)

panel <- read.csv(
  panel_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

time_results <- read.csv(
  time_results_file,
  stringsAsFactors = FALSE
)

txi <- readRDS(txi_file)

if (is.null(txi$abundance)) {
  stop("The tximport object does not contain a gene-level TPM matrix.")
}

# ==============================================================================
# 2. Labels and metadata
# ==============================================================================

display_name <- panel$gene_name

display_name[display_name == "F3primeHa"] <- "F3'Ha"
display_name[display_name == "F3primeHb"] <- "F3'Hb"
display_name[display_name == "F3prime5primeHa"] <- "F3'5'Ha"
display_name[display_name == "F3prime5primeHb"] <- "F3'5'Hb"
display_name[display_name == "MYB123_TT2"] <- "MYB123/TT2"
display_name[display_name == "bHLH42_TT8"] <- "bHLH42/TT8"

panel$display_name <- display_name

pathway_labels <- c(
  "Phenylpropanoidweg" = "Phenylpropanoid pathway",
  "Genereller_Flavonoidweg" = "Core flavonoid pathway",
  "Hydroxylierung" = "Hydroxylation",
  "Flavonbiosynthese" = "Flavone biosynthesis",
  "Flavonolbiosynthese" = "Flavonol biosynthesis",
  "Anthocyaninbiosynthese_und_Modifikation" =
    "Anthocyanin biosynthesis/modification",
  "Proanthocyanidinbiosynthese" =
    "Proanthocyanidin biosynthesis",
  "Regulation" = "Regulation"
)

missing_pathway_labels <- setdiff(
  unique(panel$pathway_section),
  names(pathway_labels)
)

if (length(missing_pathway_labels) > 0) {
  stop(
    "No heatmap label defined for pathway section(s): ",
    paste(missing_pathway_labels, collapse = ", ")
  )
}

panel$pathway_display <- unname(pathway_labels[panel$pathway_section])

metadata <- as.data.frame(colData(dds))
metadata$sample_id <- rownames(metadata)
metadata$genotype <- factor(metadata$genotype, levels = c("AA", "Aa", "aa"))
metadata$timepoint <- factor(metadata$timepoint, levels = c("t0", "t24"))

sample_order <- with(
  metadata,
  sample_id[order(genotype, timepoint, plant_index)]
)

cat("\nSample order used for the final heatmap:\n")
print(
  data.frame(
    Sample = sample_order,
    Genotype = metadata[sample_order, "genotype"],
    Timepoint = metadata[sample_order, "timepoint"],
    stringsAsFactors = FALSE
  )
)

# ==============================================================================
# 3. Annotations and colours
# ==============================================================================

row_annotation <- data.frame(
  "Pathway section" = panel$pathway_display,
  check.names = FALSE,
  row.names = panel$display_name
)

column_annotation_A <- metadata[
  sample_order,
  c("genotype", "timepoint"),
  drop = FALSE
]
colnames(column_annotation_A) <- c("Genotype", "Timepoint")

genotype_colors <- c(
  "AA" = "#B1859F",
  "Aa" = "#D38BB4",
  "aa" = "#B4B5B6"
)

timepoint_colors <- c(
  "t0" = "#75A4A4",
  "t24" = "#D08554"
)

pathway_colors <- c(
  "Phenylpropanoid pathway" = "#C0DFDA",
  "Core flavonoid pathway" = "#BFC9D6",
  "Hydroxylation" = "#BFC9D6",
  "Flavone biosynthesis" = "#CFCECC",
  "Flavonol biosynthesis" = "#E0E0C7",
  "Anthocyanin biosynthesis/modification" = "#C8C1DA",
  "Proanthocyanidin biosynthesis" = "#A29D99",
  "Regulation" = "#95A0BC"
)

annotation_colors <- list(
  "Genotype" = genotype_colors,
  "Timepoint" = timepoint_colors,
  "Pathway section" = pathway_colors
)

heatmap_colors <- colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(101)

# ==============================================================================
# 4. Panel A: sample-wise expression heatmap
# ==============================================================================

missing_tpm_genes <- setdiff(panel$gene_id, rownames(txi$abundance))
missing_tpm_samples <- setdiff(sample_order, colnames(txi$abundance))

if (length(missing_tpm_genes) > 0) {
  stop(
    "The following panel genes are missing from the TPM matrix: ",
    paste(missing_tpm_genes, collapse = ", ")
  )
}

if (length(missing_tpm_samples) > 0) {
  stop(
    "The following samples are missing from the TPM matrix: ",
    paste(missing_tpm_samples, collapse = ", ")
  )
}

tpm_matrix <- txi$abundance[
  panel$gene_id,
  sample_order,
  drop = FALSE
]
rownames(tpm_matrix) <- panel$display_name

format_tpm <- function(current_value) {
  if (is.na(current_value)) {
    return("NA")
  }
  if (current_value < 0.005) {
    return("0")
  }
  if (current_value < 1) {
    return(formatC(current_value, format = "f", digits = 2))
  }
  if (current_value < 100) {
    return(formatC(current_value, format = "f", digits = 1))
  }
  formatC(current_value, format = "f", digits = 0)
}

tpm_display <- matrix(
  vapply(as.vector(tpm_matrix), format_tpm, character(1)),
  nrow = nrow(tpm_matrix),
  ncol = ncol(tpm_matrix),
  dimnames = dimnames(tpm_matrix)
)

vst_object <- vst(dds, blind = FALSE)

vst_matrix <- assay(vst_object)[
  panel$gene_id,
  sample_order,
  drop = FALSE
]
rownames(vst_matrix) <- panel$display_name

row_zscore <- t(
  apply(
    vst_matrix,
    1,
    function(current_row) {
      current_sd <- sd(current_row)

      if (is.na(current_sd) || current_sd == 0) {
        return(rep(NA_real_, length(current_row)))
      }

      (current_row - mean(current_row)) / current_sd
    }
  )
)

rownames(row_zscore) <- panel$display_name
colnames(row_zscore) <- sample_order

finite_z <- abs(row_zscore[is.finite(row_zscore)])
z_limit <- max(2, ceiling(max(finite_z, na.rm = TRUE)))
z_breaks <- seq(
  -z_limit,
  z_limit,
  length.out = length(heatmap_colors) + 1
)

# ==============================================================================
# 5. Panel B: genotype-specific t24-vs-t0 log2 fold changes
# ==============================================================================

comparison_to_genotype <- c(
  "AA_t24_vs_t0" = "AA",
  "Aa_t24_vs_t0" = "Aa",
  "aa_t24_vs_t0" = "aa"
)

lfc_matrix <- matrix(
  NA_real_,
  nrow = nrow(panel),
  ncol = 3,
  dimnames = list(panel$display_name, c("AA", "Aa", "aa"))
)

significance_matrix <- matrix(
  "",
  nrow = nrow(panel),
  ncol = 3,
  dimnames = list(panel$display_name, c("AA", "Aa", "aa"))
)

for (current_row in seq_len(nrow(time_results))) {
  current_gene <- time_results$gene_name[current_row]
  panel_row <- match(current_gene, panel$gene_name)

  current_genotype <- unname(
    comparison_to_genotype[time_results$comparison[current_row]]
  )

  if (is.na(panel_row) || is.na(current_genotype)) {
    next
  }

  current_display_name <- panel$display_name[panel_row]

  if (isTRUE(time_results$eligible_for_inference[current_row])) {
    lfc_matrix[current_display_name, current_genotype] <-
      time_results$log2FoldChange[current_row]

    if (
      !is.na(time_results$padj_genomewide[current_row]) &&
        time_results$padj_genomewide[current_row] < 0.05
    ) {
      significance_matrix[current_display_name, current_genotype] <- "**"
    } else if (
      !is.na(time_results$padj_goi[current_row]) &&
        time_results$padj_goi[current_row] < 0.05
    ) {
      significance_matrix[current_display_name, current_genotype] <- "*"
    }
  }
}

# Add significance stars above TPM values in t24 cells.
tpm_display_with_stars <- tpm_display

for (current_genotype in c("AA", "Aa", "aa")) {
  current_t24_samples <- sample_order[
    metadata[sample_order, "genotype"] == current_genotype &
      metadata[sample_order, "timepoint"] == "t24"
  ]

  current_sig <- significance_matrix[, current_genotype]

  for (current_sample in current_t24_samples) {
    tpm_display_with_stars[, current_sample] <- ifelse(
      current_sig != "",
      paste0(current_sig, "\n", tpm_display[, current_sample]),
      tpm_display[, current_sample]
    )
  }
}

finite_lfc <- lfc_matrix[is.finite(lfc_matrix)]
lfc_limit <- max(2, ceiling(max(abs(finite_lfc), na.rm = TRUE)))
lfc_breaks <- seq(
  -lfc_limit,
  lfc_limit,
  length.out = length(heatmap_colors) + 1
)

# ==============================================================================
# 6. Plot builders
# ==============================================================================

make_panel_A <- function(output_file = NA_character_, silent = TRUE) {
  pheatmap(
    row_zscore,
    filename = if (is.na(output_file)) NA_character_ else output_file,
    silent = silent,
    width = 20,
    height = 19,
    color = heatmap_colors,
    breaks = z_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = column_annotation_A,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    annotation_names_row = FALSE,
    border_color = "white",
    na_col = "#D9D9D9",
    display_numbers = tpm_display_with_stars,
    number_color = "black",
    fontsize_number = 5.8,
    angle_col = 45,
    fontsize = 10,
    fontsize_row = 9.5,
    fontsize_col = 9.2,
    gaps_col = c(3, 6, 9, 12, 15),
    legend = FALSE,
    annotation_legend = FALSE,
    main = paste0(
      "(A) Expression Heatmap\n",
      "Colors: row-wise z-scored VST expression; values: TPM\n",
      "Stars in t24 cells indicate the genotype-level t24-vs-t0 contrast"
    )
  )
}

make_panel_B <- function(output_file = NA_character_, silent = TRUE) {
  pheatmap(
    lfc_matrix,
    filename = if (is.na(output_file)) NA_character_ else output_file,
    silent = silent,
    width = 6,
    height = 19,
    color = heatmap_colors,
    breaks = lfc_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    border_color = "white",
    display_numbers = significance_matrix,
    number_color = "black",
    fontsize_number = 16,
    fontsize = 10,
    fontsize_col = 9.5,
    angle_col = 0,
    na_col = "#D9D9D9",
    main = paste0(
      "(B) log2 Fold Change\n",
      "t24 vs. t0\n",
      "* panel FDR < 0.05; ** genome-wide FDR < 0.05"
    )
  )
}

# ==============================================================================
# 7. Custom label backgrounds
# ==============================================================================

add_row_label_backgrounds <- function(pheatmap_object, fill_colors) {
  idx <- which(pheatmap_object$gtable$layout$name == "row_names")

  if (length(idx) != 1) {
    stop("Could not uniquely identify the row_names grob in Panel A.")
  }

  old_labels <- pheatmap_object$gtable$grobs[[idx]]
  n_labels <- length(old_labels$label)

  if (length(fill_colors) != n_labels) {
    stop(
      "Number of row-label colours does not match the number of gene labels."
    )
  }

  old_labels$gp$fontface <- "italic"

  background <- rectGrob(
    x = unit(rep(0.5, n_labels), "npc"),
    y = old_labels$y,
    width = unit(rep(1, n_labels), "npc"),
    height = unit(rep(0.96 / n_labels, n_labels), "npc"),
    gp = gpar(fill = fill_colors, col = "white", lwd = 0.6)
  )

  pheatmap_object$gtable$grobs[[idx]] <- grobTree(background, old_labels)
  pheatmap_object
}

add_panelA_col_label_backgrounds <- function(pheatmap_object, fill_colors) {
  idx <- which(pheatmap_object$gtable$layout$name == "col_names")

  if (length(idx) != 1) {
    stop("Could not uniquely identify the col_names grob in Panel A.")
  }

  old_labels <- pheatmap_object$gtable$grobs[[idx]]
  n_labels <- length(old_labels$label)

  if (length(fill_colors) != n_labels) {
    stop(
      paste0(
        "Number of column-label colours does not match the number of ",
        "sample labels."
      )
    )
  }

  background <- rectGrob(
    x = old_labels$x,
    y = unit(rep(0.52, n_labels), "npc"),
    width = unit(rep(0.97 / n_labels, n_labels), "npc"),
    height = unit(rep(0.98, n_labels), "npc"),
    gp = gpar(fill = fill_colors, col = "white", lwd = 0.6)
  )

  old_labels$gp$col <- "black"
  old_labels$gp$fontsize <- 9.2
  old_labels$gp$fontface <- "bold"

  pheatmap_object$gtable$grobs[[idx]] <- grobTree(background, old_labels)
  pheatmap_object
}

add_panelB_col_label_backgrounds <- function(pheatmap_object, fill_colors) {
  idx <- which(pheatmap_object$gtable$layout$name == "col_names")

  if (length(idx) != 1) {
    stop("Could not uniquely identify the col_names grob in Panel B.")
  }

  old_labels <- pheatmap_object$gtable$grobs[[idx]]
  n_labels <- length(old_labels$label)

  if (length(fill_colors) != n_labels) {
    stop(
      paste0(
        "Number of column-label colours does not match the number of ",
        "genotype labels."
      )
    )
  }

  background <- rectGrob(
    x = old_labels$x,
    y = unit(rep(0.5, n_labels), "npc"),
    width = unit(rep(0.96 / n_labels, n_labels), "npc"),
    height = unit(rep(0.95, n_labels), "npc"),
    gp = gpar(fill = fill_colors, col = "white", lwd = 0.8)
  )

  foreground <- textGrob(
    label = old_labels$label,
    x = old_labels$x,
    y = unit(rep(0.5, n_labels), "npc"),
    hjust = 0.5,
    vjust = 0.5,
    gp = gpar(fontsize = 9.5, fontface = "bold", col = "black")
  )

  pheatmap_object$gtable$grobs[[idx]] <- grobTree(background, foreground)
  pheatmap_object
}

# ==============================================================================
# 8. Save function
# ==============================================================================

save_grob_multi <- function(
    grob_obj,
    png_file,
    pdf_file,
    svg_file,
    width,
    height
) {
  png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  grid.newpage()
  grid.draw(grob_obj)
  dev.off()

  pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  grid.newpage()
  grid.draw(grob_obj)
  dev.off()

  svg(
    filename = svg_file,
    width = width,
    height = height,
    pointsize = 12
  )
  grid.newpage()
  grid.draw(grob_obj)
  dev.off()
}

# ==============================================================================
# 9. Build final figure
# ==============================================================================

panel_A_obj <- make_panel_A(NA_character_, silent = TRUE)
panel_B_obj <- make_panel_B(NA_character_, silent = TRUE)

row_label_colors <- unname(pathway_colors[panel$pathway_display])
panel_A_obj <- add_row_label_backgrounds(
  panel_A_obj,
  row_label_colors
)

panel_A_sample_label_colors <- unname(
  genotype_colors[metadata[sample_order, "genotype"]]
)
panel_A_obj <- add_panelA_col_label_backgrounds(
  panel_A_obj,
  panel_A_sample_label_colors
)

panel_B_label_colors <- unname(genotype_colors[colnames(lfc_matrix)])
panel_B_obj <- add_panelB_col_label_backgrounds(
  panel_B_obj,
  panel_B_label_colors
)

common_n_heights <- min(
  length(panel_A_obj$gtable$heights),
  length(panel_B_obj$gtable$heights)
)

panel_B_obj$gtable$heights[seq_len(common_n_heights)] <-
  panel_A_obj$gtable$heights[seq_len(common_n_heights)]

figure_landscape <- arrangeGrob(
  grobs = list(panel_A_obj$gtable, panel_B_obj$gtable),
  ncol = 2,
  widths = c(5.4, 1.7)
)

# ==============================================================================
# 10. Annotation legend
# ==============================================================================

make_legend_item <- function(label, color, fontsize = 8) {
  arrangeGrob(
    rectGrob(
      width = unit(0.16, "in"),
      height = unit(0.16, "in"),
      gp = gpar(fill = color, col = NA)
    ),
    textGrob(
      label,
      x = unit(0, "npc"),
      just = "left",
      gp = gpar(fontsize = fontsize)
    ),
    ncol = 2,
    widths = unit.c(unit(0.21, "in"), unit(1, "null"))
  )
}

make_legend_group <- function(
    title,
    labels,
    colors,
    ncol = length(labels),
    fontsize = 8
) {
  items <- lapply(
    seq_along(labels),
    function(i) {
      make_legend_item(labels[i], colors[i], fontsize = fontsize)
    }
  )

  body <- arrangeGrob(grobs = items, ncol = ncol)

  arrangeGrob(
    textGrob(
      title,
      x = unit(0, "npc"),
      just = "left",
      gp = gpar(fontsize = 9, fontface = "bold")
    ),
    body,
    ncol = 1,
    heights = c(0.35, 1)
  )
}

legend_timepoint <- make_legend_group(
  "Timepoint",
  labels = names(timepoint_colors),
  colors = unname(timepoint_colors),
  ncol = 1,
  fontsize = 8.8
)

legend_genotype <- make_legend_group(
  "Genotype",
  labels = names(genotype_colors),
  colors = unname(genotype_colors),
  ncol = 1,
  fontsize = 8.8
)

legend_pathway <- make_legend_group(
  "Pathway section",
  labels = names(pathway_colors),
  colors = unname(pathway_colors),
  ncol = 2,
  fontsize = 8.2
)

annotation_legend <- arrangeGrob(
  textGrob(
    "Annotation legend",
    x = unit(0, "npc"),
    just = "left",
    gp = gpar(fontsize = 10.5, fontface = "bold")
  ),
  arrangeGrob(
    legend_timepoint,
    legend_genotype,
    legend_pathway,
    ncol = 3,
    widths = c(0.75, 0.85, 2.6)
  ),
  ncol = 1,
  heights = c(0.28, 1)
)

figure_with_legend <- arrangeGrob(
  figure_landscape,
  annotation_legend,
  ncol = 1,
  heights = c(16, 2.8)
)

# ==============================================================================
# 11. Save final outputs
# ==============================================================================

output_prefix <- "full_flavonoid_expression_heatmap"

png_file <- file.path(plots_dir, paste0(output_prefix, ".png"))
pdf_file <- file.path(plots_dir, paste0(output_prefix, ".pdf"))
svg_file <- file.path(final_plots_dir, paste0(output_prefix, ".svg"))

save_grob_multi(
  grob_obj = figure_with_legend,
  png_file = png_file,
  pdf_file = pdf_file,
  svg_file = svg_file,
  width = 20,
  height = 19
)

cat("\nSUCCESS: flavonoid expression heatmap created.\n")
cat("PNG: ", png_file, "\n", sep = "")
cat("PDF: ", pdf_file, "\n", sep = "")
cat("SVG: ", svg_file, "\n", sep = "")
