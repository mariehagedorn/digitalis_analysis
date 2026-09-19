############################################################
# Digitalis purpurea RNA-seq
# DESeq2 quality control
#
# PCA and sample-distance heatmap based on
# variance-stabilized gene-level expression data.
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readr)
})


############################################################
# 2. Define paths
#
# The script assumes that it is run from the root directory
# of the GitHub repository.
############################################################

project_dir <- "RNAseq"

metadata_dir <- file.path(
  project_dir,
  "metadata"
)

results_dir <- file.path(
  project_dir,
  "results"
)

plots_dir <- file.path(
  project_dir,
  "plots"
)

final_plots_dir <- file.path(
  plots_dir,
  "final"
)

# Create output directories if they do not already exist
dir.create(
  metadata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  final_plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Load tximport data and sample metadata
############################################################

txi <- readRDS(
  file.path(
    results_dir,
    "tximport_gene_level.rds"
  )
)

samples <- readRDS(
  file.path(
    metadata_dir,
    "samples.rds"
  )
)


############################################################
# 4. Validate metadata and sample order
############################################################

required_columns <- c(
  "sample_id",
  "raw_id",
  "plant_id",
  "genotype",
  "timepoint",
  "kallisto_file"
)

stopifnot(
  all(
    required_columns %in% colnames(samples)
  )
)

stopifnot(
  nrow(samples) == 18
)

stopifnot(
  identical(
    colnames(txi$counts),
    rownames(samples)
  )
)


############################################################
# 5. Define factor order
############################################################

samples$genotype <- factor(
  samples$genotype,
  levels = c(
    "AA",
    "Aa",
    "aa"
  )
)

samples$timepoint <- factor(
  samples$timepoint,
  levels = c(
    "t0",
    "t24"
  )
)

samples$plant_id <- factor(
  samples$plant_id
)


############################################################
# 6. Create DESeq2 dataset for quality control
#
# This design is used only for creating the DESeqDataSet.
# The VST below is calculated with blind = TRUE.
############################################################

dds_qc <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ genotype * timepoint
)


############################################################
# 7. Filter genes with very low counts
#
# Keep genes with at least 10 estimated counts
# in at least three samples.
############################################################

keep <- rowSums(
  counts(dds_qc) >= 10
) >= 3

cat(
  "Genes before filtering:",
  nrow(dds_qc),
  "\n"
)

cat(
  "Genes after filtering:",
  sum(keep),
  "\n"
)

dds_qc <- dds_qc[
  keep,
]


############################################################
# 8. Variance-stabilizing transformation
############################################################

vsd <- vst(
  dds_qc,
  blind = TRUE
)


############################################################
# 9. Calculate PCA
############################################################

pca_data <- plotPCA(
  vsd,
  intgroup = c(
    "genotype",
    "timepoint"
  ),
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(
    pca_data,
    "percentVar"
  ),
  digits = 1
)

pca_data$sample_id <- rownames(
  pca_data
)


############################################################
# 10. Save PCA coordinates
############################################################

write_csv(
  pca_data,
  file.path(
    results_dir,
    "PCA_coordinates_v4.csv"
  )
)


############################################################
# 11. Create PCA plot
############################################################

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = genotype,
    shape = timepoint,
    label = sample_id
  )
) +
  geom_point(
    size = 4
  ) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "AA" = "#C00000",
      "Aa" = "#F7B6D2",
      "aa" = "grey50"
    )
  ) +
  labs(
    title = "PCA of D. purpurea RNA-seq samples",
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    ),
    color = "Genotype",
    shape = "Timepoint"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


############################################################
# 12. Save PCA as PNG
############################################################

ggsave(
  filename = file.path(
    plots_dir,
    "PCA_DESeq2_v4.png"
  ),
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)


############################################################
# 13. Save PCA as PDF
############################################################

ggsave(
  filename = file.path(
    plots_dir,
    "PCA_DESeq2_v4.pdf"
  ),
  plot = pca_plot,
  width = 9,
  height = 7
)


############################################################
# 13b. Save final PCA as SVG
############################################################

grDevices::svg(
  filename = file.path(
    final_plots_dir,
    "PCA_DESeq2_final.svg"
  ),
  width = 9,
  height = 7,
  pointsize = 12
)

print(pca_plot)

grDevices::dev.off()


############################################################
# 14. Calculate VST-based sample distances
############################################################

sample_distances <- dist(
  t(
    assay(vsd)
  )
)

sample_distance_matrix <- as.matrix(
  sample_distances
)


############################################################
# 15. Prepare heatmap annotations
############################################################

annotation <- data.frame(
  Genotype = samples$genotype,
  Timepoint = samples$timepoint
)

rownames(annotation) <- rownames(
  samples
)

annotation <- annotation[
  colnames(sample_distance_matrix),
  ,
  drop = FALSE
]


############################################################
# 16. Define annotation colors
############################################################

annotation_colors <- list(

  Genotype = c(
    "AA" = "#C00000",
    "Aa" = "#F7B6D2",
    "aa" = "grey50"
  ),

  Timepoint = c(
    "t0" = "#377EB8",
    "t24" = "#FFB000"
  )

)


############################################################
# 17. Define heatmap colors
#
# This reproduces the default pheatmap color palette:
# blue = low distance / high similarity
# red  = high distance / low similarity
############################################################

distance_colors <- colorRampPalette(
  rev(
    RColorBrewer::brewer.pal(
      7,
      "RdYlBu"
    )
  )
)(
  100
)


############################################################
# 18. Function for drawing the heatmap
#
# Dendrograms are calculated directly from the
# displayed VST-based sample distances.
############################################################

draw_sample_distance_heatmap <- function() {

  pheatmap(

    sample_distance_matrix,

    color = distance_colors,

    clustering_distance_rows = sample_distances,

    clustering_distance_cols = sample_distances,

    annotation_col = annotation,

    annotation_row = annotation,

    annotation_colors = annotation_colors,

    main = "Sample distances based on VST counts",

    border_color = NA

  )

}


############################################################
# 19. Save heatmap as PNG
############################################################

png(
  filename = file.path(
    plots_dir,
    "Sample_distance_heatmap_v4.png"
  ),
  width = 2400,
  height = 2100,
  res = 300
)

draw_sample_distance_heatmap()

dev.off()


############################################################
# 20. Save heatmap as PDF
############################################################

pdf(
  file = file.path(
    plots_dir,
    "Sample_distance_heatmap_v4.pdf"
  ),
  width = 8,
  height = 7
)

draw_sample_distance_heatmap()

dev.off()


############################################################
# 21. Save processed QC objects
############################################################

saveRDS(
  dds_qc,
  file.path(
    results_dir,
    "dds_quality_control_v4.rds"
  )
)

saveRDS(
  vsd,
  file.path(
    results_dir,
    "vst_quality_control_v4.rds"
  )
)


############################################################
# 22. Completion message
############################################################

cat(
  "\nQuality control completed successfully.\n"
)

cat(
  "Genes retained:",
  nrow(dds_qc),
  "\n"
)

cat(
  "PC1:",
  percent_variance[1],
  "%\n"
)

cat(
  "PC2:",
  percent_variance[2],
  "%\n"
)

cat(
  "\nCreated files:\n"
)

cat(
  file.path(
    plots_dir,
    "PCA_DESeq2_v4.png"
  ),
  "\n"
)

cat(
  file.path(
    plots_dir,
    "PCA_DESeq2_v4.pdf"
  ),
  "\n"
)

cat(
  file.path(
    plots_dir,
    "Sample_distance_heatmap_v4.png"
  ),
  "\n"
)

cat(
  file.path(
    plots_dir,
    "Sample_distance_heatmap_v4.pdf"
  ),
  "\n"
)
