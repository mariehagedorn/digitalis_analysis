############################################################
# Digitalis purpurea RNA-seq
# DESeq2 quality control
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

stopifnot(
  identical(
    colnames(txi$counts),
    rownames(samples)
  )
)


############################################################
# 4. Create DESeq2 dataset for quality control
#
# This design is used here only for construction of the
# DESeqDataSet. The VST below is calculated with blind = TRUE.
############################################################

dds_qc <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ genotype * timepoint
)


############################################################
# 5. Remove genes with very low counts
#
# A gene is retained if it has at least 10 estimated counts
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
# 6. Variance-stabilizing transformation
############################################################

vsd <- vst(
  dds_qc,
  blind = TRUE
)


############################################################
# 7. PCA calculation
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

write_csv(
  pca_data,
  file.path(
    results_dir,
    "PCA_coordinates.csv"
  )
)


############################################################
# 8. PCA plot
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

ggsave(
  filename = file.path(
    plots_dir,
    "PCA_DESeq2.png"
  ),
  plot = pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(
    plots_dir,
    "PCA_DESeq2.pdf"
  ),
  plot = pca_plot,
  width = 9,
  height = 7
)


############################################################
# 9. Sample-distance heatmap
############################################################

sample_distances <- dist(
  t(
    assay(vsd)
  )
)

sample_distance_matrix <- as.matrix(
  sample_distances
)

annotation <- data.frame(
  Genotype = samples$genotype,
  Timepoint = samples$timepoint
)

rownames(annotation) <- rownames(
  samples
)

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

png(
  filename = file.path(
    plots_dir,
    "Sample_distance_heatmap.png"
  ),
  width = 2400,
  height = 2100,
  res = 300
)

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = annotation_colors,
  main = "Sample distances based on VST counts",
  border_color = NA
)

dev.off()


############################################################
# Save final sample-distance heatmap as SVG
############################################################

grDevices::svg(
  filename = file.path(
    final_plots_dir,
    "Sample_distance_heatmap_final.svg"
  ),
  width = 8,
  height = 7,
  pointsize = 12
)

pheatmap(
  sample_distance_matrix,
  annotation_col = annotation,
  annotation_row = annotation,
  annotation_colors = annotation_colors,
  main = "Sample distances based on VST counts",
  border_color = NA
)

grDevices::dev.off()


############################################################
# 10. Save processed objects
############################################################

saveRDS(
  dds_qc,
  file.path(
    results_dir,
    "dds_quality_control.rds"
  )
)

saveRDS(
  vsd,
  file.path(
    results_dir,
    "vst_quality_control.rds"
  )
)


############################################################
# 11. Completion message
############################################################

cat(
  "\nQuality control completed successfully.\n"
)

cat(
  "PCA variance:\n"
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
