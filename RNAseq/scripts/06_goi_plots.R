#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea RNA-seq
# Visualization of genes of interest
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})


############################################################
# 2. Define paths
#
# The script assumes that it is run from the root directory
# of the GitHub repository.
############################################################

base_dir <- "RNAseq"

results_dir <- file.path(
  base_dir,
  "results"
)

metadata_dir <- file.path(
  base_dir,
  "metadata"
)

plots_dir <- file.path(
  base_dir,
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
# 3. Load data
############################################################

dds <- readRDS(
  file.path(
    results_dir,
    "dds_paired_model.rds"
  )
)

normalized_counts <- readRDS(
  file.path(
    results_dir,
    "normalized_counts_all_genes.rds"
  )
)

goi <- read.csv(
  file.path(
    metadata_dir,
    "genes_of_interest.csv"
  ),
  stringsAsFactors = FALSE
)

goi_filter <- read.csv(
  file.path(
    results_dir,
    "GOI_filter_status.csv"
  ),
  stringsAsFactors = FALSE
)

time_results <- read.csv(
  file.path(
    results_dir,
    "GOI_time_effects.csv"
  ),
  stringsAsFactors = FALSE
)


############################################################
# 4. Prepare sample metadata
############################################################

metadata <- as.data.frame(
  colData(dds)
)

metadata$sample_id <- rownames(
  metadata
)

metadata$genotype <- factor(
  metadata$genotype,
  levels = c(
    "AA",
    "Aa",
    "aa"
  )
)

metadata$timepoint <- factor(
  metadata$timepoint,
  levels = c(
    "t0",
    "t24"
  )
)

sample_order <- with(
  metadata,
  sample_id[
    order(
      genotype,
      plant_index,
      timepoint
    )
  ]
)


############################################################
# 5. Determine expression status of genes of interest
############################################################

filter_index <- match(
  goi$gene_id,
  goi_filter$gene_id
)

expression_status <- ifelse(
  goi_filter$total_count[
    filter_index
  ] == 0,
  "not_detected",
  ifelse(
    goi_filter$passes_standard_filter[
      filter_index
    ],
    "passes_filter",
    "low_expression"
  )
)

row_annotation <- data.frame(
  Category = goi$category,
  Expression = expression_status,
  row.names = goi$gene_name
)


############################################################
# 6. Convert normalized counts to long format
############################################################

normalized_goi <- normalized_counts[
  goi$gene_id,
  ,
  drop = FALSE
]

long_data <- as.data.frame(
  as.table(
    normalized_goi
  ),
  responseName = "normalized_count",
  stringsAsFactors = FALSE
)

colnames(
  long_data
)[1:2] <- c(
  "gene_id",
  "sample_id"
)

long_data$gene_name <- goi$gene_name[
  match(
    long_data$gene_id,
    goi$gene_id
  )
]

long_data$category <- goi$category[
  match(
    long_data$gene_id,
    goi$gene_id
  )
]

metadata_index <- match(
  long_data$sample_id,
  metadata$sample_id
)

long_data$plant_id <- as.character(
  metadata$plant_id[
    metadata_index
  ]
)

long_data$genotype <- factor(
  metadata$genotype[
    metadata_index
  ],
  levels = c(
    "AA",
    "Aa",
    "aa"
  )
)

long_data$timepoint <- factor(
  metadata$timepoint[
    metadata_index
  ],
  levels = c(
    "t0",
    "t24"
  )
)

long_data$log2_normalized_count <- log2(
  long_data$normalized_count + 1
)

long_data$gene_name <- factor(
  long_data$gene_name,
  levels = goi$gene_name
)

write.csv(
  long_data,
  file.path(
    results_dir,
    "GOI_normalized_counts_long.csv"
  ),
  row.names = FALSE
)


############################################################
# 7. Calculate mean, SD, and SE per genotype and timepoint
############################################################

split_groups <- split(
  long_data,
  interaction(
    long_data$gene_id,
    long_data$genotype,
    long_data$timepoint,
    drop = TRUE
  )
)

expression_summary <- do.call(
  rbind,
  lapply(
    split_groups,
    function(current_data) {

      current_values <-
        current_data$normalized_count

      data.frame(
        gene_id =
          current_data$gene_id[1],

        gene_name =
          as.character(
            current_data$gene_name[1]
          ),

        category =
          current_data$category[1],

        genotype =
          as.character(
            current_data$genotype[1]
          ),

        timepoint =
          as.character(
            current_data$timepoint[1]
          ),

        n =
          length(
            current_values
          ),

        mean_normalized_count =
          mean(
            current_values
          ),

        sd_normalized_count =
          sd(
            current_values
          ),

        se_normalized_count =
          sd(
            current_values
          ) /
          sqrt(
            length(
              current_values
            )
          )
      )
    }
  )
)

rownames(
  expression_summary
) <- NULL

expression_summary <- expression_summary[
  order(
    match(
      expression_summary$gene_id,
      goi$gene_id
    ),
    match(
      expression_summary$genotype,
      c(
        "AA",
        "Aa",
        "aa"
      )
    ),
    match(
      expression_summary$timepoint,
      c(
        "t0",
        "t24"
      )
    )
  ),
]

write.csv(
  expression_summary,
  file.path(
    results_dir,
    "GOI_expression_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 8. Define colors
############################################################

genotype_colors <- c(
  "AA" = "#C00000",
  "Aa" = "#F3A6C8",
  "aa" = "#7F7F7F"
)

timepoint_colors <- c(
  "t0" = "#C7B5E3",
  "t24" = "#5B2A86"
)

annotation_colors <- list(

  Genotype = genotype_colors,

  Timepoint = timepoint_colors,

  Category = c(
    "structural" = "#4C78A8",
    "transcription_factor" = "#B279A2"
  ),

  Expression = c(
    "passes_filter" = "#2E8B57",
    "low_expression" = "#F0A202",
    "not_detected" = "#6C757D"
  )
)

heatmap_colors <- colorRampPalette(
  c(
    "#2166AC",
    "#F7F7F7",
    "#B2182B"
  )
)(101)


############################################################
# 9. Paired expression plots
#
# Each line represents one individual plant.
# The y-axis is scaled independently for each gene.
############################################################

paired_plot <- ggplot(
  long_data,
  aes(
    x = timepoint,
    y = log2_normalized_count,
    group = plant_id,
    color = genotype
  )
) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(
    aes(
      shape = timepoint
    ),
    size = 2.5
  ) +
  facet_wrap(
    ~ gene_name,
    ncol = 4,
    scales = "free_y",
    drop = FALSE
  ) +
  scale_color_manual(
    values = genotype_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "t0" = 16,
      "t24" = 17
    )
  ) +
  labs(
    title =
      "Expression of anthocyanin-related genes",
    subtitle =
      "Paired D. purpurea leaf samples before and after 24 h",
    x = "Timepoint",
    y = expression(
      log[2] * "(normalized count + 1)"
    ),
    color = "Genotype",
    shape = "Timepoint"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


############################################################
# 10. Save paired expression plots
############################################################

ggsave(
  filename = file.path(
    plots_dir,
    "GOI_paired_expression.png"
  ),
  plot = paired_plot,
  width = 15,
  height = 12,
  dpi = 300
)

ggsave(
  filename = file.path(
    plots_dir,
    "GOI_paired_expression.pdf"
  ),
  plot = paired_plot,
  width = 15,
  height = 12
)

grDevices::svg(
  filename = file.path(
    final_plots_dir,
    "GOI_paired_expression_final.svg"
  ),
  width = 15,
  height = 12,
  pointsize = 12
)

print(
  paired_plot
)

grDevices::dev.off()


############################################################
# 11. VST heatmap of genes of interest
#
# Expression values are standardized separately for each
# gene using row-wise z-scores.
#
# The heatmap therefore visualizes expression patterns
# rather than absolute differences in expression between
# genes.
############################################################

vst_object <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vst_object
)[
  goi$gene_id,
  sample_order,
  drop = FALSE
]

rownames(
  vst_matrix
) <- goi$gene_name

row_zscore <- t(
  apply(
    vst_matrix,
    1,
    function(
      current_row
    ) {

      current_sd <- sd(
        current_row
      )

      if (
        is.na(
          current_sd
        ) ||
        current_sd == 0
      ) {
        return(
          rep(
            0,
            length(
              current_row
            )
          )
        )
      }

      (
        current_row -
          mean(
            current_row
          )
      ) /
        current_sd
    }
  )
)

rownames(
  row_zscore
) <- goi$gene_name

colnames(
  row_zscore
) <- sample_order

column_annotation <- metadata[
  sample_order,
  c(
    "genotype",
    "timepoint"
  ),
  drop = FALSE
]

colnames(
  column_annotation
) <- c(
  "Genotype",
  "Timepoint"
)


############################################################
# 12. Save VST heatmap
############################################################

pheatmap(
  row_zscore,
  filename = file.path(
    plots_dir,
    "GOI_VST_heatmap.png"
  ),
  width = 14,
  height = 9,
  color = heatmap_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = column_annotation,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  border_color = NA,
  angle_col = 45,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 8,
  main =
    "Horz et al. anthocyanin-related genes: row z-scores"
)

pheatmap(
  row_zscore,
  filename = file.path(
    plots_dir,
    "GOI_VST_heatmap.pdf"
  ),
  width = 14,
  height = 9,
  color = heatmap_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = column_annotation,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  border_color = NA,
  angle_col = 45,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 8,
  main =
    "Horz et al. anthocyanin-related genes: row z-scores"
)


############################################################
# 13. Create log2-fold-change matrix for time responses
############################################################

comparison_to_genotype <- c(
  "AA_t24_vs_t0" = "AA",
  "Aa_t24_vs_t0" = "Aa",
  "aa_t24_vs_t0" = "aa"
)

lfc_matrix <- matrix(
  NA_real_,
  nrow = nrow(
    goi
  ),
  ncol = 3,
  dimnames = list(
    goi$gene_name,
    c(
      "AA",
      "Aa",
      "aa"
    )
  )
)

significance_matrix <- matrix(
  "",
  nrow = nrow(
    goi
  ),
  ncol = 3,
  dimnames = list(
    goi$gene_name,
    c(
      "AA",
      "Aa",
      "aa"
    )
  )
)

for (
  current_row in seq_len(
    nrow(
      time_results
    )
  )
) {

  current_gene <-
    time_results$gene_name[
      current_row
    ]

  current_genotype <-
    comparison_to_genotype[
      time_results$comparison[
        current_row
      ]
    ]

  lfc_matrix[
    current_gene,
    current_genotype
  ] <-
    time_results$log2FoldChange[
      current_row
    ]

  if (
    !is.na(
      time_results$padj_genomewide[
        current_row
      ]
    ) &&
      time_results$padj_genomewide[
        current_row
      ] < 0.05
  ) {

    significance_matrix[
      current_gene,
      current_genotype
    ] <- "**"

  } else if (
    !is.na(
      time_results$padj_goi[
        current_row
      ]
    ) &&
      time_results$padj_goi[
        current_row
      ] < 0.05
  ) {

    significance_matrix[
      current_gene,
      current_genotype
    ] <- "*"
  }
}

lfc_limit <- max(
  2,
  ceiling(
    max(
      abs(
        lfc_matrix
      ),
      na.rm = TRUE
    )
  )
)

lfc_breaks <- seq(
  -lfc_limit,
  lfc_limit,
  length.out = 102
)


############################################################
# 14. Save log2-fold-change heatmap
############################################################

pheatmap(
  lfc_matrix,
  filename = file.path(
    plots_dir,
    "GOI_time_log2FC_heatmap.png"
  ),
  width = 7,
  height = 10,
  color = heatmap_colors,
  breaks = lfc_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  border_color = "white",
  display_numbers = significance_matrix,
  number_color = "black",
  fontsize_number = 13,
  fontsize = 11,
  na_col = "#D9D9D9",
  main =
    "t24 vs. t0: * GOI FDR, ** genome-wide FDR"
)

pheatmap(
  lfc_matrix,
  filename = file.path(
    plots_dir,
    "GOI_time_log2FC_heatmap.pdf"
  ),
  width = 7,
  height = 10,
  color = heatmap_colors,
  breaks = lfc_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  border_color = "white",
  display_numbers = significance_matrix,
  number_color = "black",
  fontsize_number = 13,
  fontsize = 11,
  na_col = "#D9D9D9",
  main =
    "t24 vs. t0: * GOI FDR, ** genome-wide FDR"
)


############################################################
# 15. Completion message
############################################################

cat(
  "\nGene-of-interest visualizations were created successfully:\n"
)

cat(
  "- GOI_paired_expression.png and .pdf\n"
)

cat(
  "- GOI_VST_heatmap.png and .pdf\n"
)

cat(
  "- GOI_time_log2FC_heatmap.png and .pdf\n"
)

cat(
  "- GOI_expression_summary.csv\n"
)
