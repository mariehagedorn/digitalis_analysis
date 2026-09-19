############################################################
# Digitalis purpurea RNA-seq
# Initial Kallisto quality-control PCA
#
# Input:
# Individual Kallisto transcript-abundance files
#
# PCA basis:
# log2(TPM + 1)
# 500 transcripts with the highest variance
#
# Output:
# PNG and SVG
############################################################


############################################################
# 1. Load packages
############################################################

library(ggplot2)
library(ggrepel)


############################################################
# 2. Define input and output files
############################################################

metadata_file <- paste0(
  "/vol/data/digitalis_rnaseq/deseq2/",
  "metadata/samples.csv"
)

png_file <- paste0(
  "/vol/data/digitalis_rnaseq/deseq2/",
  "plots/PCA_initial_Kallisto_QC.png"
)

svg_file <- paste0(
  "/vol/data/digitalis_rnaseq/deseq2/",
  "plots/PCA_initial_Kallisto_QC.svg"
)


############################################################
# 3. Import sample metadata
############################################################

samples <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample_id",
  "genotype",
  "timepoint",
  "kallisto_file"
)

missing_columns <- setdiff(
  required_columns,
  colnames(samples)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(samples$sample_id)) {
  stop("Duplicate sample IDs were found.")
}

if (!all(file.exists(samples$kallisto_file))) {
  stop(
    "At least one Kallisto file listed in samples.csv ",
    "could not be found."
  )
}

samples$genotype <- factor(
  samples$genotype,
  levels = c("AA", "Aa", "aa")
)

samples$timepoint <- factor(
  samples$timepoint,
  levels = c("t0", "t24")
)


############################################################
# 4. Function for reading Kallisto TPM values
############################################################

read_kallisto_tpm <- function(file_path) {

  abundance <- read.delim(
    gzfile(file_path),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_kallisto_columns <- c(
    "target_id",
    "tpm"
  )

  if (!all(required_kallisto_columns %in% colnames(abundance))) {
    stop(
      "The required columns target_id and tpm are missing in: ",
      file_path
    )
  }

  if (anyDuplicated(abundance$target_id)) {
    stop(
      "Duplicate transcript IDs were found in: ",
      file_path
    )
  }

  tpm_values <- abundance$tpm
  names(tpm_values) <- abundance$target_id

  return(tpm_values)
}


############################################################
# 5. Create transcript-level TPM matrix
############################################################

first_sample <- read_kallisto_tpm(
  samples$kallisto_file[1]
)

transcript_ids <- names(first_sample)

tpm_matrix <- matrix(
  NA_real_,
  nrow = length(transcript_ids),
  ncol = nrow(samples),
  dimnames = list(
    transcript_ids,
    samples$sample_id
  )
)

tpm_matrix[, 1] <- first_sample[transcript_ids]

if (nrow(samples) > 1) {

  for (i in 2:nrow(samples)) {

    sample_tpm <- read_kallisto_tpm(
      samples$kallisto_file[i]
    )

    if (!setequal(names(sample_tpm), transcript_ids)) {
      stop(
        "The transcript IDs differ in sample: ",
        samples$sample_id[i]
      )
    }

    tpm_matrix[, i] <- sample_tpm[transcript_ids]
  }
}

if (anyNA(tpm_matrix)) {
  stop("Missing TPM values were found in the TPM matrix.")
}


############################################################
# 6. Filter weakly expressed transcripts
#
# Retain transcripts with TPM >= 1
# in at least three samples
############################################################

keep_transcript <- rowSums(
  tpm_matrix >= 1
) >= 3

filtered_tpm <- tpm_matrix[
  keep_transcript,
  ,
  drop = FALSE
]

if (nrow(filtered_tpm) < 2) {
  stop("Too few transcripts remained after filtering.")
}


############################################################
# 7. Log2 transformation
############################################################

log_tpm <- log2(
  filtered_tpm + 1
)


############################################################
# 8. Select the 500 most variable transcripts
############################################################

transcript_variance <- apply(
  log_tpm,
  1,
  var
)

transcript_variance[
  !is.finite(transcript_variance)
] <- 0

number_top_transcripts <- min(
  500,
  nrow(log_tpm)
)

top_transcripts <- names(
  sort(
    transcript_variance,
    decreasing = TRUE
  )
)[seq_len(number_top_transcripts)]

pca_input <- t(
  log_tpm[
    top_transcripts,
    ,
    drop = FALSE
  ]
)


############################################################
# 9. Calculate PCA
############################################################

pca_result <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

explained_variance <- (
  pca_result$sdev^2 /
  sum(pca_result$sdev^2)
) * 100


############################################################
# 10. Prepare PCA coordinates
############################################################

pca_data <- data.frame(
  sample_id = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  stringsAsFactors = FALSE
)

pca_data <- merge(
  samples,
  pca_data,
  by = "sample_id",
  sort = FALSE
)

pca_data <- pca_data[
  match(samples$sample_id, pca_data$sample_id),
]

stopifnot(
  identical(
    as.character(pca_data$sample_id),
    as.character(samples$sample_id)
  )
)


############################################################
# 11. Create PCA plot
############################################################

genotype_colors <- c(
  "AA" = "#CC0000",
  "Aa" = "#F5ACC8",
  "aa" = "#888888"
)

timepoint_shapes <- c(
  "t0" = 16,
  "t24" = 17
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = genotype,
    shape = timepoint
  )
) +

  geom_point(
    size = 4.5
  ) +

  geom_text_repel(
    aes(label = sample_id),
    size = 4.8,
    show.legend = FALSE,
    max.overlaps = Inf,
    seed = 42
  ) +

  scale_color_manual(
    values = genotype_colors,
    name = "Genotype"
  ) +

  scale_shape_manual(
    values = timepoint_shapes,
    name = "Timepoint"
  ) +

  labs(
    title = paste0(
      "Initial quality-control PCA of ",
      "D. purpurea RNA-seq samples"
    ),
    subtitle = "Kallisto transcript abundances: log2(TPM + 1)",
    x = paste0(
      "PC1: ",
      round(explained_variance[1], 1),
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      round(explained_variance[2], 1),
      "% variance"
    )
  ) +

  theme_classic(
    base_size = 18
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 23
    ),
    plot.subtitle = element_text(
      size = 17
    ),
    axis.title = element_text(
      size = 18
    ),
    axis.text = element_text(
      size = 14
    ),
    legend.title = element_text(
      face = "bold",
      size = 18
    ),
    legend.text = element_text(
      size = 16
    )
  )


############################################################
# 12. Save PCA as PNG
############################################################

ggsave(
  filename = png_file,
  plot = pca_plot,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 300,
  bg = "white"
)


############################################################
# 13. Save PCA as SVG
############################################################

grDevices::svg(
  filename = svg_file,
  width = 14,
  height = 10.5,
  pointsize = 12
)

print(pca_plot)

grDevices::dev.off()


############################################################
# 14. Print summary
############################################################

cat(
  "\nInitial Kallisto QC PCA completed successfully.\n\n"
)

cat(
  "Samples:",
  ncol(tpm_matrix),
  "\n"
)

cat(
  "Transcripts before filtering:",
  nrow(tpm_matrix),
  "\n"
)

cat(
  "Transcripts after filtering:",
  nrow(filtered_tpm),
  "\n"
)

cat(
  "Transcripts used for PCA:",
  ncol(pca_input),
  "\n"
)

cat(
  "PC1:",
  round(explained_variance[1], 1),
  "%\n"
)

cat(
  "PC2:",
  round(explained_variance[2], 1),
  "%\n"
)

cat(
  "\nPNG saved to:\n",
  png_file,
  "\n"
)

cat(
  "\nSVG saved to:\n",
  svg_file,
  "\n"
)
