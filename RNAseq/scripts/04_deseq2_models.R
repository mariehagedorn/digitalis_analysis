#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea RNA-seq
# DESeq2 models for differential expression analysis
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
})


############################################################
# 2. Define paths
#
# The script assumes that it is run from the root directory
# of the GitHub repository.
############################################################

base_dir <- "RNAseq"

metadata_dir <- file.path(
  base_dir,
  "metadata"
)

results_dir <- file.path(
  base_dir,
  "results"
)

txi_file <- file.path(
  results_dir,
  "tximport_gene_level.rds"
)

samples_file <- file.path(
  metadata_dir,
  "samples.csv"
)

goi_file <- file.path(
  metadata_dir,
  "genes_of_interest.csv"
)

# Create results directory if it does not already exist
dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Load data
############################################################

txi <- readRDS(
  txi_file
)

samples <- read.csv(
  samples_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

goi <- read.csv(
  goi_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 4. Validate metadata
############################################################

required_columns <- c(
  "sample_id",
  "plant_id",
  "genotype",
  "timepoint"
)

missing_columns <- setdiff(
  required_columns,
  colnames(samples)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in samples.csv: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(samples$sample_id)) {
  stop(
    "Duplicate sample IDs were found in samples.csv."
  )
}

rownames(samples) <- samples$sample_id

if (is.null(colnames(txi$counts))) {
  stop(
    "The tximport count matrix does not contain sample names."
  )
}

if (!setequal(
  colnames(txi$counts),
  rownames(samples)
)) {
  stop(
    "Sample names in the count matrix and metadata do not match."
  )
}

samples <- samples[
  colnames(txi$counts),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(txi$counts),
    rownames(samples)
  )
)


############################################################
# 5. Define experimental factors
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
# 6. Number plants within genotypes
#
# This is required for the paired model because each plant
# was sampled at both t0 and t24.
############################################################

plant_map <- unique(
  samples[
    ,
    c(
      "genotype",
      "plant_id"
    )
  ]
)

plant_map$plant_index <- ave(
  seq_len(nrow(plant_map)),
  plant_map$genotype,
  FUN = seq_along
)

map_key <- paste(
  plant_map$genotype,
  plant_map$plant_id,
  sep = "__"
)

sample_key <- paste(
  samples$genotype,
  samples$plant_id,
  sep = "__"
)

samples$plant_index <- factor(
  plant_map$plant_index[
    match(
      sample_key,
      map_key
    )
  ]
)

if (any(is.na(samples$plant_index))) {
  stop(
    "Plant assignment within genotypes failed."
  )
}

cat(
  "Experimental design:\n"
)

print(
  with(
    samples,
    table(
      genotype,
      timepoint
    )
  )
)

cat(
  "\nPlant assignment:\n"
)

print(
  unique(
    samples[
      ,
      c(
        "plant_id",
        "genotype",
        "plant_index"
      )
    ]
  )
)


############################################################
# 7. Create base DESeq2 dataset
############################################################

dds_base <- DESeqDataSetFromTximport(
  txi = txi,
  colData = samples,
  design = ~ 1
)

cat(
  "\nGenes before expression filtering:",
  nrow(dds_base),
  "\n"
)


############################################################
# 8. Expression filtering
#
# Standard criterion:
# at least 10 estimated counts in at least three samples.
#
# Predefined genes of interest are retained independently
# of whether they pass the standard expression filter.
############################################################

standard_keep <- rowSums(
  counts(dds_base) >= 10
) >= 3

is_goi <- rownames(dds_base) %in% goi$gene_id

goi_filter_status <- goi

goi_filter_status$total_count <- rowSums(
  counts(dds_base)
)[goi$gene_id]

goi_filter_status$passes_standard_filter <- standard_keep[
  goi$gene_id
]

write.csv(
  goi_filter_status,
  file.path(
    results_dir,
    "GOI_filter_status.csv"
  ),
  row.names = FALSE
)

keep <- standard_keep | is_goi

dds_base <- dds_base[
  keep,
]

cat(
  "Genes after expression filtering:",
  nrow(dds_base),
  "\n"
)

cat(
  "Genes of interest retained after filtering:",
  sum(
    rownames(dds_base) %in% goi$gene_id
  ),
  "of",
  nrow(goi),
  "\n"
)


############################################################
# Helper function: validate model matrix
############################################################

check_model_matrix <- function(
  dds,
  model_name
) {

  model_matrix <- model.matrix(
    design(dds),
    data = as.data.frame(
      colData(dds)
    )
  )

  model_rank <- qr(
    model_matrix
  )$rank

  cat(
    "\n",
    model_name,
    ": ",
    ncol(model_matrix),
    " model parameters, rank ",
    model_rank,
    "\n",
    sep = ""
  )

  if (
    model_rank <
      ncol(model_matrix)
  ) {
    stop(
      "The model matrix for ",
      model_name,
      " does not have full rank."
    )
  }
}


############################################################
# 9. Paired model
#
# This model is used to investigate:
# - t24 vs. t0 within each genotype
# - differences in the temporal response between genotypes
#
# Genotype main effects are not interpreted from this model
# because plant identity is nested within genotype.
############################################################

dds_paired <- dds_base

design(dds_paired) <-
  ~ genotype +
  genotype:plant_index +
  genotype:timepoint

check_model_matrix(
  dds_paired,
  "Paired model"
)

dds_paired <- DESeq(
  dds_paired
)

saveRDS(
  dds_paired,
  file.path(
    results_dir,
    "dds_paired_model.rds"
  )
)


############################################################
# 10. Genotype differences at t0
############################################################

dds_t0 <- dds_base[
  ,
  dds_base$timepoint == "t0"
]

dds_t0$genotype <- droplevels(
  dds_t0$genotype
)

design(dds_t0) <- ~ genotype

check_model_matrix(
  dds_t0,
  "t0 model"
)

dds_t0 <- DESeq(
  dds_t0
)

saveRDS(
  dds_t0,
  file.path(
    results_dir,
    "dds_t0_model.rds"
  )
)


############################################################
# 11. Genotype differences at t24
############################################################

dds_t24 <- dds_base[
  ,
  dds_base$timepoint == "t24"
]

dds_t24$genotype <- droplevels(
  dds_t24$genotype
)

design(dds_t24) <- ~ genotype

check_model_matrix(
  dds_t24,
  "t24 model"
)

dds_t24 <- DESeq(
  dds_t24
)

saveRDS(
  dds_t24,
  file.path(
    results_dir,
    "dds_t24_model.rds"
  )
)


############################################################
# 12. Save normalized counts
#
# Normalization is based on all retained genes.
# Gene-of-interest counts are additionally exported
# together with their annotation.
############################################################

normalized_counts <- counts(
  dds_paired,
  normalized = TRUE
)

saveRDS(
  normalized_counts,
  file.path(
    results_dir,
    "normalized_counts_all_genes.rds"
  )
)

goi_normalized_counts <- cbind(
  goi,
  as.data.frame(
    normalized_counts[
      goi$gene_id,
      ,
      drop = FALSE
    ]
  )
)

write.csv(
  goi_normalized_counts,
  file.path(
    results_dir,
    "GOI_normalized_counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    colData(dds_paired)
  ),
  file.path(
    results_dir,
    "model_metadata.csv"
  ),
  row.names = TRUE
)


############################################################
# 13. Print model coefficients
############################################################

cat(
  "\nCoefficients of the paired model:\n"
)

print(
  resultsNames(
    dds_paired
  )
)

cat(
  "\nCoefficients of the t0 model:\n"
)

print(
  resultsNames(
    dds_t0
  )
)

cat(
  "\nCoefficients of the t24 model:\n"
)

print(
  resultsNames(
    dds_t24
  )
)

cat(
  "\nAll three DESeq2 models were calculated successfully.\n"
)
