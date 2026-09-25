#!/usr/bin/env Rscript

# ============================================================
# 10_flavonoid_candidate_check.R
#
# Verify presence and expression of the predefined
# flavonoid-related candidate panel using gene-level
# tximport counts and TPM estimates.
#
# The standard expression filter used in this project is:
#   >= 10 estimated counts in at least 3 samples
#
# Expected project structure:
#   <project>/deseq2/results/tximport_gene_level.rds
#   <project>/deseq2/results/dds_paired_model.rds
#   <project>/deseq2/metadata/full_flavonoid_candidate_panel_v3.csv
#
# By default:
#   <project> = /vol/data/digitalis_rnaseq
#
# To use a different project directory, set:
#   DIGITALIS_RNASEQ_DIR=/path/to/project
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
})


# ------------------------------------------------------------
# 1. Paths and input files
# ------------------------------------------------------------

project_dir <- Sys.getenv(
  "DIGITALIS_RNASEQ_DIR",
  unset = "/vol/data/digitalis_rnaseq"
)

deseq2_dir <- file.path(
  project_dir,
  "deseq2"
)

results_dir <- file.path(
  deseq2_dir,
  "results"
)

metadata_dir <- file.path(
  deseq2_dir,
  "metadata"
)

candidate_file_options <- c(
  file.path(
    metadata_dir,
    "full_flavonoid_candidate_panel_v3.csv"
  ),
  file.path(
    getwd(),
    "full_flavonoid_candidate_panel_v3.csv"
  )
)

candidate_file <- candidate_file_options[
  file.exists(candidate_file_options)
][1]

if (
  length(candidate_file) == 0 ||
    is.na(candidate_file)
) {
  stop(
    paste0(
      "Could not find full_flavonoid_candidate_panel_v3.csv. ",
      "Expected it in ",
      metadata_dir,
      " or in the current working directory."
    )
  )
}

tximport_file <- file.path(
  results_dir,
  "tximport_gene_level.rds"
)

dds_file <- file.path(
  results_dir,
  "dds_paired_model.rds"
)

if (!file.exists(tximport_file)) {
  stop(
    "Missing input file: ",
    tximport_file
  )
}

if (!file.exists(dds_file)) {
  stop(
    "Missing input file: ",
    dds_file
  )
}

txi <- readRDS(
  tximport_file
)

dds <- readRDS(
  dds_file
)

panel <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Validate candidate panel and tximport object
# ------------------------------------------------------------

required_columns <- c(
  "pathway_section",
  "enzyme",
  "gene_name",
  "transcript_id",
  "gene_id",
  "category",
  "source_status"
)

missing_columns <- setdiff(
  required_columns,
  colnames(panel)
)

if (length(missing_columns) > 0) {
  stop(
    "Candidate panel is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

if (anyDuplicated(panel$gene_id)) {
  stop(
    "Candidate panel contains duplicated DESeq2 gene IDs."
  )
}

if (anyDuplicated(panel$gene_name)) {
  stop(
    "Candidate panel contains duplicated candidate gene names."
  )
}

if (
  is.null(txi$counts) ||
    is.null(txi$abundance)
) {
  stop(
    paste0(
      "The tximport object must contain gene-level ",
      "estimated counts and TPM abundance values."
    )
  )
}


# ------------------------------------------------------------
# 3. Match samples and metadata
# ------------------------------------------------------------

common_samples <- Reduce(
  intersect,
  list(
    colnames(txi$counts),
    colnames(txi$abundance),
    colnames(dds)
  )
)

if (length(common_samples) != 18) {
  stop(
    "Expected 18 common RNA-seq samples, found ",
    length(common_samples),
    "."
  )
}

metadata <- as.data.frame(
  colData(dds)[
    common_samples,
    ,
    drop = FALSE
  ]
)

metadata$sample_id <- rownames(
  metadata
)

if (!"plant_id" %in% colnames(metadata)) {
  metadata$plant_id <- sub(
    "_(t0|t24)$",
    "",
    metadata$sample_id
  )
}

required_metadata <- c(
  "genotype",
  "timepoint",
  "plant_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(metadata)
)

if (length(missing_metadata) > 0) {
  stop(
    "DESeq2 metadata are missing required columns: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. Candidate expression status
#
# Standard filter:
# >= 10 estimated counts in at least 3 samples
# ------------------------------------------------------------

counts_matrix <- txi$counts[
  ,
  common_samples,
  drop = FALSE
]

tpm_matrix <- txi$abundance[
  ,
  common_samples,
  drop = FALSE
]

candidate_summary <- panel

candidate_summary$present_in_tximport <-
  panel$gene_id %in%
  rownames(counts_matrix)

candidate_summary$total_counts <- NA_real_
candidate_summary$mean_TPM <- NA_real_
candidate_summary$max_TPM <- NA_real_
candidate_summary$samples_TPM_above_1 <- NA_integer_
candidate_summary$samples_counts_at_least_10 <- NA_integer_
candidate_summary$passes_standard_filter <- FALSE
candidate_summary$expression_status <- "missing_from_tximport"


for (i in seq_len(nrow(candidate_summary))) {

  gene_id <- candidate_summary$gene_id[i]

  if (
    !candidate_summary$present_in_tximport[i]
  ) {
    next
  }

  gene_counts <- as.numeric(
    counts_matrix[
      gene_id,
    ]
  )

  gene_tpm <- as.numeric(
    tpm_matrix[
      gene_id,
    ]
  )

  candidate_summary$total_counts[i] <- sum(
    gene_counts,
    na.rm = TRUE
  )

  candidate_summary$mean_TPM[i] <- mean(
    gene_tpm,
    na.rm = TRUE
  )

  candidate_summary$max_TPM[i] <- max(
    gene_tpm,
    na.rm = TRUE
  )

  candidate_summary$samples_TPM_above_1[i] <- sum(
    gene_tpm > 1,
    na.rm = TRUE
  )

  candidate_summary$samples_counts_at_least_10[i] <- sum(
    gene_counts >= 10,
    na.rm = TRUE
  )

  candidate_summary$passes_standard_filter[i] <-
    candidate_summary$samples_counts_at_least_10[i] >= 3

  candidate_summary$expression_status[i] <-
    if (candidate_summary$total_counts[i] == 0) {
      "not_detected"
    } else if (
      candidate_summary$passes_standard_filter[i]
    ) {
      "passes_filter"
    } else {
      "low_expression"
    }
}


# ------------------------------------------------------------
# 5. TPM and count values for each sample
# ------------------------------------------------------------

per_sample_list <- vector(
  mode = "list",
  length = nrow(panel)
)

for (i in seq_len(nrow(panel))) {

  gene_id <- panel$gene_id[i]

  if (
    !gene_id %in%
      rownames(counts_matrix)
  ) {
    next
  }

  per_sample_list[[i]] <- data.frame(
    pathway_section = panel$pathway_section[i],
    enzyme = panel$enzyme[i],
    gene_name = panel$gene_name[i],
    gene_id = gene_id,
    sample_id = common_samples,
    genotype = as.character(
      metadata$genotype
    ),
    timepoint = as.character(
      metadata$timepoint
    ),
    plant_id = as.character(
      metadata$plant_id
    ),
    estimated_counts = as.numeric(
      counts_matrix[
        gene_id,
      ]
    ),
    TPM = as.numeric(
      tpm_matrix[
        gene_id,
      ]
    ),
    stringsAsFactors = FALSE
  )
}

non_empty_entries <- !vapply(
  per_sample_list,
  is.null,
  logical(1)
)

if (!any(non_empty_entries)) {
  stop(
    "None of the candidate genes were found in the tximport object."
  )
}

per_sample <- do.call(
  rbind,
  per_sample_list[
    non_empty_entries
  ]
)


# ------------------------------------------------------------
# 6. Mean expression by genotype and time point
# ------------------------------------------------------------

group_means <- aggregate(
  cbind(
    estimated_counts,
    TPM
  ) ~
    pathway_section +
    enzyme +
    gene_name +
    gene_id +
    genotype +
    timepoint,
  data = per_sample,
  FUN = mean
)


# ------------------------------------------------------------
# 7. Export results
# ------------------------------------------------------------

candidate_summary_file <- file.path(
  results_dir,
  "full_flavonoid_candidate_expression_check_v3.csv"
)

per_sample_file <- file.path(
  results_dir,
  "full_flavonoid_TPM_counts_per_sample_v3.csv"
)

group_means_file <- file.path(
  results_dir,
  "full_flavonoid_group_means_v3.csv"
)

write.csv(
  candidate_summary,
  candidate_summary_file,
  row.names = FALSE
)

write.csv(
  per_sample,
  per_sample_file,
  row.names = FALSE
)

write.csv(
  group_means,
  group_means_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------

cat(
  "\nFlavonoid candidate expression check completed.\n\n"
)

cat(
  "Candidate panel: ",
  nrow(panel),
  " genes\n",
  sep = ""
)

cat(
  "Present in tximport: ",
  sum(
    candidate_summary$present_in_tximport
  ),
  "\n",
  sep = ""
)

cat(
  "Passed standard filter: ",
  sum(
    candidate_summary$passes_standard_filter
  ),
  "\n\n",
  sep = ""
)

cat(
  "Expression-status summary:\n"
)

print(
  table(
    candidate_summary$expression_status
  )
)

cat(
  "\nOutput files:\n",
  "- ",
  candidate_summary_file,
  "\n- ",
  per_sample_file,
  "\n- ",
  group_means_file,
  "\n",
  sep = ""
)
