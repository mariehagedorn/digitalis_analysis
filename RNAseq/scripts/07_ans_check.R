#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea RNA-seq
# ANS abundance check
#
# Compares ANS abundance estimates from the original
# Kallisto output with gene-level DESeq2 counts.
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

# Create results directory if it does not already exist
dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Load sample metadata and DESeq2 object
############################################################

samples <- read.csv(
  file.path(
    metadata_dir,
    "samples.csv"
  ),
  stringsAsFactors = FALSE
)

dds <- readRDS(
  file.path(
    results_dir,
    "dds_paired_model.rds"
  )
)

ans_gene_id <- "DP129852"


############################################################
# 4. Read ANS abundance directly from each Kallisto file
############################################################

read_ans <- function(
  kallisto_file
) {

  abundance <- read.delim(
    gzfile(
      kallisto_file
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_columns <- c(
    "target_id",
    "est_counts",
    "tpm"
  )

  if (!all(
    required_columns %in%
      colnames(
        abundance
      )
  )) {
    stop(
      "Required columns are missing in: ",
      kallisto_file
    )
  }

  transcript_gene_id <- sub(
    "\\.[0-9]+$",
    "",
    abundance$target_id
  )

  ans_rows <- abundance[
    transcript_gene_id ==
      ans_gene_id,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      ans_rows
    ) == 0
  ) {

    return(
      data.frame(
        detected_transcripts = 0,
        transcript_ids = NA_character_,
        kallisto_est_counts = 0,
        kallisto_tpm = 0
      )
    )
  }

  data.frame(
    detected_transcripts =
      nrow(
        ans_rows
      ),

    transcript_ids =
      paste(
        ans_rows$target_id,
        collapse = ";"
      ),

    kallisto_est_counts =
      sum(
        ans_rows$est_counts,
        na.rm = TRUE
      ),

    kallisto_tpm =
      sum(
        ans_rows$tpm,
        na.rm = TRUE
      )
  )
}


############################################################
# 5. Extract ANS values from all Kallisto samples
############################################################

ans_kallisto <- do.call(
  rbind,
  lapply(
    samples$kallisto_file,
    read_ans
  )
)


############################################################
# 6. Compare with gene-level DESeq2 counts
############################################################

gene_counts <- counts(
  dds,
  normalized = FALSE
)

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

if (!ans_gene_id %in%
    rownames(
      gene_counts
    )) {
  stop(
    "ANS was not found in the DESeq2 object."
  )
}

sample_index <- match(
  samples$sample_id,
  colnames(
    gene_counts
  )
)

if (
  any(
    is.na(
      sample_index
    )
  )
) {
  stop(
    "At least one sample could not be matched."
  )
}


############################################################
# 7. Combine ANS values
############################################################

ans_results <- data.frame(

  sample_id =
    samples$sample_id,

  plant_id =
    samples$plant_id,

  genotype =
    samples$genotype,

  timepoint =
    samples$timepoint,

  ans_kallisto,

  deseq_gene_count =
    gene_counts[
      ans_gene_id,
      sample_index
    ],

  normalized_gene_count =
    normalized_counts[
      ans_gene_id,
      sample_index
    ],

  stringsAsFactors = FALSE
)

ans_results$log2_normalized_count_plus1 <- log2(
  ans_results$normalized_gene_count + 1
)

ans_results$genotype <- factor(
  ans_results$genotype,
  levels = c(
    "AA",
    "Aa",
    "aa"
  )
)

ans_results$timepoint <- factor(
  ans_results$timepoint,
  levels = c(
    "t0",
    "t24"
  )
)

ans_results <- ans_results[
  order(
    ans_results$genotype,
    ans_results$plant_id,
    ans_results$timepoint
  ),
]


############################################################
# 8. Save results
############################################################

write.csv(
  ans_results,
  file.path(
    results_dir,
    "ANS_sample_level_check.csv"
  ),
  row.names = FALSE
)


############################################################
# 9. Prepare concise console output
############################################################

print_output <- ans_results[
  ,
  c(
    "sample_id",
    "genotype",
    "timepoint",
    "transcript_ids",
    "kallisto_est_counts",
    "kallisto_tpm",
    "normalized_gene_count"
  )
]

print_output$kallisto_est_counts <- round(
  print_output$kallisto_est_counts,
  3
)

print_output$kallisto_tpm <- round(
  print_output$kallisto_tpm,
  3
)

print_output$normalized_gene_count <- round(
  print_output$normalized_gene_count,
  3
)


############################################################
# 10. Completion message
############################################################

cat(
  "\nANS values per sample:\n\n"
)

print(
  print_output,
  row.names = FALSE
)

cat(
  "\nANS control table saved as ANS_sample_level_check.csv.\n"
)
