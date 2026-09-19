#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea RNA-seq
# Gene-of-interest differential expression results
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

results_dir <- file.path(
  base_dir,
  "results"
)

metadata_dir <- file.path(
  base_dir,
  "metadata"
)

# Create results directory if it does not already exist
dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 3. Load DESeq2 models and gene-of-interest information
############################################################

dds_paired <- readRDS(
  file.path(
    results_dir,
    "dds_paired_model.rds"
  )
)

dds_t0 <- readRDS(
  file.path(
    results_dir,
    "dds_t0_model.rds"
  )
)

dds_t24 <- readRDS(
  file.path(
    results_dir,
    "dds_t24_model.rds"
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


############################################################
# 4. Helper function for extracting genes of interest
############################################################

extract_goi <- function(
  res,
  comparison,
  analysis_block
) {

  res_df <- as.data.frame(
    res
  )

  res_df$gene_id <- rownames(
    res_df
  )

  res_goi <- res_df[
    match(
      goi$gene_id,
      res_df$gene_id
    ),
    ,
    drop = FALSE
  ]

  if (!all(
    res_goi$gene_id == goi$gene_id
  )) {
    stop(
      "The order of the gene-of-interest IDs could not be assigned correctly."
    )
  }

  filter_index <- match(
    goi$gene_id,
    goi_filter$gene_id
  )

  output <- data.frame(
    analysis_block = analysis_block,
    comparison = comparison,
    gene_id = goi$gene_id,
    gene_name = goi$gene_name,
    category = goi$category,
    total_count =
      goi_filter$total_count[
        filter_index
      ],
    passes_standard_filter =
      goi_filter$passes_standard_filter[
        filter_index
      ],
    baseMean = res_goi$baseMean,
    log2FoldChange =
      res_goi$log2FoldChange,
    lfcSE = res_goi$lfcSE,
    stat = res_goi$stat,
    pvalue = res_goi$pvalue,
    padj_genomewide = res_goi$padj,
    stringsAsFactors = FALSE
  )

  output$expression_status <- ifelse(
    output$total_count == 0,
    "not_detected",
    ifelse(
      output$passes_standard_filter,
      "passes_filter",
      "low_expression"
    )
  )

  output$effect_direction <- ifelse(
    is.na(
      output$log2FoldChange
    ),
    "not_testable",
    ifelse(
      output$log2FoldChange > 0,
      "positive",
      ifelse(
        output$log2FoldChange < 0,
        "negative",
        "no_change"
      )
    )
  )


  ##########################################################
  # GOI-specific Benjamini-Hochberg correction
  #
  # Correction is performed separately for each comparison
  # across the testable, predefined genes of interest.
  ##########################################################

  output$padj_goi <- NA_real_

  testable <- !is.na(
    output$pvalue
  )

  output$padj_goi[
    testable
  ] <- p.adjust(
    output$pvalue[
      testable
    ],
    method = "BH"
  )

  output$significant_goi_FDR <- (
    !is.na(
      output$padj_goi
    ) &
      output$padj_goi < 0.05
  )

  output$significant_genomewide_FDR <- (
    !is.na(
      output$padj_genomewide
    ) &
      output$padj_genomewide < 0.05
  )

  return(
    output
  )
}


############################################################
# 5. Helper function for genotype-by-time interactions
#
# Positive log2FC:
# The first-listed genotype shows a stronger change
# from t0 to t24 than the second-listed genotype.
############################################################

interaction_result <- function(
  positive_coefficient,
  negative_coefficient
) {

  coefficient_names <- resultsNames(
    dds_paired
  )

  required <- c(
    positive_coefficient,
    negative_coefficient
  )

  if (!all(
    required %in% coefficient_names
  )) {
    stop(
      "At least one interaction coefficient was not found: ",
      paste(
        setdiff(
          required,
          coefficient_names
        ),
        collapse = ", "
      )
    )
  }

  numeric_contrast <- setNames(
    rep(
      0,
      length(
        coefficient_names
      )
    ),
    coefficient_names
  )

  numeric_contrast[
    positive_coefficient
  ] <- 1

  numeric_contrast[
    negative_coefficient
  ] <- -1

  results(
    dds_paired,
    contrast = numeric_contrast,
    alpha = 0.05
  )
}


############################################################
# 6. Time effect: t24 vs. t0 within each genotype
#
# Positive log2FC = higher expression at t24.
############################################################

time_AA <- results(
  dds_paired,
  name = "genotypeAA.timepointt24",
  alpha = 0.05
)

time_Aa <- results(
  dds_paired,
  name = "genotypeAa.timepointt24",
  alpha = 0.05
)

time_aa <- results(
  dds_paired,
  name = "genotypeaa.timepointt24",
  alpha = 0.05
)


############################################################
# 7. Genotype differences at t0
#
# Positive log2FC = higher expression in the first-listed
# genotype.
############################################################

t0_Aa_vs_AA <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "Aa",
    "AA"
  ),
  alpha = 0.05
)

t0_aa_vs_AA <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "aa",
    "AA"
  ),
  alpha = 0.05
)

t0_aa_vs_Aa <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "aa",
    "Aa"
  ),
  alpha = 0.05
)


############################################################
# 8. Genotype differences at t24
############################################################

t24_Aa_vs_AA <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "Aa",
    "AA"
  ),
  alpha = 0.05
)

t24_aa_vs_AA <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "aa",
    "AA"
  ),
  alpha = 0.05
)

t24_aa_vs_Aa <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "aa",
    "Aa"
  ),
  alpha = 0.05
)


############################################################
# 9. Genotype-by-time interactions
#
# Example:
# interaction_Aa_vs_AA =
# temporal change in Aa minus temporal change in AA.
############################################################

interaction_Aa_vs_AA <- interaction_result(
  "genotypeAa.timepointt24",
  "genotypeAA.timepointt24"
)

interaction_aa_vs_AA <- interaction_result(
  "genotypeaa.timepointt24",
  "genotypeAA.timepointt24"
)

interaction_aa_vs_Aa <- interaction_result(
  "genotypeaa.timepointt24",
  "genotypeAa.timepointt24"
)


############################################################
# 10. Combine results for all genes of interest
############################################################

all_results <- rbind(

  extract_goi(
    time_AA,
    "AA_t24_vs_t0",
    "time_effect"
  ),

  extract_goi(
    time_Aa,
    "Aa_t24_vs_t0",
    "time_effect"
  ),

  extract_goi(
    time_aa,
    "aa_t24_vs_t0",
    "time_effect"
  ),

  extract_goi(
    t0_Aa_vs_AA,
    "t0_Aa_vs_AA",
    "genotype_at_t0"
  ),

  extract_goi(
    t0_aa_vs_AA,
    "t0_aa_vs_AA",
    "genotype_at_t0"
  ),

  extract_goi(
    t0_aa_vs_Aa,
    "t0_aa_vs_Aa",
    "genotype_at_t0"
  ),

  extract_goi(
    t24_Aa_vs_AA,
    "t24_Aa_vs_AA",
    "genotype_at_t24"
  ),

  extract_goi(
    t24_aa_vs_AA,
    "t24_aa_vs_AA",
    "genotype_at_t24"
  ),

  extract_goi(
    t24_aa_vs_Aa,
    "t24_aa_vs_Aa",
    "genotype_at_t24"
  ),

  extract_goi(
    interaction_Aa_vs_AA,
    "interaction_Aa_vs_AA",
    "interaction"
  ),

  extract_goi(
    interaction_aa_vs_AA,
    "interaction_aa_vs_AA",
    "interaction"
  ),

  extract_goi(
    interaction_aa_vs_Aa,
    "interaction_aa_vs_Aa",
    "interaction"
  )
)


############################################################
# 11. Save result tables
############################################################

write.csv(
  all_results,
  file.path(
    results_dir,
    "GOI_all_contrasts.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results[
    all_results$analysis_block ==
      "time_effect",
  ],
  file.path(
    results_dir,
    "GOI_time_effects.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results[
    all_results$analysis_block ==
      "genotype_at_t0",
  ],
  file.path(
    results_dir,
    "GOI_genotype_comparisons_t0.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results[
    all_results$analysis_block ==
      "genotype_at_t24",
  ],
  file.path(
    results_dir,
    "GOI_genotype_comparisons_t24.csv"
  ),
  row.names = FALSE
)

write.csv(
  all_results[
    all_results$analysis_block ==
      "interaction",
  ],
  file.path(
    results_dir,
    "GOI_interactions.csv"
  ),
  row.names = FALSE
)

significant_results <- all_results[
  all_results$significant_goi_FDR |
    all_results$significant_genomewide_FDR,
]

write.csv(
  significant_results,
  file.path(
    results_dir,
    "GOI_significant_results.csv"
  ),
  row.names = FALSE
)


############################################################
# 12. Create summary table for each comparison
############################################################

comparison_order <- unique(
  all_results$comparison
)

summary_table <- do.call(
  rbind,
  lapply(
    comparison_order,
    function(
      current_comparison
    ) {

      current <- all_results[
        all_results$comparison ==
          current_comparison,
      ]

      data.frame(
        comparison =
          current_comparison,

        testable_genes =
          sum(
            !is.na(
              current$pvalue
            )
          ),

        nominal_p_below_0.05 =
          sum(
            current$pvalue < 0.05,
            na.rm = TRUE
          ),

        significant_goi_FDR =
          sum(
            current$significant_goi_FDR,
            na.rm = TRUE
          ),

        significant_genomewide_FDR =
          sum(
            current$significant_genomewide_FDR,
            na.rm = TRUE
          )
      )
    }
  )
)

rownames(
  summary_table
) <- NULL

write.csv(
  summary_table,
  file.path(
    results_dir,
    "GOI_results_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 13. Completion message
############################################################

cat(
  "\nSummary of gene-of-interest analysis:\n\n"
)

print(
  summary_table,
  row.names = FALSE
)

cat(
  "\nGene-of-interest results were saved successfully.\n"
)

cat(
  "MYB75b has no counts and is therefore not statistically testable.\n"
)
