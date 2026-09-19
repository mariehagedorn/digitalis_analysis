#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea RNA-seq
# Global genotype-by-time interaction test
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
# 3. Load DESeq2 model and gene-of-interest information
############################################################

dds_lrt <- readRDS(
  file.path(
    results_dir,
    "dds_paired_model.rds"
  )
)

goi <- read.csv(
  file.path(
    metadata_dir,
    "genes_of_interest.csv"
  ),
  stringsAsFactors = FALSE
)


############################################################
# 4. Global genotype-by-time interaction test
#
# Full model:
# genotype-specific temporal responses
#
# Reduced model:
# one common temporal response across genotypes
############################################################

reduced_design <-
  ~ genotype +
  genotype:plant_index +
  timepoint

dds_lrt <- DESeq(
  dds_lrt,
  test = "LRT",
  reduced = reduced_design
)

saveRDS(
  dds_lrt,
  file.path(
    results_dir,
    "dds_global_interaction_LRT.rds"
  )
)


############################################################
# 5. Extract LRT results
############################################################

lrt_results <- results(
  dds_lrt,
  alpha = 0.05
)

lrt_df <- as.data.frame(
  lrt_results
)

lrt_df$gene_id <- rownames(
  lrt_df
)

goi_lrt <- lrt_df[
  match(
    goi$gene_id,
    lrt_df$gene_id
  ),
  ,
  drop = FALSE
]

output <- data.frame(
  gene_id = goi$gene_id,
  gene_name = goi$gene_name,
  category = goi$category,
  baseMean = goi_lrt$baseMean,
  LRT_statistic = goi_lrt$stat,
  pvalue = goi_lrt$pvalue,
  padj_genomewide = goi_lrt$padj,
  stringsAsFactors = FALSE
)


############################################################
# 6. GOI-specific Benjamini-Hochberg correction
############################################################

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


############################################################
# 7. Save results
############################################################

write.csv(
  output,
  file.path(
    results_dir,
    "GOI_global_interaction_LRT.csv"
  ),
  row.names = FALSE
)

significant_output <- output[
  output$significant_goi_FDR |
    output$significant_genomewide_FDR,
]


############################################################
# 8. Print summary
############################################################

cat(
  "\nGlobal genotype-by-time interaction test:\n\n"
)

cat(
  "Testable genes of interest:",
  sum(
    testable
  ),
  "\n"
)

cat(
  "Significant at GOI-specific FDR:",
  sum(
    output$significant_goi_FDR
  ),
  "\n"
)

cat(
  "Significant at genome-wide FDR:",
  sum(
    output$significant_genomewide_FDR
  ),
  "\n\n"
)

if (
  nrow(
    significant_output
  ) == 0
) {

  cat(
    "No gene of interest showed a significant global interaction.\n"
  )

} else {

  print(
    significant_output,
    row.names = FALSE
  )
}

cat(
  "\nGlobal interaction analysis completed successfully.\n"
)
