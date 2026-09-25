#!/usr/bin/env Rscript

# ============================================================
# 11_persistent_genotype_effects.R
#
# Identify persistent aa-associated genotype effects across
# t0 and t24, perform apeglm shrinkage, rank the stringent
# candidate set, and export the full set plus the Top 10.
#
# Expected input:
#   <project>/deseq2/results/dds_t0_model.rds
#   <project>/deseq2/results/dds_t24_model.rds
#
# By default, <project> is /vol/data/digitalis_rnaseq.
# This can be changed by setting the environment variable:
#   DIGITALIS_RNASEQ_DIR=/path/to/project
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_dir <- Sys.getenv(
  "DIGITALIS_RNASEQ_DIR",
  unset = "/vol/data/digitalis_rnaseq"
)

results_dir <- file.path(
  project_dir,
  "deseq2",
  "results"
)

outdir <- file.path(
  results_dir,
  "genotype_background_analysis"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dds_t0_file <- file.path(
  results_dir,
  "dds_t0_model.rds"
)

dds_t24_file <- file.path(
  results_dir,
  "dds_t24_model.rds"
)

if (!file.exists(dds_t0_file)) {
  stop("Missing input file: ", dds_t0_file)
}

if (!file.exists(dds_t24_file)) {
  stop("Missing input file: ", dds_t24_file)
}


# ------------------------------------------------------------
# 2. Load DESeq2 objects
# ------------------------------------------------------------

dds_t0 <- readRDS(dds_t0_file)
dds_t24 <- readRDS(dds_t24_file)

required_names <- c(
  "Intercept",
  "genotype_Aa_vs_AA",
  "genotype_aa_vs_AA"
)

if (!all(required_names %in% resultsNames(dds_t0))) {
  stop(
    "Unexpected coefficient names in dds_t0. Found: ",
    paste(resultsNames(dds_t0), collapse = ", ")
  )
}

if (!all(required_names %in% resultsNames(dds_t24))) {
  stop(
    "Unexpected coefficient names in dds_t24. Found: ",
    paste(resultsNames(dds_t24), collapse = ", ")
  )
}


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

result_to_df <- function(res, prefix) {

  data.frame(
    gene_id = rownames(res),
    setNames(
      data.frame(
        baseMean = res$baseMean,
        log2FC = res$log2FoldChange,
        lfcSE = res$lfcSE,
        pvalue = res$pvalue,
        padj = res$padj,
        check.names = FALSE
      ),
      paste0(
        c(
          "baseMean_",
          "log2FC_",
          "lfcSE_",
          "pvalue_",
          "padj_"
        ),
        prefix
      )
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


is_sig <- function(x, alpha = 0.05) {
  !is.na(x) & x < alpha
}


shrink_direct_aa_vs_AA <- function(dds) {

  coef_name <- "genotype_aa_vs_AA"

  if (!coef_name %in% resultsNames(dds)) {
    stop("Coefficient not found: ", coef_name)
  }

  lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
}


shrink_aa_vs_Aa <- function(dds) {

  dds_relevel <- dds

  dds_relevel$genotype <- relevel(
    droplevels(dds_relevel$genotype),
    ref = "Aa"
  )

  # Dispersions are retained; only the Wald model is refitted
  # to obtain aa vs Aa as a named coefficient for apeglm.
  dds_relevel <- nbinomWaldTest(dds_relevel)

  coef_name <- "genotype_aa_vs_Aa"

  if (!coef_name %in% resultsNames(dds_relevel)) {
    stop(
      "Coefficient ",
      coef_name,
      " not found after releveling. Found: ",
      paste(resultsNames(dds_relevel), collapse = ", ")
    )
  }

  lfcShrink(
    dds_relevel,
    coef = coef_name,
    type = "apeglm"
  )
}


# ------------------------------------------------------------
# 4. Genome-wide genotype contrasts at t0 and t24
# ------------------------------------------------------------

res_t0_Aa_vs_AA <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "Aa",
    "AA"
  ),
  alpha = 0.05
)

res_t0_aa_vs_AA <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "aa",
    "AA"
  ),
  alpha = 0.05
)

res_t0_aa_vs_Aa <- results(
  dds_t0,
  contrast = c(
    "genotype",
    "aa",
    "Aa"
  ),
  alpha = 0.05
)

res_t24_Aa_vs_AA <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "Aa",
    "AA"
  ),
  alpha = 0.05
)

res_t24_aa_vs_AA <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "aa",
    "AA"
  ),
  alpha = 0.05
)

res_t24_aa_vs_Aa <- results(
  dds_t24,
  contrast = c(
    "genotype",
    "aa",
    "Aa"
  ),
  alpha = 0.05
)


# ------------------------------------------------------------
# 5. Merge contrasts into one gene-level table
# ------------------------------------------------------------

comparison_tables <- list(
  result_to_df(
    res_t0_Aa_vs_AA,
    "Aa_vs_AA_t0"
  ),
  result_to_df(
    res_t24_Aa_vs_AA,
    "Aa_vs_AA_t24"
  ),
  result_to_df(
    res_t0_aa_vs_AA,
    "aa_vs_AA_t0"
  ),
  result_to_df(
    res_t24_aa_vs_AA,
    "aa_vs_AA_t24"
  ),
  result_to_df(
    res_t0_aa_vs_Aa,
    "aa_vs_Aa_t0"
  ),
  result_to_df(
    res_t24_aa_vs_Aa,
    "aa_vs_Aa_t24"
  )
)

all_genotype_results <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "gene_id",
      all = TRUE,
      sort = FALSE
    )
  },
  comparison_tables
)


# ------------------------------------------------------------
# 6. Define persistent aa-associated genes
#
# Initial set:
# - aa vs AA significant at t0 and t24
# - aa vs Aa significant at t0 and t24
# - same log2FC direction in all four comparisons
#
# Stringent set:
# - additionally, no significant AA vs Aa difference at t0
#   or t24
# ------------------------------------------------------------

sig_four <- with(
  all_genotype_results,
  is_sig(padj_aa_vs_AA_t0) &
    is_sig(padj_aa_vs_AA_t24) &
    is_sig(padj_aa_vs_Aa_t0) &
    is_sig(padj_aa_vs_Aa_t24)
)

same_direction <- with(
  all_genotype_results,
  (
    log2FC_aa_vs_AA_t0 > 0 &
      log2FC_aa_vs_AA_t24 > 0 &
      log2FC_aa_vs_Aa_t0 > 0 &
      log2FC_aa_vs_Aa_t24 > 0
  ) |
    (
      log2FC_aa_vs_AA_t0 < 0 &
        log2FC_aa_vs_AA_t24 < 0 &
        log2FC_aa_vs_Aa_t0 < 0 &
        log2FC_aa_vs_Aa_t24 < 0
    )
)

initial_persistent <- all_genotype_results[
  sig_four & same_direction,
]

red_genotype_difference <- with(
  initial_persistent,
  is_sig(padj_Aa_vs_AA_t0) |
    is_sig(padj_Aa_vs_AA_t24)
)

strict_aa <- initial_persistent[
  !red_genotype_difference,
]


# ------------------------------------------------------------
# 7. Raw-effect score
#
# This is retained only to preserve a stable table/figure
# numbering order. Biological prioritization below is based on
# apeglm-shrunken effect sizes.
# ------------------------------------------------------------

strict_aa$min_abs_rawLFC_all4 <- with(
  strict_aa,
  pmin(
    abs(log2FC_aa_vs_AA_t0),
    abs(log2FC_aa_vs_AA_t24),
    abs(log2FC_aa_vs_Aa_t0),
    abs(log2FC_aa_vs_Aa_t24)
  )
)

strict_aa <- strict_aa[
  order(
    strict_aa$min_abs_rawLFC_all4,
    decreasing = TRUE
  ),
]

strict_aa$raw_rank <- seq_len(
  nrow(strict_aa)
)

# Stable lookup numbering used by the supplementary figure
strict_aa$table_no <- seq_len(
  nrow(strict_aa)
)

strict_aa$table_no_fmt <- sprintf(
  "%02d",
  strict_aa$table_no
)


# ------------------------------------------------------------
# 8. apeglm shrinkage
# ------------------------------------------------------------

shr_t0_aa_vs_AA <- shrink_direct_aa_vs_AA(
  dds_t0
)

shr_t24_aa_vs_AA <- shrink_direct_aa_vs_AA(
  dds_t24
)

shr_t0_aa_vs_Aa <- shrink_aa_vs_Aa(
  dds_t0
)

shr_t24_aa_vs_Aa <- shrink_aa_vs_Aa(
  dds_t24
)


get_shrunk_lfc <- function(res, ids) {

  res$log2FoldChange[
    match(
      ids,
      rownames(res)
    )
  ]
}


strict_aa$shrLFC_t0_vs_AA <- get_shrunk_lfc(
  shr_t0_aa_vs_AA,
  strict_aa$gene_id
)

strict_aa$shrLFC_t24_vs_AA <- get_shrunk_lfc(
  shr_t24_aa_vs_AA,
  strict_aa$gene_id
)

strict_aa$shrLFC_t0_vs_Aa <- get_shrunk_lfc(
  shr_t0_aa_vs_Aa,
  strict_aa$gene_id
)

strict_aa$shrLFC_t24_vs_Aa <- get_shrunk_lfc(
  shr_t24_aa_vs_Aa,
  strict_aa$gene_id
)

if (
  any(
    is.na(
      strict_aa[
        ,
        c(
          "shrLFC_t0_vs_AA",
          "shrLFC_t24_vs_AA",
          "shrLFC_t0_vs_Aa",
          "shrLFC_t24_vs_Aa"
        )
      ]
    )
  )
) {
  warning(
    "At least one stringent gene has an NA shrunken log2FC."
  )
}


# ------------------------------------------------------------
# 9. Conservative shrinkage-based ranking
#
# Ranking score = smallest absolute shrunken log2FC across
# all four aa-vs-red-genotype comparisons.
# ------------------------------------------------------------

strict_aa$min_abs_shrLFC_all4 <- with(
  strict_aa,
  pmin(
    abs(shrLFC_t0_vs_AA),
    abs(shrLFC_t24_vs_AA),
    abs(shrLFC_t0_vs_Aa),
    abs(shrLFC_t24_vs_Aa)
  )
)

strict_aa$shrink_rank <- rank(
  -strict_aa$min_abs_shrLFC_all4,
  ties.method = "first"
)

top10 <- strict_aa[
  order(
    strict_aa$shrink_rank
  ),
]

top10 <- head(
  top10,
  10
)


# ------------------------------------------------------------
# 10. Export screening summary
# ------------------------------------------------------------

screening_summary <- data.frame(
  step = c(
    "Significant aa vs AA and aa vs Aa at t0 and t24, same direction",
    "Stringent set after excluding genes significant between AA and Aa",
    "Top genes retained for main-text table"
  ),
  n_genes = c(
    nrow(initial_persistent),
    nrow(strict_aa),
    nrow(top10)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  screening_summary,
  file.path(
    outdir,
    "persistent_genotype_screening_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Export full stringent set and Top 10
# ------------------------------------------------------------

saveRDS(
  strict_aa,
  file.path(
    outdir,
    "strict_aa_signature_57_genes_FINAL_with_shrinkage.rds"
  )
)

write.csv(
  strict_aa,
  file.path(
    outdir,
    "strict_aa_signature_57_genes_FINAL_with_shrinkage.csv"
  ),
  row.names = FALSE
)

top10_export <- top10[
  ,
  c(
    "shrink_rank",
    "gene_id",
    "shrLFC_t0_vs_AA",
    "shrLFC_t24_vs_AA",
    "shrLFC_t0_vs_Aa",
    "shrLFC_t24_vs_Aa",
    "min_abs_shrLFC_all4"
  )
]

names(top10_export)[1] <- "rank"

write.csv(
  top10_export,
  file.path(
    outdir,
    "top10_persistent_aa_genotype_effects_shrunken.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. ANS diagnostic output
#
# ANS is not forced into the stringent set. This block exports
# its genotype-comparison statistics so that its exclusion can
# be traced transparently.
# ------------------------------------------------------------

ans_id <- "DP129852"

if (
  ans_id %in% rownames(res_t0_aa_vs_AA)
) {

  ans_summary <- data.frame(
    gene_id = ans_id,
    comparison = c(
      "aa_vs_AA_t0",
      "aa_vs_AA_t24",
      "aa_vs_Aa_t0",
      "aa_vs_Aa_t24",
      "Aa_vs_AA_t0",
      "Aa_vs_AA_t24"
    ),
    log2FC = c(
      res_t0_aa_vs_AA[ans_id, "log2FoldChange"],
      res_t24_aa_vs_AA[ans_id, "log2FoldChange"],
      res_t0_aa_vs_Aa[ans_id, "log2FoldChange"],
      res_t24_aa_vs_Aa[ans_id, "log2FoldChange"],
      res_t0_Aa_vs_AA[ans_id, "log2FoldChange"],
      res_t24_Aa_vs_AA[ans_id, "log2FoldChange"]
    ),
    pvalue = c(
      res_t0_aa_vs_AA[ans_id, "pvalue"],
      res_t24_aa_vs_AA[ans_id, "pvalue"],
      res_t0_aa_vs_Aa[ans_id, "pvalue"],
      res_t24_aa_vs_Aa[ans_id, "pvalue"],
      res_t0_Aa_vs_AA[ans_id, "pvalue"],
      res_t24_Aa_vs_AA[ans_id, "pvalue"]
    ),
    padj = c(
      res_t0_aa_vs_AA[ans_id, "padj"],
      res_t24_aa_vs_AA[ans_id, "padj"],
      res_t0_aa_vs_Aa[ans_id, "padj"],
      res_t24_aa_vs_Aa[ans_id, "padj"],
      res_t0_Aa_vs_AA[ans_id, "padj"],
      res_t24_Aa_vs_AA[ans_id, "padj"]
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    ans_summary,
    file.path(
      outdir,
      "ANS_genotype_comparison_diagnostics.csv"
    ),
    row.names = FALSE
  )

  if (
    "cooks" %in% assayNames(dds_t0) &&
      ans_id %in% rownames(dds_t0)
  ) {

    ans_cooks <- data.frame(
      sample = colnames(dds_t0),
      genotype = as.character(
        colData(dds_t0)$genotype
      ),
      normalized_count = as.numeric(
        counts(
          dds_t0,
          normalized = TRUE
        )[ans_id, ]
      ),
      cooks_distance = as.numeric(
        assays(dds_t0)[["cooks"]][ans_id, ]
      ),
      stringsAsFactors = FALSE
    )

    write.csv(
      ans_cooks,
      file.path(
        outdir,
        "ANS_t0_counts_and_cooks_distance.csv"
      ),
      row.names = FALSE
    )
  }
}


# ------------------------------------------------------------
# 13. Console summary
# ------------------------------------------------------------

cat(
  "\nPersistent genotype analysis completed.\n",
  "Initial persistent set: ",
  nrow(initial_persistent),
  " genes\n",
  "Stringent aa-associated set: ",
  nrow(strict_aa),
  " genes\n",
  "Top 10 table exported: ",
  nrow(top10),
  " genes\n",
  "Output directory: ",
  outdir,
  "\n",
  sep = ""
)
