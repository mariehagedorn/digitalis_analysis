# Import Kallisto abundance estimates and summarize expression at gene level
############################################################
# Digitalis purpurea RNA-seq
# Import Kallisto estimates with tximport
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(tximport)
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

sample_file <- file.path(
  metadata_dir,
  "samples.csv"
)

# Create directories if they do not already exist
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


############################################################
# 3. Import sample metadata
############################################################

samples <- read_csv(
  sample_file,
  show_col_types = FALSE
)

samples <- as.data.frame(samples)

stopifnot(
  nrow(samples) == 18
)

stopifnot(
  !anyDuplicated(samples$sample_id)
)

stopifnot(
  all(file.exists(samples$kallisto_file))
)


############################################################
# 4. Define experimental factors
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

rownames(samples) <- samples$sample_id


############################################################
# 5. Define Kallisto input files
############################################################

kallisto_files <- samples$kallisto_file

names(kallisto_files) <- samples$sample_id


############################################################
# 6. Read transcript IDs
############################################################

first_quantification <- read_tsv(
  unname(kallisto_files[1]),
  show_col_types = FALSE,
  progress = FALSE
)

transcript_ids <- first_quantification$target_id

stopifnot(
  length(transcript_ids) == 89013
)

stopifnot(
  !anyDuplicated(transcript_ids)
)


############################################################
# 7. Create transcript-to-gene mapping
#
# Example:
# DP100003.1 and DP100003.2 -> DP100003
############################################################

tx2gene <- data.frame(
  transcript_id = transcript_ids,
  gene_id = sub(
    "\\.[0-9]+$",
    "",
    transcript_ids
  )
)

number_transcripts <- nrow(tx2gene)

number_genes <- length(
  unique(tx2gene$gene_id)
)

cat(
  "Number of transcripts:",
  number_transcripts,
  "\n"
)

cat(
  "Number of genes:",
  number_genes,
  "\n"
)

stopifnot(
  number_transcripts == 89013
)

stopifnot(
  number_genes == 35916
)

write_tsv(
  tx2gene,
  file.path(
    metadata_dir,
    "tx2gene.tsv"
  )
)


############################################################
# 8. Import Kallisto estimates with tximport
############################################################

txi <- tximport(
  files = kallisto_files,
  type = "kallisto",
  tx2gene = tx2gene
)


############################################################
# 9. Validate imported data
############################################################

stopifnot(
  nrow(txi$counts) == 35916
)

stopifnot(
  ncol(txi$counts) == 18
)

stopifnot(
  identical(
    colnames(txi$counts),
    samples$sample_id
  )
)

cat(
  "Gene-level count matrix:",
  nrow(txi$counts),
  "genes x",
  ncol(txi$counts),
  "samples\n"
)


############################################################
# 10. Check Horz et al. transcription factors
############################################################

horz_tf_ids <- c(
  "DP112203",
  "DP103545",
  "DP109418",
  "DP102989",
  "DP105472",
  "DP111282"
)

tf_check <- data.frame(
  gene_id = horz_tf_ids,
  present = horz_tf_ids %in% rownames(txi$counts)
)

print(tf_check)

write_tsv(
  tf_check,
  file.path(
    results_dir,
    "Horz_TF_reference_check.tsv"
  )
)

stopifnot(
  all(tf_check$present)
)


############################################################
# 11. Save imported data
############################################################

saveRDS(
  txi,
  file.path(
    results_dir,
    "tximport_gene_level.rds"
  )
)

saveRDS(
  samples,
  file.path(
    metadata_dir,
    "samples.rds"
  )
)


############################################################
# 12. Completion message
############################################################

cat(
  "\ntximport completed successfully.\n"
)

cat(
  "No original Kallisto files were modified.\n"
)
