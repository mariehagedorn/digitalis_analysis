#!/usr/bin/env Rscript

# ============================================================
# 12_persistent_genotype_annotation_positions.R
#
# Add functional annotation and genomic positions to the
# stringent persistent aa-associated gene set.
#
# Expected input:
#   <project>/deseq2/results/genotype_background_analysis/
#       strict_aa_signature_57_genes_FINAL_with_shrinkage.rds
#   <project>/reference/DR1_v1.annoV2.anno.txt.gz
#   <project>/reference/DR1_v1.annoV2.anno.gff.gz
#
# By default, <project> is /vol/data/digitalis_rnaseq.
# ============================================================


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

reference_dir <- file.path(
  project_dir,
  "reference"
)

strict_file <- file.path(
  outdir,
  "strict_aa_signature_57_genes_FINAL_with_shrinkage.rds"
)

annotation_file <- file.path(
  reference_dir,
  "DR1_v1.annoV2.anno.txt.gz"
)

gff_file <- file.path(
  reference_dir,
  "DR1_v1.annoV2.anno.gff.gz"
)

for (f in c(
  strict_file,
  annotation_file,
  gff_file
)) {
  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }
}


# ------------------------------------------------------------
# 2. Load stringent gene set
# ------------------------------------------------------------

strict_aa <- readRDS(
  strict_file
)


# ------------------------------------------------------------
# 3. Read functional annotation
# ------------------------------------------------------------

anno_raw <- tryCatch(
  read.delim(
    gzfile(annotation_file),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  error = function(e) NULL
)

if (
  is.null(anno_raw) ||
    !all(
      c(
        "ID",
        "Anno"
      ) %in% names(anno_raw)
    )
) {

  anno_raw <- read.delim(
    gzfile(annotation_file),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )

  if (ncol(anno_raw) < 2) {
    stop(
      "Functional annotation file has fewer than two columns."
    )
  }

  names(anno_raw)[1:2] <- c(
    "ID",
    "Anno"
  )
}

anno_raw <- anno_raw[
  !is.na(anno_raw$ID) &
    anno_raw$ID != "",
]

anno_raw$gene_id <- sub(
  "\\.[0-9]+$",
  "",
  anno_raw$ID
)


# ------------------------------------------------------------
# 4. Parse Arabidopsis-based annotation
#
# The annotation file may contain:
#   AT-ID.description
# or
#   AT-ID.symbol.description
#
# The parser retains the original annotation and creates
# conservative display fields without assuming orthology.
# ------------------------------------------------------------

parse_one_annotation <- function(x) {

  if (
    is.na(x) ||
      trimws(x) == ""
  ) {
    return(
      c(
        at_homolog = NA_character_,
        function_text = NA_character_
      )
    )
  }

  x <- trimws(x)

  parts <- strsplit(
    x,
    "\\."
  )[[1]]

  parts <- trimws(parts)

  at_idx <- grep(
    "^AT[1-5CM]G[0-9]{5}$",
    parts,
    ignore.case = FALSE
  )

  if (length(at_idx) == 0) {

    at_match <- regmatches(
      x,
      regexpr(
        "AT[1-5CM]G[0-9]{5}",
        x
      )
    )

    if (length(at_match) == 0 || at_match == "") {
      return(
        c(
          at_homolog = NA_character_,
          function_text = x
        )
      )
    }

    return(
      c(
        at_homolog = at_match,
        function_text = x
      )
    )
  }

  i <- at_idx[1]
  at_id <- parts[i]

  following <- parts[
    seq.int(
      i + 1,
      length(parts)
    )
  ]

  following <- following[
    following != ""
  ]

  if (length(following) == 0) {
    return(
      c(
        at_homolog = at_id,
        function_text = NA_character_
      )
    )
  }

  first_after <- following[1]

  looks_like_symbol <- grepl(
    "^[A-Za-z0-9_-]{1,20}$",
    first_after
  )

  if (
    looks_like_symbol &&
      length(following) >= 2
  ) {

    homolog_label <- paste0(
      first_after,
      " (",
      at_id,
      ")"
    )

    function_text <- paste(
      following[-1],
      collapse = "."
    )

  } else {

    homolog_label <- at_id

    function_text <- paste(
      following,
      collapse = "."
    )
  }

  c(
    at_homolog = homolog_label,
    function_text = function_text
  )
}


parsed <- t(
  vapply(
    anno_raw$Anno,
    parse_one_annotation,
    FUN.VALUE = c(
      at_homolog = "",
      function_text = ""
    )
  )
)

anno_raw$putative_At_homolog <- parsed[
  ,
  "at_homolog"
]

anno_raw$functional_annotation <- parsed[
  ,
  "function_text"
]


collapse_unique <- function(x) {

  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]

  x <- unique(
    trimws(x)
  )

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = " / "
  )
}


anno_gene_ids <- unique(
  anno_raw$gene_id
)

anno_gene <- data.frame(
  gene_id = anno_gene_ids,
  putative_At_homolog = vapply(
    anno_gene_ids,
    function(id) {
      collapse_unique(
        anno_raw$putative_At_homolog[
          anno_raw$gene_id == id
        ]
      )
    },
    FUN.VALUE = character(1)
  ),
  functional_annotation = vapply(
    anno_gene_ids,
    function(id) {
      collapse_unique(
        anno_raw$functional_annotation[
          anno_raw$gene_id == id
        ]
      )
    },
    FUN.VALUE = character(1)
  ),
  annotation_raw = vapply(
    anno_gene_ids,
    function(id) {
      collapse_unique(
        anno_raw$Anno[
          anno_raw$gene_id == id
        ]
      )
    },
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Read GFF
# ------------------------------------------------------------

gff_raw <- read.delim(
  gzfile(gff_file),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE,
  fill = TRUE
)

if (ncol(gff_raw) < 9) {
  stop("GFF file does not contain 9 columns.")
}

gff_raw <- gff_raw[
  ,
  1:9
]

names(gff_raw) <- c(
  "seqid",
  "source",
  "type",
  "start",
  "end",
  "score",
  "strand",
  "phase",
  "attributes"
)


extract_attribute <- function(
  x,
  key
) {

  pattern <- paste0(
    "(?:^|;)",
    key,
    "=([^;]+)"
  )

  hit <- regexec(
    pattern,
    x,
    perl = TRUE
  )

  matches <- regmatches(
    x,
    hit
  )

  vapply(
    matches,
    function(z) {
      if (length(z) >= 2) {
        z[2]
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1)
  )
}


gff_genes <- gff_raw[
  gff_raw$type == "gene",
]

gff_genes$gene_id <- extract_attribute(
  gff_genes$attributes,
  "ID"
)

gff_genes$gene_id <- sub(
  "^gene:",
  "",
  gff_genes$gene_id
)

gff_genes$gene_id <- sub(
  "\\.[0-9]+$",
  "",
  gff_genes$gene_id
)

gff_genes <- gff_genes[
  !is.na(gff_genes$gene_id),
]

gff_genes <- gff_genes[
  !duplicated(gff_genes$gene_id),
]


# ------------------------------------------------------------
# 6. Determine contig lengths
#
# Preferred source:
#   ##sequence-region directives in the GFF header
#
# Fallback:
#   maximum annotated coordinate on each contig
# ------------------------------------------------------------

header_lines <- readLines(
  gzfile(gff_file),
  n = 5000,
  warn = FALSE
)

seq_region_lines <- grep(
  "^##sequence-region[[:space:]]+",
  header_lines,
  value = TRUE
)

if (length(seq_region_lines) > 0) {

  seq_parts <- strsplit(
    seq_region_lines,
    "[[:space:]]+"
  )

  contig_lengths <- data.frame(
    seqid = vapply(
      seq_parts,
      function(x) x[2],
      character(1)
    ),
    start = as.numeric(
      vapply(
        seq_parts,
        function(x) x[3],
        character(1)
      )
    ),
    end = as.numeric(
      vapply(
        seq_parts,
        function(x) x[4],
        character(1)
      )
    ),
    length_source = "GFF_sequence_region",
    stringsAsFactors = FALSE
  )

} else {

  warning(
    paste0(
      "No ##sequence-region directives were found in the first ",
      "5000 GFF header lines. Contig display lengths will use ",
      "the maximum annotated coordinate as a fallback."
    )
  )

  contig_lengths <- aggregate(
    end ~ seqid,
    data = gff_raw,
    FUN = max
  )

  contig_lengths$start <- 1

  contig_lengths$length_source <-
    "maximum_annotated_coordinate"

  contig_lengths <- contig_lengths[
    ,
    c(
      "seqid",
      "start",
      "end",
      "length_source"
    )
  ]
}


# ------------------------------------------------------------
# 7. Add annotations and genomic positions while preserving
#    the original stringent-table order
# ------------------------------------------------------------

strict_aa_annotated <- strict_aa

anno_match <- match(
  strict_aa_annotated$gene_id,
  anno_gene$gene_id
)

strict_aa_annotated$putative_At_homolog <-
  anno_gene$putative_At_homolog[
    anno_match
  ]

strict_aa_annotated$functional_annotation <-
  anno_gene$functional_annotation[
    anno_match
  ]

strict_aa_annotated$annotation_raw <-
  anno_gene$annotation_raw[
    anno_match
  ]

gff_match <- match(
  strict_aa_annotated$gene_id,
  gff_genes$gene_id
)

strict_aa_annotated$seqid <-
  gff_genes$seqid[
    gff_match
  ]

strict_aa_annotated$start <-
  gff_genes$start[
    gff_match
  ]

strict_aa_annotated$end <-
  gff_genes$end[
    gff_match
  ]

strict_aa_annotated$strand <-
  gff_genes$strand[
    gff_match
  ]


# ------------------------------------------------------------
# 8. ANS genomic position
# ------------------------------------------------------------

ans_id <- "DP129852"

ans_position <- gff_genes[
  gff_genes$gene_id == ans_id,
  c(
    "gene_id",
    "seqid",
    "start",
    "end",
    "strand"
  )
]

if (nrow(ans_position) == 0) {
  warning(
    "ANS (DP129852) was not found among GFF gene features."
  )
}


# ------------------------------------------------------------
# 9. Export
# ------------------------------------------------------------

saveRDS(
  strict_aa_annotated,
  file.path(
    outdir,
    "strict_aa_signature_57_genes_annotated_positions.rds"
  )
)

write.csv(
  strict_aa_annotated,
  file.path(
    outdir,
    "strict_aa_signature_57_genes_annotated_positions.csv"
  ),
  row.names = FALSE
)

saveRDS(
  gff_genes,
  file.path(
    outdir,
    "gff_genes.rds"
  )
)

saveRDS(
  contig_lengths,
  file.path(
    outdir,
    "contig_lengths.rds"
  )
)

write.csv(
  contig_lengths,
  file.path(
    outdir,
    "contig_lengths.csv"
  ),
  row.names = FALSE
)

write.csv(
  ans_position,
  file.path(
    outdir,
    "ANS_genomic_position.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

cat(
  "\nAnnotation and genomic-position analysis completed.\n",
  "Stringent genes: ",
  nrow(strict_aa_annotated),
  "\nGenes with functional annotation: ",
  sum(
    !is.na(
      strict_aa_annotated$functional_annotation
    )
  ),
  "\nGenes with genomic position: ",
  sum(
    !is.na(
      strict_aa_annotated$seqid
    )
  ),
  "\nANS position rows: ",
  nrow(ans_position),
  "\nContig-length source(s): ",
  paste(
    unique(
      contig_lengths$length_source
    ),
    collapse = ", "
  ),
  "\nOutput directory: ",
  outdir,
  "\n",
  sep = ""
)
