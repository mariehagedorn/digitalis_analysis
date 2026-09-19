#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 1: Image Import and ID Recognition
#
# Purpose:
# - Locate input images
# - Extract sample IDs from filenames
# - Match ruler and ruler-free images
# - Check completeness of image pairs
# - Confirm availability of negative-control images
#
# Raw image data are not included in this repository.
############################################################


############################################################
# 1. Define input directories
#
# Expected local directory structure:
#
# data/
# ├── 02_ruler_decolorized/
# ├── 04_decolorized/
# └── 05_negative_control/
############################################################

data_dir <- file.path(
  "DAB-Image_Quantification",
  "data"
)

ruler_dir <- file.path(
  data_dir,
  "02_ruler_decolorized"
)

decolorized_dir <- file.path(
  data_dir,
  "04_decolorized"
)

negative_dir <- file.path(
  data_dir,
  "05_negative_control"
)


############################################################
# 2. Locate image files
############################################################

ruler_files <- list.files(
  ruler_dir,
  pattern = "\\.JPG$",
  full.names = TRUE
)

decolorized_files <- list.files(
  decolorized_dir,
  pattern = "\\.JPG$",
  full.names = TRUE
)

negative_files <- list.files(
  negative_dir,
  pattern = "\\.JPG$",
  full.names = TRUE
)


############################################################
# 3. Extract sample IDs
############################################################

extract_MH_ID <- function(files) {

  sub(
    "^(MH[0-9]+).*",
    "\\1",
    basename(files)
  )
}


ruler_IDs <- extract_MH_ID(
  ruler_files
)

decolorized_IDs <- extract_MH_ID(
  decolorized_files
)


############################################################
# 4. Match ruler and ruler-free images
############################################################

ruler_table <- data.frame(
  MH_ID = ruler_IDs,
  ruler_file = ruler_files,
  stringsAsFactors = FALSE
)

decolorized_table <- data.frame(
  MH_ID = decolorized_IDs,
  decolorized_file = decolorized_files,
  stringsAsFactors = FALSE
)


sample_overview <- merge(
  ruler_table,
  decolorized_table,
  by = "MH_ID",
  all = TRUE
)


sample_overview$ruler_exists <-
  !is.na(sample_overview$ruler_file)

sample_overview$decolorized_exists <-
  !is.na(sample_overview$decolorized_file)


sample_overview <- sample_overview[
  order(sample_overview$MH_ID),
]


############################################################
# 5. Check image-pair completeness
############################################################

missing_ruler <- sample_overview$MH_ID[
  !sample_overview$ruler_exists
]

missing_decolorized <- sample_overview$MH_ID[
  !sample_overview$decolorized_exists
]


cat(
  "\nSample overview\n"
)

print(
  sample_overview
)


cat(
  "\nMissing ruler images:\n"
)

if (length(missing_ruler) == 0) {

  cat(
    "None\n"
  )

} else {

  print(
    missing_ruler
  )
}


cat(
  "\nMissing ruler-free images:\n"
)

if (length(missing_decolorized) == 0) {

  cat(
    "None\n"
  )

} else {

  print(
    missing_decolorized
  )
}


############################################################
# 6. Check negative-control images
############################################################

negative_check <- data.frame(

  negative_control_image = file.exists(
    file.path(
      negative_dir,
      "negative_control_histo_front.JPG"
    )
  ),

  negative_control_ruler = file.exists(
    file.path(
      negative_dir,
      "negative_control_histo_ruler_front.JPG"
    )
  )

)


cat(
  "\nNegative-control check:\n"
)

print(
  negative_check
)


############################################################
# 7. Completion message
############################################################

cat(
  "\nStep 1 completed successfully.\n"
)
