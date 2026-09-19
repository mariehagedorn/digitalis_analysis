#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 3: Color Channel Normalization
#
# Purpose:
# - Detect neutral white background pixels near image edges
# - Use the background as an internal white-balance reference
# - Normalize RGB channels to a target value of 245
# - Apply the correction to full-resolution images
# - Generate reference masks and before/after QC images
#
# This step corrects global color differences only.
# Spatial illumination gradients are corrected in Step 4.
#
# Raw images and generated outputs are not included in
# this repository.
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(magick)
  library(png)
})


############################################################
# 2. Load image information from Step 1
############################################################

source(
  file.path(
    "DAB-Image_Quantification",
    "scripts",
    "step1_image_import.R"
  )
)

if (!exists("sample_overview")) {
  stop(
    "Step 1 did not create 'sample_overview'."
  )
}


############################################################
# 3. Define normalization parameters
############################################################

target_white <- 245

preview_geometry <- "1400x1400>"

border_fraction <- 0.10

primary_brightness_threshold <- 210

primary_chroma_threshold <- 30

fallback_chroma_threshold <- 45

minimum_reference_pixels <- 1000

gain_limits <- c(
  0.75,
  1.35
)


############################################################
# 4. Prepare input files
############################################################

sample_input <- data.frame(
  MH_ID = sample_overview$MH_ID,
  image_role = "sample",
  input_file = sample_overview$decolorized_file,
  stringsAsFactors = FALSE
)

sample_input <- sample_input[
  sample_overview$decolorized_exists,
]


negative_control_file <- file.path(
  negative_dir,
  "negative_control_histo_front.JPG"
)

if (!file.exists(negative_control_file)) {
  stop(
    "Negative-control image was not found."
  )
}


negative_control_input <- data.frame(
  MH_ID = "negative_control",
  image_role = "negative_control",
  input_file = negative_control_file,
  stringsAsFactors = FALSE
)


normalization_input <- rbind(
  sample_input,
  negative_control_input
)


############################################################
# 5. Define output directories
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step3_color_normalization"
)

normalized_dir <- file.path(
  output_dir,
  "normalized_images"
)

mask_dir <- file.path(
  output_dir,
  "white_reference_masks"
)

qc_dir <- file.path(
  output_dir,
  "QC"
)


for (directory in c(
  output_dir,
  normalized_dir,
  mask_dir,
  qc_dir
)) {

  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
# 6. Helper function: extract RGB matrices
############################################################

extract_rgb <- function(image) {

  image_array <- image_data(
    image,
    channels = "rgb"
  )

  image_array <- aperm(
    image_array,
    c(3, 2, 1)
  )

  list(
    red = matrix(
      as.integer(image_array[, , 1]),
      nrow = dim(image_array)[1]
    ),

    green = matrix(
      as.integer(image_array[, , 2]),
      nrow = dim(image_array)[1]
    ),

    blue = matrix(
      as.integer(image_array[, , 3]),
      nrow = dim(image_array)[1]
    )
  )
}


############################################################
# 7. Detect white-background reference pixels
############################################################

detect_white_reference <- function(rgb) {

  red <- rgb$red
  green <- rgb$green
  blue <- rgb$blue

  height <- nrow(red)
  width <- ncol(red)


  brightness <- (
    red +
    green +
    blue
  ) / 3


  chroma <- (
    pmax(
      red,
      green,
      blue
    ) -
    pmin(
      red,
      green,
      blue
    )
  )


  border_height <- max(
    1,
    round(
      height * border_fraction
    )
  )

  border_width <- max(
    1,
    round(
      width * border_fraction
    )
  )


  border_mask <- matrix(
    FALSE,
    nrow = height,
    ncol = width
  )

  border_mask[
    seq_len(border_height),
  ] <- TRUE

  border_mask[
    (height - border_height + 1):height,
  ] <- TRUE

  border_mask[
    ,
    seq_len(border_width)
  ] <- TRUE

  border_mask[
    ,
    (width - border_width + 1):width
  ] <- TRUE


  ##########################################################
  # Primary selection:
  # bright, low-chroma pixels near image edges
  ##########################################################

  reference_mask <- (
    border_mask &
    brightness >= primary_brightness_threshold &
    chroma <= primary_chroma_threshold
  )


  ##########################################################
  # Adaptive fallback if too few reference pixels are found
  ##########################################################

  if (
    sum(reference_mask) <
    minimum_reference_pixels
  ) {

    adaptive_threshold <- quantile(
      brightness[border_mask],
      0.85,
      na.rm = TRUE
    )

    reference_mask <- (
      border_mask &
      brightness >= adaptive_threshold &
      chroma <= fallback_chroma_threshold
    )
  }


  ##########################################################
  # Final fallback: brightest border pixels
  ##########################################################

  if (
    sum(reference_mask) <
    minimum_reference_pixels
  ) {

    brightest_threshold <- quantile(
      brightness[border_mask],
      0.95,
      na.rm = TRUE
    )

    reference_mask <- (
      border_mask &
      brightness >= brightest_threshold
    )
  }


  if (
    sum(reference_mask) < 50
  ) {

    stop(
      "Too few white-reference pixels were detected."
    )
  }


  reference_mask
}


############################################################
# 8. Apply RGB white balance
############################################################

apply_white_balance <- function(
  image,
  gains,
  output_path
) {

  rgb_raw <- image_data(
    image,
    channels = "rgb"
  )

  rgb_array <- aperm(
    rgb_raw,
    c(3, 2, 1)
  )

  corrected <- array(
    as.numeric(rgb_array) / 255,
    dim = dim(rgb_array)
  )


  for (channel in seq_len(3)) {

    corrected[, , channel] <-
      corrected[, , channel] *
      gains[channel]
  }


  corrected[
    corrected < 0
  ] <- 0

  corrected[
    corrected > 1
  ] <- 1


  png::writePNG(
    corrected,
    output_path
  )
}


############################################################
# 9. Create QC overlay
############################################################

create_reference_overlay <- function(
  image,
  reference_mask,
  output_path
) {

  overlay <- array(
    0,
    dim = c(
      nrow(reference_mask),
      ncol(reference_mask),
      4
    )
  )

  overlay[, , 1] <- 1

  overlay[, , 4] <- ifelse(
    reference_mask,
    0.45,
    0
  )


  temporary_file <- tempfile(
    fileext = ".png"
  )

  png::writePNG(
    overlay,
    temporary_file
  )


  overlay_image <- image_read(
    temporary_file
  )

  unlink(
    temporary_file
  )


  result <- image_composite(
    image,
    overlay_image,
    operator = "over"
  )


  image_write(
    result,
    output_path,
    format = "png"
  )
}


############################################################
# 10. Normalize all images
############################################################

normalization_results <- vector(
  "list",
  nrow(normalization_input)
)


for (
  i in seq_len(
    nrow(normalization_input)
  )
) {

  sample_id <-
    normalization_input$MH_ID[i]

  current_file <-
    normalization_input$input_file[i]


  cat(
    "Normalizing ",
    sample_id,
    "...\n",
    sep = ""
  )


  normalization_results[[i]] <- tryCatch(

    {

      ######################################################
      # Load and orient image
      ######################################################

      image <- image_read(
        current_file
      )

      image <- image_orient(
        image
      )

      image <- image_convert(
        image,
        colorspace = "sRGB"
      )


      ######################################################
      # Detect white reference on reduced preview
      ######################################################

      preview <- image_scale(
        image,
        preview_geometry
      )

      preview_rgb <- extract_rgb(
        preview
      )

      reference_mask <-
        detect_white_reference(
          preview_rgb
        )


      ######################################################
      # Median RGB values of white reference
      ######################################################

      reference_rgb <- c(

        red = median(
          preview_rgb$red[
            reference_mask
          ],
          na.rm = TRUE
        ),

        green = median(
          preview_rgb$green[
            reference_mask
          ],
          na.rm = TRUE
        ),

        blue = median(
          preview_rgb$blue[
            reference_mask
          ],
          na.rm = TRUE
        )
      )


      ######################################################
      # Calculate channel-specific correction factors
      ######################################################

      gains <- (
        target_white /
        reference_rgb
      )


      gains <- pmax(
        gain_limits[1],
        pmin(
          gain_limits[2],
          gains
        )
      )


      ######################################################
      # Normalize full-resolution image
      ######################################################

      normalized_file <- file.path(
        normalized_dir,
        paste0(
          sample_id,
          "_color_normalized.png"
        )
      )


      apply_white_balance(
        image = image,
        gains = gains,
        output_path = normalized_file
      )


      ######################################################
      # Save white-reference mask
      ######################################################

      mask_file <- file.path(
        mask_dir,
        paste0(
          sample_id,
          "_white_reference_mask.png"
        )
      )

      png::writePNG(
        reference_mask * 1,
        mask_file
      )


      ######################################################
      # Save visual QC overlay
      ######################################################

      overlay_file <- file.path(
        qc_dir,
        paste0(
          sample_id,
          "_white_reference_overlay.png"
        )
      )


      create_reference_overlay(
        image = preview,
        reference_mask = reference_mask,
        output_path = overlay_file
      )


      ######################################################
      # Return processing information
      ######################################################

      data.frame(
        MH_ID = sample_id,

        image_role =
          normalization_input$image_role[i],

        reference_pixels =
          sum(reference_mask),

        reference_red =
          reference_rgb["red"],

        reference_green =
          reference_rgb["green"],

        reference_blue =
          reference_rgb["blue"],

        gain_red =
          gains["red"],

        gain_green =
          gains["green"],

        gain_blue =
          gains["blue"],

        normalization_success =
          TRUE,

        normalized_file =
          normalized_file,

        stringsAsFactors = FALSE
      )
    },


    error = function(e) {

      warning(
        "Color normalization failed for ",
        sample_id,
        ": ",
        conditionMessage(e)
      )


      data.frame(
        MH_ID = sample_id,
        image_role =
          normalization_input$image_role[i],
        reference_pixels = NA,
        reference_red = NA,
        reference_green = NA,
        reference_blue = NA,
        gain_red = NA,
        gain_green = NA,
        gain_blue = NA,
        normalization_success = FALSE,
        normalized_file = NA,
        stringsAsFactors = FALSE
      )
    }
  )
}


############################################################
# 11. Combine results
############################################################

color_normalization_table <- do.call(
  rbind,
  normalization_results
)

rownames(
  color_normalization_table
) <- NULL


############################################################
# 12. Save processing summary
############################################################

write.csv(
  color_normalization_table,
  file.path(
    output_dir,
    "step3_color_normalization_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 13. Completion message
############################################################

cat(
  "\nStep 3 completed successfully.\n"
)
