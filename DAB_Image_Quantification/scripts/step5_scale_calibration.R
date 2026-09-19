#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 5: Scale Calibration and Leaf Area Determination
#
# Purpose:
# - Use ruler-containing images identified in Step 1
# - Manually select two ruler marks exactly 1 cm apart
# - Calculate image-specific pixel-to-cm conversion factors
# - Automatically segment the leaf
# - Determine leaf area in cm² from the final leaf mask
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
  library(EBImage)
})


############################################################
# 2. Run Step 1
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
# 3. Parameters
############################################################

processing_long_side_px <- 2000

display_long_side_px <- 2400

reference_length_cm <- 1

safety_border_fraction <- 0.10

threshold_min <- 12
threshold_max <- 70

opening_size <- 3
closing_size <- 11
final_closing_size <- 7

minimum_component_fraction <- 0.0015
minimum_interior_fraction <- 0.0008


############################################################
# 4. Prepare ruler-image input
############################################################

step5_input <- data.frame(
  MH_ID = sample_overview$MH_ID,
  image_role = "sample",
  input_file = sample_overview$ruler_file,
  stringsAsFactors = FALSE
)

step5_input <- step5_input[
  sample_overview$ruler_exists,
]


negative_control_file <- file.path(
  negative_dir,
  "negative_control_histo_ruler_front.JPG"
)

if (file.exists(negative_control_file)) {

  step5_input <- rbind(
    step5_input,

    data.frame(
      MH_ID = "negative_control",
      image_role = "negative_control",
      input_file = negative_control_file,
      stringsAsFactors = FALSE
    )
  )
}


############################################################
# 5. Output directories
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step5_scale_calibration"
)

mask_dir <- file.path(
  output_dir,
  "leaf_masks"
)

qc_dir <- file.path(
  output_dir,
  "QC"
)


dir.create(
  mask_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 6. Read image and create processing/display versions
############################################################

prepare_image <- function(file_path) {

  image <- image_read(
    file_path
  )

  image <- image_orient(
    image
  )

  image <- image_convert(
    image,
    colorspace = "sRGB"
  )


  processing_image <- image_resize(
    image,
    sprintf(
      "%dx%d>",
      processing_long_side_px,
      processing_long_side_px
    ),
    filter = "Lanczos"
  )


  display_image <- image_resize(
    image,
    sprintf(
      "%dx%d>",
      display_long_side_px,
      display_long_side_px
    ),
    filter = "Lanczos"
  )


  processing_info <- image_info(
    processing_image
  )

  display_info <- image_info(
    display_image
  )


  processing_file <- tempfile(
    fileext = ".png"
  )

  image_write(
    processing_image,
    processing_file,
    format = "png"
  )

  processing_rgb <- png::readPNG(
    processing_file
  )

  unlink(
    processing_file
  )


  if (dim(processing_rgb)[3] > 3) {
    processing_rgb <- processing_rgb[, , 1:3]
  }


  list(
    processing_rgb = processing_rgb,

    display_image = display_image,

    processing_width =
      processing_info$width[1],

    processing_height =
      processing_info$height[1],

    display_width =
      display_info$width[1],

    display_height =
      display_info$height[1]
  )
}


############################################################
# 7. Manual ruler calibration
############################################################

manual_calibration <- function(
  display_image,
  display_width,
  display_height,
  processing_width,
  processing_height,
  sample_id
) {

  if (!interactive()) {

    stop(
      "Manual ruler calibration requires an interactive R session."
    )
  }


  plot(
    display_image,
    main = paste(
      sample_id,
      "- select two ruler marks exactly 1 cm apart"
    )
  )


  points_clicked <- locator(
    n = 2,
    type = "p",
    pch = 16,
    col = "red"
  )


  if (
    is.null(points_clicked) ||
    length(points_clicked$x) != 2
  ) {

    stop(
      "Ruler calibration was not completed."
    )
  }


  ##########################################################
  # Convert clicked coordinates from display image
  # to processing-image coordinates
  ##########################################################

  scale_x <- (
    processing_width /
    display_width
  )

  scale_y <- (
    processing_height /
    display_height
  )


  x_processing <-
    points_clicked$x *
    scale_x

  y_processing <-
    points_clicked$y *
    scale_y


  distance_px <- sqrt(
    (
      x_processing[2] -
      x_processing[1]
    )^2 +
    (
      y_processing[2] -
      y_processing[1]
    )^2
  )


  pixels_per_cm <- (
    distance_px /
    reference_length_cm
  )


  cm_per_pixel <- (
    reference_length_cm /
    distance_px
  )


  cm2_per_pixel <- (
    cm_per_pixel^2
  )


  list(
    reference_length_px = distance_px,
    pixels_per_cm = pixels_per_cm,
    cm_per_pixel = cm_per_pixel,
    cm2_per_pixel = cm2_per_pixel
  )
}


############################################################
# 8. Otsu threshold
############################################################

otsu_threshold <- function(values) {

  values <- values[
    is.finite(values)
  ]

  values <- round(
    pmax(
      0,
      pmin(
        255,
        values
      )
    )
  )


  histogram_counts <- tabulate(
    values + 1,
    nbins = 256
  )


  probabilities <- (
    histogram_counts /
    sum(histogram_counts)
  )


  intensity <- 0:255

  cumulative_probability <-
    cumsum(
      probabilities
    )

  cumulative_mean <-
    cumsum(
      probabilities *
      intensity
    )


  total_mean <-
    cumulative_mean[256]


  denominator <- (
    cumulative_probability *
    (
      1 -
      cumulative_probability
    )
  )


  between_class_variance <-
    rep(
      -Inf,
      256
    )


  valid <- (
    denominator > 0
  )


  between_class_variance[valid] <- (
    total_mean *
    cumulative_probability[valid] -
    cumulative_mean[valid]
  )^2 /
    denominator[valid]


  which.max(
    between_class_variance
  ) - 1
}


############################################################
# 9. Morphological mask cleaning
############################################################

clean_mask <- function(
  mask,
  opening,
  closing
) {

  mask_image <- EBImage::Image(
    mask * 1,
    colormode = EBImage::Grayscale
  )


  if (opening > 1) {

    mask_image <- EBImage::opening(
      mask_image,
      EBImage::makeBrush(
        opening,
        shape = "disc"
      )
    )
  }


  if (closing > 1) {

    mask_image <- EBImage::closing(
      mask_image,
      EBImage::makeBrush(
        closing,
        shape = "disc"
      )
    )
  }


  mask_image <- EBImage::fillHull(
    mask_image
  )


  as.array(
    mask_image
  ) > 0.5
}


############################################################
# 10. Automatic leaf segmentation
############################################################

segment_leaf <- function(rgb) {

  height <- dim(rgb)[1]
  width <- dim(rgb)[2]

  total_pixels <- (
    height *
    width
  )


  red <- rgb[, , 1]
  green <- rgb[, , 2]
  blue <- rgb[, , 3]


  brightness <- (
    red +
    green +
    blue
  ) / 3


  maximum_channel <- pmax(
    red,
    green,
    blue
  )

  minimum_channel <- pmin(
    red,
    green,
    blue
  )

  chroma <- (
    maximum_channel -
    minimum_channel
  )

  saturation <- (
    chroma /
    pmax(
      maximum_channel,
      1 / 255
    )
  )


  ##########################################################
  # Outer 10% safety region for background estimation
  ##########################################################

  border_height <- max(
    1,
    round(
      height *
      safety_border_fraction
    )
  )

  border_width <- max(
    1,
    round(
      width *
      safety_border_fraction
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


  interior_mask <- (
    !border_mask
  )


  ##########################################################
  # Estimate white-background RGB reference
  ##########################################################

  brightness_cutoff <- quantile(
    brightness[border_mask],
    0.60,
    na.rm = TRUE
  )


  saturation_cutoff <- quantile(
    saturation[border_mask],
    0.75,
    na.rm = TRUE
  )


  saturation_cutoff <- min(
    0.15,
    max(
      0.03,
      saturation_cutoff
    )
  )


  background_mask <- (
    border_mask &
    brightness >= brightness_cutoff &
    saturation <= saturation_cutoff
  )


  if (
    sum(background_mask) < 500
  ) {

    brightness_cutoff <- quantile(
      brightness[border_mask],
      0.75,
      na.rm = TRUE
    )

    background_mask <- (
      border_mask &
      brightness >= brightness_cutoff
    )
  }


  if (
    sum(background_mask) < 100
  ) {

    stop(
      "Too few background pixels detected."
    )
  }


  background_rgb <- c(
    red = median(
      red[background_mask]
    ),

    green = median(
      green[background_mask]
    ),

    blue = median(
      blue[background_mask]
    )
  )


  ##########################################################
  # RGB distance from white background
  ##########################################################

  foreground_score <- sqrt(
    (
      red -
      background_rgb["red"]
    )^2 +
    (
      green -
      background_rgb["green"]
    )^2 +
    (
      blue -
      background_rgb["blue"]
    )^2
  ) /
    sqrt(3) *
    255


  raw_threshold <- otsu_threshold(
    foreground_score
  )


  applied_threshold <- max(
    threshold_min,
    min(
      threshold_max,
      raw_threshold
    )
  )


  candidate_mask <- (
    foreground_score >=
    applied_threshold
  )


  candidate_mask <- clean_mask(
    candidate_mask,
    opening_size,
    closing_size
  )


  ##########################################################
  # Connected components
  ##########################################################

  labels <- EBImage::bwlabel(
    EBImage::Image(
      candidate_mask * 1,
      colormode = EBImage::Grayscale
    )
  )

  labels <- round(
    as.array(labels)
  )


  component_count <- max(
    labels
  )


  if (
    component_count < 1
  ) {

    stop(
      "No leaf component detected."
    )
  }


  component_areas <- tabulate(
    labels[
      labels > 0
    ],
    nbins = component_count
  )


  interior_labels <- labels[
    interior_mask &
    labels > 0
  ]


  interior_pixels <- tabulate(
    interior_labels,
    nbins = component_count
  )


  minimum_component_pixels <- max(
    100,
    round(
      total_pixels *
      minimum_component_fraction
    )
  )


  minimum_interior_pixels <- max(
    50,
    round(
      total_pixels *
      minimum_interior_fraction
    )
  )


  eligible <- which(
    component_areas >=
      minimum_component_pixels &
    interior_pixels >=
      minimum_interior_pixels
  )


  if (
    length(eligible) == 0
  ) {

    stop(
      "No plausible leaf component detected."
    )
  }


  ##########################################################
  # Prefer component with greatest presence in image interior
  ##########################################################

  selected_component <- eligible[
    which.max(
      interior_pixels[
        eligible
      ]
    )
  ]


  final_mask <- (
    labels ==
    selected_component
  )


  final_mask <- clean_mask(
    final_mask,
    1,
    final_closing_size
  )


  list(
    mask = final_mask,
    leaf_pixels = sum(
      final_mask
    ),
    threshold = applied_threshold
  )
}


############################################################
# 11. Process all ruler images
############################################################

results <- vector(
  "list",
  nrow(step5_input)
)


for (
  i in seq_len(
    nrow(step5_input)
  )
) {

  sample_id <-
    step5_input$MH_ID[i]

  current_file <-
    step5_input$input_file[i]


  cat(
    "\nProcessing ",
    sample_id,
    "...\n",
    sep = ""
  )


  results[[i]] <- tryCatch(

    {

      ######################################################
      # Prepare image
      ######################################################

      image <- prepare_image(
        current_file
      )


      ######################################################
      # Manual 1-cm ruler calibration
      ######################################################

      calibration <- manual_calibration(

        display_image =
          image$display_image,

        display_width =
          image$display_width,

        display_height =
          image$display_height,

        processing_width =
          image$processing_width,

        processing_height =
          image$processing_height,

        sample_id =
          sample_id
      )


      ######################################################
      # Leaf segmentation
      ######################################################

      segmentation <- segment_leaf(
        image$processing_rgb
      )


      ######################################################
      # Convert leaf pixels to cm²
      ######################################################

      leaf_area_cm2 <- (
        segmentation$leaf_pixels *
        calibration$cm2_per_pixel
      )


      ######################################################
      # Save final leaf mask
      ######################################################

      mask_file <- file.path(
        mask_dir,
        paste0(
          sample_id,
          "_leaf_mask.png"
        )
      )


      png::writePNG(
        segmentation$mask * 1,
        mask_file
      )


      ######################################################
      # Store leaf-level measurements
      ######################################################

      data.frame(

        MH_ID =
          sample_id,

        image_role =
          step5_input$image_role[i],

        reference_length_px =
          calibration$reference_length_px,

        pixels_per_cm =
          calibration$pixels_per_cm,

        cm2_per_pixel =
          calibration$cm2_per_pixel,

        leaf_area_pixels =
          segmentation$leaf_pixels,

        leaf_area_cm2 =
          leaf_area_cm2,

        segmentation_threshold =
          segmentation$threshold,

        leaf_mask_file =
          mask_file,

        success =
          TRUE,

        stringsAsFactors = FALSE
      )
    },


    error = function(e) {

      warning(
        "Step 5 failed for ",
        sample_id,
        ": ",
        conditionMessage(e)
      )


      data.frame(

        MH_ID =
          sample_id,

        image_role =
          step5_input$image_role[i],

        reference_length_px = NA,
        pixels_per_cm = NA,
        cm2_per_pixel = NA,
        leaf_area_pixels = NA,
        leaf_area_cm2 = NA,
        segmentation_threshold = NA,
        leaf_mask_file = NA,
        success = FALSE,

        stringsAsFactors = FALSE
      )
    }
  )
}


############################################################
# 12. Combine results
############################################################

leaf_area_table <- do.call(
  rbind,
  results
)

rownames(
  leaf_area_table
) <- NULL


############################################################
# 13. Save leaf-area table
############################################################

write.csv(
  leaf_area_table,
  file.path(
    output_dir,
    "step5_leaf_area.csv"
  ),
  row.names = FALSE
)


saveRDS(
  leaf_area_table,
  file.path(
    output_dir,
    "step5_leaf_area.rds"
  )
)


############################################################
# 14. Completion message
############################################################

cat(
  "\nStep 5 completed successfully.\n"
)
