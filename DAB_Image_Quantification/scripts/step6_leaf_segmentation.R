#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 6: Final Leaf Mask Generation
#
# Purpose:
# - Use illumination-corrected analysis images from Step 4
# - Segment the leaf directly in the analysis image
# - Generate a binary leaf mask matching the full-resolution
#   analysis image
# - Use Step-5 scale information only for QC comparison
# - Prepare the image/mask pairs required for DAB
#   quantification in Step 7
#
# The Step-5 ruler-image mask is not transferred to the
# analysis image.
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
# 2. Load outputs from Steps 4 and 5
############################################################

source(
  file.path(
    "DAB-Image_Quantification",
    "scripts",
    "step4_illumination_correction.R"
  )
)

source(
  file.path(
    "DAB-Image_Quantification",
    "scripts",
    "step5_scale_calibration.R"
  )
)


if (!exists("illumination_corrected_image_table")) {
  stop(
    "Step 4 did not create 'illumination_corrected_image_table'."
  )
}

if (!exists("leaf_area_table")) {
  stop(
    "Step 5 did not create 'leaf_area_table'."
  )
}


############################################################
# 3. Parameters
############################################################

processing_long_side_px <- 2000

safety_border_fraction <- 0.10

threshold_min <- 12
threshold_max <- 70

opening_size <- 3
closing_size <- 11
final_closing_size <- 7

minimum_component_fraction <- 0.0015
minimum_interior_fraction <- 0.0008

area_difference_warning_percent <- 15


############################################################
# 4. Prepare input table
############################################################

analysis_images <-
  illumination_corrected_image_table[
    ,
    c(
      "MH_ID",
      "image_role",
      "corrected_file"
    )
  ]


step6_input <- merge(
  analysis_images,
  leaf_area_table[
    ,
    c(
      "MH_ID",
      "pixels_per_cm",
      "cm2_per_pixel",
      "leaf_area_cm2"
    )
  ],
  by = "MH_ID",
  all.x = TRUE,
  sort = FALSE
)


step6_input <- step6_input[
  file.exists(step6_input$corrected_file) &
  is.finite(step6_input$cm2_per_pixel),
]


############################################################
# 5. Output directories
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step6_leaf_segmentation"
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
# 6. Read analysis image
############################################################

read_analysis_image <- function(file_path) {

  full_image <- image_read(
    file_path
  )

  full_image <- image_convert(
    full_image,
    colorspace = "sRGB"
  )


  full_info <- image_info(
    full_image
  )


  processing_image <- image_resize(
    full_image,
    sprintf(
      "%dx%d>",
      processing_long_side_px,
      processing_long_side_px
    ),
    filter = "Lanczos"
  )


  processing_info <- image_info(
    processing_image
  )


  temporary_file <- tempfile(
    fileext = ".png"
  )

  image_write(
    processing_image,
    temporary_file,
    format = "png"
  )


  rgb <- png::readPNG(
    temporary_file
  )

  unlink(
    temporary_file
  )


  if (
    length(dim(rgb)) == 3 &&
    dim(rgb)[3] > 3
  ) {
    rgb <- rgb[, , 1:3]
  }


  list(
    rgb = rgb,

    full_width =
      full_info$width[1],

    full_height =
      full_info$height[1],

    processing_width =
      processing_info$width[1],

    processing_height =
      processing_info$height[1]
  )
}


############################################################
# 7. Otsu threshold
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


  counts <- tabulate(
    values + 1,
    nbins = 256
  )

  probabilities <-
    counts /
    sum(counts)

  intensity <- 0:255

  cumulative_probability <-
    cumsum(probabilities)

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


  valid <- denominator > 0


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
# 8. Morphological mask cleaning
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
# 9. Segment leaf in analysis image
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
  # Outer safety region
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


  safety_border <- matrix(
    FALSE,
    nrow = height,
    ncol = width
  )


  safety_border[
    seq_len(border_height),
  ] <- TRUE

  safety_border[
    (height - border_height + 1):height,
  ] <- TRUE

  safety_border[
    ,
    seq_len(border_width)
  ] <- TRUE

  safety_border[
    ,
    (width - border_width + 1):width
  ] <- TRUE


  interior_mask <- (
    !safety_border
  )


  ##########################################################
  # Background reference
  ##########################################################

  brightness_cutoff <- quantile(
    brightness[safety_border],
    0.60,
    na.rm = TRUE
  )


  saturation_cutoff <- quantile(
    saturation[safety_border],
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
    safety_border &
    brightness >= brightness_cutoff &
    saturation <= saturation_cutoff
  )


  if (
    sum(background_mask) < 500
  ) {

    brightness_cutoff <- quantile(
      brightness[safety_border],
      0.75,
      na.rm = TRUE
    )


    background_mask <- (
      safety_border &
      brightness >= brightness_cutoff
    )
  }


  if (
    sum(background_mask) < 100
  ) {

    stop(
      "Too few suitable background pixels detected."
    )
  }


  background_rgb <- c(

    red = median(
      red[background_mask],
      na.rm = TRUE
    ),

    green = median(
      green[background_mask],
      na.rm = TRUE
    ),

    blue = median(
      blue[background_mask],
      na.rm = TRUE
    )
  )


  ##########################################################
  # RGB distance from background
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


  candidate_mask[
    !is.finite(foreground_score)
  ] <- FALSE


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
    labels,
    na.rm = TRUE
  )


  if (
    component_count < 1
  ) {

    stop(
      "No connected leaf component detected."
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


  ##########################################################
  # Penalize components dominated by the actual image edge
  ##########################################################

  actual_edge <- matrix(
    FALSE,
    nrow = height,
    ncol = width
  )

  actual_edge[1, ] <- TRUE
  actual_edge[height, ] <- TRUE
  actual_edge[, 1] <- TRUE
  actual_edge[, width] <- TRUE


  edge_labels <- labels[
    actual_edge &
    labels > 0
  ]


  edge_pixels <- tabulate(
    edge_labels,
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


  interior_fraction <- (
    interior_pixels /
    pmax(
      component_areas,
      1
    )
  )


  edge_coverage <- (
    edge_pixels /
    max(
      sum(actual_edge),
      1
    )
  )


  component_score <- (
    interior_pixels *
    sqrt(
      pmax(
        interior_fraction,
        0
      )
    ) *
    pmax(
      0.15,
      1 -
        pmin(
          0.85,
          edge_coverage * 4
        )
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


  selected_component <- eligible[
    which.max(
      component_score[
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

    leaf_pixels =
      sum(final_mask),

    threshold =
      applied_threshold,

    background_red =
      background_rgb["red"] * 255,

    background_green =
      background_rgb["green"] * 255,

    background_blue =
      background_rgb["blue"] * 255
  )
}


############################################################
# 10. Save full-resolution leaf mask
#
# Nearest-neighbour interpolation preserves the binary mask.
############################################################

save_full_resolution_mask <- function(
  processing_mask,
  full_width,
  full_height,
  output_path
) {

  temporary_file <- tempfile(
    fileext = ".png"
  )


  png::writePNG(
    processing_mask * 1,
    temporary_file
  )


  mask_image <- image_read(
    temporary_file
  )


  full_mask_image <- image_resize(
    mask_image,
    sprintf(
      "%dx%d!",
      full_width,
      full_height
    ),
    filter = "point"
  )


  image_write(
    full_mask_image,
    output_path,
    format = "png"
  )


  unlink(
    temporary_file
  )


  full_mask <- png::readPNG(
    output_path
  )


  if (
    length(dim(full_mask)) == 3
  ) {

    full_mask <- full_mask[, , 1]
  }


  full_mask <- (
    full_mask > 0.5
  )


  if (
    nrow(full_mask) != full_height ||
    ncol(full_mask) != full_width
  ) {

    stop(
      "Full-resolution leaf mask dimensions do not match the analysis image."
    )
  }


  png::writePNG(
    full_mask * 1,
    output_path
  )


  full_mask
}


############################################################
# 11. Process all analysis images
############################################################

results <- vector(
  "list",
  nrow(step6_input)
)


for (
  i in seq_len(
    nrow(step6_input)
  )
) {

  sample_id <-
    step6_input$MH_ID[i]

  current_file <-
    step6_input$corrected_file[i]


  cat(
    "\nSegmenting ",
    sample_id,
    "...\n",
    sep = ""
  )


  results[[i]] <- tryCatch(

    {

      ######################################################
      # Read analysis image
      ######################################################

      image <- read_analysis_image(
        current_file
      )


      ######################################################
      # Segment leaf
      ######################################################

      segmentation <- segment_leaf(
        image$rgb
      )


      ######################################################
      # QC comparison with Step-5 leaf area
      ######################################################

      estimated_area_cm2 <- (
        segmentation$leaf_pixels *
        step6_input$cm2_per_pixel[i]
      )


      step5_area_cm2 <-
        step6_input$leaf_area_cm2[i]


      area_difference_percent <- (
        abs(
          estimated_area_cm2 -
          step5_area_cm2
        ) /
        max(
          step5_area_cm2,
          .Machine$double.eps
        ) *
        100
      )


      qc_warning <- (
        area_difference_percent >
        area_difference_warning_percent
      )


      ######################################################
      # Create full-resolution mask
      ######################################################

      mask_file <- file.path(
        mask_dir,
        paste0(
          sample_id,
          "_analysis_leaf_mask.png"
        )
      )


      full_mask <- save_full_resolution_mask(
        segmentation$mask,
        image$full_width,
        image$full_height,
        mask_file
      )


      ######################################################
      # Return compact analysis information
      ######################################################

      data.frame(

        MH_ID =
          sample_id,

        image_role =
          step6_input$image_role[i],

        analysis_image =
          current_file,

        leaf_mask_file =
          mask_file,

        image_width_px =
          image$full_width,

        image_height_px =
          image$full_height,

        leaf_pixels_full_resolution =
          sum(full_mask),

        step5_leaf_area_cm2 =
          step5_area_cm2,

        step6_estimated_leaf_area_cm2 =
          estimated_area_cm2,

        area_difference_percent =
          area_difference_percent,

        qc_warning =
          qc_warning,

        success =
          TRUE,

        stringsAsFactors = FALSE
      )
    },


    error = function(e) {

      warning(
        "Step 6 failed for ",
        sample_id,
        ": ",
        conditionMessage(e)
      )


      data.frame(

        MH_ID =
          sample_id,

        image_role =
          step6_input$image_role[i],

        analysis_image =
          current_file,

        leaf_mask_file = NA,

        image_width_px = NA,

        image_height_px = NA,

        leaf_pixels_full_resolution = NA,

        step5_leaf_area_cm2 =
          step6_input$leaf_area_cm2[i],

        step6_estimated_leaf_area_cm2 = NA,

        area_difference_percent = NA,

        qc_warning = NA,

        success =
          FALSE,

        stringsAsFactors = FALSE
      )
    }
  )
}


############################################################
# 12. Combine results
############################################################

leaf_segmentation_table <- do.call(
  rbind,
  results
)

rownames(
  leaf_segmentation_table
) <- NULL


############################################################
# 13. Prepare Step-7 input
############################################################

step7_input_table <- leaf_segmentation_table[
  leaf_segmentation_table$success,
  c(
    "MH_ID",
    "image_role",
    "analysis_image",
    "leaf_mask_file",
    "image_width_px",
    "image_height_px",
    "leaf_pixels_full_resolution",
    "step5_leaf_area_cm2"
  )
]


############################################################
# 14. Save tables
############################################################

write.csv(
  leaf_segmentation_table,
  file.path(
    output_dir,
    "step6_leaf_segmentation.csv"
  ),
  row.names = FALSE
)


saveRDS(
  step7_input_table,
  file.path(
    output_dir,
    "step6_DAB_analysis_input.rds"
  )
)


############################################################
# 15. Completion message
############################################################

cat(
  "\nStep 6 completed successfully.\n"
)
