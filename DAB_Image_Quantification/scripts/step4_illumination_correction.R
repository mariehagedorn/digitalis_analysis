#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 4: Spatial Illumination Correction
#
# Purpose:
# - Use color-normalized images generated in Step 3
# - Detect neutral, border-connected background regions
# - Sample background brightness across a spatial grid
# - Estimate a smooth two-dimensional illumination surface
# - Correct position-dependent brightness variation
# - Apply the same local gain to all RGB channels
# - Generate background masks and diagnostic maps
#
# Step 3 corrects global color differences.
# Step 4 corrects spatial brightness gradients.
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
# 2. Run Step 3
############################################################

source(
  file.path(
    "DAB-Image_Quantification",
    "scripts",
    "step3_color_normalization.R"
  )
)

if (!exists("color_normalization_table")) {
  stop(
    "Step 3 did not create 'color_normalization_table'."
  )
}


############################################################
# 3. Define illumination-correction parameters
############################################################

preview_geometry <- "1200x1200>"

target_background <- target_white

border_fraction <- 0.22

primary_brightness_threshold <- 165
primary_chroma_threshold <- 25

fallback_brightness_threshold <- 120
fallback_chroma_threshold <- 40

adaptive_brightness_quantile <- 0.60

grid_long_side_cells <- 16
grid_short_side_cells <- 12

minimum_pixels_per_cell <- 80
minimum_grid_cells <- 12

minimum_background_component_pixels <- 200

outlier_mad_multiplier <- 3.5
minimum_residual_threshold <- 2

gain_lower_limit <- 0.80
gain_upper_limit <- 1.25


############################################################
# 4. Prepare Step-4 input
############################################################

step4_input <- color_normalization_table[
  color_normalization_table$normalization_success,
  c(
    "MH_ID",
    "image_role",
    "normalized_file"
  )
]

step4_input <- step4_input[
  file.exists(step4_input$normalized_file),
]

rownames(step4_input) <- NULL


if (nrow(step4_input) == 0) {
  stop(
    "No normalized Step-3 images are available."
  )
}


############################################################
# 5. Define output directories
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step4_illumination_correction"
)

corrected_dir <- file.path(
  output_dir,
  "corrected_images"
)

mask_dir <- file.path(
  output_dir,
  "background_masks"
)

map_dir <- file.path(
  output_dir,
  "illumination_maps"
)

qc_dir <- file.path(
  output_dir,
  "QC"
)


for (directory in c(
  output_dir,
  corrected_dir,
  mask_dir,
  map_dir,
  qc_dir
)) {

  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
# 6. Extract RGB matrices
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

  h <- dim(image_array)[1]
  w <- dim(image_array)[2]


  list(
    red = matrix(
      as.integer(image_array[, , 1]),
      nrow = h,
      ncol = w
    ),

    green = matrix(
      as.integer(image_array[, , 2]),
      nrow = h,
      ncol = w
    ),

    blue = matrix(
      as.integer(image_array[, , 3]),
      nrow = h,
      ncol = w
    )
  )
}


############################################################
# 7. Calculate brightness and chroma
############################################################

brightness_chroma <- function(rgb) {

  brightness <- (
    rgb$red +
    rgb$green +
    rgb$blue
  ) / 3


  chroma <- (
    pmax(
      rgb$red,
      rgb$green,
      rgb$blue
    ) -
    pmin(
      rgb$red,
      rgb$green,
      rgb$blue
    )
  )


  list(
    brightness = brightness,
    chroma = chroma
  )
}


############################################################
# 8. Create outer-border mask
############################################################

create_border_mask <- function(
  height,
  width,
  fraction
) {

  border_rows <- max(
    1,
    floor(height * fraction)
  )

  border_columns <- max(
    1,
    floor(width * fraction)
  )


  mask <- matrix(
    FALSE,
    nrow = height,
    ncol = width
  )


  mask[
    seq_len(border_rows),
  ] <- TRUE

  mask[
    (height - border_rows + 1):height,
  ] <- TRUE

  mask[
    ,
    seq_len(border_columns)
  ] <- TRUE

  mask[
    ,
    (width - border_columns + 1):width
  ] <- TRUE


  mask
}


############################################################
# 9. Retain large border-connected components
#
# Four-neighbour connectivity is used.
# Components that do not touch the outer image border or
# contain too few pixels are removed.
############################################################

keep_border_components <- function(
  candidate_mask,
  minimum_pixels
) {

  candidate_mask[is.na(candidate_mask)] <- FALSE

  height <- nrow(candidate_mask)
  width <- ncol(candidate_mask)

  candidate_vector <- as.vector(
    candidate_mask
  )

  retained <- rep(
    FALSE,
    length(candidate_vector)
  )

  visited <- rep(
    FALSE,
    length(candidate_vector)
  )


  ##########################################################
  # Identify selected pixels touching the outer image edge
  ##########################################################

  top <- which(
    candidate_mask[1, ]
  )

  bottom <- which(
    candidate_mask[height, ]
  )

  left <- which(
    candidate_mask[, 1]
  )

  right <- which(
    candidate_mask[, width]
  )


  border_seeds <- unique(
    c(
      1 + (top - 1) * height,

      height +
        (bottom - 1) * height,

      left,

      right +
        (width - 1) * height
    )
  )


  if (length(border_seeds) == 0) {

    return(
      matrix(
        FALSE,
        nrow = height,
        ncol = width
      )
    )
  }


  ##########################################################
  # Connected-component search
  ##########################################################

  for (seed in border_seeds) {

    if (
      !candidate_vector[seed] ||
      visited[seed]
    ) {
      next
    }


    queue <- seed

    visited[seed] <- TRUE

    component <- integer(0)


    while (length(queue) > 0) {

      current <- queue[1]

      queue <- queue[-1]

      component <- c(
        component,
        current
      )


      row <- (
        (current - 1) %% height
      ) + 1

      column <- (
        (current - 1) %/% height
      ) + 1


      neighbours <- integer(0)


      if (row > 1) {
        neighbours <- c(
          neighbours,
          current - 1
        )
      }

      if (row < height) {
        neighbours <- c(
          neighbours,
          current + 1
        )
      }

      if (column > 1) {
        neighbours <- c(
          neighbours,
          current - height
        )
      }

      if (column < width) {
        neighbours <- c(
          neighbours,
          current + height
        )
      }


      neighbours <- neighbours[
        candidate_vector[neighbours] &
        !visited[neighbours]
      ]


      if (length(neighbours) > 0) {

        visited[neighbours] <- TRUE

        queue <- c(
          queue,
          neighbours
        )
      }
    }


    if (
      length(component) >=
      minimum_pixels
    ) {

      retained[component] <- TRUE
    }
  }


  matrix(
    retained,
    nrow = height,
    ncol = width
  )
}


############################################################
# 10. Determine grid dimensions
############################################################

get_grid_dimensions <- function(
  height,
  width
) {

  if (width >= height) {

    list(
      rows = grid_short_side_cells,
      columns = grid_long_side_cells
    )

  } else {

    list(
      rows = grid_long_side_cells,
      columns = grid_short_side_cells
    )
  }
}


############################################################
# 11. Summarize background pixels within grid cells
############################################################

summarize_grid <- function(
  rgb,
  brightness,
  background_mask,
  grid_rows,
  grid_columns
) {

  height <- nrow(
    brightness
  )

  width <- ncol(
    brightness
  )


  row_breaks <- round(
    seq(
      1,
      height + 1,
      length.out = grid_rows + 1
    )
  )


  column_breaks <- round(
    seq(
      1,
      width + 1,
      length.out = grid_columns + 1
    )
  )


  results <- list()


  for (
    grid_row in seq_len(grid_rows)
  ) {

    row_start <-
      row_breaks[grid_row]

    row_end <-
      row_breaks[grid_row + 1] - 1


    for (
      grid_column in seq_len(grid_columns)
    ) {

      column_start <-
        column_breaks[grid_column]

      column_end <-
        column_breaks[grid_column + 1] - 1


      local_mask <- background_mask[
        row_start:row_end,
        column_start:column_end,
        drop = FALSE
      ]


      n_pixels <- sum(
        local_mask
      )


      if (
        n_pixels <
        minimum_pixels_per_cell
      ) {
        next
      }


      local_brightness <- brightness[
        row_start:row_end,
        column_start:column_end,
        drop = FALSE
      ]


      results[[length(results) + 1]] <-
        data.frame(

          x = (
            mean(
              c(
                column_start,
                column_end
              )
            ) - 1
          ) /
            max(
              1,
              width - 1
            ),

          y = (
            mean(
              c(
                row_start,
                row_end
              )
            ) - 1
          ) /
            max(
              1,
              height - 1
            ),

          background_pixels =
            n_pixels,

          median_luminance =
            median(
              local_brightness[
                local_mask
              ],
              na.rm = TRUE
            )
        )
    }
  }


  if (length(results) == 0) {
    return(
      data.frame()
    )
  }


  do.call(
    rbind,
    results
  )
}


############################################################
# 12. Detect background pixels
############################################################

detect_background <- function(rgb) {

  bc <- brightness_chroma(
    rgb
  )

  brightness <- bc$brightness
  chroma <- bc$chroma


  height <- nrow(
    brightness
  )

  width <- ncol(
    brightness
  )


  border <- create_border_mask(
    height,
    width,
    border_fraction
  )


  grid_dimensions <- get_grid_dimensions(
    height,
    width
  )


  ##########################################################
  # Internal function for testing one threshold combination
  ##########################################################

  test_mask <- function(
    brightness_threshold,
    chroma_threshold
  ) {

    candidate <- (
      border &
      brightness >= brightness_threshold &
      chroma <= chroma_threshold
    )


    cleaned <- keep_border_components(
      candidate,
      minimum_background_component_pixels
    )


    grid <- summarize_grid(
      rgb,
      brightness,
      cleaned,
      grid_dimensions$rows,
      grid_dimensions$columns
    )


    list(
      mask = cleaned,
      grid = grid
    )
  }


  ##########################################################
  # Primary thresholds
  ##########################################################

  primary <- test_mask(
    primary_brightness_threshold,
    primary_chroma_threshold
  )


  if (
    nrow(primary$grid) >=
    minimum_grid_cells
  ) {

    return(
      list(
        mask = primary$mask,
        grid = primary$grid,
        brightness = brightness
      )
    )
  }


  ##########################################################
  # Fallback thresholds
  ##########################################################

  fallback <- test_mask(
    fallback_brightness_threshold,
    fallback_chroma_threshold
  )


  if (
    nrow(fallback$grid) >=
    minimum_grid_cells
  ) {

    return(
      list(
        mask = fallback$mask,
        grid = fallback$grid,
        brightness = brightness
      )
    )
  }


  ##########################################################
  # Adaptive threshold
  ##########################################################

  neutral_border <- (
    border &
    chroma <= fallback_chroma_threshold
  )


  adaptive_threshold <- quantile(
    brightness[neutral_border],
    adaptive_brightness_quantile,
    na.rm = TRUE
  )


  adaptive_threshold <- max(
    140,
    adaptive_threshold
  )


  adaptive_candidate <- (
    neutral_border &
    brightness >= adaptive_threshold
  )


  adaptive_mask <- keep_border_components(
    adaptive_candidate,
    minimum_background_component_pixels
  )


  adaptive_grid <- summarize_grid(
    rgb,
    brightness,
    adaptive_mask,
    grid_dimensions$rows,
    grid_dimensions$columns
  )


  if (
    nrow(adaptive_grid) <
    minimum_grid_cells
  ) {

    stop(
      "Insufficient background grid cells for illumination fitting."
    )
  }


  list(
    mask = adaptive_mask,
    grid = adaptive_grid,
    brightness = brightness
  )
}


############################################################
# 13. Fit smooth two-dimensional illumination surface
#
# A quadratic polynomial is used to capture broad spatial
# brightness gradients without following small image features.
############################################################

fit_illumination_surface <- function(
  grid
) {

  model_data <- grid

  model_data$x2 <-
    model_data$x^2

  model_data$y2 <-
    model_data$y^2

  model_data$xy <-
    model_data$x *
    model_data$y


  initial_model <- lm(
    median_luminance ~
      x +
      y +
      x2 +
      y2 +
      xy,
    data = model_data,
    weights = background_pixels
  )


  residual_values <- residuals(
    initial_model
  )


  residual_mad <- mad(
    residual_values,
    center = median(
      residual_values
    ),
    constant = 1.4826
  )


  residual_limit <- max(
    minimum_residual_threshold,
    outlier_mad_multiplier *
      residual_mad
  )


  keep <- (
    abs(residual_values) <=
    residual_limit
  )


  if (
    sum(keep) >= minimum_grid_cells &&
    any(!keep)
  ) {

    final_model <- lm(
      median_luminance ~
        x +
        y +
        x2 +
        y2 +
        xy,
      data = model_data[
        keep,
      ],
      weights =
        background_pixels[
          keep
        ]
    )

  } else {

    final_model <-
      initial_model
  }


  final_model
}


############################################################
# 14. Predict illumination surface
############################################################

predict_surface <- function(
  model,
  height,
  width
) {

  coefficients <- coef(
    model
  )


  x <- seq(
    0,
    1,
    length.out = width
  )

  y <- seq(
    0,
    1,
    length.out = height
  )


  surface <- matrix(
    NA_real_,
    nrow = height,
    ncol = width
  )


  for (
    row_index in seq_len(height)
  ) {

    surface[row_index, ] <-

      coefficients["(Intercept)"] +

      coefficients["x"] * x +

      coefficients["y"] *
        y[row_index] +

      coefficients["x2"] *
        x^2 +

      coefficients["y2"] *
        y[row_index]^2 +

      coefficients["xy"] *
        x *
        y[row_index]
  }


  surface
}


############################################################
# 15. Convert illumination surface to local gains
############################################################

calculate_gain <- function(
  illumination_surface
) {

  gain <- (
    target_background /
    illumination_surface
  )


  gain[
    gain < gain_lower_limit
  ] <- gain_lower_limit

  gain[
    gain > gain_upper_limit
  ] <- gain_upper_limit


  gain
}


############################################################
# 16. Apply the same local gain to all RGB channels
############################################################

apply_illumination_correction <- function(
  image,
  gain_matrix,
  output_path
) {

  rgb_raw <- image_data(
    image,
    channels = "rgb"
  )


  width <- dim(
    rgb_raw
  )[2]

  height <- dim(
    rgb_raw
  )[3]


  corrected <- array(
    as.raw(255),
    dim = c(
      4,
      width,
      height
    )
  )


  corrected[1:3, , ] <-
    rgb_raw


  gain_transposed <- t(
    gain_matrix
  )


  for (channel in 1:3) {

    values <- matrix(
      as.integer(
        corrected[
          channel,
          ,
        ]
      ),
      nrow = width,
      ncol = height
    )


    values <- round(
      values *
      gain_transposed
    )


    values[
      values < 0
    ] <- 0

    values[
      values > 255
    ] <- 255


    corrected[
      channel,
      ,
    ] <- as.raw(
      values
    )
  }


  png::writePNG(
    corrected,
    output_path
  )
}


############################################################
# 17. Process all normalized images
############################################################

illumination_results <- vector(
  "list",
  nrow(step4_input)
)


for (
  i in seq_len(
    nrow(step4_input)
  )
) {

  sample_id <-
    step4_input$MH_ID[i]

  current_file <-
    step4_input$normalized_file[i]


  cat(
    "Correcting illumination for ",
    sample_id,
    "...\n",
    sep = ""
  )


  illumination_results[[i]] <- tryCatch(

    {

      ######################################################
      # Load normalized image
      ######################################################

      image <- image_read(
        current_file
      )

      image <- image_convert(
        image,
        colorspace = "sRGB"
      )


      image_info_current <-
        image_info(
          image
        )


      full_width <-
        image_info_current$width[1]

      full_height <-
        image_info_current$height[1]


      ######################################################
      # Background detection on preview
      ######################################################

      preview <- image_scale(
        image,
        preview_geometry
      )

      preview_rgb <- extract_rgb(
        preview
      )


      background <- detect_background(
        preview_rgb
      )


      ######################################################
      # Fit smooth illumination model
      ######################################################

      surface_model <-
        fit_illumination_surface(
          background$grid
        )


      ######################################################
      # Predict full-resolution illumination surface
      ######################################################

      full_surface <- predict_surface(
        surface_model,
        full_height,
        full_width
      )


      full_gain <- calculate_gain(
        full_surface
      )


      ######################################################
      # Apply illumination correction
      ######################################################

      corrected_file <- file.path(
        corrected_dir,
        paste0(
          sample_id,
          "_illumination_corrected.png"
        )
      )


      apply_illumination_correction(
        image,
        full_gain,
        corrected_file
      )


      ######################################################
      # Save background mask
      ######################################################

      mask_file <- file.path(
        mask_dir,
        paste0(
          sample_id,
          "_background_mask.png"
        )
      )


      png::writePNG(
        background$mask * 1,
        mask_file
      )


      ######################################################
      # Save illumination map
      ######################################################

      preview_height <-
        nrow(
          preview_rgb$red
        )

      preview_width <-
        ncol(
          preview_rgb$red
        )


      preview_surface <- predict_surface(
        surface_model,
        preview_height,
        preview_width
      )


      illumination_map <- (
        preview_surface / 255
      )

      illumination_map[
        illumination_map < 0
      ] <- 0

      illumination_map[
        illumination_map > 1
      ] <- 1


      illumination_map_file <- file.path(
        map_dir,
        paste0(
          sample_id,
          "_illumination_map.png"
        )
      )


      png::writePNG(
        illumination_map,
        illumination_map_file
      )


      ######################################################
      # Save local gain map
      ######################################################

      preview_gain <- calculate_gain(
        preview_surface
      )


      gain_map <- (
        preview_gain -
        gain_lower_limit
      ) /
        (
          gain_upper_limit -
          gain_lower_limit
        )


      gain_map_file <- file.path(
        map_dir,
        paste0(
          sample_id,
          "_gain_map.png"
        )
      )


      png::writePNG(
        gain_map,
        gain_map_file
      )


      ######################################################
      # Return processing information
      ######################################################

      data.frame(
        MH_ID = sample_id,

        image_role =
          step4_input$image_role[i],

        background_pixels =
          sum(
            background$mask
          ),

        grid_cells =
          nrow(
            background$grid
          ),

        corrected_file =
          corrected_file,

        correction_success =
          TRUE,

        stringsAsFactors = FALSE
      )
    },


    error = function(e) {

      warning(
        "Illumination correction failed for ",
        sample_id,
        ": ",
        conditionMessage(e)
      )


      data.frame(
        MH_ID = sample_id,

        image_role =
          step4_input$image_role[i],

        background_pixels = NA,

        grid_cells = NA,

        corrected_file = NA,

        correction_success = FALSE,

        stringsAsFactors = FALSE
      )
    }
  )
}


############################################################
# 18. Combine processing results
############################################################

illumination_correction_table <- do.call(
  rbind,
  illumination_results
)

rownames(
  illumination_correction_table
) <- NULL


############################################################
# 19. Create compact output table for subsequent steps
############################################################

illumination_corrected_image_table <-
  illumination_correction_table[
    illumination_correction_table$correction_success,
    c(
      "MH_ID",
      "image_role",
      "corrected_file"
    )
  ]


############################################################
# 20. Save processing summary
############################################################

write.csv(
  illumination_correction_table,
  file.path(
    output_dir,
    "step4_illumination_correction_summary.csv"
  ),
  row.names = FALSE
)


############################################################
# 21. Completion message
############################################################

cat(
  "\nStep 4 completed successfully.\n"
)
