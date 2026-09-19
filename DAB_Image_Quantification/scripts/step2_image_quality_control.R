#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 2: Technical Image Quality Control
#
# Purpose:
# - Load ruler and ruler-free images identified in Step 1
# - Record technical image properties
# - Calculate RGB channel statistics
# - Generate diagnostic channel panels and histograms
# - Export QC summaries
#
# Raw images and generated QC outputs are not included in
# this repository.
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(magick)
  library(openxlsx)
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
# 3. Prepare image input table
############################################################

ruler_input <- data.frame(
  MH_ID = sample_overview$MH_ID,
  image_type = "ruler",
  file_path = sample_overview$ruler_file,
  stringsAsFactors = FALSE
)

decolorized_input <- data.frame(
  MH_ID = sample_overview$MH_ID,
  image_type = "ruler_free",
  file_path = sample_overview$decolorized_file,
  stringsAsFactors = FALSE
)

image_input <- rbind(
  ruler_input,
  decolorized_input
)

image_input <- image_input[
  order(
    image_input$MH_ID,
    image_input$image_type
  ),
]

rownames(image_input) <- NULL


############################################################
# 4. Define output directory
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step2_image_quality_control"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

panel_dir <- file.path(
  output_dir,
  "channel_panels"
)

histogram_dir <- file.path(
  output_dir,
  "histograms"
)

dir.create(
  panel_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  histogram_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 5. Helper functions
############################################################

channel_statistics <- function(values) {

  data.frame(
    mean = mean(
      values,
      na.rm = TRUE
    ),
    median = median(
      values,
      na.rm = TRUE
    ),
    sd = sd(
      values,
      na.rm = TRUE
    )
  )
}


add_prefix <- function(
  x,
  prefix
) {

  colnames(x) <- paste0(
    prefix,
    "_",
    colnames(x)
  )

  x
}


create_qc_panel <- function(
  image,
  sample_id,
  image_type,
  output_path
) {

  preview <- image_scale(
    image,
    "1400x1400>"
  )

  red <- image_channel(
    preview,
    "red"
  )

  green <- image_channel(
    preview,
    "green"
  )

  blue <- image_channel(
    preview,
    "blue"
  )

  label <- paste(
    sample_id,
    image_type,
    sep = " - "
  )

  images <- c(
    image_annotate(
      preview,
      paste(label, "Original"),
      gravity = "northwest",
      location = "+20+20",
      size = 30,
      boxcolor = "white"
    ),

    image_annotate(
      red,
      paste(label, "Red"),
      gravity = "northwest",
      location = "+20+20",
      size = 30,
      boxcolor = "white"
    ),

    image_annotate(
      green,
      paste(label, "Green"),
      gravity = "northwest",
      location = "+20+20",
      size = 30,
      boxcolor = "white"
    ),

    image_annotate(
      blue,
      paste(label, "Blue"),
      gravity = "northwest",
      location = "+20+20",
      size = 30,
      boxcolor = "white"
    )
  )

  montage <- image_montage(
    images,
    tile = "2x2",
    geometry = "+20+20"
  )

  image_write(
    montage,
    path = output_path,
    format = "png"
  )
}


create_rgb_histograms <- function(
  red,
  green,
  blue,
  sample_id,
  image_type,
  output_path
) {

  png(
    output_path,
    width = 1800,
    height = 1400,
    res = 180
  )

  par(
    mfrow = c(3, 1),
    mar = c(4, 4, 3, 1)
  )

  hist(
    red,
    breaks = seq(-0.5, 255.5, 1),
    main = paste(
      sample_id,
      image_type,
      "- Red channel"
    ),
    xlab = "Pixel intensity",
    xlim = c(0, 255)
  )

  hist(
    green,
    breaks = seq(-0.5, 255.5, 1),
    main = paste(
      sample_id,
      image_type,
      "- Green channel"
    ),
    xlab = "Pixel intensity",
    xlim = c(0, 255)
  )

  hist(
    blue,
    breaks = seq(-0.5, 255.5, 1),
    main = paste(
      sample_id,
      image_type,
      "- Blue channel"
    ),
    xlab = "Pixel intensity",
    xlim = c(0, 255)
  )

  dev.off()
}


############################################################
# 6. Analyze images
############################################################

qc_results <- vector(
  "list",
  nrow(image_input)
)


for (i in seq_len(nrow(image_input))) {

  sample_id <- image_input$MH_ID[i]
  image_type <- image_input$image_type[i]
  current_file <- image_input$file_path[i]

  if (
    is.na(current_file) ||
    !file.exists(current_file)
  ) {

    qc_results[[i]] <- data.frame(
      MH_ID = sample_id,
      image_type = image_type,
      width_px = NA,
      height_px = NA,
      total_pixels = NA,
      format = NA,
      colorspace = NA,
      red_mean = NA,
      red_median = NA,
      red_sd = NA,
      green_mean = NA,
      green_median = NA,
      green_sd = NA,
      blue_mean = NA,
      blue_median = NA,
      blue_sd = NA,
      read_success = FALSE
    )

    next
  }


  qc_results[[i]] <- tryCatch(

    {

      image <- image_read(
        current_file
      )

      info <- image_info(
        image
      )

      pixels <- image_data(
        image,
        channels = "rgb"
      )

      red_values <- as.integer(
        pixels[1, , ]
      )

      green_values <- as.integer(
        pixels[2, , ]
      )

      blue_values <- as.integer(
        pixels[3, , ]
      )


      red_stats <- add_prefix(
        channel_statistics(red_values),
        "red"
      )

      green_stats <- add_prefix(
        channel_statistics(green_values),
        "green"
      )

      blue_stats <- add_prefix(
        channel_statistics(blue_values),
        "blue"
      )


      create_qc_panel(
        image,
        sample_id,
        image_type,
        file.path(
          panel_dir,
          paste0(
            sample_id,
            "_",
            image_type,
            "_channels.png"
          )
        )
      )


      create_rgb_histograms(
        red_values,
        green_values,
        blue_values,
        sample_id,
        image_type,
        file.path(
          histogram_dir,
          paste0(
            sample_id,
            "_",
            image_type,
            "_histograms.png"
          )
        )
      )


      data.frame(
        MH_ID = sample_id,
        image_type = image_type,

        width_px = info$width[1],
        height_px = info$height[1],

        total_pixels =
          info$width[1] *
          info$height[1],

        format =
          as.character(
            info$format[1]
          ),

        colorspace =
          as.character(
            info$colorspace[1]
          ),

        red_stats,
        green_stats,
        blue_stats,

        read_success = TRUE,

        stringsAsFactors = FALSE
      )
    },


    error = function(e) {

      warning(
        "Image processing failed for ",
        sample_id,
        ": ",
        conditionMessage(e)
      )

      data.frame(
        MH_ID = sample_id,
        image_type = image_type,
        width_px = NA,
        height_px = NA,
        total_pixels = NA,
        format = NA,
        colorspace = NA,
        red_mean = NA,
        red_median = NA,
        red_sd = NA,
        green_mean = NA,
        green_median = NA,
        green_sd = NA,
        blue_mean = NA,
        blue_median = NA,
        blue_sd = NA,
        read_success = FALSE
      )
    }
  )
}


############################################################
# 7. Combine QC results
############################################################

qc_table <- do.call(
  rbind,
  qc_results
)

rownames(qc_table) <- NULL


############################################################
# 8. Export QC table
############################################################

excel_output <- file.path(
  output_dir,
  "step2_image_quality_control.xlsx"
)

write.xlsx(
  qc_table,
  excel_output,
  overwrite = TRUE
)


############################################################
# 9. Export text summary
############################################################

text_output <- file.path(
  output_dir,
  "step2_image_quality_control.txt"
)

summary_lines <- c(
  "DAB Image Quantification Workflow",
  "Step 2: Technical Image Quality Control",
  "",
  paste(
    "Images evaluated:",
    nrow(qc_table)
  ),
  paste(
    "Successfully processed:",
    sum(qc_table$read_success)
  ),
  paste(
    "Processing failures:",
    sum(!qc_table$read_success)
  )
)

writeLines(
  summary_lines,
  text_output
)


############################################################
# 10. Completion message
############################################################

cat(
  "\nStep 2 completed successfully.\n"
)
