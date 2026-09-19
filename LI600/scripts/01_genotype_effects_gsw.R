#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# LI-600 analysis
#
# Genotype effects on stomatal conductance (gsw)
#
# Conditions:
# - Control
# - Light stress
# - Salinity stress
#
# For each condition, morning and evening measurements are
# analyzed separately using linear mixed-effects models.
#
# Model:
# gsw ~ genotype * day + (1 | plantID)
#
# Statistical analysis:
# - Type III ANOVA
# - Pairwise genotype comparisons using emmeans
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})


############################################################
# 2. Import data
#
# Select the corresponding LI-600 data file for each
# experimental condition.
############################################################

message("Select CONTROL dataset")
control <- read_excel(
  file.choose()
)

message("Select LIGHT-STRESS dataset")
light <- read_excel(
  file.choose()
)

message("Select SALINITY-STRESS dataset")
salinity <- read_excel(
  file.choose()
)


############################################################
# 3. Data preparation function
############################################################

prepare_data <- function(data) {

  data <- data %>%
    mutate(
      plantID = factor(
        plantID
      ),

      genotype = factor(
        genotype,
        levels = c(
          "r",
          "h",
          "w"
        )
      ),

      time = factor(
        time,
        levels = c(
          "a.m.",
          "p.m."
        ),
        labels = c(
          "Morning",
          "Evening"
        )
      ),

      day = as.numeric(
        day
      ),

      gsw = as.numeric(
        gsw
      ),

      phips2 = as.numeric(
        phips2
      )
    )

  # Retain the same missing-value filtering used in the
  # original analysis scripts.
  data <- data %>%
    filter(
      !is.na(gsw),
      !is.na(phips2)
    )

  # Correct LI-600 overflow values where required.
  data$gsw <- ifelse(
    data$gsw > 10,
    data$gsw / 1000,
    data$gsw
  )

  return(
    data
  )
}


############################################################
# 4. Prepare datasets
############################################################

control <- prepare_data(
  control
)

light <- prepare_data(
  light
)

salinity <- prepare_data(
  salinity
)


############################################################
# 5. Genotype colors
############################################################

genotype_colors <- c(
  "r" = "#C00000",
  "h" = "#F7B6D2",
  "w" = "grey50"
)

genotype_labels <- c(
  "r" = "Homozygous red",
  "h" = "Heterozygous",
  "w" = "Homozygous white"
)


############################################################
# 6. Visualization function
#
# Mean gsw ± standard error is displayed across experimental
# days for all three genotypes.
############################################################

plot_condition <- function(
  data,
  condition_name
) {

  plot_data <- data %>%
    group_by(
      time,
      genotype,
      day
    ) %>%
    summarise(
      mean_gsw = mean(
        gsw,
        na.rm = TRUE
      ),

      se_gsw = sd(
        gsw,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            !is.na(gsw)
          )
        ),

      .groups = "drop"
    )

  p <- ggplot(
    plot_data,
    aes(
      x = day,
      y = mean_gsw,
      color = genotype,
      fill = genotype,
      group = genotype
    )
  ) +

    geom_ribbon(
      aes(
        ymin = mean_gsw - se_gsw,
        ymax = mean_gsw + se_gsw
      ),
      alpha = 0.25,
      color = NA
    ) +

    geom_line(
      linewidth = 1.2
    ) +

    geom_point(
      size = 2.5
    ) +

    geom_errorbar(
      aes(
        ymin = mean_gsw - se_gsw,
        ymax = mean_gsw + se_gsw
      ),
      width = 0.2
    ) +

    facet_wrap(
      ~ time,
      ncol = 1
    ) +

    scale_color_manual(
      values = genotype_colors,
      labels = genotype_labels,
      name = "Genotype"
    ) +

    scale_fill_manual(
      values = genotype_colors,
      labels = genotype_labels,
      name = "Genotype"
    ) +

    scale_x_continuous(
      breaks = sort(
        unique(
          data$day
        )
      )
    ) +

    labs(
      title = paste(
        condition_name,
        "- stomatal conductance"
      ),
      x = "Day",
      y = expression(
        g[sw]~"(mol m"^{-2}~s^{-1}*")"
      )
    ) +

    theme_classic(
      base_size = 14
    ) +

    theme(
      legend.position = "top",
      strip.text = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        color = "black"
      )
    )

  return(
    p
  )
}


############################################################
# 7. Statistical analysis function
#
# Morning and evening measurements are analyzed separately.
#
# Fixed effects:
# - genotype
# - day
# - genotype × day
#
# Random effect:
# - plantID
############################################################

analyze_genotype_effects <- function(
  data,
  condition_name
) {

  cat(
    "\n\n========================================\n"
  )

  cat(
    condition_name,
    "\n"
  )

  cat(
    "========================================\n"
  )


  ##########################################################
  # Morning
  ##########################################################

  data_morning <- data %>%
    filter(
      time == "Morning"
    )

  model_morning <- lmer(
    gsw ~
      genotype * day +
      (1 | plantID),
    data = data_morning
  )

  cat(
    "\nMorning - Type III ANOVA\n\n"
  )

  anova_morning <- as.data.frame(
    anova(
      model_morning,
      type = 3
    )
  )

  print(
    anova_morning
  )

  cat(
    "\nMorning - pairwise genotype comparisons\n\n"
  )

  posthoc_morning <- emmeans(
    model_morning,
    pairwise ~ genotype
  )

  print(
    posthoc_morning
  )


  ##########################################################
  # Evening
  ##########################################################

  data_evening <- data %>%
    filter(
      time == "Evening"
    )

  model_evening <- lmer(
    gsw ~
      genotype * day +
      (1 | plantID),
    data = data_evening
  )

  cat(
    "\nEvening - Type III ANOVA\n\n"
  )

  anova_evening <- as.data.frame(
    anova(
      model_evening,
      type = 3
    )
  )

  print(
    anova_evening
  )

  cat(
    "\nEvening - pairwise genotype comparisons\n\n"
  )

  posthoc_evening <- emmeans(
    model_evening,
    pairwise ~ genotype
  )

  print(
    posthoc_evening
  )


  ##########################################################
  # Return analysis objects
  ##########################################################

  invisible(
    list(
      morning_model =
        model_morning,

      morning_anova =
        anova_morning,

      morning_posthoc =
        posthoc_morning,

      evening_model =
        model_evening,

      evening_anova =
        anova_evening,

      evening_posthoc =
        posthoc_evening
    )
  )
}


############################################################
# 8. Create exploratory plots
############################################################

plot_control <- plot_condition(
  control,
  "Control"
)

plot_light <- plot_condition(
  light,
  "Light stress"
)

plot_salinity <- plot_condition(
  salinity,
  "Salinity stress"
)

print(
  plot_control
)

print(
  plot_light
)

print(
  plot_salinity
)


############################################################
# 9. Run genotype-effect analyses
############################################################

results_control <- analyze_genotype_effects(
  control,
  "Control"
)

results_light <- analyze_genotype_effects(
  light,
  "Light stress"
)

results_salinity <- analyze_genotype_effects(
  salinity,
  "Salinity stress"
)


############################################################
# 10. Completion message
############################################################

cat(
  "\nGenotype-effect analyses for gsw completed successfully.\n"
)
