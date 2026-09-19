#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# LI-600 analysis
#
# Treatment effects on PSII efficiency (PhiPSII)
#
# Comparisons:
# - Control vs. light stress
# - Control vs. salinity stress
#
# Two complementary analyses are performed:
#
# 1. Genotype-specific treatment effects
#    Separate analyses for each genotype and measurement phase.
#
#    Model:
#    phips2 ~ treatment * day + (1 | plantID)
#
# 2. Genotype-adjusted treatment effects
#    Treatment effects across all genotypes, with genotype
#    included as an additive model term.
#
#    Model:
#    phips2 ~ treatment * day_factor + genotype + (1 | plantID)
#
# Statistical analysis:
# - Type III ANOVA
# - Estimated marginal means using emmeans
#
# Note:
# The original genotype-specific PhiPSII script was no longer
# available. This section was reconstructed using the same
# analysis structure as the corresponding gsw workflow.
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(lmerTest)
  library(emmeans)
})


############################################################
# 2. Import datasets
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
# 3. Data preparation
############################################################

prepare_data <- function(data, treatment_label) {

  data %>%
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

      treatment = treatment_label,

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

      phips2 = as.numeric(
        phips2
      )
    ) %>%

    filter(
      !is.na(plantID),
      !is.na(genotype),
      !is.na(time),
      !is.na(day),
      !is.na(phips2)
    )
}


control <- prepare_data(
  control,
  "control"
)

light <- prepare_data(
  light,
  "lightstress"
)

salinity <- prepare_data(
  salinity,
  "salinitystress"
)


############################################################
# 4. Create treatment-comparison datasets
############################################################

light_data <- bind_rows(
  control,
  light
) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "control",
        "lightstress"
      )
    )
  )


salinity_data <- bind_rows(
  control,
  salinity
) %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "control",
        "salinitystress"
      )
    )
  )


############################################################
# 5. Genotype-specific treatment effects
#
# Separate models for:
# - each genotype
# - Morning
# - Evening
#
# Model:
# phips2 ~ treatment * day + (1 | plantID)
############################################################

analyze_by_genotype <- function(
  data,
  comparison_name
) {

  results <- list()

  for (geno in c("r", "h", "w")) {

    for (measurement_time in c(
      "Morning",
      "Evening"
    )) {

      df <- data %>%
        filter(
          genotype == geno,
          time == measurement_time
        ) %>%
        droplevels()


      model <- lmer(
        phips2 ~
          treatment * day +
          (1 | plantID),
        data = df
      )


      anova_result <- as.data.frame(
        anova(
          model,
          type = 3
        )
      )


      posthoc <- emmeans(
        model,
        pairwise ~ treatment
      )


      result_name <- paste(
        comparison_name,
        geno,
        measurement_time,
        sep = "_"
      )


      results[[result_name]] <- list(
        model = model,
        anova = anova_result,
        posthoc = posthoc
      )


      cat(
        "\n\n========================================\n",
        comparison_name,
        "\nGenotype: ",
        geno,
        "\nMeasurement: ",
        measurement_time,
        "\n========================================\n",
        sep = ""
      )

      cat(
        "\nType III ANOVA\n\n"
      )

      print(
        anova_result
      )

      cat(
        "\nTreatment comparison\n\n"
      )

      print(
        posthoc
      )
    }
  }

  invisible(
    results
  )
}


############################################################
# 6. Run genotype-specific analyses
############################################################

genotype_results_light <- analyze_by_genotype(
  light_data,
  "Control vs. light stress"
)

genotype_results_salinity <- analyze_by_genotype(
  salinity_data,
  "Control vs. salinity stress"
)


############################################################
# 7. Genotype-adjusted treatment effects
#
# Morning and evening are analyzed separately.
#
# Day is treated as a categorical factor.
#
# Model:
# phips2 ~ treatment * day_factor + genotype + (1 | plantID)
############################################################

analyze_overall_treatment <- function(
  data,
  comparison_name
) {

  data <- data %>%
    mutate(
      day_factor = factor(
        day,
        levels = sort(
          unique(day)
        )
      )
    )


  results <- list()


  for (measurement_time in c(
    "Morning",
    "Evening"
  )) {

    df <- data %>%
      filter(
        time == measurement_time
      ) %>%
      droplevels()


    model <- lmer(
      phips2 ~
        treatment * day_factor +
        genotype +
        (1 | plantID),
      data = df
    )


    anova_result <- as.data.frame(
      anova(
        model,
        type = 3
      )
    )


    ########################################################
    # Genotype-adjusted estimated marginal means
    ########################################################

    emm <- emmeans(
      model,
      ~ treatment | day_factor,
      weights = "equal"
    )


    emm_ci <- as.data.frame(
      confint(
        emm,
        level = 0.95
      )
    )


    results[[measurement_time]] <- list(
      model = model,
      anova = anova_result,
      emmeans = emm,
      confidence_intervals = emm_ci
    )


    cat(
      "\n\n========================================\n",
      comparison_name,
      "\nMeasurement: ",
      measurement_time,
      "\n========================================\n",
      sep = ""
    )

    cat(
      "\nType III ANOVA\n\n"
    )

    print(
      anova_result
    )
  }


  invisible(
    results
  )
}


############################################################
# 8. Run genotype-adjusted analyses
############################################################

overall_results_light <- analyze_overall_treatment(
  light_data,
  "Control vs. light stress"
)

overall_results_salinity <- analyze_overall_treatment(
  salinity_data,
  "Control vs. salinity stress"
)


############################################################
# 9. Completion message
############################################################

cat(
  "\nTreatment-effect analyses for PhiPSII completed successfully.\n"
)
