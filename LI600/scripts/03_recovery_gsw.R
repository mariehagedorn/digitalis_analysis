#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# LI-600 analysis
#
# Overnight recovery of stomatal conductance (gsw)
#
# Recovery:
# Delta gsw = Morning - Evening
#
# Conditions:
# - Control
# - Light stress
# - Salinity stress
#
# Model:
# recovery_gsw ~ genotype * day + (1 | plantID)
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
  library(lme4)
  library(lmerTest)
  library(emmeans)
})


############################################################
# 2. Import datasets
############################################################

message("Select CONTROL dataset")
control <- read_excel(file.choose())

message("Select LIGHT-STRESS dataset")
light <- read_excel(file.choose())

message("Select SALINITY-STRESS dataset")
salinity <- read_excel(file.choose())


############################################################
# 3. Recovery analysis function
############################################################

analyze_recovery_gsw <- function(data, condition_name) {

  data <- data %>%
    mutate(
      plantID = factor(plantID),
      genotype = factor(
        genotype,
        levels = c("r", "h", "w")
      ),
      treatment = factor(treatment),
      time = factor(
        time,
        levels = c("a.m.", "p.m."),
        labels = c("Morning", "Evening")
      ),
      day = as.numeric(day),
      gsw = as.numeric(gsw)
    ) %>%
    filter(
      !is.na(gsw)
    )


  ##########################################################
  # Calculate Morning - Evening recovery
  ##########################################################

  morning <- data %>%
    filter(
      time == "Morning"
    ) %>%
    select(
      plantID,
      genotype,
      treatment,
      day,
      gsw
    ) %>%
    rename(
      gsw_morning = gsw
    )

  evening <- data %>%
    filter(
      time == "Evening"
    ) %>%
    select(
      plantID,
      genotype,
      treatment,
      day,
      gsw
    ) %>%
    rename(
      gsw_evening = gsw
    )

  recovery <- morning %>%
    inner_join(
      evening,
      by = c(
        "plantID",
        "genotype",
        "treatment",
        "day"
      )
    ) %>%
    mutate(
      recovery_gsw =
        gsw_morning - gsw_evening
    )


  ##########################################################
  # Linear mixed-effects model
  ##########################################################

  model <- lmer(
    recovery_gsw ~
      genotype * day +
      (1 | plantID),
    data = recovery
  )


  ##########################################################
  # Type III ANOVA
  ##########################################################

  anova_result <- as.data.frame(
    anova(
      model,
      type = 3
    )
  )


  ##########################################################
  # Pairwise genotype comparisons
  ##########################################################

  posthoc <- emmeans(
    model,
    pairwise ~ genotype
  )


  ##########################################################
  # Output
  ##########################################################

  cat(
    "\n\n========================================\n",
    condition_name,
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
    "\nPairwise genotype comparisons\n\n"
  )

  print(
    posthoc
  )


  invisible(
    list(
      recovery = recovery,
      model = model,
      anova = anova_result,
      posthoc = posthoc
    )
  )
}


############################################################
# 4. Run analyses
############################################################

results_control <- analyze_recovery_gsw(
  control,
  "Control"
)

results_light <- analyze_recovery_gsw(
  light,
  "Light stress"
)

results_salinity <- analyze_recovery_gsw(
  salinity,
  "Salinity stress"
)


############################################################
# 5. Completion message
############################################################

cat(
  "\nRecovery analyses for gsw completed successfully.\n"
)
