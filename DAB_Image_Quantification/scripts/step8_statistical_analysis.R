#!/usr/bin/env Rscript


############################################################
# Digitalis purpurea
# DAB image quantification workflow
#
# Step 8: Statistical Analysis
#
# Purpose:
# - Import leaf-level DAB burden values generated in Step 7
# - Assign genotype and treatment information
# - Test genotype, treatment, and genotype × treatment effects
# - Compare genotypes within each treatment
# - Compare Light stress and Salinity stress with Control
# - Evaluate basic model diagnostics
#
# Biological unit:
# One leaf from one plant.
# Individual pixels are not treated as biological replicates.
#
# Raw data and generated results are not included in this
# repository.
############################################################


############################################################
# 1. Load packages
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(car)
  library(sandwich)
})


############################################################
# 2. Select Step-7 result file
############################################################

cat(
  "Select the Step-7 handoff RDS file.\n"
)

input_file <- file.choose()

input_object <- readRDS(
  input_file
)


############################################################
# 3. Extract leaf-level results
############################################################

if (
  !is.null(
    input_object$analysis_data
  )
) {

  analysis_data <-
    input_object$analysis_data

} else if (
  !is.null(
    input_object$leaf_level_results
  )
) {

  analysis_data <-
    input_object$leaf_level_results

} else if (
  !is.null(
    input_object$results
  )
) {

  analysis_data <-
    input_object$results

} else {

  stop(
    paste(
      "Input RDS must contain",
      "'analysis_data', 'leaf_level_results',",
      "or 'results'."
    )
  )
}


if (!is.data.frame(analysis_data)) {

  stop(
    "Leaf-level results were not stored as a data frame."
  )
}


############################################################
# 4. Identify primary DAB-burden variable
#
# The original analysis used:
# combined_normalized_total_burden
#
# The simplified Step-7 documentation script stores the same
# primary metric under the shorter name:
# DAB_burden
############################################################

if (
  "combined_normalized_total_burden" %in%
  names(analysis_data)
) {

  primary_endpoint <-
    "combined_normalized_total_burden"

} else if (
  "DAB_burden" %in%
  names(analysis_data)
) {

  primary_endpoint <-
    "DAB_burden"

} else {

  stop(
    "No recognized DAB-burden variable was found."
  )
}


############################################################
# 5. Exclude unsuccessful rows and negative control
############################################################

if (
  "success" %in%
  names(analysis_data)
) {

  analysis_data <- analysis_data %>%
    filter(
      success %in% TRUE
    )
}


if (
  "is_negative_control" %in%
  names(analysis_data)
) {

  analysis_data <- analysis_data %>%
    filter(
      !(is_negative_control %in% TRUE)
    )
}


if (
  "image_role" %in%
  names(analysis_data)
) {

  analysis_data <- analysis_data %>%
    filter(
      !tolower(
        as.character(image_role)
      ) %in%
        c(
          "negative_control",
          "negative control",
          "nc"
        )
    )
}


############################################################
# 6. Assign plant number, genotype, and treatment
############################################################

analysis_data <- analysis_data %>%
  mutate(

    MH_ID =
      trimws(
        as.character(MH_ID)
      ),

    MH_number =
      suppressWarnings(
        as.integer(
          sub(
            "^MH0*",
            "",
            toupper(MH_ID)
          )
        )
      ),


    genotype = case_when(

      between(
        MH_number,
        1L,
        27L
      ) ~ "Homozygous red",

      between(
        MH_number,
        28L,
        54L
      ) ~ "Homozygous white",

      between(
        MH_number,
        55L,
        81L
      ) ~ "Heterozygous",

      TRUE ~ NA_character_
    ),


    treatment = case_when(

      between(
        MH_number,
        1L,
        9L
      ) ~ "Control",

      between(
        MH_number,
        10L,
        18L
      ) ~ "Light stress",

      between(
        MH_number,
        19L,
        27L
      ) ~ "Salinity stress",


      between(
        MH_number,
        28L,
        36L
      ) ~ "Control",

      between(
        MH_number,
        37L,
        45L
      ) ~ "Light stress",

      between(
        MH_number,
        46L,
        54L
      ) ~ "Salinity stress",


      between(
        MH_number,
        55L,
        63L
      ) ~ "Control",

      between(
        MH_number,
        64L,
        72L
      ) ~ "Light stress",

      between(
        MH_number,
        73L,
        81L
      ) ~ "Salinity stress",

      TRUE ~ NA_character_
    )
  )


if (
  anyNA(
    analysis_data$genotype
  ) ||
  anyNA(
    analysis_data$treatment
  )
) {

  stop(
    paste(
      "At least one biological sample could not be",
      "assigned to a genotype or treatment."
    )
  )
}


############################################################
# 7. Prepare factors
############################################################

analysis_data <- analysis_data %>%
  mutate(

    genotype = factor(
      genotype,
      levels = c(
        "Homozygous red",
        "Homozygous white",
        "Heterozygous"
      )
    ),

    treatment = factor(
      treatment,
      levels = c(
        "Control",
        "Light stress",
        "Salinity stress"
      )
    )

  ) %>%
  filter(
    is.finite(
      .data[[primary_endpoint]]
    )
  )


############################################################
# 8. Verify experimental groups
############################################################

group_counts <- analysis_data %>%
  count(
    genotype,
    treatment,
    name = "n"
  )


if (
  nrow(group_counts) != 9 ||
  any(
    group_counts$n == 0
  )
) {

  stop(
    paste(
      "At least one genotype × treatment group",
      "is missing."
    )
  )
}


############################################################
# 9. Use sum-to-zero contrasts for Type-III tests
############################################################

old_contrasts <- options(
  "contrasts"
)

options(
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
)


############################################################
# 10. Fit primary factorial model
############################################################

primary_formula <- as.formula(
  paste0(
    "`",
    primary_endpoint,
    "` ~ genotype * treatment"
  )
)


primary_model <- lm(
  primary_formula,
  data = analysis_data,
  na.action = na.exclude
)


############################################################
# 11. HC3-robust covariance matrix
############################################################

robust_vcov <- sandwich::vcovHC(
  primary_model,
  type = "HC3"
)


############################################################
# 12. HC3-robust Type-III tests
#
# Tests:
# - genotype
# - treatment
# - genotype × treatment
############################################################

anova_HC3 <- car::Anova(
  primary_model,
  type = 3,
  white.adjust = "hc3"
)


cat(
  "\nHC3-robust Type-III ANOVA:\n"
)

print(
  anova_HC3
)


############################################################
# 13. Genotype comparisons within each treatment
#
# All three genotype pairs are compared separately within
# Control, Light stress, and Salinity stress.
#
# Tukey adjustment controls the family-wise error rate for
# the pairwise genotype comparisons.
############################################################

genotype_emmeans <- emmeans(
  primary_model,
  ~ genotype | treatment,
  vcov. = robust_vcov
)


genotype_comparisons <- as.data.frame(
  summary(
    pairs(
      genotype_emmeans,
      adjust = "tukey"
    ),
    infer = c(
      TRUE,
      TRUE
    )
  )
)


cat(
  "\nGenotype comparisons within treatments:\n"
)

print(
  genotype_comparisons
)


############################################################
# 14. Overall treatment means
#
# Genotypes are weighted equally so that differences in
# sample size do not give one genotype greater weight.
############################################################

treatment_emmeans <- emmeans(
  primary_model,
  ~ treatment,
  vcov. = robust_vcov,
  weights = "equal"
)


treatment_means <- as.data.frame(
  summary(
    treatment_emmeans,
    infer = c(
      TRUE,
      TRUE
    )
  )
)


############################################################
# 15. Planned stress-versus-Control comparisons
#
# Two comparisons:
# - Light stress vs Control
# - Salinity stress vs Control
#
# Holm adjustment is applied across the two planned tests.
############################################################

stress_vs_control <- as.data.frame(
  summary(
    contrast(
      treatment_emmeans,

      method = list(

        "Light stress - Control" =
          c(
            -1,
            1,
            0
          ),

        "Salinity stress - Control" =
          c(
            -1,
            0,
            1
          )
      ),

      adjust = "holm"
    ),

    infer = c(
      TRUE,
      TRUE
    )
  )
)


cat(
  "\nStress-versus-Control comparisons:\n"
)

print(
  stress_vs_control
)


############################################################
# 16. Genotype × treatment estimated means
############################################################

cell_means <- as.data.frame(
  summary(
    emmeans(
      primary_model,
      ~ genotype * treatment,
      vcov. = robust_vcov
    ),
    infer = c(
      TRUE,
      TRUE
    )
  )
)


############################################################
# 17. Model diagnostics
############################################################

model_residuals <- residuals(
  primary_model
)

standardized_residuals <- rstandard(
  primary_model
)


############################################################
# 17A. Residual normality
############################################################

shapiro_result <- shapiro.test(
  model_residuals
)


############################################################
# 17B. Variance homogeneity
############################################################

levene_data <- analysis_data %>%
  mutate(
    experimental_group = interaction(
      genotype,
      treatment,
      drop = TRUE
    )
  )


levene_formula <- as.formula(
  paste0(
    "`",
    primary_endpoint,
    "` ~ experimental_group"
  )
)


levene_result <- car::leveneTest(
  levene_formula,
  data = levene_data,
  center = median
)


############################################################
# 17C. Influential observations
############################################################

cooks_distance <- cooks.distance(
  primary_model
)


cook_threshold <- (
  4 /
  nrow(
    analysis_data
  )
)


cook_table <- data.frame(

  MH_ID =
    analysis_data$MH_ID,

  genotype =
    analysis_data$genotype,

  treatment =
    analysis_data$treatment,

  cooks_distance =
    as.numeric(
      cooks_distance
    ),

  above_4_over_n =
    as.numeric(
      cooks_distance
    ) >
    cook_threshold
)


cat(
  "\nShapiro-Wilk residual test:\n"
)

print(
  shapiro_result
)


cat(
  "\nLevene test:\n"
)

print(
  levene_result
)


cat(
  "\nCook's-distance screening threshold: ",
  cook_threshold,
  "\n",
  sep = ""
)


############################################################
# 18. Leaf-level primary endpoint table
############################################################

DAB_burden_per_leaf <- analysis_data %>%
  transmute(

    MH_ID =
      MH_ID,

    genotype =
      genotype,

    treatment =
      treatment,

    DAB_burden =
      as.numeric(
        .data[[primary_endpoint]]
      )
  )


############################################################
# 19. Define output directory
############################################################

output_dir <- file.path(
  "DAB-Image_Quantification",
  "results",
  "step8_statistical_analysis"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 20. Save statistical tables
############################################################

write.csv(
  DAB_burden_per_leaf,
  file.path(
    output_dir,
    "DAB_burden_per_leaf.csv"
  ),
  row.names = FALSE
)


write.csv(
  group_counts,
  file.path(
    output_dir,
    "group_counts.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    anova_HC3
  ),
  file.path(
    output_dir,
    "HC3_type3_ANOVA.csv"
  )
)


write.csv(
  genotype_comparisons,
  file.path(
    output_dir,
    "genotype_within_treatment.csv"
  ),
  row.names = FALSE
)


write.csv(
  treatment_means,
  file.path(
    output_dir,
    "treatment_marginal_means.csv"
  ),
  row.names = FALSE
)


write.csv(
  stress_vs_control,
  file.path(
    output_dir,
    "stress_vs_control.csv"
  ),
  row.names = FALSE
)


write.csv(
  cell_means,
  file.path(
    output_dir,
    "genotype_treatment_means.csv"
  ),
  row.names = FALSE
)


write.csv(
  cook_table,
  file.path(
    output_dir,
    "cooks_distance.csv"
  ),
  row.names = FALSE
)


############################################################
# 21. Save complete statistical object
############################################################

saveRDS(

  list(

    primary_endpoint =
      primary_endpoint,

    analysis_data =
      analysis_data,

    model =
      primary_model,

    robust_vcov =
      robust_vcov,

    anova_HC3 =
      anova_HC3,

    genotype_comparisons =
      genotype_comparisons,

    treatment_means =
      treatment_means,

    stress_vs_control =
      stress_vs_control,

    genotype_treatment_means =
      cell_means,

    shapiro =
      shapiro_result,

    levene =
      levene_result,

    cooks_distance =
      cook_table
  ),

  file.path(
    output_dir,
    "step8_statistical_results.rds"
  ),

  compress = "gzip"
)


############################################################
# 22. Restore contrast settings
############################################################

options(
  contrasts =
    old_contrasts$contrasts
)


############################################################
# 23. Completion message
############################################################

cat(
  "\nStep 8 completed successfully.\n"
)
