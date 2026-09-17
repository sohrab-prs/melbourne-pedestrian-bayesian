# ============================================================
# 20_frequentist_comparison.R
# Melbourne Pedestrian Bayesian Project
#
# Bonus frequentist comparison to Bayesian M4.
#
# Frequentist analogue:
#   Negative Binomial mixed model
#   + sensor random intercept
#   + date random intercept
#
# mgcv::bam() is used for computational efficiency.
#
# Parametric coefficient inference is extracted directly from
# coef() and vcov(), avoiding summary.gam() because the fitted
# NB BAM generated numerical deviance-residual warnings inside
# summary().
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(mgcv)
library(ggplot2)

# ============================================================
# PART A — DATA
# ============================================================

dat <- read_csv(
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  ),
  col_types = cols(
    date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

stopifnot(
  nrow(dat) == 18671,
  n_distinct(dat$sensor_index) == 53,
  n_distinct(dat$date) == 364,
  sum(is.na(dat)) == 0
)

dat <- dat %>%
  mutate(
    
    day_of_week =
      factor(
        dow_index,
        levels = 1:7,
        labels = c(
          "Monday",
          "Tuesday",
          "Wednesday",
          "Thursday",
          "Friday",
          "Saturday",
          "Sunday"
        )
      ),
    
    day_of_week =
      relevel(
        day_of_week,
        ref = "Monday"
      ),
    
    sensor_factor =
      factor(
        sensor_index
      ),
    
    date_factor =
      factor(
        date
      )
  )

# ============================================================
# PART B — FREQUENTIST MODEL
# ============================================================

cat(
  "\nFitting frequentist Negative Binomial mixed model...\n"
)

n_threads <- min(
  4L,
  max(
    1L,
    parallel::detectCores(
      logical = FALSE
    )
  )
)

freq_model <- bam(
  
  daily_ped_count ~
    
    day_of_week +
    
    public_holiday_num +
    
    z_mean_temp +
    
    z_log1p_precipitation +
    
    z_business_count +
    
    z_resident_population +
    
    z_built_coverage +
    
    z_structure_height +
    
    s(
      sensor_factor,
      bs = "re"
    ) +
    
    s(
      date_factor,
      bs = "re"
    ),
  
  data =
    dat,
  
  family =
    nb(
      link = "log"
    ),
  
  method =
    "fREML",
  
  discrete =
    TRUE,
  
  nthreads =
    n_threads
)

cat(
  "Frequentist model completed.\n"
)

# ------------------------------------------------------------
# Basic numerical QA
# ------------------------------------------------------------

fitted_mu <- fitted(
  freq_model
)

stopifnot(
  length(fitted_mu) == nrow(dat),
  all(is.finite(fitted_mu)),
  all(fitted_mu > 0)
)

model_converged <-
  if (
    !is.null(
      freq_model$converged
    )
  ) {
    isTRUE(
      freq_model$converged
    )
  } else {
    NA
  }

cat(
  "Finite positive fitted means: PASS\n"
)

cat(
  "Fitted mean range:",
  min(fitted_mu),
  "to",
  max(fitted_mu),
  "\n"
)

cat(
  "Model convergence flag:",
  model_converged,
  "\n"
)

# ============================================================
# PART C — PARAMETRIC COEFFICIENTS
# ============================================================

# Mapping between mgcv coefficient names and Bayesian M4
# parameter names.

parameter_map <- c(
  
  "day_of_weekTuesday" =
    "beta_dow[2]",
  
  "day_of_weekWednesday" =
    "beta_dow[3]",
  
  "day_of_weekThursday" =
    "beta_dow[4]",
  
  "day_of_weekFriday" =
    "beta_dow[5]",
  
  "day_of_weekSaturday" =
    "beta_dow[6]",
  
  "day_of_weekSunday" =
    "beta_dow[7]",
  
  "public_holiday_num" =
    "beta_holiday",
  
  "z_mean_temp" =
    "beta_temp",
  
  "z_log1p_precipitation" =
    "beta_rain",
  
  "z_business_count" =
    "beta_business",
  
  "z_resident_population" =
    "beta_population",
  
  "z_built_coverage" =
    "beta_coverage",
  
  "z_structure_height" =
    "beta_height"
)

# ------------------------------------------------------------
# Extract estimates directly from fitted BAM.
# This avoids summary.gam(), which produced the numerical
# deviance-residual warning in this model.
# ------------------------------------------------------------

all_coefficients <-
  coef(
    freq_model
  )

all_vcov <-
  vcov(
    freq_model
  )

required_parametric_names <- c(
  "(Intercept)",
  names(
    parameter_map
  )
)

missing_parameters <- setdiff(
  required_parametric_names,
  names(
    all_coefficients
  )
)

if (
  length(
    missing_parameters
  ) > 0
) {
  
  stop(
    paste(
      "Expected frequentist coefficients were not found:",
      paste(
        missing_parameters,
        collapse = ", "
      )
    )
  )
}

parametric_estimates <-
  all_coefficients[
    required_parametric_names
  ]

parametric_se <-
  sqrt(
    diag(
      all_vcov
    )[
      required_parametric_names
    ]
  )

stopifnot(
  all(
    is.finite(
      parametric_estimates
    )
  ),
  all(
    is.finite(
      parametric_se
    )
  ),
  all(
    parametric_se > 0
  )
)

parametric_z <-
  parametric_estimates /
  parametric_se

parametric_p <-
  2 *
  pnorm(
    abs(
      parametric_z
    ),
    lower.tail = FALSE
  )

p_table <- tibble(
  
  frequentist_parameter =
    required_parametric_names,
  
  estimate =
    as.numeric(
      parametric_estimates
    ),
  
  standard_error =
    as.numeric(
      parametric_se
    ),
  
  z_value =
    as.numeric(
      parametric_z
    ),
  
  p_value =
    as.numeric(
      parametric_p
    )
)

freq_effects <- p_table %>%
  
  filter(
    frequentist_parameter !=
      "(Intercept)"
  ) %>%
  
  mutate(
    
    bayesian_parameter =
      unname(
        parameter_map[
          frequentist_parameter
        ]
      ),
    
    frequentist_rate_ratio =
      exp(
        estimate
      ),
    
    frequentist_percent_change =
      100 *
      (
        exp(
          estimate
        ) -
          1
      ),
    
    frequentist_ci_q025 =
      estimate -
      1.96 *
      standard_error,
    
    frequentist_ci_q975 =
      estimate +
      1.96 *
      standard_error,
    
    frequentist_percent_q025 =
      100 *
      (
        exp(
          frequentist_ci_q025
        ) -
          1
      ),
    
    frequentist_percent_q975 =
      100 *
      (
        exp(
          frequentist_ci_q975
        ) -
          1
      )
  )

stopifnot(
  all(
    !is.na(
      freq_effects$bayesian_parameter
    )
  )
)

# ============================================================
# PART D — BAYESIAN M4 RESULTS
# ============================================================

bayes_effects <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_final_effect_estimates.csv"
  ),
  show_col_types = FALSE
) %>%
  
  select(
    parameter,
    effect,
    percent_change_median,
    percent_change_q025,
    percent_change_q975,
    directional_probability
  )

comparison <- freq_effects %>%
  
  left_join(
    bayes_effects,
    by = c(
      "bayesian_parameter" =
        "parameter"
    )
  ) %>%
  
  transmute(
    
    parameter =
      bayesian_parameter,
    
    effect =
      effect,
    
    bayesian_percent_change =
      percent_change_median,
    
    bayesian_q025 =
      percent_change_q025,
    
    bayesian_q975 =
      percent_change_q975,
    
    bayesian_directional_probability =
      directional_probability,
    
    frequentist_percent_change =
      frequentist_percent_change,
    
    frequentist_q025 =
      frequentist_percent_q025,
    
    frequentist_q975 =
      frequentist_percent_q975,
    
    frequentist_p_value =
      p_value,
    
    same_direction =
      sign(
        bayesian_percent_change
      ) ==
      sign(
        frequentist_percent_change
      )
  )

stopifnot(
  nrow(comparison) == 13,
  all(!is.na(comparison$effect))
)

write_csv(
  comparison,
  here(
    "outputs",
    "tables",
    "bayesian_frequentist_comparison.csv"
  )
)

# ============================================================
# PART E — NEGATIVE BINOMIAL DISPERSION
# ============================================================

freq_theta <-
  freq_model$family$getTheta(
    TRUE
  )

stopifnot(
  length(freq_theta) == 1,
  is.finite(freq_theta),
  freq_theta > 0
)

bayes_variance <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_variance_parameters.csv"
  ),
  show_col_types = FALSE
)

bayes_phi <- bayes_variance %>%
  filter(
    parameter == "phi"
  )

stopifnot(
  nrow(bayes_phi) == 1
)

dispersion_comparison <- tibble(
  
  quantity =
    "Negative Binomial dispersion",
  
  bayesian_median =
    bayes_phi$median,
  
  bayesian_q025 =
    bayes_phi$q025,
  
  bayesian_q975 =
    bayes_phi$q975,
  
  frequentist_estimate =
    freq_theta
)

write_csv(
  dispersion_comparison,
  here(
    "outputs",
    "tables",
    "bayesian_frequentist_dispersion.csv"
  )
)

# ============================================================
# PART F — RANDOM-EFFECT VARIANCE COMPONENTS
# ============================================================

# gam.vcomp() may return:
#   - a numeric vector,
#   - a matrix, or
#   - a list containing the variance-component object in $vc.
#
# For likelihood-based smoothness estimation, $vc commonly
# contains standard deviations plus confidence limits.

variance_components_raw <-
  suppressWarnings(
    gam.vcomp(
      freq_model
    )
  )

if (
  is.list(
    variance_components_raw
  ) &&
  !is.null(
    variance_components_raw$vc
  )
) {
  
  vc <-
    variance_components_raw$vc
  
  if (
    is.matrix(vc) ||
    is.data.frame(vc)
  ) {
    
    variance_components <-
      as.data.frame(
        vc
      ) %>%
      rownames_to_column(
        "component"
      )
    
  } else if (
    is.numeric(vc)
  ) {
    
    variance_components <-
      tibble(
        component =
          names(vc),
        
        std.dev =
          as.numeric(vc)
      )
    
  } else {
    
    stop(
      "gam.vcomp() returned an unsupported $vc object."
    )
  }
  
} else if (
  is.matrix(
    variance_components_raw
  ) ||
  is.data.frame(
    variance_components_raw
  )
) {
  
  variance_components <-
    as.data.frame(
      variance_components_raw
    ) %>%
    rownames_to_column(
      "component"
    )
  
} else if (
  is.numeric(
    variance_components_raw
  )
) {
  
  variance_components <-
    tibble(
      component =
        names(
          variance_components_raw
        ),
      
      std.dev =
        as.numeric(
          variance_components_raw
        )
    )
  
} else {
  
  stop(
    "gam.vcomp() returned an unsupported object."
  )
}

write_csv(
  variance_components,
  here(
    "outputs",
    "tables",
    "frequentist_variance_components.csv"
  )
)

# ============================================================
# PART G — COMPARISON FIGURE
# ============================================================

# Construct long plotting data explicitly rather than relying
# on pivot parsing.

plot_data <- bind_rows(
  
  comparison %>%
    transmute(
      
      effect =
        effect,
      
      method =
        "Bayesian M4",
      
      estimate =
        bayesian_percent_change,
      
      q025 =
        bayesian_q025,
      
      q975 =
        bayesian_q975
    ),
  
  comparison %>%
    transmute(
      
      effect =
        effect,
      
      method =
        "Frequentist NB mixed model",
      
      estimate =
        frequentist_percent_change,
      
      q025 =
        frequentist_q025,
      
      q975 =
        frequentist_q975
    )
) %>%
  
  mutate(
    
    effect =
      factor(
        effect,
        levels =
          rev(
            unique(
              comparison$effect
            )
          )
      ),
    
    method =
      factor(
        method,
        levels = c(
          "Bayesian M4",
          "Frequentist NB mixed model"
        )
      )
  )

p <- ggplot(
  plot_data,
  aes(
    x =
      estimate,
    y =
      effect,
    shape =
      method
  )
) +
  
  geom_vline(
    xintercept =
      0,
    linetype =
      "dashed"
  ) +
  
  geom_errorbar(
    aes(
      xmin =
        q025,
      xmax =
        q975
    ),
    orientation =
      "y",
    position =
      position_dodge(
        width =
          0.5
      ),
    width =
      0.15
  ) +
  
  geom_point(
    position =
      position_dodge(
        width =
          0.5
      ),
    size =
      2
  ) +
  
  labs(
    title =
      "Bayesian and Frequentist Effect Estimates",
    subtitle =
      "Percentage change in expected daily pedestrian activity",
    x =
      "Percentage change",
    y =
      NULL,
    shape =
      "Method"
  ) +
  
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "bayesian_frequentist_comparison.png"
  ),
  p,
  width =
    10,
  height =
    9,
  dpi =
    300
)

# ============================================================
# PART H — QA
# ============================================================

cat(
  "\n========== STEP 20 FREQUENTIST COMPARISON ==========\n"
)

cat(
  "\nModel numerical QA:\n"
)

cat(
  "  Observations:",
  nrow(dat),
  "\n"
)

cat(
  "  Physical CPU threads used:",
  n_threads,
  "\n"
)

cat(
  "  Convergence flag:",
  model_converged,
  "\n"
)

cat(
  "  Fitted mean range:",
  min(fitted_mu),
  "to",
  max(fitted_mu),
  "\n"
)

cat(
  "\nEffect comparison:\n"
)

print(
  comparison,
  n = Inf,
  width = Inf
)

cat(
  "\nNegative Binomial dispersion:\n"
)

print(
  dispersion_comparison,
  width = Inf
)

cat(
  "\nFrequentist variance components:\n"
)

print(
  as_tibble(
    variance_components
  ),
  n = Inf,
  width = Inf
)

cat(
  "\nEffects with same estimated direction:",
  sum(
    comparison$same_direction
  ),
  "/",
  nrow(
    comparison
  ),
  "\n"
)

cat(
  "====================================================\n"
)

message(
  "Step 20 frequentist comparison completed."
)