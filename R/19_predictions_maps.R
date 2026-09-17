# ============================================================
# 19_predictions_maps.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Final inference and prediction from accepted M4:
#
#   1. Fixed-effect posterior summaries
#   2. Rate ratios / percentage changes
#   3. Directional posterior probabilities
#   4. Representative scenario predictions
#   5. Existing-sensor baseline predictions
#   6. Prediction/effect figures and map
#
# No model is refitted.
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(coda)
library(sf)
library(ggplot2)

set.seed(20260917)

# ============================================================
# PART A — INPUTS
# ============================================================

sample_file <- here(
  "outputs",
  "mcmc",
  "m4_sensor_day_negbin_samples.rds"
)

data_file <- here(
  "data",
  "processed",
  "modelling_dataset_selected_2025.csv"
)

sensor_file <- here(
  "data",
  "interim",
  "pedestrian_sensors_2025.gpkg"
)

stopifnot(
  file.exists(sample_file),
  file.exists(data_file),
  file.exists(sensor_file)
)

dir.create(
  here("outputs", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here("outputs", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PART B — LOAD POSTERIOR
# ============================================================

m4_samples <- readRDS(
  sample_file
)

stopifnot(
  is.mcmc.list(m4_samples),
  length(m4_samples) == 4
)

posterior_matrix <- as.matrix(
  m4_samples
)

# Use up to 10,000 posterior draws for prediction.
# This is sufficient for stable prediction summaries while
# keeping post-processing light.

n_prediction_draws <- min(
  10000L,
  nrow(posterior_matrix)
)

prediction_draw_ids <- sample(
  seq_len(
    nrow(posterior_matrix)
  ),
  n_prediction_draws,
  replace = FALSE
)

pred_matrix <- posterior_matrix[
  prediction_draw_ids,
  ,
  drop = FALSE
]

# ============================================================
# PART C — FINAL FIXED-EFFECT SUMMARIES
# ============================================================

effect_parameters <- c(
  "beta_dow[2]",
  "beta_dow[3]",
  "beta_dow[4]",
  "beta_dow[5]",
  "beta_dow[6]",
  "beta_dow[7]",
  "beta_holiday",
  "beta_temp",
  "beta_rain",
  "beta_business",
  "beta_population",
  "beta_coverage",
  "beta_height"
)

effect_labels <- c(
  
  "beta_dow[2]" =
    "Tuesday vs Monday",
  
  "beta_dow[3]" =
    "Wednesday vs Monday",
  
  "beta_dow[4]" =
    "Thursday vs Monday",
  
  "beta_dow[5]" =
    "Friday vs Monday",
  
  "beta_dow[6]" =
    "Saturday vs Monday",
  
  "beta_dow[7]" =
    "Sunday vs Monday",
  
  "beta_holiday" =
    "Public holiday",
  
  "beta_temp" =
    "Mean temperature (+1 SD)",
  
  "beta_rain" =
    "log(1 + precipitation) (+1 SD)",
  
  "beta_business" =
    "Businesses within 500 m (+1 SD)",
  
  "beta_population" =
    "Resident population within 500 m (+1 SD)",
  
  "beta_coverage" =
    "Built coverage within 500 m (+1 SD)",
  
  "beta_height" =
    "Median structure height within 500 m (+1 SD)"
)

effect_groups <- c(
  
  "beta_dow[2]" = "Day of week",
  "beta_dow[3]" = "Day of week",
  "beta_dow[4]" = "Day of week",
  "beta_dow[5]" = "Day of week",
  "beta_dow[6]" = "Day of week",
  "beta_dow[7]" = "Day of week",
  
  "beta_holiday" = "Calendar",
  "beta_temp" = "Weather",
  "beta_rain" = "Weather",
  
  "beta_business" = "Urban environment",
  "beta_population" = "Urban environment",
  "beta_coverage" = "Urban environment",
  "beta_height" = "Urban environment"
)

effect_estimates <- lapply(
  effect_parameters,
  function(parameter_name) {
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    beta_q025 <-
      unname(
        quantile(
          x,
          0.025
        )
      )
    
    beta_median <-
      median(x)
    
    beta_q975 <-
      unname(
        quantile(
          x,
          0.975
        )
      )
    
    p_positive <-
      mean(
        x > 0
      )
    
    p_negative <-
      mean(
        x < 0
      )
    
    tibble(
      
      parameter =
        parameter_name,
      
      effect =
        unname(
          effect_labels[
            parameter_name
          ]
        ),
      
      group =
        unname(
          effect_groups[
            parameter_name
          ]
        ),
      
      posterior_mean_beta =
        mean(x),
      
      posterior_sd_beta =
        sd(x),
      
      beta_q025 =
        beta_q025,
      
      beta_median =
        beta_median,
      
      beta_q975 =
        beta_q975,
      
      rate_ratio_median =
        exp(
          beta_median
        ),
      
      rate_ratio_q025 =
        exp(
          beta_q025
        ),
      
      rate_ratio_q975 =
        exp(
          beta_q975
        ),
      
      percent_change_median =
        100 *
        (
          exp(
            beta_median
          ) -
            1
        ),
      
      percent_change_q025 =
        100 *
        (
          exp(
            beta_q025
          ) -
            1
        ),
      
      percent_change_q975 =
        100 *
        (
          exp(
            beta_q975
          ) -
            1
        ),
      
      probability_positive =
        p_positive,
      
      probability_negative =
        p_negative,
      
      posterior_direction =
        ifelse(
          p_positive >=
            p_negative,
          "Positive",
          "Negative"
        ),
      
      directional_probability =
        max(
          p_positive,
          p_negative
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  effect_estimates,
  here(
    "outputs",
    "tables",
    "m4_final_effect_estimates.csv"
  )
)

# ============================================================
# PART D — VARIANCE / DISPERSION PARAMETERS
# ============================================================

variance_parameters <- c(
  "phi",
  "sigma_sensor",
  "sigma_day"
)

variance_labels <- c(
  
  "phi" =
    "Negative Binomial dispersion",
  
  "sigma_sensor" =
    "Sensor random-effect SD",
  
  "sigma_day" =
    "Day random-effect SD"
)

variance_summary <- lapply(
  variance_parameters,
  function(parameter_name) {
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    tibble(
      
      parameter =
        parameter_name,
      
      label =
        unname(
          variance_labels[
            parameter_name
          ]
        ),
      
      posterior_mean =
        mean(x),
      
      posterior_sd =
        sd(x),
      
      q025 =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      
      median =
        median(x),
      
      q975 =
        unname(
          quantile(
            x,
            0.975
          )
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  variance_summary,
  here(
    "outputs",
    "tables",
    "m4_variance_parameters.csv"
  )
)

# ============================================================
# PART E — REPRESENTATIVE PREDICTION SCENARIOS
# ============================================================

scenario_table <- tibble(
  
  scenario = c(
    "Reference Monday",
    "Tuesday",
    "Wednesday",
    "Thursday",
    "Friday",
    "Saturday",
    "Sunday",
    "Public holiday Monday",
    "Temperature +1 SD",
    "Precipitation +1 SD",
    "Business count +1 SD",
    "Resident population +1 SD",
    "Built coverage +1 SD",
    "Structure height +1 SD"
  ),
  
  dow = c(
    1, 2, 3, 4, 5, 6, 7,
    1, 1, 1, 1, 1, 1, 1
  ),
  
  holiday = c(
    0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0
  ),
  
  temp = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0
  ),
  
  rain = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0
  ),
  
  business = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0
  ),
  
  population = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0
  ),
  
  coverage = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0
  ),
  
  height = c(
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1
  )
)

# ------------------------------------------------------------
# Build day-of-week coefficient matrix
# ------------------------------------------------------------

dow_effect_matrix <- matrix(
  0,
  nrow = n_prediction_draws,
  ncol = 7
)

for (
  k in 2:7
) {
  
  dow_effect_matrix[
    ,
    k
  ] <-
    pred_matrix[
      ,
      paste0(
        "beta_dow[",
        k,
        "]"
      )
    ]
}

# ------------------------------------------------------------
# Generate scenario inference + prediction
# ------------------------------------------------------------

scenario_results <- lapply(
  seq_len(
    nrow(
      scenario_table
    )
  ),
  function(i) {
    
    s <-
      scenario_table[
        i,
      ]
    
    eta_fixed <-
      pred_matrix[
        ,
        "alpha"
      ] +
      dow_effect_matrix[
        ,
        s$dow
      ] +
      pred_matrix[
        ,
        "beta_holiday"
      ] *
      s$holiday +
      pred_matrix[
        ,
        "beta_temp"
      ] *
      s$temp +
      pred_matrix[
        ,
        "beta_rain"
      ] *
      s$rain +
      pred_matrix[
        ,
        "beta_business"
      ] *
      s$business +
      pred_matrix[
        ,
        "beta_population"
      ] *
      s$population +
      pred_matrix[
        ,
        "beta_coverage"
      ] *
      s$coverage +
      pred_matrix[
        ,
        "beta_height"
      ] *
      s$height
    
    # Conditional expected count:
    # random sensor and day effects fixed at zero.
    expected_count <-
      exp(
        eta_fixed
      )
    
    # Posterior prediction for a new sensor-day:
    # integrate over sensor heterogeneity, day heterogeneity,
    # and observation-level Negative Binomial variation.
    new_sensor_effect <-
      rnorm(
        n_prediction_draws,
        mean = 0,
        sd =
          pred_matrix[
            ,
            "sigma_sensor"
          ]
      )
    
    new_day_effect <-
      rnorm(
        n_prediction_draws,
        mean = 0,
        sd =
          pred_matrix[
            ,
            "sigma_day"
          ]
      )
    
    predictive_mu <-
      exp(
        eta_fixed +
          new_sensor_effect +
          new_day_effect
      )
    
    predictive_count <-
      rnbinom(
        n_prediction_draws,
        size =
          pred_matrix[
            ,
            "phi"
          ],
        mu =
          predictive_mu
      )
    
    tibble(
      
      scenario =
        s$scenario,
      
      expected_count_q025 =
        unname(
          quantile(
            expected_count,
            0.025
          )
        ),
      
      expected_count_median =
        median(
          expected_count
        ),
      
      expected_count_q975 =
        unname(
          quantile(
            expected_count,
            0.975
          )
        ),
      
      predictive_count_q025 =
        unname(
          quantile(
            predictive_count,
            0.025
          )
        ),
      
      predictive_count_median =
        median(
          predictive_count
        ),
      
      predictive_count_q975 =
        unname(
          quantile(
            predictive_count,
            0.975
          )
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  scenario_results,
  here(
    "outputs",
    "tables",
    "m4_scenario_predictions.csv"
  )
)

# ============================================================
# PART F — EXISTING-SENSOR BASELINE PREDICTIONS
# ============================================================

dat <- read_csv(
  data_file,
  col_types = cols(
    date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

sensor_level <- dat %>%
  select(
    sensor_index,
    location_id,
    z_business_count,
    z_resident_population,
    z_built_coverage,
    z_structure_height
  ) %>%
  distinct() %>%
  mutate(
    sensor_index =
      as.integer(
        sensor_index
      )
  ) %>%
  arrange(
    sensor_index
  )

stopifnot(
  nrow(
    sensor_level
  ) == 53,
  all(
    sensor_level$sensor_index ==
      seq_len(53)
  )
)

sensor_prediction_results <- lapply(
  seq_len(
    nrow(
      sensor_level
    )
  ),
  function(j) {
    
    s <-
      sensor_level[
        j,
      ]
    
    u_name <-
      paste0(
        "u_sensor[",
        j,
        "]"
      )
    
    eta_sensor <-
      pred_matrix[
        ,
        "alpha"
      ] +
      pred_matrix[
        ,
        "beta_business"
      ] *
      s$z_business_count +
      pred_matrix[
        ,
        "beta_population"
      ] *
      s$z_resident_population +
      pred_matrix[
        ,
        "beta_coverage"
      ] *
      s$z_built_coverage +
      pred_matrix[
        ,
        "beta_height"
      ] *
      s$z_structure_height +
      pred_matrix[
        ,
        u_name
      ]
    
    mu_sensor <-
      exp(
        eta_sensor
      )
    
    tibble(
      
      sensor_index =
        j,
      
      location_id =
        as.character(
          s$location_id
        ),
      
      reference_expected_q025 =
        unname(
          quantile(
            mu_sensor,
            0.025
          )
        ),
      
      reference_expected_median =
        median(
          mu_sensor
        ),
      
      reference_expected_q975 =
        unname(
          quantile(
            mu_sensor,
            0.975
          )
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  sensor_prediction_results,
  here(
    "outputs",
    "tables",
    "m4_sensor_reference_predictions.csv"
  )
)

# ============================================================
# PART G — SENSOR PREDICTION MAP
# ============================================================

sensor_sf <- st_read(
  sensor_file,
  quiet = TRUE
)

possible_id_names <- c(
  "location_id",
  "Location_ID",
  "LOCATION_ID",
  "sensor_id",
  "Sensor_ID"
)

spatial_id_name <-
  possible_id_names[
    possible_id_names %in%
      names(
        sensor_sf
      )
  ][1]

if (
  is.na(
    spatial_id_name
  )
) {
  
  stop(
    paste(
      "Could not identify sensor ID field. Available columns:",
      paste(
        names(
          sensor_sf
        ),
        collapse = ", "
      )
    )
  )
}

sensor_sf <- sensor_sf %>%
  mutate(
    location_id_join =
      as.character(
        .data[[spatial_id_name]]
      )
  )

sensor_prediction_map <- sensor_sf %>%
  select(
    location_id_join
  ) %>%
  inner_join(
    sensor_prediction_results %>%
      mutate(
        location_id_join =
          as.character(
            location_id
          )
      ),
    by =
      "location_id_join"
  )

stopifnot(
  inherits(
    sensor_prediction_map,
    "sf"
  ),
  nrow(
    sensor_prediction_map
  ) == 53
)

sensor_prediction_map <-
  st_transform(
    sensor_prediction_map,
    7855
  )

p_map <- ggplot(
  sensor_prediction_map
) +
  geom_sf(
    aes(
      colour =
        reference_expected_median,
      size =
        reference_expected_median
    ),
    alpha =
      0.85
  ) +
  labs(
    title =
      "M4 Predicted Baseline Pedestrian Activity by Sensor",
    subtitle =
      paste(
        "Posterior median expected daily count:",
        "Monday, non-holiday, average weather;"
      ),
    colour =
      "Expected count",
    size =
      "Expected count"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_sensor_prediction_map.png"
  ),
  p_map,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# PART H — EFFECT FIGURE
# ============================================================

effect_plot_data <- effect_estimates %>%
  mutate(
    effect =
      factor(
        effect,
        levels =
          rev(
            effect
          )
      )
  )

p_effect <- ggplot(
  effect_plot_data,
  aes(
    x =
      percent_change_median,
    y =
      effect
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
      xmin = percent_change_q025,
      xmax = percent_change_q975
    ),
    orientation = "y",
    width = 0.18
  ) +
  geom_point(
    size =
      2
  ) +
  facet_grid(
    group ~ .,
    scales =
      "free_y",
    space =
      "free_y"
  ) +
  labs(
    title =
      "M4 Posterior Effects on Expected Pedestrian Activity",
    subtitle =
      "Posterior median percentage change and 95% credible interval",
    x =
      "Percentage change in expected daily pedestrian count",
    y =
      NULL
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_final_effect_estimates.png"
  ),
  p_effect,
  width = 10,
  height = 10,
  dpi = 300
)

# ============================================================
# PART I — SCENARIO FIGURE
# ============================================================

scenario_plot_data <- scenario_results %>%
  mutate(
    scenario =
      factor(
        scenario,
        levels =
          rev(
            scenario
          )
      )
  )

p_scenarios <- ggplot(
  scenario_plot_data,
  aes(
    x =
      expected_count_median,
    y =
      scenario
  )
) +
  geom_errorbar(
    aes(
      xmin = expected_count_q025,
      xmax = expected_count_q975
    ),
    orientation = "y",
    width = 0.18
  ) +
  geom_point(
    size =
      2
  ) +
  labs(
    title =
      "M4 Population-Level Scenario Predictions",
    subtitle =
      "Expected count conditional on sensor and day random effects equal to zero",
    x =
      "Expected daily pedestrian count",
    y =
      NULL
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_scenario_predictions.png"
  ),
  p_scenarios,
  width = 9,
  height = 8,
  dpi = 300
)

# ============================================================
# PART J — QA
# ============================================================

cat(
  "\n========== STEP 19 FINAL M4 INFERENCE ==========\n"
)

cat(
  "\nFixed-effect estimates:\n"
)

effect_estimates %>%
  select(
    effect,
    percent_change_median,
    percent_change_q025,
    percent_change_q975,
    probability_positive,
    probability_negative,
    posterior_direction,
    directional_probability
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat(
  "\nVariance / dispersion parameters:\n"
)

print(
  variance_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nRepresentative scenario predictions:\n"
)

print(
  scenario_results,
  n = Inf,
  width = Inf
)

cat(
  "\nExisting-sensor reference predictions:\n"
)

cat(
  "  Sensors:",
  nrow(
    sensor_prediction_results
  ),
  "\n"
)

cat(
  "  Minimum posterior median expected count:",
  min(
    sensor_prediction_results$reference_expected_median
  ),
  "\n"
)

cat(
  "  Median across sensors:",
  median(
    sensor_prediction_results$reference_expected_median
  ),
  "\n"
)

cat(
  "  Maximum posterior median expected count:",
  max(
    sensor_prediction_results$reference_expected_median
  ),
  "\n"
)

cat(
  "\nPosterior draws used for predictions:",
  n_prediction_draws,
  "\n"
)

cat(
  "================================================\n"
)

message(
  "Step 19 final M4 inference and prediction completed."
)