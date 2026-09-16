# ============================================================
# 08_eda.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Perform exploratory analysis at the appropriate statistical
# levels before selecting predictors for the Bayesian models.
#
# Analysis levels:
#   - response:       sensor-day (18,671 observations)
#   - spatial:        sensor (53 observations)
#   - weather:        calendar day (365 observations)
#
# Outputs:
#   outputs/tables/
#   outputs/figures/
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ------------------------------------------------------------
# 1. Create output directories
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 2. Load datasets
# ------------------------------------------------------------

model_data <- read_csv(
  here(
    "data",
    "processed",
    "modelling_dataset_raw_2025.csv"
  ),
  col_types = cols(
    date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

spatial_data <- read_csv(
  here(
    "data",
    "processed",
    "spatial_covariates_500m.csv"
  ),
  show_col_types = FALSE
)

daily_context <- read_csv(
  here(
    "data",
    "processed",
    "daily_context_2025.csv"
  ),
  col_types = cols(
    date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 3. Structural QA
# ------------------------------------------------------------

stopifnot(
  nrow(model_data) == 18671,
  n_distinct(model_data$location_id) == 53,
  nrow(spatial_data) == 53,
  nrow(daily_context) == 365
)

# holiday_name is intentionally NA on non-public-holidays.
# All other variables should be complete.

substantive_variables <- setdiff(
  names(model_data),
  "holiday_name"
)

substantive_missing <- colSums(
  is.na(
    model_data[
      substantive_variables
    ]
  )
)

stopifnot(
  all(
    substantive_missing == 0
  )
)

cat(
  "\n========== EDA QA ==========\n"
)

cat(
  "Sensor-day observations:",
  nrow(model_data),
  "\n"
)

cat(
  "Unique sensors:",
  n_distinct(model_data$location_id),
  "\n"
)

cat(
  "Daily weather/calendar records:",
  nrow(daily_context),
  "\n"
)

cat(
  "Substantive modelling-variable missing values:",
  sum(substantive_missing),
  "\n"
)

cat(
  "holiday_name missing values (expected):",
  sum(is.na(model_data$holiday_name)),
  "\n"
)

# ============================================================
# PART A — RESPONSE DISTRIBUTION
# ============================================================

# ------------------------------------------------------------
# 4. Overall pedestrian response summary
# ------------------------------------------------------------

response_summary <- model_data %>%
  summarise(
    observations = n(),
    
    sensors =
      n_distinct(location_id),
    
    mean =
      mean(daily_ped_count),
    
    sd =
      sd(daily_ped_count),
    
    median =
      median(daily_ped_count),
    
    q25 =
      quantile(
        daily_ped_count,
        0.25
      ),
    
    q75 =
      quantile(
        daily_ped_count,
        0.75
      ),
    
    minimum =
      min(daily_ped_count),
    
    maximum =
      max(daily_ped_count),
    
    variance =
      var(daily_ped_count),
    
    variance_to_mean =
      var(daily_ped_count) /
      mean(daily_ped_count)
  )

write_csv(
  response_summary,
  here(
    "outputs",
    "tables",
    "eda_response_summary.csv"
  )
)

cat(
  "\nOverall pedestrian response:\n"
)

print(
  response_summary
)

# ------------------------------------------------------------
# 5. Raw response histogram
# ------------------------------------------------------------

p_response <- ggplot(
  model_data,
  aes(
    x = daily_ped_count
  )
) +
  geom_histogram(
    bins = 50
  ) +
  labs(
    title = "Distribution of Daily Pedestrian Counts",
    x = "Daily pedestrian count",
    y = "Sensor-days"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "eda_daily_pedestrian_distribution.png"
  ),
  p_response,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 6. Log-scale response distribution
# ------------------------------------------------------------

p_response_log <- ggplot(
  model_data,
  aes(
    x = daily_ped_count
  )
) +
  geom_histogram(
    bins = 50
  ) +
  scale_x_log10() +
  labs(
    title = "Distribution of Daily Pedestrian Counts: Log Scale",
    x = "Daily pedestrian count (log10 scale)",
    y = "Sensor-days"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "eda_daily_pedestrian_distribution_log.png"
  ),
  p_response_log,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# PART B — SENSOR HETEROGENEITY
# ============================================================

# ------------------------------------------------------------
# 7. Pedestrian activity by sensor
# ------------------------------------------------------------

sensor_response <- model_data %>%
  group_by(
    location_id
  ) %>%
  summarise(
    n_days = n(),
    
    mean_daily_count =
      mean(daily_ped_count),
    
    median_daily_count =
      median(daily_ped_count),
    
    sd_daily_count =
      sd(daily_ped_count),
    
    minimum_daily_count =
      min(daily_ped_count),
    
    maximum_daily_count =
      max(daily_ped_count),
    
    .groups = "drop"
  )

write_csv(
  sensor_response,
  here(
    "outputs",
    "tables",
    "eda_sensor_response_summary.csv"
  )
)

cat(
  "\nSensor-level mean pedestrian counts:\n"
)

print(
  summary(
    sensor_response$mean_daily_count
  )
)

cat(
  "\nSensor-level median pedestrian counts:\n"
)

print(
  summary(
    sensor_response$median_daily_count
  )
)

# Sensor means figure
p_sensor_mean <- ggplot(
  sensor_response,
  aes(
    x = reorder(
      factor(location_id),
      mean_daily_count
    ),
    y = mean_daily_count
  )
) +
  geom_point() +
  coord_flip() +
  labs(
    title = "Mean Daily Pedestrian Activity by Sensor",
    x = "Sensor ID",
    y = "Mean daily pedestrian count"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "eda_sensor_mean_counts.png"
  ),
  p_sensor_mean,
  width = 8,
  height = 10,
  dpi = 300
)

# ============================================================
# PART C — TEMPORAL PATTERNS
# ============================================================

# ------------------------------------------------------------
# 8. Day-of-week pedestrian pattern
#
# First calculate a mean within each sensor/day-of-week.
# Then average those sensor means so sensors contribute
# approximately equally despite different numbers of valid days.
# ------------------------------------------------------------

dow_levels <- c(
  "Monday",
  "Tuesday",
  "Wednesday",
  "Thursday",
  "Friday",
  "Saturday",
  "Sunday"
)

dow_sensor <- model_data %>%
  group_by(
    location_id,
    day_of_week
  ) %>%
  summarise(
    sensor_mean =
      mean(daily_ped_count),
    
    .groups = "drop"
  )

dow_summary <- dow_sensor %>%
  group_by(
    day_of_week
  ) %>%
  summarise(
    mean_of_sensor_means =
      mean(sensor_mean),
    
    median_of_sensor_means =
      median(sensor_mean),
    
    sd_of_sensor_means =
      sd(sensor_mean),
    
    sensors =
      n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    day_of_week =
      factor(
        day_of_week,
        levels = dow_levels
      )
  ) %>%
  arrange(
    day_of_week
  )

write_csv(
  dow_summary,
  here(
    "outputs",
    "tables",
    "eda_day_of_week_summary.csv"
  )
)

cat(
  "\nBalanced day-of-week response summary:\n"
)

print(
  dow_summary
)

p_dow <- ggplot(
  dow_summary,
  aes(
    x = day_of_week,
    y = mean_of_sensor_means
  )
) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "Mean Pedestrian Activity by Day of Week",
    x = "Day of week",
    y = "Mean of sensor-specific daily means"
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  here(
    "outputs",
    "figures",
    "eda_day_of_week.png"
  ),
  p_dow,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Public-holiday pattern
# ------------------------------------------------------------

holiday_sensor <- model_data %>%
  group_by(
    location_id,
    public_holiday
  ) %>%
  summarise(
    sensor_mean =
      mean(daily_ped_count),
    
    n_days = n(),
    
    .groups = "drop"
  )

holiday_summary <- holiday_sensor %>%
  group_by(
    public_holiday
  ) %>%
  summarise(
    mean_of_sensor_means =
      mean(sensor_mean),
    
    median_of_sensor_means =
      median(sensor_mean),
    
    sensors =
      n(),
    
    .groups = "drop"
  )

write_csv(
  holiday_summary,
  here(
    "outputs",
    "tables",
    "eda_public_holiday_summary.csv"
  )
)

cat(
  "\nPublic-holiday response summary:\n"
)

print(
  holiday_summary
)

# ------------------------------------------------------------
# 10. Seasonal pedestrian pattern
# ------------------------------------------------------------

season_sensor <- model_data %>%
  group_by(
    location_id,
    season
  ) %>%
  summarise(
    sensor_mean =
      mean(daily_ped_count),
    
    .groups = "drop"
  )

season_summary <- season_sensor %>%
  group_by(
    season
  ) %>%
  summarise(
    mean_of_sensor_means =
      mean(sensor_mean),
    
    median_of_sensor_means =
      median(sensor_mean),
    
    sensors =
      n(),
    
    .groups = "drop"
  ) %>%
  arrange(
    season
  )

write_csv(
  season_summary,
  here(
    "outputs",
    "tables",
    "eda_season_summary.csv"
  )
)

cat(
  "\nSeason response summary:\n"
)

print(
  season_summary
)

# ============================================================
# PART D — WEATHER
# ============================================================

# ------------------------------------------------------------
# 11. Weather summaries
# ------------------------------------------------------------

weather_summary <- daily_context %>%
  summarise(
    mean_temp_mean =
      mean(mean_temp_c),
    
    mean_temp_sd =
      sd(mean_temp_c),
    
    mean_temp_min =
      min(mean_temp_c),
    
    mean_temp_max =
      max(mean_temp_c),
    
    max_temp_mean =
      mean(max_temp_c),
    
    precipitation_mean =
      mean(precipitation_mm),
    
    precipitation_median =
      median(precipitation_mm),
    
    precipitation_max =
      max(precipitation_mm),
    
    rainy_days =
      sum(precipitation_mm > 0)
  )

write_csv(
  weather_summary,
  here(
    "outputs",
    "tables",
    "eda_weather_summary.csv"
  )
)

cat(
  "\nWeather summary:\n"
)

print(
  weather_summary
)

# ------------------------------------------------------------
# 12. Weather internal correlations
# ------------------------------------------------------------

weather_cor <- cor(
  daily_context %>%
    select(
      mean_temp_c,
      max_temp_c,
      precipitation_mm
    ),
  method = "spearman",
  use = "complete.obs"
)

cat(
  "\nWeather Spearman correlation matrix:\n"
)

print(
  round(
    weather_cor,
    3
  )
)

write_csv(
  as.data.frame(
    weather_cor
  ) %>%
    rownames_to_column(
      "variable"
    ),
  here(
    "outputs",
    "tables",
    "eda_weather_spearman_correlations.csv"
  )
)

# ============================================================
# PART E — SPATIAL PREDICTOR CORRELATIONS
# ============================================================

# ------------------------------------------------------------
# 13. Define ONLY the seven substantive spatial candidates
#
# Auxiliary variables such as structure_count_500m,
# built_area_m2 and business_category_count_500m are retained
# for QA/descriptive purposes but are not candidate model
# predictors at this stage.
# ------------------------------------------------------------

spatial_vars <- c(
  "business_count_500m",
  "business_mix_entropy_500m",
  "dwelling_count_500m",
  "resident_population_500m",
  "income_context_500m",
  "built_coverage_500m",
  "median_structure_height_500m"
)

stopifnot(
  all(
    spatial_vars %in%
      names(spatial_data)
  )
)

# ------------------------------------------------------------
# 14. Pearson spatial correlation matrix
# ------------------------------------------------------------

spatial_cor_pearson <- cor(
  spatial_data %>%
    select(
      all_of(spatial_vars)
    ),
  method = "pearson",
  use = "complete.obs"
)

write_csv(
  as.data.frame(
    spatial_cor_pearson
  ) %>%
    rownames_to_column(
      "variable"
    ),
  here(
    "outputs",
    "tables",
    "eda_spatial_pearson_correlations.csv"
  )
)

# ------------------------------------------------------------
# 15. Spearman spatial correlation matrix
#
# This is particularly useful because several intensity
# variables may be skewed.
# ------------------------------------------------------------

spatial_cor_spearman <- cor(
  spatial_data %>%
    select(
      all_of(spatial_vars)
    ),
  method = "spearman",
  use = "complete.obs"
)

write_csv(
  as.data.frame(
    spatial_cor_spearman
  ) %>%
    rownames_to_column(
      "variable"
    ),
  here(
    "outputs",
    "tables",
    "eda_spatial_spearman_correlations.csv"
  )
)

cat(
  "\nSpatial predictor Spearman correlations (53 sensors):\n"
)

print(
  round(
    spatial_cor_spearman,
    2
  )
)

# ------------------------------------------------------------
# 16. Spatial correlation heatmap
# ------------------------------------------------------------

spatial_cor_long <- as.data.frame(
  spatial_cor_spearman
) %>%
  rownames_to_column(
    "variable_1"
  ) %>%
  pivot_longer(
    cols = -variable_1,
    names_to = "variable_2",
    values_to = "rho"
  )

p_spatial_cor <- ggplot(
  spatial_cor_long,
  aes(
    x = variable_1,
    y = variable_2,
    fill = rho
  )
) +
  geom_tile() +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        rho
      )
    ),
    size = 3
  ) +
  labs(
    title = "Spearman Correlations Among Candidate Spatial Predictors",
    x = NULL,
    y = NULL,
    fill = "rho"
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )

ggsave(
  here(
    "outputs",
    "figures",
    "eda_spatial_correlation_heatmap.png"
  ),
  p_spatial_cor,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================
# PART F — SPATIAL PREDICTORS VS PEDESTRIAN ACTIVITY
# ============================================================

# ------------------------------------------------------------
# 17. Create one row per sensor
# ------------------------------------------------------------

sensor_level <- spatial_data %>%
  select(
    location_id,
    all_of(
      spatial_vars
    )
  ) %>%
  inner_join(
    sensor_response,
    by = "location_id"
  )

stopifnot(
  nrow(sensor_level) == 53
)

# ------------------------------------------------------------
# 18. Correlation with sensor-level pedestrian activity
# ------------------------------------------------------------

spatial_outcome_cor <- tibble(
  predictor =
    spatial_vars,
  
  pearson_with_mean =
    vapply(
      spatial_vars,
      function(v) {
        cor(
          sensor_level[[v]],
          sensor_level$mean_daily_count,
          method = "pearson",
          use = "complete.obs"
        )
      },
      numeric(1)
    ),
  
  spearman_with_mean =
    vapply(
      spatial_vars,
      function(v) {
        cor(
          sensor_level[[v]],
          sensor_level$mean_daily_count,
          method = "spearman",
          use = "complete.obs"
        )
      },
      numeric(1)
    ),
  
  pearson_with_median =
    vapply(
      spatial_vars,
      function(v) {
        cor(
          sensor_level[[v]],
          sensor_level$median_daily_count,
          method = "pearson",
          use = "complete.obs"
        )
      },
      numeric(1)
    ),
  
  spearman_with_median =
    vapply(
      spatial_vars,
      function(v) {
        cor(
          sensor_level[[v]],
          sensor_level$median_daily_count,
          method = "spearman",
          use = "complete.obs"
        )
      },
      numeric(1)
    )
) %>%
  arrange(
    desc(
      abs(
        spearman_with_mean
      )
    )
  )

write_csv(
  spatial_outcome_cor,
  here(
    "outputs",
    "tables",
    "eda_spatial_predictor_outcome_correlations.csv"
  )
)

cat(
  "\nSpatial predictors vs sensor-level pedestrian activity:\n"
)

print(
  spatial_outcome_cor
)

# ============================================================
# PART G — DISTRIBUTIONAL / TRANSFORMATION CHECKS
# ============================================================

# ------------------------------------------------------------
# 19. Skewness helper
# ------------------------------------------------------------

sample_skewness <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (
    length(x) < 3 ||
    sd(x) == 0
  ) {
    return(
      NA_real_
    )
  }
  
  mean(
    (
      (
        x - mean(x)
      ) /
        sd(x)
    )^3
  )
}

# ------------------------------------------------------------
# 20. Spatial predictor spread
# ------------------------------------------------------------

spatial_spread <- tibble(
  predictor =
    spatial_vars,
  
  minimum =
    vapply(
      spatial_vars,
      function(v) {
        min(
          spatial_data[[v]]
        )
      },
      numeric(1)
    ),
  
  maximum =
    vapply(
      spatial_vars,
      function(v) {
        max(
          spatial_data[[v]]
        )
      },
      numeric(1)
    ),
  
  mean =
    vapply(
      spatial_vars,
      function(v) {
        mean(
          spatial_data[[v]]
        )
      },
      numeric(1)
    ),
  
  sd =
    vapply(
      spatial_vars,
      function(v) {
        sd(
          spatial_data[[v]]
        )
      },
      numeric(1)
    ),
  
  skewness =
    vapply(
      spatial_vars,
      function(v) {
        sample_skewness(
          spatial_data[[v]]
        )
      },
      numeric(1)
    )
) %>%
  mutate(
    coefficient_of_variation =
      sd / mean
  )

write_csv(
  spatial_spread,
  here(
    "outputs",
    "tables",
    "eda_spatial_predictor_spread.csv"
  )
)

cat(
  "\nSpatial predictor spread/skewness:\n"
)

print(
  spatial_spread
)

# ------------------------------------------------------------
# 21. Compare raw vs log1p transformations for count/intensity
# predictors
# ------------------------------------------------------------

log_candidate_vars <- c(
  "business_count_500m",
  "dwelling_count_500m",
  "resident_population_500m"
)

transformation_diagnostics <- tibble(
  predictor =
    log_candidate_vars,
  
  raw_skewness =
    vapply(
      log_candidate_vars,
      function(v) {
        sample_skewness(
          spatial_data[[v]]
        )
      },
      numeric(1)
    ),
  
  log1p_skewness =
    vapply(
      log_candidate_vars,
      function(v) {
        sample_skewness(
          log1p(
            spatial_data[[v]]
          )
        )
      },
      numeric(1)
    )
)

# Add precipitation separately because it varies by day
precip_transform <- tibble(
  predictor =
    "precipitation_mm",
  
  raw_skewness =
    sample_skewness(
      daily_context$precipitation_mm
    ),
  
  log1p_skewness =
    sample_skewness(
      log1p(
        daily_context$precipitation_mm
      )
    )
)

transformation_diagnostics <-
  bind_rows(
    transformation_diagnostics,
    precip_transform
  )

write_csv(
  transformation_diagnostics,
  here(
    "outputs",
    "tables",
    "eda_transformation_diagnostics.csv"
  )
)

cat(
  "\nRaw versus log1p skewness:\n"
)

print(
  transformation_diagnostics
)

# ============================================================
# PART H — TEMPERATURE REDUNDANCY
# ============================================================

# ------------------------------------------------------------
# 22. Explicit mean-vs-maximum temperature check
# ------------------------------------------------------------

temperature_pearson <- cor(
  daily_context$mean_temp_c,
  daily_context$max_temp_c,
  method = "pearson"
)

temperature_spearman <- cor(
  daily_context$mean_temp_c,
  daily_context$max_temp_c,
  method = "spearman"
)

temperature_redundancy <- tibble(
  comparison =
    "mean_temp_c vs max_temp_c",
  
  pearson =
    temperature_pearson,
  
  spearman =
    temperature_spearman
)

write_csv(
  temperature_redundancy,
  here(
    "outputs",
    "tables",
    "eda_temperature_redundancy.csv"
  )
)

cat(
  "\nMean vs maximum temperature redundancy:\n"
)

print(
  temperature_redundancy
)

# ============================================================
# PART I — WEATHER VS CITY-WIDE PEDESTRIAN PATTERN
# ============================================================

# ------------------------------------------------------------
# 23. Calculate city-wide daily pedestrian response
#
# This is descriptive only. It is NOT substituted for the
# hierarchical sensor-day response in the Bayesian models.
# ------------------------------------------------------------

daily_response <- model_data %>%
  group_by(
    date
  ) %>%
  summarise(
    sensors_reporting =
      n_distinct(location_id),
    
    mean_sensor_count =
      mean(daily_ped_count),
    
    median_sensor_count =
      median(daily_ped_count),
    
    .groups = "drop"
  )

daily_weather_response <- daily_context %>%
  select(
    date,
    mean_temp_c,
    max_temp_c,
    precipitation_mm
  ) %>%
  inner_join(
    daily_response,
    by = "date"
  )

weather_outcome_cor <- tibble(
  predictor =
    c(
      "mean_temp_c",
      "max_temp_c",
      "precipitation_mm"
    ),
  
  spearman_with_daily_mean =
    c(
      cor(
        daily_weather_response$mean_temp_c,
        daily_weather_response$mean_sensor_count,
        method = "spearman"
      ),
      
      cor(
        daily_weather_response$max_temp_c,
        daily_weather_response$mean_sensor_count,
        method = "spearman"
      ),
      
      cor(
        daily_weather_response$precipitation_mm,
        daily_weather_response$mean_sensor_count,
        method = "spearman"
      )
    )
)

write_csv(
  weather_outcome_cor,
  here(
    "outputs",
    "tables",
    "eda_weather_outcome_correlations.csv"
  )
)

cat(
  "\nWeather vs descriptive city-wide pedestrian activity:\n"
)

print(
  weather_outcome_cor
)

# ------------------------------------------------------------
# 24. Precipitation distribution
# ------------------------------------------------------------

p_precip <- ggplot(
  daily_context,
  aes(
    x = precipitation_mm
  )
) +
  geom_histogram(
    bins = 40
  ) +
  labs(
    title = "Distribution of Daily Precipitation",
    x = "Daily precipitation (mm)",
    y = "Days"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "eda_precipitation_distribution.png"
  ),
  p_precip,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 25. Temperature distribution
# ------------------------------------------------------------

p_temperature <- ggplot(
  daily_context,
  aes(
    x = mean_temp_c
  )
) +
  geom_histogram(
    bins = 35
  ) +
  labs(
    title = "Distribution of Daily Mean Temperature",
    x = "Daily mean temperature (C)",
    y = "Days"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "eda_mean_temperature_distribution.png"
  ),
  p_temperature,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# PART J — FINAL EDA SUMMARY
# ============================================================

cat(
  "\n========== EDA COMPLETE ==========\n"
)

cat(
  "EDA tables saved to outputs/tables/\n"
)

cat(
  "EDA figures saved to outputs/figures/\n"
)

cat(
  "No predictors have been removed at this stage.\n"
)

cat(
  "Predictor selection will be performed in the next step.\n"
)

cat(
  "==================================\n"
)