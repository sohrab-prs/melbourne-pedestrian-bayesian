# ============================================================
# 09_select_predictors.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Freeze the primary predictor set and prepare variables for
# Bayesian modelling.
#
# Primary predictors:
#
# Temporal:
#   day_of_week
#   public_holiday
#
# Weather:
#   mean_temp_c
#   log1p(precipitation_mm)
#
# Spatial:
#   business_count_500m
#   resident_population_500m
#   built_coverage_500m
#   median_structure_height_500m
#
# Continuous predictors are standardised at the level at
# which they actually vary:
#
#   - spatial predictors: 53 unique sensors
#   - weather predictors: 365 unique calendar days
#
# Outputs:
#   data/processed/modelling_dataset_selected_2025.csv
#   data/processed/predictor_scaling_2025.csv
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)

# ============================================================
# PART A — LOAD DATA
# ============================================================

# ------------------------------------------------------------
# 1. Load complete sensor-day modelling data
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

# ------------------------------------------------------------
# 2. Load complete 365-day weather/calendar table
# ------------------------------------------------------------

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
# 3. Structural validation
# ------------------------------------------------------------

stopifnot(
  nrow(model_data) == 18671
)

stopifnot(
  n_distinct(
    model_data$location_id
  ) == 53
)

stopifnot(
  nrow(daily_context) == 365
)

stopifnot(
  anyDuplicated(
    daily_context$date
  ) == 0
)

stopifnot(
  min(
    daily_context$date
  ) == as.Date("2025-01-01")
)

stopifnot(
  max(
    daily_context$date
  ) == as.Date("2025-12-31")
)

# ============================================================
# PART B — DEFINE PRIMARY PREDICTORS
# ============================================================

# ------------------------------------------------------------
# 4. Predictor sets
# ------------------------------------------------------------

spatial_predictors <- c(
  "business_count_500m",
  "resident_population_500m",
  "built_coverage_500m",
  "median_structure_height_500m"
)

weather_predictors <- c(
  "mean_temp_c",
  "precipitation_mm"
)

stopifnot(
  all(
    spatial_predictors %in%
      names(model_data)
  )
)

stopifnot(
  all(
    weather_predictors %in%
      names(daily_context)
  )
)

# ============================================================
# PART C — STANDARDISE SPATIAL PREDICTORS AT SENSOR LEVEL
# ============================================================

# ------------------------------------------------------------
# 5. Create one unique row per sensor
# ------------------------------------------------------------

sensor_predictors <- model_data %>%
  select(
    location_id,
    all_of(
      spatial_predictors
    )
  ) %>%
  distinct()

stopifnot(
  nrow(sensor_predictors) == 53
)

stopifnot(
  anyDuplicated(
    sensor_predictors$location_id
  ) == 0
)

# ------------------------------------------------------------
# 6. Calculate spatial scaling parameters
# ------------------------------------------------------------

spatial_scaling <- tibble(
  predictor =
    spatial_predictors,
  
  mean =
    vapply(
      spatial_predictors,
      function(v) {
        mean(
          sensor_predictors[[v]]
        )
      },
      numeric(1)
    ),
  
  sd =
    vapply(
      spatial_predictors,
      function(v) {
        sd(
          sensor_predictors[[v]]
        )
      },
      numeric(1)
    ),
  
  level =
    "sensor"
)

stopifnot(
  all(
    spatial_scaling$sd > 0
  )
)

# ------------------------------------------------------------
# 7. Standardise spatial predictors
# ------------------------------------------------------------

sensor_predictors_selected <- sensor_predictors %>%
  mutate(
    z_business_count =
      (
        business_count_500m -
          mean(
            business_count_500m
          )
      ) /
      sd(
        business_count_500m
      ),
    
    z_resident_population =
      (
        resident_population_500m -
          mean(
            resident_population_500m
          )
      ) /
      sd(
        resident_population_500m
      ),
    
    z_built_coverage =
      (
        built_coverage_500m -
          mean(
            built_coverage_500m
          )
      ) /
      sd(
        built_coverage_500m
      ),
    
    z_structure_height =
      (
        median_structure_height_500m -
          mean(
            median_structure_height_500m
          )
      ) /
      sd(
        median_structure_height_500m
      )
  ) %>%
  select(
    location_id,
    z_business_count,
    z_resident_population,
    z_built_coverage,
    z_structure_height
  )

# ============================================================
# PART D — STANDARDISE WEATHER AT DATE LEVEL
# ============================================================

# ------------------------------------------------------------
# 8. Use the complete validated 365-day weather table
# ------------------------------------------------------------

daily_predictors <- daily_context %>%
  select(
    date,
    mean_temp_c,
    precipitation_mm
  ) %>%
  arrange(
    date
  )

stopifnot(
  nrow(daily_predictors) == 365
)

stopifnot(
  anyDuplicated(
    daily_predictors$date
  ) == 0
)

stopifnot(
  !anyNA(
    daily_predictors$mean_temp_c
  )
)

stopifnot(
  !anyNA(
    daily_predictors$precipitation_mm
  )
)

# ------------------------------------------------------------
# 9. Transform precipitation
#
# log1p:
#   - retains dry days at exactly zero;
#   - reduces strong right skew;
#   - allows diminishing marginal effects of heavier rainfall.
# ------------------------------------------------------------

daily_predictors <- daily_predictors %>%
  mutate(
    log1p_precipitation =
      log1p(
        precipitation_mm
      )
  )

# ------------------------------------------------------------
# 10. Calculate weather scaling parameters
# ------------------------------------------------------------

weather_scaling <- tibble(
  predictor =
    c(
      "mean_temp_c",
      "log1p_precipitation"
    ),
  
  mean =
    c(
      mean(
        daily_predictors$mean_temp_c
      ),
      
      mean(
        daily_predictors$log1p_precipitation
      )
    ),
  
  sd =
    c(
      sd(
        daily_predictors$mean_temp_c
      ),
      
      sd(
        daily_predictors$log1p_precipitation
      )
    ),
  
  level =
    "date"
)

stopifnot(
  all(
    weather_scaling$sd > 0
  )
)

# ------------------------------------------------------------
# 11. Standardise weather predictors
# ------------------------------------------------------------

daily_predictors_selected <- daily_predictors %>%
  mutate(
    z_mean_temp =
      (
        mean_temp_c -
          mean(
            mean_temp_c
          )
      ) /
      sd(
        mean_temp_c
      ),
    
    z_log1p_precipitation =
      (
        log1p_precipitation -
          mean(
            log1p_precipitation
          )
      ) /
      sd(
        log1p_precipitation
      )
  ) %>%
  select(
    date,
    z_mean_temp,
    z_log1p_precipitation
  )

# ============================================================
# PART E — PREPARE CATEGORICAL PREDICTORS
# ============================================================

# ------------------------------------------------------------
# 12. Day-of-week coding
#
# Monday = reference category.
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

model_selected <- model_data %>%
  mutate(
    day_of_week =
      factor(
        day_of_week,
        levels = dow_levels
      ),
    
    dow_index =
      as.integer(
        day_of_week
      ),
    
    public_holiday_num =
      as.integer(
        public_holiday
      )
  )

stopifnot(
  all(
    model_selected$dow_index %in%
      1:7
  )
)

stopifnot(
  all(
    model_selected$public_holiday_num %in%
      c(
        0,
        1
      )
  )
)

# ============================================================
# PART F — JOIN STANDARDISED PREDICTORS
# ============================================================

# ------------------------------------------------------------
# 13. Join sensor-level spatial z scores
# ------------------------------------------------------------

model_selected <- model_selected %>%
  left_join(
    sensor_predictors_selected,
    by = "location_id"
  )

# ------------------------------------------------------------
# 14. Join date-level weather z scores
# ------------------------------------------------------------

model_selected <- model_selected %>%
  left_join(
    daily_predictors_selected,
    by = "date"
  )

# Join must not change observation count
stopifnot(
  nrow(model_selected) ==
    nrow(model_data)
)

# ============================================================
# PART G — SELECT FINAL MODELLING VARIABLES
# ============================================================

# ------------------------------------------------------------
# 15. Retain only variables needed for primary models
# ------------------------------------------------------------

model_selected <- model_selected %>%
  select(
    # Identifiers
    location_id,
    sensor_index,
    date,
    
    # Response
    daily_ped_count,
    
    # Temporal predictors
    day_of_week,
    dow_index,
    public_holiday_num,
    
    # Weather
    z_mean_temp,
    z_log1p_precipitation,
    
    # Spatial context
    z_business_count,
    z_resident_population,
    z_built_coverage,
    z_structure_height
  ) %>%
  arrange(
    sensor_index,
    date
  )

# ============================================================
# PART H — FINAL VALIDATION
# ============================================================

# ------------------------------------------------------------
# 16. Structural checks
# ------------------------------------------------------------

stopifnot(
  nrow(model_selected) == 18671
)

stopifnot(
  n_distinct(
    model_selected$location_id
  ) == 53
)

stopifnot(
  n_distinct(
    model_selected$sensor_index
  ) == 53
)

stopifnot(
  anyDuplicated(
    model_selected[
      c(
        "location_id",
        "date"
      )
    ]
  ) == 0
)

stopifnot(
  sum(
    is.na(
      model_selected
    )
  ) == 0
)

# ------------------------------------------------------------
# 17. Verify sensor-level z scores
# ------------------------------------------------------------

sensor_z_check <- model_selected %>%
  select(
    location_id,
    z_business_count,
    z_resident_population,
    z_built_coverage,
    z_structure_height
  ) %>%
  distinct()

stopifnot(
  nrow(sensor_z_check) == 53
)

# Means approximately zero
stopifnot(
  all(
    abs(
      colMeans(
        sensor_z_check[
          -1
        ]
      )
    ) < 1e-10
  )
)

# SD approximately one
stopifnot(
  all(
    abs(
      vapply(
        sensor_z_check[
          -1
        ],
        sd,
        numeric(1)
      ) - 1
    ) < 1e-10
  )
)

# ------------------------------------------------------------
# 18. Verify date-level z scores
# ------------------------------------------------------------

date_z_check <- daily_predictors_selected

stopifnot(
  nrow(date_z_check) == 365
)

stopifnot(
  abs(
    mean(
      date_z_check$z_mean_temp
    )
  ) < 1e-10
)

stopifnot(
  abs(
    sd(
      date_z_check$z_mean_temp
    ) - 1
  ) < 1e-10
)

stopifnot(
  abs(
    mean(
      date_z_check$z_log1p_precipitation
    )
  ) < 1e-10
)

stopifnot(
  abs(
    sd(
      date_z_check$z_log1p_precipitation
    ) - 1
  ) < 1e-10
)

# ------------------------------------------------------------
# 19. Check whether all 365 dates appear in response data
#
# This is diagnostic only.
# The modelling table does not need every calendar date if
# no eligible sensor had a complete 24-hour observation.
# ------------------------------------------------------------

response_dates <- sort(
  unique(
    model_selected$date
  )
)

missing_response_dates <- setdiff(
  daily_context$date,
  response_dates
)

# ============================================================
# PART I — PRINT QA
# ============================================================

cat(
  "\n========== PREDICTOR SELECTION QA ==========\n"
)

cat(
  "Observations:",
  nrow(model_selected),
  "\n"
)

cat(
  "Sensors:",
  n_distinct(
    model_selected$location_id
  ),
  "\n"
)

cat(
  "Calendar days available in weather data:",
  nrow(daily_context),
  "\n"
)

cat(
  "Calendar dates represented in pedestrian response:",
  n_distinct(
    model_selected$date
  ),
  "\n"
)

cat(
  "Calendar dates with no complete eligible sensor-day:",
  length(
    missing_response_dates
  ),
  "\n"
)

if (
  length(
    missing_response_dates
  ) > 0
) {
  
  cat(
    "\nDates absent from pedestrian modelling observations:\n"
  )
  
  print(
    missing_response_dates
  )
}

cat(
  "\nPrimary fixed-effect predictors:\n"
)

cat(
  "  day_of_week\n",
  "  public_holiday\n",
  "  mean_temp_c\n",
  "  log1p(precipitation_mm)\n",
  "  business_count_500m\n",
  "  resident_population_500m\n",
  "  built_coverage_500m\n",
  "  median_structure_height_500m\n",
  sep = ""
)

cat(
  "\nSensor-level standardised predictor summaries:\n"
)

print(
  summary(
    sensor_z_check[
      -1
    ]
  )
)

cat(
  "\nDate-level standardised predictor summaries:\n"
)

print(
  summary(
    date_z_check[
      c(
        "z_mean_temp",
        "z_log1p_precipitation"
      )
    ]
  )
)

cat(
  "\nDay-of-week frequencies:\n"
)

print(
  table(
    model_selected$day_of_week
  )
)

cat(
  "\nPublic-holiday frequencies:\n"
)

print(
  table(
    model_selected$public_holiday_num
  )
)

cat(
  "=============================================\n"
)

# ============================================================
# PART J — SAVE OUTPUTS
# ============================================================

# ------------------------------------------------------------
# 20. Combine scaling metadata
# ------------------------------------------------------------

predictor_scaling <- bind_rows(
  spatial_scaling,
  weather_scaling
)

write_csv(
  predictor_scaling,
  here(
    "data",
    "processed",
    "predictor_scaling_2025.csv"
  )
)

write_csv(
  model_selected,
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  )
)

message(
  "Primary predictor set frozen and modelling variables prepared successfully."
)