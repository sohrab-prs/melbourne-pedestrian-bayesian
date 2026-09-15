# ============================================================
# 07_build_modelling_dataset.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Assemble the complete sensor-day modelling dataset by joining:
#
#   - daily pedestrian counts (2025)
#   - calendar variables (2025)
#   - weather variables (2025)
#   - sensor-level 500 m spatial covariates
#
# Inputs:
#   data/processed/pedestrian_daily_2025.csv
#   data/processed/daily_context_2025.csv
#   data/processed/spatial_covariates_500m.csv
#
# Outputs:
#   data/processed/modelling_dataset_raw_2025.csv
#   data/processed/sensor_index_2025.csv
# ============================================================

library(here)
library(readr)
library(dplyr)

# ------------------------------------------------------------
# 1. Load processed component datasets
# ------------------------------------------------------------

pedestrian_daily <- read_csv(
  here(
    "data",
    "processed",
    "pedestrian_daily_2025.csv"
  ),
  col_types = cols(
    location_id = col_double(),
    date = col_date(),
    daily_ped_count = col_double()
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

spatial_covariates <- read_csv(
  here(
    "data",
    "processed",
    "spatial_covariates_500m.csv"
  ),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Validate component datasets
# ------------------------------------------------------------

# Pedestrian outcome data
stopifnot(
  nrow(pedestrian_daily) == 18671
)

stopifnot(
  n_distinct(
    pedestrian_daily$location_id
  ) == 53
)

stopifnot(
  anyDuplicated(
    pedestrian_daily[
      c(
        "location_id",
        "date"
      )
    ]
  ) == 0
)

stopifnot(
  !anyNA(
    pedestrian_daily$daily_ped_count
  )
)

stopifnot(
  all(
    pedestrian_daily$daily_ped_count >= 0
  )
)

# Daily temporal/weather data
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

# Sensor-level spatial covariates
stopifnot(
  nrow(spatial_covariates) == 53
)

stopifnot(
  anyDuplicated(
    spatial_covariates$location_id
  ) == 0
)

# ------------------------------------------------------------
# 3. Verify sensor IDs match exactly
# ------------------------------------------------------------

pedestrian_sensor_ids <- sort(
  unique(
    pedestrian_daily$location_id
  )
)

spatial_sensor_ids <- sort(
  unique(
    spatial_covariates$location_id
  )
)

stopifnot(
  identical(
    pedestrian_sensor_ids,
    spatial_sensor_ids
  )
)

# ------------------------------------------------------------
# 4. Join daily calendar/weather context
# ------------------------------------------------------------

model_data <- pedestrian_daily %>%
  left_join(
    daily_context,
    by = "date"
  )

# Join must not change the observation count
stopifnot(
  nrow(model_data) ==
    nrow(pedestrian_daily)
)

# ------------------------------------------------------------
# 5. Join sensor-level spatial covariates
# ------------------------------------------------------------

model_data <- model_data %>%
  left_join(
    spatial_covariates,
    by = "location_id"
  )

stopifnot(
  nrow(model_data) ==
    nrow(pedestrian_daily)
)

# ------------------------------------------------------------
# 6. Create reproducible sensor index for hierarchical models
# ------------------------------------------------------------

sensor_index <- tibble(
  location_id =
    sort(
      unique(
        model_data$location_id
      )
    )
) %>%
  mutate(
    sensor_index =
      row_number()
  )

stopifnot(
  nrow(sensor_index) == 53
)

stopifnot(
  all(
    sensor_index$sensor_index ==
      seq_len(53)
  )
)

model_data <- model_data %>%
  left_join(
    sensor_index,
    by = "location_id"
  )

# ------------------------------------------------------------
# 7. Add useful deterministic time indices
#
# These are retained for EDA and possible diagnostics.
# They are not automatically included as final model predictors.
# ------------------------------------------------------------

model_data <- model_data %>%
  mutate(
    day_of_year =
      as.integer(
        format(
          date,
          "%j"
        )
      )
  )

# ------------------------------------------------------------
# 8. Arrange columns clearly
# ------------------------------------------------------------

model_data <- model_data %>%
  select(
    # Observation identifiers
    location_id,
    sensor_index,
    date,
    day_of_year,
    
    # Response
    daily_ped_count,
    
    # Calendar
    day_of_week,
    weekend,
    month,
    season,
    holiday_name,
    public_holiday,
    
    # Weather
    mean_temp_c,
    max_temp_c,
    precipitation_mm,
    
    # Economic / activity context
    business_count_500m,
    business_mix_entropy_500m,
    business_category_count_500m,
    
    # Residential context
    dwelling_count_500m,
    resident_population_500m,
    income_context_500m,
    
    # Morphology
    built_area_m2,
    built_coverage_500m,
    median_structure_height_500m,
    structure_count_500m
  ) %>%
  arrange(
    sensor_index,
    date
  )

# ------------------------------------------------------------
# 9. Final structural validation
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
  n_distinct(
    model_data$sensor_index
  ) == 53
)

stopifnot(
  anyDuplicated(
    model_data[
      c(
        "location_id",
        "date"
      )
    ]
  ) == 0
)

# ------------------------------------------------------------
# 10. Check missingness
# ------------------------------------------------------------

missing_summary <- model_data %>%
  summarise(
    across(
      everything(),
      ~ sum(
        is.na(.x)
      )
    )
  )

cat(
  "\n========== MODELLING DATA QA ==========\n"
)

cat(
  "Observations:",
  nrow(model_data),
  "\n"
)

cat(
  "Sensors:",
  n_distinct(
    model_data$location_id
  ),
  "\n"
)

cat(
  "Date range:",
  as.character(
    min(model_data$date)
  ),
  "to",
  as.character(
    max(model_data$date)
  ),
  "\n"
)

cat(
  "\nMissing values by variable:\n"
)

print(
  missing_summary
)

# ------------------------------------------------------------
# 11. Response summary
# ------------------------------------------------------------

response_summary <- model_data %>%
  summarise(
    n = n(),
    
    mean =
      mean(
        daily_ped_count
      ),
    
    sd =
      sd(
        daily_ped_count
      ),
    
    median =
      median(
        daily_ped_count
      ),
    
    minimum =
      min(
        daily_ped_count
      ),
    
    maximum =
      max(
        daily_ped_count
      ),
    
    variance =
      var(
        daily_ped_count
      ),
    
    variance_to_mean =
      var(
        daily_ped_count
      ) /
      mean(
        daily_ped_count
      )
  )

cat(
  "\nResponse summary:\n"
)

print(
  response_summary
)

# ------------------------------------------------------------
# 12. Observation counts per sensor
# ------------------------------------------------------------

sensor_observation_summary <- model_data %>%
  count(
    location_id,
    sensor_index,
    name = "n_days"
  ) %>%
  arrange(
    n_days
  )

cat(
  "\nComplete observations per sensor:\n"
)

print(
  summary(
    sensor_observation_summary$n_days
  )
)

# ------------------------------------------------------------
# 13. Check repeated spatial values within sensors
#
# Spatial covariates should be constant across dates for each
# sensor because they are sensor-level characteristics.
# ------------------------------------------------------------

spatial_consistency <- model_data %>%
  group_by(
    location_id
  ) %>%
  summarise(
    n_business =
      n_distinct(
        business_count_500m
      ),
    
    n_business_mix =
      n_distinct(
        business_mix_entropy_500m
      ),
    
    n_dwellings =
      n_distinct(
        dwelling_count_500m
      ),
    
    n_population =
      n_distinct(
        resident_population_500m
      ),
    
    n_income =
      n_distinct(
        income_context_500m
      ),
    
    n_coverage =
      n_distinct(
        built_coverage_500m
      ),
    
    n_height =
      n_distinct(
        median_structure_height_500m
      ),
    
    .groups = "drop"
  )

stopifnot(
  all(
    spatial_consistency[
      , -1
    ] == 1
  )
)

# ------------------------------------------------------------
# 14. Check temporal variables vary only by date
# ------------------------------------------------------------

date_context_check <- model_data %>%
  group_by(
    date
  ) %>%
  summarise(
    n_mean_temp =
      n_distinct(
        mean_temp_c
      ),
    
    n_max_temp =
      n_distinct(
        max_temp_c
      ),
    
    n_precip =
      n_distinct(
        precipitation_mm
      ),
    
    n_public_holiday =
      n_distinct(
        public_holiday
      ),
    
    .groups = "drop"
  )

stopifnot(
  all(
    date_context_check[
      , -1
    ] == 1
  )
)

# ------------------------------------------------------------
# 15. Save complete raw modelling dataset
# ------------------------------------------------------------

write_csv(
  model_data,
  here(
    "data",
    "processed",
    "modelling_dataset_raw_2025.csv"
  )
)

write_csv(
  sensor_index,
  here(
    "data",
    "processed",
    "sensor_index_2025.csv"
  )
)

cat(
  "\nModelling dataset saved successfully.\n"
)

cat(
  "=======================================\n"
)