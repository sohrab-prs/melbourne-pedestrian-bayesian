# ============================================================
# 02_clean_pedestrian.R
# Melbourne Pedestrian Bayesian Project
#
# Construct complete daily pedestrian observations and
# determine the final eligible sensor sample.
# ============================================================

library(here)
library(readr)
library(dplyr)

# ------------------------------------------------------------
# 1. Load raw hourly pedestrian data
# ------------------------------------------------------------

counts_raw <- read_csv(
  here(
    "data",
    "raw",
    "pedestrian_counts_hourly_2025.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    sensing_date = as.Date(sensing_date)
  )

# ------------------------------------------------------------
# 2. Basic integrity checks
# ------------------------------------------------------------

stopifnot(
  min(counts_raw$sensing_date, na.rm = TRUE) ==
    as.Date("2025-01-01"),
  
  max(counts_raw$sensing_date, na.rm = TRUE) ==
    as.Date("2025-12-31")
)

expected_hours <- sort(
  unique(
    counts_raw$hourday[
      !is.na(counts_raw$hourday)
    ]
  )
)

cat(
  "Hourly values found:",
  paste(expected_hours, collapse = ", "),
  "\n"
)

stopifnot(
  length(expected_hours) == 24
)

duplicate_sensor_hours <- counts_raw %>%
  count(
    location_id,
    sensing_date,
    hourday,
    name = "n"
  ) %>%
  filter(n > 1)

stopifnot(
  nrow(duplicate_sensor_hours) == 0
)

stopifnot(
  sum(
    counts_raw$pedestriancount < 0,
    na.rm = TRUE
  ) == 0
)

# ------------------------------------------------------------
# 3. Check spatial-location stability
# ------------------------------------------------------------

location_check <- counts_raw %>%
  filter(
    !is.na(location)
  ) %>%
  distinct(
    location_id,
    location
  ) %>%
  count(
    location_id,
    name = "n_locations"
  )

# ------------------------------------------------------------
# 4. Compute hourly coverage by sensor
# ------------------------------------------------------------

sensor_hourly_coverage <- counts_raw %>%
  group_by(location_id) %>%
  summarise(
    n_hourly_records = n(),
    
    hourly_completeness_pct =
      100 * n_hourly_records / (365 * 24),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5. Preliminary sensor eligibility
# ------------------------------------------------------------

preliminary_sensors <- sensor_hourly_coverage %>%
  inner_join(
    location_check,
    by = "location_id"
  ) %>%
  filter(
    n_locations == 1,
    hourly_completeness_pct >= 90
  )

cat(
  "Sensors passing location + 90% hourly coverage:",
  nrow(preliminary_sensors),
  "\n"
)

# ------------------------------------------------------------
# 6. Construct sensor-day summaries
# ------------------------------------------------------------

sensor_days_all <- counts_raw %>%
  semi_join(
    preliminary_sensors,
    by = "location_id"
  ) %>%
  group_by(
    location_id,
    sensing_date
  ) %>%
  summarise(
    n_records = n(),
    
    n_distinct_hours =
      n_distinct(hourday),
    
    n_nonmissing_counts =
      sum(!is.na(pedestriancount)),
    
    total_count_observed =
      sum(
        pedestriancount,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    complete_day =
      n_records == 24 &
      n_distinct_hours == 24 &
      n_nonmissing_counts == 24,
    
    daily_ped_count =
      if_else(
        complete_day,
        as.numeric(total_count_observed),
        NA_real_
      )
  )

# ------------------------------------------------------------
# 7. Complete-day coverage by sensor
# ------------------------------------------------------------

sensor_day_coverage <- sensor_days_all %>%
  group_by(location_id) %>%
  summarise(
    days_with_any_records = n(),
    
    complete_days =
      sum(complete_day),
    
    incomplete_days =
      sum(!complete_day),
    
    complete_day_pct =
      100 * complete_days / 365,
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 8. Final sensor eligibility
# ------------------------------------------------------------

sensor_eligibility <- preliminary_sensors %>%
  left_join(
    sensor_day_coverage,
    by = "location_id"
  ) %>%
  mutate(
    eligible_final =
      complete_days >= 300
  )

final_sensor_ids <- sensor_eligibility %>%
  filter(
    eligible_final
  ) %>%
  pull(location_id)

cat(
  "Final eligible sensors:",
  length(final_sensor_ids),
  "\n"
)

# ------------------------------------------------------------
# 9. Final daily pedestrian modelling data
# ------------------------------------------------------------

pedestrian_daily_2025 <- sensor_days_all %>%
  filter(
    location_id %in% final_sensor_ids,
    complete_day
  ) %>%
  transmute(
    location_id,
    date = sensing_date,
    daily_ped_count
  ) %>%
  arrange(
    location_id,
    date
  )

# ------------------------------------------------------------
# 10. Final validation
# ------------------------------------------------------------

stopifnot(
  !anyNA(
    pedestrian_daily_2025$daily_ped_count
  )
)

stopifnot(
  all(
    pedestrian_daily_2025$daily_ped_count >= 0
  )
)

stopifnot(
  anyDuplicated(
    pedestrian_daily_2025[
      c("location_id", "date")
    ]
  ) == 0
)

# ------------------------------------------------------------
# 11. Save derived outputs
# ------------------------------------------------------------

dir.create(
  here(
    "data",
    "processed"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  pedestrian_daily_2025,
  here(
    "data",
    "processed",
    "pedestrian_daily_2025.csv"
  )
)

write_csv(
  sensor_eligibility,
  here(
    "data",
    "processed",
    "sensor_eligibility_2025.csv"
  )
)

message(
  "Step 4 pedestrian cleaning completed."
)