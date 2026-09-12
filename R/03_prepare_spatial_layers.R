# ============================================================
# 03_prepare_spatial_layers.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Construct the final spatial point layer for eligible
# pedestrian sensors using their observed 2025 locations.
#
# Input:
#   data/raw/pedestrian_counts_hourly_2025.csv
#   data/raw/pedestrian_sensor_locations_current.csv
#   data/processed/sensor_eligibility_2025.csv
#
# Output:
#   data/interim/pedestrian_sensors_2025.gpkg
# ============================================================

library(here)
library(readr)
library(dplyr)
library(stringr)
library(sf)
library(tibble)

# ------------------------------------------------------------
# 1. Create interim directory if required
# ------------------------------------------------------------

dir.create(
  here("data", "interim"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load required data
# ------------------------------------------------------------

counts_raw <- read_csv(
  here(
    "data",
    "raw",
    "pedestrian_counts_hourly_2025.csv"
  ),
  show_col_types = FALSE
)

sensor_eligibility <- read_csv(
  here(
    "data",
    "processed",
    "sensor_eligibility_2025.csv"
  ),
  show_col_types = FALSE
)

sensors_current <- read_csv(
  here(
    "data",
    "raw",
    "pedestrian_sensor_locations_current.csv"
  ),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 3. Identify final eligible sensor IDs
# ------------------------------------------------------------

final_ids <- sensor_eligibility %>%
  filter(
    eligible_final
  ) %>%
  pull(location_id)

cat(
  "\nFinal eligible sensor IDs:",
  length(final_ids),
  "\n"
)

# Our Step 4 result should contain 53 final sensors
stopifnot(
  length(final_ids) == 53
)

stopifnot(
  anyDuplicated(final_ids) == 0
)

# ------------------------------------------------------------
# 4. Extract one unique 2025 location per eligible sensor
# ------------------------------------------------------------

sensor_location_2025 <- counts_raw %>%
  filter(
    location_id %in% final_ids,
    !is.na(location)
  ) %>%
  distinct(
    location_id,
    location
  )

# Number of spatial locations associated with each sensor ID
location_counts <- sensor_location_2025 %>%
  count(
    location_id,
    name = "n_locations"
  )

# All final sensors must have exactly one recorded location
stopifnot(
  all(location_counts$n_locations == 1)
)

# All 53 final sensors must have a usable location
stopifnot(
  n_distinct(sensor_location_2025$location_id) ==
    length(final_ids)
)

# ------------------------------------------------------------
# 5. Attach one representative sensor name
# ------------------------------------------------------------

sensor_names_2025 <- counts_raw %>%
  filter(
    location_id %in% final_ids,
    !is.na(sensor_name)
  ) %>%
  group_by(location_id) %>%
  summarise(
    sensor_name = first(sensor_name),
    .groups = "drop"
  )

sensor_location_2025 <- sensor_location_2025 %>%
  left_join(
    sensor_names_2025,
    by = "location_id"
  )

# ------------------------------------------------------------
# 6. Parse longitude and latitude from location field
# ------------------------------------------------------------

extract_melbourne_coordinates <- function(x) {
  
  nums <- as.numeric(
    str_extract_all(
      x,
      "-?[0-9]+(?:\\.[0-9]+)?"
    )[[1]]
  )
  
  latitude_candidates <- nums[
    nums >= -40 &
      nums <= -35
  ]
  
  longitude_candidates <- nums[
    nums >= 140 &
      nums <= 150
  ]
  
  latitude <- if (
    length(latitude_candidates) >= 1
  ) {
    latitude_candidates[[1]]
  } else {
    NA_real_
  }
  
  longitude <- if (
    length(longitude_candidates) >= 1
  ) {
    longitude_candidates[[1]]
  } else {
    NA_real_
  }
  
  c(
    longitude = longitude,
    latitude = latitude
  )
}

parsed_coordinates <- t(
  vapply(
    sensor_location_2025$location,
    extract_melbourne_coordinates,
    numeric(2)
  )
)

sensor_location_2025 <- sensor_location_2025 %>%
  mutate(
    longitude =
      parsed_coordinates[, "longitude"],
    
    latitude =
      parsed_coordinates[, "latitude"]
  )

# ------------------------------------------------------------
# 7. Validate parsed coordinates
# ------------------------------------------------------------

stopifnot(
  !anyNA(
    sensor_location_2025$longitude
  )
)

stopifnot(
  !anyNA(
    sensor_location_2025$latitude
  )
)

stopifnot(
  all(
    sensor_location_2025$longitude >= 140 &
      sensor_location_2025$longitude <= 150
  )
)

stopifnot(
  all(
    sensor_location_2025$latitude >= -40 &
      sensor_location_2025$latitude <= -35
  )
)

stopifnot(
  nrow(sensor_location_2025) == 53
)

# ------------------------------------------------------------
# 8. Create WGS84 spatial sensor layer
# ------------------------------------------------------------

sensors_2025_wgs84 <- st_as_sf(
  sensor_location_2025,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
)

# ------------------------------------------------------------
# 9. Transform to GDA2020 / MGA Zone 55
#
# EPSG:7855 uses metres and is appropriate for Melbourne.
# All future buffering and area calculations will use this CRS.
# ------------------------------------------------------------

sensors_2025 <- st_transform(
  sensors_2025_wgs84,
  crs = 7855
)

stopifnot(
  st_crs(sensors_2025)$epsg == 7855
)

stopifnot(
  nrow(sensors_2025) == 53
)

stopifnot(
  sum(st_is_empty(sensors_2025)) == 0
)

stopifnot(
  sum(!st_is_valid(sensors_2025)) == 0
)

# ------------------------------------------------------------
# 10. Cross-check against current sensor register
#
# Historical 2025 coordinates remain authoritative.
# This comparison is only a QA check.
# ------------------------------------------------------------

current_check <- sensors_current %>%
  filter(
    location_id %in% final_ids,
    !is.na(longitude),
    !is.na(latitude)
  ) %>%
  select(
    location_id,
    current_longitude = longitude,
    current_latitude = latitude
  ) %>%
  distinct(
    location_id,
    .keep_all = TRUE
  )

current_check_sf <- st_as_sf(
  current_check,
  coords = c(
    "current_longitude",
    "current_latitude"
  ),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(
    crs = 7855
  )

# Restrict both datasets to common sensor IDs
common_ids <- intersect(
  sensors_2025$location_id,
  current_check_sf$location_id
)

historical_check <- sensors_2025 %>%
  filter(
    location_id %in% common_ids
  ) %>%
  arrange(location_id)

current_check_sf <- current_check_sf %>%
  filter(
    location_id %in% common_ids
  ) %>%
  arrange(location_id)

stopifnot(
  identical(
    as.character(historical_check$location_id),
    as.character(current_check_sf$location_id)
  )
)

coordinate_difference_m <- as.numeric(
  st_distance(
    historical_check,
    current_check_sf,
    by_element = TRUE
  )
)

coordinate_check <- tibble(
  location_id =
    historical_check$location_id,
  
  coordinate_difference_m =
    coordinate_difference_m
)

# ------------------------------------------------------------
# 11. Print QA summary
# ------------------------------------------------------------

cat(
  "\n========== SENSOR SPATIAL QA ==========\n"
)

cat(
  "Final sensors:",
  nrow(sensors_2025),
  "\n"
)

cat(
  "CRS EPSG:",
  st_crs(sensors_2025)$epsg,
  "\n"
)

cat(
  "Sensors cross-checked against current register:",
  nrow(coordinate_check),
  "\n"
)

cat(
  "\nCurrent-vs-2025 coordinate differences (metres):\n"
)

print(
  summary(
    coordinate_check$coordinate_difference_m
  )
)

cat(
  "\nTen largest coordinate differences:\n"
)

print(
  coordinate_check %>%
    arrange(
      desc(coordinate_difference_m)
    ) %>%
    head(10)
)

cat(
  "\n2025 sensor bounding box:\n"
)

print(
  st_bbox(sensors_2025)
)

cat(
  "=======================================\n"
)

# ------------------------------------------------------------
# 12. Save final 2025 sensor spatial layer
# ------------------------------------------------------------

sensor_output <- here(
  "data",
  "interim",
  "pedestrian_sensors_2025.gpkg"
)

if (file.exists(sensor_output)) {
  file.remove(sensor_output)
}

st_write(
  sensors_2025,
  sensor_output,
  quiet = TRUE
)

message(
  "Final 2025 pedestrian sensor spatial layer created successfully."
)