# ============================================================
# 04_build_sensor_catchments.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Construct the primary 500 m circular pedestrian-accessibility
# catchment around each eligible 2025 pedestrian sensor.
#
# Input:
#   data/interim/pedestrian_sensors_2025.gpkg
#
# Output:
#   data/interim/sensor_catchments_500m.gpkg
# ============================================================

library(here)
library(sf)
library(dplyr)

# ------------------------------------------------------------
# 1. Load final eligible 2025 sensor points
# ------------------------------------------------------------

sensors_2025 <- st_read(
  here(
    "data",
    "interim",
    "pedestrian_sensors_2025.gpkg"
  ),
  quiet = TRUE
)

# ------------------------------------------------------------
# 2. Validate sensor layer
# ------------------------------------------------------------

stopifnot(
  nrow(sensors_2025) == 53
)

stopifnot(
  st_crs(sensors_2025)$epsg == 7855
)

stopifnot(
  sum(st_is_empty(sensors_2025)) == 0
)

stopifnot(
  sum(!st_is_valid(sensors_2025)) == 0
)

stopifnot(
  anyDuplicated(
    sensors_2025$location_id
  ) == 0
)

# ------------------------------------------------------------
# 3. Create 500 m circular catchments
# ------------------------------------------------------------

sensor_catchments_500m <- st_buffer(
  sensors_2025,
  dist = 500
)

# ------------------------------------------------------------
# 4. Calculate catchment areas
# ------------------------------------------------------------

catchment_area_m2 <- as.numeric(
  st_area(
    sensor_catchments_500m
  )
)

# Theoretical area of a 500 m radius circle
theoretical_area_m2 <- pi * 500^2

sensor_catchments_500m <-
  sensor_catchments_500m %>%
  mutate(
    catchment_area_m2 =
      catchment_area_m2,
    
    area_difference_pct =
      100 *
      (
        catchment_area_m2 -
          theoretical_area_m2
      ) /
      theoretical_area_m2
  )

# ------------------------------------------------------------
# 5. Validate catchments
# ------------------------------------------------------------

stopifnot(
  nrow(sensor_catchments_500m) == 53
)

stopifnot(
  st_crs(sensor_catchments_500m)$epsg == 7855
)

stopifnot(
  sum(
    st_is_empty(
      sensor_catchments_500m
    )
  ) == 0
)

stopifnot(
  sum(
    !st_is_valid(
      sensor_catchments_500m
    )
  ) == 0
)

stopifnot(
  all(
    sensor_catchments_500m$catchment_area_m2 > 0
  )
)

stopifnot(
  anyDuplicated(
    sensor_catchments_500m$location_id
  ) == 0
)

# Area should be very close to the theoretical value
stopifnot(
  max(
    abs(
      sensor_catchments_500m$area_difference_pct
    )
  ) < 1
)

# ------------------------------------------------------------
# 6. Print QA summary
# ------------------------------------------------------------

cat(
  "\n========== 500 M CATCHMENT QA ==========\n"
)

cat(
  "Number of catchments:",
  nrow(sensor_catchments_500m),
  "\n"
)

cat(
  "CRS EPSG:",
  st_crs(sensor_catchments_500m)$epsg,
  "\n"
)

cat(
  "Theoretical 500 m circle area:",
  round(
    theoretical_area_m2,
    2
  ),
  "m2\n"
)

cat(
  "\nObserved catchment areas:\n"
)

print(
  summary(
    sensor_catchments_500m$catchment_area_m2
  )
)

cat(
  "\nPercent difference from theoretical area:\n"
)

print(
  summary(
    sensor_catchments_500m$area_difference_pct
  )
)

cat(
  "\nInvalid geometries:",
  sum(
    !st_is_valid(
      sensor_catchments_500m
    )
  ),
  "\n"
)

cat(
  "Empty geometries:",
  sum(
    st_is_empty(
      sensor_catchments_500m
    )
  ),
  "\n"
)

cat(
  "\nCatchment bounding box:\n"
)

print(
  st_bbox(
    sensor_catchments_500m
  )
)

cat(
  "========================================\n"
)

# ------------------------------------------------------------
# 7. Save 500 m catchment layer
# ------------------------------------------------------------

catchment_output <- here(
  "data",
  "interim",
  "sensor_catchments_500m.gpkg"
)

if (
  file.exists(
    catchment_output
  )
) {
  file.remove(
    catchment_output
  )
}

st_write(
  sensor_catchments_500m,
  catchment_output,
  quiet = TRUE
)

message(
  "500 m pedestrian sensor catchments created successfully."
)