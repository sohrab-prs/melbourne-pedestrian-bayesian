# ============================================================
# 05_engineer_spatial_covariates.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Engineer candidate urban/geospatial explanatory variables
# within 500 m of each final pedestrian sensor.
#
# Candidate spatial predictors:
#   1. resident_population_500m
#   2. income_context_500m
#   3. dwelling_count_500m
#   4. business_count_500m
#   5. business_mix_entropy_500m
#   6. built_coverage_500m
#   7. median_structure_height_500m
#
# Inputs:
#   data/interim/sensor_catchments_500m.gpkg
#   data/interim/pedestrian_sensors_2025.gpkg
#   data/raw/clue_business_establishments_2024.csv
#   data/raw/clue_residential_dwellings_2024.csv
#   data/raw/building_footprints_2023.geojson
#   ABS 2021 G01 / G02 / SA1 boundary files
#
# Outputs:
#   data/processed/spatial_covariates_500m.csv
#   data/interim/sensor_spatial_covariates_500m.gpkg
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(sf)

# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------

safe_max <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  max(
    x,
    na.rm = TRUE
  )
}

safe_min <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  min(
    x,
    na.rm = TRUE
  )
}

# Normalised Shannon entropy.
# 0 = no diversity / one category
# 1 = categories evenly represented.
normalised_entropy <- function(n) {
  
  n <- n[
    !is.na(n) &
      n > 0
  ]
  
  if (length(n) <= 1) {
    return(0)
  }
  
  p <- n / sum(n)
  
  entropy <-
    -sum(
      p * log(p)
    )
  
  entropy /
    log(length(n))
}

# ------------------------------------------------------------
# 2. Load final sensors and 500 m catchments
# ------------------------------------------------------------

sensors_2025 <- st_read(
  here(
    "data",
    "interim",
    "pedestrian_sensors_2025.gpkg"
  ),
  quiet = TRUE
)

catchments <- st_read(
  here(
    "data",
    "interim",
    "sensor_catchments_500m.gpkg"
  ),
  quiet = TRUE
)

stopifnot(
  nrow(sensors_2025) == 53,
  nrow(catchments) == 53,
  st_crs(sensors_2025)$epsg == 7855,
  st_crs(catchments)$epsg == 7855,
  anyDuplicated(catchments$location_id) == 0
)

catchment_base <- catchments %>%
  st_drop_geometry() %>%
  select(
    location_id,
    catchment_area_m2
  )

# ============================================================
# PART A — BUSINESS ACTIVITY
# ============================================================

# ------------------------------------------------------------
# 3. Load 2024 business establishments
# ------------------------------------------------------------

business_raw <- read_csv(
  here(
    "data",
    "raw",
    "clue_business_establishments_2024.csv"
  ),
  show_col_types = FALSE
)

stopifnot(
  all(business_raw$census_year == 2024),
  sum(is.na(business_raw$longitude)) == 0,
  sum(is.na(business_raw$latitude)) == 0
)

business_sf <- st_as_sf(
  business_raw,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(
    7855
  )

# ------------------------------------------------------------
# 4. Assign establishments to overlapping 500 m catchments
# ------------------------------------------------------------

business_in_catchments <- st_join(
  business_sf %>%
    select(
      industry_anzsic4_description
    ),
  catchments %>%
    select(
      location_id
    ),
  join = st_intersects,
  left = FALSE
) %>%
  st_drop_geometry()

# Total establishments
business_counts <- business_in_catchments %>%
  count(
    location_id,
    name = "business_count_500m"
  )

# ------------------------------------------------------------
# 5. Economic activity diversity
# ------------------------------------------------------------

business_category_counts <-
  business_in_catchments %>%
  filter(
    !is.na(
      industry_anzsic4_description
    ),
    industry_anzsic4_description != ""
  ) %>%
  count(
    location_id,
    industry_anzsic4_description,
    name = "n_businesses"
  )

business_diversity <-
  business_category_counts %>%
  group_by(location_id) %>%
  summarise(
    business_mix_entropy_500m =
      normalised_entropy(
        n_businesses
      ),
    
    business_category_count_500m =
      n(),
    
    .groups = "drop"
  )

business_covariates <-
  catchment_base %>%
  select(location_id) %>%
  left_join(
    business_counts,
    by = "location_id"
  ) %>%
  left_join(
    business_diversity,
    by = "location_id"
  ) %>%
  mutate(
    business_count_500m =
      replace_na(
        business_count_500m,
        0L
      ),
    
    business_category_count_500m =
      replace_na(
        business_category_count_500m,
        0L
      ),
    
    business_mix_entropy_500m =
      replace_na(
        business_mix_entropy_500m,
        0
      )
  )

# ============================================================
# PART B — RESIDENTIAL DWELLINGS
# ============================================================

# ------------------------------------------------------------
# 6. Load 2024 residential dwellings
# ------------------------------------------------------------

dwellings_raw <- read_csv(
  here(
    "data",
    "raw",
    "clue_residential_dwellings_2024.csv"
  ),
  show_col_types = FALSE
)

stopifnot(
  all(dwellings_raw$census_year == 2024),
  sum(is.na(dwellings_raw$longitude)) == 0,
  sum(is.na(dwellings_raw$latitude)) == 0
)

# Ensure dwelling count is numeric
dwellings_raw <- dwellings_raw %>%
  mutate(
    dwelling_number =
      as.numeric(
        dwelling_number
      )
  )

dwellings_sf <- st_as_sf(
  dwellings_raw,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(
    7855
  )

# ------------------------------------------------------------
# 7. Aggregate dwelling numbers within catchments
# ------------------------------------------------------------

dwellings_in_catchments <- st_join(
  dwellings_sf %>%
    select(
      dwelling_number
    ),
  catchments %>%
    select(
      location_id
    ),
  join = st_intersects,
  left = FALSE
) %>%
  st_drop_geometry()

dwelling_covariates <-
  dwellings_in_catchments %>%
  group_by(location_id) %>%
  summarise(
    dwelling_count_500m =
      sum(
        dwelling_number,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

dwelling_covariates <-
  catchment_base %>%
  select(location_id) %>%
  left_join(
    dwelling_covariates,
    by = "location_id"
  ) %>%
  mutate(
    dwelling_count_500m =
      replace_na(
        dwelling_count_500m,
        0
      )
  )

# ============================================================
# PART C — ABS CENSUS POPULATION AND INCOME
# ============================================================

# ------------------------------------------------------------
# 8. Locate ABS files
# ------------------------------------------------------------

gcp_root <- here(
  "data",
  "raw",
  "abs_gcp_sa1_vic_2021"
)

boundary_root <- here(
  "data",
  "raw",
  "abs_sa1_boundaries_2021"
)

g01_file <- list.files(
  gcp_root,
  pattern = "G01.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

g02_file <- list.files(
  gcp_root,
  pattern = "G02.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

sa1_file <- list.files(
  boundary_root,
  pattern = "\\.shp$",
  recursive = TRUE,
  full.names = TRUE
)

stopifnot(
  length(g01_file) == 1,
  length(g02_file) == 1,
  length(sa1_file) == 1
)

# ------------------------------------------------------------
# 9. Load required Census variables only
# ------------------------------------------------------------

g01 <- read_csv(
  g01_file[[1]],
  col_types = cols(
    SA1_CODE_2021 =
      col_character(),
    
    .default =
      col_guess()
  ),
  show_col_types = FALSE
) %>%
  select(
    sa1_code =
      SA1_CODE_2021,
    
    resident_population =
      Tot_P_P
  )

g02 <- read_csv(
  g02_file[[1]],
  col_types = cols(
    SA1_CODE_2021 =
      col_character(),
    
    .default =
      col_guess()
  ),
  show_col_types = FALSE
) %>%
  select(
    sa1_code =
      SA1_CODE_2021,
    
    median_hh_income_weekly =
      Median_tot_hhd_inc_weekly
  )

census_context <- g01 %>%
  left_join(
    g02,
    by = "sa1_code"
  )

stopifnot(
  anyDuplicated(
    census_context$sa1_code
  ) == 0
)

# ------------------------------------------------------------
# 10. Load Victorian SA1 geography
# ------------------------------------------------------------

sa1_raw <- st_read(
  sa1_file[[1]],
  quiet = TRUE
)

sa1_vic <- sa1_raw %>%
  filter(
    STE_NAME21 == "Victoria",
    !st_is_empty(geometry)
  ) %>%
  mutate(
    sa1_code =
      as.character(
        SA1_CODE21
      )
  ) %>%
  select(
    sa1_code,
    geometry
  ) %>%
  left_join(
    census_context,
    by = "sa1_code"
  ) %>%
  st_transform(
    7855
  )

stopifnot(
  sum(
    is.na(
      sa1_vic$resident_population
    )
  ) == 0
)

stopifnot(
  sum(
    is.na(
      sa1_vic$median_hh_income_weekly
    )
  ) == 0
)

# ------------------------------------------------------------
# 11. Calculate full SA1 areas
# ------------------------------------------------------------

sa1_vic <- sa1_vic %>%
  mutate(
    sa1_area_m2 =
      as.numeric(
        st_area(
          sa1_vic
        )
      )
  )

# Keep only SA1s intersecting at least one catchment
sa1_intersects_any <-
  lengths(
    st_intersects(
      sa1_vic,
      catchments
    )
  ) > 0

sa1_relevant <- sa1_vic[
  sa1_intersects_any,
]

# ------------------------------------------------------------
# 12. Intersect Census areas with sensor catchments
# ------------------------------------------------------------

census_intersections <-
  suppressWarnings(
    st_intersection(
      sa1_relevant,
      catchments %>%
        select(
          location_id
        )
    )
  )

census_intersections <-
  census_intersections %>%
  mutate(
    intersection_area_m2 =
      as.numeric(
        st_area(
          census_intersections
        )
      ),
    
    sa1_area_share =
      pmin(
        intersection_area_m2 /
          sa1_area_m2,
        1
      ),
    
    estimated_population =
      resident_population *
      sa1_area_share
  )

# ------------------------------------------------------------
# 13. Aggregate Census context to sensor catchments
# ------------------------------------------------------------

census_covariates <-
  census_intersections %>%
  st_drop_geometry() %>%
  group_by(location_id) %>%
  summarise(
    resident_population_500m =
      sum(
        estimated_population,
        na.rm = TRUE
      ),
    
    income_context_500m =
      if (
        sum(
          estimated_population,
          na.rm = TRUE
        ) > 0
      ) {
        
        weighted.mean(
          median_hh_income_weekly,
          estimated_population,
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    .groups = "drop"
  )

# ============================================================
# PART D — BUILDING MORPHOLOGY
# ============================================================

# ------------------------------------------------------------
# 14. Load building footprints
# ------------------------------------------------------------

buildings <- st_read(
  here(
    "data",
    "raw",
    "building_footprints_2023.geojson"
  ),
  quiet = TRUE
) %>%
  st_transform(
    7855
  )

stopifnot(
  sum(
    !st_is_valid(buildings)
  ) == 0
)

# Keep only buildings intersecting at least one catchment
building_intersects_any <-
  lengths(
    st_intersects(
      buildings,
      catchments
    )
  ) > 0

buildings_relevant <-
  buildings[
    building_intersects_any,
  ]

# ------------------------------------------------------------
# 15. Calculate unioned built coverage within each catchment
#
# st_union removes overlapping footprint polygons so stacked
# podium/tower polygons do not inflate ground coverage.
# ------------------------------------------------------------

building_intersections <-
  suppressWarnings(
    st_intersection(
      buildings_relevant %>%
        select(
          objectid,
          structure_id
        ),
      catchments %>%
        select(
          location_id
        )
    )
  )

built_union <- building_intersections %>%
  group_by(location_id) %>%
  summarise(
    do_union = TRUE,
    .groups = "drop"
  )

built_coverage <-
  tibble(
    location_id =
      built_union$location_id,
    
    built_area_m2 =
      as.numeric(
        st_area(
          built_union
        )
      )
  ) %>%
  left_join(
    catchment_base,
    by = "location_id"
  ) %>%
  mutate(
    built_coverage_500m =
      built_area_m2 /
      catchment_area_m2
  )

# ------------------------------------------------------------
# 16. Derive one height value per unique structure
# ------------------------------------------------------------

structures <- buildings_relevant %>%
  filter(
    !is.na(structure_id)
  ) %>%
  group_by(
    structure_id
  ) %>%
  summarise(
    structure_max =
      safe_max(
        structure_max_elevation
      ),
    
    structure_min =
      safe_min(
        structure_min_elevation
      ),
    
    do_union = TRUE,
    
    .groups = "drop"
  ) %>%
  mutate(
    structure_height_m =
      structure_max -
      structure_min
  ) %>%
  filter(
    !is.na(structure_height_m),
    structure_height_m >= 0
  )

# Use one representative point per structure to avoid
# double-counting vertically stacked footprint polygons.
structure_points <-
  suppressWarnings(
    st_point_on_surface(
      structures
    )
  )

# ------------------------------------------------------------
# 17. Assign structures to catchments
# ------------------------------------------------------------

structures_in_catchments <- st_join(
  structure_points %>%
    select(
      structure_id,
      structure_height_m
    ),
  catchments %>%
    select(
      location_id
    ),
  join = st_intersects,
  left = FALSE
) %>%
  st_drop_geometry()

height_covariates <-
  structures_in_catchments %>%
  group_by(
    location_id
  ) %>%
  summarise(
    median_structure_height_500m =
      median(
        structure_height_m,
        na.rm = TRUE
      ),
    
    structure_count_500m =
      n_distinct(
        structure_id
      ),
    
    .groups = "drop"
  )

# ============================================================
# PART E — COMBINE SPATIAL COVARIATES
# ============================================================

# ------------------------------------------------------------
# 18. Build final sensor-level spatial covariate table
# ------------------------------------------------------------

spatial_covariates_500m <-
  catchment_base %>%
  select(
    location_id
  ) %>%
  
  left_join(
    business_covariates,
    by = "location_id"
  ) %>%
  
  left_join(
    dwelling_covariates,
    by = "location_id"
  ) %>%
  
  left_join(
    census_covariates,
    by = "location_id"
  ) %>%
  
  left_join(
    built_coverage %>%
      select(
        location_id,
        built_area_m2,
        built_coverage_500m
      ),
    by = "location_id"
  ) %>%
  
  left_join(
    height_covariates,
    by = "location_id"
  ) %>%
  
  arrange(
    location_id
  )

# ------------------------------------------------------------
# 19. Validate final covariate table
# ------------------------------------------------------------

stopifnot(
  nrow(
    spatial_covariates_500m
  ) == 53
)

stopifnot(
  anyDuplicated(
    spatial_covariates_500m$location_id
  ) == 0
)

stopifnot(
  all(
    spatial_covariates_500m$
      business_count_500m >= 0
  )
)

stopifnot(
  all(
    spatial_covariates_500m$
      dwelling_count_500m >= 0
  )
)

stopifnot(
  all(
    spatial_covariates_500m$
      resident_population_500m >= 0
  )
)

stopifnot(
  all(
    spatial_covariates_500m$
      built_coverage_500m >= 0 &
      spatial_covariates_500m$
      built_coverage_500m <= 1
  )
)

stopifnot(
  all(
    spatial_covariates_500m$
      business_mix_entropy_500m >= 0 &
      spatial_covariates_500m$
      business_mix_entropy_500m <= 1
  )
)

# ------------------------------------------------------------
# 20. Print QA summary
# ------------------------------------------------------------

cat(
  "\n========== SPATIAL COVARIATE QA ==========\n"
)

cat(
  "Sensors:",
  nrow(
    spatial_covariates_500m
  ),
  "\n"
)

cat(
  "\nMissing values by variable:\n"
)

print(
  spatial_covariates_500m %>%
    summarise(
      across(
        everything(),
        ~ sum(is.na(.x))
      )
    )
)

cat(
  "\nSpatial covariate summaries:\n"
)

print(
  summary(
    spatial_covariates_500m %>%
      select(
        business_count_500m,
        business_mix_entropy_500m,
        dwelling_count_500m,
        resident_population_500m,
        income_context_500m,
        built_coverage_500m,
        median_structure_height_500m
      )
  )
)

cat(
  "\nBusiness establishments assigned to catchments:",
  nrow(business_in_catchments),
  "\n"
)

cat(
  "Dwelling records assigned to catchments:",
  nrow(dwellings_in_catchments),
  "\n"
)

cat(
  "Relevant SA1 polygons:",
  nrow(sa1_relevant),
  "\n"
)

cat(
  "Relevant building-footprint records:",
  nrow(buildings_relevant),
  "\n"
)

cat(
  "Unique relevant structures:",
  nrow(structures),
  "\n"
)

cat(
  "==========================================\n"
)

# ------------------------------------------------------------
# 21. Save modelling table
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
  spatial_covariates_500m,
  here(
    "data",
    "processed",
    "spatial_covariates_500m.csv"
  )
)

# ------------------------------------------------------------
# 22. Save a spatial version for maps and GIS outputs
# ------------------------------------------------------------

sensor_covariates_sf <-
  sensors_2025 %>%
  left_join(
    spatial_covariates_500m,
    by = "location_id"
  )

spatial_output <- here(
  "data",
  "interim",
  "sensor_spatial_covariates_500m.gpkg"
)

if (
  file.exists(
    spatial_output
  )
) {
  file.remove(
    spatial_output
  )
}

st_write(
  sensor_covariates_sf,
  spatial_output,
  quiet = TRUE
)

message(
  "Step 7 spatial covariate engineering completed successfully."
)