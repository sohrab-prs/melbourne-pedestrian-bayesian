# ============================================================
# 01_download_data.R
# Melbourne Pedestrian Bayesian Project
#
# Step 3A-3B:
# Download pedestrian, sensor-location, business-establishment,
# and residential-dwelling datasets
# ============================================================

library(here)
library(readr)
library(dplyr)

# Allow sufficient time for large City of Melbourne downloads
options(timeout = max(1200, getOption("timeout")))

# ------------------------------------------------------------
# 1. Create required directories
# ------------------------------------------------------------

dir.create(
  here("data", "raw"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here("data", "metadata"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Dataset identifiers
# ------------------------------------------------------------

api_root <- paste0(
  "https://data.melbourne.vic.gov.au/",
  "api/explore/v2.1/catalog/datasets/"
)

counts_id <- "pedestrian-counting-system-monthly-counts-per-hour"

sensors_id <- "pedestrian-counting-system-sensor-locations"

business_id <-
  "business-establishments-with-address-and-industry-classification"

dwellings_id <-
  "residential-dwellings"

buildings_id <-
  "2023-building-footprints"

weather_id <- 
  "open-meteo-historical-weather-melbourne-2025"

# ------------------------------------------------------------
# 3. Construct API request for 2025 pedestrian counts
# ------------------------------------------------------------

where_2025 <- "year(sensing_date) = 2025"

counts_url <- paste0(
  api_root,
  counts_id,
  "/exports/csv?",
  "where=",
  URLencode(where_2025, reserved = TRUE),
  "&delimiter=%2C",
  "&timezone=Australia%2FMelbourne",
  "&with_bom=false"
)

# ------------------------------------------------------------
# Historical weather: central Melbourne, calendar year 2025
# ------------------------------------------------------------

weather_url <- paste0(
  "https://archive-api.open-meteo.com/v1/archive?",
  "latitude=-37.83",
  "&longitude=144.98",
  "&start_date=2025-01-01",
  "&end_date=2025-12-31",
  "&daily=temperature_2m_mean,temperature_2m_max,precipitation_sum",
  "&timezone=Australia%2FMelbourne",
  "&models=era5",
  "&format=csv"
)

# ------------------------------------------------------------
# 4. Construct API requests for sensor and 2024 CLUE data
# ------------------------------------------------------------

# Sensor-location register
sensors_url <- paste0(
  api_root,
  sensors_id,
  "/exports/csv?",
  "delimiter=%2C",
  "&timezone=Australia%2FMelbourne",
  "&with_bom=false"
)

# CLUE contextual datasets
where_2024 <- "year(census_year) = 2024"

business_url <- paste0(
  api_root,
  business_id,
  "/exports/csv?",
  "where=",
  URLencode(where_2024, reserved = TRUE),
  "&delimiter=%2C",
  "&timezone=Australia%2FMelbourne",
  "&with_bom=false"
)

dwellings_url <- paste0(
  api_root,
  dwellings_id,
  "/exports/csv?",
  "where=",
  URLencode(where_2024, reserved = TRUE),
  "&delimiter=%2C",
  "&timezone=Australia%2FMelbourne",
  "&with_bom=false"
)

# 2023 building-footprint polygons
buildings_url <- paste0(
  api_root,
  buildings_id,
  "/exports/geojson?",
  "timezone=Australia%2FMelbourne"
)

# ------------------------------------------------------------
# 5. Define raw-data filenames
# ------------------------------------------------------------

counts_file <- here(
  "data",
  "raw",
  "pedestrian_counts_hourly_2025.csv"
)

sensors_file <- here(
  "data",
  "raw",
  "pedestrian_sensor_locations_current.csv"
)

business_file <- here(
  "data",
  "raw",
  "clue_business_establishments_2024.csv"
)

dwellings_file <- here(
  "data",
  "raw",
  "clue_residential_dwellings_2024.csv"
)

buildings_file <- here(
  "data",
  "raw",
  "building_footprints_2023.geojson"
)

weather_file <- here(
  "data",
  "raw",
  "weather_melbourne_2025.csv"
)

# ------------------------------------------------------------
# 6. Download data
# ------------------------------------------------------------

download_if_needed <- function(url, destfile, min_size = 1000) {
  
  needs_download <-
    !file.exists(destfile) ||
    is.na(file.info(destfile)$size) ||
    file.info(destfile)$size < min_size
  
  if (needs_download) {
    
    message("Downloading: ", basename(destfile))
    
    download.file(
      url = url,
      destfile = destfile,
      mode = "wb",
      method = "libcurl"
    )
    
    # Validate downloaded file size
    downloaded_size <- file.info(destfile)$size
    
    if (
      is.na(downloaded_size) ||
      downloaded_size < min_size
    ) {
      stop(
        "Downloaded file appears incomplete or invalid: ",
        basename(destfile)
      )
    }
    
  } else {
    
    message("Using existing file: ", basename(destfile))
  }
}

download_if_needed(
  counts_url,
  counts_file,
  min_size = 1e6
)

download_if_needed(
  sensors_url,
  sensors_file,
  min_size = 1000
)

download_if_needed(
  business_url,
  business_file,
  min_size = 10000
)

download_if_needed(
  dwellings_url,
  dwellings_file,
  min_size = 10000
)

download_if_needed(
  buildings_url,
  buildings_file,
  min_size = 1e6
)

download_if_needed(
  weather_url,
  weather_file,
  min_size = 5000
)

message("Download checks completed.")

# ------------------------------------------------------------
# 7. Check that files exist
# ------------------------------------------------------------

required_files <- c(
  counts_file,
  sensors_file,
  business_file,
  dwellings_file,
  buildings_file,
  weather_file
)

stopifnot(
  all(file.exists(required_files))
)

cat("\nDownloaded/raw files:\n")

print(
  file.info(required_files)[
    , c("size", "mtime")
  ]
)

# ------------------------------------------------------------
# 8. Update data provenance register
# ------------------------------------------------------------

metadata_file <- here(
  "data",
  "metadata",
  "data_sources.csv"
)

new_sources <- tibble(
  dataset_id = c(
    counts_id,
    sensors_id,
    business_id,
    dwellings_id,
    buildings_id,
    weather_id
  ),
  
  dataset_name = c(
    "Pedestrian Counting System (counts per hour)",
    "Pedestrian Counting System - Sensor Locations",
    "Business establishments location and industry classification",
    "Residential dwellings",
    "2023 Building Footprints",
    "Open-Meteo historical weather - central Melbourne"
  ),
  
  provider = c(
    rep("City of Melbourne", 5),
    "Open-Meteo"
  ),
  
  source_url = c(
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      counts_id, "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      sensors_id, "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      business_id, "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      dwellings_id, "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      buildings_id, "/"
    ),
    "https://open-meteo.com/en/docs/historical-weather-api"
  ),
  
  reference_year = c(
    "2025",
    "Current register; 2025 location derived from hourly records",
    "2024",
    "2024",
    "2023",
    "2025"
  ),
  
  download_date = rep(
    Sys.Date(),
    6
  ),
  
  spatial_unit = c(
    "Pedestrian sensor-hour",
    "Pedestrian sensor",
    "Business establishment",
    "Residential building/location",
    "Building footprint polygon",
    "Day at central Melbourne reference coordinate"
  ),
  
  geometry_type = c(
    "Point",
    "Point",
    "Point",
    "Point",
    "Polygon/MultiPolygon",
    "Point"
  ),
  
  crs_original = c(
    rep("WGS84 geographic coordinates", 4),
    "GeoJSON spatial reference",
    "WGS84 geographic coordinates"
  ),
  
  license = c(
    rep("CC BY", 5),
    "CC BY 4.0"
  ),
  
  raw_filename = c(
    basename(counts_file),
    basename(sensors_file),
    basename(business_file),
    basename(dwellings_file),
    basename(buildings_file),
    basename(weather_file)
  ),
  
  notes = c(
    "API-filtered to calendar year 2025. Raw file must not be manually edited.",
    
    paste(
      "Current sensor register used for sensor metadata.",
      "Historical 2025 spatial position is derived from hourly observations."
    ),
    
    paste(
      "API-filtered to CLUE census year 2024.",
      "Point-level business establishment data."
    ),
    
    paste(
      "API-filtered to CLUE census year 2024.",
      "Point-level residential dwelling data."
    ),
    
    paste(
      "City of Melbourne 2023 building-footprint dataset.",
      "Individual features contain date_captured values from multiple",
      "capture campaigns. Overlapping footprints will be unioned",
      "before ground-coverage calculations."
    ),
    
    paste(
      "ERA5 historical reanalysis at latitude -37.83,",
      "longitude 144.98. Daily mean/max temperature and",
      "precipitation for calendar year 2025."
    )
  )
)

# Read existing register if it already contains rows
if (
  file.exists(metadata_file) &&
  file.info(metadata_file)$size > 0
) {
  
  existing_sources <- suppressMessages(
    read_csv(
      metadata_file,
      col_types = cols(
        download_date = col_date(),
        .default = col_character()
      ),
      show_col_types = FALSE
    )
  )
  
} else {
  
  existing_sources <- tibble()
}

# Avoid duplicate entries when script is rerun
if (nrow(existing_sources) > 0) {
  
  existing_sources <- existing_sources %>%
    filter(
      !dataset_id %in% new_sources$dataset_id
    )
}

if (nrow(existing_sources) > 0) {
  
  existing_sources <- existing_sources %>%
    mutate(
      download_date = as.Date(download_date)
    )
}

updated_sources <- bind_rows(
  existing_sources,
  new_sources
)

write_csv(
  updated_sources,
  metadata_file
)

message("Data provenance register updated.")