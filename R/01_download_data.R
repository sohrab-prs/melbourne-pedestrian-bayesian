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

message("Download checks completed.")

# ------------------------------------------------------------
# 7. Check that files exist
# ------------------------------------------------------------

required_files <- c(
  counts_file,
  sensors_file,
  business_file,
  dwellings_file
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
    dwellings_id
  ),
  
  dataset_name = c(
    "Pedestrian Counting System (counts per hour)",
    "Pedestrian Counting System - Sensor Locations",
    "Business establishments location and industry classification",
    "Residential dwellings"
  ),
  
  provider = rep(
    "City of Melbourne",
    4
  ),
  
  source_url = c(
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      counts_id,
      "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      sensors_id,
      "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      business_id,
      "/"
    ),
    paste0(
      "https://data.melbourne.vic.gov.au/explore/dataset/",
      dwellings_id,
      "/"
    )
  ),
  
  reference_year = c(
    "2025",
    "Current register; 2025 location derived from hourly records",
    "2024",
    "2024"
  ),
  
  download_date = rep(
    Sys.Date(),
    4
  ),
  
  spatial_unit = c(
    "Pedestrian sensor-hour",
    "Pedestrian sensor",
    "Business establishment",
    "Residential building/location"
  ),
  
  geometry_type = rep(
    "Point",
    4
  ),
  
  crs_original = rep(
    "WGS84 geographic coordinates",
    4
  ),
  
  license = rep(
    "CC BY",
    4
  ),
  
  raw_filename = c(
    basename(counts_file),
    basename(sensors_file),
    basename(business_file),
    basename(dwellings_file)
  ),
  
  notes = c(
    paste(
      "API-filtered to calendar year 2025.",
      "Raw file must not be manually edited."
    ),
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