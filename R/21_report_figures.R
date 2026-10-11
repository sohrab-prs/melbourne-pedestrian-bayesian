# ============================================================
# 21_report_figures.R
# Report-only visualisation layer for final_report.qmd
#
# This script redraws final analytical outputs in a consistent report style,
# reproduces PPC visual summaries from saved posterior samples, and
# builds cartographic maps with an official City of Melbourne basemap
# and municipal boundary.
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(coda)

if (!requireNamespace("png", quietly = TRUE)) {
  stop(
    "Package 'png' is required for the report basemap. ",
    "Install it once with install.packages('png')."
  )
}

set.seed(20260918)

dir.create(
  here("outputs", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here("data", "metadata", "report_cartography"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# COMMON STYLE
# ============================================================

report_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      colour = "#1F2933"
    ),
    plot.subtitle = element_text(
      size = 11.5,
      colour = "#4B5563"
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 10.5,
      colour = "#30343B"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11.5
    ),
    strip.background = element_rect(
      fill = "#EEF3F7",
      colour = "#CBD5DD",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.35,
      colour = "#E5E9EC"
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    ),
    plot.caption = element_text(
      size = 9,
      colour = "#687078",
      hjust = 0
    )
  )

cols_model <- c(
  "M3" = "#0072B2",
  "M4" = "#D55E00"
)

cols4 <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7"
)

cols7 <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#D55E00",
  "#CC79A7",
  "#56B4E9",
  "#5B5B5B"
)

save_plot <- function(
    plot,
    filename,
    width,
    height
) {

  ggsave(
    here(
      "outputs",
      "figures",
      filename
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

# ============================================================
# DATA
# ============================================================

mod_dat <- read_csv(
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  ),
  show_col_types = FALSE
)

mod_dat <- mod_dat %>%
  mutate(
    location_id =
      as.character(
        location_id
      ),
    sensor_index =
      as.integer(
        sensor_index
      ),
    date =
      as.Date(
        date
      )
  )

sensor_lookup <- mod_dat %>%
  distinct(
    location_id,
    sensor_index
  ) %>%
  arrange(
    sensor_index
  )

date_lookup <- tibble(
  date =
    sort(
      unique(
        mod_dat$date
      )
    ),
  day_index =
    seq_along(
      sort(
        unique(
          mod_dat$date
        )
      )
    )
)

mod_dat <- mod_dat %>%
  left_join(
    date_lookup,
    by = "date"
  )

N <- nrow(mod_dat)
J <- n_distinct(
  mod_dat$sensor_index
)
D <- n_distinct(
  mod_dat$day_index
)

stopifnot(
  N == 18671,
  J == 53,
  D == 364
)

# ============================================================
# CARTOGRAPHY HELPERS
# ============================================================

sensor_file <- here(
  "data",
  "interim",
  "pedestrian_sensors_2025.gpkg"
)

buffer_file <- here(
  "data",
  "interim",
  "sensor_catchments_500m.gpkg"
)

stopifnot(
  file.exists(
    sensor_file
  ),
  file.exists(
    buffer_file
  )
)

sensors_sf <- st_read(
  sensor_file,
  quiet = TRUE
)

buffers_sf <- st_read(
  buffer_file,
  quiet = TRUE
)

rename_location_id <- function(x) {

  if (
    "location_id" %in%
      names(x)
  ) {
    return(x)
  }

  candidates <- c(
    "Location_ID",
    "LOCATION_ID",
    "locationid",
    "sensor_id",
    "Sensor_ID"
  )

  id_col <- candidates[
    candidates %in%
      names(x)
  ][1]

  if (
    is.na(id_col)
  ) {
    stop(
      "Could not identify the sensor location ID column."
    )
  }

  names(x)[
    names(x) ==
      id_col
  ] <- "location_id"

  x
}

sensors_sf <- rename_location_id(
  sensors_sf
)

buffers_sf <- rename_location_id(
  buffers_sf
)

sensors_sf <- sensors_sf %>%
  mutate(
    location_id =
      as.character(
        location_id
      )
  ) %>%
  filter(
    location_id %in%
      sensor_lookup$location_id
  ) %>%
  left_join(
    sensor_lookup,
    by =
      "location_id"
  ) %>%
  arrange(
    sensor_index
  )

buffers_sf <- buffers_sf %>%
  mutate(
    location_id =
      as.character(
        location_id
      )
  ) %>%
  filter(
    location_id %in%
      sensor_lookup$location_id
  )

stopifnot(
  nrow(
    sensors_sf
  ) == J
)

# ------------------------------------------------------------
# Official City of Melbourne municipal boundary
# ------------------------------------------------------------

carto_dir <- here(
  "data",
  "metadata",
  "report_cartography"
)

boundary_file <- file.path(
  carto_dir,
  "city_of_melbourne_municipal_boundary.geojson"
)

boundary_urls <- c(
  paste0(
    "https://data.melbourne.vic.gov.au/api/explore/v2.1/",
    "catalog/datasets/municipal-boundary/exports/geojson",
    "?lang=en&timezone=Australia%2FMelbourne"
  ),
  paste0(
    "https://gisags.melbourne.vic.gov.au/server_wa/rest/services/",
    "BasemapsWGS84/Basemap_Light_WMS/MapServer/76/query",
    "?where=1%3D1&outFields=*&returnGeometry=true&outSR=4326&f=geojson"
  )
)

if (
  !file.exists(
    boundary_file
  )
) {

  downloaded <- FALSE

  for (
    u in boundary_urls
  ) {

    downloaded <- tryCatch(
      {
        suppressWarnings(
          download.file(
            u,
            boundary_file,
            mode = "wb",
            quiet = TRUE,
            method = "libcurl"
          )
        )

        b_test <- st_read(
          boundary_file,
          quiet = TRUE
        )

        nrow(
          b_test
        ) > 0
      },
      error =
        function(e) {
          FALSE
        }
    )

    if (
      downloaded
    ) {
      break
    }
  }

  if (
    !downloaded
  ) {
    stop(
      "Could not download the City of Melbourne municipal boundary. ",
      "Check the internet connection and rerun R/21_report_figures.R."
    )
  }
}

boundary_sf <- st_read(
  boundary_file,
  quiet = TRUE
) %>%
  st_make_valid() %>%
  st_transform(
    4326
  )

# Use a single municipal polygon for plotting.
boundary_sf <- st_sf(
  geometry =
    st_union(
      boundary_sf
    )
)

city_bbox <- st_bbox(
  boundary_sf
)

x_pad <-
  as.numeric(
    city_bbox["xmax"] -
      city_bbox["xmin"]
  ) *
  0.025

y_pad <-
  as.numeric(
    city_bbox["ymax"] -
      city_bbox["ymin"]
  ) *
  0.025

map_bbox <- c(
  xmin =
    as.numeric(
      city_bbox["xmin"]
    ) -
    x_pad,
  ymin =
    as.numeric(
      city_bbox["ymin"]
    ) -
    y_pad,
  xmax =
    as.numeric(
      city_bbox["xmax"]
    ) +
    x_pad,
  ymax =
    as.numeric(
      city_bbox["ymax"]
    ) +
    y_pad
)

# ------------------------------------------------------------
# Official City of Melbourne light basemap.
# Cached after the first successful download.
# ------------------------------------------------------------

basemap_file <- file.path(
  carto_dir,
  "city_of_melbourne_light_basemap.png"
)

if (
  !file.exists(
    basemap_file
  )
) {

  basemap_url <- paste0(
    "https://gisags.melbourne.vic.gov.au/server_wa/rest/services/",
    "BasemapsWGS84/Basemap_Light_WMS/MapServer/export",
    "?bbox=",
    paste(
      format(
        map_bbox[
          c(
            "xmin",
            "ymin",
            "xmax",
            "ymax"
          )
        ],
        scientific = FALSE,
        trim = TRUE,
        digits = 12
      ),
      collapse = ","
    ),
    "&bboxSR=4326",
    "&imageSR=4326",
    "&size=1800,1400",
    "&format=png32",
    "&transparent=false",
    "&f=image"
  )

  try(
    suppressWarnings(
      download.file(
        basemap_url,
        basemap_file,
        mode = "wb",
        quiet = TRUE,
        method = "libcurl"
      )
    ),
    silent = TRUE
  )
}

basemap_image <- NULL

if (
  file.exists(
    basemap_file
  )
) {

  basemap_image <- tryCatch(
    png::readPNG(
      basemap_file
    ),
    error =
      function(e) {
        NULL
      }
  )

  if (
    !is.null(
      basemap_image
    )
  ) {

    # Fade the basemap toward white so analytical layers dominate.
    rgb_channels <-
      seq_len(
        min(
          3,
          dim(
            basemap_image
          )[3]
        )
      )

    basemap_image[
      ,
      ,
      rgb_channels
    ] <-
      1 -
      0.42 *
      (
        1 -
          basemap_image[
            ,
            ,
            rgb_channels
          ]
      )

    if (
      dim(
        basemap_image
      )[3] ==
        4
    ) {
      basemap_image[
        ,
        ,
        4
      ] <- 1
    }
  }
}

sensors_wgs <- st_transform(
  sensors_sf,
  4326
)

buffers_wgs <- st_transform(
  buffers_sf,
  4326
)

add_city_background <- function(
    p
) {

  if (
    !is.null(
      basemap_image
    )
  ) {

    p <- p +
      annotation_raster(
        basemap_image,
        xmin =
          map_bbox["xmin"],
        xmax =
          map_bbox["xmax"],
        ymin =
          map_bbox["ymin"],
        ymax =
          map_bbox["ymax"]
      )
  } else {

    p <- p +
      geom_sf(
        data =
          boundary_sf,
        fill =
          "#F4F6F7",
        colour =
          NA
      )
  }

  p +
    coord_sf(
      xlim =
        map_bbox[
          c(
            "xmin",
            "xmax"
          )
        ],
      ylim =
        map_bbox[
          c(
            "ymin",
            "ymax"
          )
        ],
      expand =
        FALSE,
      datum =
        NA
    )
}

map_theme <- theme_void(
  base_size =
    14
) +
  theme(
    plot.title =
      element_text(
        face =
          "bold",
        size =
          15,
        colour =
          "#1F2933"
      ),
    plot.subtitle =
      element_text(
        size =
          11.5,
        colour =
          "#4B5563"
      ),
    legend.position =
      "right",
    legend.title =
      element_text(
        face =
          "bold"
      ),
    plot.caption =
      element_text(
        size =
          8.8,
        colour =
          "#687078",
        hjust =
          0
      )
  )

map_caption <-
  paste(
    "Background and municipal boundary:",
    "City of Melbourne GIS / Open Data."
  )

# ============================================================
# FIGURE 1 — SENSOR NETWORK AND 500 m CATCHMENTS
# ============================================================

p_map <- ggplot()

p_map <- add_city_background(
  p_map
) +
  geom_sf(
    data =
      buffers_wgs,
    fill =
      "#56B4E9",
    colour =
      "#0072B2",
    linewidth =
      0.25,
    alpha =
      0.11
  ) +
  geom_sf(
    data =
      sensors_wgs,
    colour =
      "#005A9C",
    fill =
      "#FFFFFF",
    shape =
      21,
    size =
      2.45,
    stroke =
      0.75
  ) +
  geom_sf(
    data = boundary_sf,
    fill = NA,
    colour = "#263640",
    linewidth = 1.15
  ) +
  labs(
    title =
      "Retained pedestrian sensors and 500 m analytical catchments",
    subtitle =
      "Final 53-sensor modelling network within the City of Melbourne",
    caption =
      map_caption
  ) +
  map_theme

save_plot(
  p_map,
  "report_sensor_catchments_map.png",
  10.5,
  8.0
)

# ============================================================
# FIGURE 2 — SPATIAL PREDICTOR CORRELATION HEATMAP
# ============================================================

spatial_cov <- read_csv(
  here(
    "data",
    "processed",
    "spatial_covariates_500m.csv"
  ),
  show_col_types =
    FALSE
)

candidate_vars <- c(
  "business_count_500m",
  "business_mix_entropy_500m",
  "dwelling_count_500m",
  "resident_population_500m",
  "income_context_500m",
  "built_coverage_500m",
  "median_structure_height_500m"
)

candidate_vars <- candidate_vars[
  candidate_vars %in%
    names(
      spatial_cov
    )
]

pretty_names <- c(
  business_count_500m =
    "Business count",
  business_mix_entropy_500m =
    "Business mix entropy",
  dwelling_count_500m =
    "Dwelling count",
  resident_population_500m =
    "Resident population",
  income_context_500m =
    "Income context",
  built_coverage_500m =
    "Built coverage",
  median_structure_height_500m =
    "Median structure height"
)

cor_mat <- cor(
  spatial_cov[
    candidate_vars
  ],
  use =
    "pairwise.complete.obs",
  method =
    "spearman"
)

cor_long <- as.data.frame(
  as.table(
    cor_mat
  )
)

names(
  cor_long
) <- c(
  "Var1",
  "Var2",
  "correlation"
)

cor_long <- cor_long %>%
  mutate(
    Var1 =
      factor(
        Var1,
        levels =
          candidate_vars,
        labels =
          unname(
            pretty_names[
              candidate_vars
            ]
          )
      ),
    Var2 =
      factor(
        Var2,
        levels =
          candidate_vars,
        labels =
          unname(
            pretty_names[
              candidate_vars
            ]
          )
      ),
    text_colour =
      abs(
        correlation
      ) >
      0.55
  )

p_heat <- ggplot(
  cor_long,
  aes(
    x =
      Var1,
    y =
      Var2,
    fill =
      correlation
  )
) +
  geom_tile(
    colour =
      "white",
    linewidth =
      0.65
  ) +
  geom_text(
    aes(
      label =
        sprintf(
          "%.2f",
          correlation
        ),
      colour =
        text_colour
    ),
    size =
      4.0,
    fontface =
      "bold"
  ) +
  scale_fill_gradient2(
    low =
      "#B2182B",
    mid =
      "#F7F7F7",
    high =
      "#2166AC",
    midpoint =
      0,
    limits =
      c(
        -1,
        1
      ),
    name =
      "Spearman\ncorrelation"
  ) +
  scale_colour_manual(
    values =
      c(
        "FALSE" =
          "#1F2933",
        "TRUE" =
          "#FFFFFF"
      ),
    guide =
      "none"
  ) +
  coord_equal() +
  labs(
    title =
      "Correlation among candidate spatial predictors",
    x =
      NULL,
    y =
      NULL
  ) +
  report_theme +
  theme(
    axis.text.x =
      element_text(
        angle =
          38,
        hjust =
          1,
        vjust =
          1,
        size =
          10
      ),
    axis.text.y =
      element_text(
        size =
          10
      ),
    panel.grid =
      element_blank()
  )

save_plot(
  p_heat,
  "report_spatial_correlation_heatmap.png",
  10.2,
  8.2
)

# ============================================================
# FIGURE 3 — DAILY COUNT DISTRIBUTION
# ============================================================

p_eda <- ggplot(
  mod_dat,
  aes(
    x =
      daily_ped_count
  )
) +
  geom_histogram(
    bins =
      45,
    fill =
      "#2C7FB8",
    colour =
      "white",
    linewidth =
      0.2
  ) +
  geom_vline(
    xintercept =
      median(
        mod_dat$daily_ped_count
      ),
    colour =
      "#D7301F",
    linewidth =
      0.95,
    linetype =
      "dashed"
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Distribution of daily pedestrian counts",
    subtitle =
      "Dashed line denotes the observed median",
    x =
      "Daily pedestrian count",
    y =
      "Sensor-days"
  ) +
  report_theme +
  theme(
    legend.position =
      "none"
  )

save_plot(
  p_eda,
  "report_eda_daily_counts.png",
  10.0,
  6.0
)

# ============================================================
# FIGURE 4 — DAY-OF-WEEK EDA
# ============================================================

weekday_labels <- c(
  "Monday",
  "Tuesday",
  "Wednesday",
  "Thursday",
  "Friday",
  "Saturday",
  "Sunday"
)

weekday_cols <- c(
  "#4477AA",
  "#66CCEE",
  "#228833",
  "#CCBB44",
  "#EE6677",
  "#AA3377",
  "#6C6C6C"
)

weekday_summary <- mod_dat %>%
  group_by(
    dow_index
  ) %>%
  summarise(
    mean_count =
      mean(
        daily_ped_count
      ),
    .groups =
      "drop"
  ) %>%
  mutate(
    day_of_week =
      factor(
        dow_index,
        levels =
          1:7,
        labels =
          weekday_labels
      )
  )

p_weekday <- ggplot(
  weekday_summary,
  aes(
    x =
      day_of_week,
    y =
      mean_count,
    group =
      1
  )
) +
  geom_line(
    colour =
      "#75808A",
    linewidth =
      0.9
  ) +
  geom_point(
    aes(
      colour =
        day_of_week
    ),
    size =
      3.9
  ) +
  scale_colour_manual(
    values =
      setNames(
        weekday_cols,
        weekday_labels
      ),
    guide =
      "none"
  ) +
  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Mean pedestrian activity by day of week",
    x =
      "Day of week",
    y =
      "Mean daily pedestrian count"
  ) +
  report_theme +
  theme(
    axis.text.x =
      element_text(
        angle =
          28,
        hjust =
          1
      )
  )

save_plot(
  p_weekday,
  "report_eda_day_of_week.png",
  9.6,
  5.9
)

# ============================================================
# FIGURE 5 — SENSOR MEAN ACTIVITY
# ============================================================

sensor_summary <- mod_dat %>%
  group_by(
    sensor_index
  ) %>%
  summarise(
    mean_count =
      mean(
        daily_ped_count
      ),
    .groups =
      "drop"
  ) %>%
  arrange(
    mean_count
  ) %>%
  mutate(
    sensor_label =
      factor(
        paste0(
          "S",
          sensor_index
        ),
        levels =
          paste0(
            "S",
            sensor_index
          )
      )
  )

p_sensor_mean <- ggplot(
  sensor_summary,
  aes(
    x =
      mean_count,
    y =
      sensor_label,
    colour =
      mean_count
  )
) +
  geom_segment(
    aes(
      x =
        0,
      xend =
        mean_count,
      yend =
        sensor_label
    ),
    colour =
      "#D8DEE3",
    linewidth =
      0.55
  ) +
  geom_point(
    size =
      2.6
  ) +
  scale_colour_gradientn(
    colours =
      c(
        "#2C7BB6",
        "#ABD9E9",
        "#FFFFBF",
        "#FDAE61",
        "#D7191C"
      ),
    guide =
      "none"
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Mean daily pedestrian activity by sensor",
    subtitle =
      "Sensors ordered from lower to higher mean activity",
    x =
      "Mean daily pedestrian count",
    y =
      "Sensor"
  ) +
  report_theme +
  theme(
    axis.text.y =
      element_text(
        size =
          7.5
      ),
    panel.grid.major.y =
      element_blank()
  )

save_plot(
  p_sensor_mean,
  "report_eda_sensor_means.png",
  9.5,
  10.5
)

# ============================================================
# FIGURE 6 — FINAL M4 PRIOR-PREDICTIVE CHECK
# ============================================================

prior_draws <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_prior_predictive_draws.csv"
  ),
  show_col_types =
    FALSE
)

prior_summary <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_prior_predictive_summary.csv"
  ),
  show_col_types =
    FALSE
)

prior_stats <- c(
  "mean_count",
  "q95_count",
  "maximum_count",
  "variance_to_mean"
)

prior_labels <- c(
  mean_count =
    "Mean daily count",
  q95_count =
    "95th percentile",
  maximum_count =
    "Maximum count",
  variance_to_mean =
    "Variance-to-mean ratio"
)

prior_long <- prior_draws %>%
  select(
    all_of(
      prior_stats
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to =
      "statistic",
    values_to =
      "value"
  ) %>%
  mutate(
    statistic =
      factor(
        statistic,
        levels =
          prior_stats,
        labels =
          unname(
            prior_labels[
              prior_stats
            ]
          )
      )
  )

prior_obs <- prior_summary %>%
  filter(
    statistic %in%
      prior_stats
  ) %>%
  transmute(
    statistic =
      factor(
        statistic,
        levels =
          prior_stats,
        labels =
          unname(
            prior_labels[
              prior_stats
            ]
          )
      ),
    observed =
      observed
  )

prior_cols <- setNames(
  c(
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#B279A2"
  ),
  unname(
    prior_labels[
      prior_stats
    ]
  )
)

p_prior <- ggplot(
  prior_long,
  aes(
    x =
      value,
    fill =
      statistic
  )
) +
  geom_histogram(
    bins =
      35,
    colour =
      "white",
    linewidth =
      0.18,
    alpha =
      0.92
  ) +
  geom_vline(
    data =
      prior_obs,
    aes(
      xintercept =
        observed
    ),
    inherit.aes =
      FALSE,
    colour =
      "#B2182B",
    linewidth =
      1.05,
    linetype =
      "dashed"
  ) +
  facet_wrap(
    ~ statistic,
    scales =
      "free_x",
    ncol =
      2
  ) +
  scale_fill_manual(
    values =
      prior_cols
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark =
          ",",
        accuracy =
          1
      )
  ) +
  labs(
    title =
      "M4 prior-predictive distributions",
    subtitle =
      "Dashed red lines show the corresponding observed statistics",
    x =
      "Prior-predictive statistic",
    y =
      "Prior draws"
  ) +
  report_theme +
  theme(
    legend.position =
      "none"
  )

save_plot(
  p_prior,
  "report_m4_prior_predictive.png",
  10.8,
  7.7
)

# ============================================================
# FIGURE 7 — M4 MCMC TRACE PLOT
# ============================================================

m4_rds <- here(
  "outputs",
  "mcmc",
  "m4_sensor_day_negbin_samples.rds"
)

stopifnot(
  file.exists(
    m4_rds
  )
)

m4_samples <- readRDS(
  m4_rds
)

trace_parameters <- c(
  "alpha",
  "beta_holiday",
  "beta_temp",
  "beta_rain",
  "phi",
  "sigma_sensor",
  "sigma_day"
)

trace_df <- bind_rows(
  lapply(
    seq_along(
      m4_samples
    ),
    function(
      ch
    ) {

      mm <- as.matrix(
        m4_samples[[ch]]
      )

      available <- intersect(
        trace_parameters,
        colnames(
          mm
        )
      )

      out <- as.data.frame(
        mm[
          ,
          available,
          drop =
            FALSE
        ]
      )

      out$iteration <-
        seq_len(
          nrow(
            out
          )
        )

      out$chain <-
        factor(
          ch
        )

      out
    }
  )
) %>%
  filter(
    iteration ==
      1L |
      iteration %%
        4L ==
        0L
  ) %>%
  pivot_longer(
    cols =
      any_of(
        trace_parameters
      ),
    names_to =
      "parameter",
    values_to =
      "value"
  )

trace_labels <- c(
  alpha =
    "Intercept (alpha)",
  beta_holiday =
    "Public holiday",
  beta_temp =
    "Mean temperature",
  beta_rain =
    "Precipitation",
  phi =
    "NB dispersion (phi)",
  sigma_sensor =
    "Sensor SD",
  sigma_day =
    "Day SD"
)

trace_df <- trace_df %>%
  mutate(
    parameter =
      factor(
        parameter,
        levels =
          trace_parameters,
        labels =
          unname(
            trace_labels[
              trace_parameters
            ]
          )
      )
  )

p_trace <- ggplot(
  trace_df,
  aes(
    x =
      iteration,
    y =
      value,
    colour =
      chain
  )
) +
  geom_line(
    linewidth =
      0.52,
    alpha =
      0.84
  ) +
  facet_wrap(
    ~ parameter,
    scales =
      "free_y",
    ncol =
      1
  ) +
  scale_colour_manual(
    values =
      cols4,
    name =
      "Chain"
  ) +
  labs(
    title =
      "M4 representative MCMC trace plots",
    x =
      "Retained iteration",
    y =
      "Parameter value"
  ) +
  report_theme +
  theme(
    legend.position =
      "top",
    panel.spacing.y =
      grid::unit(
        0.65,
        "lines"
      )
  )

save_plot(
  p_trace,
  "report_mcmc_trace_m4.png",
  10.5,
  12.8
)

# ============================================================
# FIGURE 8 — CONFIRMATION PARAMETER RECOVERY
# ============================================================

core_rec <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_core.csv"
  ),
  show_col_types =
    FALSE
)

param_label <- c(
  alpha =
    "Intercept",
  "beta_dow[2]" =
    "Tuesday",
  "beta_dow[3]" =
    "Wednesday",
  "beta_dow[4]" =
    "Thursday",
  "beta_dow[5]" =
    "Friday",
  "beta_dow[6]" =
    "Saturday",
  "beta_dow[7]" =
    "Sunday",
  beta_holiday =
    "Public holiday",
  beta_temp =
    "Temperature",
  beta_rain =
    "Precipitation",
  beta_business =
    "Business count",
  beta_population =
    "Resident population",
  beta_coverage =
    "Built coverage",
  beta_height =
    "Structure height",
  phi =
    "NB dispersion",
  sigma_sensor =
    "Sensor SD",
  sigma_day =
    "Day SD"
)

core_rec <- core_rec %>%
  mutate(
    group =
      case_when(
        parameter ==
          "alpha" ~
          "Intercept",
        grepl(
          "^beta_dow",
          parameter
        ) |
          parameter %in%
          c(
            "beta_holiday",
            "beta_temp",
            "beta_rain"
          ) ~
          "Temporal / weather",
        parameter %in%
          c(
            "beta_business",
            "beta_population",
            "beta_coverage",
            "beta_height"
          ) ~
          "Urban context",
        parameter ==
          "phi" ~
          "NB dispersion",
        TRUE ~
          "Random-effect SDs"
      ),
    label =
      unname(
        param_label[
          parameter
        ]
      )
  )

group_cols <- c(
  "Intercept" =
    "#4D4D4D",
  "Temporal / weather" =
    "#0072B2",
  "Urban context" =
    "#009E73",
  "NB dispersion" =
    "#D55E00",
  "Random-effect SDs" =
    "#CC79A7"
)

p_core <- ggplot(
  core_rec,
  aes(
    y =
      reorder(
        label,
        posterior_mean
      ),
    x =
      posterior_mean,
    colour =
      group
  )
) +
  geom_errorbar(
    aes(
      xmin =
        q025,
      xmax =
        q975
    ),
    orientation =
      "y",
    width =
      0.18,
    linewidth =
      1.0
  ) +
  geom_point(
    size =
      3.0
  ) +
  geom_point(
    aes(
      x =
        true_value
    ),
    inherit.aes =
      TRUE,
    shape =
      23,
    fill =
      "white",
    colour =
      "#111111",
    size =
      2.6,
    stroke =
      0.75
  ) +
  facet_wrap(
    ~ group,
    scales =
      "free_x",
    ncol =
      2
  ) +
  scale_colour_manual(
    values =
      group_cols,
    guide =
      "none"
  ) +
  labs(
    title =
      "M4 confirmation parameter recovery",
    subtitle =
      "Coloured points/intervals are posterior summaries; hollow diamonds are true simulated values",
    x =
      "Parameter value",
    y =
      NULL
  ) +
  report_theme +
  theme(
    panel.grid.major.y =
      element_blank()
  )

save_plot(
  p_core,
  "report_m4_recovery_core.png",
  11.0,
  8.6
)

# ============================================================
# FIGURES 9 AND 10 — RANDOM-EFFECT RECOVERY
# ============================================================

sensor_rec <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_sensor_effects.csv"
  ),
  show_col_types =
    FALSE
)

p_sensor_rec <- ggplot(
  sensor_rec,
  aes(
    x =
      true_effect,
    y =
      posterior_mean
  )
) +
  geom_abline(
    slope =
      1,
    intercept =
      0,
    colour =
      "#41484D",
    linewidth =
      0.9,
    linetype =
      "dashed"
  ) +
  geom_point(
    aes(
      colour =
        true_effect
    ),
    size =
      2.6,
    alpha =
      0.9
  ) +
  scale_colour_gradient2(
    low =
      "#2166AC",
    mid =
      "#F7F7F7",
    high =
      "#B2182B",
    midpoint =
      0,
    name =
      "True effect"
  ) +
  coord_equal() +
  labs(
    title =
      "Sensor random-effect recovery",
    x =
      "True sensor effect",
    y =
      "Posterior mean"
  ) +
  report_theme

save_plot(
  p_sensor_rec,
  "report_m4_recovery_sensor.png",
  7.8,
  6.6
)

day_rec <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_day_effects.csv"
  ),
  show_col_types =
    FALSE
)

weekday_lookup <- mod_dat %>%
  distinct(
    date,
    dow_index
  ) %>%
  mutate(
    day_of_week =
      factor(
        dow_index,
        levels =
          1:7,
        labels =
          weekday_labels
      ),
    date =
      as.character(
        date
      )
  ) %>%
  select(
    date,
    day_of_week
  )

day_rec <- day_rec %>%
  mutate(
    date =
      as.character(
        date
      )
  ) %>%
  left_join(
    weekday_lookup,
    by =
      "date"
  )

p_day_rec <- ggplot(
  day_rec,
  aes(
    x =
      true_effect,
    y =
      posterior_mean,
    colour =
      day_of_week
  )
) +
  geom_abline(
    slope =
      1,
    intercept =
      0,
    colour =
      "#41484D",
    linewidth =
      0.9,
    linetype =
      "dashed"
  ) +
  geom_point(
    size =
      2.05,
    alpha =
      0.80
  ) +
  scale_colour_manual(
    values =
      setNames(
        cols7,
        weekday_labels
      ),
    name =
      "Day of week"
  ) +
  coord_equal() +
  labs(
    title =
      "Day random-effect recovery",
    x =
      "True day effect",
    y =
      "Posterior mean"
  ) +
  report_theme

save_plot(
  p_day_rec,
  "report_m4_recovery_day.png",
  8.6,
  7.0
)

# ============================================================
# FIGURE 11 — CORRECT DIC COMPARISON
# ============================================================

model_comp <- read_csv(
  here(
    "outputs",
    "tables",
    "model_comparison_final.csv"
  ),
  show_col_types =
    FALSE
)

model_comp <- model_comp %>%
  mutate(
    model_short =
      sub(
        "\\s.*$",
        "",
        model
      ),
    model_short =
      factor(
        model_short,
        levels =
          c(
            "M1",
            "M2",
            "M3",
            "M4"
          )
      ),
    dic_label =
      format(
        round(
          DIC
        ),
        big.mark =
          ",",
        scientific =
          FALSE
      )
  )

model_cols <- c(
  "M1" =
    "#7A7A7A",
  "M2" =
    "#56B4E9",
  "M3" =
    "#009E73",
  "M4" =
    "#D55E00"
)

p_dic <- ggplot(
  model_comp,
  aes(
    x =
      DIC,
    y =
      model_short,
    colour =
      model_short
  )
) +
  geom_segment(
    aes(
      x =
        300000,
      xend =
        DIC,
      yend =
        model_short
    ),
    linewidth =
      1.1,
    alpha =
      0.55
  ) +
  geom_point(
    size =
      4.2
  ) +
  geom_text(
    aes(
      label =
        dic_label
    ),
    nudge_y =
      0.18,
    hjust =
      0,
    size =
      3.8,
    fontface =
      "bold",
    show.legend =
      FALSE
  ) +
  scale_x_log10(
    labels =
      scales::label_number(
        big.mark =
          ",",
        accuracy =
          1
      ),
    expand =
      expansion(
        mult =
          c(
            0.05,
            0.30
          )
      )
  ) +
  scale_colour_manual(
    values =
      model_cols,
    guide =
      "none"
  ) +
  labs(
    title =
      "DIC across the four Bayesian count models",
    subtitle =
      "Logarithmic x-axis; lower DIC indicates better relative fit",
    x =
      "DIC (log scale)",
    y =
      NULL
  ) +
  report_theme +
  theme(
    panel.grid.major.y =
      element_blank()
  )

save_plot(
  p_dic,
  "report_model_comparison_dic.png",
  9.3,
  5.8
)

# ============================================================
# FIGURES 12-15 — CONSISTENT M3 / M4 PPC VISUALS
# ============================================================

mcmc_dir <- here(
  "outputs",
  "mcmc"
)

find_mcmc_file <- function(
    pattern
) {

  candidates <- list.files(
    mcmc_dir,
    pattern =
      pattern,
    full.names =
      TRUE,
    ignore.case =
      TRUE
  )

  candidates <- candidates[
    !grepl(
      "recover|parameter",
      basename(
        candidates
      ),
      ignore.case =
        TRUE
    )
  ]

  if (
    length(
      candidates
    ) ==
      0
  ) {
    stop(
      paste(
        "Could not locate MCMC file matching:",
        pattern
      )
    )
  }

  candidates[
    1
  ]
}

m3_file <- find_mcmc_file(
  "m3.*\\.rds$"
)

m4_file <- if (
  file.exists(
    m4_rds
  )
) {
  m4_rds
} else {
  find_mcmc_file(
    "m4.*\\.rds$"
  )
}

m3_mat <- as.matrix(
  readRDS(
    m3_file
  )
)

m4_mat <- as.matrix(
  readRDS(
    m4_file
  )
)

extract_indexed <- function(
    draw,
    prefix,
    n
) {

  out <- numeric(
    n
  )

  for (
    i in seq_len(
      n
    )
  ) {

    nm <- paste0(
      prefix,
      "[",
      i,
      "]"
    )

    if (
      nm %in%
        names(
          draw
        )
    ) {
      out[
        i
      ] <-
        draw[
          nm
        ]
    }
  }

  out
}

get_draw_value <- function(
    draw,
    name,
    default =
      0
) {

  if (
    name %in%
      names(
        draw
      )
  ) {
    as.numeric(
      draw[
        name
      ]
    )
  } else {
    default
  }
}

observed_y <-
  mod_dat$daily_ped_count

sensor_idx <-
  mod_dat$sensor_index

day_idx <-
  mod_dat$day_index

dow_idx <-
  mod_dat$dow_index

observed_sensor_mean <- as.numeric(
  tapply(
    observed_y,
    sensor_idx,
    mean
  )
)

observed_date_mean <- as.numeric(
  tapply(
    observed_y,
    day_idx,
    mean
  )
)

quantile_probs <- c(
  0.01,
  0.05,
  0.10,
  0.25,
  0.50,
  0.75,
  0.90,
  0.95,
  0.99
)

observed_quantiles <- as.numeric(
  quantile(
    observed_y,
    quantile_probs,
    names =
      FALSE
  )
)

calculate_global_stats <- function(
    y
) {

  sensor_means <- as.numeric(
    tapply(
      y,
      sensor_idx,
      mean
    )
  )

  date_means <- as.numeric(
    tapply(
      y,
      day_idx,
      mean
    )
  )

  c(
    mean_count =
      mean(
        y
      ),
    median_count =
      median(
        y
      ),
    sd_count =
      sd(
        y
      ),
    variance_to_mean =
      var(
        y
      ) /
      mean(
        y
      ),
    q95_count =
      as.numeric(
        quantile(
          y,
          0.95
        )
      ),
    maximum_count =
      max(
        y
      ),
    proportion_gt_50000 =
      mean(
        y >
          50000
      ),
    sd_sensor_mean =
      sd(
        sensor_means
      ),
    max_sensor_mean =
      max(
        sensor_means
      ),
    sd_date_mean =
      sd(
        date_means
      )
  )
}

simulate_ppc <- function(
    mat,
    method,
    include_day,
    n_draws =
      200
) {

  draw_ids <- sample(
    seq_len(
      nrow(
        mat
      )
    ),
    size =
      min(
        n_draws,
        nrow(
          mat
        )
      ),
    replace =
      FALSE
  )

  n_use <- length(
    draw_ids
  )

  global_list <- vector(
    "list",
    n_use
  )

  sensor_matrix <- matrix(
    NA_real_,
    nrow =
      n_use,
    ncol =
      J
  )

  date_matrix <- matrix(
    NA_real_,
    nrow =
      n_use,
    ncol =
      D
  )

  quantile_matrix <- matrix(
    NA_real_,
    nrow =
      n_use,
    ncol =
      length(
        quantile_probs
      )
  )

  for (
    ii in seq_along(
      draw_ids
    )
  ) {

    draw <-
      mat[
        draw_ids[
          ii
        ],
        ,
        drop =
          TRUE
      ]

    eta <-
      rep(
        get_draw_value(
          draw,
          "alpha"
        ),
        N
      )

    eta <-
      eta +
      get_draw_value(
        draw,
        "beta_holiday"
      ) *
      mod_dat$public_holiday_num +
      get_draw_value(
        draw,
        "beta_temp"
      ) *
      mod_dat$z_mean_temp +
      get_draw_value(
        draw,
        "beta_rain"
      ) *
      mod_dat$z_log1p_precipitation +
      get_draw_value(
        draw,
        "beta_business"
      ) *
      mod_dat$z_business_count +
      get_draw_value(
        draw,
        "beta_population"
      ) *
      mod_dat$z_resident_population +
      get_draw_value(
        draw,
        "beta_coverage"
      ) *
      mod_dat$z_built_coverage +
      get_draw_value(
        draw,
        "beta_height"
      ) *
      mod_dat$z_structure_height

    for (
      d in 2:7
    ) {

      eta[
        dow_idx ==
          d
      ] <-
        eta[
          dow_idx ==
            d
        ] +
        get_draw_value(
          draw,
          paste0(
            "beta_dow[",
            d,
            "]"
          )
        )
    }

    u_sensor <-
      extract_indexed(
        draw,
        "u_sensor",
        J
      )

    eta <-
      eta +
      u_sensor[
        sensor_idx
      ]

    if (
      include_day
    ) {

      day_effect <-
        extract_indexed(
          draw,
          "day_effect",
          D
        )

      eta <-
        eta +
        day_effect[
          day_idx
        ]
    }

    mu <-
      exp(
        eta
      )

    phi <-
      get_draw_value(
        draw,
        "phi"
      )

    y_rep <- rnbinom(
      N,
      size =
        phi,
      mu =
        mu
    )

    global_list[[ii]] <-
      tibble(
        draw =
          ii,
        method =
          method,
        statistic =
          names(
            calculate_global_stats(
              y_rep
            )
          ),
        value =
          as.numeric(
            calculate_global_stats(
              y_rep
            )
          )
      )

    sensor_matrix[
      ii,
    ] <-
      as.numeric(
        tapply(
          y_rep,
          sensor_idx,
          mean
        )
      )

    date_matrix[
      ii,
    ] <-
      as.numeric(
        tapply(
          y_rep,
          day_idx,
          mean
        )
      )

    quantile_matrix[
      ii,
    ] <-
      as.numeric(
        quantile(
          y_rep,
          quantile_probs,
          names =
            FALSE
        )
      )
  }

  list(
    global =
      bind_rows(
        global_list
      ),
    sensor =
      sensor_matrix,
    date =
      date_matrix,
    quantile =
      quantile_matrix
  )
}

set.seed(
  20260918
)

ppc_m3 <- simulate_ppc(
  m3_mat,
  "M3",
  FALSE,
  200
)

set.seed(
  20260919
)

ppc_m4 <- simulate_ppc(
  m4_mat,
  "M4",
  TRUE,
  200
)

# ------------------------------------------------------------
# Figure 12 — Global PPC statistics
# ------------------------------------------------------------

global_ppc <- bind_rows(
  ppc_m3$global,
  ppc_m4$global
)

global_obs_vec <- calculate_global_stats(
  observed_y
)

global_labels <- c(
  mean_count =
    "Mean",
  median_count =
    "Median",
  sd_count =
    "SD",
  variance_to_mean =
    "Variance / mean",
  q95_count =
    "95th percentile",
  maximum_count =
    "Maximum",
  proportion_gt_50000 =
    "Proportion > 50,000",
  sd_sensor_mean =
    "SD of sensor means",
  max_sensor_mean =
    "Maximum sensor mean",
  sd_date_mean =
    "SD of date means"
)

global_summary <- global_ppc %>%
  group_by(
    method,
    statistic
  ) %>%
  summarise(
    q025 =
      quantile(
        value,
        0.025
      ),
    median =
      median(
        value
      ),
    q975 =
      quantile(
        value,
        0.975
      ),
    .groups =
      "drop"
  ) %>%
  mutate(
    label =
      factor(
        statistic,
        levels =
          names(
            global_labels
          ),
        labels =
          unname(
            global_labels
          )
      ),
    method =
      factor(
        method,
        levels =
          c(
            "M3",
            "M4"
          )
      )
  )

global_obs <- tibble(
  statistic =
    names(
      global_obs_vec
    ),
  observed =
    as.numeric(
      global_obs_vec
    )
) %>%
  mutate(
    label =
      factor(
        statistic,
        levels =
          names(
            global_labels
          ),
        labels =
          unname(
            global_labels
          )
      )
  )

p_ppc_global <- ggplot(
  global_summary,
  aes(
    x =
      median,
    y =
      method,
    colour =
      method
  )
) +
  geom_vline(
    data =
      global_obs,
    aes(
      xintercept =
        observed
    ),
    inherit.aes =
      FALSE,
    colour =
      "#222222",
    linewidth =
      0.75,
    linetype =
      "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin =
        q025,
      xmax =
        q975
    ),
    orientation =
      "y",
    width =
      0.20,
    linewidth =
      0.9
  ) +
  geom_point(
    size =
      2.8
  ) +
  facet_wrap(
    ~ label,
    scales =
      "free_x",
    ncol =
      2
  ) +
  scale_colour_manual(
    values =
      cols_model,
    guide =
      "none"
  ) +
  labs(
    title =
      "M3 and M4 posterior predictive checks of global statistics",
    subtitle =
      "Coloured intervals are posterior-predictive 95% intervals; dashed lines are observed statistics",
    x =
      "Statistic value",
    y =
      NULL
  ) +
  report_theme

save_plot(
  p_ppc_global,
  "report_ppc_m3_m4_statistics.png",
  11.2,
  9.6
)

# ------------------------------------------------------------
# Figure 13 — Sensor means
# ------------------------------------------------------------

summarise_matrix <- function(
    x,
    method,
    unit_name
) {

  bind_rows(
    lapply(
      seq_len(
        ncol(
          x
        )
      ),
      function(
        j
      ) {

        tibble(
          unit =
            j,
          method =
            method,
          q025 =
            quantile(
              x[
                ,
                j
              ],
              0.025
            ),
          median =
            median(
              x[
                ,
                j
              ]
            ),
          q975 =
            quantile(
              x[
                ,
                j
              ],
              0.975
            )
        )
      }
    )
  ) %>%
    rename(
      !!unit_name :=
        unit
    )
}

sensor_ppc <- bind_rows(
  summarise_matrix(
    ppc_m3$sensor,
    "M3",
    "sensor_index"
  ),
  summarise_matrix(
    ppc_m4$sensor,
    "M4",
    "sensor_index"
  )
) %>%
  mutate(
    observed =
      observed_sensor_mean[
        sensor_index
      ],
    method =
      factor(
        method,
        levels =
          c(
            "M3",
            "M4"
          )
      )
  )

sensor_order <- order(
  observed_sensor_mean
)

sensor_rank <- integer(
  J
)

sensor_rank[
  sensor_order
] <-
  seq_len(
    J
  )

sensor_ppc <- sensor_ppc %>%
  mutate(
    sensor_rank =
      sensor_rank[
        sensor_index
      ]
  )

p_ppc_sensor <- ggplot(
  sensor_ppc,
  aes(
    x =
      sensor_rank,
    colour =
      method
  )
) +
  geom_linerange(
    aes(
      ymin =
        q025,
      ymax =
        q975
    ),
    linewidth =
      0.62,
    alpha =
      0.78
  ) +
  geom_point(
    aes(
      y =
        median
    ),
    size =
      1.45
  ) +
  geom_point(
    aes(
      y =
        observed
    ),
    colour =
      "#111111",
    shape =
      4,
    size =
      1.7,
    stroke =
      0.65
  ) +
  facet_wrap(
    ~ method,
    ncol =
      1
  ) +
  scale_colour_manual(
    values =
      cols_model,
    guide =
      "none"
  ) +
  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Posterior predictive checks of sensor-level mean activity",
    subtitle =
      "Sensors are ordered by observed mean; black crosses are observed means",
    x =
      "Sensors ordered by observed mean",
    y =
      "Mean daily pedestrian count"
  ) +
  report_theme

save_plot(
  p_ppc_sensor,
  "report_ppc_m3_m4_sensor_means.png",
  11.0,
  8.6
)

# ------------------------------------------------------------
# Figure 14 — Date means
# ------------------------------------------------------------

date_ppc <- bind_rows(
  summarise_matrix(
    ppc_m3$date,
    "M3",
    "day_index"
  ),
  summarise_matrix(
    ppc_m4$date,
    "M4",
    "day_index"
  )
) %>%
  mutate(
    date =
      date_lookup$date[
        day_index
      ],
    observed =
      observed_date_mean[
        day_index
      ],
    method =
      factor(
        method,
        levels =
          c(
            "M3",
            "M4"
          )
      )
  )

p_ppc_date <- ggplot(
  date_ppc,
  aes(
    x =
      date
  )
) +
  geom_ribbon(
    aes(
      ymin =
        q025,
      ymax =
        q975,
      fill =
        method
    ),
    alpha =
      0.18,
    colour =
      NA
  ) +
  geom_line(
    aes(
      y =
        median,
      colour =
        method
    ),
    linewidth =
      0.7
  ) +
  geom_line(
    aes(
      y =
        observed
    ),
    colour =
      "#202020",
    linewidth =
      0.55,
    alpha =
      0.85
  ) +
  facet_wrap(
    ~ method,
    ncol =
      1
  ) +
  scale_colour_manual(
    values =
      cols_model,
    guide =
      "none"
  ) +
  scale_fill_manual(
    values =
      cols_model,
    guide =
      "none"
  ) +
  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Posterior predictive checks of shared date-level activity",
    subtitle =
      "Black line is observed; coloured line and band are posterior-predictive median and 95% interval",
    x =
      "Date",
    y =
      "Mean pedestrian count across sensors"
  ) +
  report_theme

save_plot(
  p_ppc_date,
  "report_ppc_m3_m4_date_means.png",
  11.2,
  8.1
)

# ------------------------------------------------------------
# Figure 15 — Quantiles
# ------------------------------------------------------------

quantile_summary <- bind_rows(
  lapply(
    c(
      "M3",
      "M4"
    ),
    function(
      method
    ) {

      mat <- if (
        method ==
          "M3"
      ) {
        ppc_m3$quantile
      } else {
        ppc_m4$quantile
      }

      bind_rows(
        lapply(
          seq_along(
            quantile_probs
          ),
          function(
            j
          ) {

            tibble(
              method =
                method,
              probability =
                quantile_probs[
                  j
                ],
              observed =
                observed_quantiles[
                  j
                ],
              q025 =
                quantile(
                  mat[
                    ,
                    j
                  ],
                  0.025
                ),
              median =
                median(
                  mat[
                    ,
                    j
                  ]
                ),
              q975 =
                quantile(
                  mat[
                    ,
                    j
                  ],
                  0.975
                )
            )
          }
        )
      )
    }
  )
) %>%
  mutate(
    method =
      factor(
        method,
        levels =
          c(
            "M3",
            "M4"
          )
      )
  )

p_ppc_quant <- ggplot(
  quantile_summary,
  aes(
    x =
      observed,
    y =
      median,
    colour =
      method
  )
) +
  geom_abline(
    slope =
      1,
    intercept =
      0,
    colour =
      "#3D4247",
    linewidth =
      0.8,
    linetype =
      "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin =
        q025,
      ymax =
        q975
    ),
    width =
      0,
    linewidth =
      0.8
  ) +
  geom_point(
    size =
      2.8
  ) +
  facet_wrap(
    ~ method,
    ncol =
      1
  ) +
  scale_colour_manual(
    values =
      cols_model,
    guide =
      "none"
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  scale_y_continuous(
    labels =
      scales::label_number(
        big.mark =
          ","
      )
  ) +
  labs(
    title =
      "Posterior predictive comparison of count quantiles",
    subtitle =
      "Dashed identity line indicates perfect agreement between observed and replicated quantiles",
    x =
      "Observed quantile",
    y =
      "Posterior-predictive quantile"
  ) +
  report_theme

save_plot(
  p_ppc_quant,
  "report_ppc_m3_m4_quantiles.png",
  9.5,
  8.2
)

# ============================================================
# FIGURE 16 — FINAL M4 EFFECT ESTIMATES
# ============================================================

effects <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_final_effect_estimates.csv"
  ),
  show_col_types =
    FALSE
)

effects <- effects %>%
  mutate(
    group =
      case_when(
        grepl(
          "Monday|holiday",
          effect,
          ignore.case =
            TRUE
        ) ~
          "Calendar",
        grepl(
          "temperature|precipitation",
          effect,
          ignore.case =
            TRUE
        ) ~
          "Weather",
        TRUE ~
          "Urban context"
      ),
    effect =
      factor(
        effect,
        levels =
          rev(
            effect
          )
      )
  )

effect_cols <- c(
  "Calendar" =
    "#0072B2",
  "Weather" =
    "#D55E00",
  "Urban context" =
    "#009E73"
)

p_effects <- ggplot(
  effects,
  aes(
    x =
      percent_change_median,
    y =
      effect,
    colour =
      group
  )
) +
  geom_vline(
    xintercept =
      0,
    colour =
      "#666666",
    linewidth =
      0.75,
    linetype =
      "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin =
        percent_change_q025,
      xmax =
        percent_change_q975
    ),
    orientation =
      "y",
    width =
      0.18,
    linewidth =
      0.9
  ) +
  geom_point(
    size =
      2.9
  ) +
  scale_colour_manual(
    values =
      effect_cols,
    name =
      NULL
  ) +
  labs(
    title =
      "M4 posterior associations with expected pedestrian activity",
    subtitle =
      "Posterior median percentage change and 95% credible interval",
    x =
      "Percentage change in conditional expected daily count",
    y =
      NULL
  ) +
  report_theme

save_plot(
  p_effects,
  "report_m4_effect_estimates.png",
  10.5,
  7.5
)

# ============================================================
# FIGURE 17 — CONDITIONAL SCENARIO PREDICTIONS
# ============================================================

scenarios <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_scenario_predictions.csv"
  ),
  show_col_types =
    FALSE
)

scenarios <- scenarios %>%
  mutate(
    group =
      case_when(
        scenario ==
          "Reference Monday" ~
          "Reference",
        grepl(
          "Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday|holiday",
          scenario,
          ignore.case =
            TRUE
        ) ~
          "Calendar",
        grepl(
          "Temperature|Precipitation",
          scenario,
          ignore.case =
            TRUE
        ) ~
          "Weather",
        TRUE ~
          "Urban context"
      ),
    scenario =
      factor(
        scenario,
        levels =
          rev(
            scenario
          )
      )
  )

scenario_cols <- c(
  "Reference" =
    "#4D4D4D",
  "Calendar" =
    "#0072B2",
  "Weather" =
    "#D55E00",
  "Urban context" =
    "#009E73"
)

p_scenario <- ggplot(
  scenarios,
  aes(
    x =
      expected_count_median,
    y =
      scenario,
    colour =
      group
  )
) +
  geom_errorbar(
    aes(
      xmin =
        expected_count_q025,
      xmax =
        expected_count_q975
    ),
    orientation =
      "y",
    width =
      0.18,
    linewidth =
      0.9
  ) +
  geom_point(
    size =
      2.9
  ) +
  scale_colour_manual(
    values =
      scenario_cols,
    name =
      NULL
  ) +
  scale_x_continuous(
    labels =
      scales::label_number(
        big.mark =
          ",",
        accuracy =
          1
      )
  ) +
  labs(
    title =
      "Conditional scenario expectations for a reference sensor-day",
    subtitle =
      "Sensor and day random effects fixed at zero",
    x =
      "Expected daily pedestrian count",
    y =
      NULL
  ) +
  report_theme

save_plot(
  p_scenario,
  "report_m4_scenario_predictions.png",
  10.3,
  7.5
)

# ============================================================
# FIGURE 18 — SENSOR PREDICTION MAP WITH CITY BASEMAP
# ============================================================

sensor_pred <- read_csv(
  here(
    "outputs",
    "tables",
    "m4_sensor_reference_predictions.csv"
  ),
  show_col_types =
    FALSE
)

# Detect the posterior median expected-count column robustly.
# Step 19 may use slightly different naming conventions depending on
# the version of the inference script, so prefer known names and then
# fall back to an unambiguous "median" column.

known_prediction_cols <- c(
  "expected_count_median",
  "reference_expected_count_median",
  "posterior_median_expected_count",
  "posterior_expected_count_median",
  "median_expected_count",
  "expected_median",
  "posterior_median",
  "median"
)

prediction_col <- known_prediction_cols[
  known_prediction_cols %in%
    names(
      sensor_pred
    )
][1]

if (
  is.na(
    prediction_col
  )
) {

  median_candidates <- names(
    sensor_pred
  )[
    grepl(
      "median",
      names(
        sensor_pred
      ),
      ignore.case = TRUE
    )
  ]

  # Prefer a median column whose name also indicates an expected or
  # predicted count. This avoids accidentally selecting an unrelated
  # median covariate if such a column is ever added later.
  preferred_candidates <- median_candidates[
    grepl(
      "expected|predict|count|reference",
      median_candidates,
      ignore.case = TRUE
    )
  ]

  if (
    length(
      preferred_candidates
    ) ==
      1
  ) {

    prediction_col <- preferred_candidates[1]

  } else if (
    length(
      median_candidates
    ) ==
      1
  ) {

    prediction_col <- median_candidates[1]

  } else {

    stop(
      paste0(
        "Could not identify a unique posterior median expected-count ",
        "column in m4_sensor_reference_predictions.csv. Available columns: ",
        paste(
          names(
            sensor_pred
          ),
          collapse = ", "
        )
      )
    )
  }
}

message(
  "Figure 18 using sensor prediction column: ",
  prediction_col
)

if (
  "location_id" %in%
    names(
      sensor_pred
    )
) {

  sensor_pred <- sensor_pred %>%
    mutate(
      location_id =
        as.character(
          location_id
        )
    )

  pred_sf <- sensors_wgs %>%
    left_join(
      sensor_pred,
      by =
        "location_id"
    )
} else {

  pred_sf <- sensors_wgs %>%
    left_join(
      sensor_pred,
      by =
        "sensor_index"
    )
}

names(
  pred_sf
)[
  names(
    pred_sf
  ) ==
    prediction_col
] <- "prediction_median"

p_pred_map <- ggplot()

p_pred_map <- add_city_background(
  p_pred_map
) +
  geom_sf(
    data =
      pred_sf,
    aes(
      colour =
        prediction_median,
      size =
        prediction_median
    ),
    alpha =
      0.88
  ) +
  geom_sf(
    data = boundary_sf,
    fill = NA,
    colour = "#263640",
    linewidth = 1.15,
    inherit.aes = FALSE
  ) +
  scale_colour_gradientn(
    colours =
      c(
        "#2C7BB6",
        "#ABD9E9",
        "#FFFFBF",
        "#FDAE61",
        "#D7191C"
      ),
    labels =
      scales::label_number(
        big.mark =
          ","
      ),
    name =
      "Expected\ncount"
  ) +
  scale_size_continuous(
    range =
      c(
        2.1,
        5.4
      ),
    labels =
      scales::label_number(
        big.mark =
          ","
      ),
    name =
      "Expected\ncount"
  ) +
  labs(
    title =
      "M4 predicted baseline pedestrian activity by sensor",
    subtitle =
      "Posterior median expected daily count: Monday, non-holiday, average weather",
    caption =
      map_caption
  ) +
  map_theme

save_plot(
  p_pred_map,
  "report_m4_sensor_prediction_map.png",
  10.5,
  8.0
)

# ============================================================
# FIGURE 19 — MORAN'S I POSTERIOR DISTRIBUTIONS
# ============================================================

m4_matrix <- as.matrix(
  m4_samples
)

u_names <- grep(
  "^u_sensor\\[",
  colnames(
    m4_matrix
  ),
  value =
    TRUE
)

u_index <- as.integer(
  sub(
    ".*\\[([0-9]+)\\].*",
    "\\1",
    u_names
  )
)

u_names <- u_names[
  order(
    u_index
  )
]

u_draws <- m4_matrix[
  ,
  u_names,
  drop =
    FALSE
]

sensor_projected <- st_transform(
  sensors_sf,
  7855
) %>%
  arrange(
    sensor_index
  )

coords <- st_coordinates(
  sensor_projected
)

make_knn_weights <- function(
    coordinates,
    k
) {

  dist_mat <- as.matrix(
    dist(
      coordinates
    )
  )

  diag(
    dist_mat
  ) <- Inf

  n <- nrow(
    dist_mat
  )

  A <- matrix(
    0,
    nrow =
      n,
    ncol =
      n
  )

  for (
    i in seq_len(
      n
    )
  ) {

    neighbours <-
      order(
        dist_mat[
          i,
        ]
      )[
        seq_len(
          k
        )
      ]

    A[
      i,
      neighbours
    ] <-
      1
  }

  A <-
    (
      A +
        t(
          A
        )
    ) >
    0

  W <-
    A /
    rowSums(
      A
    )

  W
}

moran_i <- function(
    x,
    W
) {

  z <-
    x -
    mean(
      x
    )

  n <-
    length(
      z
    )

  s0 <-
    sum(
      W
    )

  (
    n /
      s0
  ) *
    (
      sum(
        W *
          outer(
            z,
            z
          )
      ) /
        sum(
          z^2
        )
    )
}

set.seed(
  20260918
)

moran_draw_ids <- sample(
  seq_len(
    nrow(
      u_draws
    )
  ),
  size =
    min(
      5000,
      nrow(
        u_draws
      )
    ),
  replace =
    FALSE
)

W4 <- make_knn_weights(
  coords,
  4
)

W6 <- make_knn_weights(
  coords,
  6
)

moran_4 <- apply(
  u_draws[
    moran_draw_ids,
    ,
    drop =
      FALSE
  ],
  1,
  moran_i,
  W =
    W4
)

moran_6 <- apply(
  u_draws[
    moran_draw_ids,
    ,
    drop =
      FALSE
  ],
  1,
  moran_i,
  W =
    W6
)

moran_df <- bind_rows(
  tibble(
    neighbourhood =
      "4-NN",
    moran =
      moran_4
  ),
  tibble(
    neighbourhood =
      "6-NN",
    moran =
      moran_6
  )
) %>%
  mutate(
    neighbourhood =
      factor(
        neighbourhood,
        levels =
          c(
            "4-NN",
            "6-NN"
          )
      )
  )

observed_moran <- tibble(
  neighbourhood =
    factor(
      c(
        "4-NN",
        "6-NN"
      ),
      levels =
        c(
          "4-NN",
          "6-NN"
        )
    ),
  observed =
    c(
      moran_i(
        colMeans(
          u_draws
        ),
        W4
      ),
      moran_i(
        colMeans(
          u_draws
        ),
        W6
      )
    )
)

moran_cols <- c(
  "4-NN" =
    "#0072B2",
  "6-NN" =
    "#D55E00"
)

p_moran <- ggplot(
  moran_df,
  aes(
    x =
      moran,
    fill =
      neighbourhood
  )
) +
  geom_histogram(
    bins =
      38,
    colour =
      "white",
    linewidth =
      0.2,
    alpha =
      0.92
  ) +
  geom_vline(
    xintercept =
      0,
    colour =
      "#303030",
    linewidth =
      0.75,
    linetype =
      "dashed"
  ) +
  geom_vline(
    data =
      observed_moran,
    aes(
      xintercept =
        observed
    ),
    inherit.aes =
      FALSE,
    colour =
      "#202020",
    linewidth =
      1.0
  ) +
  facet_wrap(
    ~ neighbourhood,
    ncol =
      1
  ) +
  scale_fill_manual(
    values =
      moran_cols,
    guide =
      "none"
  ) +
  labs(
    title =
      "Posterior distribution of Moran's I for M4 sensor effects",
    subtitle =
      "Dashed line marks zero; solid black line is Moran's I of the posterior-mean sensor effects",
    x =
      "Moran's I",
    y =
      "Posterior draws"
  ) +
  report_theme

save_plot(
  p_moran,
  "report_m4_residual_moran_posterior.png",
  9.4,
  7.0
)

# ============================================================
# FIGURE 20 — SENSOR RESIDUAL MAP WITH CITY BASEMAP
# ============================================================

u_mean <- colMeans(
  u_draws
)

residual_df <- tibble(
  sensor_index =
    seq_len(
      J
    ),
  posterior_mean =
    u_mean
)

residual_sf <- sensors_wgs %>%
  left_join(
    residual_df,
    by =
      "sensor_index"
  )

max_abs_resid <- max(
  abs(
    residual_sf$posterior_mean
  ),
  na.rm =
    TRUE
)

p_resid_map <- ggplot()

p_resid_map <- add_city_background(
  p_resid_map
) +
  geom_sf(
    data =
      residual_sf,
    aes(
      colour =
        posterior_mean,
      size =
        abs(
          posterior_mean
        )
    ),
    alpha =
      0.90
  ) +
  geom_sf(
    data = boundary_sf,
    fill = NA,
    colour = "#263640",
    linewidth = 1.15,
    inherit.aes = FALSE
  ) +
  scale_colour_gradient2(
    low =
      "#2166AC",
    mid =
      "#F7F7F7",
    high =
      "#B2182B",
    midpoint =
      0,
    limits =
      c(
        -max_abs_resid,
        max_abs_resid
      ),
    name =
      "Posterior mean\nsensor effect"
  ) +
  scale_size_continuous(
    range =
      c(
        2.0,
        5.0
      ),
    guide =
      "none"
  ) +
  labs(
    title =
      "M4 persistent sensor residual effects",
    subtitle =
      "Posterior mean sensor effects after temporal and built-environment adjustment",
    caption =
      map_caption
  ) +
  map_theme

save_plot(
  p_resid_map,
  "report_m4_sensor_residual_effect_map.png",
  10.5,
  8.0
)

# ============================================================
# FIGURE 21 — BAYESIAN / FREQUENTIST COMPARISON
# ============================================================

bf <- read_csv(
  here(
    "outputs",
    "tables",
    "bayesian_frequentist_comparison.csv"
  ),
  show_col_types =
    FALSE
)

bf_plot <- bind_rows(
  bf %>%
    transmute(
      effect,
      method =
        "Bayesian M4",
      estimate =
        bayesian_percent_change,
      q025 =
        bayesian_q025,
      q975 =
        bayesian_q975
    ),
  bf %>%
    transmute(
      effect,
      method =
        "Frequentist NB mixed model",
      estimate =
        frequentist_percent_change,
      q025 =
        frequentist_q025,
      q975 =
        frequentist_q975
    )
) %>%
  mutate(
    effect =
      factor(
        effect,
        levels =
          rev(
            unique(
              bf$effect
            )
          )
      ),
    method =
      factor(
        method,
        levels =
          c(
            "Bayesian M4",
            "Frequentist NB mixed model"
          )
      )
  )

method_cols <- c(
  "Bayesian M4" =
    "#0072B2",
  "Frequentist NB mixed model" =
    "#D55E00"
)

p_bf <- ggplot(
  bf_plot,
  aes(
    x =
      estimate,
    y =
      effect,
    colour =
      method,
    shape =
      method
  )
) +
  geom_vline(
    xintercept =
      0,
    colour =
      "#666666",
    linewidth =
      0.75,
    linetype =
      "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin =
        q025,
      xmax =
        q975
    ),
    orientation =
      "y",
    position =
      position_dodge(
        width =
          0.52
      ),
    width =
      0.16,
    linewidth =
      0.82
  ) +
  geom_point(
    position =
      position_dodge(
        width =
          0.52
      ),
    size =
      2.7
  ) +
  scale_colour_manual(
    values =
      method_cols,
    name =
      NULL
  ) +
  scale_shape_manual(
    values =
      c(
        16,
        17
      ),
    name =
      NULL
  ) +
  labs(
    title =
      "Bayesian and frequentist effect estimates",
    subtitle =
      "Percentage change in expected daily pedestrian activity",
    x =
      "Percentage change",
    y =
      NULL
  ) +
  report_theme

save_plot(
  p_bf,
  "report_bayesian_frequentist_comparison.png",
  10.5,
  7.8
)

cat(
  "\nAll report figures were successfully generated in outputs/figures/.\n"
)
