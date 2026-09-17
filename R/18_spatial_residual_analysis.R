# ============================================================
# 18_spatial_residual_analysis.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Assess residual spatial dependence in the accepted M4 model.
#
# Target:
# Persistent sensor-specific residual effects u_sensor[j].
#
# Spatial diagnostics:
#   - Moran's I using symmetrised k-nearest-neighbour weights
#   - k = 4 primary
#   - k = 6 sensitivity
#   - 999 permutation tests on posterior mean sensor effects
#   - posterior distribution of Moran's I using 5,000 MCMC draws
#
# No model is refitted.
#
# Outputs:
#   outputs/tables/m4_spatial_residual_moran.csv
#   outputs/tables/m4_spatial_residual_sensor_effects.csv
#   outputs/figures/m4_sensor_residual_effect_map.png
#   outputs/figures/m4_residual_moran_posterior.png
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(sf)
library(coda)
library(ggplot2)

set.seed(20260917)

# ============================================================
# PART A — INPUTS
# ============================================================

sensor_file <- here(
  "data",
  "interim",
  "pedestrian_sensors_2025.gpkg"
)

effect_file <- here(
  "outputs",
  "tables",
  "m4_sensor_day_negbin_sensor_effects.csv"
)

sample_file <- here(
  "outputs",
  "mcmc",
  "m4_sensor_day_negbin_samples.rds"
)

required_files <- c(
  sensor_file,
  effect_file,
  sample_file
)

stopifnot(
  all(
    file.exists(
      required_files
    )
  )
)

dir.create(
  here(
    "outputs",
    "tables"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here(
    "outputs",
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PART B — READ SENSOR LOCATIONS
# ============================================================

sensor_sf <- st_read(
  sensor_file,
  quiet = TRUE
)

sensor_effects <- read_csv(
  effect_file,
  show_col_types = FALSE
)

stopifnot(
  nrow(sensor_effects) == 53,
  "location_id" %in% names(sensor_effects)
)

# ------------------------------------------------------------
# Find location identifier in spatial file
# ------------------------------------------------------------

possible_id_names <- c(
  "location_id",
  "Location_ID",
  "LOCATION_ID",
  "sensor_id",
  "Sensor_ID"
)

spatial_id_name <-
  possible_id_names[
    possible_id_names %in%
      names(sensor_sf)
  ][1]

if (
  is.na(
    spatial_id_name
  )
) {
  stop(
    paste(
      "Could not identify a location ID column in",
      sensor_file,
      "\nAvailable columns:",
      paste(
        names(sensor_sf),
        collapse = ", "
      )
    )
  )
}

sensor_sf <- sensor_sf %>%
  mutate(
    location_id_join =
      as.character(
        .data[[spatial_id_name]]
      )
  )

sensor_effects <- sensor_effects %>%
  mutate(
    location_id_join =
      as.character(
        location_id
      )
  )

# ------------------------------------------------------------
# Keep exactly the 53 modelling sensors.
# Start with the sf object so geometry is retained safely.
# ------------------------------------------------------------

residual_sf <- sensor_sf %>%
  select(
    location_id_join
  ) %>%
  inner_join(
    sensor_effects,
    by = "location_id_join"
  )

stopifnot(
  inherits(
    residual_sf,
    "sf"
  ),
  nrow(
    residual_sf
  ) == 53,
  sum(
    st_is_empty(
      residual_sf
    )
  ) == 0,
  all(
    !is.na(
      residual_sf$posterior_mean
    )
  )
)

# Ensure one geometry per modelling sensor.
if (
  anyDuplicated(
    residual_sf$sensor_index
  ) > 0
) {
  stop(
    "Duplicate modelling sensor geometries detected."
  )
}

residual_sf <- residual_sf %>%
  mutate(
    sensor_index =
      as.integer(
        sensor_index
      )
  ) %>%
  arrange(
    sensor_index
  )

stopifnot(
  nrow(residual_sf) == 53,
  anyDuplicated(residual_sf$sensor_index) == 0,
  all(
    residual_sf$sensor_index ==
      seq_len(53)
  )
)

# ============================================================
# PART C — PROJECTED COORDINATES
# ============================================================

# GDA2020 / MGA Zone 55.
# Appropriate for Melbourne metric distance calculations.

residual_sf_projected <- st_transform(
  residual_sf,
  7855
)

coords <- st_coordinates(
  residual_sf_projected
)

stopifnot(
  nrow(coords) == 53
)

# Pairwise Euclidean distances in metres.
distance_matrix <- as.matrix(
  dist(
    coords[
      ,
      c(
        "X",
        "Y"
      )
    ]
  )
)

diag(
  distance_matrix
) <- Inf

# ============================================================
# PART D — SPATIAL WEIGHT MATRICES
# ============================================================

make_knn_weights <- function(
    distance_matrix,
    k
) {
  
  n <- nrow(
    distance_matrix
  )
  
  adjacency <- matrix(
    0,
    nrow = n,
    ncol = n
  )
  
  for (i in seq_len(n)) {
    
    neighbours <-
      order(
        distance_matrix[
          i,
        ]
      )[
        seq_len(k)
      ]
    
    adjacency[
      i,
      neighbours
    ] <- 1
  }
  
  # Symmetrise:
  # sensors are neighbours if either sensor identifies the
  # other among its k nearest neighbours.
  adjacency <-
    1 * (
      adjacency +
        t(
          adjacency
        ) >
        0
    )
  
  diag(
    adjacency
  ) <- 0
  
  row_totals <-
    rowSums(
      adjacency
    )
  
  stopifnot(
    all(
      row_totals > 0
    )
  )
  
  # Row-standardised spatial weights.
  W <-
    adjacency /
    row_totals
  
  W
}

W4 <- make_knn_weights(
  distance_matrix,
  k = 4
)

W6 <- make_knn_weights(
  distance_matrix,
  k = 6
)

# ============================================================
# PART E — MORAN'S I FUNCTIONS
# ============================================================

morans_i <- function(
    x,
    W
) {
  
  x <- as.numeric(
    x
  )
  
  n <- length(
    x
  )
  
  z <-
    x -
    mean(
      x
    )
  
  denominator <-
    sum(
      z^2
    )
  
  if (
    denominator == 0
  ) {
    return(
      NA_real_
    )
  }
  
  S0 <-
    sum(
      W
    )
  
  numerator <-
    as.numeric(
      t(z) %*%
        W %*%
        z
    )
  
  (
    n /
      S0
  ) *
    (
      numerator /
        denominator
    )
}

permutation_moran_test <- function(
    x,
    W,
    permutations = 999
) {
  
  observed_I <-
    morans_i(
      x,
      W
    )
  
  n <- length(
    x
  )
  
  expected_I <-
    -1 /
    (
      n - 1
    )
  
  permutation_I <- replicate(
    permutations,
    morans_i(
      sample(
        x,
        replace = FALSE
      ),
      W
    )
  )
  
  # Positive spatial autocorrelation test.
  p_positive <-
    (
      1 +
        sum(
          permutation_I >=
            observed_I
        )
    ) /
    (
      permutations +
        1
    )
  
  # Two-sided permutation test centred on the randomisation
  # expectation E[I] = -1/(n-1).
  p_two_sided <-
    (
      1 +
        sum(
          abs(
            permutation_I -
              expected_I
          ) >=
            abs(
              observed_I -
                expected_I
            )
        )
    ) /
    (
      permutations +
        1
    )
  
  list(
    observed_I =
      observed_I,
    
    expected_I =
      expected_I,
    
    permutation_mean =
      mean(
        permutation_I
      ),
    
    permutation_q025 =
      unname(
        quantile(
          permutation_I,
          0.025
        )
      ),
    
    permutation_q975 =
      unname(
        quantile(
          permutation_I,
          0.975
        )
      ),
    
    p_positive =
      p_positive,
    
    p_two_sided =
      p_two_sided
  )
}

# ============================================================
# PART F — CLASSICAL TEST ON POSTERIOR MEAN EFFECTS
# ============================================================

u_mean <-
  residual_sf_projected$posterior_mean

moran4_test <- permutation_moran_test(
  u_mean,
  W4,
  permutations = 999
)

moran6_test <- permutation_moran_test(
  u_mean,
  W6,
  permutations = 999
)

# ============================================================
# PART G — POSTERIOR DISTRIBUTION OF MORAN'S I
# ============================================================

m4_samples <- readRDS(
  sample_file
)

stopifnot(
  is.mcmc.list(
    m4_samples
  )
)

posterior_matrix <- as.matrix(
  m4_samples
)

u_names <- paste0(
  "u_sensor[",
  seq_len(53),
  "]"
)

stopifnot(
  all(
    u_names %in%
      colnames(
        posterior_matrix
      )
  )
)

u_matrix <-
  posterior_matrix[
    ,
    u_names,
    drop = FALSE
  ]

# ------------------------------------------------------------
# Use 5,000 posterior draws for spatial diagnostic.
# This is more than sufficient for a stable summary while
# avoiding unnecessary computation.
# ------------------------------------------------------------

n_spatial_draws <-
  min(
    5000L,
    nrow(
      u_matrix
    )
  )

spatial_draw_ids <- sample(
  seq_len(
    nrow(
      u_matrix
    )
  ),
  n_spatial_draws,
  replace = FALSE
)

u_sub <-
  u_matrix[
    spatial_draw_ids,
    ,
    drop = FALSE
  ]

posterior_moran <- function(
    X,
    W
) {
  
  # Centre each posterior draw across sensors.
  X_centered <-
    X -
    rowMeans(
      X
    )
  
  numerator <-
    rowSums(
      (
        X_centered %*%
          W
      ) *
        X_centered
    )
  
  denominator <-
    rowSums(
      X_centered^2
    )
  
  n <- ncol(
    X_centered
  )
  
  S0 <- sum(
    W
  )
  
  (
    n /
      S0
  ) *
    numerator /
    denominator
}

moran4_posterior <- posterior_moran(
  u_sub,
  W4
)

moran6_posterior <- posterior_moran(
  u_sub,
  W6
)

# ============================================================
# PART H — SUMMARISE SPATIAL DIAGNOSTICS
# ============================================================

summarise_moran <- function(
    k,
    test_result,
    posterior_values
) {
  
  tibble(
    
    neighbours_k =
      k,
    
    observed_moran_I =
      test_result$observed_I,
    
    randomisation_expected_I =
      test_result$expected_I,
    
    permutation_mean =
      test_result$permutation_mean,
    
    permutation_q025 =
      test_result$permutation_q025,
    
    permutation_q975 =
      test_result$permutation_q975,
    
    permutation_p_positive =
      test_result$p_positive,
    
    permutation_p_two_sided =
      test_result$p_two_sided,
    
    posterior_moran_mean =
      mean(
        posterior_values
      ),
    
    posterior_moran_median =
      median(
        posterior_values
      ),
    
    posterior_moran_q025 =
      unname(
        quantile(
          posterior_values,
          0.025
        )
      ),
    
    posterior_moran_q975 =
      unname(
        quantile(
          posterior_values,
          0.975
        )
      ),
    
    posterior_probability_I_gt_0 =
      mean(
        posterior_values > 0
      )
  )
}

moran_summary <- bind_rows(
  
  summarise_moran(
    4,
    moran4_test,
    moran4_posterior
  ),
  
  summarise_moran(
    6,
    moran6_test,
    moran6_posterior
  )
)

write_csv(
  moran_summary,
  here(
    "outputs",
    "tables",
    "m4_spatial_residual_moran.csv"
  )
)

# ============================================================
# PART I — SENSOR TABLE
# ============================================================

sensor_output <- residual_sf_projected %>%
  mutate(
    
    x_m =
      coords[
        ,
        "X"
      ],
    
    y_m =
      coords[
        ,
        "Y"
      ]
  ) %>%
  st_drop_geometry() %>%
  select(
    sensor_index,
    location_id,
    posterior_mean,
    q025,
    median,
    q975,
    x_m,
    y_m
  )

write_csv(
  sensor_output,
  here(
    "outputs",
    "tables",
    "m4_spatial_residual_sensor_effects.csv"
  )
)

# ============================================================
# PART J — FIGURES
# ============================================================

# ------------------------------------------------------------
# 1. Map of posterior mean residual sensor effect
# ------------------------------------------------------------

p_map <- ggplot(
  residual_sf_projected
) +
  geom_sf(
    aes(
      colour =
        posterior_mean,
      size =
        abs(
          posterior_mean
        )
    ),
    alpha =
      0.85
  ) +
  scale_colour_gradient2(
    midpoint =
      0
  ) +
  labs(
    title =
      "M4 Persistent Sensor Residual Effects",
    subtitle =
      "Posterior mean sensor random effects after temporal and built-environment adjustment",
    colour =
      "Posterior mean",
    size =
      "|Effect|"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_sensor_residual_effect_map.png"
  ),
  p_map,
  width = 8,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 2. Posterior Moran's I distributions
# ------------------------------------------------------------

posterior_moran_plot_data <- bind_rows(
  
  tibble(
    neighbours =
      "4-nearest neighbours",
    moran_I =
      moran4_posterior
  ),
  
  tibble(
    neighbours =
      "6-nearest neighbours",
    moran_I =
      moran6_posterior
  )
)

p_moran <- ggplot(
  posterior_moran_plot_data,
  aes(
    x =
      moran_I
  )
) +
  geom_histogram(
    bins = 45
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ neighbours,
    ncol = 1
  ) +
  labs(
    title =
      "Posterior Distribution of Moran's I for M4 Sensor Effects",
    subtitle =
      "Positive values indicate residual spatial clustering",
    x =
      "Moran's I",
    y =
      "Posterior draws"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_residual_moran_posterior.png"
  ),
  p_moran,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# PART K — QA
# ============================================================

cat(
  "\n========== STEP 18 SPATIAL RESIDUAL ANALYSIS ==========\n"
)

cat(
  "\nSensors analysed:",
  nrow(
    residual_sf_projected
  ),
  "\n"
)

cat(
  "Spatial CRS:",
  st_crs(
    residual_sf_projected
  )$input,
  "\n"
)

cat(
  "Posterior draws used:",
  n_spatial_draws,
  "\n"
)

cat(
  "\nMoran's I results:\n"
)

print(
  moran_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nInterpretation flags:\n"
)

for (
  i in seq_len(
    nrow(
      moran_summary
    )
  )
) {
  
  cat(
    "  k =",
    moran_summary$neighbours_k[i],
    ": permutation p(+)=",
    round(
      moran_summary$permutation_p_positive[i],
      4
    ),
    ", posterior P(I>0)=",
    round(
      moran_summary$posterior_probability_I_gt_0[i],
      4
    ),
    "\n"
  )
}

cat(
  "=========================================================\n"
)

message(
  "Step 18 spatial residual analysis completed."
)