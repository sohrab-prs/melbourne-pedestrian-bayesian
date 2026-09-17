# ============================================================
# 17c_compare_m3_m4_ppc.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Targeted posterior predictive comparison of M3 and M4.
#
# M4 was motivated by the M3 PPC finding that:
#   - sensor-level heterogeneity was reproduced well;
#   - date-level variation was underestimated;
#   - overall dispersion / upper tail were overpredicted.
#
# No model is refitted.
#
# Outputs:
#   outputs/tables/ppc_m3_m4_statistics.csv
#   outputs/tables/ppc_m3_m4_quantiles.csv
#   outputs/tables/ppc_m3_m4_sensor_means.csv
#   outputs/tables/ppc_m3_m4_date_means.csv
#   outputs/figures/ppc_m3_m4_statistics.png
#   outputs/figures/ppc_m3_m4_quantiles.png
#   outputs/figures/ppc_m3_m4_sensor_means.png
#   outputs/figures/ppc_m3_m4_date_means.png
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(coda)
library(ggplot2)

set.seed(20260917)

# ============================================================
# PART A — INPUTS
# ============================================================

required_files <- c(
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  ),
  
  here(
    "outputs",
    "mcmc",
    "m3_hierarchical_negbin_samples.rds"
  ),
  
  here(
    "outputs",
    "mcmc",
    "m4_sensor_day_negbin_samples.rds"
  )
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

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------

dat <- read_csv(
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  ),
  col_types = cols(
    date = col_date(),
    .default = col_guess()
  ),
  show_col_types = FALSE
)

stopifnot(
  nrow(dat) == 18671,
  n_distinct(dat$sensor_index) == 53,
  n_distinct(dat$date) == 364,
  sum(is.na(dat)) == 0
)

N <- nrow(dat)

J <- 53L

date_lookup <- tibble(
  date =
    sort(
      unique(
        dat$date
      )
    ),
  day_index =
    seq_along(
      sort(
        unique(
          dat$date
        )
      )
    )
)

D <- nrow(
  date_lookup
)

dat <- dat %>%
  left_join(
    date_lookup,
    by = "date"
  )

stopifnot(
  D == 364,
  sum(is.na(dat$day_index)) == 0
)

# ------------------------------------------------------------
# Group indices
# ------------------------------------------------------------

sensor_index <-
  as.integer(
    dat$sensor_index
  )

day_index <-
  as.integer(
    dat$day_index
  )

sensor_n <-
  as.numeric(
    table(
      factor(
        sensor_index,
        levels = seq_len(J)
      )
    )
  )

day_n <-
  as.numeric(
    table(
      factor(
        day_index,
        levels = seq_len(D)
      )
    )
  )

# ------------------------------------------------------------
# Predictors
# ------------------------------------------------------------

dow <-
  as.integer(
    dat$dow_index
  )

holiday <-
  as.numeric(
    dat$public_holiday_num
  )

temp <-
  as.numeric(
    dat$z_mean_temp
  )

rain <-
  as.numeric(
    dat$z_log1p_precipitation
  )

business <-
  as.numeric(
    dat$z_business_count
  )

population <-
  as.numeric(
    dat$z_resident_population
  )

coverage <-
  as.numeric(
    dat$z_built_coverage
  )

height <-
  as.numeric(
    dat$z_structure_height
  )

y_observed <-
  as.integer(
    dat$daily_ped_count
  )

# ------------------------------------------------------------
# Posterior draws
# ------------------------------------------------------------

m3_matrix <-
  as.matrix(
    readRDS(
      here(
        "outputs",
        "mcmc",
        "m3_hierarchical_negbin_samples.rds"
      )
    )
  )

m4_matrix <-
  as.matrix(
    readRDS(
      here(
        "outputs",
        "mcmc",
        "m4_sensor_day_negbin_samples.rds"
      )
    )
  )

n_ppc_draws <- 200L

m3_draw_ids <- sample(
  seq_len(
    nrow(
      m3_matrix
    )
  ),
  n_ppc_draws,
  replace = FALSE
)

m4_draw_ids <- sample(
  seq_len(
    nrow(
      m4_matrix
    )
  ),
  n_ppc_draws,
  replace = FALSE
)

# ============================================================
# PART B — HELPERS
# ============================================================

get_dow_effect <- function(
    posterior_row
) {
  
  effect <- rep(
    0,
    7
  )
  
  for (k in 2:7) {
    
    effect[k] <-
      posterior_row[
        paste0(
          "beta_dow[",
          k,
          "]"
        )
      ]
  }
  
  effect
}

extract_random_effects <- function(
    posterior_row,
    prefix,
    n_effects
) {
  
  vapply(
    seq_len(
      n_effects
    ),
    function(i) {
      
      posterior_row[
        paste0(
          prefix,
          "[",
          i,
          "]"
        )
      ]
    },
    numeric(1)
  )
}

calculate_group_means <- function(
    y,
    group,
    group_sizes
) {
  
  group_sums <-
    rowsum(
      matrix(
        y,
        ncol = 1
      ),
      group =
        group,
      reorder =
        TRUE
    )[
      ,
      1
    ]
  
  as.numeric(
    group_sums
  ) /
    group_sizes
}

calculate_statistics <- function(
    y
) {
  
  sensor_means <-
    calculate_group_means(
      y,
      sensor_index,
      sensor_n
    )
  
  day_means <-
    calculate_group_means(
      y,
      day_index,
      day_n
    )
  
  tibble(
    
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
      unname(
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
        y > 50000
      ),
    
    sd_sensor_mean =
      sd(
        sensor_means
      ),
    
    maximum_sensor_mean =
      max(
        sensor_means
      ),
    
    sd_date_mean =
      sd(
        day_means
      )
  )
}

quantile_probabilities <- c(
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

calculate_quantiles <- function(
    y
) {
  
  as.numeric(
    quantile(
      y,
      probs =
        quantile_probabilities,
      names =
        FALSE
    )
  )
}

# ============================================================
# PART C — OBSERVED FEATURES
# ============================================================

observed_statistics <-
  calculate_statistics(
    y_observed
  )

observed_quantiles <-
  calculate_quantiles(
    y_observed
  )

observed_sensor_means <-
  calculate_group_means(
    y_observed,
    sensor_index,
    sensor_n
  )

observed_day_means <-
  calculate_group_means(
    y_observed,
    day_index,
    day_n
  )

# ============================================================
# PART D — STORAGE
# ============================================================

stat_results <- list()

quantile_results <- list()

sensor_results <- list()

day_results <- list()

counter <- 1L

# ============================================================
# PART E — M3 PPC
# ============================================================

cat(
  "\n[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting targeted M3 PPC...\n",
  sep = ""
)

for (d in seq_len(n_ppc_draws)) {
  
  row <-
    m3_matrix[
      m3_draw_ids[d],
      ,
      drop = TRUE
    ]
  
  dow_effect <-
    get_dow_effect(
      row
    )
  
  u_sensor <-
    extract_random_effects(
      row,
      "u_sensor",
      J
    )
  
  eta <-
    row["alpha"] +
    dow_effect[dow] +
    row["beta_holiday"] *
    holiday +
    row["beta_temp"] *
    temp +
    row["beta_rain"] *
    rain +
    row["beta_business"] *
    business +
    row["beta_population"] *
    population +
    row["beta_coverage"] *
    coverage +
    row["beta_height"] *
    height +
    u_sensor[
      sensor_index
    ]
  
  mu <-
    exp(
      eta
    )
  
  y_rep <-
    rnbinom(
      N,
      size =
        row["phi"],
      mu =
        mu
    )
  
  stat_results[[counter]] <-
    calculate_statistics(
      y_rep
    ) %>%
    mutate(
      model =
        "M3 Hierarchical Negative Binomial",
      draw =
        d,
      .before = 1
    )
  
  quantile_results[[counter]] <-
    tibble(
      model =
        "M3 Hierarchical Negative Binomial",
      draw =
        d,
      probability =
        quantile_probabilities,
      replicated_quantile =
        calculate_quantiles(
          y_rep
        )
    )
  
  sensor_means_rep <-
    calculate_group_means(
      y_rep,
      sensor_index,
      sensor_n
    )
  
  day_means_rep <-
    calculate_group_means(
      y_rep,
      day_index,
      day_n
    )
  
  sensor_results[[counter]] <-
    tibble(
      model =
        "M3 Hierarchical Negative Binomial",
      draw =
        d,
      sensor_index =
        seq_len(J),
      replicated_sensor_mean =
        sensor_means_rep
    )
  
  day_results[[counter]] <-
    tibble(
      model =
        "M3 Hierarchical Negative Binomial",
      draw =
        d,
      day_index =
        seq_len(D),
      replicated_day_mean =
        day_means_rep
    )
  
  counter <-
    counter +
    1L
}

cat(
  "M3 targeted PPC completed.\n"
)

# ============================================================
# PART F — M4 PPC
# ============================================================

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting targeted M4 PPC...\n",
  sep = ""
)

for (d in seq_len(n_ppc_draws)) {
  
  row <-
    m4_matrix[
      m4_draw_ids[d],
      ,
      drop = TRUE
    ]
  
  dow_effect <-
    get_dow_effect(
      row
    )
  
  u_sensor <-
    extract_random_effects(
      row,
      "u_sensor",
      J
    )
  
  day_effect <-
    extract_random_effects(
      row,
      "day_effect",
      D
    )
  
  eta <-
    row["alpha"] +
    dow_effect[dow] +
    row["beta_holiday"] *
    holiday +
    row["beta_temp"] *
    temp +
    row["beta_rain"] *
    rain +
    row["beta_business"] *
    business +
    row["beta_population"] *
    population +
    row["beta_coverage"] *
    coverage +
    row["beta_height"] *
    height +
    u_sensor[
      sensor_index
    ] +
    day_effect[
      day_index
    ]
  
  mu <-
    exp(
      eta
    )
  
  y_rep <-
    rnbinom(
      N,
      size =
        row["phi"],
      mu =
        mu
    )
  
  stat_results[[counter]] <-
    calculate_statistics(
      y_rep
    ) %>%
    mutate(
      model =
        "M4 Sensor-Day Hierarchical Negative Binomial",
      draw =
        d,
      .before = 1
    )
  
  quantile_results[[counter]] <-
    tibble(
      model =
        "M4 Sensor-Day Hierarchical Negative Binomial",
      draw =
        d,
      probability =
        quantile_probabilities,
      replicated_quantile =
        calculate_quantiles(
          y_rep
        )
    )
  
  sensor_means_rep <-
    calculate_group_means(
      y_rep,
      sensor_index,
      sensor_n
    )
  
  day_means_rep <-
    calculate_group_means(
      y_rep,
      day_index,
      day_n
    )
  
  sensor_results[[counter]] <-
    tibble(
      model =
        "M4 Sensor-Day Hierarchical Negative Binomial",
      draw =
        d,
      sensor_index =
        seq_len(J),
      replicated_sensor_mean =
        sensor_means_rep
    )
  
  day_results[[counter]] <-
    tibble(
      model =
        "M4 Sensor-Day Hierarchical Negative Binomial",
      draw =
        d,
      day_index =
        seq_len(D),
      replicated_day_mean =
        day_means_rep
    )
  
  counter <-
    counter +
    1L
}

cat(
  "M4 targeted PPC completed.\n"
)

# ============================================================
# PART G — GLOBAL STATISTICS
# ============================================================

stat_draws <-
  bind_rows(
    stat_results
  )

observed_statistics_long <-
  observed_statistics %>%
  pivot_longer(
    cols =
      everything(),
    names_to =
      "statistic",
    values_to =
      "observed"
  )

ppc_statistics <-
  stat_draws %>%
  
  pivot_longer(
    cols =
      -c(
        model,
        draw
      ),
    names_to =
      "statistic",
    values_to =
      "replicated"
  ) %>%
  
  left_join(
    observed_statistics_long,
    by =
      "statistic"
  ) %>%
  
  group_by(
    model,
    statistic
  ) %>%
  
  summarise(
    
    observed =
      first(
        observed
      ),
    
    replicated_q025 =
      unname(
        quantile(
          replicated,
          0.025
        )
      ),
    
    replicated_median =
      median(
        replicated
      ),
    
    replicated_q975 =
      unname(
        quantile(
          replicated,
          0.975
        )
      ),
    
    posterior_predictive_p_upper =
      mean(
        replicated >=
          first(
            observed
          )
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    observed_inside_95 =
      observed >=
      replicated_q025 &
      observed <=
      replicated_q975
  )

write_csv(
  ppc_statistics,
  here(
    "outputs",
    "tables",
    "ppc_m3_m4_statistics.csv"
  )
)

# ============================================================
# PART H — QUANTILES
# ============================================================

quantile_draws <-
  bind_rows(
    quantile_results
  )

observed_quantile_table <-
  tibble(
    probability =
      quantile_probabilities,
    observed_quantile =
      observed_quantiles
  )

ppc_quantiles <-
  quantile_draws %>%
  
  group_by(
    model,
    probability
  ) %>%
  
  summarise(
    
    replicated_q025 =
      unname(
        quantile(
          replicated_quantile,
          0.025
        )
      ),
    
    replicated_median =
      median(
        replicated_quantile
      ),
    
    replicated_q975 =
      unname(
        quantile(
          replicated_quantile,
          0.975
        )
      ),
    
    .groups =
      "drop"
  ) %>%
  
  left_join(
    observed_quantile_table,
    by =
      "probability"
  )

write_csv(
  ppc_quantiles,
  here(
    "outputs",
    "tables",
    "ppc_m3_m4_quantiles.csv"
  )
)

# ============================================================
# PART I — SENSOR-MEAN PPC
# ============================================================

sensor_draws <-
  bind_rows(
    sensor_results
  )

ppc_sensor_means <-
  sensor_draws %>%
  
  group_by(
    model,
    sensor_index
  ) %>%
  
  summarise(
    
    replicated_q025 =
      unname(
        quantile(
          replicated_sensor_mean,
          0.025
        )
      ),
    
    replicated_median =
      median(
        replicated_sensor_mean
      ),
    
    replicated_q975 =
      unname(
        quantile(
          replicated_sensor_mean,
          0.975
        )
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    
    observed_sensor_mean =
      observed_sensor_means[
        sensor_index
      ],
    
    observed_inside_95 =
      observed_sensor_mean >=
      replicated_q025 &
      observed_sensor_mean <=
      replicated_q975
  )

write_csv(
  ppc_sensor_means,
  here(
    "outputs",
    "tables",
    "ppc_m3_m4_sensor_means.csv"
  )
)

# ============================================================
# PART J — DATE-MEAN PPC
# ============================================================

day_draws <-
  bind_rows(
    day_results
  )

ppc_day_means <-
  day_draws %>%
  
  group_by(
    model,
    day_index
  ) %>%
  
  summarise(
    
    replicated_q025 =
      unname(
        quantile(
          replicated_day_mean,
          0.025
        )
      ),
    
    replicated_median =
      median(
        replicated_day_mean
      ),
    
    replicated_q975 =
      unname(
        quantile(
          replicated_day_mean,
          0.975
        )
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    
    date =
      date_lookup$date[
        day_index
      ],
    
    observed_day_mean =
      observed_day_means[
        day_index
      ],
    
    observed_inside_95 =
      observed_day_mean >=
      replicated_q025 &
      observed_day_mean <=
      replicated_q975
  )

write_csv(
  ppc_day_means,
  here(
    "outputs",
    "tables",
    "ppc_m3_m4_date_means.csv"
  )
)

# ============================================================
# PART K — FIGURES
# ============================================================

plot_statistics <- c(
  "mean_count",
  "median_count",
  "sd_count",
  "variance_to_mean",
  "q95_count",
  "maximum_count",
  "proportion_gt_50000",
  "sd_sensor_mean",
  "maximum_sensor_mean",
  "sd_date_mean"
)

p_statistics <- ppc_statistics %>%
  filter(
    statistic %in%
      plot_statistics
  ) %>%
  ggplot(
    aes(
      x =
        model,
      y =
        replicated_median
    )
  ) +
  geom_errorbar(
    aes(
      ymin =
        replicated_q025,
      ymax =
        replicated_q975
    ),
    width =
      0.15
  ) +
  geom_point(
    size =
      2
  ) +
  geom_point(
    aes(
      y =
        observed
    ),
    shape =
      4,
    size =
      3
  ) +
  facet_wrap(
    ~ statistic,
    scales =
      "free_y",
    ncol =
      2
  ) +
  labs(
    title =
      "M3 vs M4 Posterior Predictive Checks",
    subtitle =
      "Cross = observed; point and interval = posterior predictive median and 95% interval",
    x =
      NULL,
    y =
      NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x =
      element_text(
        angle = 20,
        hjust = 1
      )
  )

ggsave(
  here(
    "outputs",
    "figures",
    "ppc_m3_m4_statistics.png"
  ),
  p_statistics,
  width = 11,
  height = 11,
  dpi = 300
)

p_quantiles <- ggplot(
  ppc_quantiles,
  aes(
    x =
      observed_quantile,
    y =
      replicated_median
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin =
        replicated_q025,
      ymax =
        replicated_q975
    ),
    width =
      0
  ) +
  geom_point(
    size =
      2
  ) +
  facet_wrap(
    ~ model,
    ncol =
      1,
    scales =
      "free"
  ) +
  labs(
    title =
      "M3 vs M4 Posterior Predictive Quantiles",
    x =
      "Observed count quantile",
    y =
      "Posterior predictive count quantile"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "ppc_m3_m4_quantiles.png"
  ),
  p_quantiles,
  width = 8,
  height = 8,
  dpi = 300
)

sensor_order <-
  order(
    observed_sensor_means
  )

p_sensor <- ppc_sensor_means %>%
  mutate(
    sensor_rank =
      match(
        sensor_index,
        sensor_order
      )
  ) %>%
  ggplot(
    aes(
      x =
        sensor_rank,
      y =
        replicated_median
    )
  ) +
  geom_linerange(
    aes(
      ymin =
        replicated_q025,
      ymax =
        replicated_q975
    )
  ) +
  geom_point(
    size =
      1
  ) +
  geom_point(
    aes(
      y =
        observed_sensor_mean
    ),
    shape =
      4,
    size =
      1.5
  ) +
  facet_wrap(
    ~ model,
    ncol =
      1
  ) +
  labs(
    title =
      "M3 vs M4 Sensor-Mean Posterior Predictive Checks",
    x =
      "Sensor rank by observed mean",
    y =
      "Daily pedestrian count"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "ppc_m3_m4_sensor_means.png"
  ),
  p_sensor,
  width = 10,
  height = 8,
  dpi = 300
)

p_day <- ppc_day_means %>%
  ggplot(
    aes(
      x =
        date,
      y =
        replicated_median
    )
  ) +
  geom_ribbon(
    aes(
      ymin =
        replicated_q025,
      ymax =
        replicated_q975
    ),
    alpha =
      0.25
  ) +
  geom_line() +
  geom_point(
    aes(
      y =
        observed_day_mean
    ),
    shape =
      4,
    size =
      0.8
  ) +
  facet_wrap(
    ~ model,
    ncol =
      1
  ) +
  labs(
    title =
      "M3 vs M4 Date-Mean Posterior Predictive Checks",
    subtitle =
      "Cross = observed daily mean across available sensors",
    x =
      "Date",
    y =
      "Mean pedestrian count"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "ppc_m3_m4_date_means.png"
  ),
  p_day,
  width = 11,
  height = 8,
  dpi = 300
)

# ============================================================
# PART L — QA
# ============================================================

global_coverage <- ppc_statistics %>%
  group_by(
    model
  ) %>%
  summarise(
    statistics_inside_95 =
      sum(
        observed_inside_95
      ),
    
    statistics_checked =
      n(),
    
    proportion_inside_95 =
      mean(
        observed_inside_95
      ),
    
    .groups =
      "drop"
  )

sensor_coverage <- ppc_sensor_means %>%
  group_by(
    model
  ) %>%
  summarise(
    sensor_mean_95_coverage =
      mean(
        observed_inside_95
      ),
    .groups =
      "drop"
  )

day_coverage <- ppc_day_means %>%
  group_by(
    model
  ) %>%
  summarise(
    date_mean_95_coverage =
      mean(
        observed_inside_95
      ),
    .groups =
      "drop"
  )

cat(
  "\n========== TARGETED M3 VS M4 PPC ==========\n"
)

cat(
  "\nGlobal PPC statistics:\n"
)

print(
  ppc_statistics,
  n = Inf,
  width = Inf
)

cat(
  "\nGlobal coverage:\n"
)

print(
  global_coverage,
  width = Inf
)

cat(
  "\nSensor-mean coverage:\n"
)

print(
  sensor_coverage,
  width = Inf
)

cat(
  "\nDate-mean coverage:\n"
)

print(
  day_coverage,
  width = Inf
)

cat(
  "\nPosterior predictive draws per model:",
  n_ppc_draws,
  "\n"
)

cat(
  "===========================================\n"
)

message(
  "Targeted M3 vs M4 posterior predictive comparison completed."
)