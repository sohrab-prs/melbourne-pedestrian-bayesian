# ============================================================
# 17_posterior_predictive_checks.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Posterior predictive checking for M1, M2 and M3.
#
# No models are refitted.
# Existing posterior MCMC draws are used to generate replicated
# sensor-day datasets under each fitted model.
#
# Computational design:
#   200 posterior draws per model
#   full N = 18,671 observations retained
#
# Outputs:
#   outputs/tables/ppc_statistics.csv
#   outputs/tables/ppc_quantiles.csv
#   outputs/tables/ppc_sensor_means.csv
#   outputs/figures/ppc_statistics.png
#   outputs/figures/ppc_quantiles.png
#   outputs/figures/ppc_sensor_means.png
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
    "m1_poisson_samples.rds"
  ),
  
  here(
    "outputs",
    "mcmc",
    "m2_negbin_samples.rds"
  ),
  
  here(
    "outputs",
    "mcmc",
    "m3_hierarchical_negbin_samples.rds"
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
# 1. Modelling data
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
  sum(is.na(dat)) == 0
)

N <- nrow(dat)

J <- 53L

# ------------------------------------------------------------
# 2. Indices
# ------------------------------------------------------------

sensor_index <-
  as.integer(
    dat$sensor_index
  )

date_index <-
  as.integer(
    factor(
      dat$date,
      levels = sort(
        unique(
          dat$date
        )
      )
    )
  )

n_dates <-
  max(
    date_index
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

date_n <-
  as.numeric(
    table(
      factor(
        date_index,
        levels = seq_len(n_dates)
      )
    )
  )

# ------------------------------------------------------------
# 3. Posterior samples
# ------------------------------------------------------------

m1_samples <- readRDS(
  here(
    "outputs",
    "mcmc",
    "m1_poisson_samples.rds"
  )
)

m2_samples <- readRDS(
  here(
    "outputs",
    "mcmc",
    "m2_negbin_samples.rds"
  )
)

m3_samples <- readRDS(
  here(
    "outputs",
    "mcmc",
    "m3_hierarchical_negbin_samples.rds"
  )
)

m1_matrix <-
  as.matrix(
    m1_samples
  )

m2_matrix <-
  as.matrix(
    m2_samples
  )

m3_matrix <-
  as.matrix(
    m3_samples
  )

# ============================================================
# PART B — PPC SETTINGS
# ============================================================

n_ppc_draws <- 200L

m1_draw_ids <- sample(
  seq_len(
    nrow(
      m1_matrix
    )
  ),
  n_ppc_draws,
  replace = FALSE
)

m2_draw_ids <- sample(
  seq_len(
    nrow(
      m2_matrix
    )
  ),
  n_ppc_draws,
  replace = FALSE
)

m3_draw_ids <- sample(
  seq_len(
    nrow(
      m3_matrix
    )
  ),
  n_ppc_draws,
  replace = FALSE
)

# ============================================================
# PART C — DESIGN COMPONENTS
# ============================================================

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

dow <-
  as.integer(
    dat$dow_index
  )

# ============================================================
# PART D — HELPER FUNCTIONS
# ============================================================

get_dow_effect <- function(
    posterior_row
) {
  
  effect <- rep(
    0,
    7
  )
  
  for (k in 2:7) {
    
    parameter_name <-
      paste0(
        "beta_dow[",
        k,
        "]"
      )
    
    effect[k] <-
      posterior_row[
        parameter_name
      ]
  }
  
  effect
}

calculate_group_means <- function(
    y,
    group,
    group_sizes
) {
  
  sums <-
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
    sums
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
  
  date_means <-
    calculate_group_means(
      y,
      date_index,
      date_n
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
        date_means
      )
  )
}

# ------------------------------------------------------------
# Quantiles used for distributional PPC
# ------------------------------------------------------------

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
# PART E — OBSERVED FEATURES
# ============================================================

y_observed <-
  as.integer(
    dat$daily_ped_count
  )

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

cat(
  "\nObserved data features prepared.\n"
)

# ============================================================
# PART F — SIMULATION STORAGE
# ============================================================

model_names <- c(
  "M1 Poisson",
  "M2 Negative Binomial",
  "M3 Hierarchical Negative Binomial"
)

stat_results <- vector(
  "list",
  length(
    model_names
  ) *
    n_ppc_draws
)

quantile_results <- vector(
  "list",
  length(
    model_names
  ) *
    n_ppc_draws
)

sensor_results <- vector(
  "list",
  length(
    model_names
  ) *
    n_ppc_draws
)

result_counter <- 1L

# ============================================================
# PART G — M1 POSTERIOR PREDICTIVE
# ============================================================

cat(
  "[",
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "] Starting M1 posterior predictive simulation...\n",
  sep = ""
)

for (d in seq_len(n_ppc_draws)) {
  
  row <-
    m1_matrix[
      m1_draw_ids[d],
      ,
      drop = TRUE
    ]
  
  dow_effect <-
    get_dow_effect(
      row
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
    height
  
  mu <-
    exp(
      eta
    )
  
  y_rep <-
    rpois(
      N,
      lambda =
        mu
    )
  
  stat_results[[result_counter]] <-
    calculate_statistics(
      y_rep
    ) %>%
    mutate(
      model =
        "M1 Poisson",
      draw =
        d,
      .before = 1
    )
  
  quantile_results[[result_counter]] <-
    tibble(
      model =
        "M1 Poisson",
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
  
  sensor_results[[result_counter]] <-
    tibble(
      model =
        "M1 Poisson",
      draw =
        d,
      sensor_index =
        seq_len(J),
      replicated_sensor_mean =
        sensor_means_rep
    )
  
  result_counter <-
    result_counter +
    1L
}

cat(
  "M1 posterior predictive simulation completed.\n"
)

# ============================================================
# PART H — M2 POSTERIOR PREDICTIVE
# ============================================================

cat(
  "[",
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "] Starting M2 posterior predictive simulation...\n",
  sep = ""
)

for (d in seq_len(n_ppc_draws)) {
  
  row <-
    m2_matrix[
      m2_draw_ids[d],
      ,
      drop = TRUE
    ]
  
  dow_effect <-
    get_dow_effect(
      row
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
    height
  
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
  
  stat_results[[result_counter]] <-
    calculate_statistics(
      y_rep
    ) %>%
    mutate(
      model =
        "M2 Negative Binomial",
      draw =
        d,
      .before = 1
    )
  
  quantile_results[[result_counter]] <-
    tibble(
      model =
        "M2 Negative Binomial",
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
  
  sensor_results[[result_counter]] <-
    tibble(
      model =
        "M2 Negative Binomial",
      draw =
        d,
      sensor_index =
        seq_len(J),
      replicated_sensor_mean =
        sensor_means_rep
    )
  
  result_counter <-
    result_counter +
    1L
}

cat(
  "M2 posterior predictive simulation completed.\n"
)

# ============================================================
# PART I — M3 POSTERIOR PREDICTIVE
# ============================================================

cat(
  "[",
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "] Starting M3 posterior predictive simulation...\n",
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
  
  u_sensor_draw <- numeric(
    J
  )
  
  for (j in seq_len(J)) {
    
    u_sensor_draw[j] <-
      row[
        paste0(
          "u_sensor[",
          j,
          "]"
        )
      ]
  }
  
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
    u_sensor_draw[
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
  
  stat_results[[result_counter]] <-
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
  
  quantile_results[[result_counter]] <-
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
  
  sensor_results[[result_counter]] <-
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
  
  result_counter <-
    result_counter +
    1L
}

cat(
  "M3 posterior predictive simulation completed.\n"
)

# ============================================================
# PART J — SUMMARISE PPC STATISTICS
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

stat_draws_long <-
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
  )

ppc_statistics <-
  stat_draws_long %>%
  
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
    "ppc_statistics.csv"
  )
)

# ============================================================
# PART K — QUANTILE PPC
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
    "ppc_quantiles.csv"
  )
)

# ============================================================
# PART L — SENSOR-MEAN PPC
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
    "ppc_sensor_means.csv"
  )
)

# ============================================================
# PART M — FIGURES
# ============================================================

# ------------------------------------------------------------
# 1. Global PPC statistics
# ------------------------------------------------------------

plot_statistics <- c(
  "mean_count",
  "sd_count",
  "variance_to_mean",
  "q95_count",
  "maximum_count",
  "sd_sensor_mean",
  "sd_date_mean"
)

plot_stat_data <-
  ppc_statistics %>%
  filter(
    statistic %in%
      plot_statistics
  )

p_statistics <- ggplot(
  plot_stat_data,
  aes(
    x = model,
    y = replicated_median
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
      3,
    stroke =
      1
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
      "Posterior Predictive Checks: Observed vs Replicated Statistics",
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
        angle = 25,
        hjust = 1
      )
  )

ggsave(
  here(
    "outputs",
    "figures",
    "ppc_statistics.png"
  ),
  p_statistics,
  width =
    11,
  height =
    10,
  dpi =
    300
)

# ------------------------------------------------------------
# 2. Distribution quantiles
# ------------------------------------------------------------

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
    scales =
      "free",
    ncol =
      1
  ) +
  labs(
    title =
      "Posterior Predictive Quantile Comparison",
    subtitle =
      "Closer to the 45-degree line indicates better reproduction of the observed count distribution",
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
    "ppc_quantiles.png"
  ),
  p_quantiles,
  width =
    8,
  height =
    11,
  dpi =
    300
)

# ------------------------------------------------------------
# 3. Sensor means
# ------------------------------------------------------------

sensor_order <-
  order(
    observed_sensor_means
  )

sensor_plot_data <-
  ppc_sensor_means %>%
  
  mutate(
    sensor_order =
      match(
        sensor_index,
        sensor_order
      )
  )

p_sensor <- ggplot(
  sensor_plot_data,
  aes(
    x =
      sensor_order,
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
      1.2
  ) +
  geom_point(
    aes(
      y =
        observed_sensor_mean
    ),
    shape =
      4,
    size =
      1.7
  ) +
  facet_wrap(
    ~ model,
    ncol =
      1,
    scales =
      "free_y"
  ) +
  labs(
    title =
      "Posterior Predictive Checks for Sensor-Level Mean Activity",
    subtitle =
      "Sensors ordered by observed mean; cross = observed sensor mean",
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
    "ppc_sensor_means.png"
  ),
  p_sensor,
  width =
    10,
  height =
    11,
  dpi =
    300
)

# ============================================================
# PART N — QA SUMMARY
# ============================================================

sensor_coverage_summary <-
  ppc_sensor_means %>%
  
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

global_coverage_summary <-
  ppc_statistics %>%
  
  group_by(
    model
  ) %>%
  
  summarise(
    
    statistics_checked =
      n(),
    
    statistics_inside_95 =
      sum(
        observed_inside_95
      ),
    
    proportion_inside_95 =
      mean(
        observed_inside_95
      ),
    
    .groups =
      "drop"
  )

cat(
  "\n========== STEP 17 POSTERIOR PREDICTIVE CHECKS ==========\n"
)

cat(
  "\nObserved statistics:\n"
)

print(
  observed_statistics,
  width = Inf
)

cat(
  "\nGlobal PPC statistics:\n"
)

print(
  ppc_statistics,
  width = Inf
)

cat(
  "\nGlobal PPC coverage summary:\n"
)

print(
  global_coverage_summary,
  width = Inf
)

cat(
  "\nSensor-mean PPC coverage:\n"
)

print(
  sensor_coverage_summary,
  width = Inf
)

cat(
  "\nPosterior predictive draws per model:",
  n_ppc_draws,
  "\n"
)

cat(
  "==========================================================\n"
)

message(
  "Step 17 posterior predictive checks completed."
)