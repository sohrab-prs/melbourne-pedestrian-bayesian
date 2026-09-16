# ============================================================
# 10_prior_predictive.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Perform prior predictive checks for the three planned
# Bayesian count models before fitting the observed response.
#
# Models:
#   M1 Poisson
#   M2 Negative Binomial
#   M3 Hierarchical Negative Binomial
#
# No observed pedestrian counts are used to generate prior
# predictive values. Only the observed predictor structure
# (days, sensors, covariates) is used.
#
# Outputs:
#   outputs/tables/prior_predictive_summary.csv
#   outputs/tables/prior_parameter_summary.csv
#   outputs/figures/prior_predictive_*.png
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

set.seed(20260916)

# ------------------------------------------------------------
# 1. Output directories
# ------------------------------------------------------------

dir.create(
  here("outputs", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here("outputs", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load selected modelling structure
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
  n_distinct(dat$location_id) == 53,
  sum(is.na(dat)) == 0
)

N <- nrow(dat)
J <- n_distinct(dat$sensor_index)

# ------------------------------------------------------------
# 3. Prior settings
# ------------------------------------------------------------

n_prior_draws <- 300

# Baseline daily count:
# median prior expectation = approximately 10,000 pedestrians/day.
#
# SD = 0.5 still allows substantial baseline uncertainty while
# avoiding implausibly extreme average sensor volumes.
alpha_mean <- log(10000)
alpha_sd   <- 0.5

# Standardised continuous predictors and categorical contrasts.
#
# 95% prior interval is approximately +/- 0.49 on the log scale,
# corresponding to a multiplicative effect of roughly
# 0.61 to 1.63 for a one-SD predictor change.
beta_sd <- 0.25

# Negative Binomial dispersion.
#
# Median phi = 10, with substantial but more controlled
# uncertainty than the original prior.
log_phi_mean <- log(10)
log_phi_sd   <- 1.0

# Sensor-level heterogeneity.
#
# Half-Normal(0, 0.5) permits meaningful unexplained
# between-sensor variability without routinely generating
# extreme sensor effects.
sigma_sensor_prior_sd <- 0.5

# ------------------------------------------------------------
# 4. Draw parameters from priors
# ------------------------------------------------------------

alpha_draws <- rnorm(
  n_prior_draws,
  mean = alpha_mean,
  sd = alpha_sd
)

# Monday is reference = 0.
# Six coefficients represent Tuesday-Sunday.
beta_dow <- matrix(
  rnorm(
    n_prior_draws * 6,
    mean = 0,
    sd = beta_sd
  ),
  nrow = n_prior_draws,
  ncol = 6
)

colnames(beta_dow) <- c(
  "Tuesday",
  "Wednesday",
  "Thursday",
  "Friday",
  "Saturday",
  "Sunday"
)

beta_holiday <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_temp <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_rain <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_business <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_population <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_coverage <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

beta_height <- rnorm(
  n_prior_draws,
  0,
  beta_sd
)

log_phi <- rnorm(
  n_prior_draws,
  mean = log_phi_mean,
  sd = log_phi_sd
)

phi_draws <- exp(
  log_phi
)

# Half-normal prior
sigma_sensor <- abs(
  rnorm(
    n_prior_draws,
    mean = 0,
    sd = sigma_sensor_prior_sd
  )
)

# ------------------------------------------------------------
# 5. Summarise parameter priors
# ------------------------------------------------------------

prior_parameter_summary <- tibble(
  parameter = c(
    "alpha",
    "beta_generic",
    "phi",
    "sigma_sensor"
  ),
  
  prior = c(
    "Normal(log(10000), 0.5^2)",
    "Normal(0, 0.25^2)",
    "log(phi) ~ Normal(log(10), 1^2)",
    "Half-Normal(0, 0.5)"
  ),
  
  q025 = c(
    quantile(alpha_draws, 0.025),
    qnorm(0.025, 0, beta_sd),
    quantile(phi_draws, 0.025),
    quantile(sigma_sensor, 0.025)
  ),
  
  median = c(
    median(alpha_draws),
    0,
    median(phi_draws),
    median(sigma_sensor)
  ),
  
  q975 = c(
    quantile(alpha_draws, 0.975),
    qnorm(0.975, 0, beta_sd),
    quantile(phi_draws, 0.975),
    quantile(sigma_sensor, 0.975)
  )
)

write_csv(
  prior_parameter_summary,
  here(
    "outputs",
    "tables",
    "prior_parameter_summary.csv"
  )
)

# ============================================================
# PART A — PREPARE DESIGN INFORMATION
# ============================================================

# ------------------------------------------------------------
# 6. Day-of-week coefficient lookup
# ------------------------------------------------------------

# dow_index:
# 1 Monday
# 2 Tuesday
# ...
# 7 Sunday
#
# Prior coefficient for Monday is fixed at zero.

dow_design <- model.matrix(
  ~ day_of_week,
  data = dat
)

# Drop intercept and Monday reference.
dow_design <- dow_design[
  ,
  -1,
  drop = FALSE
]

stopifnot(
  ncol(dow_design) == 6
)

# ------------------------------------------------------------
# 7. Continuous/categorical predictor matrix
# ------------------------------------------------------------

X <- as.matrix(
  dat %>%
    select(
      public_holiday_num,
      z_mean_temp,
      z_log1p_precipitation,
      z_business_count,
      z_resident_population,
      z_built_coverage,
      z_structure_height
    )
)

stopifnot(
  ncol(X) == 7
)

# ============================================================
# PART B — PRIOR PREDICTIVE SIMULATION
# ============================================================

# ------------------------------------------------------------
# 8. Storage
#
# We keep only summary statistics and a manageable random
# sample of predicted counts rather than saving every
# simulated count from every model.
# ------------------------------------------------------------

summary_list <- vector(
  "list",
  n_prior_draws * 3
)

sample_list <- vector(
  "list",
  n_prior_draws * 3
)

summary_counter <- 1
sample_counter <- 1

# Number of predicted observations retained from each draw
n_sample_per_draw <- 200

sample_indices <- sample(
  seq_len(N),
  size = n_sample_per_draw,
  replace = FALSE
)

# ------------------------------------------------------------
# 9. Loop over prior draws
# ------------------------------------------------------------

for (s in seq_len(n_prior_draws)) {
  
  # ----------------------------------------------------------
  # Common non-hierarchical linear predictor
  # ----------------------------------------------------------
  
  dow_effect <-
    as.numeric(
      dow_design %*%
        beta_dow[
          s,
          ,
          drop = TRUE
        ]
    )
  
  eta_common <-
    alpha_draws[s] +
    dow_effect +
    beta_holiday[s] *
    X[, "public_holiday_num"] +
    beta_temp[s] *
    X[, "z_mean_temp"] +
    beta_rain[s] *
    X[, "z_log1p_precipitation"] +
    beta_business[s] *
    X[, "z_business_count"] +
    beta_population[s] *
    X[, "z_resident_population"] +
    beta_coverage[s] *
    X[, "z_built_coverage"] +
    beta_height[s] *
    X[, "z_structure_height"]
  
  # Protect prior simulations from numerical overflow.
  # Values beyond exp(20) are already implausibly huge for
  # this application, so this cap is only numerical safety.
  eta_common_safe <- pmin(
    eta_common,
    20
  )
  
  mu_common <- exp(
    eta_common_safe
  )
  
  # ----------------------------------------------------------
  # M1: Poisson
  # ----------------------------------------------------------
  
  y_poisson <- rpois(
    N,
    lambda = mu_common
  )
  
  summary_list[[summary_counter]] <- tibble(
    model = "M1 Poisson",
    draw = s,
    mean_count = mean(y_poisson),
    median_count = median(y_poisson),
    q95_count = unname(
      quantile(y_poisson, 0.95)
    ),
    maximum_count = max(y_poisson),
    proportion_gt_100k =
      mean(y_poisson > 100000),
    proportion_gt_1m =
      mean(y_poisson > 1000000)
  )
  
  summary_counter <- summary_counter + 1
  
  sample_list[[sample_counter]] <- tibble(
    model = "M1 Poisson",
    draw = s,
    count =
      y_poisson[
        sample_indices
      ]
  )
  
  sample_counter <- sample_counter + 1
  
  # ----------------------------------------------------------
  # M2: Negative Binomial
  #
  # R parameterisation:
  # Var(Y) = mu + mu^2 / size
  #
  # size = phi
  # ----------------------------------------------------------
  
  y_nb <- rnbinom(
    N,
    size = phi_draws[s],
    mu = mu_common
  )
  
  summary_list[[summary_counter]] <- tibble(
    model = "M2 Negative Binomial",
    draw = s,
    mean_count = mean(y_nb),
    median_count = median(y_nb),
    q95_count = unname(
      quantile(y_nb, 0.95)
    ),
    maximum_count = max(y_nb),
    proportion_gt_100k =
      mean(y_nb > 100000),
    proportion_gt_1m =
      mean(y_nb > 1000000)
  )
  
  summary_counter <- summary_counter + 1
  
  sample_list[[sample_counter]] <- tibble(
    model = "M2 Negative Binomial",
    draw = s,
    count =
      y_nb[
        sample_indices
      ]
  )
  
  sample_counter <- sample_counter + 1
  
  # ----------------------------------------------------------
  # M3: Hierarchical Negative Binomial
  # ----------------------------------------------------------
  
  sensor_effect <- rnorm(
    J,
    mean = 0,
    sd = sigma_sensor[s]
  )
  
  eta_hierarchical <-
    eta_common +
    sensor_effect[
      dat$sensor_index
    ]
  
  eta_hierarchical_safe <- pmin(
    eta_hierarchical,
    20
  )
  
  mu_hierarchical <- exp(
    eta_hierarchical_safe
  )
  
  y_hnb <- rnbinom(
    N,
    size = phi_draws[s],
    mu = mu_hierarchical
  )
  
  summary_list[[summary_counter]] <- tibble(
    model = "M3 Hierarchical Negative Binomial",
    draw = s,
    mean_count = mean(y_hnb),
    median_count = median(y_hnb),
    q95_count = unname(
      quantile(y_hnb, 0.95)
    ),
    maximum_count = max(y_hnb),
    proportion_gt_100k =
      mean(y_hnb > 100000),
    proportion_gt_1m =
      mean(y_hnb > 1000000)
  )
  
  summary_counter <- summary_counter + 1
  
  sample_list[[sample_counter]] <- tibble(
    model = "M3 Hierarchical Negative Binomial",
    draw = s,
    count =
      y_hnb[
        sample_indices
      ]
  )
  
  sample_counter <- sample_counter + 1
}

# ------------------------------------------------------------
# 10. Combine simulation results
# ------------------------------------------------------------

prior_draw_summary <- bind_rows(
  summary_list
)

prior_predictive_sample <- bind_rows(
  sample_list
)

# ============================================================
# PART C — SUMMARISE PRIOR PREDICTIVE DISTRIBUTIONS
# ============================================================

# ------------------------------------------------------------
# 11. Across-draw summaries
# ------------------------------------------------------------

prior_predictive_summary <-
  prior_draw_summary %>%
  group_by(
    model
  ) %>%
  summarise(
    median_of_draw_means =
      median(mean_count),
    
    draw_mean_q025 =
      quantile(
        mean_count,
        0.025
      ),
    
    draw_mean_q975 =
      quantile(
        mean_count,
        0.975
      ),
    
    median_of_draw_medians =
      median(median_count),
    
    median_q95 =
      median(q95_count),
    
    median_draw_maximum =
      median(maximum_count),
    
    q95_draw_maximum =
      quantile(
        maximum_count,
        0.95
      ),
    
    median_prop_gt_100k =
      median(
        proportion_gt_100k
      ),
    
    median_prop_gt_1m =
      median(
        proportion_gt_1m
      ),
    
    .groups = "drop"
  )

write_csv(
  prior_predictive_summary,
  here(
    "outputs",
    "tables",
    "prior_predictive_summary.csv"
  )
)

write_csv(
  prior_draw_summary,
  here(
    "outputs",
    "tables",
    "prior_predictive_draw_summary.csv"
  )
)

# ============================================================
# PART D — FIGURES
# ============================================================

# ------------------------------------------------------------
# 12. Prior predictive counts
#
# Plot on log10 scale because priors deliberately allow a
# wide range of count magnitudes.
# ------------------------------------------------------------

p_counts <- ggplot(
  prior_predictive_sample,
  aes(
    x = count + 1
  )
) +
  geom_histogram(
    bins = 60
  ) +
  scale_x_log10() +
  facet_wrap(
    ~ model,
    scales = "free_y"
  ) +
  labs(
    title = "Prior Predictive Distribution of Daily Pedestrian Counts",
    x = "Prior predictive daily count + 1 (log10 scale)",
    y = "Simulated observations"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "prior_predictive_counts.png"
  ),
  p_counts,
  width = 11,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Distribution of mean counts across prior draws
# ------------------------------------------------------------

p_means <- ggplot(
  prior_draw_summary,
  aes(
    x = mean_count + 1
  )
) +
  geom_histogram(
    bins = 50
  ) +
  scale_x_log10() +
  facet_wrap(
    ~ model,
    scales = "free_y"
  ) +
  labs(
    title = "Prior Predictive Mean Counts Across Prior Draws",
    x = "Mean daily count + 1 (log10 scale)",
    y = "Prior draws"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "prior_predictive_draw_means.png"
  ),
  p_means,
  width = 11,
  height = 6,
  dpi = 300
)

# ============================================================
# PART E — QA OUTPUT
# ============================================================

cat(
  "\n========== PRIOR PREDICTIVE QA ==========\n"
)

cat(
  "Prior draws:",
  n_prior_draws,
  "\n"
)

cat(
  "Observed predictor rows used:",
  N,
  "\n"
)

cat(
  "Sensors represented:",
  J,
  "\n"
)

cat(
  "\nPrior parameter summary:\n"
)

print(
  prior_parameter_summary,
  width = Inf
)

cat(
  "\nPrior predictive summary:\n"
)

print(
  prior_predictive_summary,
  width = Inf
)

cat(
  "\nInterpretation thresholds:\n"
)

cat(
  "  >100,000 pedestrians/day: high but conceivable at very busy sensors\n"
)

cat(
  "  >1,000,000 pedestrians/day: implausible for an individual sensor\n"
)

cat(
  "\nNo observed pedestrian counts were used to generate these simulations.\n"
)

cat(
  "=========================================\n"
)

message(
  "Prior predictive simulation completed successfully."
)