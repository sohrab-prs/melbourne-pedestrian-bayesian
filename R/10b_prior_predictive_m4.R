# ============================================================
# 10b_prior_predictive_m4.R
# Melbourne Pedestrian Bayesian Project
#
# Final prior-predictive validation for M4:
#   Negative Binomial
#   + sensor random intercept
#   + day random effect
#
# Uses the actual final modelling design:
#   N = 18,671 sensor-days
#   J = 53 sensors
#   D = 364 represented dates
#
# No posterior information is used.
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

set.seed(20260917)

# ============================================================
# PART A — DATA
# ============================================================

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
  date = sort(unique(dat$date)),
  day_index = seq_along(sort(unique(dat$date)))
)

D <- nrow(date_lookup)

dat <- dat %>%
  left_join(
    date_lookup,
    by = "date"
  )

stopifnot(
  D == 364,
  sum(is.na(dat$day_index)) == 0
)

sensor <- as.integer(dat$sensor_index)
day <- as.integer(dat$day_index)
dow <- as.integer(dat$dow_index)

holiday <- as.numeric(dat$public_holiday_num)
temp <- as.numeric(dat$z_mean_temp)
rain <- as.numeric(dat$z_log1p_precipitation)

business <- as.numeric(dat$z_business_count)
population <- as.numeric(dat$z_resident_population)
coverage <- as.numeric(dat$z_built_coverage)
height <- as.numeric(dat$z_structure_height)

y_observed <- as.numeric(dat$daily_ped_count)

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

# ============================================================
# PART B — OBSERVED REFERENCE FEATURES
# ============================================================

observed_summary <- tibble(
  
  statistic = c(
    "mean_count",
    "median_count",
    "q95_count",
    "maximum_count",
    "variance_to_mean",
    "proportion_gt_100000",
    "proportion_gt_1000000"
  ),
  
  observed = c(
    mean(y_observed),
    median(y_observed),
    unname(quantile(y_observed, 0.95)),
    max(y_observed),
    var(y_observed) / mean(y_observed),
    mean(y_observed > 100000),
    mean(y_observed > 1000000)
  )
)

# ============================================================
# PART C — PRIORS
# ============================================================

# These are exactly the priors used in the final M4 model:
#
# alpha ~ Normal(log(10000), 0.5^2)
# beta  ~ Normal(0, 0.25^2)
# log(phi) ~ Normal(log(10), 1^2)
# sigma_sensor ~ Half-Normal(0, 0.5)
# sigma_day    ~ Half-Normal(0, 0.5)

n_prior_draws <- 300L

prior_draw_results <- vector(
  "list",
  n_prior_draws
)

cat(
  "\nGenerating M4 prior-predictive datasets...\n"
)

for (d in seq_len(n_prior_draws)) {
  
  alpha <- rnorm(
    1,
    mean = log(10000),
    sd = 0.5
  )
  
  beta_dow <- c(
    0,
    rnorm(
      6,
      mean = 0,
      sd = 0.25
    )
  )
  
  beta_holiday <- rnorm(1, 0, 0.25)
  beta_temp <- rnorm(1, 0, 0.25)
  beta_rain <- rnorm(1, 0, 0.25)
  
  beta_business <- rnorm(1, 0, 0.25)
  beta_population <- rnorm(1, 0, 0.25)
  beta_coverage <- rnorm(1, 0, 0.25)
  beta_height <- rnorm(1, 0, 0.25)
  
  phi <- exp(
    rnorm(
      1,
      mean = log(10),
      sd = 1
    )
  )
  
  sigma_sensor <- abs(
    rnorm(
      1,
      mean = 0,
      sd = 0.5
    )
  )
  
  sigma_day <- abs(
    rnorm(
      1,
      mean = 0,
      sd = 0.5
    )
  )
  
  u_sensor <- rnorm(
    J,
    mean = 0,
    sd = sigma_sensor
  )
  
  u_day <- rnorm(
    D,
    mean = 0,
    sd = sigma_day
  )
  
  eta <-
    alpha +
    beta_dow[dow] +
    beta_holiday * holiday +
    beta_temp * temp +
    beta_rain * rain +
    beta_business * business +
    beta_population * population +
    beta_coverage * coverage +
    beta_height * height +
    u_sensor[sensor] +
    u_day[day]
  
  mu <- exp(
    eta
  )
  
  stopifnot(
    all(is.finite(mu)),
    all(mu > 0)
  )
  
  y_rep <- rnbinom(
    N,
    size = phi,
    mu = mu
  )
  
  prior_draw_results[[d]] <- tibble(
    
    draw =
      d,
    
    alpha =
      alpha,
    
    phi =
      phi,
    
    sigma_sensor =
      sigma_sensor,
    
    sigma_day =
      sigma_day,
    
    mean_count =
      mean(y_rep),
    
    median_count =
      median(y_rep),
    
    q95_count =
      unname(
        quantile(
          y_rep,
          0.95
        )
      ),
    
    maximum_count =
      max(y_rep),
    
    variance_to_mean =
      var(y_rep) /
      mean(y_rep),
    
    proportion_gt_100000 =
      mean(
        y_rep > 100000
      ),
    
    proportion_gt_1000000 =
      mean(
        y_rep > 1000000
      )
  )
}

prior_draws <- bind_rows(
  prior_draw_results
)

write_csv(
  prior_draws,
  here(
    "outputs",
    "tables",
    "m4_prior_predictive_draws.csv"
  )
)

# ============================================================
# PART D — SUMMARY
# ============================================================

summary_variables <- c(
  "mean_count",
  "median_count",
  "q95_count",
  "maximum_count",
  "variance_to_mean",
  "proportion_gt_100000",
  "proportion_gt_1000000"
)

prior_summary <- prior_draws %>%
  
  select(
    all_of(
      summary_variables
    )
  ) %>%
  
  pivot_longer(
    cols = everything(),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  
  group_by(
    statistic
  ) %>%
  
  summarise(
    
    prior_q025 =
      unname(
        quantile(
          value,
          0.025
        )
      ),
    
    prior_median =
      median(
        value
      ),
    
    prior_q975 =
      unname(
        quantile(
          value,
          0.975
        )
      ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    observed_summary,
    by = "statistic"
  )

write_csv(
  prior_summary,
  here(
    "outputs",
    "tables",
    "m4_prior_predictive_summary.csv"
  )
)

# ============================================================
# PART E — FIGURE
# ============================================================

plot_data <- prior_draws %>%
  
  select(
    mean_count,
    q95_count,
    maximum_count,
    variance_to_mean
  ) %>%
  
  pivot_longer(
    cols = everything(),
    names_to = "statistic",
    values_to = "value"
  )

p <- ggplot(
  plot_data,
  aes(
    x = value
  )
) +
  geom_histogram(
    bins = 40
  ) +
  facet_wrap(
    ~ statistic,
    scales = "free",
    ncol = 2
  ) +
  labs(
    title =
      "M4 Prior Predictive Distributions",
    subtitle =
      "Final M4 priors evaluated using the actual sensor-day design",
    x =
      "Prior predictive statistic",
    y =
      "Prior draws"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_prior_predictive_checks.png"
  ),
  p,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================
# PART F — QA
# ============================================================

cat(
  "\n========== FINAL M4 PRIOR PREDICTIVE CHECK ==========\n"
)

cat(
  "\nPrior predictive draws:",
  n_prior_draws,
  "\n"
)

cat(
  "Observations per replicated dataset:",
  N,
  "\n"
)

cat(
  "Sensors:",
  J,
  "\n"
)

cat(
  "Dates:",
  D,
  "\n"
)

cat(
  "\nObserved vs prior-predictive summaries:\n"
)

print(
  prior_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nAdditional prior diagnostics:\n"
)

cat(
  "Median sigma_sensor:",
  median(
    prior_draws$sigma_sensor
  ),
  "\n"
)

cat(
  "Median sigma_day:",
  median(
    prior_draws$sigma_day
  ),
  "\n"
)

cat(
  "Median phi:",
  median(
    prior_draws$phi
  ),
  "\n"
)

cat(
  "Fraction of prior datasets with any count > 1,000,000:",
  mean(
    prior_draws$maximum_count >
      1000000
  ),
  "\n"
)

cat(
  "=====================================================\n"
)

message(
  "Final M4 prior-predictive validation completed."
)