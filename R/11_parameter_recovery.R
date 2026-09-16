# ============================================================
# 11_parameter_recovery.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Validate implementation of the hierarchical Negative
# Binomial model by simulating data from known parameters
# and assessing whether JAGS/MCMC recovers those parameters.
#
# The real pedestrian response is NOT used to generate the
# synthetic outcome.
#
# The real predictor structure is retained so the recovery
# experiment reflects the actual model-design problem.
#
# Outputs:
#   outputs/tables/parameter_recovery_results.csv
#   outputs/tables/parameter_recovery_convergence.csv
#   outputs/tables/parameter_recovery_sensor_effects.csv
#   outputs/figures/parameter_recovery_fixed_effects.png
#   outputs/figures/parameter_recovery_sensor_effects.png
#   outputs/mcmc/parameter_recovery_samples.rds
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(rjags)
library(coda)
library(ggplot2)

set.seed(20260916)

# ------------------------------------------------------------
# 1. Create output directories
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

dir.create(
  here("outputs", "mcmc"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PART A — LOAD PREDICTOR STRUCTURE
# ============================================================

# ------------------------------------------------------------
# 2. Load selected modelling data
# ------------------------------------------------------------

dat_full <- read_csv(
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
  nrow(dat_full) == 18671
)

stopifnot(
  n_distinct(
    dat_full$sensor_index
  ) == 53
)

stopifnot(
  sum(is.na(dat_full)) == 0
)

# ------------------------------------------------------------
# 3. Construct a manageable but realistic recovery design
#
# Use 80 observations per sensor.
#
# This preserves:
#   - all 53 sensors;
#   - realistic temporal predictor combinations;
#   - actual spatial predictor correlations;
#
# while avoiding unnecessary computation for a validation
# experiment.
# ------------------------------------------------------------

recovery_dat <- dat_full %>%
  group_by(
    sensor_index
  ) %>%
  slice_sample(
    n = 80,
    replace = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    sensor_index,
    date
  )

stopifnot(
  nrow(recovery_dat) ==
    53 * 80
)

stopifnot(
  n_distinct(
    recovery_dat$sensor_index
  ) == 53
)

N <- nrow(recovery_dat)

J <- 53

# ------------------------------------------------------------
# Sensor-level predictor table
#
# Spatial predictors vary only between sensors, so the
# hierarchically centred JAGS model receives one value per
# sensor rather than repeating them for every sensor-day.
# ------------------------------------------------------------

sensor_level_jags <- recovery_dat %>%
  select(
    sensor_index,
    z_business_count,
    z_resident_population,
    z_built_coverage,
    z_structure_height
  ) %>%
  distinct() %>%
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
  nrow(sensor_level_jags) == J
)

stopifnot(
  anyDuplicated(
    sensor_level_jags$sensor_index
  ) == 0
)

stopifnot(
  identical(
    sensor_level_jags$sensor_index,
    seq_len(J)
  )
)

# ============================================================
# PART B — DEFINE KNOWN SYNTHETIC TRUTH
# ============================================================

# ------------------------------------------------------------
# 4. Fixed-effect truth
#
# These values are chosen to be:
#   - substantively plausible;
#   - well inside the prior distributions;
#   - non-zero enough to test parameter recovery.
#
# They are NOT estimates from the observed pedestrian counts.
# ------------------------------------------------------------

true_alpha <- log(10000)

# Monday is reference = 0
true_beta_dow <- c(
  Monday    = 0.00,
  Tuesday   = 0.05,
  Wednesday = 0.10,
  Thursday  = 0.15,
  Friday    = 0.18,
  Saturday  = 0.16,
  Sunday    = -0.05
)

true_beta_holiday <- -0.15

true_beta_temp <- 0.10

true_beta_rain <- -0.20

true_beta_business <- 0.20

true_beta_population <- -0.10

true_beta_coverage <- 0.15

true_beta_height <- 0.12

# Negative Binomial dispersion
true_phi <- 12

# Sensor heterogeneity
true_sigma_sensor <- 0.35

# ------------------------------------------------------------
# 5. Simulate true sensor random effects
# ------------------------------------------------------------

true_u_sensor <- rnorm(
  J,
  mean = 0,
  sd = true_sigma_sensor
)

# Centre realised random effects around zero.
#
# This avoids a finite-sample random intercept shift being
# absorbed into the synthetic intercept.
true_u_sensor <-
  true_u_sensor -
  mean(true_u_sensor)

# ============================================================
# PART C — SIMULATE SYNTHETIC RESPONSE
# ============================================================

# ------------------------------------------------------------
# 6. Construct true linear predictor
# ------------------------------------------------------------

eta_true <-
  true_alpha +
  true_beta_dow[
    recovery_dat$dow_index
  ] +
  true_beta_holiday *
  recovery_dat$public_holiday_num +
  true_beta_temp *
  recovery_dat$z_mean_temp +
  true_beta_rain *
  recovery_dat$z_log1p_precipitation +
  true_beta_business *
  recovery_dat$z_business_count +
  true_beta_population *
  recovery_dat$z_resident_population +
  true_beta_coverage *
  recovery_dat$z_built_coverage +
  true_beta_height *
  recovery_dat$z_structure_height +
  true_u_sensor[
    recovery_dat$sensor_index
  ]

mu_true <- exp(
  eta_true
)

# ------------------------------------------------------------
# 7. Generate synthetic Negative Binomial counts
# ------------------------------------------------------------

y_sim <- rnbinom(
  N,
  size = true_phi,
  mu = mu_true
)

stopifnot(
  length(y_sim) == N
)

stopifnot(
  all(
    y_sim >= 0
  )
)

cat(
  "\nSynthetic response summary:\n"
)

print(
  summary(y_sim)
)

cat(
  "\nSynthetic variance-to-mean ratio:",
  var(y_sim) / mean(y_sim),
  "\n"
)

# ============================================================
# PART D — PREPARE JAGS DATA
# ============================================================

# ------------------------------------------------------------
# 8. JAGS data list
# ------------------------------------------------------------

jags_data <- list(
  
  N = N,
  
  J = J,
  
  y =
    as.integer(
      y_sim
    ),
  
  sensor =
    as.integer(
      recovery_dat$sensor_index
    ),
  
  dow =
    as.integer(
      recovery_dat$dow_index
    ),
  
  holiday =
    as.numeric(
      recovery_dat$public_holiday_num
    ),
  
  temp =
    as.numeric(
      recovery_dat$z_mean_temp
    ),
  
  rain =
    as.numeric(
      recovery_dat$z_log1p_precipitation
    ),
  
  # ----------------------------------------------------------
  # Sensor-level spatial predictors
  # ----------------------------------------------------------
  
  business_sensor =
    as.numeric(
      sensor_level_jags$z_business_count
    ),
  
  population_sensor =
    as.numeric(
      sensor_level_jags$z_resident_population
    ),
  
  coverage_sensor =
    as.numeric(
      sensor_level_jags$z_built_coverage
    ),
  
  height_sensor =
    as.numeric(
      sensor_level_jags$z_structure_height
    ),
  
  # Spatial coefficient prior
  mu_spatial =
    rep(
      0,
      4
    ),
  
  Tau_spatial =
    diag(
      rep(
        16,
        4
      )
    )
)

# ------------------------------------------------------------
# 9. Chain initial values
# ------------------------------------------------------------

make_inits <- function(
    seed,
    alpha_init,
    phi_init,
    sigma_init
) {
  
  set.seed(seed)
  
  beta_spatial_init <- rnorm(
    4,
    mean = 0,
    sd = 0.10
  )
  
  sensor_mean_init <-
    alpha_init +
    beta_spatial_init[1] *
    sensor_level_jags$z_business_count +
    beta_spatial_init[2] *
    sensor_level_jags$z_resident_population +
    beta_spatial_init[3] *
    sensor_level_jags$z_built_coverage +
    beta_spatial_init[4] *
    sensor_level_jags$z_structure_height
  
  list(
    
    alpha =
      alpha_init,
    
    beta_spatial =
      beta_spatial_init,
    
    log_phi =
      log(
        phi_init
      ),
    
    sigma_sensor =
      sigma_init,
    
    sensor_intercept =
      sensor_mean_init +
      rnorm(
        J,
        mean = 0,
        sd = sigma_init
      ),
    
    .RNG.name =
      "base::Mersenne-Twister",
    
    .RNG.seed =
      seed
  )
}

inits <- list(
  
  make_inits(
    seed = 10101,
    alpha_init = log(7000),
    phi_init = 6,
    sigma_init = 0.20
  ),
  
  make_inits(
    seed = 20202,
    alpha_init = log(9000),
    phi_init = 10,
    sigma_init = 0.35
  ),
  
  make_inits(
    seed = 30303,
    alpha_init = log(12000),
    phi_init = 18,
    sigma_init = 0.50
  ),
  
  make_inits(
    seed = 40404,
    alpha_init = log(15000),
    phi_init = 30,
    sigma_init = 0.70
  )
)

# ============================================================
# PART E — FIT SYNTHETIC DATA WITH JAGS/MCMC
# ============================================================

# ------------------------------------------------------------
# 10. Compile model
# ------------------------------------------------------------

model_file <- here(
  "models",
  "jags",
  "hierarchical_nb_recovery.jags"
)

stopifnot(
  file.exists(model_file)
)

recovery_model <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits,
  n.chains = 4,
  n.adapt = 5000,
  quiet = TRUE
)

# ------------------------------------------------------------
# 11. Burn-in
# ------------------------------------------------------------

update(
  recovery_model,
  n.iter = 5000,
  progress.bar = "none"
)

# ------------------------------------------------------------
# 12. Posterior sampling
# ------------------------------------------------------------

parameters_to_monitor <- c(
  "alpha",
  "beta_dow",
  "beta_holiday",
  "beta_temp",
  "beta_rain",
  "beta_business",
  "beta_population",
  "beta_coverage",
  "beta_height",
  "phi",
  "sigma_sensor",
  "u_sensor"
)

recovery_samples <- coda.samples(
  recovery_model,
  variable.names =
    parameters_to_monitor,
  n.iter = 10000,
  thin = 1,
  progress.bar = "none"
)

saveRDS(
  recovery_samples,
  here(
    "outputs",
    "mcmc",
    "parameter_recovery_samples.rds"
  )
)

# ============================================================
# PART F — FIXED-PARAMETER RECOVERY
# ============================================================

# ------------------------------------------------------------
# 13. Define true values using JAGS parameter names
# ------------------------------------------------------------

truth <- tibble(
  parameter = c(
    "alpha",
    "beta_dow[2]",
    "beta_dow[3]",
    "beta_dow[4]",
    "beta_dow[5]",
    "beta_dow[6]",
    "beta_dow[7]",
    "beta_holiday",
    "beta_temp",
    "beta_rain",
    "beta_business",
    "beta_population",
    "beta_coverage",
    "beta_height",
    "phi",
    "sigma_sensor"
  ),
  
  true_value = c(
    true_alpha,
    true_beta_dow["Tuesday"],
    true_beta_dow["Wednesday"],
    true_beta_dow["Thursday"],
    true_beta_dow["Friday"],
    true_beta_dow["Saturday"],
    true_beta_dow["Sunday"],
    true_beta_holiday,
    true_beta_temp,
    true_beta_rain,
    true_beta_business,
    true_beta_population,
    true_beta_coverage,
    true_beta_height,
    true_phi,
    true_sigma_sensor
  )
)

# ------------------------------------------------------------
# 14. Convert MCMC samples to matrix
# ------------------------------------------------------------

posterior_matrix <- as.matrix(
  recovery_samples
)

# Ensure all required parameters exist
stopifnot(
  all(
    truth$parameter %in%
      colnames(
        posterior_matrix
      )
  )
)

# ------------------------------------------------------------
# 15. Posterior summaries for known parameters
# ------------------------------------------------------------

recovery_results <- truth %>%
  rowwise() %>%
  mutate(
    posterior_mean =
      mean(
        posterior_matrix[
          ,
          parameter
        ]
      ),
    
    posterior_sd =
      sd(
        posterior_matrix[
          ,
          parameter
        ]
      ),
    
    q025 =
      unname(
        quantile(
          posterior_matrix[
            ,
            parameter
          ],
          0.025
        )
      ),
    
    median =
      unname(
        quantile(
          posterior_matrix[
            ,
            parameter
          ],
          0.50
        )
      ),
    
    q975 =
      unname(
        quantile(
          posterior_matrix[
            ,
            parameter
          ],
          0.975
        )
      ),
    
    recovered_95 =
      true_value >= q025 &
      true_value <= q975,
    
    standardised_error =
      (
        posterior_mean -
          true_value
      ) /
      posterior_sd
  ) %>%
  ungroup()

# ============================================================
# PART G — CONVERGENCE DIAGNOSTICS
# ============================================================

# ------------------------------------------------------------
# 16. Core-parameter MCMC object
# ------------------------------------------------------------

core_parameters <- truth$parameter

core_samples <-
  recovery_samples[
    ,
    core_parameters
  ]

# ------------------------------------------------------------
# 17. Gelman-Rubin R-hat
# ------------------------------------------------------------

rhat_values <-
  gelman.diag(
    core_samples,
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[
    ,
    1
  ]

# ------------------------------------------------------------
# 18. Effective sample size
# ------------------------------------------------------------

ess_values <- effectiveSize(
  core_samples
)

convergence_results <- tibble(
  parameter =
    names(rhat_values),
  
  rhat =
    as.numeric(
      rhat_values
    ),
  
  effective_sample_size =
    as.numeric(
      ess_values[
        names(rhat_values)
      ]
    )
)

# ============================================================
# PART H — SENSOR RANDOM-EFFECT RECOVERY
# ============================================================

# ------------------------------------------------------------
# 19. Extract posterior random effects
# ------------------------------------------------------------

sensor_effect_results <- lapply(
  seq_len(J),
  function(j) {
    
    parameter_name <-
      paste0(
        "u_sensor[",
        j,
        "]"
      )
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    tibble(
      sensor_index = j,
      
      true_effect =
        true_u_sensor[j],
      
      posterior_mean =
        mean(x),
      
      q025 =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      
      q975 =
        unname(
          quantile(
            x,
            0.975
          )
        ),
      
      recovered_95 =
        true_u_sensor[j] >=
        quantile(
          x,
          0.025
        ) &
        true_u_sensor[j] <=
        quantile(
          x,
          0.975
        )
    )
  }
) %>%
  bind_rows()

sensor_effect_correlation <- cor(
  sensor_effect_results$true_effect,
  sensor_effect_results$posterior_mean,
  method = "pearson"
)

sensor_effect_coverage <- mean(
  sensor_effect_results$recovered_95
)

# ============================================================
# PART I — SAVE TABLES
# ============================================================

write_csv(
  recovery_results,
  here(
    "outputs",
    "tables",
    "parameter_recovery_results.csv"
  )
)

write_csv(
  convergence_results,
  here(
    "outputs",
    "tables",
    "parameter_recovery_convergence.csv"
  )
)

write_csv(
  sensor_effect_results,
  here(
    "outputs",
    "tables",
    "parameter_recovery_sensor_effects.csv"
  )
)

# ============================================================
# PART J — FIGURES
# ============================================================

# ------------------------------------------------------------
# 20. Fixed-parameter recovery plot
# ------------------------------------------------------------

fixed_plot_data <- recovery_results %>%
  filter(
    !parameter %in%
      c(
        "phi",
        "sigma_sensor",
        "alpha"
      )
  )

p_fixed <- ggplot(
  fixed_plot_data,
  aes(
    x = true_value,
    y = posterior_mean
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = q025,
      ymax = q975
    ),
    width = 0
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_text(
    aes(
      label = parameter
    ),
    nudge_x = 0.01,
    nudge_y = 0.01,
    check_overlap = TRUE
  ) +
  labs(
    title = "Parameter Recovery: Fixed Effects",
    x = "True simulated value",
    y = "Posterior mean and 95% credible interval"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "parameter_recovery_fixed_effects.png"
  ),
  p_fixed,
  width = 9,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 21. Sensor-effect recovery plot
# ------------------------------------------------------------

p_sensor <- ggplot(
  sensor_effect_results,
  aes(
    x = true_effect,
    y = posterior_mean
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = q025,
      ymax = q975
    ),
    width = 0
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    title = "Parameter Recovery: Sensor Random Effects",
    subtitle =
      paste0(
        "Correlation = ",
        round(
          sensor_effect_correlation,
          3
        )
      ),
    x = "True simulated sensor effect",
    y = "Posterior mean and 95% credible interval"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "parameter_recovery_sensor_effects.png"
  ),
  p_sensor,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# PART K — QA OUTPUT
# ============================================================

cat(
  "\n========== PARAMETER RECOVERY QA ==========\n"
)

cat(
  "Synthetic observations:",
  N,
  "\n"
)

cat(
  "Synthetic sensors:",
  J,
  "\n"
)

cat(
  "\nKnown-parameter recovery:\n"
)

print(
  recovery_results,
  width = Inf
)

cat(
  "\nConvergence diagnostics:\n"
)

print(
  convergence_results,
  width = Inf
)

cat(
  "\nMaximum R-hat:",
  max(
    convergence_results$rhat
  ),
  "\n"
)

cat(
  "Minimum effective sample size:",
  min(
    convergence_results$effective_sample_size
  ),
  "\n"
)

cat(
  "\nProportion of core true parameters inside 95% credible intervals:",
  mean(
    recovery_results$recovered_95
  ),
  "\n"
)

cat(
  "\nSensor random-effect recovery correlation:",
  sensor_effect_correlation,
  "\n"
)

cat(
  "Sensor random-effect 95% interval coverage:",
  sensor_effect_coverage,
  "\n"
)

cat(
  "===========================================\n"
)


# ------------------------------------------------------------
# Formal recovery QA criteria
# ------------------------------------------------------------

max_rhat <-
  max(
    convergence_results$rhat
  )

min_ess <-
  min(
    convergence_results$effective_sample_size
  )

core_coverage <-
  mean(
    recovery_results$recovered_95
  )

max_abs_standardised_error <-
  max(
    abs(
      recovery_results$standardised_error
    )
  )

recovery_pass <-
  max_rhat < 1.05 &
  min_ess >= 400 &
  core_coverage >= 0.90 &
  sensor_effect_correlation >= 0.90

cat(
  "\nFormal recovery criteria:\n"
)

cat(
  "  Maximum R-hat < 1.05:",
  max_rhat < 1.05,
  "\n"
)

cat(
  "  Minimum ESS >= 400:",
  min_ess >= 400,
  "\n"
)

cat(
  "  Core 95% interval coverage >= 90%:",
  core_coverage >= 0.90,
  "\n"
)

cat(
  "  Sensor-effect correlation >= 0.90:",
  sensor_effect_correlation >= 0.90,
  "\n"
)

cat(
  "  Maximum absolute standardised recovery error:",
  max_abs_standardised_error,
  "\n"
)

cat(
  "\nOVERALL PARAMETER RECOVERY PASS:",
  recovery_pass,
  "\n"
)

message(
  "Hierarchical Negative Binomial parameter recovery completed."
)