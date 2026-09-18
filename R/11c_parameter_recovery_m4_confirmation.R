# ============================================================
# 11b_parameter_recovery_m4.R
# Melbourne Pedestrian Bayesian Project
#
# Final parameter-recovery experiment for M4.
#
# Uses:
#   - actual 18,671-row observation design
#   - actual 53 sensors
#   - actual 364 represented dates
#   - actual standardized predictor values
#
# Synthetic response is generated from known parameters.
# The unchanged final M4 JAGS model is then fitted.
#
# Recovery assessed for:
#   - 17 core model parameters
#   - 53 sensor random effects
#   - 364 day random effects
#
# Computational strategy:
#   - 4 independent JAGS chains in parallel
#   - same 7,000 retained iterations used for accepted M4
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(rjags)
library(coda)
library(parallel)
library(ggplot2)

set.seed(20260918)

# ============================================================
# PART A — ACTUAL DESIGN
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

sensor_level <- dat %>%
  
  select(
    sensor_index,
    location_id,
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
  nrow(sensor_level) == J,
  all(
    sensor_level$sensor_index ==
      seq_len(J)
  )
)

# ============================================================
# PART B — TRUE SYNTHETIC PARAMETERS
# ============================================================

truth_core <- c(
  
  alpha =
    log(10000),
  
  "beta_dow[2]" =
    0.05,
  
  "beta_dow[3]" =
    0.10,
  
  "beta_dow[4]" =
    0.15,
  
  "beta_dow[5]" =
    0.18,
  
  "beta_dow[6]" =
    0.16,
  
  "beta_dow[7]" =
    -0.05,
  
  beta_holiday =
    -0.15,
  
  beta_temp =
    0.10,
  
  beta_rain =
    -0.20,
  
  beta_business =
    0.20,
  
  beta_population =
    -0.10,
  
  beta_coverage =
    0.15,
  
  beta_height =
    0.12,
  
  phi =
    12,
  
  sigma_sensor =
    0.35,
  
  sigma_day =
    0.10
)

# ------------------------------------------------------------
# Generate latent effects independently from the fitted model.
# ------------------------------------------------------------

true_u_sensor <- rnorm(
  J,
  mean = 0,
  sd = truth_core[
    "sigma_sensor"
  ]
)

true_day_effect <- rnorm(
  D,
  mean = 0,
  sd = truth_core[
    "sigma_day"
  ]
)

# ============================================================
# PART C — GENERATE SYNTHETIC RESPONSE
# ============================================================

beta_dow_true <- c(
  0,
  truth_core[
    paste0(
      "beta_dow[",
      2:7,
      "]"
    )
  ]
)

sensor <- as.integer(
  dat$sensor_index
)

day <- as.integer(
  dat$day_index
)

dow <- as.integer(
  dat$dow_index
)

eta_true <-
  
  truth_core[
    "alpha"
  ] +
  
  beta_dow_true[
    dow
  ] +
  
  truth_core[
    "beta_holiday"
  ] *
  dat$public_holiday_num +
  
  truth_core[
    "beta_temp"
  ] *
  dat$z_mean_temp +
  
  truth_core[
    "beta_rain"
  ] *
  dat$z_log1p_precipitation +
  
  truth_core[
    "beta_business"
  ] *
  dat$z_business_count +
  
  truth_core[
    "beta_population"
  ] *
  dat$z_resident_population +
  
  truth_core[
    "beta_coverage"
  ] *
  dat$z_built_coverage +
  
  truth_core[
    "beta_height"
  ] *
  dat$z_structure_height +
  
  true_u_sensor[
    sensor
  ] +
  
  true_day_effect[
    day
  ]

mu_true <- exp(
  eta_true
)

stopifnot(
  all(
    is.finite(
      mu_true
    )
  ),
  all(
    mu_true > 0
  )
)

y_sim <- rnbinom(
  N,
  size =
    truth_core[
      "phi"
    ],
  mu =
    mu_true
)

stopifnot(
  all(
    is.finite(
      y_sim
    )
  )
)

# ============================================================
# PART D — JAGS DATA
# ============================================================

jags_data <- list(
  
  N = N,
  
  J = J,
  
  D = D,
  
  y =
    as.integer(
      y_sim
    ),
  
  sensor =
    sensor,
  
  day =
    day,
  
  dow =
    dow,
  
  holiday =
    as.numeric(
      dat$public_holiday_num
    ),
  
  temp =
    as.numeric(
      dat$z_mean_temp
    ),
  
  rain =
    as.numeric(
      dat$z_log1p_precipitation
    ),
  
  business_sensor =
    as.numeric(
      sensor_level$z_business_count
    ),
  
  population_sensor =
    as.numeric(
      sensor_level$z_resident_population
    ),
  
  coverage_sensor =
    as.numeric(
      sensor_level$z_built_coverage
    ),
  
  height_sensor =
    as.numeric(
      sensor_level$z_structure_height
    ),
  
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

model_file <- here(
  "models",
  "jags",
  "m4_sensor_day_negbin.jags"
)

stopifnot(
  file.exists(
    model_file
  )
)

# ============================================================
# PART E — INITIAL VALUES
# ============================================================

make_inits <- function(
    seed,
    alpha_init,
    phi_init,
    sigma_sensor_init,
    sigma_day_init
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
    sensor_level$z_business_count +
    
    beta_spatial_init[2] *
    sensor_level$z_resident_population +
    
    beta_spatial_init[3] *
    sensor_level$z_built_coverage +
    
    beta_spatial_init[4] *
    sensor_level$z_structure_height
  
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
      sigma_sensor_init,
    
    sigma_day =
      sigma_day_init,
    
    sensor_intercept =
      sensor_mean_init +
      rnorm(
        J,
        0,
        sigma_sensor_init
      ),
    
    day_effect =
      rnorm(
        D,
        0,
        sigma_day_init
      ),
    
    .RNG.name =
      "base::Mersenne-Twister",
    
    .RNG.seed =
      seed
  )
}

inits_list <- list(
  
  make_inits(
    91011,
    log(7000),
    6,
    0.20,
    0.05
  ),
  
  make_inits(
    92022,
    log(9000),
    10,
    0.35,
    0.10
  ),
  
  make_inits(
    93033,
    log(12000),
    18,
    0.55,
    0.20
  ),
  
  make_inits(
    94044,
    log(15000),
    30,
    0.75,
    0.30
  )
)

# ============================================================
# PART F — PARAMETERS / SETTINGS
# ============================================================

core_parameters <- names(
  truth_core
)

sensor_parameters <- paste0(
  "u_sensor[",
  seq_len(J),
  "]"
)

day_parameters <- paste0(
  "day_effect[",
  seq_len(D),
  "]"
)

parameters_to_monitor <- c(
  core_parameters,
  "u_sensor",
  "day_effect"
)

diagnostic_parameters <- c(
  core_parameters,
  sensor_parameters,
  day_parameters
)

n_chains <- 4L

n_adapt <- 1000L
n_burn <- 1500L
n_post <- 7000L

# ============================================================
# PART G — PARALLEL FIT
# ============================================================

cat(
  "\n[",
  format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  ),
  "] Starting full-design M4 parameter recovery...\n",
  sep = ""
)

cl <- makePSOCKcluster(
  n_chains,
  outfile = ""
)

on.exit(
  try(
    stopCluster(cl),
    silent = TRUE
  ),
  add = TRUE
)

clusterEvalQ(
  cl,
  {
    library(rjags)
    library(coda)
    NULL
  }
)

clusterExport(
  cl,
  varlist = c(
    "jags_data",
    "model_file",
    "inits_list",
    "parameters_to_monitor",
    "n_adapt",
    "n_burn",
    "n_post"
  ),
  envir = environment()
)

chain_results <- parLapply(
  cl,
  seq_len(n_chains),
  function(chain_id) {
    
    cat(
      "[",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      "] Recovery chain ",
      chain_id,
      " starting\n",
      sep = ""
    )
    
    model <- jags.model(
      
      file =
        model_file,
      
      data =
        jags_data,
      
      inits =
        inits_list[[chain_id]],
      
      n.chains =
        1,
      
      n.adapt =
        n_adapt,
      
      quiet =
        TRUE
    )
    
    update(
      model,
      n.iter =
        n_burn,
      progress.bar =
        "none"
    )
    
    samples <- coda.samples(
      
      model,
      
      variable.names =
        parameters_to_monitor,
      
      n.iter =
        n_post,
      
      thin =
        1,
      
      progress.bar =
        "none"
    )
    
    cat(
      "[",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      "] Recovery chain ",
      chain_id,
      " completed\n",
      sep = ""
    )
    
    samples[[1]]
  }
)

stopCluster(
  cl
)

recovery_samples <- mcmc.list(
  chain_results
)

dir.create(
  here("outputs", "mcmc"),
  recursive = TRUE,
  showWarnings = FALSE
)

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

saveRDS(
  recovery_samples,
  here(
    "outputs",
    "mcmc",
    "m4_parameter_recovery_confirmation_samples.rds"
  )
)

posterior_matrix <- as.matrix(
  recovery_samples
)

# ============================================================
# PART H — CORE PARAMETER RECOVERY
# ============================================================

core_summary <- lapply(
  core_parameters,
  function(parameter_name) {
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    true_value <-
      unname(
        truth_core[
          parameter_name
        ]
      )
    
    posterior_mean <-
      mean(x)
    
    posterior_sd <-
      sd(x)
    
    q025 <-
      unname(
        quantile(
          x,
          0.025
        )
      )
    
    q975 <-
      unname(
        quantile(
          x,
          0.975
        )
      )
    
    tibble(
      
      parameter =
        parameter_name,
      
      true_value =
        true_value,
      
      posterior_mean =
        posterior_mean,
      
      posterior_sd =
        posterior_sd,
      
      q025 =
        q025,
      
      median =
        median(x),
      
      q975 =
        q975,
      
      standardized_error =
        (
          posterior_mean -
            true_value
        ) /
        posterior_sd,
      
      covered_95 =
        true_value >=
        q025 &
        true_value <=
        q975
    )
  }
) %>%
  bind_rows()

write_csv(
  core_summary,
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_core.csv"
  )
)

# ============================================================
# PART I — SENSOR EFFECT RECOVERY
# ============================================================

sensor_recovery <- lapply(
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
    
    q025 <-
      unname(
        quantile(
          x,
          0.025
        )
      )
    
    q975 <-
      unname(
        quantile(
          x,
          0.975
        )
      )
    
    tibble(
      
      sensor_index =
        j,
      
      true_effect =
        true_u_sensor[j],
      
      posterior_mean =
        mean(x),
      
      q025 =
        q025,
      
      median =
        median(x),
      
      q975 =
        q975,
      
      covered_95 =
        true_u_sensor[j] >=
        q025 &
        true_u_sensor[j] <=
        q975
    )
  }
) %>%
  bind_rows()

write_csv(
  sensor_recovery,
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_sensor_effects.csv"
  )
)

# ============================================================
# PART J — DAY EFFECT RECOVERY
# ============================================================

day_recovery <- lapply(
  seq_len(D),
  function(d) {
    
    parameter_name <-
      paste0(
        "day_effect[",
        d,
        "]"
      )
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    q025 <-
      unname(
        quantile(
          x,
          0.025
        )
      )
    
    q975 <-
      unname(
        quantile(
          x,
          0.975
        )
      )
    
    tibble(
      
      day_index =
        d,
      
      date =
        date_lookup$date[d],
      
      true_effect =
        true_day_effect[d],
      
      posterior_mean =
        mean(x),
      
      q025 =
        q025,
      
      median =
        median(x),
      
      q975 =
        q975,
      
      covered_95 =
        true_day_effect[d] >=
        q025 &
        true_day_effect[d] <=
        q975
    )
  }
) %>%
  bind_rows()

write_csv(
  day_recovery,
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_day_effects.csv"
  )
)

# ============================================================
# PART K — CONVERGENCE
# ============================================================

diagnostic_samples <-
  recovery_samples[
    ,
    diagnostic_parameters
  ]

rhat_values <-
  gelman.diag(
    diagnostic_samples,
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[
    ,
    1
  ]

ess_values <-
  effectiveSize(
    diagnostic_samples
  )

convergence_results <- tibble(
  
  parameter =
    names(
      rhat_values
    ),
  
  rhat =
    as.numeric(
      rhat_values
    ),
  
  effective_sample_size =
    as.numeric(
      ess_values[
        names(
          rhat_values
        )
      ]
    )
)

write_csv(
  convergence_results,
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_convergence.csv"
  )
)

# ============================================================
# PART L — RECOVERY METRICS
# ============================================================

max_rhat <- max(
  convergence_results$rhat
)

min_ess <- min(
  convergence_results$effective_sample_size
)

core_coverage <- mean(
  core_summary$covered_95
)

core_parameters_covered <- sum(
  core_summary$covered_95
)

max_abs_standardized_error <- max(
  abs(
    core_summary$standardized_error
  )
)

critical_parameters <- c(
  "phi",
  "sigma_sensor",
  "sigma_day"
)

critical_parameters_covered <- all(
  core_summary$covered_95[
    core_summary$parameter %in%
      critical_parameters
  ]
)

sensor_correlation <- cor(
  sensor_recovery$true_effect,
  sensor_recovery$posterior_mean
)

sensor_coverage <- mean(
  sensor_recovery$covered_95
)

day_correlation <- cor(
  day_recovery$true_effect,
  day_recovery$posterior_mean
)

day_coverage <- mean(
  day_recovery$covered_95
)

recovery_metrics <- tibble(
  
  metric = c(
    "maximum_rhat",
    "minimum_ess",
    "core_95_coverage",
    "core_parameters_covered",
    "maximum_absolute_standardized_error",
    "critical_variance_parameters_all_covered",
    "sensor_effect_correlation",
    "sensor_effect_95_coverage",
    "day_effect_correlation",
    "day_effect_95_coverage"
  ),
  
  value = c(
    max_rhat,
    min_ess,
    core_coverage,
    core_parameters_covered,
    max_abs_standardized_error,
    as.numeric(
      critical_parameters_covered
    ),
    sensor_correlation,
    sensor_coverage,
    day_correlation,
    day_coverage
  )
)

write_csv(
  recovery_metrics,
  here(
    "outputs",
    "tables",
    "m4_parameter_recovery_confirmation_metrics.csv"
  )
)

# ============================================================
# PART M — PRE-SPECIFIED ACCEPTANCE CRITERIA
# ============================================================

convergence_pass <-
  max_rhat < 1.05 &
  min_ess >= 400

# With one simulation and 17 core parameters, requiring
# 100% nominal 95% coverage would be unnecessarily strict.
# At least 15/17 corresponds to >= 88.2%.

core_recovery_pass <-
  core_parameters_covered >= 15 &
  max_abs_standardized_error <= 2.5 &
  critical_parameters_covered

latent_recovery_pass <-
  sensor_correlation >= 0.95 &
  sensor_coverage >= 0.90 &
  day_correlation >= 0.90 &
  day_coverage >= 0.90

overall_pass <-
  convergence_pass &
  core_recovery_pass &
  latent_recovery_pass

# ============================================================
# PART N — FIGURES
# ============================================================

core_plot_data <- core_summary %>%
  mutate(
    parameter =
      factor(
        parameter,
        levels =
          rev(
            parameter
          )
      )
  )

p_core <- ggplot(
  core_plot_data,
  aes(
    y = parameter,
    x = posterior_mean
  )
) +
  geom_errorbar(
    aes(
      xmin = q025,
      xmax = q975
    ),
    orientation = "y",
    width = 0.15
  ) +
  geom_point(
    size = 2
  ) +
  geom_point(
    aes(
      x = true_value
    ),
    shape = 4,
    size = 3
  ) +
  labs(
    title =
      "M4 Full-Design Parameter Recovery",
    subtitle =
      "Cross = true simulated value; point and interval = posterior mean and 95% credible interval",
    x =
      "Parameter value",
    y =
      NULL
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_parameter_recovery_confirmation_core.png"
  ),
  p_core,
  width = 9,
  height = 10,
  dpi = 300
)

p_sensor <- ggplot(
  sensor_recovery,
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
  geom_point() +
  labs(
    title =
      "M4 Sensor Random-Effect Recovery",
    x =
      "True sensor effect",
    y =
      "Posterior mean"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_parameter_recovery_confirmation_sensor.png"
  ),
  p_sensor,
  width = 7,
  height = 6,
  dpi = 300
)

p_day <- ggplot(
  day_recovery,
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
  geom_point(
    alpha = 0.6
  ) +
  labs(
    title =
      "M4 Day Random-Effect Recovery",
    x =
      "True day effect",
    y =
      "Posterior mean"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "m4_parameter_recovery_confirmation_day.png"
  ),
  p_day,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# PART O — QA
# ============================================================

cat(
  "\n========== FINAL M4 PARAMETER RECOVERY — CONFIRMATION ==========\n"
)

cat(
  "\nSynthetic data:\n"
)

cat(
  "  Observations:",
  N,
  "\n"
)

cat(
  "  Sensors:",
  J,
  "\n"
)

cat(
  "  Dates:",
  D,
  "\n"
)

cat(
  "  Simulated mean count:",
  mean(y_sim),
  "\n"
)

cat(
  "  Simulated maximum count:",
  max(y_sim),
  "\n"
)

cat(
  "\nCore parameter recovery:\n"
)

print(
  core_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nRecovery metrics:\n"
)

print(
  recovery_metrics,
  n = Inf,
  width = Inf
)

cat(
  "\nConvergence:\n"
)

cat(
  "  Maximum R-hat:",
  max_rhat,
  "\n"
)

cat(
  "  Minimum ESS:",
  min_ess,
  "\n"
)

cat(
  "\nAcceptance criteria:\n"
)

cat(
  "  Convergence PASS:",
  convergence_pass,
  "\n"
)

cat(
  "  Core recovery PASS:",
  core_recovery_pass,
  "\n"
)

cat(
  "  Latent-effect recovery PASS:",
  latent_recovery_pass,
  "\n"
)

cat(
  "\nOVERALL M4 PARAMETER RECOVERY PASS:",
  overall_pass,
  "\n"
)

cat(
  "=================================================\n"
)

message(
  "Final M4 confirmation parameter-recovery experiment completed."
)