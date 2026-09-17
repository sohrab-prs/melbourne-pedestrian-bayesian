# ============================================================
# 12_fit_poisson.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Fit M1, the Bayesian Poisson regression, to the real
# Melbourne pedestrian sensor-day data using JAGS/MCMC.
#
# M1 is the benchmark count model.
#
# Inputs:
#   data/processed/modelling_dataset_selected_2025.csv
#   models/jags/m1_poisson.jags
#
# Outputs:
#   outputs/mcmc/m1_poisson_samples.rds
#   outputs/tables/m1_poisson_posterior_summary.csv
#   outputs/tables/m1_poisson_convergence.csv
#   outputs/tables/m1_poisson_fit_metrics.csv
#   outputs/tables/m1_poisson_dic.txt
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(rjags)
library(coda)

# ------------------------------------------------------------
# 1. Output directories
# ------------------------------------------------------------

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

# ============================================================
# PART A — LOAD AND VALIDATE DATA
# ============================================================

# ------------------------------------------------------------
# 2. Load selected modelling dataset
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
  nrow(dat) == 18671
)

stopifnot(
  n_distinct(
    dat$location_id
  ) == 53
)

stopifnot(
  n_distinct(
    dat$sensor_index
  ) == 53
)

stopifnot(
  sum(
    is.na(dat)
  ) == 0
)

stopifnot(
  all(
    dat$daily_ped_count >= 0
  )
)

N <- nrow(dat)

# ============================================================
# PART B — PREPARE JAGS DATA
# ============================================================

# ------------------------------------------------------------
# 3. JAGS input
# ------------------------------------------------------------

jags_data <- list(
  
  N = N,
  
  y =
    as.integer(
      dat$daily_ped_count
    ),
  
  dow =
    as.integer(
      dat$dow_index
    ),
  
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
  
  business =
    as.numeric(
      dat$z_business_count
    ),
  
  population =
    as.numeric(
      dat$z_resident_population
    ),
  
  coverage =
    as.numeric(
      dat$z_built_coverage
    ),
  
  height =
    as.numeric(
      dat$z_structure_height
    )
)

# ============================================================
# PART C — MCMC INITIAL VALUES
# ============================================================

# ------------------------------------------------------------
# 4. Reproducible dispersed chain initialisation
# ------------------------------------------------------------

make_inits <- function(
    seed,
    alpha_init
) {
  
  set.seed(seed)
  
  list(
    
    alpha =
      alpha_init,
    
    beta_holiday =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_temp =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_rain =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_business =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_population =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_coverage =
      rnorm(
        1,
        0,
        0.05
      ),
    
    beta_height =
      rnorm(
        1,
        0,
        0.05
      ),
    
    .RNG.name =
      "base::Mersenne-Twister",
    
    .RNG.seed =
      seed
  )
}

inits <- list(
  
  make_inits(
    seed = 11011,
    alpha_init = log(7000)
  ),
  
  make_inits(
    seed = 22022,
    alpha_init = log(9000)
  ),
  
  make_inits(
    seed = 33033,
    alpha_init = log(12000)
  ),
  
  make_inits(
    seed = 44044,
    alpha_init = log(15000)
  )
)

# ============================================================
# PART D — COMPILE AND RUN JAGS
# ============================================================

# ------------------------------------------------------------
# 5. Model file
# ------------------------------------------------------------

model_file <- here(
  "models",
  "jags",
  "m1_poisson.jags"
)

stopifnot(
  file.exists(
    model_file
  )
)

# ------------------------------------------------------------
# 6. Compile and adapt
# ------------------------------------------------------------

cat(
  "\n[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Compiling and adapting M1...\n",
  sep = ""
)

m1_model <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits,
  n.chains = 4,
  n.adapt = 1000,
  quiet = TRUE
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Adaptation completed. Starting burn-in...\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Burn-in
# ------------------------------------------------------------

update(
  m1_model,
  n.iter = 2000,
  progress.bar = "none"
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Burn-in completed. Starting posterior sampling...\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. Parameters to monitor
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
  "beta_height"
)

core_parameters <- c(
  "alpha",
  paste0(
    "beta_dow[",
    2:7,
    "]"
  ),
  "beta_holiday",
  "beta_temp",
  "beta_rain",
  "beta_business",
  "beta_population",
  "beta_coverage",
  "beta_height"
)

# ------------------------------------------------------------
# 9. Posterior sampling
# ------------------------------------------------------------

m1_samples <- coda.samples(
  m1_model,
  variable.names =
    parameters_to_monitor,
  n.iter = 4000,
  thin = 1,
  progress.bar = "none"
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Posterior sampling completed.\n",
  sep = ""
)

saveRDS(
  m1_samples,
  here(
    "outputs",
    "mcmc",
    "m1_poisson_samples.rds"
  )
)

# ============================================================
# PART E — POSTERIOR SUMMARY
# ============================================================

# ------------------------------------------------------------
# 10. Convert to matrix
# ------------------------------------------------------------

posterior_matrix <- as.matrix(
  m1_samples
)

stopifnot(
  all(
    core_parameters %in%
      colnames(
        posterior_matrix
      )
  )
)

# ------------------------------------------------------------
# 11. Posterior summaries
# ------------------------------------------------------------

posterior_summary <- lapply(
  core_parameters,
  function(parameter_name) {
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    tibble(
      parameter =
        parameter_name,
      
      posterior_mean =
        mean(x),
      
      posterior_sd =
        sd(x),
      
      q025 =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      
      median =
        unname(
          quantile(
            x,
            0.50
          )
        ),
      
      q975 =
        unname(
          quantile(
            x,
            0.975
          )
        ),
      
      exp_median =
        exp(
          median(x)
        ),
      
      exp_q025 =
        exp(
          quantile(
            x,
            0.025
          )
        ),
      
      exp_q975 =
        exp(
          quantile(
            x,
            0.975
          )
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  posterior_summary,
  here(
    "outputs",
    "tables",
    "m1_poisson_posterior_summary.csv"
  )
)

# ============================================================
# PART F — MCMC CONVERGENCE
# ============================================================

# ------------------------------------------------------------
# 12. Restrict convergence diagnostics to stochastic
#     parameters of substantive interest
# ------------------------------------------------------------

core_samples <-
  m1_samples[
    ,
    core_parameters
  ]

# ------------------------------------------------------------
# 13. R-hat
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
# 14. Effective sample size
# ------------------------------------------------------------

ess_values <-
  effectiveSize(
    core_samples
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
    "m1_poisson_convergence.csv"
  )
)

# ============================================================
# PART G — DIC
# ============================================================

# ------------------------------------------------------------
# 15. DIC using pD penalty
#
# DIC is required later for comparison with M2/M3.
# ------------------------------------------------------------

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting DIC sampling...\n",
  sep = ""
)

m1_dic <- dic.samples(
  m1_model,
  n.iter = 2000,
  thin = 1,
  type = "pD",
  progress.bar = "none"
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] DIC sampling completed.\n",
  sep = ""
)

mean_deviance <-
  sum(
    m1_dic$deviance
  )

effective_parameters <-
  sum(
    m1_dic$penalty
  )

dic_value <-
  mean_deviance +
  effective_parameters

fit_metrics <- tibble(
  model =
    "M1 Poisson",
  
  observations =
    N,
  
  mean_deviance =
    mean_deviance,
  
  effective_parameters_pD =
    effective_parameters,
  
  DIC =
    dic_value
)

write_csv(
  fit_metrics,
  here(
    "outputs",
    "tables",
    "m1_poisson_fit_metrics.csv"
  )
)

writeLines(
  capture.output(
    print(
      m1_dic
    )
  ),
  here(
    "outputs",
    "tables",
    "m1_poisson_dic.txt"
  )
)

# ============================================================
# PART H — FORMAL MCMC QA
# ============================================================

# ------------------------------------------------------------
# 16. Predefined convergence criteria
# ------------------------------------------------------------

max_rhat <-
  max(
    convergence_results$rhat
  )

min_ess <-
  min(
    convergence_results$effective_sample_size
  )

mcmc_pass <-
  max_rhat < 1.05 &
  min_ess >= 400

# ============================================================
# PART I — CONSOLE SUMMARY
# ============================================================

cat(
  "\n========== M1 POISSON QA ==========\n"
)

cat(
  "Observations:",
  N,
  "\n"
)

cat(
  "Sensors represented:",
  n_distinct(
    dat$sensor_index
  ),
  "\n"
)

cat(
  "Chains: 4\n"
)

cat(
  "Adaptation iterations per chain: 1000\n"
)

cat(
  "Burn-in iterations per chain: 2000\n"
)

cat(
  "Retained posterior iterations per chain: 4000\n"
)

cat(
  "Thinning: 1\n"
)

cat(
  "\nPosterior summary:\n"
)

print(
  posterior_summary,
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
  max_rhat,
  "\n"
)

cat(
  "Minimum effective sample size:",
  min_ess,
  "\n"
)

cat(
  "\nDIC metrics:\n"
)

print(
  fit_metrics,
  width = Inf
)

cat(
  "\nFormal convergence criteria:\n"
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
  "\nOVERALL M1 MCMC PASS:",
  mcmc_pass,
  "\n"
)

cat(
  "===================================\n"
)

message(
  "M1 Bayesian Poisson fit completed."
)