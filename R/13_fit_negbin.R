# ============================================================
# 13_fit_negbin.R
# Melbourne Pedestrian Bayesian Project
#
# M2: Bayesian Negative Binomial regression
#
# Computational strategy:
#   - 4 independent posterior chains
#   - chains executed in parallel on 4 CPU cores
#   - separate shorter 2-chain JAGS run for DIC
#
# Outputs:
#   outputs/mcmc/m2_negbin_samples.rds
#   outputs/tables/m2_negbin_posterior_summary.csv
#   outputs/tables/m2_negbin_convergence.csv
#   outputs/tables/m2_negbin_fit_metrics.csv
#   outputs/tables/m2_negbin_dic.txt
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(rjags)
library(coda)
library(parallel)

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
  sum(is.na(dat)) == 0,
  all(dat$daily_ped_count >= 0)
)

N <- nrow(dat)

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

model_file <- here(
  "models",
  "jags",
  "m2_negbin.jags"
)

stopifnot(
  file.exists(model_file)
)

# ============================================================
# PART B — INITIAL VALUES
# ============================================================

make_inits <- function(
    seed,
    alpha_init,
    phi_init
) {
  
  set.seed(seed)
  
  list(
    
    alpha =
      alpha_init,
    
    beta_holiday =
      rnorm(1, 0, 0.05),
    
    beta_temp =
      rnorm(1, 0, 0.05),
    
    beta_rain =
      rnorm(1, 0, 0.05),
    
    beta_business =
      rnorm(1, 0, 0.05),
    
    beta_population =
      rnorm(1, 0, 0.05),
    
    beta_coverage =
      rnorm(1, 0, 0.05),
    
    beta_height =
      rnorm(1, 0, 0.05),
    
    log_phi =
      log(phi_init),
    
    .RNG.name =
      "base::Mersenne-Twister",
    
    .RNG.seed =
      seed
  )
}

inits <- list(
  
  make_inits(
    51011,
    log(7000),
    5
  ),
  
  make_inits(
    52022,
    log(9000),
    10
  ),
  
  make_inits(
    53033,
    log(12000),
    20
  ),
  
  make_inits(
    54044,
    log(15000),
    40
  )
)

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
  "phi"
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
  "beta_height",
  "phi"
)

# ============================================================
# PART C — PARALLEL POSTERIOR CHAINS
# ============================================================

n_chains <- 4

n_adapt <- 1000

n_burn <- 1500

n_post <- 4000

cat(
  "\n[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting 4 parallel M2 chains on 4 CPU cores...\n",
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
    "inits",
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
    
    model <- jags.model(
      file = model_file,
      data = jags_data,
      inits = inits[[chain_id]],
      n.chains = 1,
      n.adapt = n_adapt,
      quiet = TRUE
    )
    
    update(
      model,
      n.iter = n_burn,
      progress.bar = "none"
    )
    
    samples <- coda.samples(
      model,
      variable.names =
        parameters_to_monitor,
      n.iter = n_post,
      thin = 1,
      progress.bar = "none"
    )
    
    samples[[1]]
  }
)

stopCluster(cl)

m2_samples <- mcmc.list(
  chain_results
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Parallel posterior sampling completed.\n",
  sep = ""
)

saveRDS(
  m2_samples,
  here(
    "outputs",
    "mcmc",
    "m2_negbin_samples.rds"
  )
)

# ============================================================
# PART D — POSTERIOR SUMMARY
# ============================================================

posterior_matrix <- as.matrix(
  m2_samples
)

stopifnot(
  all(
    core_parameters %in%
      colnames(
        posterior_matrix
      )
  )
)

posterior_summary <- lapply(
  core_parameters,
  function(parameter_name) {
    
    x <-
      posterior_matrix[
        ,
        parameter_name
      ]
    
    exponentiate <-
      parameter_name != "phi"
    
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
        if (exponentiate) {
          exp(
            median(x)
          )
        } else {
          NA_real_
        },
      
      exp_q025 =
        if (exponentiate) {
          exp(
            quantile(
              x,
              0.025
            )
          )
        } else {
          NA_real_
        },
      
      exp_q975 =
        if (exponentiate) {
          exp(
            quantile(
              x,
              0.975
            )
          )
        } else {
          NA_real_
        }
    )
  }
) %>%
  bind_rows()

write_csv(
  posterior_summary,
  here(
    "outputs",
    "tables",
    "m2_negbin_posterior_summary.csv"
  )
)

# ============================================================
# PART E — CONVERGENCE
# ============================================================

core_samples <-
  m2_samples[
    ,
    core_parameters
  ]

rhat_values <-
  gelman.diag(
    core_samples,
    autoburnin = FALSE,
    multivariate = FALSE
  )$psrf[
    ,
    1
  ]

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
    "m2_negbin_convergence.csv"
  )
)

# ============================================================
# PART F — DIC
# ============================================================

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting separate 2-chain DIC run...\n",
  sep = ""
)

dic_model <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits[1:2],
  n.chains = 2,
  n.adapt = 1000,
  quiet = TRUE
)

update(
  dic_model,
  n.iter = 1500,
  progress.bar = "none"
)

m2_dic <- dic.samples(
  dic_model,
  n.iter = 2000,
  thin = 1,
  type = "pD",
  progress.bar = "none"
)

mean_deviance <-
  sum(
    m2_dic$deviance
  )

effective_parameters <-
  sum(
    m2_dic$penalty
  )

dic_value <-
  mean_deviance +
  effective_parameters

fit_metrics <- tibble(
  
  model =
    "M2 Negative Binomial",
  
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
    "m2_negbin_fit_metrics.csv"
  )
)

writeLines(
  capture.output(
    print(
      m2_dic
    )
  ),
  here(
    "outputs",
    "tables",
    "m2_negbin_dic.txt"
  )
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] DIC run completed.\n",
  sep = ""
)

# ============================================================
# PART G — FORMAL QA
# ============================================================

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

cat(
  "\n========== M2 NEGATIVE BINOMIAL QA ==========\n"
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
  "Posterior chains: 4\n"
)

cat(
  "Parallel CPU workers: 4\n"
)

cat(
  "Adaptation per posterior chain: 1000\n"
)

cat(
  "Burn-in per posterior chain: 1500\n"
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
  "\nOVERALL M2 MCMC PASS:",
  mcmc_pass,
  "\n"
)

cat(
  "=============================================\n"
)

message(
  "M2 Bayesian Negative Binomial fit completed."
)