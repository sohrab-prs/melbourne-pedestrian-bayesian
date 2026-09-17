# ============================================================
# 14_fit_hierarchical_negbin.R
# Melbourne Pedestrian Bayesian Project
#
# M3:
# Hierarchical Negative Binomial model with sensor-specific
# random intercepts.
#
# Hierarchical centering was validated through the Step 11
# parameter-recovery experiment.
#
# Computational strategy:
#   - 4 independent posterior chains on 4 CPU workers
#   - DIC job runs concurrently on a fifth worker
#   - no thinning
#
# Outputs:
#   outputs/mcmc/m3_hierarchical_negbin_samples.rds
#   outputs/tables/m3_hierarchical_negbin_posterior_summary.csv
#   outputs/tables/m3_hierarchical_negbin_convergence.csv
#   outputs/tables/m3_hierarchical_negbin_sensor_effects.csv
#   outputs/tables/m3_hierarchical_negbin_fit_metrics.csv
#   outputs/tables/m3_hierarchical_negbin_dic.txt
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(rjags)
library(coda)
library(parallel)

# ------------------------------------------------------------
# Output directories
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

J <- 53L

# ------------------------------------------------------------
# Sensor-level spatial predictor table
# ------------------------------------------------------------

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
  anyDuplicated(sensor_level$sensor_index) == 0,
  identical(
    sensor_level$sensor_index,
    seq_len(J)
  )
)

# ------------------------------------------------------------
# JAGS data
# ------------------------------------------------------------

jags_data <- list(
  
  N = N,
  
  J = J,
  
  y =
    as.integer(
      dat$daily_ped_count
    ),
  
  sensor =
    as.integer(
      dat$sensor_index
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
  "m3_hierarchical_negbin.jags"
)

stopifnot(
  file.exists(
    model_file
  )
)

# ============================================================
# PART B — INITIAL VALUES
# ============================================================

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

posterior_inits <- list(
  
  make_inits(
    71011,
    log(7000),
    5,
    0.20
  ),
  
  make_inits(
    72022,
    log(9000),
    10,
    0.35
  ),
  
  make_inits(
    73033,
    log(12000),
    20,
    0.50
  ),
  
  make_inits(
    74044,
    log(15000),
    40,
    0.70
  )
)

dic_inits <- list(
  
  make_inits(
    75055,
    log(9000),
    10,
    0.30
  ),
  
  make_inits(
    76066,
    log(13000),
    25,
    0.60
  )
)

# ============================================================
# PART C — PARAMETERS
# ============================================================

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
  "phi",
  "sigma_sensor"
)

sensor_parameters <- paste0(
  "u_sensor[",
  seq_len(J),
  "]"
)

parameters_to_monitor <- c(
  core_parameters,
  "u_sensor"
)

diagnostic_parameters <- c(
  core_parameters,
  sensor_parameters
)

# ============================================================
# PART D — COMPUTATION SETTINGS
# ============================================================

# Lean but defensible settings.
# If convergence criteria fail, only then increase iterations.

n_adapt <- 1500

n_burn <- 2000

n_post <- 3500

dic_adapt <- 1000

dic_burn <- 1500

dic_iter <- 1500

# ============================================================
# PART E — PARALLEL POSTERIOR + DIC
# ============================================================

cat(
  "\n[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting M3: 4 posterior chains + 1 concurrent DIC job...\n",
  sep = ""
)

n_workers <- 5

cl <- makePSOCKcluster(
  n_workers,
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
    "posterior_inits",
    "dic_inits",
    "parameters_to_monitor",
    "n_adapt",
    "n_burn",
    "n_post",
    "dic_adapt",
    "dic_burn",
    "dic_iter"
  ),
  envir = environment()
)

jobs <- parLapply(
  cl,
  seq_len(n_workers),
  function(job_id) {
    
    if (job_id <= 4) {
      
      cat(
        "[",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] Posterior chain ",
        job_id,
        " starting\n",
        sep = ""
      )
      
      model <- jags.model(
        file = model_file,
        data = jags_data,
        inits = posterior_inits[[job_id]],
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
      
      cat(
        "[",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] Posterior chain ",
        job_id,
        " completed\n",
        sep = ""
      )
      
      return(
        list(
          type = "posterior",
          chain_id = job_id,
          samples = samples[[1]]
        )
      )
    }
    
    cat(
      "[",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "] Concurrent DIC job starting\n",
      sep = ""
    )
    
    dic_model <- jags.model(
      file = model_file,
      data = jags_data,
      inits = dic_inits,
      n.chains = 2,
      n.adapt = dic_adapt,
      quiet = TRUE
    )
    
    update(
      dic_model,
      n.iter = dic_burn,
      progress.bar = "none"
    )
    
    dic_object <- dic.samples(
      dic_model,
      n.iter = dic_iter,
      thin = 1,
      type = "pD",
      progress.bar = "none"
    )
    
    mean_deviance <-
      sum(
        dic_object$deviance
      )
    
    effective_parameters <-
      sum(
        dic_object$penalty
      )
    
    dic_value <-
      mean_deviance +
      effective_parameters
    
    cat(
      "[",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "] Concurrent DIC job completed\n",
      sep = ""
    )
    
    list(
      type = "dic",
      mean_deviance =
        mean_deviance,
      effective_parameters =
        effective_parameters,
      DIC =
        dic_value,
      dic_text =
        capture.output(
          print(
            dic_object
          )
        )
    )
  }
)

stopCluster(cl)

# ============================================================
# PART F — COLLECT POSTERIOR CHAINS
# ============================================================

posterior_jobs <- jobs[
  vapply(
    jobs,
    function(x) {
      identical(
        x$type,
        "posterior"
      )
    },
    logical(1)
  )
]

stopifnot(
  length(
    posterior_jobs
  ) == 4
)

posterior_jobs <- posterior_jobs[
  order(
    vapply(
      posterior_jobs,
      function(x) {
        x$chain_id
      },
      integer(1)
    )
  )
]

m3_samples <- mcmc.list(
  lapply(
    posterior_jobs,
    function(x) {
      x$samples
    }
  )
)

saveRDS(
  m3_samples,
  here(
    "outputs",
    "mcmc",
    "m3_hierarchical_negbin_samples.rds"
  )
)

cat(
  "[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] All M3 computational jobs collected.\n",
  sep = ""
)

# ============================================================
# PART G — POSTERIOR SUMMARY
# ============================================================

posterior_matrix <- as.matrix(
  m3_samples
)

stopifnot(
  all(
    diagnostic_parameters %in%
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
      !parameter_name %in%
      c(
        "phi",
        "sigma_sensor"
      )
    
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
    "m3_hierarchical_negbin_posterior_summary.csv"
  )
)

# ============================================================
# PART H — SENSOR RANDOM EFFECTS
# ============================================================

sensor_effects <- lapply(
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
      sensor_index =
        j,
      
      location_id =
        sensor_level$location_id[
          j
        ],
      
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
        )
    )
  }
) %>%
  bind_rows()

write_csv(
  sensor_effects,
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_sensor_effects.csv"
  )
)

# ============================================================
# PART I — CONVERGENCE
# ============================================================

diagnostic_samples <-
  m3_samples[
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
    "m3_hierarchical_negbin_convergence.csv"
  )
)

# ============================================================
# PART J — DIC
# ============================================================

dic_jobs <- jobs[
  vapply(
    jobs,
    function(x) {
      identical(
        x$type,
        "dic"
      )
    },
    logical(1)
  )
]

stopifnot(
  length(
    dic_jobs
  ) == 1
)

dic_result <-
  dic_jobs[[1]]

fit_metrics <- tibble(
  
  model =
    "M3 Hierarchical Negative Binomial",
  
  observations =
    N,
  
  sensors =
    J,
  
  mean_deviance =
    dic_result$mean_deviance,
  
  effective_parameters_pD =
    dic_result$effective_parameters,
  
  DIC =
    dic_result$DIC
)

write_csv(
  fit_metrics,
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_fit_metrics.csv"
  )
)

writeLines(
  dic_result$dic_text,
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_dic.txt"
  )
)

# ============================================================
# PART K — FORMAL QA
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

# ============================================================
# PART L — OUTPUT
# ============================================================

cat(
  "\n========== M3 HIERARCHICAL NB QA ==========\n"
)

cat(
  "Observations:",
  N,
  "\n"
)

cat(
  "Sensors:",
  J,
  "\n"
)

cat(
  "Posterior chains: 4\n"
)

cat(
  "Concurrent compute workers: 5\n"
)

cat(
  "Posterior adaptation per chain:",
  n_adapt,
  "\n"
)

cat(
  "Posterior burn-in per chain:",
  n_burn,
  "\n"
)

cat(
  "Retained posterior iterations per chain:",
  n_post,
  "\n"
)

cat(
  "Thinning: 1\n"
)

cat(
  "\nCore posterior summary:\n"
)

print(
  posterior_summary,
  width = Inf
)

cat(
  "\nCore + sensor-effect convergence diagnostics:\n"
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
  "\nOVERALL M3 MCMC PASS:",
  mcmc_pass,
  "\n"
)

cat(
  "===========================================\n"
)

message(
  "M3 Hierarchical Negative Binomial fit completed."
)