# ============================================================
# 17b_fit_sensor_day_negbin.R
# Melbourne Pedestrian Bayesian Project
#
# M4:
# Hierarchical Negative Binomial with:
#   - sensor random intercepts
#   - shared calendar-day random effects
#
# Motivated by M3 posterior predictive checking:
# M3 reproduced sensor-level heterogeneity but underestimated
# variation in date-level mean pedestrian activity.
#
# Computational strategy:
#   - 4 independent posterior chains
#   - 1 concurrent DIC job
#   - 5 workstation processes
#   - no thinning
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
  sum(is.na(dat)) == 0
)

N <- nrow(dat)

J <- 53L

# ------------------------------------------------------------
# Date indexing
# ------------------------------------------------------------

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

stopifnot(
  D == 364
)

dat <- dat %>%
  left_join(
    date_lookup,
    by = "date"
  )

stopifnot(
  sum(
    is.na(
      dat$day_index
    )
  ) == 0
)

# ------------------------------------------------------------
# Sensor-level predictors
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
  
  D = D,
  
  y =
    as.integer(
      dat$daily_ped_count
    ),
  
  sensor =
    as.integer(
      dat$sensor_index
    ),
  
  day =
    as.integer(
      dat$day_index
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
  "m4_sensor_day_negbin.jags"
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
    sigma_sensor_init,
    sigma_day_init
) {
  
  set.seed(seed)
  
  beta_spatial_init <- rnorm(
    4,
    0,
    0.10
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

posterior_inits <- list(
  
  make_inits(
    81011,
    log(7000),
    8,
    0.25,
    0.10
  ),
  
  make_inits(
    82022,
    log(9000),
    12,
    0.40,
    0.20
  ),
  
  make_inits(
    83033,
    log(12000),
    20,
    0.60,
    0.30
  ),
  
  make_inits(
    84044,
    log(15000),
    30,
    0.80,
    0.40
  )
)

dic_inits <- list(
  
  make_inits(
    85055,
    log(9000),
    12,
    0.35,
    0.15
  ),
  
  make_inits(
    86066,
    log(13000),
    25,
    0.65,
    0.35
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
  "sigma_sensor",
  "sigma_day"
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

# ============================================================
# PART D — COMPUTE SETTINGS
# ============================================================

# Efficient first-pass settings.
# Increase only if formal convergence criteria fail.

n_adapt <- 1000

n_burn <- 1500

n_post <- 7000

dic_adapt <- 800

dic_burn <- 1000

dic_iter <- 1200

# ============================================================
# PART E — PARALLEL FIT
# ============================================================

cat(
  "\n[",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  "] Starting M4: 4 posterior chains + concurrent DIC...\n",
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
        "] M4 posterior chain ",
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
        "] M4 posterior chain ",
        job_id,
        " completed\n",
        sep = ""
      )
      
      return(
        list(
          type =
            "posterior",
          chain_id =
            job_id,
          samples =
            samples[[1]]
        )
      )
    }
    
    cat(
      "[",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "] M4 concurrent DIC job starting\n",
      sep = ""
    )
    
    model <- jags.model(
      file = model_file,
      data = jags_data,
      inits = dic_inits,
      n.chains = 2,
      n.adapt = dic_adapt,
      quiet = TRUE
    )
    
    update(
      model,
      n.iter = dic_burn,
      progress.bar = "none"
    )
    
    dic_object <- dic.samples(
      model,
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
      "] M4 concurrent DIC job completed\n",
      sep = ""
    )
    
    list(
      type =
        "dic",
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
# PART F — COMBINE POSTERIOR CHAINS
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

stopifnot(
  length(
    posterior_jobs
  ) == 4
)

m4_samples <- mcmc.list(
  lapply(
    posterior_jobs,
    function(x) {
      x$samples
    }
  )
)

saveRDS(
  m4_samples,
  here(
    "outputs",
    "mcmc",
    "m4_sensor_day_negbin_samples.rds"
  )
)

posterior_matrix <- as.matrix(
  m4_samples
)

# ============================================================
# PART G — CORE POSTERIOR SUMMARY
# ============================================================

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
        "sigma_sensor",
        "sigma_day"
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
    "m4_sensor_day_negbin_posterior_summary.csv"
  )
)

# ============================================================
# PART H — RANDOM EFFECT SUMMARIES
# ============================================================

sensor_effects <- lapply(
  seq_len(J),
  function(j) {
    
    x <-
      posterior_matrix[
        ,
        paste0(
          "u_sensor[",
          j,
          "]"
        )
      ]
    
    tibble(
      sensor_index =
        j,
      location_id =
        sensor_level$location_id[j],
      posterior_mean =
        mean(x),
      q025 =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      median =
        median(x),
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
    "m4_sensor_day_negbin_sensor_effects.csv"
  )
)

day_effects <- lapply(
  seq_len(D),
  function(d) {
    
    x <-
      posterior_matrix[
        ,
        paste0(
          "day_effect[",
          d,
          "]"
        )
      ]
    
    tibble(
      day_index =
        d,
      date =
        date_lookup$date[d],
      posterior_mean =
        mean(x),
      q025 =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      median =
        median(x),
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
  day_effects,
  here(
    "outputs",
    "tables",
    "m4_sensor_day_negbin_day_effects.csv"
  )
)

# ============================================================
# PART I — CONVERGENCE
# ============================================================

diagnostic_samples <-
  m4_samples[
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

lowest_ess_parameters <- convergence_results %>%
  arrange(
    effective_sample_size
  ) %>%
  slice_head(
    n = 10
  )

write_csv(
  convergence_results,
  here(
    "outputs",
    "tables",
    "m4_sensor_day_negbin_convergence.csv"
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
    "M4 Sensor-Day Hierarchical Negative Binomial",
  
  observations =
    N,
  
  sensors =
    J,
  
  dates =
    D,
  
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
    "m4_sensor_day_negbin_fit_metrics.csv"
  )
)

writeLines(
  dic_result$dic_text,
  here(
    "outputs",
    "tables",
    "m4_sensor_day_negbin_dic.txt"
  )
)

# ============================================================
# PART K — QA
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
  "\n========== M4 SENSOR + DAY NB QA ==========\n"
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
  "Represented dates:",
  D,
  "\n"
)

cat(
  "Concurrent workers: 5\n"
)

cat(
  "\nCore posterior summary:\n"
)

print(
  posterior_summary,
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
  "\nLowest-ESS parameters:\n"
)

print(
  lowest_ess_parameters,
  n = 10,
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
  "\nOVERALL M4 MCMC PASS:",
  mcmc_pass,
  "\n"
)

cat(
  "===========================================\n"
)

message(
  "M4 sensor-day hierarchical Negative Binomial completed."
)