# ============================================================
# 15_mcmc_diagnostics.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Consolidate MCMC convergence diagnostics for M1, M2 and M3
# and produce trace/autocorrelation diagnostics suitable for
# the final report.
#
# This is a light post-processing step and should run locally.
#
# Outputs:
#   outputs/tables/mcmc_diagnostics_summary.csv
#   outputs/tables/mcmc_selected_autocorrelation.csv
#   outputs/figures/mcmc_trace_m1.png
#   outputs/figures/mcmc_trace_m2.png
#   outputs/figures/mcmc_trace_m3.png
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(coda)
library(ggplot2)

# ------------------------------------------------------------
# 1. Validate required files
# ------------------------------------------------------------

required_files <- c(
  
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
  ),
  
  here(
    "outputs",
    "tables",
    "m1_poisson_convergence.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m2_negbin_convergence.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_convergence.csv"
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
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here(
    "outputs",
    "tables"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PART A — LOAD SAMPLES
# ============================================================

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

stopifnot(
  is.mcmc.list(m1_samples),
  is.mcmc.list(m2_samples),
  is.mcmc.list(m3_samples)
)

stopifnot(
  length(m1_samples) == 4,
  length(m2_samples) == 4,
  length(m3_samples) == 4
)

# ============================================================
# PART B — CONSOLIDATED RHAT / ESS
# ============================================================

read_convergence <- function(
    filename,
    model_name
) {
  
  read_csv(
    filename,
    show_col_types = FALSE
  ) %>%
    mutate(
      model =
        model_name,
      .before = 1
    )
}

m1_convergence <- read_convergence(
  here(
    "outputs",
    "tables",
    "m1_poisson_convergence.csv"
  ),
  "M1 Poisson"
)

m2_convergence <- read_convergence(
  here(
    "outputs",
    "tables",
    "m2_negbin_convergence.csv"
  ),
  "M2 Negative Binomial"
)

m3_convergence <- read_convergence(
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_convergence.csv"
  ),
  "M3 Hierarchical Negative Binomial"
)

all_convergence <- bind_rows(
  m1_convergence,
  m2_convergence,
  m3_convergence
)

mcmc_summary <- all_convergence %>%
  group_by(
    model
  ) %>%
  summarise(
    
    parameters_checked =
      n(),
    
    max_rhat =
      max(
        rhat
      ),
    
    min_ess =
      min(
        effective_sample_size
      ),
    
    median_ess =
      median(
        effective_sample_size
      ),
    
    all_rhat_below_1_05 =
      all(
        rhat < 1.05
      ),
    
    all_ess_at_least_400 =
      all(
        effective_sample_size >= 400
      ),
    
    convergence_pass =
      all_rhat_below_1_05 &
      all_ess_at_least_400,
    
    .groups = "drop"
  )

write_csv(
  mcmc_summary,
  here(
    "outputs",
    "tables",
    "mcmc_diagnostics_summary.csv"
  )
)

# ============================================================
# PART C — TRACE PLOT HELPERS
# ============================================================

extract_trace <- function(
    samples,
    model_name,
    parameters
) {
  
  available <-
    colnames(
      samples[[1]]
    )
  
  stopifnot(
    all(
      parameters %in%
        available
    )
  )
  
  result <- list()
  
  counter <- 1L
  
  for (chain_id in seq_along(samples)) {
    
    chain_matrix <-
      as.matrix(
        samples[[chain_id]]
      )
    
    for (parameter_name in parameters) {
      
      result[[counter]] <-
        tibble(
          
          model =
            model_name,
          
          chain =
            factor(
              chain_id
            ),
          
          iteration =
            seq_len(
              nrow(
                chain_matrix
              )
            ),
          
          parameter =
            parameter_name,
          
          value =
            chain_matrix[
              ,
              parameter_name
            ]
        )
      
      counter <-
        counter +
        1L
    }
  }
  
  bind_rows(
    result
  )
}

# ------------------------------------------------------------
# Representative parameters selected to show:
#
# alpha      = intercept
# beta_temp  = temporal/weather effect
# beta_business = sensor-level built-environment effect
# phi        = NB dispersion
# sigma_sensor = hierarchical variation
# ------------------------------------------------------------

m1_trace_parameters <- c(
  "alpha",
  "beta_temp",
  "beta_business",
  "beta_coverage"
)

m2_trace_parameters <- c(
  "alpha",
  "beta_temp",
  "beta_business",
  "phi"
)

m3_trace_parameters <- c(
  "alpha",
  "beta_temp",
  "beta_business",
  "phi",
  "sigma_sensor"
)

m1_trace <- extract_trace(
  m1_samples,
  "M1 Poisson",
  m1_trace_parameters
)

m2_trace <- extract_trace(
  m2_samples,
  "M2 Negative Binomial",
  m2_trace_parameters
)

m3_trace <- extract_trace(
  m3_samples,
  "M3 Hierarchical Negative Binomial",
  m3_trace_parameters
)

# ============================================================
# PART D — TRACE FIGURES
# ============================================================

make_trace_plot <- function(
    trace_data,
    title_text
) {
  
  ggplot(
    trace_data,
    aes(
      x = iteration,
      y = value,
      group = chain,
      colour = chain
    )
  ) +
    geom_line(
      linewidth = 0.25,
      alpha = 0.75
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free_y",
      ncol = 1
    ) +
    labs(
      title =
        title_text,
      x =
        "Retained MCMC iteration",
      y =
        "Parameter value",
      colour =
        "Chain"
    ) +
    theme_minimal() +
    theme(
      legend.position =
        "bottom"
    )
}

p_m1 <- make_trace_plot(
  m1_trace,
  "M1 Poisson: Representative MCMC Trace Plots"
)

p_m2 <- make_trace_plot(
  m2_trace,
  "M2 Negative Binomial: Representative MCMC Trace Plots"
)

p_m3 <- make_trace_plot(
  m3_trace,
  "M3 Hierarchical Negative Binomial: Representative MCMC Trace Plots"
)

ggsave(
  here(
    "outputs",
    "figures",
    "mcmc_trace_m1.png"
  ),
  p_m1,
  width = 9,
  height = 9,
  dpi = 300
)

ggsave(
  here(
    "outputs",
    "figures",
    "mcmc_trace_m2.png"
  ),
  p_m2,
  width = 9,
  height = 9,
  dpi = 300
)

ggsave(
  here(
    "outputs",
    "figures",
    "mcmc_trace_m3.png"
  ),
  p_m3,
  width = 9,
  height = 10,
  dpi = 300
)

# ============================================================
# PART E — AUTOCORRELATION
# ============================================================

chain_autocorrelation <- function(
    x,
    lag_value
) {
  
  acf_object <- acf(
    x,
    lag.max = lag_value,
    plot = FALSE,
    demean = TRUE
  )
  
  as.numeric(
    acf_object$acf[
      lag_value + 1
    ]
  )
}

calculate_autocorrelation <- function(
    samples,
    model_name,
    parameters
) {
  
  rows <- list()
  
  counter <- 1L
  
  for (parameter_name in parameters) {
    
    for (chain_id in seq_along(samples)) {
      
      x <-
        as.matrix(
          samples[[chain_id]]
        )[
          ,
          parameter_name
        ]
      
      rows[[counter]] <-
        tibble(
          
          model =
            model_name,
          
          parameter =
            parameter_name,
          
          chain =
            chain_id,
          
          autocorrelation_lag1 =
            chain_autocorrelation(
              x,
              1
            ),
          
          autocorrelation_lag10 =
            chain_autocorrelation(
              x,
              10
            )
        )
      
      counter <-
        counter +
        1L
    }
  }
  
  bind_rows(
    rows
  )
}

autocorrelation_results <- bind_rows(
  
  calculate_autocorrelation(
    m1_samples,
    "M1 Poisson",
    m1_trace_parameters
  ),
  
  calculate_autocorrelation(
    m2_samples,
    "M2 Negative Binomial",
    m2_trace_parameters
  ),
  
  calculate_autocorrelation(
    m3_samples,
    "M3 Hierarchical Negative Binomial",
    m3_trace_parameters
  )
) %>%
  
  group_by(
    model,
    parameter
  ) %>%
  
  summarise(
    
    mean_abs_autocorrelation_lag1 =
      mean(
        abs(
          autocorrelation_lag1
        )
      ),
    
    max_abs_autocorrelation_lag1 =
      max(
        abs(
          autocorrelation_lag1
        )
      ),
    
    mean_abs_autocorrelation_lag10 =
      mean(
        abs(
          autocorrelation_lag10
        )
      ),
    
    max_abs_autocorrelation_lag10 =
      max(
        abs(
          autocorrelation_lag10
        )
      ),
    
    .groups = "drop"
  )

write_csv(
  autocorrelation_results,
  here(
    "outputs",
    "tables",
    "mcmc_selected_autocorrelation.csv"
  )
)

# ============================================================
# PART F — QA OUTPUT
# ============================================================

cat(
  "\n========== STEP 15 MCMC DIAGNOSTICS ==========\n"
)

cat(
  "\nModel-level convergence summary:\n"
)

print(
  mcmc_summary,
  width = Inf
)

cat(
  "\nSelected-parameter autocorrelation:\n"
)

print(
  autocorrelation_results,
  width = Inf
)

cat(
  "\nAll three models convergence pass:",
  all(
    mcmc_summary$convergence_pass
  ),
  "\n"
)

cat(
  "==============================================\n"
)

message(
  "Step 15 MCMC diagnostics completed."
)