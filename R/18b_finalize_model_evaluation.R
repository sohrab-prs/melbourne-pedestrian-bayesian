# ============================================================
# 18b_finalize_model_evaluation.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Produce final convergence and DIC comparison tables after M4
# was accepted as the primary model.
#
# Also creates an M4 representative trace plot.
#
# No models are refitted.
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(coda)
library(ggplot2)

# ============================================================
# PART A — FILES
# ============================================================

required_files <- c(
  
  here("outputs", "tables", "m1_poisson_fit_metrics.csv"),
  here("outputs", "tables", "m2_negbin_fit_metrics.csv"),
  here("outputs", "tables", "m3_hierarchical_negbin_fit_metrics.csv"),
  here("outputs", "tables", "m4_sensor_day_negbin_fit_metrics.csv"),
  
  here("outputs", "tables", "m1_poisson_convergence.csv"),
  here("outputs", "tables", "m2_negbin_convergence.csv"),
  here("outputs", "tables", "m3_hierarchical_negbin_convergence.csv"),
  here("outputs", "tables", "m4_sensor_day_negbin_convergence.csv"),
  
  here("outputs", "mcmc", "m4_sensor_day_negbin_samples.rds")
)

stopifnot(
  all(
    file.exists(required_files)
  )
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

# ============================================================
# PART B — FINAL CONVERGENCE SUMMARY
# ============================================================

read_convergence <- function(
    path,
    model_name
) {
  
  read_csv(
    path,
    show_col_types = FALSE
  ) %>%
    mutate(
      model = model_name,
      .before = 1
    )
}

all_convergence <- bind_rows(
  
  read_convergence(
    here(
      "outputs",
      "tables",
      "m1_poisson_convergence.csv"
    ),
    "M1 Poisson"
  ),
  
  read_convergence(
    here(
      "outputs",
      "tables",
      "m2_negbin_convergence.csv"
    ),
    "M2 Negative Binomial"
  ),
  
  read_convergence(
    here(
      "outputs",
      "tables",
      "m3_hierarchical_negbin_convergence.csv"
    ),
    "M3 Hierarchical Negative Binomial"
  ),
  
  read_convergence(
    here(
      "outputs",
      "tables",
      "m4_sensor_day_negbin_convergence.csv"
    ),
    "M4 Sensor-Day Hierarchical Negative Binomial"
  )
)

convergence_summary <- all_convergence %>%
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
  convergence_summary,
  here(
    "outputs",
    "tables",
    "mcmc_diagnostics_summary_final.csv"
  )
)

# ============================================================
# PART C — FINAL DIC COMPARISON
# ============================================================

extract_fit <- function(
    path,
    model_name,
    model_order
) {
  
  x <- read_csv(
    path,
    show_col_types = FALSE
  )
  
  tibble(
    
    model =
      model_name,
    
    model_order =
      model_order,
    
    observations =
      x$observations[1],
    
    mean_deviance =
      x$mean_deviance[1],
    
    effective_parameters_pD =
      x$effective_parameters_pD[1],
    
    DIC =
      x$DIC[1]
  )
}

model_comparison <- bind_rows(
  
  extract_fit(
    here(
      "outputs",
      "tables",
      "m1_poisson_fit_metrics.csv"
    ),
    "M1 Poisson",
    1
  ),
  
  extract_fit(
    here(
      "outputs",
      "tables",
      "m2_negbin_fit_metrics.csv"
    ),
    "M2 Negative Binomial",
    2
  ),
  
  extract_fit(
    here(
      "outputs",
      "tables",
      "m3_hierarchical_negbin_fit_metrics.csv"
    ),
    "M3 Hierarchical Negative Binomial",
    3
  ),
  
  extract_fit(
    here(
      "outputs",
      "tables",
      "m4_sensor_day_negbin_fit_metrics.csv"
    ),
    "M4 Sensor-Day Hierarchical Negative Binomial",
    4
  )
) %>%
  arrange(
    model_order
  ) %>%
  mutate(
    
    improvement_from_previous =
      lag(DIC) -
      DIC,
    
    delta_DIC_from_best =
      DIC -
      min(DIC)
  )

best_model <-
  model_comparison$model[
    which.min(
      model_comparison$DIC
    )
  ]

model_comparison <- model_comparison %>%
  left_join(
    convergence_summary %>%
      select(
        model,
        max_rhat,
        min_ess,
        convergence_pass
      ),
    by = "model"
  )

write_csv(
  model_comparison,
  here(
    "outputs",
    "tables",
    "model_comparison_final.csv"
  )
)

# ============================================================
# PART D — FINAL DIC FIGURE
# ============================================================

plot_data <- model_comparison %>%
  mutate(
    model_short =
      factor(
        c(
          "M1 Poisson",
          "M2 NB",
          "M3 Sensor NB",
          "M4 Sensor + Day NB"
        ),
        levels = c(
          "M1 Poisson",
          "M2 NB",
          "M3 Sensor NB",
          "M4 Sensor + Day NB"
        )
      )
  )

p_dic <- ggplot(
  plot_data,
  aes(
    x = model_short,
    y = DIC
  )
) +
  geom_col() +
  geom_text(
    aes(
      label =
        format(
          round(DIC),
          big.mark = ",",
          scientific = FALSE
        )
    ),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_log10() +
  labs(
    title =
      "Final Bayesian Model Comparison",
    subtitle =
      "DIC shown on a logarithmic axis; lower values indicate better relative fit",
    x =
      NULL,
    y =
      "DIC (log scale)"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "model_comparison_dic_final.png"
  ),
  p_dic,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================
# PART E — M4 TRACE PLOT
# ============================================================

m4_samples <- readRDS(
  here(
    "outputs",
    "mcmc",
    "m4_sensor_day_negbin_samples.rds"
  )
)

stopifnot(
  is.mcmc.list(m4_samples),
  length(m4_samples) == 4
)

trace_parameters <- c(
  "alpha",
  "beta_temp",
  "beta_business",
  "phi",
  "sigma_sensor",
  "sigma_day"
)

stopifnot(
  all(
    trace_parameters %in%
      colnames(
        m4_samples[[1]]
      )
  )
)

trace_rows <- list()

counter <- 1L

for (
  chain_id in seq_along(
    m4_samples
  )
) {
  
  chain_matrix <-
    as.matrix(
      m4_samples[[chain_id]]
    )
  
  for (
    parameter_name in trace_parameters
  ) {
    
    trace_rows[[counter]] <-
      tibble(
        
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

trace_data <-
  bind_rows(
    trace_rows
  )

p_trace <- ggplot(
  trace_data,
  aes(
    x =
      iteration,
    y =
      value,
    group =
      chain,
    colour =
      chain
  )
) +
  geom_line(
    linewidth = 0.22,
    alpha = 0.7
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title =
      "M4 Representative MCMC Trace Plots",
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

ggsave(
  here(
    "outputs",
    "figures",
    "mcmc_trace_m4.png"
  ),
  p_trace,
  width = 9,
  height = 11,
  dpi = 300
)

# ============================================================
# PART F — QA
# ============================================================

cat(
  "\n========== FINAL M1-M4 MODEL EVALUATION ==========\n"
)

cat(
  "\nConvergence:\n"
)

print(
  convergence_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nModel comparison:\n"
)

print(
  model_comparison,
  n = Inf,
  width = Inf
)

cat(
  "\nSequential DIC improvements:\n"
)

model_comparison %>%
  select(
    model,
    improvement_from_previous
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat(
  "\nBest DIC model:",
  best_model,
  "\n"
)

cat(
  "All models pass convergence:",
  all(
    convergence_summary$convergence_pass
  ),
  "\n"
)

cat(
  "==================================================\n"
)

message(
  "Final M1-M4 model evaluation completed."
)