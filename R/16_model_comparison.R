# ============================================================
# 16_model_comparison.R
# Melbourne Pedestrian Bayesian Project
#
# Purpose:
# Compare M1, M2 and M3 using:
#   - DIC
#   - effective parameter count
#   - convergence diagnostics
#   - key posterior dispersion / heterogeneity parameters
#
# This is a light post-processing step and should run locally.
#
# Outputs:
#   outputs/tables/model_comparison_dic.csv
#   outputs/tables/model_comparison_summary.csv
#   outputs/tables/model_comparison_key_parameters.csv
#   outputs/figures/model_comparison_dic.png
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

# ------------------------------------------------------------
# 1. Validate required files
# ------------------------------------------------------------

required_files <- c(
  
  here(
    "outputs",
    "tables",
    "m1_poisson_fit_metrics.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m2_negbin_fit_metrics.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_fit_metrics.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "mcmc_diagnostics_summary.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m1_poisson_posterior_summary.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m2_negbin_posterior_summary.csv"
  ),
  
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_posterior_summary.csv"
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
    "tables"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here(
    "outputs",
    "figures"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# PART A — DIC
# ============================================================

m1_fit <- read_csv(
  here(
    "outputs",
    "tables",
    "m1_poisson_fit_metrics.csv"
  ),
  show_col_types = FALSE
)

m2_fit <- read_csv(
  here(
    "outputs",
    "tables",
    "m2_negbin_fit_metrics.csv"
  ),
  show_col_types = FALSE
)

m3_fit <- read_csv(
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_fit_metrics.csv"
  ),
  show_col_types = FALSE
)

# Standardise columns before binding.
m1_fit_standard <- m1_fit %>%
  transmute(
    model =
      model,
    
    observations =
      observations,
    
    sensors =
      NA_integer_,
    
    mean_deviance =
      mean_deviance,
    
    effective_parameters_pD =
      effective_parameters_pD,
    
    DIC =
      DIC
  )

m2_fit_standard <- m2_fit %>%
  transmute(
    model =
      model,
    
    observations =
      observations,
    
    sensors =
      NA_integer_,
    
    mean_deviance =
      mean_deviance,
    
    effective_parameters_pD =
      effective_parameters_pD,
    
    DIC =
      DIC
  )

m3_fit_standard <- m3_fit %>%
  transmute(
    model =
      model,
    
    observations =
      observations,
    
    sensors =
      sensors,
    
    mean_deviance =
      mean_deviance,
    
    effective_parameters_pD =
      effective_parameters_pD,
    
    DIC =
      DIC
  )

dic_comparison <- bind_rows(
  m1_fit_standard,
  m2_fit_standard,
  m3_fit_standard
) %>%
  arrange(
    DIC
  ) %>%
  mutate(
    
    dic_rank =
      row_number(),
    
    delta_DIC =
      DIC -
      min(
        DIC
      ),

  )

# Because arranging by DIC gives best -> worst,
# improvement_vs_previous is not the sequential M1 -> M2 -> M3
# comparison we want for interpretation.
#
# Add explicit pairwise improvements separately.

m1_dic <- dic_comparison$DIC[
  dic_comparison$model ==
    "M1 Poisson"
]

m2_dic <- dic_comparison$DIC[
  dic_comparison$model ==
    "M2 Negative Binomial"
]

m3_dic <- dic_comparison$DIC[
  dic_comparison$model ==
    "M3 Hierarchical Negative Binomial"
]

stopifnot(
  length(m1_dic) == 1,
  length(m2_dic) == 1,
  length(m3_dic) == 1
)

dic_comparison <- dic_comparison %>%
  mutate(
    
    DIC_improvement_from_M1 =
      case_when(
        
        model ==
          "M1 Poisson" ~
          0,
        
        TRUE ~
          m1_dic -
          DIC
      ),
    
    DIC_improvement_from_previous_model =
      case_when(
        
        model ==
          "M1 Poisson" ~
          NA_real_,
        
        model ==
          "M2 Negative Binomial" ~
          m1_dic -
          m2_dic,
        
        model ==
          "M3 Hierarchical Negative Binomial" ~
          m2_dic -
          m3_dic,
        
        TRUE ~
          NA_real_
      )
  )

write_csv(
  dic_comparison,
  here(
    "outputs",
    "tables",
    "model_comparison_dic.csv"
  )
)

# ============================================================
# PART B — CONVERGENCE + DIC SUMMARY
# ============================================================

mcmc_summary <- read_csv(
  here(
    "outputs",
    "tables",
    "mcmc_diagnostics_summary.csv"
  ),
  show_col_types = FALSE
)

model_comparison_summary <- dic_comparison %>%
  
  select(
    model,
    observations,
    mean_deviance,
    effective_parameters_pD,
    DIC,
    delta_DIC,
    DIC_improvement_from_M1,
    DIC_improvement_from_previous_model,
    dic_rank
  ) %>%
  
  left_join(
    mcmc_summary,
    by = "model"
  ) %>%
  
  arrange(
    dic_rank
  )

write_csv(
  model_comparison_summary,
  here(
    "outputs",
    "tables",
    "model_comparison_summary.csv"
  )
)

# ============================================================
# PART C — KEY MODEL-SPECIFIC PARAMETERS
# ============================================================

m1_post <- read_csv(
  here(
    "outputs",
    "tables",
    "m1_poisson_posterior_summary.csv"
  ),
  show_col_types = FALSE
)

m2_post <- read_csv(
  here(
    "outputs",
    "tables",
    "m2_negbin_posterior_summary.csv"
  ),
  show_col_types = FALSE
)

m3_post <- read_csv(
  here(
    "outputs",
    "tables",
    "m3_hierarchical_negbin_posterior_summary.csv"
  ),
  show_col_types = FALSE
)

phi_m2 <- m2_post %>%
  filter(
    parameter == "phi"
  ) %>%
  transmute(
    model =
      "M2 Negative Binomial",
    
    parameter =
      "phi",
    
    posterior_mean,
    q025,
    median,
    q975
  )

phi_m3 <- m3_post %>%
  filter(
    parameter == "phi"
  ) %>%
  transmute(
    model =
      "M3 Hierarchical Negative Binomial",
    
    parameter =
      "phi",
    
    posterior_mean,
    q025,
    median,
    q975
  )

sigma_m3 <- m3_post %>%
  filter(
    parameter == "sigma_sensor"
  ) %>%
  transmute(
    model =
      "M3 Hierarchical Negative Binomial",
    
    parameter =
      "sigma_sensor",
    
    posterior_mean,
    q025,
    median,
    q975
  )

key_parameters <- bind_rows(
  phi_m2,
  phi_m3,
  sigma_m3
)

write_csv(
  key_parameters,
  here(
    "outputs",
    "tables",
    "model_comparison_key_parameters.csv"
  )
)

# ============================================================
# PART D — DIC FIGURE
# ============================================================

dic_plot_data <- dic_comparison %>%
  mutate(
    model_short =
      factor(
        case_when(
          model == "M1 Poisson" ~ "M1 Poisson",
          model == "M2 Negative Binomial" ~ "M2 NB",
          model == "M3 Hierarchical Negative Binomial" ~ "M3 Hierarchical NB"
        ),
        levels = c(
          "M1 Poisson",
          "M2 NB",
          "M3 Hierarchical NB"
        )
      )
  )

p_dic <- ggplot(
  dic_plot_data,
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
          round(
            DIC
          ),
          big.mark = ",",
          scientific = FALSE
        )
    ),
    vjust = -0.4
  ) +
  labs(
    title =
      "Bayesian Model Comparison by DIC",
    subtitle =
      "Lower DIC indicates better relative fit after complexity penalty",
    x =
      NULL,
    y =
      "DIC"
  ) +
  theme_minimal()

ggsave(
  here(
    "outputs",
    "figures",
    "model_comparison_dic.png"
  ),
  p_dic,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# PART E — QA
# ============================================================

cat(
  "\n========== STEP 16 MODEL COMPARISON ==========\n"
)

cat(
  "\nDIC comparison:\n"
)

print(
  dic_comparison,
  width = Inf
)

cat(
  "\nCombined model comparison:\n"
)

print(
  model_comparison_summary,
  width = Inf
)

cat(
  "\nKey model-specific parameters:\n"
)

print(
  key_parameters,
  width = Inf
)

cat(
  "\nDIC improvement M1 -> M2:",
  m1_dic - m2_dic,
  "\n"
)

cat(
  "DIC improvement M2 -> M3:",
  m2_dic - m3_dic,
  "\n"
)

cat(
  "Best DIC model:",
  dic_comparison$model[
    which.min(
      dic_comparison$DIC
    )
  ],
  "\n"
)

cat(
  "==============================================\n"
)

message(
  "Step 16 model comparison completed."
)