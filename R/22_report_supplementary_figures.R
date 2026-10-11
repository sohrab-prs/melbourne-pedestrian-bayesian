# ============================================================
# 22_report_supplementary_figures.R
# Additional report-only EDA and sensor-effect summaries.
#
# This script uses the locked modelling dataset and saved M4 posterior samples.
# ============================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)

dir.create(
  here("outputs", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here("outputs", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

report_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      colour = "#1F2933"
    ),
    plot.subtitle = element_text(
      size = 11.5,
      colour = "#4B5563"
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      size = 10.5,
      colour = "#30343B"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11.5
    ),
    strip.background = element_rect(
      fill = "#EEF3F7",
      colour = "#CBD5DD",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.35,
      colour = "#E5E9EC"
    ),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

save_plot <- function(plot, filename, width, height) {
  ggsave(
    here("outputs", "figures", filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

# ============================================================
# DATA
# ============================================================

mod_dat <- read_csv(
  here(
    "data",
    "processed",
    "modelling_dataset_selected_2025.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    location_id = as.character(location_id),
    sensor_index = as.integer(sensor_index),
    date = as.Date(date)
  )

stopifnot(
  nrow(mod_dat) == 18671,
  dplyr::n_distinct(mod_dat$sensor_index) == 53,
  dplyr::n_distinct(mod_dat$date) == 364
)

# ============================================================
# FIGURE A — WITHIN-YEAR NETWORK TEMPORAL TREND
# ============================================================

daily_network <- mod_dat %>%
  group_by(date) %>%
  summarise(
    mean_count = mean(daily_ped_count),
    n_sensors = n_distinct(sensor_index),
    .groups = "drop"
  ) %>%
  arrange(date)

# Centered 14-observation moving average.
daily_network$rolling_14 <- as.numeric(
  stats::filter(
    daily_network$mean_count,
    filter = rep(1 / 14, 14),
    sides = 2
  )
)

p_temporal <- ggplot(
  daily_network,
  aes(x = date, y = mean_count)
) +
  geom_line(
    colour = "#8EC9E8",
    linewidth = 0.48,
    alpha = 0.70
  ) +
  geom_line(
    aes(y = rolling_14),
    colour = "#0072B2",
    linewidth = 1.05,
    na.rm = TRUE
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    title = "Network-wide pedestrian activity through 2025",
    subtitle = "Light line: daily mean across available retained sensors; dark line: centered 14-observation moving average",
    x = "Date",
    y = "Mean daily pedestrian count"
  ) +
  report_theme

save_plot(
  p_temporal,
  "report_eda_temporal_trend.png",
  10.0,
  5.6
)

# ============================================================
# FIGURE B — BIVARIATE SENSOR-LEVEL SPATIAL RELATIONSHIPS
# ============================================================

sensor_level <- mod_dat %>%
  group_by(sensor_index, location_id) %>%
  summarise(
    mean_count = mean(daily_ped_count),
    z_business_count = first(z_business_count),
    z_resident_population = first(z_resident_population),
    z_built_coverage = first(z_built_coverage),
    z_structure_height = first(z_structure_height),
    .groups = "drop"
  )

spatial_long <- sensor_level %>%
  pivot_longer(
    cols = c(
      z_business_count,
      z_resident_population,
      z_built_coverage,
      z_structure_height
    ),
    names_to = "predictor",
    values_to = "z_value"
  ) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = c(
        "z_business_count",
        "z_resident_population",
        "z_built_coverage",
        "z_structure_height"
      ),
      labels = c(
        "Business count",
        "Resident population",
        "Built coverage",
        "Median structure height"
      )
    )
  )

spatial_cols <- c(
  "Business count" = "#0072B2",
  "Resident population" = "#009E73",
  "Built coverage" = "#D55E00",
  "Median structure height" = "#CC79A7"
)

p_bivariate <- ggplot(
  spatial_long,
  aes(
    x = z_value,
    y = mean_count,
    colour = predictor
  )
) +
  geom_point(
    size = 2.4,
    alpha = 0.82
  ) +
  geom_smooth(
    aes(fill = predictor),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.85,
    alpha = 0.13
  ) +
  facet_wrap(
    ~ predictor,
    scales = "free_x",
    ncol = 2
  ) +
  scale_colour_manual(
    values = spatial_cols,
    guide = "none"
  ) +
  scale_fill_manual(
    values = spatial_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ",")
  ) +
  labs(
    title = "Descriptive sensor-level relationships with urban context",
    subtitle = "Each point is one retained sensor; lines are unadjusted linear summaries and are not causal estimates",
    x = "Standardized spatial predictor",
    y = "Mean daily pedestrian count"
  ) +
  report_theme

save_plot(
  p_bivariate,
  "report_eda_spatial_bivariate.png",
  10.0,
  7.0
)

# ============================================================
# FIGURE C + TABLE — FINAL M4 SENSOR RANDOM EFFECTS
# ============================================================

m4_rds <- here(
  "outputs",
  "mcmc",
  "m4_sensor_day_negbin_samples.rds"
)

if (!file.exists(m4_rds)) {
  stop(
    "Could not find outputs/mcmc/m4_sensor_day_negbin_samples.rds. ",
    "The supplementary sensor-effect figure requires the saved final M4 posterior."
  )
}

m4_samples <- readRDS(m4_rds)
m4_matrix <- as.matrix(m4_samples)

u_names <- grep(
  "^u_sensor\\[",
  colnames(m4_matrix),
  value = TRUE
)

if (length(u_names) != 53) {
  stop(
    "Expected 53 u_sensor[...] columns in the final M4 posterior, but found ",
    length(u_names),
    "."
  )
}

u_index <- as.integer(
  sub(
    ".*\\[([0-9]+)\\].*",
    "\\1",
    u_names
  )
)

u_names <- u_names[order(u_index)]
u_index <- sort(u_index)

u_draws <- m4_matrix[
  ,
  u_names,
  drop = FALSE
]

u_summary <- tibble(
  sensor_index = u_index,
  posterior_median = apply(
    u_draws,
    2,
    median
  ),
  q025 = apply(
    u_draws,
    2,
    quantile,
    probs = 0.025
  ),
  q975 = apply(
    u_draws,
    2,
    quantile,
    probs = 0.975
  )
)

sensor_lookup <- mod_dat %>%
  distinct(
    sensor_index,
    location_id
  )

u_summary <- u_summary %>%
  left_join(
    sensor_lookup,
    by = "sensor_index"
  ) %>%
  arrange(
    posterior_median
  ) %>%
  mutate(
    ordered_sensor = row_number(),
    direction = if_else(
      posterior_median >= 0,
      "Positive",
      "Negative"
    )
  )

p_sensor_effects <- ggplot(
  u_summary,
  aes(
    x = ordered_sensor,
    y = posterior_median,
    colour = direction
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "#5F6368",
    linewidth = 0.75,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = q025,
      ymax = q975
    ),
    width = 0,
    linewidth = 0.70,
    alpha = 0.76
  ) +
  geom_point(
    size = 2.35
  ) +
  scale_colour_manual(
    values = c(
      "Negative" = "#2166AC",
      "Positive" = "#B2182B"
    ),
    name = "Posterior median"
  ) +
  scale_x_continuous(
    breaks = c(1, 10, 20, 30, 40, 50, 53)
  ) +
  labs(
    title = "Persistent M4 sensor random effects",
    subtitle = "Sensors ordered by posterior median; vertical intervals are 95% credible intervals on the log-count scale",
    x = "Sensor order (low to high posterior median effect)",
    y = "Sensor random effect"
  ) +
  report_theme +
  theme(
    panel.grid.major.x = element_blank()
  )

save_plot(
  p_sensor_effects,
  "report_m4_sensor_random_effects.png",
  10.2,
  5.9
)

n_extreme <- 5

extremes <- bind_rows(
  u_summary %>%
    slice_head(n = n_extreme) %>%
    mutate(group = "Lowest"),
  u_summary %>%
    slice_tail(n = n_extreme) %>%
    arrange(desc(posterior_median)) %>%
    mutate(group = "Highest")
) %>%
  select(
    group,
    sensor_index,
    location_id,
    posterior_median,
    q025,
    q975
  )

write_csv(
  extremes,
  here(
    "outputs",
    "tables",
    "report_m4_sensor_random_effect_extremes.csv"
  )
)

cat(
  "\nSupplementary report figures/tables generated successfully.\n",
  "Created:\n",
  "  outputs/figures/report_eda_temporal_trend.png\n",
  "  outputs/figures/report_eda_spatial_bivariate.png\n",
  "  outputs/figures/report_m4_sensor_random_effects.png\n",
  "  outputs/tables/report_m4_sensor_random_effect_extremes.csv\n"
)
