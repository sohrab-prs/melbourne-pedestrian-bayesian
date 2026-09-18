#!/usr/bin/env bash

set -u

PROJECT="/home/sohrab/projects/melbourne-pedestrian-bayes"
RUN_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "$PROJECT" || exit 90

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate melbourne-bayes

export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

date -Iseconds > "$RUN_DIR/started_at.txt"

git rev-parse HEAD > "$RUN_DIR/git_head.txt"

sha256sum \
  R/11c_parameter_recovery_m4_confirmation.R \
  models/jags/m4_sensor_day_negbin.jags \
  data/processed/modelling_dataset_selected_2025.csv \
  renv.lock \
  > "$RUN_DIR/input_sha256.txt"

Rscript --vanilla R/11c_parameter_recovery_m4_confirmation.R \
  > "$RUN_DIR/stdout.log" \
  2> "$RUN_DIR/stderr.log"

EXIT_CODE=$?

echo "$EXIT_CODE" > "$RUN_DIR/exit_code.txt"
date -Iseconds > "$RUN_DIR/finished_at.txt"

if [ "$EXIT_CODE" -eq 0 ]; then

  mkdir -p "$RUN_DIR/artifacts"

  cp outputs/mcmc/m4_parameter_recovery_confirmation_samples.rds \
     "$RUN_DIR/artifacts/"

  cp outputs/tables/m4_parameter_recovery_confirmation_core.csv \
     "$RUN_DIR/artifacts/"

  cp outputs/tables/m4_parameter_recovery_confirmation_sensor_effects.csv \
     "$RUN_DIR/artifacts/"

  cp outputs/tables/m4_parameter_recovery_confirmation_day_effects.csv \
     "$RUN_DIR/artifacts/"

  cp outputs/tables/m4_parameter_recovery_confirmation_convergence.csv \
     "$RUN_DIR/artifacts/"

  cp outputs/tables/m4_parameter_recovery_confirmation_metrics.csv \
     "$RUN_DIR/artifacts/"

  cp outputs/figures/m4_parameter_recovery_confirmation_core.png \
     "$RUN_DIR/artifacts/"

  cp outputs/figures/m4_parameter_recovery_confirmation_sensor.png \
     "$RUN_DIR/artifacts/"

  cp outputs/figures/m4_parameter_recovery_confirmation_day.png \
     "$RUN_DIR/artifacts/"

  sha256sum "$RUN_DIR"/artifacts/* \
    > "$RUN_DIR/output_sha256.txt"
fi

exit "$EXIT_CODE"
