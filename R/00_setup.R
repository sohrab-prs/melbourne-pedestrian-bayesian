# ============================================================
# Melbourne Pedestrian Bayesian Project
# Environment check
# ============================================================

library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(rjags)
library(coda)

cat("\n--- R version ---\n")
print(R.version.string)

cat("\n--- Project root ---\n")
print(here())

cat("\n--- sf version ---\n")
print(packageVersion("sf"))

cat("\n--- rjags version ---\n")
print(packageVersion("rjags"))

cat("\n--- JAGS version ---\n")
print(jags.version())

cat("\n--- Session information ---\n")
sessionInfo()

