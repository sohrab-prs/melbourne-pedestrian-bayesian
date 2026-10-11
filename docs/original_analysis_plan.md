# Original Analysis Plan

> **Historical note.** This document records the analysis plan before the final model-development stage. The initial plan reserved M4 for a possible spatial extension. During model checking, however, M3 reproduced persistent between-sensor variation well but under-reproduced shared date-level variability. This motivated the final M4 specification with both sensor and day random effects. Residual spatial diagnostics subsequently showed only weak remaining spatial dependence, so an additional spatial random-effect model was not fitted. The final methodology and results are documented in the project report.

---

PROJECT:
Bayesian Hierarchical Modelling of Pedestrian Activity in Melbourne

STUDY AREA:
City of Melbourne, Victoria, Australia

STUDY PERIOD:
1 January 2025 – 31 December 2025

OBSERVATIONAL UNIT:
Pedestrian sensor × calendar day

PRIMARY OUTCOME:
Daily total pedestrian count

PRIMARY SPATIAL SCALE:
500 m circular catchment around each sensor

CANDIDATE EXPLANATORY FAMILIES:
1. Calendar effects
2. Weather
3. Employment/economic activity
4. Residential/demographic characteristics
5. Land use/activity mix
6. Urban morphology

PRIMARY RESEARCH QUESTION:
How do temporal conditions and characteristics of the surrounding urban environment explain variation in pedestrian activity across central Melbourne, and how much unexplained heterogeneity remains between monitoring locations?

INFERENTIAL QUESTIONS:
1. Which temporal factors are associated with pedestrian activity?
2. Which urban characteristics are associated with pedestrian activity?
3. Does overdispersion make Negative Binomial preferable to Poisson?
4. How large is unexplained between-sensor heterogeneity?
5. Can the Bayesian model reproduce important observed spatial and temporal characteristics of pedestrian activity?

PLANNED MODELS:
M1: Bayesian Poisson regression
M2: Bayesian Negative Binomial regression
M3: Bayesian hierarchical Negative Binomial regression
M4: Optional spatial model, conditional on evidence of residual spatial structure

PRIMARY BAYESIAN IMPLEMENTATION:
R + JAGS

MODEL ASSESSMENT:
MCMC convergence diagnostics
DIC
Posterior predictive checking
Prediction
Parameter-recovery simulation

BONUS:
Formal Bayesian model diagnostics
Frequentist comparison

TEMPORAL ALIGNMENT OF CONTEXTUAL DATA:
Pedestrian activity and weather are measured in 2025.
CLUE land-use/employment/residential variables use 2024 data.
Building morphology uses the 2023 City of Melbourne footprint layer.
Resident demographic characteristics use the 2021 ABS Census.

These contextual variables represent relatively slowly changing features of the surrounding urban environment.

SENSOR ELIGIBILITY:

A pedestrian monitoring location is eligible for the main analysis if:

1. it has one valid recorded spatial location during 2025;
2. no within-year location change is detected;
3. at least 90% of the expected hourly observations are available; and
4. at least 300 complete 24-hour sensor-days remain after daily-data cleaning.

Daily pedestrian totals will be calculated only from days containing all 24 hourly observations. Incomplete sensor-days will not be imputed for the primary analysis.