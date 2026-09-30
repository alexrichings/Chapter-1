
# Chapter-1

#-------------------------------------------------------------------------------
  
# Vibe-coding a Bayesian latent class regression framework 
# Aim is to separate out the sub-components within the population e.g. booster and non-boosters, but account for covariates 
  
# Probability of boosting will be influenced by age, baseline titre etc. 
  
#-------------------------------------------------------------------------------
  
# 1. Install and load packages ----

library(cmdstanr)
library(dplyr)
library(posterior)
library(bayesplot)
library(loo)
library(here)

#-------------------------------------------------------------------------------
  
# 2. Load data ---- 

inputs <- readRDS(here("data", "inputs_dr.rds"))

# Scale data for stability 
inputs <- inputs %>%
  mutate(
    baseline_scaled = as.numeric(scale(dengns1_1_mfi_sero / 1000))
  ) %>%
  filter(!is.na(change_dengns1_1_scale), !is.na(dengns1_1_mfi_sero))

--------------------
  
# Covariate matrix ----

# Add covariates that may influence probability of boosting
# Start with baseline only; extend e.g. ~ baseline_scaled + age_scaled + sex
X <- model.matrix(~ baseline_scaled, data = inputs)[, -1, drop = FALSE]

--------------------
  
# Stan data ----

stan_data <- list(
  N     = nrow(inputs),
  K     = ncol(X),
  delta = inputs$change_dengns1_1_scale,
  X     = X
)

# Sanity checks
cat("N:", stan_data$N, "\n")
cat("K:", stan_data$K, "\n")
cat("NAs in delta:", sum(is.na(stan_data$delta)), "\n") # cannot have any NA values or else Stan will error 
cat("Range of delta:", range(stan_data$delta, na.rm = TRUE), "\n")
cat("Range of baseline_scaled:", range(stan_data$X, na.rm = TRUE), "\n")

--------------------
  
# Compile Stan model ----

mod <- cmdstan_model(here("scripts", "blcr", "01-gaussain_gamma_blcr.stan"))

--------------------
  
# Fit model ----

fit_denv1 <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  init            = 0.5
)

--------------------
  
# Diagnostics ----

fit_denv1$summary(c("mu1", "sigma1", "alpha", "beta_rate","theta_intercept", "theta_beta", "prop_booster"))

--------------------
  
# Extract booster probabilities ----

p_booster_draws <- fit_denv1$draws("p_booster", format = "matrix")
inputs$p_booster_denv1 <- apply(p_booster_draws, 2, mean)

# Compare posterior booster proportion to existing model estimate (~6%)
cat("Posterior mean booster proportion:", 
    mean(fit_denv1$draws("prop_booster", format = "matrix")), "\n")

