
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
library(ggplot2)

#-------------------------------------------------------------------------------
  
# 2. Load data ---- 

inputs <- readRDS(here("data", "inputs_dr.rds"))


# Scale data for stability 
inputs <- inputs %>%
  mutate(
    baseline_scaled = as.numeric(scale(dengns1_1_mfi_sero / 1000))
  ) %>%
  filter(!is.na(change_dengns1_1_scale), !is.na(dengns1_1_mfi_sero))

#-------------------------------------------------------------------------------
  
# 3. Covariate matrix ----

# Add covariates that may influence probability of boosting
# Start with baseline only; extend e.g. ~ baseline_scaled + age_scaled + sex
X <- model.matrix(~ baseline_scaled, data = inputs)[, -1, drop = FALSE]

#-------------------------------------------------------------------------------
  
# 4. Set up Stan data and compilation ----

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

  
# Compile Stan model 

mod <- cmdstan_model(here("scripts", "blcr", "01-gaussain_gamma_blcr.stan"))

#-------------------------------------------------------------------------------
  
# 5. Fit model ----

fit_denv1 <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  seed            = 42,
  init            = 0.5
)

#-------------------------------------------------------------------------------

# 6. Diagnostics ----

fit_denv1$summary(c("mu1", "sigma1", "alpha", "beta_rate","theta_intercept", "theta_beta", "prop_booster"))

#-------------------------------------------------------------------------------
  
# 7. Extract booster probabilities ----

p_booster_draws <- fit_denv1$draws("p_booster", format = "matrix")
inputs$p_booster_denv1 <- apply(p_booster_draws, 2, mean)

# Compare posterior booster proportion to existing model estimate (~6%)
cat("Posterior mean booster proportion:", 
    mean(fit_denv1$draws("prop_booster", format = "matrix")), "\n")

#-------------------------------------------------------------------------------

# 8. Plot mixture model with posterior uncertainty ----

# Extract posterior draws for mixture parameters
draws <- fit_denv1$draws(
  variables = c("mu1", "sigma1", "alpha", "beta_rate", "prop_booster"),
  format    = "matrix"
)

# Thin to 500 draws for speed
set.seed(1)
draw_idx  <- sample(nrow(draws), min(500, nrow(draws)))
draws_sub <- draws[draw_idx, ]

#-------------------------------------------------------------------------------

# 9. Build density grid over observed range ----

x_grid <- seq(
  min(inputs$change_dengns1_1_scale, na.rm = TRUE) - 1,
  max(inputs$change_dengns1_1_scale, na.rm = TRUE) + 1,
  length.out = 600
)

#-------------------------------------------------------------------------------

# 10. Compute component densities for each posterior draw ----

n  <- length(draw_idx)
ng <- length(x_grid)

nonb_mat  <- matrix(NA_real_, n, ng)
boost_mat <- matrix(NA_real_, n, ng)
mix_mat   <- matrix(NA_real_, n, ng)

for (i in seq_len(n)) {
  p    <- draws_sub[i, "prop_booster"]
  mu   <- draws_sub[i, "mu1"]
  sig  <- draws_sub[i, "sigma1"]
  alph <- draws_sub[i, "alpha"]
  rate <- draws_sub[i, "beta_rate"]
  
  nonb  <- (1 - p) * dnorm(x_grid, mu, sig)
  boost <- ifelse(x_grid > 0,
                  p * dgamma(x_grid, shape = alph, rate = rate),
                  0)
  
  nonb_mat[i, ]  <- nonb
  boost_mat[i, ] <- boost
  mix_mat[i, ]   <- nonb + boost
}

#-------------------------------------------------------------------------------

# 11. Summarise draws to median + 95% credible interval ----

summ <- function(mat, label) {
  data.frame(
    x         = x_grid,
    med       = apply(mat, 2, median),
    lo        = apply(mat, 2, quantile, 0.025),
    hi        = apply(mat, 2, quantile, 0.975),
    component = label
  )
}

plot_df <- bind_rows(
  summ(mix_mat,   "Marginal mixture"),
  summ(nonb_mat,  "Non-booster"),
  summ(boost_mat, "Booster")
) |>
  mutate(component = factor(component,
                            levels = c("Marginal mixture", "Non-booster", "Booster")))

# Posterior CI on prop_booster for subtitle
pb_med <- 100 * median(draws[, "prop_booster"])
pb_lo  <- 100 * quantile(draws[, "prop_booster"], 0.025)
pb_hi  <- 100 * quantile(draws[, "prop_booster"], 0.975)

#-------------------------------------------------------------------------------

# 12. Plot ----

cols  <- c("Marginal mixture" = "#1a1a1a",
           "Non-booster"      = "#4477AA",
           "Booster"          = "#EE6677")

fills <- c("Non-booster" = "#4477AA",
           "Booster"     = "#EE6677")

ggplot() +
  geom_histogram(
    data = inputs,
    aes(x = change_dengns1_1_scale, y = after_stat(density)),
    bins      = 80,
    fill      = "grey88",
    colour    = "white",
    linewidth = 0.2
  ) +
  geom_ribbon(
    data = plot_df |> filter(component %in% c("Non-booster", "Booster")),
    aes(x = x, ymin = lo, ymax = hi, fill = component),
    alpha = 0.20
  ) +
  geom_ribbon(
    data = plot_df |> filter(component == "Marginal mixture"),
    aes(x = x, ymin = lo, ymax = hi),
    fill  = "#1a1a1a",
    alpha = 0.10
  ) +
  geom_line(
    data = plot_df |> filter(component %in% c("Non-booster", "Booster")),
    aes(x = x, y = med, colour = component),
    linewidth = 0.9,
    linetype  = "dashed"
  ) +
  geom_line(
    data = plot_df |> filter(component == "Marginal mixture"),
    aes(x = x, y = med),
    colour    = "#1a1a1a",
    linewidth = 1.1
  ) +
  scale_colour_manual(values = cols[c("Non-booster", "Booster")], name = NULL) +
  scale_fill_manual(values   = fills,                              name = NULL) +
  coord_cartesian(xlim = c(-15, 30)) +
  labs(
    title    = "Gaussian-Gamma mixture model: dengue NS1 ΔMFI",
    subtitle = sprintf(
      "Estimated booster proportion: %.1f%% (95%% CrI: %.1f%%–%.1f%%)",
      pb_med, pb_lo, pb_hi),
    x       = "Change in MFI (raw / 1000)",
    y       = "Density",
    caption = "Shaded bands = posterior 95% CrI. Black = marginal mixture; dashed = individual components."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
  )



