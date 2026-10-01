
# Chapter-1

#-------------------------------------------------------------------------------

# Vibe-coding a Bayesian latent class regression framework 
# Aim is to separate out the sub-components within the population e.g. booster and non-boosters, but account for covariates 

# Probability of boosting will be influenced by age, baseline titre etc. 

# Using data from Fiji from 2013 - 2015 with decent separation between components 

#-------------------------------------------------------------------------------

# 1. Load and prepare Fiji data ----

# Load data
inputs_fj <- readRDS(here("data", "inputs_fj.rds"))

fiji_inputs <- inputs_fj |>                         
  filter(!is.na(change_denv1_scaled)) |>
  mutate(
    baseline_scaled = as.numeric(scale(DENV1 / 1000))  
  ) |>
  filter(!is.na(baseline_scaled))

#-------------------------------------------------------------------------------

# 2. Build design matrix ----

X_fiji <- model.matrix(~ baseline_scaled, data = fiji_inputs)[, -1, drop = FALSE]

#-------------------------------------------------------------------------------

# 3. Assemble Stan data ----

stan_data_fiji <- list(
  N     = nrow(fiji_inputs),
  K     = ncol(X_fiji),
  delta = fiji_inputs$change_denv1_scaled,
  X     = X_fiji
)

cat("N =", stan_data_fiji$N, "\n")
cat("Range of delta:", range(stan_data_fiji$delta), "\n")
cat("Proportion delta > 0:", mean(stan_data_fiji$delta > 0), "\n")

#-------------------------------------------------------------------------------

# 4. Compile model (reuse if already compiled) ----

mod <- cmdstan_model(here("scripts", "blcr", "01-gaussain_gamma_blcr.stan"))

#-------------------------------------------------------------------------------

# 5. Fit model ----

fit_fiji_denv1 <- mod$sample(
  data          = stan_data_fiji,
  chains        = 4,
  parallel_chains = 4,
  iter_warmup   = 1000,
  iter_sampling = 1000,
  seed          = 42,
  init          = 0.5
)

#-------------------------------------------------------------------------------

# 6. Check convergence ----

fit_fiji_denv1$diagnostic_summary()

fit_fiji_denv1$summary(
  variables = c("mu1", "sigma1", "alpha", "beta_rate",
                "theta_intercept", "prop_booster")
) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

#-------------------------------------------------------------------------------

# 7. Extract booster probabilities ----

p_booster_draws <- fit_fiji_denv1$draws("p_booster", format = "matrix")
fiji_inputs$p_booster_denv1 <- apply(p_booster_draws, 2, mean)

cat("Posterior mean booster proportion:",
    mean(fit_fiji_denv1$draws("prop_booster", format = "matrix")), "\n")

#-------------------------------------------------------------------------------

# 8. Extract posterior draws ----

draws_fiji <- fit_fiji_denv1$draws(
  variables = c("mu1", "sigma1", "alpha", "beta_rate", "prop_booster"),
  format    = "matrix"
)

set.seed(1)
draw_idx      <- sample(nrow(draws_fiji), min(500, nrow(draws_fiji)))
draws_fiji_sub <- draws_fiji[draw_idx, ]

#-------------------------------------------------------------------------------

# 9. Build density grid over observed range ----

x_grid <- seq(
  min(fiji_inputs$change_denv1_scaled, na.rm = TRUE) - 1,
  max(fiji_inputs$change_denv1_scaled, na.rm = TRUE) + 1,
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
  p    <- draws_fiji_sub[i, "prop_booster"]
  mu   <- draws_fiji_sub[i, "mu1"]
  sig  <- draws_fiji_sub[i, "sigma1"]
  alph <- draws_fiji_sub[i, "alpha"]
  rate <- draws_fiji_sub[i, "beta_rate"]
  
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

plot_df_fiji <- bind_rows(
  summ(mix_mat,   "Marginal mixture"),
  summ(nonb_mat,  "Non-booster"),
  summ(boost_mat, "Booster")
) |>
  mutate(component = factor(component,
                            levels = c("Marginal mixture", "Non-booster", "Booster")))

pb_med <- 100 * median(draws_fiji[, "prop_booster"])
pb_lo  <- 100 * quantile(draws_fiji[, "prop_booster"], 0.025)
pb_hi  <- 100 * quantile(draws_fiji[, "prop_booster"], 0.975)

#-------------------------------------------------------------------------------

# 12. Plot ----

cols  <- c("Marginal mixture" = "#1a1a1a",
           "Non-booster"      = "#4477AA",
           "Booster"          = "#EE6677")
fills <- c("Non-booster" = "#4477AA", "Booster" = "#EE6677")

bfmm_fiji <- ggplot() +
  geom_histogram(
    data = fiji_inputs,
    aes(x = change_denv1_scaled, y = after_stat(density)),
    bins = 80, fill = "grey88", colour = "white", linewidth = 0.2
  ) +
  geom_ribbon(
    data = plot_df_fiji |> filter(component %in% c("Non-booster", "Booster")),
    aes(x = x, ymin = lo, ymax = hi, fill = component), alpha = 0.20
  ) +
  geom_ribbon(
    data = plot_df_fiji |> filter(component == "Marginal mixture"),
    aes(x = x, ymin = lo, ymax = hi), fill = "#1a1a1a", alpha = 0.10
  ) +
  geom_line(
    data = plot_df_fiji |> filter(component %in% c("Non-booster", "Booster")),
    aes(x = x, y = med, colour = component),
    linewidth = 0.9, linetype = "dashed"
  ) +
  geom_line(
    data = plot_df_fiji |> filter(component == "Marginal mixture"),
    aes(x = x, y = med), colour = "#1a1a1a", linewidth = 1.1
  ) +
  scale_colour_manual(values = cols[c("Non-booster", "Booster")], name = NULL) +
  scale_fill_manual(values = fills, name = NULL) +
  coord_cartesian(xlim = c(-10, 30)) +
  labs(
    title    = "Gaussian-Gamma mixture model: dengue NS1 ΔMFI (Fiji)",
    subtitle = sprintf("Estimated booster proportion: %.1f%% (95%% CrI: %.1f%%–%.1f%%)",
                       pb_med, pb_lo, pb_hi),
    x = "Change in MFI (scaled)", y = "Density",
    caption = "Shaded bands = posterior 95% CrI. Black = marginal mixture; dashed = individual components."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3))

# Save plot
ggsave(here("outputs/figures/", "bfmm_fiji.png"), plot = bfmm_fiji, width = 13, height = 7, dpi = 300)



