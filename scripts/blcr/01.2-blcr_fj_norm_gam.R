
# Chapter-1

#-------------------------------------------------------------------------------

# Vibe-coding a Bayesian latent class regression framework 
# Aim is to separate out the sub-components within the population e.g. booster and non-boosters, but account for covariates 

# Probability of boosting will be influenced by age, baseline titre etc. 

# Using data from Fiji from 2013 - 2015 with decent separation between components 

#-------------------------------------------------------------------------------
# Load packages ----

library(cmdstanr)
library(dplyr)
library(posterior)
library(bayesplot)
library(loo)
library(here)
library(ggplot2)
library(tidyr)
library(patchwork)

#-------------------------------------------------------------------------------
# Section 1: FMM for one serotype ----
#-------------------------------------------------------------------------------

# 1.1. Load and prepare Fiji data ----

# Load data
inputs_fj <- readRDS(here("data", "inputs_fj.rds"))

inputs_fj$change_denv2_scaled <- (inputs_fj$DENV2B - inputs_fj$DENV2)/1000
inputs_fj$change_denv4_scaled <- (inputs_fj$DENV4B - inputs_fj$DENV4)/1000

fiji_inputs <- inputs_fj |>
  filter(
    !is.na(change_denv1_scaled),
    !is.na(change_denv2_scaled),
    !is.na(change_denv3_scaled),
    !is.na(change_denv4_scaled)
  ) |>
  mutate(
    baseline_denv1_scaled = as.numeric(scale(DENV1 / 1000)),
    baseline_denv2_scaled = as.numeric(scale(DENV2 / 1000)),
    baseline_denv3_scaled = as.numeric(scale(DENV3 / 1000)),
    baseline_denv4_scaled = as.numeric(scale(DENV4 / 1000))
  ) |>
  filter(
    !is.na(baseline_denv1_scaled),
    !is.na(baseline_denv2_scaled),
    !is.na(baseline_denv3_scaled),
    !is.na(baseline_denv4_scaled)
  )

#-------------------------------------------------------------------------------

# 1.2. Build design matrix ----

X_fiji <- model.matrix(~ baseline_denv2_scaled, data = fiji_inputs)[, -1, drop = FALSE]

#-------------------------------------------------------------------------------

# 1.3. Assemble Stan data ----

stan_data_fiji <- list(
  N     = nrow(fiji_inputs),
  K     = ncol(X_fiji),
  delta = fiji_inputs$change_denv2_scaled,
  X     = X_fiji
)

cat("N =", stan_data_fiji$N, "\n")
cat("Range of delta:", range(stan_data_fiji$delta), "\n")
cat("Proportion delta > 0:", mean(stan_data_fiji$delta > 0), "\n")

#-------------------------------------------------------------------------------

# 1.4. Compile model (reuse if already compiled) ----

mod <- cmdstan_model(here("scripts", "blcr", "01-blcr_norm_gam.stan"))

#-------------------------------------------------------------------------------

# 1.5. Fit model ----

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

# 1.6. Check convergence ----

fit_fiji_denv1$diagnostic_summary()

fit_fiji_denv1$summary(
  variables = c("mu1", "sigma1", "alpha", "beta_rate",
                "theta_intercept", "prop_booster")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

#-------------------------------------------------------------------------------

# 1.7. Extract booster probabilities ----

p_booster_draws <- fit_fiji_denv1$draws("p_booster", format = "matrix")
fiji_inputs$p_booster_denv1 <- apply(p_booster_draws, 2, mean)

cat("Posterior mean booster proportion:",
    mean(fit_fiji_denv1$draws("prop_booster", format = "matrix")), "\n")

#-------------------------------------------------------------------------------

# 1.8. Extract posterior draws ----

draws_fiji <- fit_fiji_denv1$draws(
  variables = c("mu1", "sigma1", "alpha", "beta_rate", "prop_booster"),
  format    = "matrix"
)

set.seed(1)
draw_idx      <- sample(nrow(draws_fiji), min(500, nrow(draws_fiji)))
draws_fiji_sub <- draws_fiji[draw_idx, ]

#-------------------------------------------------------------------------------

# 1.9. Build density grid over observed range ----

x_grid <- seq(
  min(fiji_inputs$change_denv1_scaled, na.rm = TRUE) - 1,
  max(fiji_inputs$change_denv1_scaled, na.rm = TRUE) + 1,
  length.out = 600
)

#-------------------------------------------------------------------------------

# 1.10. Compute component densities for each posterior draw ----

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

# 1.11. Summarise draws to median + 95% credible interval ----

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

# 1.12. Plot ----

cols  <- c("Marginal mixture" = "#1a1a1a",
           "Non-booster"      = "#4477AA",
           "Booster"          = "#EE6677")
fills <- c("Non-booster" = "#4477AA", "Booster" = "#EE6677")

bfmm_fiji <- ggplot() +
  geom_histogram(
    data = fiji_inputs,
    aes(x = change_denv2_scaled, y = after_stat(density)),
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
    title    = "Gaussian-Gamma mixture model: DENV2 ΔMFI (Fiji)",
    subtitle = sprintf("Estimated booster proportion: %.1f%% (95%% CrI: %.1f%%–%.1f%%)",
                       pb_med, pb_lo, pb_hi),
    x = "Change in MFI (scaled)", y = "Density",
    caption = "Shaded bands = posterior 95% CrI. Black = marginal mixture; dashed = individual components."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3))

bfmm_fiji

# Save plot
ggsave(here("outputs/figures/", "bfmm_fiji.png"), plot = bfmm_fiji, width = 13, height = 7, dpi = 300)



#-------------------------------------------------------------------------------
# Section 2: FMM for multiple serotypes simultaneously ----
#-------------------------------------------------------------------------------

# 2.1. Fit the model ----

# compile model 
mod <- cmdstan_model(here("scripts", "blcr", "01-blcr_norm_gam.stan"))

# create a list of serotypes with the relevant data 
serotypes <- list(
  list(name = "DENV1", delta_col = "change_denv1_scaled", baseline_col = "baseline_denv1_scaled",
       mfi_col = "DENV1", followup_col = "DENV1B", p_booster_col = "p_booster_denv1"),
  list(name = "DENV2", delta_col = "change_denv2_scaled", baseline_col = "baseline_denv2_scaled",
       mfi_col = "DENV2", followup_col = "DENV2B", p_booster_col = "p_booster_denv2"),
  list(name = "DENV3", delta_col = "change_denv3_scaled", baseline_col = "baseline_denv3_scaled",
       mfi_col = "DENV3", followup_col = "DENV3B", p_booster_col = "p_booster_denv3"),
  list(name = "DENV4", delta_col = "change_denv4_scaled", baseline_col = "baseline_denv4_scaled",
       mfi_col = "DENV4", followup_col = "DENV4B", p_booster_col = "p_booster_denv4")
)

# create lists for each output
fits_fiji    <- list()
results_fiji <- list()
plots_fiji   <- list()

# set features for plotting 
cols  <- c("Marginal mixture" = "#1a1a1a", "Non-booster" = "#4477AA", "Booster" = "#EE6677")
fills <- c("Non-booster" = "#4477AA", "Booster" = "#EE6677")

# define a common axis limit across serotypes 
y_max <- max(sapply(serotypes, function(s) {
  v <- fiji_inputs[[s$delta_col]]
  max(density(v[!is.na(v)])$y)
})) * 1.05

# fit the model and plot a FMM for each serotype 
for (s in serotypes) {
  
  delta_vec    <- fiji_inputs[[s$delta_col]]
  baseline_vec <- fiji_inputs[[s$baseline_col]]
  
  X <- model.matrix(~ baseline_vec)[, -1, drop = FALSE]
  
  stan_data <- list(
    N     = length(delta_vec),
    K     = ncol(X),
    delta = delta_vec,
    X     = X
  )
  
  fit <- mod$sample(
    data            = stan_data,
    seed            = 123,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    iter_sampling   = 1000,
    refresh         = 500
  )
  
  fits_fiji[[s$name]] <- fit
  
  draws <- fit$draws(
    c("mu1", "sigma1", "alpha", "beta_rate", "theta_intercept", "prop_booster"),
    format = "matrix"
  )
  
  results_fiji[[s$name]] <- list(
    summary      = fit$summary(),
    prop_booster = mean(draws[, "prop_booster"])
  )
  
  # Thin to 500 draws for plotting
  idx    <- sample(nrow(draws), 500)
  draws  <- draws[idx, ]
  
  x_min  <- min(fiji_inputs[[s$delta_col]], na.rm = TRUE)
  x_max  <- max(fiji_inputs[[s$delta_col]], na.rm = TRUE)
  x_grid <- seq(x_min - 0.5, x_max + 0.5, length.out = 400)
  
  plot_df <- purrr::map_dfr(seq_len(nrow(draws)), function(i) {
    mu    <- draws[i, "mu1"]
    sig   <- draws[i, "sigma1"]
    alph  <- draws[i, "alpha"]
    rate  <- draws[i, "beta_rate"]
    p     <- draws[i, "prop_booster"]
    
    nonb  <- (1 - p) * dnorm(x_grid, mean = mu, sd = sig)
    boost <- p * ifelse(x_grid > 0, dgamma(x_grid, shape = alph, rate = rate), 0)
    marg  <- nonb + boost
    
    data.frame(
      x            = x_grid,
      nonbooster   = nonb,
      booster      = boost,
      marginal     = marg,
      draw         = i
    )
  })
  
  plot_df_summary <- plot_df |>
    group_by(x) |>
    summarise(
      nonbooster_med = median(nonbooster),
      nonbooster_lo  = quantile(nonbooster, 0.025),
      nonbooster_hi  = quantile(nonbooster, 0.975),
      booster_med    = median(booster),
      booster_lo     = quantile(booster, 0.025),
      booster_hi     = quantile(booster, 0.975),
      marginal_med   = median(marginal),
      marginal_lo    = quantile(marginal, 0.025),
      marginal_hi    = quantile(marginal, 0.975),
      .groups        = "drop"
    )
  
  plots_fiji[[s$name]] <- ggplot() +
    geom_histogram(
      data    = fiji_inputs,
      aes(x = .data[[s$delta_col]], y = after_stat(density), fill = "Non-booster"),
      bins    = 60,
      alpha   = 0.4,
      colour  = NA
    ) +
    geom_ribbon(
      data  = plot_df_summary,
      aes(x = x, ymin = nonbooster_lo, ymax = nonbooster_hi),
      fill  = "#4477AA", alpha = 0.2
    ) +
    geom_ribbon(
      data  = plot_df_summary,
      aes(x = x, ymin = booster_lo, ymax = booster_hi),
      fill  = "#EE6677", alpha = 0.2
    ) +
    geom_ribbon(
      data  = plot_df_summary,
      aes(x = x, ymin = marginal_lo, ymax = marginal_hi),
      fill  = "#1a1a1a", alpha = 0.15
    ) +
    geom_line(
      data  = plot_df_summary,
      aes(x = x, y = nonbooster_med, colour = "Non-booster")
    ) +
    geom_line(
      data  = plot_df_summary,
      aes(x = x, y = booster_med, colour = "Booster")
    ) +
    geom_line(
      data  = plot_df_summary,
      aes(x = x, y = marginal_med, colour = "Marginal mixture"),
      linewidth = 1
    ) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = fills, name = NULL, guide = "none") +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title = s$name,
      x     = "ΔMFI (scaled)",
      y     = "Density"
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom")
}

# Save individual plots
ggsave(here("outputs/figures/Fiji", "fiji_denv1.png"), plots_fiji$DENV1, width = 8, height = 5, dpi = 300)
ggsave(here("outputs/figures/Fiji", "fiji_denv2.png"), plots_fiji$DENV2, width = 8, height = 5, dpi = 300)
ggsave(here("outputs/figures/Fiji", "fiji_denv3.png"), plots_fiji$DENV3, width = 8, height = 5, dpi = 300)
ggsave(here("outputs/figures/Fiji", "fiji_denv4.png"), plots_fiji$DENV4, width = 8, height = 5, dpi = 300)

# Save joint plot
joint_plot <- (plots_fiji$DENV1 + plots_fiji$DENV2) / (plots_fiji$DENV3 + plots_fiji$DENV4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("outputs/figures/Fiji", "fiji_all_serotypes.png"), joint_plot, width = 14, height = 10, dpi = 300)


#-------------------------------------------------------------------------------

# 2.2. Consider antibody boost summary ---- 

# Extract prop_booster summary for each serotype
prop_booster_summary <- purrr::map_dfr(names(fits_fiji), function(nm) {
  draws <- fits_fiji[[nm]]$draws("prop_booster", format = "matrix")
  data.frame(
    serotype = nm,
    median   = round(median(draws) * 100, 1),
    lo       = round(quantile(draws, 0.025) * 100, 1),
    hi       = round(quantile(draws, 0.975) * 100, 1)
  )
})

print(prop_booster_summary)

# Plot the antibody boost summary per serotype 
prop_booster_plot <- ggplot(prop_booster_summary, aes(x = serotype, y = median)) +
  geom_col(fill = "#2c7bb6", alpha = 0.5, colour = "black", width = 0.5) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width     = 0.15,
    linewidth = 0.8,
    colour    = "#1a1a1a"
  ) +
  geom_text(
    aes(y = hi, label = sprintf("%.1f%% (%.1f–%.1f%%)", median, lo, hi)),
    vjust = -0.8,
    size  = 3.5
  ) +
  scale_y_continuous(
    name   = "Estimated proportion boosting (%)",
    limits = c(0, 100),
    labels = function(x) paste0(x, "%")
  ) +
  labs(x = "Serotype", title = "Proportion boosting by DENV serotype (Fiji)") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(here("outputs/figures/Fiji", "fiji_prop_boosting.png"), prop_booster_plot, width = 8, height = 5, dpi = 300)


#-------------------------------------------------------------------------------
# Section 3: Post-FMM analyses - understanding the covariates 
#-------------------------------------------------------------------------------

# 3.1. Compute individual-level P(booster) for all serotypes ----

for (s in serotypes) {
  reg_draws <- fits_fiji[[s$name]]$draws(
    variables = c("theta_intercept", "theta_beta[1]"),
    format    = "matrix"
  )
  fiji_inputs[[s$p_booster_col]] <- plogis(
    median(reg_draws[, "theta_intercept"]) +
      median(reg_draws[, "theta_beta[1]"]) * fiji_inputs[[s$baseline_col]]
  )
}

#-------------------------------------------------------------------------------

# 3.2. Odds ratios for baseline titre effect ----

or_summary <- purrr::map_dfr(serotypes, function(s) {
  draws <- fits_fiji[[s$name]]$draws("theta_beta[1]", format = "matrix")
  or    <- exp(draws)
  data.frame(
    serotype = s$name,
    median   = round(median(or), 2),
    lo       = round(quantile(or, 0.025), 2),
    hi       = round(quantile(or, 0.975), 2)
  )
})

print(or_summary)

# Forest plot of ORs
or_plot <- ggplot(or_summary, aes(x = median, y = serotype)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_pointrange(
    aes(xmin = lo, xmax = hi),
    colour = "#2c7bb6", linewidth = 0.8, size = 0.6
  ) +
  geom_text(
    aes(label = sprintf("%.2f (%.2f–%.2f)", median, lo, hi)),
    hjust = -0.15, size = 3.5
  ) +
  scale_x_continuous(
    name   = "Odds ratio (per unit SD increase in baseline MFI)",
    limits = c(0, max(or_summary$hi) * 1.5)
  ) +
  labs(y = NULL, title = "Effect of baseline titre on probability of boosting (Fiji)") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major.y = element_blank())

ggsave(here("outputs/figures/Fiji", "fiji_or_baseline.png"), or_plot, width = 8, height = 5, dpi = 300)

#-------------------------------------------------------------------------------

# 3.3. Marginal effect curve: P(booster) across baseline titre range ----

reg_plots_marginal <- list()

for (s in serotypes) {
  reg_draws <- fits_fiji[[s$name]]$draws(
    variables = c("theta_intercept", "theta_beta[1]"),
    format    = "matrix"
  )
  
  baseline_grid <- seq(
    min(fiji_inputs[[s$baseline_col]], na.rm = TRUE),
    max(fiji_inputs[[s$baseline_col]], na.rm = TRUE),
    length.out = 300
  )
  
  prob_mat <- matrix(NA_real_, nrow(reg_draws), length(baseline_grid))
  for (i in seq_len(nrow(reg_draws))) {
    lp <- reg_draws[i, "theta_intercept"] + reg_draws[i, "theta_beta[1]"] * baseline_grid
    prob_mat[i, ] <- plogis(lp)
  }
  
  bs_mean <- mean(fiji_inputs[[s$mfi_col]], na.rm = TRUE)
  bs_sd   <- sd(fiji_inputs[[s$mfi_col]],   na.rm = TRUE)
  
  reg_df <- data.frame(
    baseline_mfi = baseline_grid * bs_sd + bs_mean,
    med = apply(prob_mat, 2, median),
    lo  = apply(prob_mat, 2, quantile, 0.025),
    hi  = apply(prob_mat, 2, quantile, 0.975)
  )
  
  reg_plots_marginal[[s$name]] <- ggplot(reg_df, aes(x = baseline_mfi)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#EE6677", alpha = 0.20) +
    geom_line(aes(y = med), colour = "#EE6677", linewidth = 1.0) +
    geom_rug(data = fiji_inputs, aes(x = .data[[s$mfi_col]]),
             sides = "b", alpha = 0.3, linewidth = 0.3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = s$name, x = "Baseline MFI", y = "P(booster | baseline titre)",
         caption = "Shaded band = posterior 95% CrI. Rug = observed baseline values.") +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank())
}

marginal_joint <- (reg_plots_marginal$DENV1 + reg_plots_marginal$DENV2) /
  (reg_plots_marginal$DENV3 + reg_plots_marginal$DENV4)

ggsave(here("outputs/figures/Fiji", "fiji_marginal_baseline.png"), marginal_joint, width = 14, height = 10, dpi = 300)

#-------------------------------------------------------------------------------

# 3.4. Scatter: baseline titre vs change in MFI, coloured by P(booster) ----

reg_plots_scatter <- list()

for (s in serotypes) {
  reg_plots_scatter[[s$name]] <- ggplot(
    fiji_inputs,
    aes(x = .data[[s$mfi_col]], y = .data[[s$delta_col]], colour = .data[[s$p_booster_col]])
  ) +
    geom_point(alpha = 0.5, size = 1.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_colour_gradientn(
      colours = c("#4477AA", "#FFFFFF", "#EE6677"),
      values  = scales::rescale(c(0, 0.5, 1)),
      limits  = c(0, 1),
      name    = "P(booster)",
      labels  = scales::percent_format(accuracy = 1)
    ) +
    labs(title = s$name, x = "Baseline MFI", y = "Change in MFI (scaled)",
         caption = "Dashed line = zero change.") +
    theme_bw(base_size = 16) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
}

scatter_joint <- (reg_plots_scatter$DENV1 + reg_plots_scatter$DENV2) /
  (reg_plots_scatter$DENV3 + reg_plots_scatter$DENV4)

ggsave(here("outputs/figures/Fiji", "fiji_scatter_baseline.png"), scatter_joint, width = 14, height = 10, dpi = 300)












