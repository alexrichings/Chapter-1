
# -------------------------------------------------------------------------------
# Script: 08-simulated-parameters.R
# -------------------------------------------------------------------------------

# Objective:
#   Simulate a population with known individual-level risk factors and exposure
#   parameters. Apply the FMM to recovered ΔMFI data and check whether the
#   model recovers:
#     (1) the FMM mixture parameters (λ, b0, a1, b1)
#     (2) the risk factor ORs via logistic regression on FMM-classified boost
#         status

# Output:
#   Console: parameter recovery table, AUC, sensitivity/specificity, OR table
#   Plots:fitted FMM overlay on simulated histogram; ROC curve

# -------------------------------------------------------------------------------

# Load packages ----

library(tidyverse)
library(pROC)

# -------------------------------------------------------------------------------
# Section 1: Simulation parameters
# -------------------------------------------------------------------------------

# 1.1. Global settings ----

set.seed(42)
n <- 1000

# -------------------------------------------------------------------------------

# 1.2. Risk factor prevalence ----

prev_young <- 0.40 # p(age < 20)
prev_bucket <- 0.30 # p(bucket water storage)
prev_repellant <- 0.50 # p(insect repellant usage)

# -------------------------------------------------------------------------------

# 1.3. Baseline exposure probability and risk factor log-ORs ----

p_baseline <- 0.20 # baseline exposure probability 
logit_baseline <- qlogis(p_baseline) # convert to log odds 

log_or_young     <- log(2.0) # OR = 2.0 — age < 20 vs ≥ 20
log_or_bucket    <- log(1.4) # OR = 1.4 — bucket storage vs not
log_or_repellant <- log(0.6) # OR = 0.6 — repellant use vs not (protective)

# -------------------------------------------------------------------------------
# 1.4. FMM distribution parameters ----

# Non-boosters: Normal(mean = 0, sd = b0)
# Boosters: Gamma(shape = a1, scale = b1)
# Gamma mean: a1 × b1 = 25.0
# SNR: (a1 × b1) / b0 = 25.0 / 4.5 ≈ 5.6 (signal-to-noise ratio)

# a larger SNR means the distributions separate better
# either because there is less noise or a bigger boost 

true_b0 <- 4.5 # noise SD 
true_a1 <- 5.0 # gamma shape
true_b1 <- 5.0 # gamma scale


# -------------------------------------------------------------------------------
# Section 2: Simulate population
# -------------------------------------------------------------------------------

# 2.1. Draw individual risk factors ----

# risk factors drawn from a Bernoulli distribution 
pop <- tibble(
  id = seq_len(n),
  young = rbinom(n, 1, prev_young),
  bucket = rbinom(n, 1, prev_bucket),
  repellant = rbinom(n, 1, prev_repellant)
)

# -------------------------------------------------------------------------------
# 2.2. Compute individual exposure probability via logistic model ----

# determine the log odds of exposure from the baseline exposure and drivers
pop <- pop |>
  mutate(
    logit_p = logit_baseline +
                log_or_young * young +
                log_or_bucket * bucket +
                log_or_repellant * repellant,
    p_exposed = plogis(logit_p) # convert log odds to probability 
  )

# -------------------------------------------------------------------------------
# 2.3. Draw exposure status ----

# allocate exposures in proportion to p(exposure)
pop <- pop |>
  mutate(exposed = rbinom(n, 1, p_exposed))

cat("Realised boost prevalence:         ", round(mean(pop$exposed), 3), "\n")
cat("Expected boost prevalence (E[p_i]):", round(mean(pop$p_exposed), 3), "\n")

# -------------------------------------------------------------------------------
# 2.4. Draw ΔMFI conditional on exposure status ----

# if exposed draw from the gamma distribution  
# otherwise draw from the normal distribution 
pop <- pop |>
  mutate(
    delta_mfi = if_else(
      exposed == 1,
      rgamma(n, shape = true_a1, scale = true_b1), # boosters
      rnorm(n,  mean  = 0, sd = true_b0) # non-boosters
    )
  )


# -------------------------------------------------------------------------------
# Section 3: FMM utility functions
# -------------------------------------------------------------------------------

# 3.1. Logistic transform for mixing proportion ----

transform.lambda <- function(x) { 1 / (1 + exp(-x)) }

# -------------------------------------------------------------------------------
# 3.2. Negative log-likelihood — Normal(0, b0) + Gamma(a1, b1) mixture
#
#   For val > 0: full mixture density
#   For val ≤ 0: normal component only (gamma has no support on negatives)

fitmixture <- function(param, val) {
  lambda1 <- transform.lambda(param[["lambda1"]])
  b0 <- pmax(param[["b0"]], 1e-5)
  a1 <- pmax(param[["a1"]], 1e-5)
  b1 <- pmax(param[["b1"]], 1e-5)

  non_boost <- dnorm(val, mean = 0, sd = b0)
  boost <- dgamma(pmax(val, 1e-10), shape = a1, scale = b1)

  lik <- if_else(
    val > 0,
    (1 - lambda1) * non_boost + lambda1 * boost,
    (1 - lambda1) * non_boost
  )

  -sum(log(lik + 1e-10))
}

# -------------------------------------------------------------------------------
# 3.3. Fit wrapper — L-BFGS-B optimisation ----

fit_pathogen <- function(data, start_par, lower = NULL) {
  data      <- data[!is.na(data)]
  start_par <- start_par[names(start_par) != "a0"]

  if (is.null(lower)) {
    lower <- c(b0 = 1e-5, a1 = 1e-5, b1 = 1e-5, lambda1 = -10)
  }

  optim(
    par     = start_par,
    fn      = fitmixture,
    method  = "L-BFGS-B",
    val     = data,
    lower   = lower,
    hessian = FALSE
  )
}

# -------------------------------------------------------------------------------
# 3.4. Extract parameters ----

extract_params <- function(opt_par) {
  list(
    lambda1 = transform.lambda(opt_par["lambda1"]),
    b0      = opt_par["b0"],
    a1      = opt_par["a1"],
    b1      = opt_par["b1"]
  )
}


# -------------------------------------------------------------------------------
# Section 4: Fit FMM to simulated ΔMFI
# -------------------------------------------------------------------------------

# 4.1. Fit model ----

fmm_fit <- fit_pathogen(
  pop$delta_mfi,
  start_par = c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)
)

est_par <- extract_params(fmm_fit$par)

cat("\nOptimiser convergence code (0 = converged):", fmm_fit$convergence, "\n")

# -------------------------------------------------------------------------------
# 4.2. Parameter recovery table ----

recovery <- tibble(
  parameter = c("λ (boost prevalence)", "b0 (noise SD)",
                "a1 (gamma shape)", "b1 (gamma scale)"),
  true      = c(mean(pop$exposed), true_b0, true_a1, true_b1),
  estimated = c(est_par$lambda1, est_par$b0, est_par$a1, est_par$b1)
) |>
  mutate(rel_error_pct = round((estimated - true) / true * 100, 1))

cat("\n---- Parameter recovery ----\n")
print(recovery)

# trade-off between a1 and b1
# the same gamma mean may come from different combinations of parameters 
# this may influence the distribution and thus the overlap between distributions


# -------------------------------------------------------------------------------
# Section 5: Classify individuals and evaluate
# -------------------------------------------------------------------------------

# 5.1. Posterior probability of boosting ----

post_prob <- function(delta, pars) {
  normal_comp <- (1 - pars$lambda1) * dnorm(delta, mean = 0, sd = pars$b0)
  gamma_comp <- pars$lambda1  * dgamma(pmax(delta, 1e-10), shape = pars$a1, scale = pars$b1)
  marginal <- normal_comp + gamma_comp
  gamma_comp / marginal
}

# boosters are defined at the p = 0.5 threshold 
pop <- pop |>
  mutate(
    p_boost = post_prob(delta_mfi, est_par),
    boost_class = as.integer(p_boost >= 0.5)
  )

# -------------------------------------------------------------------------------
# 5.2. AUC and threshold-level performance ----

roc_obj <- pROC::roc(pop$exposed, pop$p_boost, quiet = TRUE)

cat("\n---- Classification performance ----\n")
cat("AUC:         ", round(as.numeric(pROC::auc(roc_obj)), 3), "\n")
cat("Sensitivity (P ≥ 0.5):",
    round(sum(pop$boost_class == 1 & pop$exposed == 1) / sum(pop$exposed == 1), 3), "\n")
cat("Specificity (P ≥ 0.5):",
    round(sum(pop$boost_class == 0 & pop$exposed == 0) / sum(pop$exposed == 0), 3), "\n")


# -------------------------------------------------------------------------------
# Section 6: Risk factor recovery via logistic regression
# -------------------------------------------------------------------------------

# 6.1. Fit logistic regression ----

# Regression for model-derived exposures 
lr_fmm1 <- glm(boost_class ~ young, data = pop, family = binomial(link = "logit"))
lr_fmm2 <- glm(boost_class ~ bucket, data = pop, family = binomial(link = "logit"))
lr_fmm3 <- glm(boost_class ~ repellant, data = pop, family = binomial(link = "logit"))

# Regression for true exposures 
lr_true1 <- glm(exposed ~ young,     data = pop, family = binomial(link = "logit"))
lr_true2 <- glm(exposed ~ bucket,    data = pop, family = binomial(link = "logit"))
lr_true3 <- glm(exposed ~ repellant, data = pop, family = binomial(link = "logit"))

# -------------------------------------------------------------------------------
# 6.2. OR recovery table ----

# Use univariable regression as you are interested in the individual OR 
get_or_univar <- function(model, varname) {
  coef_val <- coef(model)[varname]
  ci       <- suppressMessages(confint(model)[varname, ])
  tibble(
    or      = round(exp(coef_val), 2),
    ci_low  = round(exp(ci[1]),    2),
    ci_high = round(exp(ci[2]),    2)
  )
}

or_recovery <- bind_rows(
  tibble(source = "True (oracle)",  risk_factor = "young",
         true_or = 2.0, get_or_univar(lr_true1, "young")),
  tibble(source = "True (oracle)",  risk_factor = "bucket",
         true_or = 1.4, get_or_univar(lr_true2, "bucket")),
  tibble(source = "True (oracle)",  risk_factor = "repellant",
         true_or = 0.6, get_or_univar(lr_true3, "repellant")),
  tibble(source = "FMM-classified", risk_factor = "young",
         true_or = 2.0, get_or_univar(lr_fmm1, "young")),
  tibble(source = "FMM-classified", risk_factor = "bucket",
         true_or = 1.4, get_or_univar(lr_fmm2, "bucket")),
  tibble(source = "FMM-classified", risk_factor = "repellant",
         true_or = 0.6, get_or_univar(lr_fmm3, "repellant"))
)

cat("\n---- OR recovery ----\n")
print(or_recovery)

# -------------------------------------------------------------------------------
# Section 7: Visualisation
# -------------------------------------------------------------------------------

# 7.1. Fitted FMM overlaid on simulated ΔMFI histogram ----

mixture_density <- function(x, pars) {
  normal_comp <- (1 - pars$lambda1) * dnorm(x, mean = 0, sd = pars$b0)
  gamma_comp  <-      pars$lambda1  * dgamma(pmax(x, 1e-10),
                                             shape = pars$a1, scale = pars$b1)
  tibble(
    x       = x,
    Mixture = normal_comp + gamma_comp,
    Normal  = normal_comp,
    Gamma   = gamma_comp
  )
}

x_grid  <- seq(min(pop$delta_mfi), max(pop$delta_mfi), length.out = 500)
dens_df <- mixture_density(x_grid, est_par) |>
  pivot_longer(-x, names_to = "component", values_to = "density")

p_fmm <- ggplot() +
  geom_histogram(
    data = pop,
    aes(x = delta_mfi, y = after_stat(density)),
    fill = "lightgrey", colour = "black", binwidth = 2
  ) +
  geom_line(
    data = dens_df |> filter(component != "Mixture"),
    aes(x = x, y = density, colour = component, linetype = component),
    linewidth = 1.2
  ) +
  geom_line(
    data = dens_df |> filter(component == "Mixture"),
    aes(x = x, y = density, colour = component, linetype = component),
    linewidth = 1.2
  ) +
  scale_colour_manual(
    breaks = c("Mixture", "Normal", "Gamma"),
    values = c(Mixture = "black", Normal = "red", Gamma = "purple"),
    labels = c("Mixture", "Normal (non-boost)", "Gamma (boost)")
  ) +
  scale_linetype_manual(
    breaks = c("Mixture", "Normal", "Gamma"),
    values = c(Mixture = "solid", Normal = "solid", Gamma = "solid"),
    labels = c("Mixture", "Normal (non-boost)", "Gamma (boost)")
  ) +
  labs(
    x        = "ΔMFi (scaled, /1000)",
    y        = "Density",
    title    = "Simulated data with fitted FMM",
    colour   = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# -------------------------------------------------------------------------------
# 7.2. ROC curve ----

roc_df <- tibble(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 1.2, colour = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  annotate(
    "text", x = 0.6, y = 0.2,
    label = paste0("AUC = ", round(as.numeric(pROC::auc(roc_obj)), 3)),
    size  = 5
  ) +
  labs(
    x     = "1 − Specificity",
    y     = "Sensitivity",
    title = "ROC curve: FMM boost classification vs true exposure"
  ) +
  theme_minimal(base_size = 14)

print(p_fmm)
print(p_roc)

cat("\n---- Done ----\n")


# -------------------------------------------------------------------------------
# Section 8: Repeated simulation: single scenario  ----
# -------------------------------------------------------------------------------

# 8.1. Single-run function ----

# Wraps the full pipeline: simulate population > fit FMM > classify > 
# fit univariable logistic regressions > return OR estimates and CIs
# Returns NA rows on optimiser failure so the loop doesn't break 

run_one <- function(seed) {
  
  set.seed(seed)
  
  # Simulate population
  pop <- tibble(
    id        = seq_len(n),
    young     = rbinom(n, 1, prev_young),
    bucket    = rbinom(n, 1, prev_bucket),
    repellant = rbinom(n, 1, prev_repellant)
  ) |>
    mutate(
      logit_p   = logit_baseline +
        log_or_young * young +
        log_or_bucket * bucket +
        log_or_repellant * repellant,
      p_exposed = plogis(logit_p),
      exposed   = rbinom(n, 1, p_exposed),
      delta_mfi = if_else(
        exposed == 1,
        rgamma(n, shape = true_a1, scale = true_b1),
        rnorm(n,  mean  = 0,       sd    = true_b0)
      )
    )
  
  # Fit FMM
  fit <- tryCatch(
    fit_pathogen(pop$delta_mfi, start_par = c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0) return(NULL)
  
  est_par <- extract_params(fit$par)
  
  # Classify
  pop <- pop |>
    mutate(
      p_boost     = post_prob(delta_mfi, est_par),
      boost_class = as.integer(p_boost >= 0.5)
    )
  
  # Univariable logistic regressions — FMM-classified and oracle
  get_or_univar <- function(model, varname) {
    coef_val <- coef(model)[varname]
    ci       <- suppressMessages(confint(model)[varname, ])
    tibble(or = exp(coef_val), ci_low = exp(ci[1]), ci_high = exp(ci[2]))
  }
  
  vars <- c("young", "bucket", "repellant")
  
  fmm_ors <- map_dfr(vars, function(v) {
    m <- glm(as.formula(paste("boost_class ~", v)), data = pop, family = binomial())
    get_or_univar(m, v) |> mutate(risk_factor = v, source = "fmm")
  })
  
  oracle_ors <- map_dfr(vars, function(v) {
    m <- glm(as.formula(paste("exposed ~", v)), data = pop, family = binomial())
    get_or_univar(m, v) |> mutate(risk_factor = v, source = "oracle")
  })
  
  bind_rows(fmm_ors, oracle_ors) |> mutate(seed = seed)
}

# -------------------------------------------------------------------------------
# 8.2. Run 500 iterations ----

n_sims <- 500
sim_results <- map_dfr(seq_len(n_sims), run_one, .progress = TRUE)

# -------------------------------------------------------------------------------
# 8.3. Compute bias, coverage and RMSE ----

true_ors <- c(young = 2.0, bucket = 1.4, repellant = 0.6)

summary_stats <- sim_results |>
  mutate(
    true_or   = true_ors[risk_factor],
    log_or    = log(or),
    log_true  = log(true_or),
    covers    = ci_low <= true_or & ci_high >= true_or,
    sq_error  = (log_or - log_true)^2
  ) |>
  group_by(source, risk_factor) |>
  summarise(
    bias     = round(mean(log_or) - first(log_true), 3),
    coverage = round(mean(covers), 3),
    rmse     = round(sqrt(mean(sq_error)), 3),
    .groups  = "drop"
  )

cat("\n---- Simulation summary (n_sims =", n_sims, ") ----\n")
print(summary_stats)


# -------------------------------------------------------------------------------
# Section 9: Grid simulation ----

# Want to compare the performance of different SNRs and boost prevalence
# -------------------------------------------------------------------------------

# 9.1. Parameter grid

# SNR = (a1 × b1) / b0; b0 and a1 held fixed, b1 varies to hit target SNR.
# λ varied by scaling the baseline logit to hit target mean exposure probability.

snr_values <- c(2, 4, 6, 8)
lambda_values <- c(0.05, 0.10, 0.20, 0.40)
n_sims <- 500

param_grid <- expand_grid(snr = snr_values, lambda = lambda_values)

# -------------------------------------------------------------------------------
# 9.2. Single-run function ----

run_one_grid <- function(seed, snr, lambda) {
  
  set.seed(seed)
  
  # Derive b1 from target SNR (b0 and a1 fixed)
  b1_run <- (snr * true_b0) / true_a1
  
  # Derive logit baseline to hit target mean exposure probability
  logit_base_run <- qlogis(lambda)
  
  # Simulate population
  pop <- tibble(
    id        = seq_len(n),
    young     = rbinom(n, 1, prev_young),
    bucket    = rbinom(n, 1, prev_bucket),
    repellant = rbinom(n, 1, prev_repellant)
  ) |>
    mutate(
      logit_p   = logit_base_run +
        log_or_young     * young +
        log_or_bucket    * bucket +
        log_or_repellant * repellant,
      p_exposed = plogis(logit_p),
      exposed   = rbinom(n, 1, p_exposed),
      delta_mfi = if_else(
        exposed == 1,
        rgamma(n, shape = true_a1, scale = b1_run),
        rnorm(n,  mean  = 0,       sd    = true_b0)
      )
    )
  
  # Fit FMM — return NULL on failure
  fit <- tryCatch(
    fit_pathogen(pop$delta_mfi, start_par = c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0) return(NULL)
  
  est_par <- extract_params(fit$par)
  
  # Classify
  pop <- pop |>
    mutate(
      p_boost     = post_prob(delta_mfi, est_par),
      boost_class = as.integer(p_boost >= 0.5)
    )
  
  # Univariable logistic regressions
  get_or_univar <- function(model, varname) {
    coef_val <- coef(model)[varname]
    ci <- tryCatch(
      suppressMessages(confint(model)[varname, ]),
      error = function(e) {
        se <- sqrt(diag(vcov(model)))[varname]
        c(coef_val - 1.96 * se, coef_val + 1.96 * se)
      }
    )
    tibble(or = exp(coef_val), ci_low = exp(ci[1]), ci_high = exp(ci[2]))
  }
  
  vars <- c("young", "bucket", "repellant")
  
  fmm_ors <- map_dfr(vars, function(v) {
    m <- glm(as.formula(paste("boost_class ~", v)), data = pop, family = binomial())
    get_or_univar(m, v) |> mutate(risk_factor = v, source = "fmm")
  })
  
  oracle_ors <- map_dfr(vars, function(v) {
    m <- glm(as.formula(paste("exposed ~", v)), data = pop, family = binomial())
    get_or_univar(m, v) |> mutate(risk_factor = v, source = "oracle")
  })
  
  bind_rows(fmm_ors, oracle_ors) |>
    mutate(seed = seed)
}

# -------------------------------------------------------------------------------
# 9.3. Run grid ----

# 500 iterations per cell
grid_results <- param_grid |>
  mutate(
    sims = map2(snr, lambda, function(s, l) {
      cat("Running SNR =", s, "| lambda =", l, "\n")
      map_dfr(seq_len(n_sims), run_one_grid, snr = s, lambda = l)
    })
  ) |>
  unnest(sims)

# save results 
saveRDS(grid_results, "outputs/simulation/grid_results.rds")

# -------------------------------------------------------------------------------
# 9.4. Compute bias, coverage and RMSE per grid cell

true_ors <- c(young = 2.0, bucket = 1.4, repellant = 0.6)

grid_summary <- grid_results |>
  mutate(
    true_or  = true_ors[risk_factor],
    log_or   = log(or),
    log_true = log(true_or),
    covers   = ci_low <= true_or & ci_high >= true_or,
    sq_error = (log_or - log_true)^2
  ) |>
  group_by(snr, lambda, source, risk_factor) |>
  summarise(
    n_converged = n(),
    bias        = round(mean(log_or, na.rm = TRUE) - first(log_true), 3),
    coverage    = round(mean(covers, na.rm = TRUE), 3),
    rmse        = round(sqrt(mean(sq_error, na.rm = TRUE)), 3),
    .groups     = "drop"
  )

cat("\n---- Grid simulation summary ----\n")
print(grid_summary, n = Inf)


# -------------------------------------------------------------------------------
# Section 10: Visualise grid results
# -------------------------------------------------------------------------------

# 10.1. Coverage heatmap ----

grid_summary |>
  filter(source == "fmm") |>
  ggplot(aes(x = factor(snr), y = factor(lambda), fill = coverage)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", coverage)), size = 2.8) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#fee090",
    high     = "#1a9850",
    midpoint = 0.95,
    limits   = c(0.6, 1.0),
    name     = "Coverage"
  ) +
  facet_wrap(~risk_factor) +
  labs(
    x     = "SNR",
    y     = "λ (proportion with antibody boost)",
    title = "95% coverage heatmap across SNR and λ "
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background  = element_rect(fill = "grey85", colour = "black"),
    panel.border      = element_rect(colour = "black", fill = NA)
  )

# -------------------------------------------------------------------------------

# 10.2. Bias heatmap — FMM only
  
grid_summary |>
  filter(source == "fmm") |>
  ggplot(aes(x = factor(snr), y = factor(lambda), fill = bias)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", bias)), size = 2.8) +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "white",
    high     = "#2166ac",
    midpoint = 0,
    name     = "Bias\n(log-OR)"
  ) +
  facet_wrap(~risk_factor) +
  labs(
    x     = "SNR",
    y     = "λ (proportion with antibody boost)",
    title = "Bias in log-OR estimates across SNR and λ"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey85", colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA)
  )

# -------------------------------------------------------------------------------
  
# 10.3. RMSE heatmap — FMM only

grid_summary |>
  filter(source == "fmm") |>
  ggplot(aes(x = factor(snr), y = factor(lambda), fill = rmse)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rmse)), size = 2.8) +
  scale_fill_distiller(palette = "Blues", direction = 1, name = "RMSE\n(log-OR)") +
  facet_wrap(~risk_factor) +
  labs(
    x     = "SNR",
    y     = "λ (proportion with antibody boost)",
    title = "RMSE of log-OR estimates across SNR and λ"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey85", colour = "black"),
    panel.border     = element_rect(colour = "black", fill = NA)
  )
