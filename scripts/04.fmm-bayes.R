
# 04. Bayesian framework for mixture modelling
------------------------------------------------------------

# 1. Install and load packages ----

library(rstan)
library(ggplot2)
library(patchwork)
library(here)
library(tidyverse)

------------------------------------------------------------
  
# 2. Load data ----

inputs <- readRDS(here("data", "inputs.rds"))

------------------------------------------------------------

# 3. Model set up ----

# Define model 

stan_code <- "

data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real<lower=0> sigma;  // standard deviation of normal component         
  real<lower=0> alpha;  // shape of gamma component   
  real<lower=0> b1;  // scale of gamma component
  real<lower=0, upper=1> lambda;  // proportion boosting in population 
}

model {
  // Priors
  sigma  ~ lognormal(1.6, 0.5);  // sd of non-boost: most mass between 0-10, sensible for scaled data
  alpha ~ gamma(10, 1);  // shape of boost: large alpha represents near symmetry in boosts 
  b1 ~ cauchy(0, 5);  // scale of boost: consistent with boosts up to around 40 
  lambda ~ beta(2, 23);  // mean around 0.08, reflect the FOI from serocatalytic model 
  
  // Likelihood
  for (i in 1:N) {
    if (y[i] > 0) {
      target += log_mix(lambda,
                        gamma_lpdf(y[i] | alpha, 1/b1),  
                        normal_lpdf(y[i] | 0, sigma));      
    } else {
      target += log(1 - lambda) + normal_lpdf(y[i] | 0, sigma);
    }
  }
}
"

# 2b. Prepare data ----
deng1_clean <- inputs$change_dengns1_1_scale[!is.na(inputs$change_dengns1_1_scale)]
deng2_clean <- inputs$change_dengns1_2_scale[!is.na(inputs$change_dengns1_2_scale)]
deng3_clean <- inputs$change_dengns1_3_scale[!is.na(inputs$change_dengns1_3_scale)]
deng4_clean <- inputs$change_dengns1_4_scale[!is.na(inputs$change_dengns1_4_scale)]

stan_deng1 <- list(N = length(deng1_clean), y = deng1_clean)
stan_deng2 <- list(N = length(deng2_clean), y = deng2_clean)
stan_deng3 <- list(N = length(deng3_clean), y = deng3_clean)
stan_deng4 <- list(N = length(deng4_clean), y = deng4_clean)

# 2c. Fit model ----
deng1_fit <- stan(model_code = stan_code, data = stan_deng1, chains = 4, iter = 10000, warmup = 1000, cores = 4)
deng2_fit <- stan(model_code = stan_code, data = stan_deng2, chains = 4, iter = 10000, warmup = 1000, cores = 4)
deng3_fit <- stan(model_code = stan_code, data = stan_deng3, chains = 4, iter = 10000, warmup = 1000, cores = 4)
deng4_fit <- stan(model_code = stan_code, data = stan_deng4, chains = 4, iter = 10000, warmup = 1000, cores = 4)

print(deng1_fit)
print(deng2_fit)
print(deng3_fit)
print(deng4_fit)

------------------------------------------------------------

# 3. Assess model output  ----

# 3a. Extract posteriors
post_deng1 <- as.data.frame(deng1_fit)
post_deng2 <- as.data.frame(deng2_fit)
post_deng3 <- as.data.frame(deng3_fit)
post_deng4 <- as.data.frame(deng4_fit)

# set seed 
set.seed(123)

# randomly sample 200 rows fromo posterior 
idx_deng1 <- sample(nrow(post_deng1), 200)
idx_deng2 <- sample(nrow(post_deng2), 200)
idx_deng3 <- sample(nrow(post_deng3), 200)
idx_deng4 <- sample(nrow(post_deng4), 200)

# sequences within full and positive bounds
x_all_deng1 <- seq(min(deng1_clean, na.rm = TRUE), max(deng1_clean, na.rm = TRUE), length.out = 1000)
x_all_deng2 <- seq(min(deng2_clean, na.rm = TRUE), max(deng2_clean, na.rm = TRUE), length.out = 1000)
x_all_deng3 <- seq(min(deng3_clean, na.rm = TRUE), max(deng3_clean, na.rm = TRUE), length.out = 1000)
x_all_deng4 <- seq(min(deng4_clean, na.rm = TRUE), max(deng4_clean, na.rm = TRUE), length.out = 1000)

x_pos_deng1 <- seq(0.001, max(deng1_clean, na.rm = TRUE), length.out = 1000)
x_pos_deng2 <- seq(0.001, max(deng2_clean, na.rm = TRUE), length.out = 1000)
x_pos_deng3 <- seq(0.001, max(deng3_clean, na.rm = TRUE), length.out = 1000)
x_pos_deng4 <- seq(0.001, max(deng4_clean, na.rm = TRUE), length.out = 1000)

# pre-define matrices for posterior samples
mix_deng1 <- matrix(NA, nrow = length(x_all_deng1), ncol = 200)
mix_deng2 <- matrix(NA, nrow = length(x_all_deng2), ncol = 200)
mix_deng3 <- matrix(NA, nrow = length(x_all_deng3), ncol = 200)
mix_deng4 <- matrix(NA, nrow = length(x_all_deng4), ncol = 200)

prob_deng1 <- matrix(NA, nrow = length(x_pos_deng1), ncol = 200)
prob_deng2 <- matrix(NA, nrow = length(x_pos_deng2), ncol = 200)
prob_deng3 <- matrix(NA, nrow = length(x_pos_deng3), ncol = 200)
prob_deng4 <- matrix(NA, nrow = length(x_pos_deng4), ncol = 200)


# estimate probability of boosting 
idx <- idx_deng1

for (j in seq_along(idx)) {
  i        <- idx[j]
  sigma_i  <- posterior$sigma[i]
  alpha_i  <- posterior$alpha[i]
  b1_i     <- posterior$b1[i]
  lambda_i <- posterior$lambda[i]
  
  # mixture density over full range 
  mixture_samples[, j] <- (1 - lambda_i) * dnorm(x_all, mean = 0, sd = sigma_i) + lambda_i * ifelse(x_all > 0, dgamma(x_all, shape = alpha_i, scale = b1_i), 0)
  
  # p(Gamma | y) over positive range
  numerator         <- lambda_i * dgamma(x_positive, shape = alpha_i, scale = b1_i)
  denominator       <- numerator + (1 - lambda_i) * dnorm(x_positive, mean = 0, sd = sigma_i)
  prob_samples[, j] <- numerator / denominator
}

# build dataframes
mixture_ribbon_df <- data.frame(
  x     = x_all,
  mean  = rowMeans(mixture_samples),
  lower = apply(mixture_samples, 1, quantile, 0.025),
  upper = apply(mixture_samples, 1, quantile, 0.975)
)

prob_df <- data.frame(
  x     = x_positive,
  mean  = rowMeans(prob_samples),
  lower = apply(prob_samples, 1, quantile, 0.025),
  upper = apply(prob_samples, 1, quantile, 0.975)
)

# scaling factor to map probability (0-1) onto density axis
scale_factor <- max(mixture_ribbon_df$upper) / 1

# combined plot
p_combined <- ggplot() +
  # histogram of observed data
  geom_histogram(
    data = inputs,
    aes(x = change_dengns1_1_scale, y = after_stat(density)),  
    binwidth = 2.5, fill = "lightgrey", colour = "white"
  ) +
  # mixture density ribbon and mean line
  geom_ribbon(
    data = mixture_ribbon_df,
    aes(x = x, ymin = lower, ymax = upper),
    fill = "blue", alpha = 0.2
  ) +
  geom_line(
    data = mixture_ribbon_df,
    aes(x = x, y = mean),
    colour = "blue", linewidth = 1
  ) +
  # posterior probability of boosting scaled to density axis
  geom_ribbon(
    data = prob_df,
    aes(x = x, ymin = lower * scale_factor, ymax = upper * scale_factor),
    fill = "red", alpha = 0.2
  ) +
  geom_line(
    data = prob_df,
    aes(x = x, y = mean * scale_factor),
    colour = "red", linewidth = 1
  ) +
  # 0.5 threshold line scaled to density axis
  geom_hline(
    yintercept = 0.5 * scale_factor,
    linetype = "dashed", colour = "red"
  ) +
  # secondary axis to show probability scale on right
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / scale_factor, name = "P(Boost | y)")
  ) +
  labs(
    x = "Antibody change (scaled)",
    title = "Mixture model fit and posterior probability of boosting"
  ) +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(colour = "red"),
    axis.text.y.right  = element_text(colour = "red")
  )

p_combined




# extract summary from stan fit
fit_summary <- summary(fit)$summary

# pull the parameters of interest and format
param_table <- fit_summary[c("sigma", "alpha", "b1", "lambda"), 
                           c("mean", "2.5%", "97.5%")] %>%
  as.data.frame() %>%
  rownames_to_column("Parameter") %>%
  rename(
    Mean   = mean,
    Lower  = `2.5%`,
    Upper  = `97.5%`
  ) %>%
  mutate(
    `95% CrI`   = paste0("[", round(Lower, 3), ", ", round(Upper, 3), "]"),
    Mean        = round(Mean, 3),
    Parameter   = c("Sigma (sd, normal)", "Alpha (shape, gamma)", 
                    "B1 (scale, gamma)", "Lambda (prop. boosting)")
  ) %>%
  select(Parameter, Mean, `95% CrI`)

# print as a clean table
knitr::kable(param_table, align = c("l", "c", "c"))

