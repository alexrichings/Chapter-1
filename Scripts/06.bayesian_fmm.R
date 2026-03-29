
# 06. Bayesian framework for mixture modelling
------------------------------------------------------------

# Load packages 
library(rstan)
library(ggplot2)
library(patchwork)

------------------------------------------------------------
# 1. Load data  ----

df_wide <- readRDS("Data/DR_surveillance/Serology/df_wide")
df_wide <- readRDS("Data/DR_surveillance/EN/df_wide")

df_wide$change_dengns1_4_scaled <- df_wide$change_dengns1_4 / 1000

------------------------------------------------------------
# 2. Model in Stan ----

# 2a. Model definition ----
stan_code <- "
data {
  int<lower=0> N;        
  vector[N] y;           
}
parameters {
  real mu;                    
  real<lower=0> sigma;             
  real<lower=0> alpha;             
  real<lower=0> theta;              
  real<lower=0, upper=1> lambda;   
}
model {
  mu     ~ normal(0, 1);          // a0 is fixed at 0
  sigma  ~ normal(4.13, 1);       // b0 = 4.128738
  alpha  ~ gamma(5.00, 1);        // a1 = 5.002289 (shape of gamma prior ~ your estimate)
  theta  ~ normal(4.63, 1);       // b1 = 4.634815
  lambda ~ beta(2, 47);           // lambda = 0.041 ~ beta(2, 47) gives mean ~0.04

  for (i in 1:N) {
    if (y[i] > 0) {
      target += log_mix(lambda,
                        gamma_lpdf(y[i] | alpha, 1/theta),
                        normal_lpdf(y[i] | mu, sigma));
    } else {
      target += log(1 - lambda) + normal_lpdf(y[i] | mu, sigma);
    }
  }
}
"

# 2b. Prepare data ----
stan_data <- list(
  N = nrow(df_wide),
  y = df_wide$change_dengns1_4_scaled
)

y_clean <- df_wide$change_dengns1_4_scaled[!is.na(df_wide$change_dengns1_4_scaled)]

stan_data <- list(
  N = length(y_clean),
  y = y_clean
)

# 2c. Fit model ----
fit <- stan(
  model_code = stan_code,  
  data = stan_data,
  chains = 4,
  iter = 10000,              
  warmup = 1000,
  cores = 4
)

print(fit)

------------------------------------------------------------
# 3. Visualise model ----

# 3a. Extract posteriors
posterior <- as.data.frame(fit)

set.seed(123)

# randomly sample 200 rows from the posterior 
idx <- sample(nrow(posterior), 200)

# sequences within full and positive bounds 
x_all      <- seq(min(df_wide$change_dengns1_4_scaled, na.rm = TRUE), max(df_wide$change_dengns1_4_scaled, na.rm = TRUE), length.out = 1000)
x_positive <- seq(0.001, max(df_wide$change_dengns1_4_scaled, na.rm = TRUE), length.out = 1000)

# pre-define a matrix of 200 cols for the random samples of the posteriors 
mixture_samples <- matrix(NA, nrow = length(x_all), ncol = 200)
prob_samples    <- matrix(NA, nrow = length(x_positive), ncol = 200)

# 
for (j in seq_along(idx)) {
  i        <- idx[j]
  mu_i     <- posterior$mu[i]
  sigma_i  <- posterior$sigma[i]
  alpha_i  <- posterior$alpha[i]
  theta_i  <- posterior$theta[i]
  lambda_i <- posterior$lambda[i]
  
  # mixture density over full range
  mixture_samples[, j] <- 
    (1 - lambda_i) * dnorm(x_all, mean = mu_i, sd = sigma_i) +
    lambda_i * ifelse(x_all > 0, dgamma(x_all, shape = alpha_i, scale = theta_i), 0)
  
  # p(Gamma | y) over positive range
  numerator        <- lambda_i * dgamma(x_positive, shape = alpha_i, scale = theta_i)
  denominator      <- numerator + (1 - lambda_i) * dnorm(x_positive, mean = mu_i, sd = sigma_i)
  prob_samples[, j] <- numerator / denominator
}

# build dataframes
mixture_lines_df <- data.frame(
  x       = rep(x_all, 200),
  density = as.vector(mixture_samples),
  sample  = rep(1:200, each = length(x_all))
)

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

# plot 1 - mixture density
p1a <- ggplot() +
  geom_histogram(data = df_wide, aes(x = change_dengns1_4_scaled, y = after_stat(density)),
                 binwidth = 4, fill = "lightgrey", colour = "white") +
  geom_line(data = mixture_lines_df, aes(x = x, y = density, group = sample),
            colour = "black", alpha = 0.05) +
  labs(x = "MIA change", y = "Density",
       title = "Mixture model fit with uncertainty") +
  theme_minimal() 

p1b <- ggplot() +
  geom_histogram(data = df_wide, aes(x = change_dengns1_4_scaled, y = after_stat(density)),
                 binwidth = 4, fill = "lightgrey", colour = "white") +
  geom_ribbon(data = mixture_ribbon_df, aes(x = x, ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(data = mixture_ribbon_df, aes(x = x, y = mean), colour = "blue", linewidth = 1) +
  labs(x = "MIA change", y = "Density",
       title = "Mixture model fit with uncertainty") +
  theme_minimal()

# plot 2 - posterior probability of boosting
p2 <- ggplot(prob_df, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean), colour = "blue", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  labs(x = "MIA", y = "P(Gamma | MIA)",
       title = "Posterior probability of boosting component") +
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0,1))

# combine
p1a / p2
p1b / p2

------------------------------------------------------------
# 4. Regression analyses ----



