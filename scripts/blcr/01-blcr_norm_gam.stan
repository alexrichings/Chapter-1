
// gaussian_gamma_blcr.stan
data {
  int<lower=0> N;
  int<lower=1> K;          // number of covariates
  vector[N] delta;
  matrix[N, K] X;          // covariate matrix: baseline + any others
}

parameters {
  real mu1;
  real<lower=0> sigma1;
  real log_alpha_m1;       // unconstrained; alpha = 1 + exp(log_alpha_m1)
  real<lower=0> beta_rate;
  real theta_intercept;
  vector[K] theta_beta;
}

transformed parameters {
  real alpha = 1.0 + exp(log_alpha_m1);   // always > 1, smooth
  vector[N] log_theta_booster;
  vector[N] log_theta_nonbooster;
  {
    vector[N] lp = theta_intercept + X * theta_beta;
    log_theta_booster    = -log1p_exp(-lp);  // log P(booster)
    log_theta_nonbooster = -log1p_exp(lp);   // log P(non-booster)
  }
}

model {
  mu1             ~ normal(0, 1);
  sigma1          ~ normal(0, 1);
  log_alpha_m1    ~ normal(0, 1);
  beta_rate ~ gamma(2, 9.0);   // prior mean rate ≈ 0.22 → Gamma mean ≈ 5/0.22 ≈ 22.5
  theta_intercept ~ normal(0, 1.5);  // weakly informative on booster proportion
  theta_beta      ~ normal(0, 1);    // one prior applies to all covariate coefficients
  
  for (i in 1:N) {
    real lp_nonbooster = log_theta_nonbooster[i] + normal_lpdf(delta[i] | mu1, sigma1);
    if (delta[i] > 0) {
      real lp_booster = log_theta_booster[i] + gamma_lpdf(delta[i] | alpha, beta_rate);
      target += log_sum_exp(lp_nonbooster, lp_booster);
    } else {
      target += lp_nonbooster;
    }
  }
}

generated quantities {
  vector[N] p_booster;
  real prop_booster;
  for (i in 1:N) {
    if (delta[i] > 0) {
      real lp_nb = log_theta_nonbooster[i] + normal_lpdf(delta[i] | mu1, sigma1);
      real lp_b  = log_theta_booster[i]    + gamma_lpdf(delta[i] | alpha, beta_rate);
      p_booster[i] = exp(lp_b - log_sum_exp(lp_nb, lp_b));
    } else {
      p_booster[i] = 0;
    }
  }
  prop_booster = mean(inv_logit(theta_intercept + X * theta_beta));
}



