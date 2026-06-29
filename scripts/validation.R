
# simulation
------------------------------------------------------------------------------------
  
# Aim is to simulate a data set with known exposures, underlying probability of infection etc.
  
library(ggplot2)
library(dplyr)
library(tidyverse)
library(flextable)

------------------------------------------------------------------------------------
# 1. Simulation Conditions ---- 
------------------------------------------------------------------------------------
  
# basing the model on exposure to mosquitoes 
# 1a. Set up underlying conditions ----

# set seed 
set.seed(123)

# set up 
n <- 1000
p_mosq <- 0.9

p0 <- 0.08 # baseline probability of boosting (no mosquito exposure)
OR <- 3 # odds ratio for boosting after exposure to mosquitoes, relative to boosting with no exposure 

# MIA data parameters 
# assume the data comes from a mixture model made from a normal and gamma distribution 
noise_mean <- 0
noise_sd <- 0.3
gamma_shape <- 3
gamma_scale <- 0.8


# 1b. Set up mosquito exposure ----

# draw n people from a binomial distribution 
# each person gets one trial 
# with enough iterations around p_mosq x n will have been exposed to mosqsuitoes 
mosquito <- rbinom(n, 1, p_mosq)

# 1c. Set up infection probabilities ---- 

# convert probability of boosting in the non-exposed group to odds of boosting in non-exposed group 
odds0 <- p0 / (1 - p0)

# determine the odds of boosting for each individual based on exposure status and OR 
odds <- odds0 * ifelse(mosquito == 1, OR, 1)

# determine the probability of infection from the odds 
p_infect <- odds / (1 + odds)

# random binomial trials to determine who is infected 
infected <- rbinom(n, 1, p_infect)


# 1d. Simulate MIA titre changes ----

# if someone is infected, they come from the gamma dsitribution 
# if not, they come from the normal distribution 
mia <- ifelse(
  infected == 1,
  rgamma(n, shape = gamma_shape, scale = gamma_scale),
  rnorm(n, mean = noise_mean, sd = noise_sd)
)

# 1e. Create a df for exposures ----
dat <- data.frame(
  mosquito = factor(mosquito, labels = c("No mosquitoes", "Mosquitoes")),
  mia = mia
)

# plot the change in titres by exposure status 
ggplot(dat, aes(x = mia)) +
  geom_histogram(
    bins = 30,
    fill = "grey70",
    color = "black"
  ) +
  facet_wrap(~ mosquito, ncol = 2) +
  labs(
    x = "Change in MFI",
    y = "Count",
    title = "Overall MIA Titre Change Distributions by Mosquito Exposure Status"
  ) +
  theme_minimal()

# plot the change in titres overal 
ggplot(dat, aes(x = mia)) +
  geom_histogram(
    binwidth = 0.25,
    fill = "lightgrey",
    color = "black"
  ) +
  labs(
    x = "Change in MFI",
    y = "Count",
    title = "Overall MIA Titre Change Distributions"
  ) +
  theme_minimal()

------------------------------------------------------------------------------------
# 2. Run Mixture Models ---- 
------------------------------------------------------------------------------------
  
  # 2a. Set up the mixture model ----

transform.lambda <- function(x){1/(1 + exp(-x))}
param = c(a0 = 0, b0=1, a1=2, b1=3, lambda1 = 1)

# mixture function 
fitmixture <- function(param, val) {
  
  lambda1 <- transform.lambda(param[["lambda1"]])
  
  # add low positive integer to ensure no issues at zero 
  a0 <- param[["a0"]]                
  b0 <- abs(param[["b0"]]) + 1e-6    
  a1 <- abs(param[["a1"]]) + 1e-6    
  b1 <- abs(param[["b1"]]) + 1e-6    
  
  mix <- lambda1 * dgamma(val, shape = a1, scale = b1) +
    (1 - lambda1) * dnorm(val, mean = a0, sd = b0)
  
  # prevent log of zero 
  mix <- pmax(mix, 1e-300)
  -sum(log(mix))
  
}

# optimisiation function 
optimise <- function(val){
  val = as.numeric(val)
  opt = optim(par = param, 
              fn = fitmixture, 
              method="L-BFGS-B", 
              val = val,
              lower = c(a0 = -Inf, b0 = 1e-6, a1 = 1e-6, b1 = 1e-6, lambda1 = -10),
              upper = c(a0 =  Inf, b0 =  Inf,  a1 =  Inf,  b1 =  Inf,  lambda1 =  10), 
              hessian=FALSE) 
  return(opt)
}

opt <- optimise(dat$mia)

# back transform the lambda value
opt$par[["lambda1"]] <- 1/(1 + exp(-opt$par[["lambda1"]]))

# Set mixture function 
mixture_denv <- function(par){
  
  gamma = dgamma(change_vals, shape = (par[["a1"]]), scale = (par[["b1"]]))
  normal = dnorm(change_vals, mean = (par[["a0"]]), sd=(par[["b0"]]))
  mixture = (par[["lambda1"]] * gamma) + ((1 - par[["lambda1"]]) * normal)
  
  return(mixture)
}

# Extrapolate the values from the mixture 
change_vals = seq(-5, 20, length.out = 1000)
mixture_denv = mixture_denv(opt$par)

# plot mixture model 
denv_df <- data.frame(x = change_vals, y = mixture_denv)

ggplot(data = denv_df, 
       aes(x = change_vals, y = mixture_denv)) + 
  geom_area(alpha = 0.7, fill = "limegreen") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Mixture model for simulated data") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-2, 6))

# plot the componantes 
comp_builder <- function(param, mixture){
  
  gamma <- param[["lambda1"]] * dgamma(change_vals, shape = param[["a1"]], scale = param[["b1"]])
  normal <- (1 - param[["lambda1"]]) * dnorm(change_vals, mean = param[["a0"]], sd = param[["b0"]])
  
  df <- data.frame(
    x = rep(change_vals, 3),
    y = c(mixture, gamma, normal),
    component = rep(c("Mixture", "Gamma", "Normal"),
                    each = length(change_vals))
  )
  
  return(df)
  
}

denv_comp <- comp_builder(opt$par, mixture_denv) 

# Components in the simulated mixture model 
ggplot() +
  geom_density(data = dat, aes(x = mia), stat = "density", alpha = 0.7, fill = "#EDEDED", linewidth = 1) + 
  #geom_histogram(data = dat, aes(x = mia, y = ..density..), binwidth = 0.25, fill = "lightgrey", color = "black") + 
  geom_line(data = denv_comp, aes(x = x, y = y, color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "mixture model for simulated data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-2, 6))

# 2b. Posterior probability ----

post_prob <- function(data, param, id){
  
  # create a vector populated with NA 
  posterior <- rep(NA_real_, length(data))
  complete <- !is.na(data)
  
  # determine the components of the model 
  gamma <- param[["lambda1"]] * dgamma(data[complete], shape = param[["a1"]], scale = param[["b1"]])
  normal <- (1 - param[["lambda1"]]) * dnorm(data[complete], mean = param[["a0"]], sd = param[["b0"]])
  
  # determine the maginal likelihood 
  marginal_lik <- normal + gamma 
  
  # determine the posterior probability for people without NA 
  posterior[complete] <- gamma / marginal_lik
  
  data.frame(ID = id, boost = data, probability = posterior)
  
}

# add a ID column to the data 
dat <- dat %>% mutate(ID = row_number())

# determine the probability of boosting given the mixture model parameters 
deng_prob <- post_prob(dat$mia, opt$par, dat$ID)

# plot posterior probability of boosting 
ymax <- max(denv_comp$y)

ggplot() +
  geom_histogram(data = dat, aes(x = mia, y = ..density..), 
                 fill = "lightgrey", color = "black", binwidth = 0.25) + 
  geom_line(data = denv_comp, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng_prob, aes(x = boost, y = probability * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "Change in MFI",
    y = "Density",
    title = "Mixture model for simulated data ") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-2, 6))

# join the probabilities to the data 
dat <- dat %>% 
  left_join(deng_prob, by = "ID") %>%
  rename(denv_prob = probability)

table(dat$mosquito, useNA = "ifany")

------------------------------------------------------------------------------------
  # 3. Regression ---- 
------------------------------------------------------------------------------------
  
# discretise the data like with the Kucharski et al. (2018) paper 

# add column to the data 
dat <- dat %>%
mutate(
  MIA_status_denv = if_else(denv_prob >= 0.5, 1L, 0L)
)

# 3a. Models ----

# logistic regression with binary outcome 
model1.1 <- glm(MIA_status_denv ~ mosquito , data = dat, family = "binomial")

# fractional logistic regression with continuous outcome bounded by 0,1
model1.2 <- glm(denv_prob ~ mosquito, family = quasibinomial(link = "logit"), data = dat)

model.list=list(model1.1, model1.2)
names.model=c("mosquito","mosquito")

data.tally = dat[,names.model]

# calculate ORs for the different models 
store.oddsTable = NULL
for (jj in seq_along(model.list)) {
  modelT <- model.list[[jj]]
  
  mf <- model.frame(modelT)
  x  <- mf[[2]]
  
  tab <- table(x, useNA = "ifany")
  counts_str <- paste(names(tab), tab, sep="=", collapse="; ")
  
  estm.D <- exp(coef(modelT))
  conf.D <- exp(confint(modelT))
  pval   <- coef(summary(modelT))[, 4][2]
  
  out_row <- cbind(
    var = names.model[jj],
    N_used = nobs(modelT),
    counts = counts_str,
    OR_CI = paste0(round(estm.D[2], 2), " (", paste(round(conf.D[2,], 2), collapse="-"), ")"),
    p = round(pval, 2)
  )
  
  store.oddsTable <- rbind(store.oddsTable, out_row)
}

store.oddsTable

# 3b. Plot the ORs ----

# extracy the data into a df 
forest_dat <- as.data.frame(store.oddsTable, stringsAsFactors = FALSE) %>%
  mutate(
    model = c("Logistic regression", "Fractional logistic"),
    OR    = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
    lower = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
    upper = as.numeric(str_extract(OR_CI, "(?<=-)[0-9.]+"))
  )

# forest plot 
ggplot(forest_dat, aes(x = OR, y = model)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  theme_minimal() + 
  geom_vline(xintercept = 3, linetype = "dashed", colour = "red")

------------------------------------------------------------------------------------
  # 4. Multi-simulation Set Up  ---- 
------------------------------------------------------------------------------------
  
  # 4a. Set up functions within a simulation ----

# transform lambda 
transform.lambda <- function(x) 1/(1 + exp(-x))

# fit mixture 
fitmixture <- function(param, val) {
  
  lambda1 <- transform.lambda(param[["lambda1"]])
  
  # add low positive integer to ensure no issues at zero 
  a0 <- param[["a0"]]                
  b0 <- abs(param[["b0"]]) + 1e-6    
  a1 <- abs(param[["a1"]]) + 1e-6    
  b1 <- abs(param[["b1"]]) + 1e-6    
  
  mix <- lambda1 * dgamma(val, shape = a1, scale = b1) +
    (1 - lambda1) * dnorm(val, mean = a0, sd = b0)
  
  # prevent log of zero 
  mix <- pmax(mix, 1e-300)
  -sum(log(mix))
  
}

# posterior probability 
post_prob_vec <- function(data, param){
  
  # create a vector populated with NA 
  posterior <- rep(NA_real_, length(data))
  complete <- !is.na(data)
  
  # determine the components of the model 
  gamma <- param[["lambda1"]] * dgamma(data[complete], shape = param[["a1"]], scale = param[["b1"]])
  normal <- (1 - param[["lambda1"]]) * dnorm(data[complete], mean = param[["a0"]], sd = param[["b0"]])
  
  # determine the maginal likelihood 
  marginal_lik <- normal + gamma 
  
  # determine the posterior probability for people without NA 
  posterior[complete] <- gamma / marginal_lik
  
  posterior
  
}


# 4b. Run a single simulation ----
one_sim <- function(seed,
                    n = 1000,
                    p_mosq = 0.9,
                    p0 = 0.08,
                    OR_true = 3,
                    noise_mean = 0,
                    noise_sd = 0.3,
                    gamma_shape = 3,
                    gamma_scale = 0.8,
                    threshold = 0.5) {
  
  set.seed(seed)
  
  # mosquito exposure ----
  mosquito <- rbinom(n, 1, p_mosq)
  
  # simulate infection status ----
  odds0 <- p0 / (1 - p0)
  odds  <- odds0 * ifelse(mosquito == 1, OR_true, 1)
  p_inf <- odds / (1 + odds)
  infected <- rbinom(n, 1, p_inf)
  
  # allocate to respective component ----
  mia <- ifelse(
    infected == 1,
    rgamma(n, shape = gamma_shape, scale = gamma_scale),
    rnorm(n, mean = noise_mean, sd = noise_sd)
  )
  
  # build exposure df ----
  dat <- data.frame(
    mosquito = factor(mosquito, labels = c("No mosquitoes", "Mosquitoes")),
    mia = mia
  )
  
  # fit mixture to MIA data ----
  #init <- c(a0 = 0, b0=1, a1=2, b1=3, lambda1 = 1)
  # Better starting values when OR is higher
  init <- c(a0 = 0, b0 = 0.3, a1 = 2.4, b1 = 2.4, lambda1 = 0.5)
  
  opt <- optim(
    par = init,
    fn = fitmixture,
    method = "L-BFGS-B",
    val = dat$mia,
    lower = c(a0 = -Inf, b0 = 1e-6, a1 = 1e-6, b1 = 1e-6, lambda1 = -10),
    upper = c(a0 =  Inf, b0 =  Inf,  a1 =  Inf,  b1 =  Inf,  lambda1 =  10),
    hessian = FALSE
  )
  
  par <- opt$par
  par[["lambda1"]] <- transform.lambda(par[["lambda1"]])
  
  # posterior probability of boosting ----
  p_boost <- post_prob_vec(dat$mia, par)
  
  # logistic regression ----
  boost_bin <- as.integer(p_boost >= threshold)
  model_LR  <- glm(boost_bin ~ mosquito, data = dat, family = binomial())
  OR_LR <- exp(coef(model_LR)[2])
  
  # fractional logistic regression ----
  model_FLR  <- glm(p_boost ~ mosquito, data = dat, family = quasibinomial(link="logit"))
  OR_FLR <- exp(coef(model_FLR)[2])
  
  data.frame(
    seed        = seed,
    OR_true     = OR_true,
    OR_LR       = OR_LR,
    lambda      = par[["lambda1"]],
    p0_est      = plogis(coef(model_LR)[1]),
    noise_mean  = par[["a0"]],
    noise_sd    = par[["b0"]],
    gamma_shape = par[["a1"]],
    gamma_scale = par[["b1"]]
  )
}


# 4c. Run multiple simulations ----

# number of simulations 
B <- 500  

# allow the function to silently overcome NULL values 
safe_one_sim <- function(seed, ...) {
  tryCatch(
    one_sim(seed, ...),
    error = function(e) NULL
  )
}

results <- lapply(1:B, function(s) safe_one_sim(seed = 1000 + s))
out <- do.call(rbind, results[!sapply(results, is.null)])

message("Succeeded: ", nrow(out), " / ", B)


# 4d. Simulation description ----

# tabulate the output parameters from the simulations 

validation_table <- data.frame(
  Parameter   = c(
    "Baseline boost probability",
    "OR (mosquito exposure)",
    "Infection rate (lambda)",
    "Normal: mean",
    "Normal: standard deviation",
    "Gamma: shape",
    "Gamma: scale"
  ),
  True_Value  = c(0.08, 3.0, round(mean(infected), 2), 0.00, 0.30, 3.0, 0.8),
  Mean_Recovered = round(c(
    mean(out$p0_est),
    mean(out$OR_LR),
    mean(out$lambda),
    mean(out$noise_mean),
    mean(out$noise_sd),
    mean(out$gamma_shape),
    mean(out$gamma_scale)
  ), 2),
  SD = round(c(
    sd(out$p0_est),
    sd(out$OR_LR),
    sd(out$lambda),
    sd(out$noise_mean),
    sd(out$noise_sd),
    sd(out$gamma_shape),
    sd(out$gamma_scale)
  ), 2)
)

vt <- validation_table %>%
  flextable() %>%
  set_header_labels(
    Parameter      = "Parameter",
    True_Value     = "True value",
    Mean_Recovered = "Mean recovered",
    SD             = "SD"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

vt

save_as_docx(vt, path = here("outputs", "validation_table.docx"))

# plot the parameters 

# reshape for plotting
plot_dat <- out %>%
  pivot_longer(
    cols = c(OR_LR, lambda, p0_est, noise_mean, noise_sd, gamma_shape, gamma_scale),
    names_to  = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(parameter = recode(parameter,
                            OR_LR       = "OR (mosquito exposure)",
                            lambda      = "Infection rate (lambda)",
                            p0_est      = "Baseline boost probability (p0)",
                            noise_mean  = "Normal mean",
                            noise_sd    = "Normal SD",
                            gamma_shape = "Gamma shape",
                            gamma_scale = "Gamma scale"
  ))

# true values for vertical lines
true_df <- data.frame(
  parameter = c(
    "OR (mosquito exposure)",
    "Infection rate (lambda)",
    "Baseline boost probability (p0)",
    "Normal mean",
    "Normal SD",
    "Gamma shape",
    "Gamma scale"
  ),
  true_value = c(3.0, mean(infected), 0.08, 0.0, 0.3, 3.0, 0.8)
)

ggplot(plot_dat, aes(x = estimate)) +
  geom_histogram(bins = 40, fill = "lightgrey", color = "black") +
  geom_vline(data = true_df, aes(xintercept = true_value),
             linetype = "dashed", color = "red") +
  facet_wrap(~parameter, scales = "free_x") +
  labs(
    x     = "Estimated value",
    y     = "Count",
    title = "Parameter recovery across 500 simulations"
  ) +
  theme_classic() +
  theme(strip.background = element_blank())

