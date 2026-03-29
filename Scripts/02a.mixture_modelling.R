# mixture modelling for DR data 

library(ggplot2)
----------------------------------------------------------------------------------------------
  # 1. Load the data ---- 
df_wide <- readRDS("Data/DR_surveillance/EN/df_wide")
df_wide <- readRDS("Data/DR_surveillance/Serology/df_wide")

----------------------------------------------------------------------------------------------
  # 2. Plot the data ----

# All DENV serotypes 
ggplot(data = df_wide) + 
  geom_density(aes(x = change_dengns1_1/1000, fill = "DENV1"), alpha = 0.5) + 
  geom_density(aes(x = change_dengns1_2/1000, fill = "DENV2"), alpha = 0.5) + 
  geom_density(aes(x = change_dengns1_3/1000, fill = "DENV3"), alpha = 0.5) + 
  geom_density(aes(x = change_dengns1_4/1000, fill = "DENV4"), alpha = 0.5) + 
  geom_density(aes(x = change_zika_ns1/1000, fill = "ZIKV"), alpha = 0.5) + 
  geom_density(aes(x = change_chik_e1/1000, fill = "CHIKV"), alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in MFI values in DR",
       fill = "Serotype") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-15, 15))

# deng1
ggplot(data = df_wide) + 
  geom_density(aes(x = change_dengns1_1/1000), fill = "red", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV1 MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, 0.03))

# deng2
ggplot(data = df_wide) + 
  geom_density(aes(x = change_dengns1_2/1000), fill = "yellow", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV2 MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, 0.03))

# deng3
ggplot(data = df_wide) + 
  geom_density(aes(x = change_dengns1_3/1000), fill = "green", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV3 MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, 0.03))

# deng4 
ggplot(data = df_wide) + 
  geom_density(aes(x = change_dengns1_4/1000), fill = "blue", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV4 MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60))

# CHIKV
ggplot(data = df_wide) + 
  geom_density(aes(x = change_chik_e1/1000), fill = "purple", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in CHIKV MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60))

# ZIKV
ggplot(data = df_wide) + 
  geom_density(aes(x = change_zika_ns1/1000), fill = "violet", alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in ZIKV MFI values in DR") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, 0.03))

----------------------------------------------------------------------------------------------
  # 3. Prep the data ----

# create a new df to edit
inputs <- df_wide

# scale the data 
inputs$change_dengns1_1_scale <- inputs$change_dengns1_1/1000
inputs$change_dengns1_2_scale <- inputs$change_dengns1_2/1000
inputs$change_dengns1_3_scale <- inputs$change_dengns1_3/1000
inputs$change_dengns1_4_scale <- inputs$change_dengns1_4/1000
inputs$change_zika_ns1_scale <- inputs$change_zika_ns1/1000
inputs$change_chik_e1_scale <- inputs$change_chik_e1/1000

inputs$dengns1_1_scale <- inputs$dengns1_1_mfi_sero/1000
inputs$dengns1_2_scale <- inputs$dengns1_1_mfi_sero/1000
inputs$dengns1_3_scale <- inputs$dengns1_1_mfi_sero/1000
inputs$dengns1_4_scale <- inputs$dengns1_1_mfi_sero/1000
inputs$zika_ns1_scale <- inputs$zika_ns1_mfi_sero/1000
inputs$chik_e1_scale <- inputs$chik_e1_mfi_sero/1000

sum(is.na(inputs$change_dengns1_1_scale))

--------------------------------------------------------------------------
  # 4. Constrain mixing proportions ----

# Mixing proportions should range from 0 - 1
transform.lambda <- function(x){1/(1 + exp(-x))}

--------------------------------------------------------------------------
  # 5. Model set up ---- 

# Set the parameters 

# a0: mean for the normal distribution 
# b0: standard deviation for the normal distribution 
# a1: shape for the gamma distribution > controls peak and skewness
# b1: scale for the gamma distribution > controls spread along the axis

# Set the function 
fitmixture <- function(param, val, a0_fixed = NULL) {
  
  lambda1 <- transform.lambda(param[["lambda1"]])
  
  # non-boost components: use fixed a0 if provided, otherwise estimate it
  a0 <- if (is.null(a0_fixed)) param[["a0"]] else a0_fixed
  b0 <- pmax(param[["b0"]], 1e-5)
  
  a1 <- pmax(param[["a1"]], 1e-5)
  b1 <- pmax(param[["b1"]], 1e-5)
  
  non_boost <- dnorm(val, mean = a0, sd = b0)
  #boost     <- dgamma(pmax(val, 1e-10), shape = a1, scale = b1)
  boost <- dgamma(val, shape = a1, scale = b1)
  
  
  epsilon <- 1e-10
  -sum(log((1 - lambda1) * non_boost + lambda1 * boost + epsilon))
}

----
fitmixture <- function(param, val, a0_fixed = NULL) {
  
  lambda1 <- transform.lambda(param[["lambda1"]])
  a0 <- if (is.null(a0_fixed)) param[["a0"]] else a0_fixed
  b0 <- pmax(param[["b0"]], 1e-5)
  a1 <- pmax(param[["a1"]], 1e-5)
  b1 <- pmax(param[["b1"]], 1e-5)
  
  non_boost <- dnorm(val, mean = a0, sd = b0)
  boost     <- dgamma(pmax(val, 1e-10), shape = a1, scale = b1)
  
  epsilon <- 1e-10
  
  # for positive values: full mixture
  # for negative values: only normal component
  lik <- ifelse(val > 0,
                (1 - lambda1) * non_boost + lambda1 * boost,
                (1 - lambda1) * non_boost)
  
  -sum(log(lik + epsilon))
}
----

# function to fit the mixture model to the data from any pathogen 
fit_pathogen <- function(
    data,
    start_par,
    fix_a0 = TRUE,     
    lower = NULL
) {
  data <- data[!is.na(data)]
  
  if (fix_a0) {
    # a0 is fixed at 0
    start_par <- start_par[names(start_par) != "a0"]
    if (is.null(lower)) lower <- c(b0 = 0, a1 = 0, b1 = 0, lambda1 = -10)
    
    opt <- optim(
      par     = start_par,
      fn      = fitmixture,
      method  = "L-BFGS-B",
      val     = data,
      lower   = lower,
      a0_fixed = 0,          
      hessian = FALSE
    )
    
  } else {
    # a0 is estimated
    if (!("a0" %in% names(start_par))) start_par <- c(a0 = 0, start_par)
    if (is.null(lower)) lower <- c(a0 = -Inf, b0 = 0, a1 = 0, b1 = 0, lambda1 = -10)
    
    opt <- optim(
      par     = start_par,
      fn      = fitmixture,
      method  = "L-BFGS-B",
      val     = data,
      lower   = lower,
      a0_fixed = NULL,      
      hessian = FALSE
    )
  }
  
  opt
}


# provide the starting parameters 
start_par <- c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)

# deng data

# deng1
deng1_fixed <- fit_pathogen(inputs$change_dengns1_1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
deng1_free <- fit_pathogen(inputs$change_dengns1_1_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)
deng1_start_fixed <- fit_pathogen(inputs$dengns1_1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)

# deng2
deng2_fixed <- fit_pathogen(inputs$change_dengns1_2_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
deng2_free <- fit_pathogen(inputs$change_dengns1_2_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)

# deng3
deng3_fixed <- fit_pathogen(inputs$change_dengns1_3_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
deng3_free <- fit_pathogen(inputs$change_dengns1_3_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)

# deng4
deng4_fixed <- fit_pathogen(inputs$change_dengns1_4_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
deng4_free <- fit_pathogen(inputs$change_dengns1_4_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)

# zika
zika_fixed <- fit_pathogen(inputs$change_zika_ns1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
zika_free <- fit_pathogen(inputs$change_zika_ns1_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)

# chikungunya
chik_fixed <- fit_pathogen(inputs$change_chik_e1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = TRUE)
chik_free <- fit_pathogen(inputs$change_chik_e1_scale, c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1), fix_a0 = FALSE)

# retrieve parameters 

deng1_fixed_par <- deng1_fixed$par
deng1_free_par <- deng1_free$par
deng1_start_par <- deng1_start_fixed$par

deng2_fixed_par <- deng2_fixed$par
deng2_free_par <- deng2_free$par

deng3_fixed_par <- deng3_fixed$par
deng3_free_par <- deng3_free$par

deng4_fixed_par <- deng4_fixed$par
deng4_free_par <- deng4_free$par

zika_fixed_par <- zika_fixed$par
zika_free_par <- zika_free$par

chik_fixed_par <- chik_fixed$par
chik_free_par <- chik_free$par

--------------------------------------------------------------------------
  # 6. Model wrangling ---- 

# function to retrieve parameters ----
extract_params <- function(opt_par, fix_a0 = TRUE) {
  list(
    lambda1 = transform.lambda(opt_par["lambda1"]),
    a0 = if (fix_a0) 0 else opt_par["a0"],
    b0 = opt_par["b0"],
    a1 = opt_par["a1"],
    b1 = opt_par["b1"]
  )
}


# retrieve parameters for all pathogens ----

# deng1 
deng1_free_par <- extract_params(deng1_free_par, fix_a0 = FALSE)
deng1_fixed_par <- extract_params(deng1_fixed_par, fix_a0 = TRUE)
deng1_start_par <- extract_params(deng1_start_par, fix_a0 = TRUE)

# deng2
deng2_free_par <- extract_params(deng2_free_par, fix_a0 = FALSE)
deng2_fixed_par <- extract_params(deng2_fixed_par, fix_a0 = TRUE)

# deng3
deng3_free_par <- extract_params(deng3_free_par, fix_a0 = FALSE)
deng3_fixed_par <- extract_params(deng3_fixed_par, fix_a0 = TRUE)

# deng4
deng4_free_par <- extract_params(deng4_free_par, fix_a0 = FALSE)
deng4_fixed_par <- extract_params(deng4_fixed_par, fix_a0 = TRUE)

# zika
zika_free_par <- extract_params(zika_free_par, fix_a0 = FALSE)
zika_fixed_par <- extract_params(zika_fixed_par, fix_a0 = TRUE)

# chik
chik_free_par <- extract_params(chik_free_par, fix_a0 = FALSE)
chik_fixed_par <- extract_params(chik_fixed_par, fix_a0 = TRUE)


# function to model the components ----
mixture_density <- function(x, pars) {
  
  normal <- dnorm(x, mean = pars$a0, sd = pars$b0)
  gamma  <- dgamma(x, shape = pars$a1, scale = pars$b1)
  
  ((1 - pars$lambda1) * normal) + (pars$lambda1 * gamma)
}

# simulate data based off the optimised parameters ----
simulate_mixture <- function(data, pars, n = 1000) {
  
  data <- data[is.finite(data)] # edit 
  x_vals <- seq(min(data), max(data), length.out = n)
  y_vals <- mixture_density(x_vals, pars)
  data.frame(x = x_vals, y = y_vals)
}

# simulate data for each pathogen ----

sim_deng1_fixed_df <- simulate_mixture(inputs$change_dengns1_1_scale, deng1_fixed_par)
sim_deng1_free_df <- simulate_mixture(inputs$change_dengns1_1_scale, deng1_free_par)
sim_deng1_start_df <- simulate_mixture(inputs$dengns1_1_scale, deng1_start_par)

sim_deng2_fixed_df <- simulate_mixture(inputs$change_dengns1_2_scale, deng2_fixed_par)
sim_deng2_free_df <- simulate_mixture(inputs$change_dengns1_2_scale, deng2_free_par)

sim_deng3_fixed_df <- simulate_mixture(inputs$change_dengns1_3_scale, deng3_fixed_par)
sim_deng3_free_df <- simulate_mixture(inputs$change_dengns1_3_scale, deng3_free_par)

sim_deng4_fixed_df <- simulate_mixture(inputs$change_dengns1_4_scale, deng4_fixed_par)
sim_deng4_free_df <- simulate_mixture(inputs$change_dengns1_4_scale, deng4_free_par)

sim_zika_fixed_df <- simulate_mixture(inputs$change_zika_ns1_scale, zika_fixed_par)
sim_zika_free_df <- simulate_mixture(inputs$change_zika_ns1_scale, zika_free_par)

sim_chik_fixed_df <- simulate_mixture(inputs$change_chik_e1_scale, chik_fixed_par)
sim_chik_free_df <- simulate_mixture(inputs$change_chik_e1_scale, chik_free_par)


# Use ggplot to plot the distribution of the mixture model ----

# deng1 ----
ggplot(data = sim_deng1_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV1 data") + 
  theme_minimal()

ggplot(data = sim_deng1_free_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV1 data") + 
  theme_minimal()

ggplot(data = sim_deng1_start_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV1 data") + 
  theme_minimal()

# deng2 ----
ggplot(data = sim_deng2_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV2 data") + 
  theme_minimal()

ggplot(data = sim_deng2_free_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV2 data") + 
  theme_minimal()

# deng3 ----
ggplot(data = sim_deng3_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV3 data") + 
  theme_minimal()

# deng4 ----
ggplot(data = sim_deng4_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆DENV4 MIA", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV4 data") + 
  theme_minimal()

# chik ----
ggplot(data = sim_chik_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI CHIKV data") + 
  theme_minimal()

# zika ----
ggplot(data = sim_zika_fixed_df, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI Zika data") + 
  theme_minimal()



# function to create a df for individual component densities ---- 
comp_density <- function(pars, simulated_data){
  
  normal <- (1 - pars$lambda1) * dnorm(simulated_data$x, mean = pars$a0, sd = pars$b0)
  gamma <- pars$lambda1 * dgamma(pmax(simulated_data$x, 1e-10), shape = pars$a1, scale = pars$b1)
  mixture <- normal + gamma 
  
  x = rep(simulated_data$x, 3)
  y = c(mixture, normal, gamma)
  component = rep(c("Mixture", "Normal", "Gamma"), each = length(simulated_data$x))
  
  data.frame(x = x, y = y, component = component)
}

# create the component df for each pathogen ----
deng1_comp_df <- comp_density(deng1_fixed_par, sim_deng1_fixed_df)
deng1_start_comp_df <- comp_density(deng1_start_par, sim_deng1_start_df)

deng2_comp_df <- comp_density(deng2_fixed_par, sim_deng2_fixed_df)
deng3_comp_df <- comp_density(deng3_fixed_par, sim_deng3_fixed_df)
deng4_comp_df <- comp_density(deng4_fixed_par, sim_deng4_fixed_df)
zika_comp_df <- comp_density(zika_fixed_par, sim_zika_fixed_df)

# plot the components of each mixture model ----

# deng1 ----
ggplot(deng1_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng1_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "Two-component mixture model for MIA DENV1 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

ggplot(deng1_start_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng1_start_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "MIA",
    y = "Density",
    title = "Two-component mixture model for MIA DENV1 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 0.1))

# deng2 ----
ggplot(deng2_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng2_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "Two-component mixture model for MIA DENV2 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# deng3 ----
ggplot(deng3_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng3_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "Two-component mixture model for MIA DENV3 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# deng4 ----
ggplot(deng4_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng4_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MFI",
    y = "Density",
    title = "Two-component mixture model for MIA DENV4 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# zika ----
ggplot(zika_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(zika_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "Two-component mixture model for MIA Zika data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

--------------------------------------------------------------------------
  # 7. Determining posterior probability of infection ---- 

# posterior probability = (likelihood x prior) / marginal likelihood 
# probability of boosting given a particular increase in MFI 

# function to determine the posterior probability of infection ----
post_prob <- function(data, pars, id){
  
  x <- data 
  normal <- (1 - pars$lambda1) * dnorm(data, mean = pars$a0, sd = pars$b0)
  gamma <- pars$lambda1 * dgamma(pmax(data, 1e-10), shape = pars$a1, scale = pars$b1)
  marginal_lik <- normal + gamma 
  posterior <- gamma / marginal_lik
  
  data.frame(id = id, x = x, y = posterior)
}

# function to determine the posterior probability of naivety----
post_naive <- function(data, pars, id){
  
  x <- data 
  naive <- (1 - pars$lambda1) * dnorm(data, mean = pars$a0, sd = pars$b0)
  expo <- pars$lambda1 * dnorm(data, mean = pars$a1, sd = pars$b1)
  marginal_lik <- expo + naive
  posterior <- naive / marginal_lik
  
  data.frame(id = id, x = x, y = posterior)
}

# determine the posterior probability of infection for each pathogen ----
deng1_inf_prob <- post_prob(inputs$change_dengns1_1_scale, deng1_fixed_par, inputs$cohort_id)
deng1_naive_prob <- post_naive(inputs$dengns1_1_scale, deng1_start_par, inputs$cohort_id)

deng2_inf_prob <- post_prob(inputs$change_dengns1_2_scale, deng2_fixed_par, inputs$cohort_id)
deng3_inf_prob <- post_prob(inputs$change_dengns1_3_scale, deng3_fixed_par, inputs$cohort_id)
deng4_inf_prob <- post_prob(inputs$change_dengns1_4_scale, deng4_fixed_par, inputs$cohort_id)
zika_inf_prob <- post_prob(inputs$change_zika_ns1_scale, zika_fixed_par, inputs$cohort_id)

# maximum probability of infection given a ∆MFI value ----
summary(deng1_inf_prob$y) # 0.04
summary(deng2_inf_prob$y) # 0.527
summary(deng3_inf_prob$y) # 0.544
summary(deng4_inf_prob$y) # 0.506
summary(zika_inf_prob$y) # 0.857

# scaling factor 
ymax <- 0.125

# deng1 ----
ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_1_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng1_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng1_inf_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "∆MIA",
    y = "Density",
    title = "DENV1") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-30, 50), ylim = c(0, ymax)) +
  theme(panel.grid = element_blank())

ggplot() +
  geom_histogram(data = inputs, aes(x = dengns1_1_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 2) + 
  geom_line(data = deng1_start_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng1_naive_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "MIA",
    y = "Density",
    title = "DENV1") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(0, 70), ylim = c(0, ymax))




# deng2 ----
ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_2_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng2_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng2_inf_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "∆MIA",
    y = "Density",
    title = "DENV2") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, ymax))

# deng3 ----
ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_3_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng3_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng3_inf_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "∆MIA",
    y = "Density",
    title = "DENV3") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, ymax))

# deng4 ----
ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_4_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng4_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng4_inf_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "∆MIA",
    y = "Density",
    title = "DENV4") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, ymax))

# zika ----
ggplot() +
  geom_histogram(data = inputs, aes(x = change_zika_ns1_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = zika_comp_df, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = zika_inf_prob, aes(x = x, y = y * ymax), size = 1, linetype = "dashed") + 
  labs(
    x = "∆MIA",
    y = "Density",
    title = "ZIKV") +
  theme_minimal() + 
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~ . / ymax, name = "Posterior Probability")) + 
  coord_cartesian(xlim = c(-60, 60), ylim = c(0, ymax))

--------------------------------------------------------------------------
  # 8. Impact of initial MFI on probability of infection  ---- 

# p(infection) = 1 - p(no infection) 
# p(infection) = 1 - p(no infection(a) * no infection(b) ... no(infection(n)))

# create a df with infection data ----

df_wide <- df_wide %>% 
  left_join(deng1_inf_prob %>% select(id, y) %>% rename(prob_deng1 = y), by = c("cohort_id" = "id")) %>%
  left_join(deng2_inf_prob %>% select(id, y) %>% rename(prob_deng2 = y), by = c("cohort_id" = "id")) %>%
  left_join(deng3_inf_prob %>% select(id, y) %>% rename(prob_deng3 = y), by = c("cohort_id" = "id")) %>%
  left_join(deng4_inf_prob %>% select(id, y) %>% rename(prob_deng4 = y), by = c("cohort_id" = "id")) %>%
  left_join(zika_inf_prob %>% select(id, y) %>% rename(prob_zika = y), by = c("cohort_id" = "id")) 


# MFI vs p(infection) curve ----

# deng1 ----
ggplot(df_wide, aes(dengns1_1_mfi_sero/1000, prob_deng1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "deng1"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,60))

# deng2 ----
ggplot(df_wide, aes(dengns1_2_mfi_sero/1000, prob_deng2)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "deng2"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,60))

# deng3 ----
ggplot(df_wide, aes(dengns1_3_mfi_sero/1000, prob_deng3)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "deng3"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,60))

# deng4 ----
ggplot(df_wide, aes(dengns1_4_mfi_sero/1000, prob_deng4)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "deng4"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,60))

# zika----
ggplot(df_wide, aes(zika_ns1_mfi_sero/1000, prob_zika)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "zika"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,60))

---------------------------------------------------------------------------------
  # 11. Account for low probability bias ---- 

subset_prob <- function(data, col, level){
  prob <- data %>% filter({{col}} > level)
  return(prob)
}

deng1_01 <- subset_prob(df_wide, prob_deng1, 0.01)
deng2_01 <- subset_prob(df_wide, prob_deng2, 0.01)
deng3_01 <- subset_prob(df_wide, prob_deng3, 0.01)
deng4_01 <- subset_prob(df_wide, prob_deng4, 0.01)
zika_01 <- subset_prob(df_wide, prob_zika, 0.01)

# deng1 

ggplot(deng1_01, aes(x = dengns1_1_mfi_sero/1000, y = prob_deng1)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = TRUE, color = "orange") +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "DENV1: (p > 0.01)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  xlim(0, 60)

# deng2
ggplot(deng2_01, aes(x = dengns1_2_mfi_sero/1000, y = prob_deng2)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = TRUE, color = "orange") +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "DENV2: (p > 0.01)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  xlim(0, 60)

# deng3
ggplot(deng3_01, aes(x = dengns1_3_mfi_sero/1000, y = prob_deng3)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = TRUE, color = "orange") +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "DENV3: (p > 0.01)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  xlim(0, 60)

# deng4
ggplot(deng4_01, aes(x = dengns1_4_mfi_sero/1000, y = prob_deng4)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = TRUE, color = "orange") +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "DENV4: (p > 0.01)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  xlim(0, 60)

# zika 
ggplot(zika_01, aes(x = zika_ns1_mfi_sero/1000, y = prob_zika)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "loess", se = TRUE, color = "orange") +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "ZIKV: (p > 0.01)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  xlim(0, 60)

---------------------------------------------------------------------------------
  # 11. Save the data ----

saveRDS(df_wide, "Data/DR_surveillance/EN/df_wide")


------
  
naive <- inputs %>% filter(inputs$dengns1_1_scale <= 10)
primed <- inputs %>% filter(inputs$dengns1_1_scale > 10)

ggplot(naive, aes(x = age_sero)) + 
  geom_histogram(fill = "lightgrey", color = "black", binwidth = 2) + 
  theme_minimal()

ggplot(primed, aes(x = age_sero)) + 
  geom_histogram(fill = "lightgrey", color = "black", binwidth = 2) + 
  theme_minimal()

summary(naive$age_sero)
summary(primed$age_sero)

ggplot(inputs, aes(x = age_sero)) +
  geom_histogram(color = "black", fill = "lightgrey", binwidth = 5) + 
  theme_minimal()


-----------
  
# e. fit mixture model for the basic initial MFI data ----

# Scalar version for optim
fitmixture_init <- function(param, val) {
  
  lambda1 <- transform.lambda(param[["lambda1"]])
  a0 <- param[["a0"]]
  b0 <- param[["b0"]]
  a1 <- param[["a1"]]
  b1 <- param[["b1"]]
  
  naive <- dnorm(val, mean = a0, sd = b0)
  exposed <- dnorm(val, mean = a1, sd = b1)
  
  lik <- (1 - lambda1) * naive + lambda1 * exposed
  
  -sum(log(lik))
}

# Run optim
opt <- optim(
  par     = c(a0 = 0, b0 = 1, a1 = 20, b1 = 5, lambda1 = 1),
  fn      = fitmixture_init,
  method  = "L-BFGS-B",
  lower  = c(a0 = -Inf, b0 = 0.01, a1 = -Inf, b1 = 0.01, lambda1 = -Inf),
  val     = inputs$dengns1_1_scale,
  hessian = FALSE
)

extract_param <- function(opt_par) {
  list(
    lambda1 = transform.lambda(opt_par["lambda1"]),
    a0 = 0,
    b0 = opt_par["b0"],
    a1 = opt_par["a1"],
    b1 = opt_par["b1"]
  )
}

deng1_init_param <- extract_param(opt$par)
deng1_init_param

mixture_den <- function(x, pars) {
  
  naive <- dnorm(x, mean = pars$a0, sd = pars$b0)
  exposed  <- dnorm(x, mean = pars$a1, sd = pars$b1)
  
  ((1 - pars$lambda1) * naive) + (pars$lambda1 * exposed)
}

# simulate data based off the optimised parameters ----
simulate_mix <- function(data, pars, n = 1000) {
  
  data <- data[is.finite(data)] # edit 
  x_vals <- seq(min(data), max(data), length.out = n)
  y_vals <- mixture_den(x_vals, pars)
  data.frame(x = x_vals, y = y_vals)
}

# simulate data for each pathogen ----

sim_deng1_init <- simulate_mix(inputs$dengns1_1_scale, deng1_init_param)

ggplot(data = sim_deng1_init, 
       aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "MFI", 
       y = "Density", 
       title = "Two component mixture model for MFI DENV1 data") + 
  theme_minimal()

comp_den <- function(pars, simulated_data){
  
  naive <- (1 - pars$lambda1) * dnorm(simulated_data$x, mean = pars$a0, sd = pars$b0)
  exposed <- pars$lambda1 * dnorm(simulated_data$x, mean = pars$a1, sd = pars$b1)
  mixture <- naive + exposed
  
  x = rep(simulated_data$x, 3)
  y = c(mixture, naive, exposed)
  component = rep(c("Mixture", "Naive", "Exposed"), each = length(simulated_data$x))
  
  data.frame(x = x, y = y, component = component)
}

# create the component df for each pathogen ----
deng1_initcomp_df <- comp_den(deng1_init_param, sim_deng1_init)

# deng1 ----
ggplot(deng1_initcomp_df, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng1_initcomp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Naive" = "red", "Exposed" = "purple")) +
  labs(
    x = "MFI",
    y = "Density",
    title = "Two-component mixture model for MIA DENV1 data"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))




