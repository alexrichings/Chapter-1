
# Frequentist finite mixture model (FMM) framework ----  
------------------------------------------------------------------------
  
# 1. Packages ---- 

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(gtsummary)
library(stringr)
library(here)
library(readr)

------------------------------------------------------------------------

# 2. Data ---- 

inputs <- read_csv(here("data", "serology_inputs.csv"))

------------------------------------------------------------------------

# 3. Descriptive Analysis ---- 

# Manage the data ----

# Determine the change in values and scale 
inputs$change_denv1_scaled <- (inputs$DENV1B - inputs$DENV1) / 1000
inputs$change_denv3_scaled <- (inputs$DENV3B - inputs$DENV3) / 1000

# Remove outliers 
inputs$change_denv1_scaled[inputs$change_denv1_scaled <= -5] <- NA
inputs$change_denv3_scaled[inputs$change_denv3_scaled <= -5] <- NA


# Plot the crude change values ----
  
# Density 
# denv3 
ggplot(data = inputs, aes(x = change_denv3_scaled)) + 
geom_density(stat = "density", alpha = 0.7, fill = "grey") + 
labs(x = "∆MFI", 
     y = "Density",
     title = "Change in DENV3 titre (2013 - 2015)") + 
theme_minimal() + 
coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.4))

ggsave(here("outputs", "denv3_density_plot.png"), width = 8, height = 6, dpi = 300)

# denv1 
ggplot(data = inputs, aes(x = change_denv1_scaled)) + 
  geom_density(stat = "density", alpha = 0.7, fill = "grey") + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV1 MIA titre (2013 - 2015)") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.4))

ggsave(here("outputs", "denv1_density_plot.png"), width = 8, height = 6, dpi = 300)


# Count 
# denv3 
ggplot(data = inputs, aes(x = change_denv3_scaled)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(
    x = "∆MFI",
    y = "Count",
    title = "Change in DENV3 MFI titre (2013 - 2015)"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 30), ylim = c(0, 120))

ggsave(here("outputs", "denv3_count_plot.png"), width = 8, height = 6, dpi = 300)

# denv1
ggplot(data = inputs, aes(x = change_denv1_scaled)) +
  geom_histogram(binwidth = 1, fill = "grey", color = "black") +
  labs(
    x = "∆MFI",
    y = "Count",
    title = "Change in DENV1 MFI titre (2013 - 2015)"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 30), ylim = c(0, 120))

ggsave(here("outputs", "denv1_count_plot.png"), width = 8, height = 6, dpi = 300)

-------------------------------------------------------------------------

# 4. FMM ---- 

# Transform lambda 
transform.lambda <- function(x){1/(1 + exp(-x))}

# Set the parameters 

# a0: mean for the normal distribution 
# b0: standard deviation for the normal distribution 
# a1: shape for the gamma distribution > controls peak and skewness
# b1: scale for the gamma distribution > controls spread along the axis

param = c(a0 = 0, b0=1, a1=2, b1=3, lambda1 = 1)

# Set the function 
fitmixture <- function(param, val) {
  
  # double square brackets retrieves just the number
  lambda1 = transform.lambda(param[["lambda1"]]) 
  
  # return the likelihood function 
  return(
    - sum(log(
      lambda1 * dgamma(val, shape = abs(param[["a1"]]), scale = abs(param[["b1"]])) +
        (1-lambda1) * dnorm(val, mean = abs(param[["a0"]]), sd=abs(param[["b0"]]))
    ))
  )
}

# store the boost values 
change_denv1 <- inputs$change_denv1_scaled
change_denv3 <- inputs$change_denv3_scaled

change_denv1 <- change_denv1[!is.na(change_denv1)]
change_denv3 <- change_denv3[!is.na(change_denv3)]

# create a list of the data you want to optimise 
boost_data <- list(
  deng1 = change_denv1,
  deng3 = change_denv3
)

# Optimise 
optimise <- function(val){
  opt = optim(par = param, 
              fn = fitmixture, 
              method="L-BFGS-B", 
              val = val,
              lower=c(rep(0,4),-10), 
              hessian=FALSE) 
  return(opt)
}

# Use lapply 
opt <- lapply(boost_data, optimise)

# lambda needs to be backtransformed 
opt$deng1$par[["lambda1"]] <- 1/(1 + exp(-opt$deng1$par[["lambda1"]]))
opt$deng3$par[["lambda1"]] <- 1/(1 + exp(-opt$deng3$par[["lambda1"]]))


# Simulate data from the parameters 
mixture_denv <- function(par){
  
  gamma = dgamma(change_vals, shape = (par[["a1"]]), scale = (par[["b1"]]))
  normal = dnorm(change_vals, mean = (par[["a0"]]), sd=(par[["b0"]]))
  mixture = (par[["lambda1"]] * gamma) + ((1 - par[["lambda1"]]) * normal)
  
  return(mixture)
}

# Length to simulate over 
change_vals = seq(-5, 20, length.out = 1000)

mixture_denv1 = mixture_denv(opt$deng1$par)
mixture_denv3 = mixture_denv(opt$deng3$par)

# Make a data frame of the x and y values to plot 
denv1_df <- data.frame(x = change_vals, y = mixture_denv1)
denv3_df <- data.frame(x = change_vals, y = mixture_denv3)

# DENV1
ggplot(data = denv1_df, 
       aes(x = change_vals, y = mixture_denv1)) + 
  geom_area(alpha = 0.7, fill = "limegreen") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for DENV1") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.4))

# DENV3
ggplot(data = denv3_df, 
       aes(x = change_vals, y = mixture_denv3)) + 
  geom_area(alpha = 0.7, fill = "limegreen") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "Two component mixture model for DENV3") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.4))

# Build a data frame for the components 
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

denv1_comp <- comp_builder(opt$deng1$par, mixture_denv1)
denv3_comp <- comp_builder(opt$deng3$par, mixture_denv3)

# Plot 
# denv1
ggplot() +
  geom_histogram(data = inputs, aes(x = change_denv1_scaled, y = ..density..), fill = "lightgrey", color = "black", binwidth = 2) + 
  geom_line(data = denv1_comp, aes(x = x, y = y, color = component), linewidth = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) +
  labs(
    x = "∆MFI",
    y = "Density",
    title = "Two component mixture model of DENV1 titre change (2013 - 2015)"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.3))

# denv3
ggplot() +
  geom_histogram(data = inputs, aes(x = change_denv3_scaled, y = ..density..), fill = "lightgrey", color = "black", binwidth = 2) + 
  geom_line(data = denv3_comp, aes(x = x, y = y, color = component), linewidth = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) +
  labs(
    x = "∆MFI",
    y = "Density",
    title = "Two component mixture model of DENV3 titre change (2013 - 2015)"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.3))

----------------------------------------------------------------------------------------
  
# 5. Probability of infection ----

# function to determine the posterior probability of infection 
post_prob <- function(data, param, id){
  
  posterior <- rep(NA_real_, length(data))
  complete <- !is.na(data)
  
  gamma <- param[["lambda1"]] * dgamma(data[complete], shape = param[["a1"]], scale = param[["b1"]])
  normal <- (1 - param[["lambda1"]]) * dnorm(data[complete], mean = param[["a0"]], sd = param[["b0"]])
  marginal_lik <- normal + gamma 
  posterior[complete] <- gamma / marginal_lik
  
  data.frame(PPID2 = id, boost = data, probability = posterior)
}

deng1_prob <- post_prob(inputs$change_denv1_scaled, opt$deng1$par, inputs$PPID2)
deng3_prob <- post_prob(inputs$change_denv3_scaled, opt$deng3$par, inputs$PPID2)

# deng1 
ymax <- 0.3

ggplot() +
  geom_histogram(data = inputs, aes(x = change_denv1_scaled, y = ..density..), fill = "lightgrey", binwidth = 2, color = "black") + 
  geom_line(data = denv1_comp, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng1_prob, aes(x = boost, y = probability * ymax), size = 1, linetype = "dashed", color = "grey40") + 
  labs(
    x = "∆MFI",
    y = "Density") +
  theme_minimal() + 
  scale_y_continuous(
    name = "density",
    sec.axis = sec_axis(~ . / ymax, name = "probability of boosting")) + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.3)) + 
  theme(
    axis.title.y.right = element_text(angle = 90, colour = "grey40"),
    axis.text.y.right  = element_text(colour = "grey40"),
    axis.ticks.y.right = element_line(colour = "grey40"),
    axis.line.y.right  = element_line(colour = "grey40")
  )

ggsave(here("outputs", "denv1_fmm_prob_plot.png"), width = 8, height = 6, dpi = 300)

# deng3
ymax <- 0.3

ggplot() +
  geom_histogram(data = inputs, aes(x = change_denv3_scaled, y = ..density..), fill = "lightgrey", binwidth = 2, color = "black") + 
  geom_line(data = denv3_comp, aes(x = x, y = y , color = component), size = 1) + 
  scale_color_manual(values = c("Mixture" = "black", "Gamma" = "purple", "Normal" = "red")) + 
  geom_line(data = deng3_prob, aes(x = boost, y = probability * ymax), size = 1, linetype = "dashed", color = "grey40") + 
  labs(
    x = "∆MFI",
    y = "Density") +
  theme_minimal() + 
  scale_y_continuous(
    name = "density",
    sec.axis = sec_axis(~ . / ymax, name = "probability of boosting")) + 
  coord_cartesian(xlim = c(-5, 20), ylim = c(0,0.3)) + 
  theme(
    axis.title.y.right = element_text(angle = 90, colour = "grey40"),
    axis.text.y.right  = element_text(colour = "grey40"),
    axis.ticks.y.right = element_line(colour = "grey40"),
    axis.line.y.right  = element_line(colour = "grey40")
  )

ggsave(here("outputs", "denv3_fmm_prob_plot.png"), width = 8, height = 6, dpi = 300)

----------------------------------------------------------------------------------------

# 5. Join the data ----

inputs <- inputs %>% 
  left_join(deng1_prob, by = "PPID2") %>%
  rename(denv1_prob = probability)

inputs <- inputs %>% 
  left_join(deng3_prob, by = "PPID2") %>%
  rename(denv3_prob = probability)

inputs_fj <- inputs %>% select(-boost.x, -boost.y)

saveRDS(inputs_fj, here("data", "inputs_fj.rds"))

