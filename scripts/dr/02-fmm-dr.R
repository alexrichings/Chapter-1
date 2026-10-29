
# Chapter 1 

-------------------------------------------------------------------------------
  
# Script: 02-fmm-dr ----  

# Purpose: 
# Run two-component finite mixture model (FMM) in a frequentist framework

# Output: 
# Individual probability of boosting for each person per dengue serotype (inputs_dr) and
# fmm parameters 

# Plots: 
# deng1_fmm_prob_plot.png
# deng2_fmm_prob_plot.png
# deng3_fmm_prob_plot.png
# deng4_fmm_prob_plot.png
# facet_fmm_plot.png

-------------------------------------------------------------------------------
  
# Load packages and data ----

# Packages 
library("ggplot2")
library("dplyr")
library("here")

# Paired data 
df_wide <- readRDS(here("data", "df_wide.rds"))

# Population ≤ 20 years 
young <- readRDS(here("data", "young.rds"))

-------------------------------------------------------------------------------
# Section 1 ----
------------------------------------------------------------------------------- 
  
# 1.1. Edits to data ----

# Create data frame to edit 
inputs <- df_wide
names(inputs)

# Scale data to equalise scale and facilitate convergence 
inputs$change_dengns1_1_scale <- inputs$change_dengns1_1/1000
inputs$change_dengns1_2_scale <- inputs$change_dengns1_2/1000
inputs$change_dengns1_3_scale <- inputs$change_dengns1_3/1000
inputs$change_dengns1_4_scale <- inputs$change_dengns1_4/1000

# Scale data to equalise scale and facilitate convergence 
young$change_dengns1_1_scale <- young$change_dengns1_1/1000
young$change_dengns1_2_scale <- young$change_dengns1_2/1000
young$change_dengns1_3_scale <- young$change_dengns1_3/1000
young$change_dengns1_4_scale <- young$change_dengns1_4/1000


-------------------------------------------------------------------------------
# Section 2 ----
------------------------------------------------------------------------------- 
  
# 2.1: Model set up ----
  
# Constrain mixing proportions within [0,1]
transform.lambda <- function(x){1/(1 + exp(-x))}

# fmm function 
fitmixture <- function(param, val) {
  lambda1 <- transform.lambda(param[["lambda1"]])
  
  # non-boost 
  a0 <- 0 # mean 
  b0 <- pmax(param[["b0"]], 1e-5) # sd 
  
  # boost 
  a1 <- pmax(param[["a1"]], 1e-5) # shape (skew)
  b1 <- pmax(param[["b1"]], 1e-5) # scale (spread)
  
  # components 
  non_boost <- dnorm(val, mean = a0, sd = b0)
  boost     <- dgamma(pmax(val, 1e-10), shape = a1, scale = b1)
  
  # for x > 0: full mixture
  # for x ≤ 0: only normal component
  lik <- ifelse(val > 0, ((1 - lambda1) * non_boost) + (lambda1 * boost), (1 - lambda1) * non_boost)
  
  # likelihood 
  epsilon <- 1e-10
  -sum(log(lik + epsilon))
}

# function to fit to different pathogens 
fit_pathogen <- function(data, start_par, lower = NULL) {
  
  data <- data[!is.na(data)]
  
  # remove a0 if it was accidentally included
  start_par <- start_par[names(start_par) != "a0"]
  
  # default lower bounds
  if (is.null(lower)) {
    lower <- c(
      b0 = 1e-5,
      a1 = 1e-5,
      b1 = 1e-5,
      lambda1 = -10
    )
  }
  
  opt <- optim(
    par     = start_par,
    fn      = fitmixture,
    method  = "L-BFGS-B",
    val     = data,
    lower   = lower,
    hessian = FALSE
  )
  
  opt
}

-------------------------------------------------------------------------------
  
# 2.2. Model fitting  ----

# Starting parameters 
start_par <- c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)

# Fit fmm

# a) whole population 
deng1_fmm <- fit_pathogen(inputs$change_dengns1_1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng2_fmm <- fit_pathogen(inputs$change_dengns1_2_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng3_fmm <- fit_pathogen(inputs$change_dengns1_3_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng4_fmm <- fit_pathogen(inputs$change_dengns1_4_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))

deng1_par <- deng1_fmm$par
deng2_par <- deng2_fmm$par
deng3_par <- deng3_fmm$par
deng4_par <- deng4_fmm$par

# b) population ≤ 20 years 
youth_deng1_fmm <- fit_pathogen(young$change_dengns1_1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
youth_deng2_fmm <- fit_pathogen(young$change_dengns1_2_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
youth_deng3_fmm <- fit_pathogen(young$change_dengns1_3_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
youth_deng4_fmm <- fit_pathogen(young$change_dengns1_4_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))

youth_deng1_par <- youth_deng1_fmm$par
youth_deng2_par <- youth_deng2_fmm$par
youth_deng3_par <- youth_deng3_fmm$par
youth_deng4_par <- youth_deng4_fmm$par


# retrieve parameters 
extract_params <- function(opt_par) {
  list(
    lambda1 = transform.lambda(opt_par["lambda1"]),
    a0 = 0,
    b0 = opt_par["b0"],
    a1 = opt_par["a1"],
    b1 = opt_par["b1"]
  )
}

# a) whole population 
deng1_par <- extract_params(deng1_par)
deng2_par <- extract_params(deng2_par)
deng3_par <- extract_params(deng3_par)
deng4_par <- extract_params(deng4_par)

# b) population ≤ 20 years 
youth_deng1_par <- extract_params(youth_deng1_par)
youth_deng2_par <- extract_params(youth_deng2_par)
youth_deng3_par <- extract_params(youth_deng3_par)
youth_deng4_par <- extract_params(youth_deng4_par)


# Specify model components 
mixture_density <- function(x, pars) {
  
  normal <- dnorm(x, mean = pars$a0, sd = pars$b0)
  gamma  <- dgamma(x, shape = pars$a1, scale = pars$b1)
  
  ((1 - pars$lambda1) * normal) + (pars$lambda1 * gamma)
}

# Simulate data based off the optimised parameters 
simulate_mixture <- function(data, pars, n = 1000) {
  
  data <- data[is.finite(data)] 
  
  x_vals <- seq(min(data), max(data), length.out = n)
  y_vals <- mixture_density(x_vals, pars)
  
  data.frame(x = x_vals, y = y_vals)
}

# a) whole population 
sim_deng1 <- simulate_mixture(inputs$change_dengns1_1_scale, deng1_par)
sim_deng2 <- simulate_mixture(inputs$change_dengns1_2_scale, deng2_par)
sim_deng3 <- simulate_mixture(inputs$change_dengns1_3_scale, deng3_par)
sim_deng4 <- simulate_mixture(inputs$change_dengns1_4_scale, deng4_par)

# b) population ≤ 20 years 
youth_sim_deng1 <- simulate_mixture(young$change_dengns1_1_scale, youth_deng1_par)
youth_sim_deng2 <- simulate_mixture(young$change_dengns1_2_scale, youth_deng2_par)
youth_sim_deng3 <- simulate_mixture(young$change_dengns1_3_scale, youth_deng3_par)
youth_sim_deng4 <- simulate_mixture(young$change_dengns1_4_scale, youth_deng4_par)


-------------------------------------------------------------------------------

# 2.3: FMM visualisation ----  

# Plot components 

# Function to plot each component 
plot_mfi <- function(data, figure_label, serotype, color) {
  ggplot(data, aes(x = x, y = y)) + 
    geom_area(alpha = 0.7, fill = color) + 
    labs(x = figure_label, y = "Density", title = serotype) + 
    theme_minimal()
}

# a) Whole population 
plot_mfi(sim_deng1, "∆MFI", "DENV1 FMM", "lightgrey")
plot_mfi(sim_deng2, "∆MFI", "DENV2 FMM", "lightgrey")
plot_mfi(sim_deng3, "∆MFI", "DENV3 FMM", "lightgrey")
plot_mfi(sim_deng4, "∆MFI", "DENV4 FMM", "lightgrey")

# b) Population ≤20 years 
plot_mfi(youth_sim_deng1, "∆MFI", "DENV1 FMM (≤ 20 years)", "lightgrey")
plot_mfi(youth_sim_deng2, "∆MFI", "DENV2 FMM (≤ 20 years)", "lightgrey")
plot_mfi(youth_sim_deng3, "∆MFI", "DENV3 FMM (≤ 20 years)", "lightgrey")
plot_mfi(youth_sim_deng4, "∆MFI", "DENV4 FMM (≤ 20 years)", "lightgrey")


# Function to create a df for individual component densities
comp_density <- function(pars, simulated_data){
  
  normal <- (1 - pars$lambda1) * dnorm(simulated_data$x, mean = pars$a0, sd = pars$b0)
  gamma <- pars$lambda1 * dgamma(pmax(simulated_data$x, 1e-10), shape = pars$a1, scale = pars$b1)
  mixture <- normal + gamma 
  
  x = rep(simulated_data$x, 3)
  y = c(mixture, normal, gamma)
  component = rep(c("Mixture", "Normal", "Gamma"), each = length(simulated_data$x))
  
  data.frame(x = x, y = y, component = component)
}

# a) Whole population 
deng1_comp <- comp_density(deng1_par, sim_deng1)
deng2_comp <- comp_density(deng2_par, sim_deng2)
deng3_comp <- comp_density(deng3_par, sim_deng3)
deng4_comp <- comp_density(deng4_par, sim_deng4)

# b) Population ≤20 years 
youth_deng1_comp <- comp_density(youth_deng1_par, youth_sim_deng1)
youth_deng2_comp <- comp_density(youth_deng2_par, youth_sim_deng2)
youth_deng3_comp <- comp_density(youth_deng3_par, youth_sim_deng3)
youth_deng4_comp <- comp_density(youth_deng4_par, youth_sim_deng4)


# Plotting 

# Function to plot each FMM
plot_fmm <- function(data, figure_label, serotype) {
  ggplot(data, aes(x = x, y = y, color = component)) + 
    geom_line(size = 1.2) + 
    geom_area(data = subset(deng1_comp, component == "Mixture"), aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
    scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
    labs(x = figure_label, 
         y = "Density", 
         title = serotype) + 
    theme_minimal()
}

# a) Whole population 
plot_fmm(deng1_comp, "∆MFI", "DENV1 FMM")
plot_fmm(deng2_comp, "∆MFI", "DENV2 FMM")
plot_fmm(deng3_comp, "∆MFI", "DENV3 FMM")
plot_fmm(deng4_comp, "∆MFI", "DENV4 FMM")

# a) Population ≤ 20 years 
plot_fmm(youth_deng1_comp, "∆MFI", "DENV1 FMM (≤ 20 years)")
plot_fmm(youth_deng2_comp, "∆MFI", "DENV2 FMM (≤ 20 years)")
plot_fmm(youth_deng3_comp, "∆MFI", "DENV3 FMM (≤ 20 years)")
plot_fmm(youth_deng4_comp, "∆MFI", "DENV4 FMM (≤ 20 years)")


-------------------------------------------------------------------------------
  
# 6. Probability of boosting ---- 

# posterior probability = (likelihood x prior) / marginal likelihood 
# probability of boosting given a particular increase in MFI 

# function to determine the probability of infection
post_prob <- function(data, pars, id){
  
  x <- data 
  
  normal <- (1 - pars$lambda1) * dnorm(data, mean = pars$a0, sd = pars$b0)
  gamma <- pars$lambda1 * dgamma(pmax(data, 1e-10), shape = pars$a1, scale = pars$b1)
  
  marginal_lik <- normal + gamma 
  posterior <- gamma / marginal_lik
  
  data.frame(serosurvey_id = id, x = x, y = posterior)
  
}

# a) Whole population 
deng1_prob <- post_prob(inputs$change_dengns1_1_scale, deng1_par, inputs$serosurvey_id)
deng2_prob <- post_prob(inputs$change_dengns1_2_scale, deng2_par, inputs$serosurvey_id)
deng3_prob <- post_prob(inputs$change_dengns1_3_scale, deng3_par, inputs$serosurvey_id)
deng4_prob <- post_prob(inputs$change_dengns1_4_scale, deng4_par, inputs$serosurvey_id)

summary(deng1_prob$y)
summary(deng2_prob$y)
summary(deng3_prob$y)
summary(deng4_prob$y)

# b) Population ≤ 20 years
youth_deng1_prob <- post_prob(young$change_dengns1_1_scale, youth_deng1_par, young$serosurvey_id)
youth_deng2_prob <- post_prob(young$change_dengns1_2_scale, youth_deng2_par, young$serosurvey_id)
youth_deng3_prob <- post_prob(young$change_dengns1_3_scale, youth_deng3_par, young$serosurvey_id)
youth_deng4_prob <- post_prob(young$change_dengns1_4_scale, youth_deng4_par, young$serosurvey_id)

summary(youth_deng1_prob$y)
summary(youth_deng2_prob$y)
summary(youth_deng3_prob$y)
summary(youth_deng4_prob$y)


# Create a fine grid for plotting 
prob_curve <- function(data, pars, n = 1000) {
  
  data <- data[is.finite(data)]
  x    <- seq(min(data), max(data), length.out = n)
  
  normal       <- (1 - pars$lambda1) * dnorm(x, mean = pars$a0, sd = pars$b0)
  gamma        <- pars$lambda1 * dgamma(pmax(x, 1e-10), shape = pars$a1, scale = pars$b1)
  marginal_lik <- normal + gamma
  posterior    <- gamma / marginal_lik
  
  data.frame(x = x, y = posterior)
}

# a) Whole population — plot curves
deng1_prob_curve <- prob_curve(inputs$change_dengns1_1_scale, deng1_par)
deng2_prob_curve <- prob_curve(inputs$change_dengns1_2_scale, deng2_par)
deng3_prob_curve <- prob_curve(inputs$change_dengns1_3_scale, deng3_par)
deng4_prob_curve <- prob_curve(inputs$change_dengns1_4_scale, deng4_par)

# b) Population ≤ 20 years — plot curves
youth_deng1_prob_curve <- prob_curve(young$change_dengns1_1_scale, youth_deng1_par)
youth_deng2_prob_curve <- prob_curve(young$change_dengns1_2_scale, youth_deng2_par)
youth_deng3_prob_curve <- prob_curve(young$change_dengns1_3_scale, youth_deng3_par)
youth_deng4_prob_curve <- prob_curve(young$change_dengns1_4_scale, youth_deng4_par)

-------------------------------------------------------------------------------
  
# 7. Visualise boosting components & probability of boosting ---- 

# Function for plotting full fmm 
plot_model <- function(
    inputs,
    input_col,
    comp,
    prob,
    title,
    ymax     = 0.135,
    xlims    = c(-15, 55),
    binwidth = 4
) {
  
  x_thresh <- prob$x[which.min(abs(prob$y - 0.5))]
  
  ggplot() +
    geom_histogram(
      data = inputs,
      aes(x = {{ input_col }}, y = ..density..),
      fill = "lightgrey", color = "black", binwidth = binwidth
    ) +
    geom_line(
      data = subset(comp, component == "Normal"),
      aes(x = x, y = y), color = "red", size = 2
    ) +
    geom_line(
      data = subset(comp, component == "Gamma"),
      aes(x = x, y = y), color = "purple", size = 2
    ) +
    geom_line(
      data = subset(comp, component == "Mixture"),
      aes(x = x, y = y), color = "black", linetype = "dashed", size = 2
    ) +
    geom_line(
      data = prob,
      aes(x = x, y = y * ymax), size = 2, linetype = "dashed", color = "#7994B2"
    ) +
    geom_vline(xintercept = x_thresh, linetype = "solid", colour = "black", size = 1) +
    labs(x = "∆MFI", y = "Density", title = title) +
    theme_minimal(base_size = 26) +
    scale_y_continuous(
      name     = "Density",
      expand   = c(0, 0),
      sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xlims, ylim = c(0, ymax * 1.05)) +
    theme(
      legend.position    = "none",
      panel.grid         = element_blank(),
      axis.line          = element_line(colour = "black"),
      axis.ticks         = element_line(colour = "black"),
      axis.title.y.right = element_text(angle = -90, colour = "#7994B2"),
      axis.text.y.right  = element_text(colour = "#7994B2", hjust = 0.5),
      axis.ticks.y.right = element_line(colour = "#7994B2"),
      axis.line.y.right  = element_line(colour = "#7994B2"),
      plot.title         = element_text(size = 26)
    )
}

# a) Whole population 
p1 <- plot_model(inputs, change_dengns1_1_scale, deng1_comp, deng1_prob_curve, "DENV1")
p2 <- plot_model(inputs, change_dengns1_2_scale, deng2_comp, deng2_prob_curve, "DENV2")
p3 <- plot_model(inputs, change_dengns1_3_scale, deng3_comp, deng3_prob_curve, "DENV3")
p4 <- plot_model(inputs, change_dengns1_4_scale, deng4_comp, deng4_prob_curve, "DENV4")

ggsave(here("outputs", "deng1_fmm_prob_plot.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(here("outputs", "deng2_fmm_prob_plot.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(here("outputs", "deng3_fmm_prob_plot.png"), plot = p3, width = 10, height = 8, dpi = 300)
ggsave(here("outputs", "deng4_fmm_prob_plot.png"), plot = p4, width = 10, height = 8, dpi = 300)

# b) Population ≤ 20 years 
young_p1 <- plot_model(young, change_dengns1_1_scale, youth_deng1_comp, youth_deng1_prob_curve, "DENV1")
young_p2 <- plot_model(young, change_dengns1_2_scale, youth_deng2_comp, youth_deng2_prob_curve, "DENV2")
young_p3 <- plot_model(young, change_dengns1_3_scale, youth_deng3_comp, youth_deng3_prob_curve, "DENV3")
young_p4 <- plot_model(young, change_dengns1_4_scale, youth_deng4_comp, youth_deng4_prob_curve, "DENV4")

ggsave(here("outputs/figures/under20/", "young_deng1_fmm_prob_plot.png"), plot = young_p1, width = 10, height = 8, dpi = 300)
ggsave(here("outputs/figures/under20/", "young_deng2_fmm_prob_plot.png"), plot = young_p2, width = 10, height = 8, dpi = 300)
ggsave(here("outputs/figures/under20/", "young_deng3_fmm_prob_plot.png"), plot = young_p3, width = 10, height = 8, dpi = 300)
ggsave(here("outputs/figures/under20/", "young_deng4_fmm_prob_plot.png"), plot = young_p4, width = 10, height = 8, dpi = 300)



# Facet all the plots 
plot_facet_fmm <- function(
    data,
    comp_list,    # named list: list(DENV1 = deng1_comp, DENV2 = deng2_comp, ...)
    prob_list,    # named list: list(DENV1 = deng1_prob_curve, ...)
    ymax     = 0.135,
    xlims    = c(-15, 55),
    binwidth = 4
) {
  
  serotypes <- c("DENV1", "DENV2", "DENV3", "DENV4")
  cols      <- c("change_dengns1_1_scale", "change_dengns1_2_scale",
                 "change_dengns1_3_scale", "change_dengns1_4_scale")
  
  comp_all <- bind_rows(mapply(
    function(df, s) mutate(df, serotype = s),
    comp_list, serotypes, SIMPLIFY = FALSE
  ))
  
  prob_all <- bind_rows(mapply(
    function(df, s) mutate(df, serotype = s),
    prob_list, serotypes, SIMPLIFY = FALSE
  ))
  
  inputs_long <- bind_rows(mapply(
    function(col, s) mutate(data, serotype = s, value = .data[[col]]),
    cols, serotypes, SIMPLIFY = FALSE
  ))
  
  thresholds <- bind_rows(mapply(
    function(df, s) data.frame(serotype = s, x = df$x[which.min(abs(df$y - 0.5))]),
    prob_list, serotypes, SIMPLIFY = FALSE
  ))
  
  # Relabel for display
  label_map <- c("DENV1" = "DENV-1", "DENV2" = "DENV-2",
                 "DENV3" = "DENV-3", "DENV4" = "DENV-4")
  
  comp_all$serotype    <- label_map[comp_all$serotype]
  prob_all$serotype    <- label_map[prob_all$serotype]
  inputs_long$serotype <- label_map[inputs_long$serotype]
  thresholds$serotype  <- label_map[thresholds$serotype]
  
  ggplot() +
    geom_histogram(
      data = inputs_long,
      aes(x = value, y = ..density..),
      fill = "lightgrey", color = "black", binwidth = binwidth
    ) +
    geom_line(
      data = subset(comp_all, component == "Normal"),
      aes(x = x, y = y), color = "red", size = 2
    ) +
    geom_line(
      data = subset(comp_all, component == "Gamma"),
      aes(x = x, y = y), color = "purple", size = 2
    ) +
    geom_line(
      data = subset(comp_all, component == "Mixture"),
      aes(x = x, y = y), color = "black", linetype = "dashed", size = 2
    ) +
    geom_line(
      data = prob_all,
      aes(x = x, y = y * ymax), color = "#7994B2", linetype = "dashed", size = 2
    ) +
    geom_vline(data = thresholds, aes(xintercept = x), color = "black", size = 1) +
    facet_wrap(~ serotype, ncol = 2) +
    scale_y_continuous(
      name     = "Density",
      expand   = c(0, 0),
      sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = xlims, ylim = c(0, ymax * 1.05)) +
    labs(x = "∆MFI", y = "Density") +
    theme_minimal(base_size = 18) +
    theme(
      legend.position    = "none",
      panel.grid         = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.text         = element_text(size = 18),
      strip.background   = element_rect(fill = "grey90", colour = "black"),
      axis.line          = element_line(colour = "black"),
      axis.ticks         = element_line(colour = "black"),
      axis.title.y.right = element_text(angle = -90, colour = "#7994B2"),
      axis.text.y.right  = element_text(colour = "#7994B2"),
      axis.ticks.y.right = element_line(colour = "#7994B2"),
      axis.line.y.right  = element_line(colour = "#7994B2")
    )
}

# a) Whole population
facet_deng <- plot_facet_fmm(
  data      = inputs,
  comp_list = list(DENV1 = deng1_comp, DENV2 = deng2_comp, DENV3 = deng3_comp, DENV4 = deng4_comp),
  prob_list = list(DENV1 = deng1_prob_curve, DENV2 = deng2_prob_curve, DENV3 = deng3_prob_curve, DENV4 = deng4_prob_curve))

ggsave(here("outputs", "facet_fmm_plot.png"), width = 14, height = 10, dpi = 300)


# b) Population (≤ 20 years)
young_facet_deng <- plot_facet_fmm(
  data      = young,
  comp_list = list(DENV1 = youth_deng1_comp, DENV2 = youth_deng2_comp, DENV3 = youth_deng3_comp, DENV4 = youth_deng4_comp),
  prob_list = list(DENV1 = youth_deng1_prob_curve, DENV2 = youth_deng2_prob_curve, DENV3 = youth_deng3_prob_curve, DENV4 = youth_deng4_prob_curve))

ggsave(here("outputs/figures/under20/", "young_facet_fmm_plot.png"), plot = young_facet_deng, width = 14, height = 10, dpi = 300)


-------------------------------------------------------------------------------

# 7. Join data ----

# need to assign the probabilities of boosting back to each individual 

# a) Whole population 
inputs <- inputs %>% 
  left_join(deng1_prob %>% select(serosurvey_id, deng1_prob = y), by = "serosurvey_id") %>% 
  left_join(deng2_prob %>% select(serosurvey_id, deng2_prob = y), by = "serosurvey_id") %>% 
  left_join(deng3_prob %>% select(serosurvey_id, deng3_prob = y), by = "serosurvey_id") %>% 
  left_join(deng4_prob %>% select(serosurvey_id, deng4_prob = y), by = "serosurvey_id") 

# b) Population ≤20 years 

-------------------------------------------------------------------------------

# 9. Save data ---- 

inputs_dr <- inputs 
# saveRDS(inputs, here("data", "inputs_dr.rds"))


