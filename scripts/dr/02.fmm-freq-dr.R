
# 02. fmm-freq ----

# finite mixture model in a frequentist framework for the Dominican Republic data 

-------------------------------------------------------------------------------

# 1. Install and load packages ----

library("ggplot2")
library("dplyr")
library("here")

-------------------------------------------------------------------------------
  
# 2. Load data ----

df_wide <- readRDS(here("data", "df_wide.rds"))

----------------------------------------------------------------------------------------------

# 3. Data wrangling ----

# create a new df to edit
inputs <- df_wide

# scale the data 

# changes in MFI 
inputs$change_dengns1_1_scale <- inputs$change_dengns1_1/1000
inputs$change_dengns1_2_scale <- inputs$change_dengns1_2/1000
inputs$change_dengns1_3_scale <- inputs$change_dengns1_3/1000
inputs$change_dengns1_4_scale <- inputs$change_dengns1_4/1000

--------------------------------------------------------------------------

# 4. Model set up ----
  
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

--------------------------------------------------------------------------
  
# 5. Model fitting  ----

# provide the starting parameters 
start_par <- c(a0 = 0, b0 = 1, a1 = 2, b1 = 3, lambda1 = 1)

# fit fmm
deng1_fmm <- fit_pathogen(inputs$change_dengns1_1_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng2_fmm <- fit_pathogen(inputs$change_dengns1_2_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng3_fmm <- fit_pathogen(inputs$change_dengns1_3_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))
deng4_fmm <- fit_pathogen(inputs$change_dengns1_4_scale, c(b0 = 1, a1 = 2, b1 = 3, lambda1 = 1))

# store parameters 
deng1_par <- deng1_fmm$par
deng2_par <- deng2_fmm$par
deng3_par <- deng3_fmm$par
deng4_par <- deng4_fmm$par

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

deng1_par <- extract_params(deng1_par)
deng2_par <- extract_params(deng2_par)
deng3_par <- extract_params(deng3_par)
deng4_par <- extract_params(deng4_par)

# model components 
mixture_density <- function(x, pars) {
  
  normal <- dnorm(x, mean = pars$a0, sd = pars$b0)
  gamma  <- dgamma(x, shape = pars$a1, scale = pars$b1)
  
  ((1 - pars$lambda1) * normal) + (pars$lambda1 * gamma)
}

# simulate data based off the optimised parameters 
simulate_mixture <- function(data, pars, n = 1000) {
  
  data <- data[is.finite(data)] 
  
  x_vals <- seq(min(data), max(data), length.out = n)
  y_vals <- mixture_density(x_vals, pars)
  
  data.frame(x = x_vals, y = y_vals)
}

sim_deng1 <- simulate_mixture(inputs$change_dengns1_1_scale, deng1_par)
sim_deng2 <- simulate_mixture(inputs$change_dengns1_2_scale, deng2_par)
sim_deng3 <- simulate_mixture(inputs$change_dengns1_3_scale, deng3_par)
sim_deng4 <- simulate_mixture(inputs$change_dengns1_4_scale, deng4_par)

# plot distribution 

# deng1 
ggplot(data = sim_deng1, aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "DENV1 FMM") + 
  theme_minimal()

# deng2 
ggplot(data = sim_deng2, aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "DENV2 FMM") + 
  theme_minimal()

# deng3 
ggplot(data = sim_deng3, aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "DENV3 FMM") + 
  theme_minimal()

# deng4 
ggplot(data = sim_deng4, aes(x = x, y = y)) + 
  geom_area(alpha = 0.7, fill = "lightgrey") + 
  labs(x = "∆MFI", 
       y = "Density", 
       title = "DENV4 FMM") + 
  theme_minimal()


# function to create a df for individual component densities
comp_density <- function(pars, simulated_data){
  
  normal <- (1 - pars$lambda1) * dnorm(simulated_data$x, mean = pars$a0, sd = pars$b0)
  gamma <- pars$lambda1 * dgamma(pmax(simulated_data$x, 1e-10), shape = pars$a1, scale = pars$b1)
  mixture <- normal + gamma 
  
  x = rep(simulated_data$x, 3)
  y = c(mixture, normal, gamma)
  component = rep(c("Mixture", "Normal", "Gamma"), each = length(simulated_data$x))
  
  data.frame(x = x, y = y, component = component)
}

deng1_comp <- comp_density(deng1_par, sim_deng1)
deng2_comp <- comp_density(deng2_par, sim_deng2)
deng3_comp <- comp_density(deng3_par, sim_deng3)
deng4_comp <- comp_density(deng4_par, sim_deng4)

# plot components 

# deng1 
ggplot(deng1_comp, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng1_comp, component == "Mixture"), aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "FMM DENV1"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# deng2 
ggplot(deng2_comp, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng2_comp, component == "Mixture"), aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "FMM DENV2"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# deng3 
ggplot(deng3_comp, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng3_comp, component == "Mixture"), aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "FMM DENV3"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

# deng4 
ggplot(deng4_comp, aes(x = x, y = y, color = component)) +
  geom_line(size = 1.2) +
  geom_area(data = subset(deng4_comp, component == "Mixture"), aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "FMM DENV4"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(-20, 50), ylim = c(0, 0.1))

--------------------------------------------------------------------------
  
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

deng1_prob <- post_prob(inputs$change_dengns1_1_scale, deng1_par, inputs$serosurvey_id)
deng2_prob <- post_prob(inputs$change_dengns1_2_scale, deng2_par, inputs$serosurvey_id)
deng3_prob <- post_prob(inputs$change_dengns1_3_scale, deng3_par, inputs$serosurvey_id)
deng4_prob <- post_prob(inputs$change_dengns1_4_scale, deng4_par, inputs$serosurvey_id)

summary(deng1_prob$y)
summary(deng2_prob$y)
summary(deng3_prob$y)
summary(deng4_prob$y)

# scaling factor 
ymax <- 0.125

# deng1 
x_deng1 <- deng1_prob$x[which.min(abs(deng1_prob$y - 0.5))]

ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_1_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = subset(deng1_comp, component == "Normal"), aes(x = x, y = y), color = "red", size = 2) +
  geom_line(data = subset(deng1_comp, component == "Gamma"), aes(x = x, y = y), color = "purple", size = 2) + 
  geom_line(data = subset(deng1_comp, component == "Mixture"), aes(x = x, y = y), color = "black", linetype = "dashed", size = 2) +
  geom_line(data = deng1_prob, aes(x = x, y = y * ymax), size = 2, linetype = "dashed", color = "#7994B2") + 
  geom_vline(xintercept = x_deng1, linetype = "solid", colour = "black", size = 1) +
  labs(x = "∆MFI", y = "Density", title = "DENV1") +
  theme_minimal(base_size = 26) + 
  scale_y_continuous(
    name = "Density",
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-15, 55), ylim = c(0, ymax * 1.05)) +
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

ggsave(here("outputs", "deng1_fmm_prob_plot.png"), width = 10, height = 8, dpi = 300)


# deng2 
x_deng2 <- deng2_prob$x[which.min(abs(deng2_prob$y - 0.5))]

ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_2_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng2_comp, aes(x = x, y = y, color = component), size = 1) + 
  geom_line(data = subset(deng2_comp, component == "Normal"), aes(x = x, y = y), color = "red", size = 2) +
  geom_line(data = subset(deng2_comp, component == "Gamma"), aes(x = x, y = y), color = "purple", size = 2) + 
  geom_line(data = subset(deng2_comp, component == "Mixture"), aes(x = x, y = y), color = "black", linetype = "dashed", size = 2) +
  geom_line(data = deng2_prob, aes(x = x, y = y * ymax), size = 2, linetype = "dashed", color = "#7994B2") + 
  geom_vline(xintercept = x_deng2, linetype = "solid", colour = "black", size = 1) +
  labs(x = "∆MFI", y = "Density", title = "DENV2") +
  theme_minimal(base_size = 26) + 
  scale_y_continuous(
    name = "Density",
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-15, 55), ylim = c(0, ymax * 1.05)) +
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

ggsave(here("outputs", "deng2_fmm_prob_plot.png"), width = 10, height = 8, dpi = 300)

# deng3 
x_deng3 <- deng3_prob$x[which.min(abs(deng3_prob$y - 0.5))]

ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_3_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng3_comp, aes(x = x, y = y, color = component), size = 1) + 
  geom_line(data = subset(deng3_comp, component == "Normal"), aes(x = x, y = y), color = "red", size = 2) +
  geom_line(data = subset(deng3_comp, component == "Gamma"), aes(x = x, y = y), color = "purple", size = 2) + 
  geom_line(data = subset(deng3_comp, component == "Mixture"), aes(x = x, y = y), color = "black", linetype = "dashed", size = 2) +
  geom_line(data = deng3_prob, aes(x = x, y = y * ymax), size = 2, linetype = "dashed", color = "#7994B2") + 
  geom_vline(xintercept = x_deng3, linetype = "solid", colour = "black", size = 1) +
  labs(x = "∆MFI", y = "Density", title = "DENV3") +
  theme_minimal(base_size = 26) + 
  scale_y_continuous(
    name = "Density",
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-15, 55), ylim = c(0, ymax * 1.05)) +
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

ggsave(here("outputs", "deng3_fmm_prob_plot.png"), width = 10, height = 8, dpi = 300)

# deng4 
x_deng4 <- deng4_prob$x[which.min(abs(deng4_prob$y - 0.5))]

ggplot() +
  geom_histogram(data = inputs, aes(x = change_dengns1_4_scale, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) + 
  geom_line(data = deng4_comp, aes(x = x, y = y, color = component), size = 1) + 
  geom_line(data = subset(deng4_comp, component == "Normal"), aes(x = x, y = y), color = "red", size = 2) +
  geom_line(data = subset(deng4_comp, component == "Gamma"), aes(x = x, y = y), color = "purple", size = 2) + 
  geom_line(data = subset(deng4_comp, component == "Mixture"), aes(x = x, y = y), color = "black", linetype = "dashed", size = 2) +
  geom_line(data = deng4_prob, aes(x = x, y = y * ymax), size = 2, linetype = "dashed", color = "#7994B2") + 
  geom_vline(xintercept = x_deng4, linetype = "solid", colour = "black", size = 1) +
  labs(x = "∆MFI", y = "Density", title = "DENV4") +
  theme_minimal(base_size = 26) + 
  scale_y_continuous(
    name = "Density",
    expand = c(0, 0),
    sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-15, 55), ylim = c(0, ymax * 1.05)) +
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

ggsave(here("outputs", "deng4_fmm_prob_plot.png"), width = 10, height = 8, dpi = 300)


# facet all the plots 

# make a df of all the component data 
comp_all <- bind_rows(
  mutate(deng1_comp, serotype = "DENV1"),
  mutate(deng2_comp, serotype = "DENV2"),
  mutate(deng3_comp, serotype = "DENV3"),
  mutate(deng4_comp, serotype = "DENV4")
)

# make a df of all the probability data 
prob_all <- bind_rows(
  mutate(deng1_prob, serotype = "DENV1"),
  mutate(deng2_prob, serotype = "DENV2"),
  mutate(deng3_prob, serotype = "DENV3"),
  mutate(deng4_prob, serotype = "DENV4")
)

# join the inputs data into a long form 
inputs_long <- bind_rows(
  mutate(inputs, serotype = "DENV1", value = change_dengns1_1_scale),
  mutate(inputs, serotype = "DENV2", value = change_dengns1_2_scale),
  mutate(inputs, serotype = "DENV3", value = change_dengns1_3_scale),
  mutate(inputs, serotype = "DENV4", value = change_dengns1_4_scale)
)

# add the threshold lines for p > 0.5
thresholds <- data.frame(
  serotype = c("DENV1", "DENV2", "DENV3", "DENV4"),
  x = c(x_deng1, x_deng2, x_deng3, x_deng4)
)

serotype_map <- c(
  "DENV1" = "DENV-1",
  "DENV2" = "DENV-2",
  "DENV3" = "DENV-3",
  "DENV4" = "DENV-4"
)

comp_all$serotype    <- serotype_map[comp_all$serotype]
prob_all$serotype    <- serotype_map[prob_all$serotype]
inputs_long$serotype <- serotype_map[inputs_long$serotype]
thresholds$serotype  <- serotype_map[thresholds$serotype]


# set max value 
ymax <- 0.10

ggplot() +
  geom_histogram(data = inputs_long, aes(x = value, y = ..density..), fill = "lightgrey", color = "black", binwidth = 4) +
  geom_line(data = subset(comp_all, component == "Normal"), aes(x = x, y = y), color = "red", size = 2) +
  geom_line(data = subset(comp_all, component == "Gamma"), aes(x = x, y = y), color = "purple", size = 2) +
  geom_line(data = subset(comp_all, component == "Mixture"), aes(x = x, y = y), color = "black", linetype = "dashed", size = 2) +
  geom_line(data = prob_all, aes(x = x, y = y * ymax), color = "#7994B2", linetype = "dashed", size = 2) +
  geom_vline(data = thresholds, aes(xintercept = x), color = "black", size = 1) +
  facet_wrap(~ serotype, ncol = 2) +
  scale_y_continuous(name = "Density", limits = c(0, ymax * 1.05), expand = c(0, 0), sec.axis = sec_axis(~ . / ymax, name = "Boosting Probability")) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-15, 55)) +
  labs(x = "∆MFI", y = "Density") +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.y.right = element_text(angle = -90, colour = "#7994B2"),
    axis.text.y.right  = element_text(colour = "#7994B2"),
    axis.ticks.y.right = element_line(colour = "#7994B2"),
    axis.line.y.right  = element_line(colour = "#7994B2"),
    plot.title = element_text(size = 18)
  )

ggsave(here("outputs", "facet_fmm_plot.png"), width = 14, height = 10, dpi = 300)


--------------------------------------------------------------------------

# 7. Join data ----

# need to assign the probabilities of boosting back to each individual 
inputs <- inputs %>% 
  left_join(deng1_prob %>% select(serosurvey_id, deng1_prob = y), by = "serosurvey_id") %>% 
  left_join(deng2_prob %>% select(serosurvey_id, deng2_prob = y), by = "serosurvey_id") %>% 
  left_join(deng3_prob %>% select(serosurvey_id, deng3_prob = y), by = "serosurvey_id") %>% 
  left_join(deng4_prob %>% select(serosurvey_id, deng4_prob = y), by = "serosurvey_id") 

# split into dengue boosters and non-boosters 

--------------------------------------------------------------------------

# 9. Save data ---- 

inputs_dr <- inputs 
saveRDS(inputs, here("data", "inputs_dr.rds"))




