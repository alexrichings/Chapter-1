
# Attic ----

# 02.fmm-freq ----

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

--------------------------------------------------------------------------
  


