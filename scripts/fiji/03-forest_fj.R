
# Forest Plots 
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
library(colorspace)

------------------------------------------------------------------------
  
# 2. Data ---- 

inputs <- readRDS(here("data", "inputs_fj.rds"))

deng1_n97_lr <- read.csv(here("outputs", "lr_deng1_n97.csv"))
deng3_n97_lr <- read.csv(here("outputs", "lr_deng3_n97.csv"))
deng1_n260_lr <- read.csv(here("outputs", "lr_deng1_n260.csv"))
deng3_n260_lr <- read.csv(here("outputs", "lr_deng3_n260.csv"))
deng1_n97_flr <- read.csv(here("outputs", "flr_deng1_n97.csv"))
deng3_n97_flr <- read.csv(here("outputs", "flr_deng3_n97.csv"))
deng1_n260_flr <- read.csv(here("outputs", "flr_deng1_n260.csv"))
deng3_n260_flr <- read.csv(here("outputs", "flr_deng3_n260.csv"))

lr_deng_n97 <- read.csv(here("outputs", "lr_deng_n97.csv"))
lr_deng_n260 <- read.csv(here("outputs", "lr_deng_n260.csv"))
flr_deng_n97 <- read.csv(here("outputs", "flr_deng_n97.csv"))
flr_deng_n260 <- read.csv(here("outputs", "flr_deng_n260.csv"))

---------------------------------------------------------------------------------------
  
# Serotype-specific forest plots ----  

# function to combine the different models 
add_or_ci <- function(df, outcome, model, sample) {
  df %>%
    mutate(
      OR  = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
      LCL = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
      UCL = as.numeric(str_extract(OR_CI, "(?<=-)[0-9.]+")),
      outcome = outcome,
      model   = model,
      sample  = sample
    )
}

# combine all models into a single df 
all_df <- bind_rows(
  add_or_ci(deng1_n97_lr,    outcome="DENV1", model="LR",  sample="Reduced"),
  add_or_ci(deng3_n97_lr,    outcome="DENV3", model="LR",  sample="Reduced"),
  add_or_ci(deng1_n260_lr,  outcome="DENV1", model="LR",  sample="Total"),
  add_or_ci(deng3_n260_lr,  outcome="DENV3", model="LR",  sample="Total"),
  add_or_ci(deng1_n97_flr,   outcome="DENV1", model="FLR", sample="Reduced"),
  add_or_ci(deng3_n97_flr,   outcome="DENV3", model="FLR", sample="Reduced"),
  add_or_ci(deng1_n260_flr, outcome="DENV1", model="FLR", sample="Total"),
  add_or_ci(deng3_n260_flr, outcome="DENV3", model="FLR", sample="Total")
) %>%
  mutate(var = factor(var, levels = rev(unique(var))))

# plot the forest plot 
forest_4p <- 
  ggplot(all_df, aes(x = OR, y = var, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 position = position_dodge(width = 0.6),
                 height = 0.2) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_x_log10() +
  facet_grid(outcome ~ sample, switch = "y") +
  labs(
    x = "Odds Ratio (log scale)",
    y = "",
    colour = "Model",
    title = "Risk factors for DENV1 and DENV3 infection in Fiji (2013 - 2015)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    
    panel.border = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    panel.spacing = unit(1.2, "lines"),
    
    strip.background = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) 

forest_4p
ggsave(here("outputs", "forest_4p.png"), width = 8, height = 6, dpi = 300)

---------------------------------------------------------------------------------------

# Forest plot without serotype ---- 

# just plot logistic regression output 
collapsed_df <- bind_rows(
  add_or_ci(lr_deng_n97 , outcome="DENV", model="LR", sample="Reduced"),
  add_or_ci(lr_deng_n260, outcome="DENV", model="LR", sample="Total")
  # ,
  # add_or_ci(flr_deng_n97, outcome="DENV", model="FLR", sample="Reduced"),
  # add_or_ci(flr_deng_n260, outcome="DENV", model="FLR", sample="Total")
) %>%
  mutate(var = factor(var, levels = rev(unique(var))))

# add outcome for facet 
collapsed_df$outcome <- "DENV"

# rename labels 
var_labels <- c(
  AGE_U_20 = "Age under 20",
  SEX = "Male",
  ETHNIC = "iTaukei ethnicity",
  I_MOS = "Mosquito exposure",
  I_TIR = "Used car tires",
  I_WAT = "Open water container(s)",
  I_AC = "Air conditioning",
  I_BLK = "Blocked drains",
  GEOG = "Urban or peri-urban",
  FEVER_2YR = "Fever (past 2 years)",
  DOC_2YR = "Doctor visit (past 2 years)",
  HH_D = "Cohabitant doctor visit (past 2 years)"
)

# forest plot 
forest_collapsed_lr_fiji <- 
  ggplot(collapsed_df, aes(x = OR, y = var, colour = sample)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), position = position_dodge(width = 0.6), height = 0.2) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_y_discrete(labels = var_labels) + 
  scale_x_log10() +
  facet_grid(. ~ outcome) +   
  labs(
    x = "Odds Ratio (log scale)",
    y = "",
    colour = "Data"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    strip.background = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(0, "lines"), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

forest_collapsed_lr_fiji
ggsave(here("outputs", "forest_collapsed_lr_fiji.png"), width = 12, height = 6, dpi = 300)

---------------------------------------------------------------------------------------

# data from kucharski et al. (2018) ----

ak_df <- data.frame(
  var = c("AGE_U_20", "SEX", "ETHNIC","I_MOS","I_TIR","I_WAT","I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR", "HH_D"),
  number = c(61, 49, 85, 90, 61, 61, 23, 53, 50, 20, 16, 9),
  OR_CI = c(
    "0.49 (0.21-1.13)",
    "0.81 (0.36-1.84)",
    "1.33 (0.39-5.32)",
    "4.19 (0.68-80.85)",
    "1.80 (0.77-4.42)",
    "1.49 (0.64-3.58)",
    "0.46 (0.15-1.26)",
    "1.04 (0.46-2.38)",
    "2.18 (0.95-5.11)",
    "2.94 (1.08-8.38)",
    "3.15 (1.06-10.13)",
    "2.08 (0.52-8.94)"
  ),
  p_value = c(0.10, 0.62, 0.66, 0.19, 0.18, 0.37, 0.15, 0.92, 0.07, 0.04, 0.04, 0.30),
  stringsAsFactors = FALSE
) 

ak_df <- add_or_ci(ak_df , outcome="DENV", model="LR", sample="Reduced") %>%
  mutate(var = factor(var, levels = rev(unique(var))))

# plot forest plot of AK data 
forest_ak <- 
  ggplot(ak_df, aes(x = OR, y = var, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 position = position_dodge(width = 0.6),
                 height = 0.2) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_x_log10() +
  facet_grid(outcome ~ sample, switch = "y") +
  scale_color_manual(values = "#00BFC4") + 
  labs(
    x = "Odds Ratio (log scale)",
    y = "",
    colour = "Model",
    title = "Risk factors for DENV infection from Kucharski et al. (2018)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    
    panel.border = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    panel.spacing = unit(1.2, "lines"),
    
    strip.background = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

forest_ak
ggsave(here("outputs", "forest_ak.png"), width = 8, height = 6, dpi = 300)

---------------------------------------------------------------------------------------
  
# Integrating different forest plots ---- 

# add a source to each set of data 

# all dengue 
collapsed_plot_df <- collapsed_df %>%
  mutate(source = "Current")

# AK plot 
ak_plot_df <- ak_df %>%
  mutate(
    source = "Kucharski et al.",
    model  = "LR"   # Kucharski is LR only, so make it explicit
  )

# combine the two 
combined_df <- bind_rows(collapsed_plot_df, ak_plot_df)
var_levels <- combined_df %>%
  distinct(var) %>%
  pull(var) %>%
  as.character()

combined_df <- combined_df %>%
  mutate(var = factor(as.character(var), levels = rev(var_levels)))

pd <- position_dodge(width = 0.75)

# plot the integrated forest plot 
forest_compare <-
  ggplot(combined_df,
         aes(x = OR, y = var,
             colour = model,
             shape  = source,
             group  = interaction(model, source))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL),
                 position = pd, height = 0.2) +
  geom_point(position = pd, size = 2.8) +
  scale_x_log10() +
  facet_grid(outcome ~ sample,
             switch = "y",
             labeller = labeller(sample = c("Reduced" = "Seronegative",
                                            "Total"   = "All participants"))) +
  labs(
    x = "Odds Ratio (log scale)",
    y = "",
    colour = "Model",
    shape  = "Study",
    title = "Risk factors for DENV infection: This study vs Kucharski et al. (2018)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    panel.spacing = unit(1.2, "lines"),
    strip.background = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) + 
  scale_y_discrete(labels = c(
    "AGE_U_20" = "AGE (U20)",
    "ETHNIC" = "ETHNICITY",
    "I_MOS" = "MOSQUITO",
    "I_TIR" = "TIRES",
    "I_WAT" = "WATER",
    "I_AC" = "AC",
    "I_BLK" = "BLOCKAGE",
    "GEOG" = "LOCATION",
    "FEVER_2YR" = "FEVER",
    "DOC_2YR" = "DOCTOR",
    "HH_D" = "HH DOCTOR"
  ))

forest_compare
ggsave(here("outputs", "forest_compare.png"), width = 13, height = 7, dpi = 300)

