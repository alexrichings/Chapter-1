# Cross-reactivity 

library(ggplot2)
-----------------------------------------------------------------------------------------------------------
# Load in the data ----
inf_prob <- readRDS("Data/DR_surveillance/Serology/inf_prob")

-----------------------------------------------------------------------------------------------------------
# 1. Correlation between boosts ----

# look at how the probability of infection changes between serotypes and viruses 

# use ghrexplore 
library(GHRexplore)

# correlation of boosts 
plot_correlation(inf_prob, 
                 var = c("deng1_boost","deng2_boost","deng3_boost","deng4_boost","zika_boost"),
                 plot_type = c("raster", "number"),
                 palette = "RdBu",
                 method = "pearson") 

# boosts in deng2/deng3 are completely correlated 
# zika shows the lowest correlation boosts

# correlation of start values 
plot_correlation(inf_prob, 
                 var = c("deng1_start","deng2_start","deng3_start","deng4_start","zika_start"),
                 plot_type = c("raster", "number"),
                 palette = "RdBu",
                 method = "pearson") 

# the correlation structure remains identical as with the boosts, which suggests that the similarities 
# between the MFI data is driven more so by cross-reactivity because heterogeneity of exposure in a set 
# time frame would disrupt correlations between serotypes/viruses 

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
  xlim(-70, 70) + 
  ylim(0, 0.03)

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
  xlim(-70, 70) + 
  ylim(0, 0.03)

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
  xlim(-70, 70) + 
  ylim(0, 0.03)

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
  xlim(-70, 70) + 
  ylim(0, 0.03)

# zika ----
ggplot(zika_comp_df, aes(x = x, y = y, color = component)) +
  geom_line(linewidth = 1.2) +
  geom_area(data = subset(zika_comp_df, component == "Mixture"), 
            aes(x = x, y = y), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
  scale_color_manual(values = c("Mixture" = "black", "Normal" = "red", "Gamma" = "purple")) +
  labs(
    x = "∆MIA",
    y = "Density",
    title = "Two-component mixture model for MIA Zika data"
  ) +
  theme_minimal() + 
  xlim(-70, 70) + 
  ylim(0, 0.03)

-----------------------------------------------------------------------------------------------------------
# 2. Components of cross-reactivity 
  
# there was no circulation of zika cirus in the dominican republic between 2021 and 2022
# yet there is signal for boost which indicates cross-reactivity during serological testing 
  
# assume that none of the dengue serotypes cross-react, only dengue-zika cross-reaction occurs 
  
# zika MFI with no circulation but XR: 
  # ∆MFI_z = µ_z + xr_boosts 

  
  
 
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  

