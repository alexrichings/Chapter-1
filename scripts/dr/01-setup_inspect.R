
# Chapter 1 

--------------------------------------------------------------------------------

# Script: 01-setup_inspect ----  

# Purpose: 
# Import raw serological data, clean variables, visual inspection of data

# Output: 
# Paired serological data file: df_wide 

--------------------------------------------------------------------------------

# Load packages ----

library("readxl")
library("here")
library("ggplot2")
library("dplyr")

-------------------------------------------------------------------------------
# Section 1: Load and manage data ----
-------------------------------------------------------------------------------
  
# 1.1. Import raw data ----

# df_0: national serological survey; df_1: follow up cohort 
df_0 <- read_excel(here("data", "Copy of DR_serosurvey_MASTER_with_MBA_UPDATED.xlsx"))
df_1 <- read_excel(here("data", "Copy of FINAL_COHORT_PLUS_MBA.xlsx"))
                   
-------------------------------------------------------------------------------

# 1.2. Generate paired data set & rename/mutate variables ---- 

# round 1: x_sero 
# round 2: x_cohort 

# Include the specific variables to carry through 
df_wide <- df_0 %>%
  select(serosurvey_id, age, gender, setting, occupation, education, province, 
         region, occupation_2, education2) %>%
  rename(age_sero = age, gender_sero = gender, setting_sero = setting, 
         occupation_sero = occupation, education_sero = education, 
         occupation2_sero = occupation_2, education2_sero = education2) %>%
  inner_join(
    df_1 %>% 
      select(
        serosurvey_id, cohort_id, date_sero, date_cohort, age, gender, occupation, date_diff, 
        dengns1_1_mfi_sero, dengns1_2_mfi_sero, dengns1_3_mfi_sero, dengns1_4_mfi_sero, zika_ns1_mfi_sero, 
        chik_e1_mfi_sero, dengns1_1_mfi_cohort, dengns1_2_mfi_cohort,dengns1_3_mfi_cohort, dengns1_4_mfi_cohort, 
        zika_ns1_mfi_cohort, chik_e1_mfi_cohort, dengns1_1_seropos_sero, dengns1_2_seropos_sero, 
        dengns1_3_seropos_sero, dengns1_4_seropos_sero, zika_ns1_seropos_sero, 
        chik_e1_seropos_sero,dengns1_1_seropos_cohort, dengns1_2_seropos_cohort, dengns1_3_seropos_cohort, 
        dengns1_4_seropos_cohort, zika_ns1_seropos_cohort, chik_e1_seropos_cohort) %>%
      rename(age_cohort = age, gender_cohort = gender, occupation_cohort = occupation), 
    by = "serosurvey_id") %>% 
  relocate(age_cohort, .after = age_sero) %>%
  mutate(
    change_dengns1_1 = dengns1_1_mfi_cohort - dengns1_1_mfi_sero,
    change_dengns1_2 = dengns1_2_mfi_cohort - dengns1_2_mfi_sero,
    change_dengns1_3 = dengns1_3_mfi_cohort - dengns1_3_mfi_sero,
    change_dengns1_4 = dengns1_4_mfi_cohort - dengns1_4_mfi_sero,
    change_zika_ns1  = zika_ns1_mfi_cohort  - zika_ns1_mfi_sero,
    change_chik_e1   = chik_e1_mfi_cohort   - chik_e1_mfi_sero
  )

# Save the data 
# saveRDS(df_wide, here("data", "df_wide.rds"))

-------------------------------------------------------------------------------
# Section 2: Inspect the data ----
-------------------------------------------------------------------------------
  
# 2.1. Raw titres ----

# Function to plot mfi
plot_mfi <- function(data, input, figure_label, serotype, color) {
  ggplot(data, aes({{ input }})) + 
    geom_density(stat = "density", alpha = 0.7, fill = color) + 
    labs(x = figure_label, y = "Density", title = serotype) + 
    theme_minimal()
}

# Cut-off values for dengue serotypes 
# DENV-1: ≥739; DENV-2: ≥1777; DENV-3: ≥451; DENV-4: ≥420 

# Plot baseline titre
plot_mfi(df_wide, dengns1_1_mfi_sero, "MFI", "DENV1", "skyblue")
plot_mfi(df_wide, dengns1_2_mfi_sero, "MFI", "DENV2", "skyblue")
plot_mfi(df_wide, dengns1_3_mfi_sero, "MFI", "DENV3", "skyblue")
plot_mfi(df_wide, dengns1_4_mfi_sero, "MFI", "DENV4", "skyblue")

# Plot follow-up 
plot_mfi(df_wide, dengns1_1_mfi_cohort, "MFI", "DENV1", "skyblue")
plot_mfi(df_wide, dengns1_2_mfi_cohort, "MFI", "DENV2", "skyblue")
plot_mfi(df_wide, dengns1_3_mfi_cohort, "MFI", "DENV3", "skyblue")
plot_mfi(df_wide, dengns1_4_mfi_cohort, "MFI", "DENV4", "skyblue")

-------------------------------------------------------------------------------
  
# 2.2. Change in titre ----

# Plot ∆MFI
plot_mfi(df_wide, change_dengns1_1, "∆MFI", "Change in DENV1", "skyblue")
plot_mfi(df_wide, change_dengns1_2, "∆MFI", "Change in DENV2", "skyblue")
plot_mfi(df_wide, change_dengns1_3, "∆MFI", "Change in DENV3", "skyblue")
plot_mfi(df_wide, change_dengns1_4, "∆MFI", "Change in DENV4", "skyblue")

-------------------------------------------------------------------------------
  
# 2.3. Inspect for under 20 years ----

# Filter by baseline age 
young <- df_wide |> filter(age_sero <= 20)

# Plot baseline titre
plot_mfi(young, dengns1_1_mfi_sero, "MFI", "DENV1", "skyblue")
plot_mfi(young, dengns1_2_mfi_sero, "MFI", "DENV2", "skyblue")
plot_mfi(young, dengns1_3_mfi_sero, "MFI", "DENV3", "skyblue")
plot_mfi(young, dengns1_4_mfi_sero, "MFI", "DENV4", "skyblue")

# Plot follow-up 
plot_mfi(young, dengns1_1_mfi_cohort, "MFI", "DENV1", "skyblue")
plot_mfi(young, dengns1_2_mfi_cohort, "MFI", "DENV2", "skyblue")
plot_mfi(young, dengns1_3_mfi_cohort, "MFI", "DENV3", "skyblue")
plot_mfi(young, dengns1_4_mfi_cohort, "MFI", "DENV4", "skyblue")

# Plot ∆MFI
plot_mfi(young, change_dengns1_1, "∆MFI", "Change in DENV1", "skyblue")
plot_mfi(young, change_dengns1_2, "∆MFI", "Change in DENV2", "skyblue")
plot_mfi(young, change_dengns1_3, "∆MFI", "Change in DENV3", "skyblue")
plot_mfi(young, change_dengns1_4, "∆MFI", "Change in DENV4", "skyblue")

