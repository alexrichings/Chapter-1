
# Chapter 1 

--------------------------------------------------------------------------------

# Script: 01-setup_inspect ----  

# Purpose: 
# Import raw serological data, clean variables, visual inspection of data

# Output: 
# Paired serological data file: df_wide 

--------------------------------------------------------------------------------

# 1. Packages ----

library("readxl")
library("here")
library("ggplot2")
library("dplyr")

-------------------------------------------------------------------------------

# 2. Import raw data ----

# df_0: national serological survey; df_1: follow up cohort 
df_0 <- read_excel(here("data", "Copy of DR_serosurvey_MASTER_with_MBA_UPDATED.xlsx"))
df_1 <- read_excel(here("data", "Copy of FINAL_COHORT_PLUS_MBA.xlsx"))
                   

-------------------------------------------------------------------------------

# 3. Generate paired data set & rename/mutate variables ---- 

# round 1: x_sero 
# round 2: x_cohort 

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

-------------------------------------------------------------------------------
  
# 4. Save data ---- 

# saveRDS(df_wide, here("data", "df_wide.rds"))

-------------------------------------------------------------------------------

# 5. Visual data inspection  ---- 

# Changes in MFI 
# No major skew in any distribution - transformation may be unwarranted 

# DENV1 
ggplot(data = df_wide, mapping = aes(change_dengns1_1/1000)) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV1") + 
  theme_minimal() + 
  xlim(c(-15,15))

# DENV2
ggplot(data = df_wide, mapping = aes(change_dengns1_2/1000)) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV2") + 
  theme_minimal() + 
  xlim(c(-15,15))

# DENV3
ggplot(data = df_wide, mapping = aes((change_dengns1_3/1000))) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV3") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-15, 15))

# DENV4
ggplot(data = df_wide, mapping = aes((change_dengns1_4/1000))) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI for DENV4",
       y = "Density",
       title = "Change in DENV4 MFI") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-15, 15))




