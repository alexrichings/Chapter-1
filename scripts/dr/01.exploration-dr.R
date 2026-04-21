
# 01.exploration ----  

# preliminary description and visualisation of data from the Dominican Republic 
-------------------------------------------------------------------------------
# 1. Install and load packages ----

library("readxl")
library("ggplot2")
library("summarytools")
library("dplyr")
library("lme4")
library("gtsummary")
library("serofoi")
library("janitor")
library("tidyr")
library("here")

-------------------------------------------------------------------------------
# 2. Load in the data ----

df_0 <- read_excel(here("data", "Copy of DR_serosurvey_MASTER_with_MBA_UPDATED.xlsx"))
df_1 <- read_excel(here("data", "Copy of FINAL_COHORT_PLUS_MBA.xlsx"))
                   
-------------------------------------------------------------------------------
# 3. Inspect the variables ----

length(unique(df_0$serosurvey_id))
length(unique(df_1$cohort_id))

-------------------------------------------------------------------------------
# 4. Basic descriptive analysis ----

# Seroprevalence 
summary(df_0$dengns1_1_seropos) # DENV1: 93% 
summary(df_0$dengns1_2_seropos) # DENV2: 92% 
summary(df_0$dengns1_3_seropos) # DENV3: 95% 
summary(df_0$dengns1_4_seropos) # DENV4: 91% 
summary(df_0$zika_ns1_seropos) # ZIKV: 75% 
summary(df_0$chik_e1_seropos) # DENV1: 74% 

-------------------------------------------------------------------------------
# 6. Accounting for pairing structure ---- 

df_wide <- df_0 %>%
  select(serosurvey_id, age, gender, setting, occupation, education, province, region, occupation_2, education2) %>%
  rename(age_sero = age, gender_sero = gender, setting_sero = setting, occupation_sero = occupation, education_sero = education, occupation2_sero = occupation_2, education2_sero = education2) %>%
  inner_join(
    df_1 %>% 
      select(
        serosurvey_id, cohort_id, date_sero, date_cohort, age, gender, occupation, date_diff, 
        dengns1_1_mfi_sero, dengns1_2_mfi_sero, dengns1_3_mfi_sero, dengns1_4_mfi_sero, zika_ns1_mfi_sero, chik_e1_mfi_sero, 
        dengns1_1_mfi_cohort, dengns1_2_mfi_cohort, dengns1_3_mfi_cohort, dengns1_4_mfi_cohort, zika_ns1_mfi_cohort, chik_e1_mfi_cohort,
        dengns1_1_seropos_sero, dengns1_2_seropos_sero, dengns1_3_seropos_sero, dengns1_4_seropos_sero, zika_ns1_seropos_sero, chik_e1_seropos_sero,
        dengns1_1_seropos_cohort, dengns1_2_seropos_cohort, dengns1_3_seropos_cohort, dengns1_4_seropos_cohort, zika_ns1_seropos_cohort, chik_e1_seropos_cohort) %>%
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
# 7. Explore pairing ---- 

# look at changes in MFI 

# deng1
ggplot(data = df_wide, mapping = aes(change_dengns1_1/1000)) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV1") + 
  theme_minimal() + 
  xlim(c(-15,15))

# deng2
ggplot(data = df_wide, mapping = aes(change_dengns1_2/1000)) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV2") + 
  theme_minimal() + 
  xlim(c(-15,15))

# deng3
ggplot(data = df_wide, mapping = aes((change_dengns1_3/1000))) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI",
       y = "Density",
       title = "DENV3") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-15, 15))

# deng4
ggplot(data = df_wide, mapping = aes((change_dengns1_4/1000))) + 
  geom_density(stat = "density", alpha = 0.7, fill = "skyblue") + 
  labs(x = "∆MFI for DENV4",
       y = "Density",
       title = "Change in DENV4 MFI") + 
  theme_minimal() + 
  coord_cartesian(xlim = c(-15, 15))

# plot all on the same axis

ggplot(data = df_wide) + 
  geom_density(aes(x = (change_dengns1_1/1000), fill = "DENV1"), alpha = 0.5) + 
  geom_density(aes(x = (change_dengns1_2/1000), fill = "DENV2"), alpha = 0.5) + 
  geom_density(aes(x = (change_dengns1_3/1000), fill = "DENV3"), alpha = 0.5) + 
  geom_density(aes(x = (change_dengns1_4/1000), fill = "DENV4"), alpha = 0.5) + 
  labs(x = "∆MFI", 
       y = "Density",
       title = "Change in DENV MFI values in DR",
       fill = "Serotype") + 
  theme_minimal() +
  coord_cartesian(xlim = c(-15, 15))

-------------------------------------------

# 8. Save data ----

# saveRDS(df_wide, here("data", "df_wide.rds"))


