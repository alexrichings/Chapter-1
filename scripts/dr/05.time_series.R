
# 05. Code to produce a time series of case data 
------------------------------------------------------------
  
# 1. Install and load packages ----

library(readxl)
library(dplyr)
library(tidyr)
library(ISOweek)
library(ggplot2)
library(here)
library(lubridate)
library(tidyverse)

------------------------------------------------------------
  
# 2. Load data ----

# data from end of 2017 until 2023 
data <- read_xlsx(here("data", "dengue-digepi_10-08-23.xlsx"))

data_2013 <- read_csv(here("data", "2013 etv.csv"))
data_2014 <- read_csv(here("data", "etv 2014.csv"))
data_2015 <- read_csv(here("data", "2015 etc.csv"))
data_2016 <- read_csv(here("data", "2016 etv.csv"))
data_2017 <- read_csv(here("data", "2017 etv.csv"))
data_2018 <- read_csv(here("data", "etv 2018.csv"))
data_2019 <- read_csv(here("data", "etv 2019.csv"))
data_2020 <- read_csv(here("data", "etv 2020.csv"))
data_2021 <- read_csv(here("data", "etv 2021.csv"))

------------------------------------------------------------
  
# 3. Process data ----

# load the data for 2013 - 2017 and save together 
df_13_17 <- list(
  data_2013,
  data_2014,
  data_2015,
  data_2016,
  data_2017
)

df_13_17 <- bind_rows(df_13_17)

# rename temporal variables for each period 

# 2013 - 2017 
df_13_17 <- df_13_17 %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# 2018 - 2023 
df_18_23 <- data %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# align cases by day 

# 2013 - 2017 
daily_cases_13_17 <- df_13_17 %>%
  mutate(date = dmy(`Fecha inicio síntomas`)) %>%
  filter(!is.na(date)) %>%
  count(date, name = "cases") %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day"),
    fill = list(cases = 0)
  )

# 2018 - 2023 
daily_cases_18_23 <- df_18_23 %>%
  mutate(date = as.Date(`Fecha inicio síntomas`)) %>%
  filter(!is.na(date)) %>%
  count(date, name = "cases") %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day"),
    fill = list(cases = 0)
  )

# align cases by month + province 

df_13_17 <- df_13_17 %>%
  mutate(date = dmy(`Fecha inicio síntomas`),
         month_date = floor_date(date, "month"))

df_18_23 <- df_18_23 %>%
  mutate(date = as.Date(`Fecha inicio síntomas`),
         month_date = floor_date(date, "month"))

# 2013 - 2017 
monthly_cases_13_17 <- df_13_17 %>%
  filter(!is.na(month_date)) %>%
  count(month_date, Provincia, name = "cases") %>%
  group_by(Provincia) %>%
  complete(
    month_date = seq(min(month_date), max(month_date), by = "month"),
    fill = list(cases = 0)
  ) %>%
  ungroup()

# 2018 - 2023 
monthly_cases_18_23 <- df_18_23 %>%
  filter(!is.na(month_date)) %>%
  count(month_date, Provincia, name = "cases") %>%
  group_by(Provincia) %>%
  complete(
    month_date = seq(min(month_date), max(month_date), by = "month"),
    fill = list(cases = 0)
  ) %>%
  ungroup()



# join the two data sets 

# add a cut-off to separate them 
cutoff <- as.Date("2017-12-25")

df_pre_cutoff <- daily_cases_13_17 %>% filter(date <= cutoff)
df_post_cutoff <- daily_cases_18_23 %>% filter(date > cutoff)

daily_cases_combined <- bind_rows(df_pre_cutoff, df_post_cutoff) %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day"),
    fill = list(cases = 0)
  )

saveRDS(daily_cases_combined, here("data", "daily_cases_combined.rds"))

# join monthly data 
monthly_cases_prov <- bind_rows(
  monthly_cases_13_17,
  monthly_cases_18_23
) %>%
  arrange(month_date, Provincia)

saveRDS(monthly_cases_prov, here("data", "monthly_cases_prov.rds"))


------------------------------------------------------------
  
# 4. Plot ----

# daily cases 
ggplot(daily_cases_combined, aes(x = date, y = cases)) +
  geom_line(colour = "#2E86C1", linewidth = 0.1, alpha = 0.7) +
  labs(
    x = "Date of symptom onset",
    y = "Daily Dengue Cases"
  ) +
  scale_x_date(
    limits = as.Date(c("2013-01-01", "2023-12-31")),
    breaks = seq(as.Date("2013-01-01"), as.Date("2023-12-31"), by = "year"),
    date_labels = "%Y"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(here("outputs", "dengue_dr_case_series2.png"), width = 12, height = 6, dpi = 300)




