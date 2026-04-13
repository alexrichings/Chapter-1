
# 05. Code to produce a time series of case data 
------------------------------------------------------------
  
# 1. Install and load packages ----

library(readxl)
library(dplyr)
library(tidyr)
library(ISOweek)
library(GHRexplore)
library(ggplot2)
library(here)

------------------------------------------------------------
  
# 2. Load data ----

data <- read_xlsx(here("data", "dengue-digepi_10-08-23.xlsx"))

------------------------------------------------------------
  
# 3. Process data ----

# rename variables 
data <- data %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# align cases by week 
weekly_cases <- data %>%
  filter(!is.na(year), !is.na(week)) %>%
  mutate(
    year = as.integer(year),
    week = as.integer(week),
    week_date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1")
    )
  ) %>%
  count(week_date, name = "cases") %>%
  arrange(week_date) %>%
  complete(
    week_date = seq(min(week_date), max(week_date), by = "week"),
    fill = list(cases = 0)
  )

# align weekly cases by province 
weekly_cases_prov <- data %>%
  filter(!is.na(year), !is.na(week), !is.na(Provincia)) %>%
  mutate(
    year = as.integer(year),
    week = as.integer(week),
    week_date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1"))
  ) %>%
  count(Provincia, week_date, name = "cases") %>%
  arrange(Provincia, week_date) %>%
  group_by(Provincia) %>%
  complete(
    week_date = seq(min(week_date), max(week_date), by = "week"),
    fill = list(cases = 0)
  ) %>%
  ungroup()

------------------------------------------------------------
  
# 4. Plot ----

# weekly cases 
ggplot(weekly_cases, aes(x = week_date, y = cases)) +
  geom_line() +
  labs(x = "Week", y = "Cases", title = "Weekly dengue cases in Dominican Republic (2018 - 2023)") + 
  theme_minimal() + 
  ylim(0, 2000) 

# weekly cases by province
ggplot(weekly_cases_prov, aes(x = week_date, y = cases, colour = Provincia)) +
  geom_line() +
  labs(x = "Week", y = "Cases", title = "Weekly dengue cases by Provincia") +
  theme_minimal()

# plot time series with GHRmodel package
plot_timeseries(weekly_cases_prov, var = "cases", type = "counts", time = "week_date", area = "Provincia")

