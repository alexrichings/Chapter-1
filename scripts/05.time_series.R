
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

# have data from different years 
# need to go back and separate each year in the data sets 

# add to a list 
# overlap between year ends 

df_list <- list(
  data_2013,
  data_2014,
  data_2015,
  data_2016,
  data_2017,
  data_2018,
  data_2019,
  data_2020,
  data_2021
)

# df for all data from 2013-2021 
df <- bind_rows(df_list)

# rename variables for 2013-2021
df_data <- df %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# rename variables for 2017-2023
df_data_mix <- data %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# align cases by day 
daily_cases_13_21 <- df_data %>%
  mutate(date = dmy(`Fecha inicio síntomas`)) %>%
  filter(!is.na(date)) %>%
  count(date, name = "cases") %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day"),
    fill = list(cases = 0)
  )

daily_cases_17_23 <- data %>%
  mutate(date = as.Date(`Fecha inicio síntomas`)) %>%
  filter(!is.na(date)) %>%
  count(date, name = "cases") %>%
  arrange(date) %>%
  complete(
    date = seq(min(date), max(date), by = "day"),
    fill = list(cases = 0)
  )


# join the two data sets 

# daily_cases up to end of 2021
df_pre22 <- daily_cases_13_21 %>%
  filter(date <= as.Date("2021-12-31"))

# daily_cases_data from 2022 onwards
df_post22 <- daily_cases_17_23 %>%
  filter(date >= as.Date("2021-06-01"))

# bind
daily_cases_combined <- bind_rows(df_pre22, df_post22) %>% arrange(date)

------------------------------------------------------------
  
# 4. Plot ----

# daily cases 
p <- ggplot(daily_cases_combined, aes(x = date, y = cases)) +
  geom_line(colour = "#2E86C1", linewidth = 0.1, alpha = 0.7) +
  labs(
    x = "Date of symptom onset",
    y = "Cases"
  ) +
  scale_x_date(
    limits = as.Date(c("2013-01-01", "2023-12-31")),
    breaks = seq(as.Date("2013-01-01"), as.Date("2023-12-31"), by = "year"),
    date_labels = "%Y"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p
