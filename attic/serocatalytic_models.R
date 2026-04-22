
# 10. Serocatalytic models ----

# load in the libaries ----
library(readxl)
library(dplyr)
library(serofoi)
library(tidyr)

# load in the data ----
df_wide <- readRDS("Data/DR_surveillance/EN/df_wide")
df_0 <- read_excel('Data/DR_surveillance/EN/Copy of DR_serosurvey_MASTER_with_MBA_UPDATED.xlsx')
df_1 <- read_excel('Data/DR_surveillance/EN/Copy of FINAL_COHORT_PLUS_MBA.xlsx')

ggplot(data = df_wide, aes(x = age_sero)) +
  geom_histogram(binwidth = 5, fill = "lightgrey", color = "black") + 
  coord_cartesian(xlim = c(0, 100)) +
  theme_minimal()

# subset the data by pathogen ----
make.df <- function(data, pathogen_serostatus, year) {
  path_df <- data %>%
    group_by(age) %>%
    summarise(
      n_sample = sum(!is.na({{pathogen_serostatus}})),
      n_seropositive = sum({{pathogen_serostatus}}, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      survey_year = year,
      age_min = age,
      age_max = age
    ) %>%
    select(survey_year, n_sample, n_seropositive, age_min, age_max)
  
  return(path_df)
}

deng1.df <- make.df(df_0, dengns1_1_seropos, 2021)
deng2.df <- make.df(df_0, dengns1_2_seropos, 2021)
deng3.df <- make.df(df_0, dengns1_3_seropos, 2021)
deng4.df <- make.df(df_0, dengns1_4_seropos, 2021)
zika.df <- make.df(df_0, zika_ns1_seropos, 2021)
chik.df <- make.df(df_0, chik_e1_seropos, 2021)

# plot serosurvey ----
deng1.plot <- plot_serosurvey(deng1.df, bin_serosurvey = TRUE, size_text = 15)
deng2.plot <- plot_serosurvey(deng2.df, bin_serosurvey = TRUE, size_text = 15)
deng3.plot <- plot_serosurvey(deng3.df, bin_serosurvey = TRUE, size_text = 15)
deng4.plot <- plot_serosurvey(deng4.df, bin_serosurvey = TRUE, size_text = 15)
zika.plot <- plot_serosurvey(zika.df, bin_serosurvey = TRUE, size_text = 15)
chik.plot <- plot_serosurvey(chik.df, bin_serosurvey = TRUE, size_text = 15)

# constant FOI 
deng1_constant <- fit_seromodel(serosurvey = deng1.df, model_type = "constant")
deng2_constant <- fit_seromodel(serosurvey = deng2.df, model_type = "constant")
deng3_constant <- fit_seromodel(serosurvey = deng3.df, model_type = "constant")
deng4_constant <- fit_seromodel(serosurvey = deng4.df, model_type = "constant")
zika_constant <- fit_seromodel(serosurvey = zika.df, model_type = "constant")
chik_constant <- fit_seromodel(serosurvey = chik.df, model_type = "constant")

deng1_constant.plot <- plot_seromodel(deng1_constant, serosurvey = deng1.df, bin_serosurvey = TRUE, size_text = 6)
deng2_constant.plot <- plot_seromodel(deng2_constant, serosurvey = deng2.df, bin_serosurvey = TRUE, size_text = 6)
deng3_constant.plot <- plot_seromodel(deng3_constant, serosurvey = deng3.df, bin_serosurvey = TRUE, size_text = 6)
deng4_constant.plot <- plot_seromodel(deng4_constant, serosurvey = deng4.df, bin_serosurvey = TRUE, size_text = 6)
zika_constant.plot <- plot_seromodel(zika_constant, serosurvey = zika.df, bin_serosurvey = TRUE, size_text = 6)
chik_constant.plot <- plot_seromodel(chik_constant, serosurvey = chik.df, bin_serosurvey = TRUE, size_text = 6)

# time-varying FOI ----

# set an age-group index 
deng1_index <- get_foi_index(serosurvey = deng1.df, group_size = 5, model_type = "time")
deng2_index <- get_foi_index(serosurvey = deng2.df, group_size = 5, model_type = "time")
deng3_index <- get_foi_index(serosurvey = deng3.df, group_size = 5, model_type = "time")
deng4_index <- get_foi_index(serosurvey = deng4.df, group_size = 5, model_type = "time")
zika_index <- get_foi_index(serosurvey = zika.df, group_size = 5, model_type = "time")
chik_index <- get_foi_index(serosurvey = chik.df, group_size = 5, model_type = "time")

# set the FOI priors 
foi_prior <- sf_normal(mean = 0.077, sd = 0.05)

deng1_time <- fit_seromodel(serosurvey = deng1.df, model_type = "time", foi_prior = foi_prior, foi_index = deng1_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
deng2_time <- fit_seromodel(serosurvey = deng2.df, model_type = "time", foi_prior = foi_prior, foi_index = deng2_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
deng3_time <- fit_seromodel(serosurvey = deng3.df, model_type = "time", foi_prior = foi_prior, foi_index = deng3_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
deng4_time <- fit_seromodel(serosurvey = deng4.df, model_type = "time", foi_prior = foi_prior, foi_index = deng4_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
zika_time <- fit_seromodel(serosurvey = zika.df, model_type = "time", foi_prior = foi_prior, foi_index = zika_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))
chik_time <- fit_seromodel(serosurvey = chik.df, model_type = "time", foi_prior = foi_prior, foi_index = chik_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))

deng1_time.plot <- plot_seromodel(deng1_time, serosurvey = deng1.df, bin_serosurvey = TRUE, size_text = 6)
deng2_time.plot <- plot_seromodel(deng2_time, serosurvey = deng2.df, bin_serosurvey = TRUE, size_text = 6)
deng3_time.plot <- plot_seromodel(deng3_time, serosurvey = deng3.df, bin_serosurvey = TRUE, size_text = 6)
deng4_time.plot <- plot_seromodel(deng4_time, serosurvey = deng4.df, bin_serosurvey = TRUE, size_text = 6)
zika_time.plot <- plot_seromodel(zika_time, serosurvey = zika.df, bin_serosurvey = TRUE, size_text = 6)
chik_time.plot <- plot_seromodel(chik_time, serosurvey = chik.df, bin_serosurvey = TRUE, size_text = 6)

# make a universal dengue metric 
deng.df <- df_0 %>% 
  mutate(
    deng_stat = ifelse(dengns1_1_seropos == 1 | dengns1_2_seropos == 1 | dengns1_3_seropos == 1 | dengns1_4_seropos == 1, 1, 0)
  )

  # plot serosurvey 
deng_total.df <- make.df(deng.df, deng_stat, 2021)
deng.plot <- plot_serosurvey(deng_total.df, bin_serosurvey = TRUE, size_text = 15)

# fit the seromodel for the combined dengue data
deng_index <- get_foi_index(serosurvey = deng_total.df, group_size = 5, model_type = "time")
deng_time <- fit_seromodel(serosurvey = deng_total.df, model_type = "time", foi_prior = foi_prior, foi_index = deng_index,
                            is_log_foi = TRUE, iter = 4000, warmup = 2000, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15))

deng_time.plot <- plot_seromodel(deng_time, serosurvey = deng_total.df, bin_serosurvey = TRUE, size_text = 6)

------------------------------------------------------------------------------------------------------------------------------------------------
# surveillance data ---- 

surv <- read_excel("Data/DR_surveillance/dengue-digepi_10-08-23.xlsx")

# recode the temporal variables 
surv <- surv %>%
  rename(
    week = `Semana inicio síntomas`,
    month = `Mes inicio síntomas`,
    year = `Año inicio síntomas`
  )

# count the weekly cases 
weekly_cases <- surv %>%
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

# plot weekly cases 
ggplot(weekly_cases, aes(x = week_date, y = cases)) +
  geom_line() +
  labs(x = "Week", y = "Cases", title = "Weekly dengue cases in Dominican Republic (2018 - 2023)") + 
  theme_minimal() + 
  ylim(0, 2000) 

# need to determine reporting rate 

# reporting rate





summary(df_wide$age_cohort, useNA = "ifany")
sum(weekly_cases$cases, na.rm = TRUE)

head(weekly_cases$week_date)

# 63,461 cases 
# 2016 days 

# 0.175 foi 
