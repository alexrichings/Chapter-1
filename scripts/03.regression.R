
# 03.regression 

# regression analyses using outputs from FMMs 

-------------------------------------------------------------------------------
  
# 1. Install and load packages ----

library("ggplot2")
library("dplyr")
library("here")

-------------------------------------------------------------------------------
  
# 2. Load data ----

inputs <- readRDS(here("data", "inputs.rds"))

-------------------------------------------------------------------------------
  
# 3. Data wrangling ----

# use p = 0.5 threshold to assign boosting status 
inputs <- inputs %>%
  mutate(deng1_boost = factor(if_else(deng1_prob < 0.5, "no", "yes"))) %>% 
  mutate(deng2_boost = factor(if_else(deng2_prob < 0.5, "no", "yes"))) %>% 
  mutate(deng3_boost = factor(if_else(deng3_prob < 0.5, "no", "yes"))) %>% 
  mutate(deng4_boost = factor(if_else(deng4_prob < 0.5, "no", "yes"))) 

summary(inputs$deng1_boost, useNA = "ifany")
summary(inputs$deng2_boost, useNA = "ifany")
summary(inputs$deng3_boost, useNA = "ifany")
summary(inputs$deng4_boost, useNA = "ifany")

# use any sign of boosting to assign a dengue boost 
inputs <- inputs %>% 
  mutate(deng_boost = factor(if_else(
    deng1_boost == "yes" | deng2_boost == "yes" | deng3_boost == "yes" | deng4_boost == "yes", "yes", "no")))

summary(inputs$deng_boost, useNA = "ifany")

-------------------------------------------------------------------------------

# 3. Descriptive Analyses ----

# Outcome is binary (yes/no)
table(inputs$deng_boost, useNA = "ifany") # 1 NA values present 

# Recategorise 

# age 
inputs <- inputs %>%
  mutate(age_group = cut(age_cohort, breaks = c(6, 20, 30, 40, 50, 60, 70, Inf), right = FALSE,
    labels = c("6–20", "20–30", "30–40", "40–50","50–60", "60–70", "70+")))

inputs <- inputs %>%
  mutate(age_u20 = cut(age_cohort, breaks = c(6, 20, Inf), right = FALSE,
                         labels = c("under 20", "over 20")))

# education 
inputs <- inputs %>%
  mutate(
    education_group = as.character(education_sero),
    education_group = na_if(education_group, "dont_know_refuse"),
    education_group = na_if(education_group, "NA"),
    education_group = factor(education_group)
  )

# Demographic characteristics 
table(inputs$age_group, useNA = "ifany")
table(inputs$age_u20, useNA = "ifany")
table(inputs$gender_cohort, useNA = "ifany")
table(inputs$setting_sero, useNA = "ifany")
table(inputs$occupation_sero, useNA = "ifany")
table(inputs$education_group, useNA = "ifany")
table(inputs$province, useNA = "ifany")
table(inputs$region, useNA = "ifany")

-------------------------------------------------------------------------------
  
# 4. Logistic regression ----

model1.1 <- glm(deng_boost ~ age_group, data = inputs, family = "binomial")
model1.2 <- glm(deng_boost ~ age_u20, data = inputs, family = "binomial")
model1.3 <- glm(deng_boost ~ gender_cohort, data = inputs, family = "binomial")
model1.4 <- glm(deng_boost ~ setting_sero, data = inputs, family = "binomial")
#model1.5 <- glm(deng_boost ~ occupation_sero, data = inputs, family = "binomial")
model1.6 <- glm(deng_boost ~ province, data = inputs, family = "binomial")
model1.7 <- glm(deng_boost ~ region, data = inputs, family = "binomial")

model.list=list(model1.1, model1.2, model1.3, model1.4, model1.6, model1.7)
names.model=c("age_group","age_u20","gender_cohort","setting_sero","province", "region")

data.tally = inputs[,names.model]

# calculate ORs for the different models 
store.oddsTable = NULL
for (jj in seq_along(model.list)) {
  modelT <- model.list[[jj]]
  
  mf <- model.frame(modelT)
  x  <- mf[[2]]
  
  tab <- table(x, useNA = "ifany")
  counts_str <- paste(names(tab), tab, sep="=", collapse="; ")
  
  estm.D <- exp(coef(modelT))
  conf.D <- exp(confint(modelT))
  pval   <- coef(summary(modelT))[, 4][2]
  
  out_row <- cbind(
    var = names.model[jj],
    N_used = nobs(modelT),
    counts = counts_str,
    OR_CI = paste0(round(estm.D[2], 2), " (", paste(round(conf.D[2,], 2), collapse="-"), ")"),
    p = round(pval, 2)
  )
  
  store.oddsTable <- rbind(store.oddsTable, out_row)
}

store.oddsTable <- NULL

for (jj in seq_along(model.list)) {
  
  modelT <- model.list[[jj]]
  varname <- names.model[jj]
  
  mf <- model.frame(modelT)
  x  <- mf[[2]]
  
  tab <- table(x, useNA = "ifany")
  counts_str <- paste(names(tab), tab, sep="=", collapse="; ")
  
  estm <- coef(modelT)
  conf <- confint(modelT)
  pval <- coef(summary(modelT))[, 4]
  
  # remove intercept
  estm <- estm[-1]
  conf <- conf[-1, , drop = FALSE]
  pval <- pval[-1]
  
  # loop through each level
  for (i in seq_along(estm)) {
    
    out_row <- data.frame(
      var = varname,
      level = names(estm)[i],
      N_used = nobs(modelT),
      counts = counts_str,
      OR_CI = paste0(
        round(exp(estm[i]), 2),
        " (",
        paste(round(exp(conf[i, ]), 2), collapse = "-"),
        ")"
      ),
      p = round(pval[i], 3),
      stringsAsFactors = FALSE
    )
    
    store.oddsTable <- rbind(store.oddsTable, out_row)
  }
}

store.oddsTable

write.csv(store.oddsTable, here("outputs", "lr_deng_fmm.csv"))
lr_deng_fmm <- read.csv(here("outputs", "lr_deng_fmm.csv"))
