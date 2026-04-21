
# Regression analyses ----
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

------------------------------------------------------------------------
  
# 2. Data ---- 

inputs <- readRDS(here("data", "inputs_fj.rds"))

------------------------------------------------------------------------
 
# 3. Descriptive Analyses ----

# Outcome is continuous and bounded by [0,1]

table(inputs$denv3_prob, useNA = "ifany") # 3 NA values present 
table(inputs$denv1_prob, useNA = "ifany") # 12 NA values present 

# subset to seronegatives 
inputs_og <- inputs %>% filter(inputs$ELISA1<=9)

# tabulate the characteristics 
dim(inputs_og)

# Demographic characteristics 
table(inputs_og$AGE_U_20, useNA = "ifany")
table(inputs_og$SEX, useNA = "ifany")
table(inputs_og$ETHNIC, useNA = "ifany")

# Environmental factors present 
table(inputs_og$I_MOS, useNA = "ifany")
table(inputs_og$I_TIR, useNA = "ifany")
table(inputs_og$I_WAT, useNA = "ifany")
table(inputs_og$I_AC, useNA = "ifany")
table(inputs_og$I_BLK, useNA = "ifany")

# Location 

# Modify variables 
inputs <- inputs %>% mutate(GEOG2 = if_else(GEOG %in% c(1, 2), 1L, 0L)) %>% relocate(GEOG2, .after = GEOG)
inputs_og <- inputs_og %>% mutate(GEOG2 = if_else(GEOG %in% c(1, 2), 1L, 0L)) %>% relocate(GEOG2, .after = GEOG)

table(inputs_og$GEOG, useNA = "ifany")
table(inputs_og$GEOG2, useNA = "ifany")

# Health seeking behaviour 
table(inputs_og$FEVER_2YR, useNA = "ifany")
table(inputs_og$DOC_2YR, useNA = "ifany")
table(inputs_og$HH_D, useNA = "ifany")

----------------------------------------------------------------------------------------
  
# 4. Logistic regression ----

# use a p = 0.5 cut off for the MIA data 
inputs <- inputs %>%
  mutate(MIA_status_denv3 = if_else(denv3_prob >= 0.5, 1L, 0L))

inputs <- inputs %>%
  mutate(MIA_status_denv1 = if_else(denv1_prob >= 0.5, 1L, 0L))

inputs_og <- inputs_og %>%
  mutate(MIA_status_denv3 = if_else(denv3_prob >= 0.5, 1L, 0L))

inputs_og <- inputs_og %>%
  mutate(MIA_status_denv1 = if_else(denv1_prob >= 0.5, 1L, 0L))


# model type A: n = 97 
# logistic regression with the reduced number of participants 

# DENV3
modelA3.1 <- glm(MIA_status_denv3 ~ AGE_U_20 , data = inputs_og,family = "binomial")
modelA3.2 <- glm(MIA_status_denv3 ~ SEX , data = inputs_og,family = "binomial")
modelA3.3 <- glm(MIA_status_denv3 ~ ETHNIC , data = inputs_og,family = "binomial")
modelA3.4 <- glm(MIA_status_denv3 ~ I_MOS , data = inputs_og,family = "binomial")
modelA3.5 <- glm(MIA_status_denv3 ~ I_TIR , data = inputs_og,family = "binomial")
modelA3.6 <- glm(MIA_status_denv3 ~ I_WAT , data = inputs_og,family = "binomial")
modelA3.7 <- glm(MIA_status_denv3 ~ I_AC , data = inputs_og,family = "binomial")
modelA3.8 <- glm(MIA_status_denv3 ~ I_BLK , data = inputs_og,family = "binomial")
modelA3.9 <- glm(MIA_status_denv3 ~ GEOG2 , data = inputs_og,family = "binomial")
modelA3.10 <- glm(MIA_status_denv3 ~ FEVER_2YR , data = inputs_og,family = "binomial")
modelA3.11 <- glm(MIA_status_denv3 ~ DOC_2YR , data = inputs_og,family = "binomial")
modelA3.12 <- glm(MIA_status_denv3 ~ HH_D , data = inputs_og,family = "binomial")

model.list=list(modelA3.1, modelA3.2, modelA3.3, modelA3.4, modelA3.5, modelA3.6, modelA3.7, modelA3.8, modelA3.9, modelA3.10, modelA3.11, modelA3.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng3_n97.csv"))
deng3_n97_lr <- read.csv(here("outputs", "lr_deng3_n97.csv"))

# DENV1
modelA1.1 <- glm(MIA_status_denv1 ~ AGE_U_20 , data = inputs_og,family = "binomial")
modelA1.2 <- glm(MIA_status_denv1 ~ SEX , data = inputs_og,family = "binomial")
modelA1.3 <- glm(MIA_status_denv1 ~ ETHNIC , data = inputs_og,family = "binomial")
modelA1.4 <- glm(MIA_status_denv1 ~ I_MOS , data = inputs_og,family = "binomial")
modelA1.5 <- glm(MIA_status_denv1 ~ I_TIR , data = inputs_og,family = "binomial")
modelA1.6 <- glm(MIA_status_denv1 ~ I_WAT , data = inputs_og,family = "binomial")
modelA1.7 <- glm(MIA_status_denv1 ~ I_AC , data = inputs_og,family = "binomial")
modelA1.8 <- glm(MIA_status_denv1 ~ I_BLK , data = inputs_og,family = "binomial")
modelA1.9 <- glm(MIA_status_denv1 ~ GEOG2 , data = inputs_og,family = "binomial")
modelA1.10 <- glm(MIA_status_denv1 ~ FEVER_2YR , data = inputs_og,family = "binomial")
modelA1.11 <- glm(MIA_status_denv1 ~ DOC_2YR , data = inputs_og,family = "binomial")
modelA1.12 <- glm(MIA_status_denv1 ~ HH_D , data = inputs_og,family = "binomial")

model.list=list(modelA1.1, modelA1.2, modelA1.3, modelA1.4, modelA1.5, modelA1.6, modelA1.7, modelA1.8, modelA1.9, modelA1.10, modelA1.11, modelA1.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng1_n97.csv"))
deng1_n97_lr <- read.csv(here("outputs", "lr_deng1_n97.csv"))



# model type B: n = 260 (total population)
# logistic regression with the full number of participants 

table(inputs$MIA_status_denv3, useNA = "ifany")

# n = 260

# DENV3
modelB3.1 <- glm(MIA_status_denv3 ~ AGE_U_20 , data = inputs,family = "binomial")
modelB3.2 <- glm(MIA_status_denv3 ~ SEX , data = inputs,family = "binomial")
modelB3.3 <- glm(MIA_status_denv3 ~ ETHNIC , data = inputs,family = "binomial")
modelB3.4 <- glm(MIA_status_denv3 ~ I_MOS , data = inputs,family = "binomial")
modelB3.5 <- glm(MIA_status_denv3 ~ I_TIR , data = inputs,family = "binomial")
modelB3.6 <- glm(MIA_status_denv3 ~ I_WAT , data = inputs,family = "binomial")
modelB3.7 <- glm(MIA_status_denv3 ~ I_AC , data = inputs,family = "binomial")
modelB3.8 <- glm(MIA_status_denv3 ~ I_BLK , data = inputs,family = "binomial")
modelB3.9 <- glm(MIA_status_denv3 ~ GEOG2 , data = inputs,family = "binomial")
modelB3.10 <- glm(MIA_status_denv3 ~ FEVER_2YR , data = inputs,family = "binomial")
modelB3.11 <- glm(MIA_status_denv3 ~ DOC_2YR , data = inputs,family = "binomial")
modelB3.12 <- glm(MIA_status_denv3 ~ HH_D , data = inputs,family = "binomial")

model.list=list(modelB3.1, modelB3.2, modelB3.3, modelB3.4, modelB3.5, modelB3.6, modelB3.7, modelB3.8, modelB3.9, modelB3.10, modelB3.11, modelB3.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng3_n260.csv"))
deng3_n260_lr <- read.csv(here("outputs", "lr_deng3_n260.csv"))


# DENV1
modelB1.1 <- glm(MIA_status_denv1 ~ AGE_U_20 , data = inputs,family = "binomial")
modelB1.2 <- glm(MIA_status_denv1 ~ SEX , data = inputs,family = "binomial")
modelB1.3 <- glm(MIA_status_denv1 ~ ETHNIC , data = inputs,family = "binomial")
modelB1.4 <- glm(MIA_status_denv1 ~ I_MOS , data = inputs,family = "binomial")
modelB1.5 <- glm(MIA_status_denv1 ~ I_TIR , data = inputs,family = "binomial")
modelB1.6 <- glm(MIA_status_denv1 ~ I_WAT , data = inputs,family = "binomial")
modelB1.7 <- glm(MIA_status_denv1 ~ I_AC , data = inputs,family = "binomial")
modelB1.8 <- glm(MIA_status_denv1 ~ I_BLK , data = inputs,family = "binomial")
modelB1.9 <- glm(MIA_status_denv1 ~ GEOG2 , data = inputs,family = "binomial")
modelB1.10 <- glm(MIA_status_denv1 ~ FEVER_2YR , data = inputs,family = "binomial")
modelB1.11 <- glm(MIA_status_denv1 ~ DOC_2YR , data = inputs,family = "binomial")
modelB1.12 <- glm(MIA_status_denv1 ~ HH_D , data = inputs,family = "binomial")

model.list=list(modelB1.1, modelB1.2, modelB1.3, modelB1.4, modelB1.5, modelB1.6, modelB1.7, modelB1.8, modelB1.9, modelB1.10, modelB1.11, modelB1.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng1_n260.csv"))
deng1_n260_lr <- read.csv(here("outputs", "lr_deng1_n260.csv"))

----------------------------------------------------------------------------------------

# 5. Fractional logistic regression ----

# model type C: n = 97 
# FLR with reduced number of participants 

# denv3
modelC3.1 <- glm(denv3_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.2 <- glm(denv3_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.3 <- glm(denv3_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.4 <- glm(denv3_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.5 <- glm(denv3_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.6 <- glm(denv3_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.7 <- glm(denv3_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.8 <- glm(denv3_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.9 <- glm(denv3_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.10 <- glm(denv3_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.11 <- glm(denv3_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC3.12 <- glm(denv3_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs_og)

model.list=list(modelC3.1, modelC3.2, modelC3.3, modelC3.4, modelC3.5, modelC3.6, modelC3.7, modelC3.8, modelC3.9, modelC3.10, modelC3.11, modelC3.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng3_n97.csv"))
deng3_n97_flr <- read.csv(here("outputs", "flr_deng3_n97.csv"))


# denv1
modelC1.1 <- glm(denv1_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.2 <- glm(denv1_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.3 <- glm(denv1_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.4 <- glm(denv1_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.5 <- glm(denv1_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.6 <- glm(denv1_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.7 <- glm(denv1_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.8 <- glm(denv1_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.9 <- glm(denv1_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.10 <- glm(denv1_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.11 <- glm(denv1_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC1.12 <- glm(denv1_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs_og)

model.list=list(modelC1.1, modelC1.2, modelC1.3, modelC1.4, modelC1.5, modelC1.6, modelC1.7, modelC1.8, modelC1.9, modelC1.10, modelC1.11, modelC1.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng1_n97.csv"))
deng1_n97_flr <- read.csv(here("outputs", "flr_deng1_n97.csv"))


# model type D: n = total 
# FLR with total number of participants (n = 260)

# denv3
modelD3.1 <- glm(denv3_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs)
modelD3.2 <- glm(denv3_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs)
modelD3.3 <- glm(denv3_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs)
modelD3.4 <- glm(denv3_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs)
modelD3.5 <- glm(denv3_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs)
modelD3.6 <- glm(denv3_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs)
modelD3.7 <- glm(denv3_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs)
modelD3.8 <- glm(denv3_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs)
modelD3.9 <- glm(denv3_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs)
modelD3.10 <- glm(denv3_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD3.11 <- glm(denv3_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD3.12 <- glm(denv3_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs)

model.list=list(modelD3.1, modelD3.2, modelD3.3, modelD3.4, modelD3.5, modelD3.6, modelD3.7, modelD3.8, modelD3.9, modelD3.10, modelD3.11, modelD3.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng3_n260.csv"))
deng3_n260_flr <- read.csv(here("outputs", "flr_deng3_n260.csv"))


# denv1
modelD1.1 <- glm(denv1_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs)
modelD1.2 <- glm(denv1_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs)
modelD1.3 <- glm(denv1_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs)
modelD1.4 <- glm(denv1_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs)
modelD1.5 <- glm(denv1_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs)
modelD1.6 <- glm(denv1_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs)
modelD1.7 <- glm(denv1_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs)
modelD1.8 <- glm(denv1_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs)
modelD1.9 <- glm(denv1_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs)
modelD1.10 <- glm(denv1_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD1.11 <- glm(denv1_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD1.12 <- glm(denv1_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs)

model.list=list(modelD1.1, modelD1.2, modelD1.3, modelD1.4, modelD1.5, modelD1.6, modelD1.7, modelD1.8, modelD1.9, modelD1.10, modelD1.11, modelD1.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng1_n260.csv")
deng1_n260_flr <- read.csv(here("outputs", "flr_deng1_n260.csv"))

---------------------------------------------------------------------------------------

# 6. Collapsing DENV serotypes  ----

# inspect the data 
table(inputs$denv1_prob, useNA = "ifany")
table(inputs$denv3_prob, useNA = "ifany")

# create a df that is the maximum probability of the above columns 
inputs <- inputs %>% mutate(max_prob = pmax(denv1_prob, denv3_prob, na.rm = TRUE))

# create the binary form of the variable 
inputs <- inputs %>% mutate(MIA_status_denv = if_else(max_prob >= 0.5, 1L, 0L))

# assign to the old data set 
inputs_og <- inputs %>% filter(inputs$ELISA1<=9)
dim(inputs_og)

---------------------------------------------------------------------------------------
  
# 7. Logistic regression (without serotype) ----


# model type A: n = 97 ----
# LR for reduced participant data set (only those seronegative)

modelA.1 <- glm(MIA_status_denv ~ AGE_U_20 , data = inputs_og,family = "binomial")
modelA.2 <- glm(MIA_status_denv ~ SEX , data = inputs_og,family = "binomial")
modelA.3 <- glm(MIA_status_denv ~ ETHNIC , data = inputs_og,family = "binomial")
modelA.4 <- glm(MIA_status_denv ~ I_MOS , data = inputs_og,family = "binomial")
modelA.5 <- glm(MIA_status_denv ~ I_TIR , data = inputs_og,family = "binomial")
modelA.6 <- glm(MIA_status_denv ~ I_WAT , data = inputs_og,family = "binomial")
modelA.7 <- glm(MIA_status_denv ~ I_AC , data = inputs_og,family = "binomial")
modelA.8 <- glm(MIA_status_denv ~ I_BLK , data = inputs_og,family = "binomial")
modelA.9 <- glm(MIA_status_denv ~ GEOG2 , data = inputs_og,family = "binomial")
modelA.10 <- glm(MIA_status_denv ~ FEVER_2YR , data = inputs_og,family = "binomial")
modelA.11 <- glm(MIA_status_denv ~ DOC_2YR , data = inputs_og,family = "binomial")
modelA.12 <- glm(MIA_status_denv ~ HH_D , data = inputs_og,family = "binomial")

model.list=list(modelA.1, modelA.2, modelA.3, modelA.4, modelA.5, modelA.6, modelA.7, modelA.8, modelA.9, modelA.10, modelA.11, modelA.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng_n97.csv"))
lr_deng_n97 <- read.csv(here("outputs", "lr_deng_n97.csv"))

# model type B: n = total ----
# LR but for all participants 

modelB.1 <- glm(MIA_status_denv ~ AGE_U_20 , data = inputs,family = "binomial")
modelB.2 <- glm(MIA_status_denv ~ SEX , data = inputs,family = "binomial")
modelB.3 <- glm(MIA_status_denv ~ ETHNIC , data = inputs,family = "binomial")
modelB.4 <- glm(MIA_status_denv ~ I_MOS , data = inputs,family = "binomial")
modelB.5 <- glm(MIA_status_denv ~ I_TIR , data = inputs,family = "binomial")
modelB.6 <- glm(MIA_status_denv ~ I_WAT , data = inputs,family = "binomial")
modelB.7 <- glm(MIA_status_denv ~ I_AC , data = inputs,family = "binomial")
modelB.8 <- glm(MIA_status_denv ~ I_BLK , data = inputs,family = "binomial")
modelB.9 <- glm(MIA_status_denv ~ GEOG2 , data = inputs,family = "binomial")
modelB.10 <- glm(MIA_status_denv ~ FEVER_2YR , data = inputs,family = "binomial")
modelB.11 <- glm(MIA_status_denv ~ DOC_2YR , data = inputs,family = "binomial")
modelB.12 <- glm(MIA_status_denv ~ HH_D , data = inputs,family = "binomial")

model.list=list(modelB.1, modelB.2, modelB.3, modelB.4, modelB.5, modelB.6, modelB.7, modelB.8, modelB.9, modelB.10, modelB.11, modelB.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "lr_deng_n260.csv"))
lr_deng_n260 <- read.csv(here("outputs", "lr_deng_n260.csv"))

---------------------------------------------------------------------------------------
  
# 8. Fractional logistic regression (without serotype) ----

# model type C: n = 97 
# FLR for the reduced data set (seronegatives only)

modelC.1 <- glm(max_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.2 <- glm(max_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.3 <- glm(max_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.4 <- glm(max_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.5 <- glm(max_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.6 <- glm(max_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.7 <- glm(max_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.8 <- glm(max_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.9 <- glm(max_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.10 <- glm(max_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.11 <- glm(max_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs_og)
modelC.12 <- glm(max_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs_og)

model.list=list(modelC.1, modelC.2, modelC.3, modelC.4, modelC.5, modelC.6, modelC.7, modelC.8, modelC.9, modelC.10, modelC.11, modelC.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

data.tally = inputs_og[,names.model]

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng_n97.csv"))
flr_deng_n97 <- read.csv(here("outputs", "flr_deng_n97.csv"))

# model type D: n = total 
# FLR for full data set 

modelD.1 <- glm(max_prob ~ AGE_U_20, family = quasibinomial(link = "logit"), data = inputs)
modelD.2 <- glm(max_prob ~ SEX, family = quasibinomial(link = "logit"), data = inputs)
modelD.3 <- glm(max_prob ~ ETHNIC, family = quasibinomial(link = "logit"), data = inputs)
modelD.4 <- glm(max_prob ~ I_MOS, family = quasibinomial(link = "logit"), data = inputs)
modelD.5 <- glm(max_prob ~ I_TIR, family = quasibinomial(link = "logit"), data = inputs)
modelD.6 <- glm(max_prob ~ I_WAT, family = quasibinomial(link = "logit"), data = inputs)
modelD.7 <- glm(max_prob ~ I_AC, family = quasibinomial(link = "logit"), data = inputs)
modelD.8 <- glm(max_prob ~ I_BLK, family = quasibinomial(link = "logit"), data = inputs)
modelD.9 <- glm(max_prob ~ GEOG2, family = quasibinomial(link = "logit"), data = inputs)
modelD.10 <- glm(max_prob ~ FEVER_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD.11 <- glm(max_prob ~ DOC_2YR, family = quasibinomial(link = "logit"), data = inputs)
modelD.12 <- glm(max_prob ~ HH_D, family = quasibinomial(link = "logit"), data = inputs)

model.list=list(modelD.1, modelD.2, modelD.3, modelD.4, modelD.5, modelD.6, modelD.7, modelD.8, modelD.9, modelD.10, modelD.11, modelD.12)
names.model=c("AGE_U_20","SEX","ETHNIC","I_MOS","I_TIR","I_WAT", "I_AC","I_BLK","GEOG","FEVER_2YR","DOC_2YR","HH_D")

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

store.oddsTable
write.csv(store.oddsTable, here("outputs", "flr_deng_n260.csv"))
flr_deng_n260 <- read.csv(here("outputs", "flr_deng_n260.csv"))







