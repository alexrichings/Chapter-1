
# Regression analyses ---- 

--------------------------------------------------------------------------------

library(lme4)  
  
# load in the data ----
df_wide <- readRDS("Data/DR_surveillance/EN/df_wide")

inputs <- df_wide
  
# 7. Logistic regression ----

# use a p = 0.5 cut off for the MIA data 
inputs <- inputs %>% mutate(MIA_status_denv1 = if_else(prob_deng1 >= 0.5, 1L, 0L))
inputs <- inputs %>% mutate(MIA_status_denv2 = if_else(prob_deng2 >= 0.5, 1L, 0L))
inputs <- inputs %>% mutate(MIA_status_denv3 = if_else(prob_deng3 >= 0.5, 1L, 0L))
inputs <- inputs %>% mutate(MIA_status_denv4 = if_else(prob_deng4 >= 0.5, 1L, 0L))
inputs <- inputs %>% mutate(MIA_status_zika = if_else(prob_zika >= 0.5, 1L, 0L))

# model type A: logistic regression ----

# DENV1
modelA1.1 <- glm(MIA_status_denv1 ~ age_cohort, data = inputs, family = "binomial")
modelA1.2 <- glm(MIA_status_denv1 ~ gender_cohort, data = inputs, family = "binomial")
modelA1.3 <- glm(MIA_status_denv1 ~ setting_sero, data = inputs, family = "binomial")
# modelA1.4 <- glm(MIA_status_denv1 ~ occupation2_sero, data = inputs, family = "binomial")
# modelA1.5 <- glm(MIA_status_denv1 ~ education2_sero, data = inputs, family = "binomial")
modelA1.6 <- glm(MIA_status_denv1 ~ province, data = inputs, family = "binomial")
modelA1.7 <- glm(MIA_status_denv1 ~ region, data = inputs, family = "binomial")

model.list=list(modelA1.1, modelA1.2, modelA1.3, modelA1.6, modelA1.7)
model.list=list(modelA1.1, modelA1.2, modelA1.3, modelA1.4, modelA1.5, modelA1.6, modelA1.7)

names.model=c("age_cohort","gender_cohort","education2_sero","province", "region")
names.model=c("age_cohort","gender_cohort","setting_sero","occupation2_sero","education2_sero","province", "region")

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
  
  coefs <- summary(modelT)$coefficients
  estm.D <- exp(coef(modelT))
  conf.D <- exp(confint.default(modelT))  # faster than confint()
  
  for (k in 2:nrow(coefs)) {
    
    out_row <- cbind(
      var = rownames(coefs)[k],
      N_used = nobs(modelT),
      OR_CI = paste0(
        round(estm.D[k], 2),
        " (",
        round(conf.D[k,1], 2), "-",
        round(conf.D[k,2], 2),
        ")"
      ),
      p = round(coefs[k,4], 3)
    )
    
    store.oddsTable <- rbind(store.oddsTable, out_row)
  }
}

store.oddsTable

write.csv(store.oddsTable,"Outputs/")


table(inputs$province, inputs$MIA_status_denv1, useNA = "always")

