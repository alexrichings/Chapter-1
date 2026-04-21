
# 03.regression 

# regression analyses using outputs from FMMs 

-------------------------------------------------------------------------------
  
# 1. Install and load packages ----

library("ggplot2")
library("dplyr")
library("here")
library("flextable")
library("officer")

-------------------------------------------------------------------------------
  
# 2. Load data ----

inputs <- readRDS(here("data", "inputs_dr.rds"))

-------------------------------------------------------------------------------
  
# 3. Data wrangling ----

# use p = 0.5 threshold to assign boosting status 

# any evidence of a boost (p ≥ 0.5) 

inputs <- inputs %>% 
  mutate(deng_boost = if_else(
    pmax(deng1_prob, deng2_prob, deng3_prob, deng4_prob, na.rm = TRUE) >= 0.5, 1L, 0L))

# 5.78% with evidence of boosting to any serotype (2021-22)
round(prop.table(table(inputs$deng_boost, useNA = "ifany")) * 100, 2)

-------------------------------------------------------------------------------

# 3. Data recategorisation ----

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

-------------------------------------------------------------------------------
  
# 4. Descriptive analyses ----

# Outcome is binary (yes/no)
table(inputs$deng_boost, useNA = "ifany") # 1 NA values present 

# Demographic characteristics 
table(inputs$age_group, useNA = "ifany")
table(inputs$age_u20, useNA = "ifany")
table(inputs$gender_cohort, useNA = "ifany")
table(inputs$setting_sero, useNA = "ifany")
table(inputs$occupation_sero, useNA = "ifany")
table(inputs$education_group, useNA = "ifany")
table(inputs$province, useNA = "ifany")
table(inputs$region, useNA = "ifany")

# Stratify by age 
age_bins   <- c(0, 10, 20, 30, 40, 50, 60, 70, Inf)
age_labels <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")

binom_ci <- function(x, n) {
  if (n == 0) return("NA")
  ht       <- binom.test(x, n, conf.level = 0.95)
  mean_pct <- signif(100 * ht$estimate, 3)
  ci1      <- signif(100 * ht$conf.int[1], 3)
  ci2      <- signif(100 * ht$conf.int[2], 3)
  paste0(mean_pct, "% (", ci1, "-", ci2, "%)")
}

# Tabulate boosts by age 
table_boost <- inputs %>%
  group_by(age_group) %>%
  summarise(
    N       = n(),
    n_boost = sum(deng_boost, na.rm = TRUE),
    propn   = binom_ci(n_boost, N),
    .groups = "drop"
  ) %>%
  bind_rows(
    inputs %>%
      summarise(
        age_group = "Total",
        N         = n(),
        n_boost   = sum(deng_boost, na.rm = TRUE),
        propn     = binom_ci(n_boost, N)
      )
  )

table_boost
write.csv(table_boost, here("outputs", "table_boost.csv"))

ft <- table_boost %>%
  flextable() %>%
  set_header_labels(
    age_group = "Age group",
    N         = "N",
    n_boost   = "N boosting",
    propn     = "Proportion boosting (95% CI)"
  ) %>%
  bold(part = "header") %>%
  bold(i = nrow(table_boost), part = "body") %>%  
  hline(i = nrow(table_boost) - 1) %>%            
  autofit()

save_as_docx(ft, path = here("outputs", "table_boost.docx"))

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
