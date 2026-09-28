
# Chapter 1 

-------------------------------------------------------------------------------
  
# Script: 03-boost-regression ----  

# Purpose: 
# Tabulate exposures (including by different covariates), tabulate boosts, 
# tabulate exposures, run regression 

# Output: 
# Individual probability of boosting for each person per dengue serotype (inputs_dr) and
# fmm parameters, regression tables 

-------------------------------------------------------------------------------


# 1. Install and load packages ----

library("ggplot2")
library("dplyr")
library("here")
library("flextable")
library("officer")
library("tidyverse")
library("officer")
library("tidyr")
library("stringr")
library("mgcv")
library("mgcViz")
library("scales")

-------------------------------------------------------------------------------
  
# 2. Load data ----

inputs <- readRDS(here("data", "inputs_dr.rds"))

-------------------------------------------------------------------------------
  
# 3. Data wrangling ----

# use p = 0.5 threshold to assign boosting status 

# any evidence of a boost (p ≥ 0.5) 

inputs <- inputs %>% 
  mutate(deng_boost = if_else(
    pmax(as.numeric(deng1_prob), 
         as.numeric(deng2_prob), 
         as.numeric(deng3_prob), 
         as.numeric(deng4_prob), na.rm = TRUE) >= 0.5,1L, 0L))

# 5.78% with evidence of boosting to any serotype (2021-22)
round(prop.table(table(inputs$deng_boost, useNA = "ifany")) * 100, 2)

-------------------------------------------------------------------------------

# 3. Data recategorisation ----

# Recategorise 

# age 
inputs <- inputs %>%
  mutate(age_group = cut(age_cohort, breaks = c(6, 20, 30, 40, 50, 60, 70, Inf), right = FALSE,
                         labels = c("6–20", "20–30", "30–40", "40–50","50–60", "60–70", "70+")))

-------------------------------------------------------------------------------
  
# 4. Tabulating by exposure  ----

# confidence intervals 

# binomial confidence intervals 
binom_ci_raw <- function(x, n) {
  if (n == 0 | is.na(x) | is.na(n)) return(NA_character_)
  ht <- binom.test(x, n, conf.level = 0.95)
  paste0(
    formatC(100 * ht$estimate, digits = 1, format = "f"), "% (",
    formatC(100 * ht$conf.int[1], digits = 1, format = "f"), "-",
    formatC(100 * ht$conf.int[2], digits = 1, format = "f"), "%)"
  )
}

# difference in proportion confidence intervals
diff_ci <- function(x1, n1, x2, n2) {
  if (any(c(n1, n2) == 0) || any(is.na(c(x1, x2, n1, n2)))) return(NA_character_)
  pt <- prop.test(c(x1, x2), c(n1, n2), conf.level = 0.95, correct = FALSE)
  
  est <- pt$estimate[1] - pt$estimate[2]
  ci1 <- max(0, pt$conf.int[1])
  ci2 <- pt$conf.int[2]
  
  paste0(
    formatC(100 * est, digits = 2, format = "f"), "% (",
    formatC(100 * ci1, digits = 2, format = "f"), "-",
    formatC(100 * ci2, digits = 2, format = "f"), "%)"
  )
}

# get complete set of data 
complete_ids <- inputs %>% 
  filter(!is.na(deng_boost), !is.na(deng1_prob), !is.na(dengns1_1_seropos_sero)) %>%
  pull(serosurvey_id)

inputs_complete <- inputs %>% filter(serosurvey_id %in% complete_ids)

# a) tabulate by boost ----

serotypes <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "DENV")
boost_vars <- c("deng1_prob", "deng2_prob", "deng3_prob", "deng4_prob", "any")

table_boost <- map2_dfr(serotypes, boost_vars, function(name, var) {
  
  if (var == "any") {
    n_total <- sum(!is.na(inputs_complete$deng_boost))
    n_boost <- sum(inputs_complete$deng_boost == 1, na.rm = TRUE)
  } else {
    prob_col <- inputs_complete[[var]]
    n_total  <- sum(!is.na(prob_col))
    n_boost  <- sum(prob_col >= 0.5, na.rm = TRUE)
  }
  
  tibble(
    Serotype = name,
    N        = n_total,
    n_boost  = n_boost,
    propn    = binom_ci_raw(n_boost, n_total)
  )
})

table_boost
write.csv(table_boost, here("outputs", "table_boost.csv"))

# save table to 
bt <- table_boost %>%
  flextable() %>%
  set_header_labels(Serotype = "Serotype", N = "N", n_boost = "N boost", propn = "Boosting %") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

bt
save_as_docx(bt, path = here("outputs", "table_boost.docx"))



# b) tabulate boosts by age ----

# stratify
age_bins   <- c(0, 10, 20, 30, 40, 50, 60, 70, Inf)
age_labels <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")

table_boost_age <- inputs_complete %>%
  group_by(age_group) %>%
  summarise(
    N       = n(),
    n_boost = sum(deng_boost, na.rm = TRUE),
    propn   = binom_ci_raw(n_boost, N),
    .groups = "drop"
  ) %>%
  bind_rows(
    inputs_complete %>%
      summarise(
        age_group = "Total",
        N         = n(),
        n_boost   = sum(deng_boost, na.rm = TRUE),
        propn     = binom_ci_raw(n_boost, N)
      )
  )

table_boost_age
write.csv(table_boost_age, here("outputs", "table_boost_age.csv"))

bat <- table_boost_age %>%
  flextable() %>%
  set_header_labels(
    age_group = "Age group",
    N         = "N",
    n_boost   = "N boosting",
    propn     = "Proportion boosting (95% CI)"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  bold(i = nrow(table_boost_age), part = "body") %>%  
  hline(i = nrow(table_boost_age) - 1) %>%            
  autofit()

save_as_docx(bat, path = here("outputs", "table_boost_age.docx"))

# plot 

table_boost_age %>%
  filter(age_group != "Total") %>%
  mutate(
    age_group = factor(age_group, levels = unique(age_group)),
    nums  = str_extract_all(propn, "[0-9]+\\.?[0-9]*"),
    pct   = as.numeric(map_chr(nums, 1)),
    lower = as.numeric(map_chr(nums, 2)),
    upper = as.numeric(map_chr(nums, 3))
  ) %>%
  ggplot(aes(x = age_group, y = pct)) +
  geom_col(fill = "#185FA540", colour = "#185FA5", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, colour = "#185FA5") +
  scale_y_continuous(limits = c(0, 20), labels = function(x) paste0(x, "%")) +
  labs(x = "Age group", y = "Proportion boosted (%)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

# c) seroconversion ----

# use the MIA cut-offs from the data 

serotypes_conv <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "DENV")
sero_vars <- c("dengns1_1", "dengns1_2", "dengns1_3", "dengns1_4", "any")

table_seropos <- map2_dfr(serotypes_conv, sero_vars, function(name, var) {
  
  if (var == "any") {
    
    year1_vars <- paste0(c("dengns1_1", "dengns1_2", "dengns1_3", "dengns1_4"), "_seropos_sero")
    year2_vars <- paste0(c("dengns1_1", "dengns1_2", "dengns1_3", "dengns1_4"), "_seropos_cohort")
   
     n_total <- sum(!is.na(inputs_complete$deng_boost))
    
    # find any seropositive 
    n_pos1 <- sum(rowSums(inputs_complete[year1_vars] == 1, na.rm = TRUE) > 0)
    n_pos2 <- sum(rowSums(inputs_complete[year2_vars] == 1, na.rm = TRUE) > 0)
    
  } else {
  
    year1_var <- paste0(var, "_seropos_sero")
    year2_var <- paste0(var, "_seropos_cohort")
    
    n_total <- sum(!is.na(inputs_complete[[year1_var]]))
    
    n_pos1  <- sum(inputs_complete[[year1_var]] == 1, na.rm = TRUE)
    n_pos2  <- sum(inputs_complete[[year2_var]] == 1, na.rm = TRUE)
  }
  
  tibble(
    Serotype   = name,
    N = n_total,
    n_2021 = n_pos1,
    sero_2021 = binom_ci_raw(n_pos1, n_total),
    n_2022 = n_pos2,
    sero_2022 = binom_ci_raw(n_pos2, n_total),
    Difference = diff_ci(n_pos2, n_total, n_pos1, n_total)  
  )
})

table_seropos
write.csv(table_seropos, here("outputs", "table_seropos.csv"))

st <- table_seropos %>%
  flextable() %>%
  set_header_labels(
    Serotype = "Serotype",
    N = "N",
    n_2021 = "Sero(+) (2021)",
    sero_2021 = "Sero(+) % (2021)",
    n_2022 = "Sero(+) (2022)",
    sero_2022 = "Sero(+) % (2022)",
    Difference = "Difference"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

st
save_as_docx(st, path = here("outputs", "table_seropos.docx"))

-------------------------------------------------------------------------------

# 5. Descriptive analysis ----

names(inputs)

# re-group 

inputs <- inputs %>%
  mutate(age_u20 = cut(age_cohort, breaks = c(6, 20, Inf), right = FALSE,
                       labels = c("under 20", "over 20"))) %>%
  mutate(age_u20 = relevel(factor(age_u20), ref = "over 20"))

# education 
inputs <- inputs %>%
  mutate(
    education_group = as.character(education_sero),
    education_group = na_if(education_group, "dont_know_refuse"),
    education_group = na_if(education_group, "NA"),
    education_group = factor(education_group)
  )

# look at low education vs some degree of education 
inputs <- inputs %>%
  mutate(
    education_group2 = case_when(
      education_group %in% c("no_formal", "primary") ~ "Primary and under",
      education_group %in% c("secondary", "technical", "university") ~ "Beyond primary",
      TRUE ~ NA_character_
    ),
    education_group2 = factor(education_group2)
  )

# "Other" category too small for regression 
inputs <- inputs %>%
  mutate(
    gender = ifelse(gender_sero == "Other", NA, gender_sero),
    gender = factor(gender)
  )



# age, gender, setting, education, province, region

table(inputs$age_u20, useNA = "ifany")
table(inputs$gender, useNA = "ifany")
table(inputs$setting_sero, useNA = "ifany")
table(inputs$education_group2, useNA = "ifany")
table(inputs$province, useNA = "ifany")
table(inputs$region, useNA = "ifany")


# GAMs ----

# initial titre against change 

# initial titre against probability of boosting 
gam_denv1_start_boost <- getViz(gam(data = inputs, formula = deng1_prob ~ s(dengns1_1_mfi_sero), family = betar(link = "logit")))
gam_denv2_start_boost <- getViz(gam(data = inputs, formula = deng2_prob ~ (dengns1_2_mfi_sero), family = betar(link = "logit")))
gam_denv3_start_boost <- getViz(gam(data = inputs, formula = deng3_prob ~ s(dengns1_3_mfi_sero), family = betar(link = "logit")))
gam_denv4_start_boost <- getViz(gam(data = inputs, formula = deng4_prob ~ s(dengns1_4_mfi_sero), family = betar(link = "logit")))

# linear relatinship between initial MFI and probability of boosting 
plot(gam_denv1_start_boost, allTerms = TRUE, seWithMean = TRUE) + 
  l_ciPoly(alpha = 0.7, fill = "lightblue") +
  l_fitLine(linetype = 1) +
  theme_minimal() +
  labs(
    title = "DENV1",
    x = "Initial MFI",
    y = "Log odds of mean probability of boosting"
  )

# initial MFI against boost (binary)
gam_denv_start_boost <- getViz(gam(data = inputs, formula = deng_boost ~ s(dengns1_1_mfi_sero), family = binomial(link = "logit")))
plot(gam_denv_start_boost, allTerms = TRUE, seWithMean = TRUE) + 
  l_ciPoly(alpha = 0.7, fill = "lightblue") +
  l_fitLine(linetype = 1) +
  theme_minimal() +
  labs(
    title = "DENV",
    x = "Initial MFI",
    y = "Log odds of boosting"
  )
                               
                               
                               
# age against boosting 
gam_age_boost_viz <- getViz(gam(data = inputs, formula = deng_boost ~ s(age_cohort), family = "binomial"))
gam_age_boost <- gam(deng_boost ~ s(age_cohort), data = inputs, family = "binomial")

plot(gam_age_boost_viz, allTerms = T, seWithMean = T) + 
  l_ciPoly(alpha=0.7) +
  l_fitLine(linetype = 1)  +
  # l_ciBar() +
  l_rug() +
  theme_bw() +
  xlab("Age") +
  ylab("Log odds of boosting") 

# prediction grid 
age_grid <- data.frame(age_cohort = seq(min(inputs$age_cohort), max(inputs$age_cohort), length.out = 200))

# predict on response scale 
pred <- predict(m, newdata = newdat, se.fit = TRUE, type = "response")

plot_dat <- newdat %>%
  mutate(
    fit = pred$fit,
    se = pred$se.fit,
    upper = fit + 1.96 * se,
    lower = fit - 1.96 * se)

ggplot(plot_dat, aes(x = age_cohort, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4C78A8", alpha = 0.2) +
  geom_line(color = "black", size = 1) +
  geom_rug(data = inputs, aes(x = age_cohort), inherit.aes = FALSE, alpha = 0.15) +
  theme_minimal(base_size = 18) +
  labs(
    x = "Age",
    y = "Probability of boosting"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(here("outputs", "age_boost_gam.png"), width = 8, height = 6, dpi = 300)

-------------------------------------------------------------------------------

# 6. Logistic regression ----

# model 1a: use all participants - univariable 

model1.1 <- glm(deng_boost ~ age_u20, data = inputs, family = "binomial")
model1.2 <- glm(deng_boost ~ gender, data = inputs, family = "binomial")
model1.3 <- glm(deng_boost ~ setting_sero, data = inputs, family = "binomial")
model1.4 <- glm(deng_boost ~ education_group2, data = inputs, family = "binomial")
#model1.5 <- glm(deng_boost ~ province, data = inputs, family = "binomial")
model1.6 <- glm(deng_boost ~ region, data = inputs, family = "binomial")

model.list=list(model1.1, model1.2, model1.3, model1.4, model1.6)
names.model=c("age_u20","gender_sero","setting_sero", "education_group2", "region")

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


# model 1b: use all participants - multivariable regression 

model1.1a <- glm(deng_boost ~ age_u20, data = inputs, family = "binomial")
model1.2a <- glm(deng_boost ~ gender, data = inputs, family = "binomial")
model1.3a <- glm(deng_boost ~ setting_sero + age_u20 + region, data = inputs, family = "binomial")
model1.4a <- glm(deng_boost ~ education_group2 + age_u20 + region + gender + setting_sero, data = inputs, family = "binomial")
#model1.5a <- glm(deng_boost ~ province, data = inputs, family = "binomial")
model1.6a <- glm(deng_boost ~ region + age_u20, data = inputs, family = "binomial")


model.list <- list(model1.1a, model1.2a, model1.3a, model1.4a, model1.6a)
exposure.vars <- c("age_u20", "gender", "setting_sero", "education_group2", "region")

store.oddsTable <- NULL

for (jj in seq_along(model.list)) {
  
  modelT  <- model.list[[jj]]
  varname <- exposure.vars[jj]
  
  estm <- coef(modelT)
  conf <- confint.default(modelT)
  pval <- coef(summary(modelT))[, 4]
  
  # pull the exposure coefficient only (by name, not position)
  expo_coef <- estm[grep(varname, names(estm))]
  expo_conf <- conf[grep(varname, rownames(conf)), , drop = FALSE]
  expo_pval <- pval[grep(varname, names(pval))]
  
  # counts and reference level
  expo_data  <- model.frame(modelT)[[varname]]
  ref_level  <- levels(factor(expo_data))[1]
  test_level <- levels(factor(expo_data))[2]
  
  out_row <- data.frame(
    var       = varname,
    reference = ref_level,
    comparison= test_level,
    N_used    = nobs(modelT),
    N_ref     = sum(expo_data == ref_level, na.rm = TRUE),
    N_comp    = sum(expo_data == test_level, na.rm = TRUE),
    OR_CI     = paste0(round(exp(expo_coef), 2), " (",
                       round(exp(expo_conf[1]), 2), "-",
                       round(exp(expo_conf[2]), 2), ")"),
    p         = round(expo_pval, 3),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  store.oddsTable <- rbind(store.oddsTable, out_row)
}

store.oddsTable

write.csv(store.oddsTable, here("outputs", "mlr_deng_fmm.csv"))
mlr_deng_fmm <- read.csv(here("outputs", "mlr_deng_fmm.csv"))

# write into tabel format 

lr  <- read.csv(here("outputs", "lr_deng_fmm.csv"))
mlr <- read.csv(here("outputs", "mlr_deng_fmm.csv"))

table_or <- data.frame(
  Variable   = c("Age", "Gender", "Setting", "Education", "Region"),
  Reference  = c("Over 20", "Female", "Rural", "Beyond primary", "B"),
  Comparison = c("Under 20", "Male", "Urban", "Primary and under", "D"),
  N_unadj    = lr$N_used,
  uOR_CI     = lr$OR_CI,
  u_p        = lr$p,
  N_adj      = mlr$N_used,
  aOR_CI     = mlr$OR_CI,
  a_p        = mlr$p
)

st <- table_or %>%
  flextable() %>%
  set_header_labels(
    Variable   = "Variable",
    Reference  = "Reference",
    Comparison = "Comparison",
    N_unadj    = "N (unadjusted)",
    uOR_CI     = "uOR (95% CI)",
    u_p        = "p",
    N_adj      = "N (adjusted)",
    aOR_CI     = "aOR (95% CI)",
    a_p        = "p"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

st

save_as_docx(st, path = here("outputs", "lr_comparison_OR.docx"))


----
  
  table_or <- data.frame(
    Variable = c(
      "Age / years", "  Over 20 (ref)", "  Under 20",
      "Gender", "  Female (ref)", "  Male",
      "Setting", "  Rural (ref)", "  Urban",
      "Education", "  Beyond primary (ref)", "  Primary or less",
      "Region", "  Southwest (ref)", "  Northeast"
    ),
    N_unadj = c(
      "", "", "1,019",
      "", "", "1,007",
      "", "", "1,019",
      "", "", "911",
      "", "", "1,019"
    ),
    uOR_CI = c(
      "", "", "2.33 (1.26-4.14)",
      "", "", "1.54 (0.89-2.61)",
      "", "", "1.05 (0.62-1.79)",
      "", "", "0.99 (0.53-1.82)",
      "", "", "1.00 (0.59-1.72)"
    ),
    u_p = c(
      "", "", "0.001",
      "", "", "0.120",
      "", "", "0.850",
      "", "", "0.980",
      "", "", "0.990"
    ),
    N_adj = c(
      "", "", "1,019",
      "", "", "1,007",
      "", "", "1,019",
      "", "", "900",
      "", "", "1,019"
    ),
    aOR_CI = c(
      "", "", "2.33 (1.29-4.21)",
      "", "", "1.54 (0.90-2.62)",
      "", "", "1.05 (0.62-1.78)",
      "", "", "0.92 (0.49-1.72)",
      "", "", "0.87 (0.51-1.50)"
    ),
    a_p = c(
      "", "", "0.005",
      "", "", "0.116",
      "", "", "0.859",
      "", "", "0.798",
      "", "", "0.623"
    )
  )

# rows that are variable-level headers (bold)
header_rows <- c(1, 4, 7, 10, 13)

st <- table_or %>%
  flextable() %>%
  set_header_labels(
    Variable = "Variable",
    N_unadj  = "N (unadjusted)",
    uOR_CI   = "uOR (95% CI)",
    u_p      = "p-value",
    N_adj    = "N (adjusted)",
    aOR_CI   = "aOR (95% CI)",
    a_p      = "p"
  ) %>%
  bold(part = "header") %>%
  bold(i = header_rows, j = "Variable", part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

st

save_as_docx(st, path = here("outputs", "lr_comparison_OR.docx"))




# model 2: use seronegatives 

# defined from the cut offs 
table(inputs$dengns1_1_seropos_sero, useNA = "ifany") # 53
table(inputs$dengns1_2_seropos_sero, useNA = "ifany") # 59 
table(inputs$dengns1_3_seropos_sero, useNA = "ifany") # 39 
table(inputs$dengns1_4_seropos_sero, useNA = "ifany") # 73 

inputs_seroneg <- inputs %>%
  dplyr::filter(
    dengns1_1_seropos_sero == 0 &
      dengns1_2_seropos_sero == 0 &
      dengns1_3_seropos_sero == 0 &
      dengns1_4_seropos_sero == 0
  )

# 26 seronegative to all serotypes 
dim(inputs_seroneg)

model2.1 <- glm(deng_boost ~ age_u20, data = inputs_seroneg, family = "binomial")
model2.2 <- glm(deng_boost ~ gender, data = inputs_seroneg, family = "binomial")
model2.3 <- glm(deng_boost ~ setting_sero, data = inputs_seroneg, family = "binomial")
model2.4 <- glm(deng_boost ~ education_group2, data = inputs_seroneg, family = "binomial")
#model2.5 <- glm(deng_boost ~ province, data = inputs_seroneg, family = "binomial")
model2.6 <- glm(deng_boost ~ region, data = inputs_seroneg, family = "binomial")

model.list=list(model2.1, model2.2, model2.3, model2.4, model2.6)
names.model=c("age_u20","gender_sero","setting_sero", "education_group2", "region")

data.tally = inputs[,names.model]

# calculate ORs for the different models 
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

write.csv(store.oddsTable, here("outputs", "lr_deng_fmm_seroneg.csv"))
lr_deng_fmm_seroneg <- read.csv(here("outputs", "lr_deng_fmm_seroneg.csv"))


-------------------------------------------------------------------------------
  
# 6. Forest plot ----

# separate CI into separate columns 
lr_deng_fmm <- lr_deng_fmm %>%
  mutate(
    OR  = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
    LCL = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
    UCL = as.numeric(str_extract(OR_CI, "[0-9.]+(?=\\))"))
  )

lr_deng_fmm_seroneg <- lr_deng_fmm_seroneg %>%
  mutate(
    OR  = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
    LCL = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
    UCL = as.numeric(str_extract(OR_CI, "[0-9.]+(?=\\))"))
  )

# label variables 
dr_vars <- c(
  "age_u20"          = "Age under 20*",
  "gender_sero"      = "Male",
  "setting_sero"     = "Urban",
  "education_group2" = "Primary education or less",
  "region"           = "Southeast region"
)

# plot

lr_deng_fmm <- lr_deng_fmm %>% mutate(outcome = "DENV") %>%
  mutate(var = fct_relevel(var, "region", "setting_sero", "education_group2", "gender_sero", "age_u20"))

forest_collapsed_lr_dr <- 
  ggplot(lr_deng_fmm, aes(x = OR, y = var)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, colour = "#5BA4A4") +
  geom_point(size = 2.8, colour = "#5BA4A4") +
  scale_y_discrete(labels = dr_vars) +
  scale_x_log10() +
  facet_grid(. ~ outcome) +
  labs(
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border       = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    strip.background   = element_rect(fill = "grey90", colour = "grey35", linewidth = 0.8),
    strip.text         = element_text(face = "bold"),
    strip.placement    = "outside",
    panel.spacing      = unit(0, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank()
  )

forest_collapsed_lr_dr
ggsave(here("outputs", "forest_collapsed_lr_dr.png"), width = 12, height = 6, dpi = 300)


# plot seronegatives and full participants 

# relevel the covariates 
lr_deng_fmm_seroneg <- lr_deng_fmm_seroneg %>% mutate(outcome = "DENV") %>%
  mutate(var = fct_relevel(var, "region", "setting_sero", "education_group2", "gender_sero", "age_u20"))

# join the two models 
lr_deng_fmm <- lr_deng_fmm %>%
  mutate(sample = "Total",
         source = "Current",
         model  = "LR")

lr_deng_fmm_seroneg <- lr_deng_fmm_seroneg %>%
  mutate(sample = "Seronegative",
         source = "Current",
         model  = "LR")

combined_fmm_df <- bind_rows(
  lr_deng_fmm,
  lr_deng_fmm_seroneg
) %>%
  mutate(SE = (log(UCL) - log(LCL)) / (2 * 1.96)) # add the standard error 

var_levels <- combined_fmm_df %>%
  distinct(var) %>%
  pull(var) %>%
  as.character()

combined_fmm_df <- combined_fmm_df %>%
  mutate(var = factor(var, levels = rev(var_levels))) %>%
  mutate(OR  = ifelse(var == "education_group2" & sample == "Seronegative", NA, OR),
         LCL = ifelse(var == "education_group2" & sample == "Seronegative", NA, LCL),
         UCL = ifelse(var == "education_group2" & sample == "Seronegative", NA, UCL))

# add the SE ratio to the data 
se_ratio_df <- combined_fmm_df %>%
  select(var, sample, SE) %>%
  pivot_wider(names_from = sample, values_from = SE) %>%
  mutate(SE_ratio = Total / Seronegative)

combined_fmm_df <- combined_fmm_df %>%
  left_join(se_ratio_df %>% select(var, SE_ratio), by = "var") %>%
  mutate(SE_ratio_label = sprintf("%.2f", SE_ratio))


pd <- position_dodge(width = 0.75)

forest_fmm_compare <-
  ggplot(combined_fmm_df,
         aes(x = OR, y = var,
             colour = sample,
             group  = sample)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), position = pd, height = 0.2) +
  geom_point(position = pd, size = 2.8) +
  scale_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1.0", "10.0", "100.0")
  ) +
  facet_grid(. ~ outcome) +
  scale_colour_manual(
    name   = "Data",
    values = c("Seronegative" = "#F08080",
               "Total"        = "#40BDB8"),
    labels = c("Seronegative" = "Seronegatives",
               "Total"        = "All participants")
  ) +
  scale_y_discrete(labels = c(
    "age_u20"          = "Age under 20*",
    "gender_sero"      = "Male",
    "education_group2" = "Primary education or less",
    "setting_sero"     = "Urban",
    "region"           = "Southeast region"
  )) +
  labs(
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.border       = element_rect(fill = NA, linewidth = 0.8, colour = "grey35"),
    panel.spacing      = unit(1.2, "lines"),
    strip.background   = element_rect(fill = "grey85", colour = "grey35", linewidth = 0.8),
    strip.text         = element_text(face = "bold"),
    strip.placement    = "outside",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "bottom", 
    legend.justification = "centre",
    legend.background    = element_rect(fill = "white", colour = "grey35", linewidth = 0.5),
  )  +
  annotate("segment", x = 0, xend = 1000, y = 1.8, yend = 1.8, linetype = "dashed", colour = "grey50", linewidth = 0.8) + 
  coord_cartesian(xlim = c(0.01, 100)) 

forest_fmm_compare

# convert to a gtable form 
gt <- ggplotGrob(forest_fmm_compare)
right_col <- 7
top_row   <- 10
bot_row   <- 10
strip_row <- 8

se_labels <- combined_fmm_df %>%
  filter(sample == "Total") %>%
  arrange(match(var, levels(combined_fmm_df$var))) %>%
  pull(SE_ratio_label) %>%
  rev()

n_rows <- length(se_labels)

header_grob <- grobTree(
  rectGrob(gp = gpar(fill = "grey85", col = "grey35", lwd = 0.8 * .pt)),
  textGrob("SE ratio", gp = gpar(fontface = "bold", fontsize = 18))
)

row_grobs <- lapply(se_labels, function(lab) {
  grobTree(
    rectGrob(gp = gpar(fill = "white", col = "grey35", lwd = 0.8 * .pt)),
    textGrob(lab, gp = gpar(fontsize = 18), x = 0.5, y = 0.5)
  )
})

se_column_grob <- frameGrob(
  layout = grid.layout(nrow = n_rows, ncol = 1,
                       heights = unit(rep(1, n_rows), "null"))
)

for (i in seq_along(row_grobs)) {
  se_column_grob <- placeGrob(se_column_grob, row_grobs[[i]], row = i, col = 1)
}

gt      <- gtable_add_cols(gt, unit(2.8, "cm"), pos = right_col)
new_col <- right_col + 1

gt <- gtable_add_grob(gt, header_grob,
                      t = strip_row, b = strip_row,
                      l = new_col,  r = new_col, name = "se_header")
gt <- gtable_add_grob(gt, se_column_grob,
                      t = top_row,  b = bot_row,
                      l = new_col,  r = new_col, name = "se_cells")

# plot
grid.newpage()
grid.draw(gt)

ggsave(here("outputs", "forest_fmm_compare.png"), plot = gt, width = 13, height = 7, dpi = 300)
