
# 03.regression 

# regression analyses using outputs from FMMs 

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
library("stringr"")

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


# c) seroconversion 

# use the MIA cut-offs from the data 

serotypes_conv <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4")
sero_vars <- c("dengns1_1", "dengns1_2", "dengns1_3", "dengns1_4")

table_seropos <- map2_dfr(serotypes_conv, sero_vars, function(name, var) {
  
  year1_var <- paste0(var, "_seropos_sero")
  year2_var <- paste0(var, "_seropos_cohort")
  
  n_total <- sum(!is.na(inputs_complete[[year1_var]]))
  
  n_pos1  <- sum(inputs_complete[[year1_var]] == 1, na.rm = TRUE)
  n_pos2  <- sum(inputs_complete[[year2_var]] == 1, na.rm = TRUE)
  
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



-------------------------------------------------------------------------------

# 6. Logistic regression ----

# model 1: use all participants 

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


-------------------------------------------------------------------------------
  
# 6. Forest plot ----

# separate CI into separate columns 
lr_deng_fmm <- lr_deng_fmm %>%
  mutate(
    OR  = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
    LCL = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
    UCL = as.numeric(str_extract(OR_CI, "[0-9.]+(?=\\))"))
  )

# label variables 
dr_vars <- c(
  "age_u20"          = "Age under 20",
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




