
# Chapter 1 ----

--------------------------------------------------------------------------------
  
# Objective: tabulate the number of people seropositive for the Fiji data and 
# those who seroconvert 
  
--------------------------------------------------------------------------------
  
# Load data and packages ----

library(here)
library(tidyverse)
library(dplyr)
library(PropCIs)
library(flextable)

inputs <- readRDS(here("data", "inputs_fj.rds"))

--------------------------------------------------------------------------------
# Section 1: Table for serostatus and seroconversion ----
--------------------------------------------------------------------------------

# 1.1. Functions for analysis  ----

# Binomial confidence intervals ----
binom_ci_raw <- function(x, n) {
  if (n == 0 | is.na(x) | is.na(n)) return(NA_character_)
  ht <- binom.test(x, n, conf.level = 0.95)
  paste0(
    formatC(100 * ht$estimate, digits = 1, format = "f"), "% (",
    formatC(100 * ht$conf.int[1], digits = 1, format = "f"), "-",
    formatC(100 * ht$conf.int[2], digits = 1, format = "f"), "%)"
  )
}

# Confidence intervals for differences in proportion ----
diff_ci_paired <- function(x1, n, x2, n_both) {
  if (any(is.na(c(x1, x2, n, n_both))) || n == 0) return(NA_character_)
  
  n10 <- x1 - n_both    # pos at t1 only
  n01 <- x2 - n_both    # pos at t2 only
  
  diff <- (x2 - x1) / n
  ci   <- diffpropci.mp(n10, n01, n, conf.level = 0.95)
  
  paste0(
    formatC(100 * diff,           digits = 2, format = "f"), "% (",
    formatC(100 * ci$conf.int[1], digits = 2, format = "f"), "-",
    formatC(100 * ci$conf.int[2], digits = 2, format = "f"), "%)"
  )
}
  
--------------------------------------------------------------------------------
  
# 1.2. Seroconversion ----

# Define names for tables
serotypes_conv <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "DENV")
sero_vars      <- c("DENV1", "DENV2", "DENV3", "DENV4", "any")

table_seropos_2013_2015 <- map2_dfr(serotypes_conv, sero_vars, function(name, var) {
  
  if (var == "any") {
    year1_vars <- c("DENV1P",  "DENV2P",  "DENV3P",  "DENV4P")
    year2_vars <- c("DENV1BP", "DENV2BP", "DENV3BP", "DENV4BP")
    
    n_total  <- sum(rowSums(!is.na(inputs[year1_vars])) > 0)
    n_pos1   <- sum(rowSums(inputs[year1_vars]  == 1, na.rm = TRUE) > 0)
    n_pos2   <- sum(rowSums(inputs[year2_vars]  == 1, na.rm = TRUE) > 0)
    n_both   <- sum(
      rowSums(inputs[year1_vars] == 1, na.rm = TRUE) > 0 &
        rowSums(inputs[year2_vars] == 1, na.rm = TRUE) > 0
    )
    
  } else {
    year1_var <- paste0(var, "P")
    year2_var <- paste0(var, "BP")
    
    n_total  <- sum(!is.na(inputs[[year1_var]]))
    n_pos1   <- sum(inputs[[year1_var]] == 1, na.rm = TRUE)
    n_pos2   <- sum(inputs[[year2_var]] == 1, na.rm = TRUE)
    n_both   <- sum(inputs[[year1_var]] == 1 & inputs[[year2_var]] == 1, na.rm = TRUE)
  }
  
  tibble(
    Serotype   = name,
    N          = n_total,
    n_2013     = n_pos1,
    sero_2013  = binom_ci_raw(n_pos1, n_total),
    n_2015     = n_pos2,
    sero_2015  = binom_ci_raw(n_pos2, n_total),
    Difference = diff_ci_paired(n_pos1, n_total, n_pos2, n_both)
  )
})

write.csv(table_seropos_2013_2015, here("outputs/tables/Fiji/", "table_seropos_2013_2015.csv"))

st_2013_2015 <- table_seropos_2013_2015 %>%
  flextable() %>%
  set_header_labels(
    Serotype   = "Serotype",
    N          = "N",
    n_2013     = "Sero(+) (2013)",
    sero_2013  = "Sero(+) % (2013)",
    n_2015     = "Sero(+) (2015)",
    sero_2015  = "Sero(+) % (2015)",
    Difference = "Difference"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

save_as_docx(st_2013_2015, path = here("outputs/tables/Fiji/", "table_seropos_2013_2015.docx"))
  
# 1.5. Seroconversion ----

# Define names for tables
serotypes_conv <- c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "DENV")
sero_vars      <- c("DENV1", "DENV2", "DENV3", "DENV4", "any")

table_seropos_2013_2015 <- map2_dfr(serotypes_conv, sero_vars, function(name, var) {
  
  if (var == "any") {
    
    year1_vars <- c("DENV1P",  "DENV2P",  "DENV3P",  "DENV4P")
    year2_vars <- c("DENV1BP", "DENV2BP", "DENV3BP", "DENV4BP")
    
    n_total <- sum(rowSums(!is.na(inputs[year1_vars])) > 0)
    
    n_pos1 <- sum(rowSums(inputs[year1_vars] == 1, na.rm = TRUE) > 0)
    n_pos2 <- sum(rowSums(inputs[year2_vars] == 1, na.rm = TRUE) > 0)
    
  } else {
    
    year1_var <- paste0(var, "P")   # e.g. DENV1P
    year2_var <- paste0(var, "BP")  # e.g. DENV1BP
    
    n_total <- sum(!is.na(inputs[[year1_var]]))
    
    n_pos1  <- sum(inputs[[year1_var]] == 1, na.rm = TRUE)
    n_pos2  <- sum(inputs[[year2_var]] == 1, na.rm = TRUE)
  }
  
  tibble(
    Serotype   = name,
    N          = n_total,
    n_2013     = n_pos1,
    sero_2013  = binom_ci_raw(n_pos1, n_total),
    n_2015     = n_pos2,
    sero_2015  = binom_ci_raw(n_pos2, n_total),
    Difference = diff_ci(n_pos2, n_total, n_pos1, n_total)
  )
})

write.csv(table_seropos_2013_2015, here("outputs", "table_seropos_2013_2015.csv"))

st_2013_2015 <- table_seropos_2013_2015 %>%
  flextable() %>%
  set_header_labels(
    Serotype   = "Serotype",
    N          = "N",
    n_2013     = "Sero(+) (2013)",
    sero_2013  = "Sero(+) % (2013)",
    n_2015     = "Sero(+) (2015)",
    sero_2015  = "Sero(+) % (2015)",
    Difference = "Difference"
  ) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  autofit()

save_as_docx(st_2013_2015, path = here("outputs", "table_seropos_2013_2015.docx"))


