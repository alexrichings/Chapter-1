
# ---- Adaptation of work by O'Driscoll et al. (2025)

library(rstan)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(mvtnorm)
library(loo)
library(qs)
library(data.table)
library(bayesplot)
library(posterior)
library(here)

source(here('R', 'RFunctions.R'))

# --- Antigen columns and labels
mfi_cols  <- c('dengns1_1_mfi','dengns1_2_mfi','dengns1_3_mfi',
               'dengns1_4_mfi','zika_ns1_mfi','chik_e1_mfi')
pathogens <- c('DENV1','DENV2','DENV3','DENV4','ZIKA','CHIK')

# --- Declare present pathogens
# All three arboviruses have documented circulation in the Dominican Republic
# Starting with DENV1-4 + CHIK as present, ZIKA as absent (cross-reactor only)
# You can change this once you've seen initial results
present <- c('DENV1','DENV2','DENV3','DENV4','CHIK')
nonpres <- pathogens[!pathogens %in% present]        # 'ZIKA'
pathogens <- c(present, nonpres)                     # present pathogens always first
mfi_cols  <- c(mfi_cols[pathogens %in% present], 
               mfi_cols[pathogens %in% nonpres])     # reorder to match

# --- Remove rows with missing age or region
df_model <- df[!is.na(df$age) & !is.na(df$region), ]
cat('N after removing missing age/region:', nrow(df_model), '\n')

# --- Location index from region
df_model$locID <- as.integer(factor(df_model$region))
loc <- levels(factor(df_model$region))              # 'a','b','c','d','e'

# --- Age groups
df_model$ageG <- cut(df_model$age,
                     breaks = c(0, 10, 20, 30, 40, 50, 60, Inf),
                     labels = FALSE,
                     right  = FALSE)
ageG_labels <- c('0-9','10-19','20-29','30-39','40-49','50-59','60+')

# --- Remove any remaining rows with missing age group
df_model <- df_model[!is.na(df_model$ageG), ]
cat('N after removing missing age group:', nrow(df_model), '\n')

# --- Log10-transform MFI (clamp negatives and zeros to 1 first)
mfi_raw       <- df_model[, mfi_cols]
mfi_clamped   <- as.data.frame(lapply(mfi_raw, function(x) pmax(x, 1)))
mfi_log       <- log10(mfi_clamped)
colnames(mfi_log) <- pathogens

# Quick check of transformed range
cat('\nLog10 MFI ranges:\n')
print(sapply(mfi_log, range))

# --- Compile Stan data list
data <- list()
data$y    <- as.matrix(mfi_log)
data$pres <- c(rep(1, length(present)), rep(0, length(nonpres)))
data$N    <- nrow(data$y)
data$nP   <- ncol(data$y)
data$nPp  <- sum(data$pres)
data$nL   <- length(loc)                             # 5 regions
data$nA   <- length(ageG_labels)                     # 7 age groups
data$ageG <- df_model$ageG
data$loc  <- df_model$locID
data$NperL  <- as.vector(table(df_model$locID))
data$NperLA <- t(table(df_model$locID, df_model$ageG))
data$ageProp <- as.matrix(table(df_model$locID, df_model$ageG) / data$NperL)

# --- Infection status combination matrix
data$infM <- inf_matrix(data$nP, pres = data$pres)
data$nC   <- nrow(data$infM)

# Indexing helpers
npos <- rowSums(data$infM)
wpos <- wneg <- matrix(0, ncol = data$nP, nrow = data$nC)
for(c in 1:nrow(data$infM)) for(p in 1:data$nP){
  if(npos[c] > 0)         wpos[c, 1:npos[c]]             <- which(data$infM[c,] == 1)
  if(npos[c] < data$nP)   wneg[c, 1:(data$nP - npos[c])] <- which(data$infM[c,] == 0)
}
data$npos <- npos
data$wpos <- wpos
data$wneg <- wneg

# --- Sanity checks
cat('\n--- Sanity checks ---\n')
cat('N individuals:              ', data$N,   '\n')
cat('N antigens:                 ', data$nP,  '\n')
cat('N present pathogens:        ', data$nPp, '\n')
cat('N locations (regions):      ', data$nL,  '\n')
cat('N age groups:               ', data$nA,  '\n')
cat('N infection status combos:  ', data$nC,  '\n')
cat('Any NAs in y?               ', anyNA(data$y), '\n')
cat('Any -Inf in y?              ', any(is.infinite(data$y)), '\n')
cat('\nIndividuals per region:\n')
print(data$NperL)
cat('\nIndividuals per region x age group:\n')
print(data$NperLA)


write.csv(df_model[, c('age', 'region', mfi_cols)],
          'data/serology_data.csv',
          row.names = FALSE)

list.files('StanModels')

dir.create('Outputs', showWarnings = FALSE)
