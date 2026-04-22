# 03. Covariate analyses 

# Load the data ----

df_wide <- readRDS("Data/DR_surveillance/EN/df_wide")

-------------------------------------------------------------------------------------------------------------------------------
# 1. Plot probability of infection against age ---- 

# deng1 ----

# probability 
ggplot(data = df_wide, aes(x = age_cohort, y = prob_deng1)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "posterior probability of boosting", title = "deng1") +
  ylim(0,1) +
  theme_minimal()

# change in MFI
ggplot(data = df_wide, aes(x = age_cohort, y = change_dengns1_1/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "change in MFI", title = "deng1") +
  theme_minimal()

# initial MFI
ggplot(data = df_wide, aes(x = age_cohort, y = dengns1_1_mfi_sero/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "initial MFI", title = "deng1") +
  theme_minimal() + 
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 60))

ggplot(data = df_wide, aes(x = age_cohort, y = dengns1_1_mfi_cohort/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "final MFI", title = "deng1") +
  theme_minimal() + 
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 60))

# initial MFI vs change in MFI
ggplot(data = df_wide, aes(x = dengns1_1_mfi_sero, y = change_dengns1_1/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "initial MFI", y = "change in MFI", title = "deng1") +
  theme_minimal()

# deng2 ----

# probability 
ggplot(data = df_wide, aes(x = age_cohort, y = prob_deng2)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "posterior probability of boosting",title = "deng2") +
  ylim(0,1) +
  theme_minimal()

# change in MFI
ggplot(data = df_wide, aes(x = age_cohort, y = change_dengns1_2/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "change in MFI", title = "deng2") +
  theme_minimal()

# initial MFI
ggplot(data = df_wide, aes(x = age_cohort, y = dengns1_2_mfi_sero/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "initial MFI", title = "deng2") +
  theme_minimal()

# initial MFI vs change in MFI
ggplot(data = df_wide, aes(x = dengns1_2_mfi_sero, y = change_dengns1_2/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "initial MFI", y = "change in MFI", title = "deng2") +
  theme_minimal()

# deng3
ggplot(data = df_wide, aes(x = age_cohort, y = prob_deng3)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "posterior probability of boosting",title = "deng3") +
  ylim(0,1) +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = change_dengns1_3/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "change in MFI", title = "deng3") +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = dengns1_3_mfi_sero/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "initial MFI", title = "deng3") +
  theme_minimal()

ggplot(data = df_wide, aes(x = dengns1_3_mfi_sero, y = change_dengns1_3/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "initial MFI", y = "change in MFI", title = "deng3") +
  theme_minimal()

# deng4
ggplot(data = df_wide, aes(x = age_cohort, y = prob_deng4)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "posterior probability of boosting", title = "deng4") +
  ylim(0,1) +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = change_dengns1_4/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "change in MFI", title = "deng4") +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = dengns1_4_mfi_sero/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "initial MFI", title = "deng4") +
  theme_minimal()

ggplot(data = df_wide, aes(x = dengns1_4_mfi_sero, y = change_dengns1_4/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "initial MFI", y = "change in MFI", title = "deng4") +
  theme_minimal()

# zika
ggplot(data = df_wide, aes(x = age_cohort, y = prob_zika)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "posterior probability of boosting", title = "zika") +
  ylim(0,1) +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = change_zika_ns1/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "change in MFI", title = "zika") +
  theme_minimal()

ggplot(data = df_wide, aes(x = age_cohort, y = zika_ns1_mfi_sero/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "age", y = "initial MFI", title = "zika") +
  theme_minimal()

ggplot(data = df_wide, aes(x = zika_ns1_mfi_sero, y = change_zika_ns1/1000)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(se = TRUE) + 
  labs(x = "initial MFI", y = "change in MFI", title = "deng1") +
  theme_minimal()

-------------------------------------------------------------------------------------------------------------------------------
# 2. Location ---- 
  
# facet wrap if for each location 

# function 
facet_cov <- function(data, x_var, color, facet, title ){
  plot <- ggplot(data = data) + 
    geom_density(aes(x = {{x_var}}), fill = color, alpha = 0.5) + 
    labs(x = "∆MFI", 
         y = "Density",
         title = title) + 
    theme_minimal() + 
    facet_wrap(vars({{ facet }}))
    
  return(plot)
}

# DENV1
facet_cov(data = df_wide, x_var = change_dengns1_1, color = "red", facet = province, 
  title = "Change in DENV1 MFI values in DR by location")

# DENV2
facet_cov(data = df_wide, x_var = change_dengns1_2, color = "yellow", facet = province, 
          title = "Change in DENV2 MFI values in DR by location")

# DENV3
facet_cov(data = df_wide, x_var = change_dengns1_3, color = "green", facet = province, 
          title = "Change in DENV3 MFI values in DR by location")

# DENV4
facet_cov(data = df_wide, x_var = change_dengns1_4, color = "blue", facet = province, 
          title = "Change in DENV4 MFI values in DR by location")

# ZIKV
facet_cov(data = df_wide, x_var = change_zika_ns1, color = "violet", facet = province, 
          title = "Change in ZIKV MFI values in DR by location")

-------------------------------------------------------------------------------------------------------------------------------
# 4. Table 4 replication ---- 

# this table looks to see how attack rates (proportion infected and proportion seroconverting) differs by age

binom.calc <- function(x,n){
  htest <- binom.test(x,n, p = 1,conf.level=0.95)
  meanA=100*htest$estimate %>% as.numeric()  %>% signif(digits=3)
  conf1=100*htest$conf.int[1] %>% signif(digits=3)
  conf2=100*htest$conf.int[2] %>% signif(digits=3)
  paste(meanA,"% (",conf1,"-",conf2,"%)",sep="") 
}

binom.distn.calc <- function(p,n,digitsA=2){
  
  xx = seq(0,n,0.0001)
  yy = pbinom(xx,size=n,prob=p)
  meanA=100*p  %>% signif(digits=digitsA)
  conf1=(100*max(xx[yy<0.025])/n) %>% signif(digits=digitsA)
  conf2=(100*max(xx[yy<=0.975])/n) %>% signif(digits=digitsA)
  paste(meanA,"% (",conf1,"-",conf2,"%)",sep="") 
  
}

inputs0 <- df_wide
age.bins <- c(seq(0, 70, 10), 120)  # 70+ as last bin

results <- NULL

for(ii in 1:(length(age.bins) - 1)){
  
  age.ID <- inputs0$age_sero >= age.bins[ii] & inputs0$age_sero < age.bins[ii+1]
  
  n.total          <- sum(age.ID, na.rm = TRUE)
  prob             <- inputs0$prob_deng1[age.ID]
  prop.inf         <- sum(prob, na.rm = TRUE) / n.total
  
  sero.neg.r1      <- sum(inputs0$dengns1_1_seropos_sero[age.ID]  == 0, na.rm = TRUE)
  sero.convert     <- sum(inputs0$dengns1_1_seropos_sero[age.ID]  == 0 & 
                            inputs0$dengns1_1_seropos_cohort[age.ID] == 1, na.rm = TRUE)
  
  ci.inf    <- binom.distn.calc(prop.inf, n.total)
  ci.sero   <- if(sero.neg.r1 > 0) binom.calc(sero.convert, sero.neg.r1) else "N/A"
  
  age.label <- if(ii == length(age.bins)-1) "70+" else paste0(age.bins[ii], "-", age.bins[ii+1]-1)
  
  results <- rbind(results, data.frame(
    Age             = age.label,
    N               = n.total,
    Propn_infected  = ci.inf,
    Seroneg         = sero.neg.r1,
    Seroconverted   = sero.convert,
    Serocon_pct     = ci.sero
  ))
}

# Total row
total.n       <- sum(results$N)
total.prop    <- sum(inputs0$prob_deng1, na.rm = TRUE) / total.n
total.seroneg <- sum(results$Seroneg)
total.serocon <- sum(results$Seroconverted)

results <- rbind(results, data.frame(
  Age            = "Total",
  N              = total.n,
  Propn_infected = binom.distn.calc(total.prop, total.n),
  Seroneg        = total.seroneg,
  Seroconverted  = total.serocon,
  Serocon_pct    = binom.calc(total.serocon, total.seroneg)
))

results$Propn_infected <- gsub("-Inf", "0", results$Propn_infected)

print(results)


sum(df_wide$dengns1_1_seropos_sero == 0, na.rm = TRUE)
sum(df_wide$dengns1_1_seropos_cohort == 0, na.rm = TRUE)

sum(df_wide$dengns1_1_seropos_sero == 1, na.rm = TRUE)
sum(df_wide$dengns1_1_seropos_cohort == 1, na.rm = TRUE)



  


