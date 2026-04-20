

# AK Suggestions ---- 

library(dplyr)
library(ggplot2)

# 1. Probility of boosting by age ----

x <- x %>%
  left_join(deng1_plat %>% mutate(id = as.integer(id)) %>% select(id, deng1_plateau), by = "id") %>%
  left_join(deng2_plat %>% mutate(id = as.integer(id)) %>% select(id, deng2_plateau), by = "id") %>%
  left_join(deng3_plat %>% mutate(id = as.integer(id)) %>% select(id, deng3_plateau), by = "id") %>%
  left_join(deng4_plat %>% mutate(id = as.integer(id)) %>% select(id, deng4_plateau), by = "id") %>%
  left_join(zika_plat %>% mutate(id = as.integer(id)) %>% select(id, zika_plateau), by = "id")

var_gam_plot <- function(data, x_var, y_var, x_lab, y_lab, title, xmin, xmax, ymin, ymax){
  plot <- 
    ggplot(data = data, aes(x = .data[[x_var]], y = .data[[y_var]])) + 
    geom_point(alpha = 0.5) + 
    geom_smooth(se = TRUE) + 
    labs(x = x_lab, y = y_lab, title = title) + 
    theme_minimal() + 
    coord_cartesian(xlim = c(xmin,xmax), ylim = c(ymin, ymax))
  
  return(plot)
}

# raw probability 

var_gam_plot(x, "age_calc", "p_deng1", "age", "probability", "deng1 boost prob (raw)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_deng2", "age", "probability", "deng2 boost prob (raw)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_deng3", "age", "probability", "deng3 boost prob (raw)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_deng4", "age", "probability", "deng4 boost prob (raw)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_zika", "age", "probability", "zika boot prob (raw)", 0, 100, 0, 1)

# plateau 
var_gam_plot(x, "age_calc", "deng1_plateau", "age", "probability", "deng1 boost prob (plateau)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "deng2_plateau", "age", "probability", "deng2 boost prob (plateau)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "deng3_plateau", "age", "probability", "deng3 boost prob (plateau)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "deng4_plateau", "age", "probability", "deng4 boost prob (plateau)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "zika_plateau", "age", "probability", "zika boost prob (plateau)", 0, 100, 0, 1)

# sigmoid
var_gam_plot(x, "age_calc", "p_corrected_deng1", "age", "probability", "deng1 boost prob (sigmoid)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_corrected_deng2", "age", "probability", "deng2 boost prob (sigmoid)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_corrected_deng3", "age", "probability", "deng3 boost prob (sigmoid)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_corrected_deng4", "age", "probability", "deng4 boost prob (sigmoid)", 0, 100, 0, 1)
var_gam_plot(x, "age_calc", "p_corrected_zika", "age", "probability", "zika boost prob (sigmoid)", 0, 100, 0, 1)


# 2. Probability of boosting by location ----

table(x$location, useNA = "ifany")

# count how many people are in each location 
x %>% group_by(location) %>% summarise(n = n())
x %>% summarise(n = n())

# maybe use the data online to work out who is from where? 

# 3. Primary and secondary infections ---- 

# primary infections 
# those with initial MFI values below a certain threshold level 

x <- x %>%
  left_join(
    inputs %>% mutate(person_id = as.integer(person_id)) %>%
      select(
        person_id, 
        dengns1_1_mfi_1, dengns1_1_mfi_2,
        dengns1_2_mfi_1, dengns1_2_mfi_2,
        dengns1_3_mfi_1, dengns1_3_mfi_2,
        dengns1_4_mfi_1, dengns1_4_mfi_2,
        zika_ns1_mfi_1, zika_ns1_mfi_2),
    by = c("id" = "person_id")
  )

# look at age against MFI 

# initial MFI
deng1_age_mfi0 <- var_gam_plot(x, "age_calc", "dengns1_1_mfi_1", "age", "initial mfi", "deng1: age vs initial mfi", 0, 100, 0, 60000)
deng2_age_mfi0 <- var_gam_plot(x, "age_calc", "dengns1_2_mfi_1", "age", "initial mfi", "deng2: age vs initial mfi", 0, 100, 0, 60000)
deng3_age_mfi0 <- var_gam_plot(x, "age_calc", "dengns1_3_mfi_1", "age", "initial mfi", "deng3: age vs initial mfi", 0, 100, 0, 60000)
deng4_age_mfi0 <- var_gam_plot(x, "age_calc", "dengns1_4_mfi_1", "age", "initial mfi", "deng4: age vs initial mfi", 0, 100, 0, 60000)
zika_age_mfi0 <- var_gam_plot(x, "age_calc", "zika_ns1_mfi_1", "age", "initial mfi", "zika: age vs initial mfi", 0, 100, 0, 60000)

# final MFI
deng1_age_mfi1 <- var_gam_plot(x, "age_calc", "dengns1_1_mfi_2", "age", "final mfi", "deng1: age vs final mfi", 0, 100, 0, 60000)
deng2_age_mfi1 <- var_gam_plot(x, "age_calc", "dengns1_2_mfi_2", "age", "final mfi", "deng2: age vs final mfi", 0, 100, 0, 60000)
deng3_age_mfi1 <- var_gam_plot(x, "age_calc", "dengns1_3_mfi_2", "age", "final mfi", "deng3: age vs final mfi", 0, 100, 0, 60000)
deng4_age_mfi1 <- var_gam_plot(x, "age_calc", "dengns1_4_mfi_2", "age", "final mfi", "deng4: age vs final mfi", 0, 100, 0, 60000)
zika_age_mfi1 <- var_gam_plot(x, "age_calc", "zika_ns1_mfi_2", "age", "final mfi", "zika: age vs final mfi", 0, 100, 0, 60000)

# change in MFI 
deng1_age_mfi <- var_gam_plot(x, "age_calc", "x_deng1", "age", "change mfi", "deng1: age vs change mfi", 0, 100, -60, 60)
deng2_age_mfi <- var_gam_plot(x, "age_calc", "x_deng2", "age", "change mfi", "deng2: age vs change mfi", 0, 100, -60, 60)
deng3_age_mfi <- var_gam_plot(x, "age_calc", "x_deng3", "age", "change mfi", "deng3: age vs change mfi", 0, 100, -60, 60)
deng4_age_mfi <- var_gam_plot(x, "age_calc", "x_deng4", "age", "change mfi", "deng4: age vs change mfi", 0, 100, -60, 60)
zika_age_mfi <- var_gam_plot(x, "age_calc", "x_zika", "age", "change mfi", "zika: age vs change mfi", 0, 100, -60, 60)

# plot with patchwork 
(deng1_age_mfi0 / deng2_age_mfi0 / deng3_age_mfi0 / deng4_age_mfi0 / zika_age_mfi0) |
(deng1_age_mfi1 / deng2_age_mfi1 / deng3_age_mfi1 / deng4_age_mfi1 / zika_age_mfi1) |
(deng1_age_mfi / deng2_age_mfi / deng3_age_mfi / deng4_age_mfi / zika_age_mfi) 

# age distribution of cases 
# infections before cases 

# 4. Spatial mapping ---- 

table(x$location, useNA = "ifany")
length(unique(x$location))

# 5. Reporting rates ---- 
# standard or age-adjusted 

# 6. Infection relative to suveillance data by age and location 

# 7. Probability of infection by initial MFI ----

# deng1 ----
ggplot(x, aes(dengns1_1_mfi_1/1000, p_deng1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Probability of infection",
    title = "deng1"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,50))

ggplot(x, aes(dengns1_1_mfi_1/1000, x_deng1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  labs(
    x = "Initial MFI",
    y = "Boost",
    title = "deng1"
  ) +
  theme_minimal() + 
  coord_cartesian(xlim = c(0,50))



