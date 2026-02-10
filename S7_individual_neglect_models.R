rm(list = ls())

library(rptR)
library(lme4)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(performance)
library(glmmTMB)

# 1. Load data
setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")
data <- read.csv("./COINr_Individual_NI.csv", sep=";", header=TRUE)
hatch <- read.csv("./Julian_Date_COINr_NI.csv", sep=";", header=TRUE)

# 2. Prepare and merge incubation data
hatch <- hatch %>%
  mutate(Season_Nest = paste(Season, Nest, sep = "_"))

data <- data %>%
  mutate(Season_Nest = paste(season, nest, sep = "_")) %>%
  select(Season_Nest, nest, season, ringno, session, sx, NI, mean_incubation_duration) %>%
  left_join(
    hatch %>% select(Season_Nest, session, Days_to_Hatch, Incubation_Duration),
    by = c("Season_Nest", "session")
  )

# 3. Add partner information
partner_data <- data %>%
  select(season, nest, session, 
         Partner_ringno = ringno, 
         Partner_Sx = sx, 
         Partner_NI = NI,
         Partner_Incubation = mean_incubation_duration)

data <- data %>%
  left_join(partner_data, by = c("season", "nest", "session")) %>%
  filter(sx != Partner_Sx) %>%
  select(-Partner_Sx)

# 4. Data cleaning
fix_num <- function(x) as.numeric(gsub(",", ".", x))

data <- data %>%
  mutate(
    NI = fix_num(NI),
    t = (31 - as.numeric(Days_to_Hatch)) / 31,
    NI_plus = NI + 0.01,
    session = as.factor(session),
    season = as.factor(season),
    ringno = as.factor(ringno),
    Partner_ringno = as.factor(Partner_ringno)
  )

# 5. Split by sex
data_female <- data %>% filter(sx == "f")
data_male   <- data %>% filter(sx == "m")

# 6. Model comparison
m_lm   <- lm(log(NI + 0.001) ~ t + season + sx, data = data)
m_lmm  <- lmer(log(NI + 0.001) ~ t + season + sx + (1 | ringno), data = data)
m_glmm <- glmer(NI_plus ~ t + season + sx + (1 | ringno), 
                family = Gamma(link = "log"), data = data)
mod_glmm <- glmmTMB(NI ~ t + season + sx + (1 | ringno), 
                    family = tweedie(link = "log"), 
                    data = data)

comparison <- compare_performance(m_lm, m_lmm, m_glmm, mod_glmm, rank = TRUE)
print(comparison)

# 7. Sex differences model (main test)
model_sex_diff <- glmmTMB(
  NI ~ t + season + sx + (1 | ringno), 
  family = tweedie(link = "log"), 
  data = data
)

summary(model_sex_diff)

res_glmm <- simulateResiduals(model_sex_diff)
plot(res_glmm)

# 8. Sex-specific models
model_neglect_female <- glmmTMB(
  NI ~ t + season + (1 | ringno), 
  family = tweedie(link = "log"), 
  data = data_female
)

model_neglect_male <- glmmTMB(
  NI ~ t + season + (1 | ringno), 
  family = tweedie(link = "log"), 
  data = data_male
)

summary(model_neglect_female)
summary(model_neglect_male)

plot(simulateResiduals(model_neglect_female))
plot(simulateResiduals(model_neglect_male))

# 9. Repeatability
rep_male <- rpt(
  log(NI + 0.001) ~ t + season + (1 | ringno), 
  grname = "ringno", 
  data = data_male, 
  datatype = "Gaussian", 
  nboot = 1000, npermut = 1000
)

rep_female <- rpt(
  log(NI + 0.001) ~ t + season + (1 | ringno), 
  grname = "ringno", 
  data = data_female, 
  datatype = "Gaussian", 
  nboot = 1000, npermut = 1000
)

print(rep_male)
print(rep_female)

# 10. Save processed data
write.csv2(data, "Ind_neglect_Julian_Date.csv", row.names = FALSE)
