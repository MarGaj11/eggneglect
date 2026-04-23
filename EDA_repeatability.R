# PURPOSE: measure repeatability of being neglectful across incubation stage
library(rptR)
library(tidyverse)

S <- readRDS("EDA_sex_neglect_wide.RDS")

# ------------------------------
# Prepare dph_c exactly as in your modelling script
# ------------------------------
S <- S %>%
  filter(!is.na(day_prior_hatch),
         day_prior_hatch >= 2,
         day_prior_hatch <= 28) %>%
  mutate(
    dph_c = day_prior_hatch - mean(day_prior_hatch, na.rm = TRUE),
    is_neglecting = ifelse(sum_neglect_sec > 0, 1, 0)
  )

# ------------------------------
# Overall repeatability
# Binary: neglecting or not
# ------------------------------
rep_binary <- rpt(
  is_neglecting ~ dph_c + (1 | ringno),
  grname = "ringno",
  data = S,
  datatype = "Binary",
  nboot = 1000,
  npermut = 1000
)

print(rep_binary)

# ------------------------------
# Repeatability by sex
# ------------------------------
S_male <- S %>% filter(sx == "m")
S_female <- S %>% filter(sx == "f")

rep_binary_male <- rpt(
  is_neglecting ~ dph_c + (1 | ringno),
  grname = "ringno",
  data = S_male,
  datatype = "Binary",
  nboot = 1000,
  npermut = 1000
)

rep_binary_female <- rpt(
  is_neglecting ~ dph_c + (1 | ringno),
  grname = "ringno",
  data = S_female,
  datatype = "Binary",
  nboot = 1000,
  npermut = 1000
)

print(rep_binary_male)
print(rep_binary_female)