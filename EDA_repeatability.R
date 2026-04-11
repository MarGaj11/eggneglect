# PURPOSE: measure repeatability of being neglectful
library(rptR)
library(tidyverse)
S <- readRDS("EDA_sex_neglect_wide.RDS")
# ------------------------------
# Overall repeatability in sum_neglect, accordingly in n_gaps - is neglecting or not (Binary)
# ------------------------------
S$is_neglecting <- ifelse(S$sum_neglect_sec > 0, 1, 0)

rep_binary <- rpt(
  is_neglecting ~ session + (1 | ringno),
  grname = "ringno",
  data = S,
  datatype = "Binary", # Rozkład dwumianowy (0/1)
  nboot = 1000, npermut = 1000
)

print(rep_binary)

# ------------------------------
# Repeatability by sex - is neglectful or not - binary
# ------------------------------
S_male <- S %>% filter(sx == "m")
S_female <- S %>% filter(sx == "f")

rep_binary_male <- rpt(
  is_neglecting ~ session + (1 | ringno),
  grname = "ringno",
  data = S_male,
  datatype = "Binary", # Rozkład dwumianowy (0/1)
  nboot = 1000, npermut = 1000
)

rep_binary_female <- rpt(
  is_neglecting ~ session + (1 | ringno),
  grname = "ringno",
  data = S_female,
  datatype = "Binary",
  nboot = 1000,
  npermut = 1000
)

print(rep_binary_male)
print(rep_binary_female)

