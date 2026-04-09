# PURPOSE: measure repeatability of being neglectful
library(rptR)
S <- readRDS("EDA_sex_neglect_wide.RDS")
# ------------------------------
# Overall repeatability in sum_neglect, accordingly in n_gaps - is neglecting or not (Binary)
# ------------------------------
S$is_neglecting <- ifelse(S$sum_neglect_sec > 0, 1, 0)

rep_binary <- rpt(
  is_neglecting ~ session + season + (1 | ringno),
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
  is_neglecting ~ session + season + (1 | ringno),
  grname = "ringno",
  data = S_male,
  datatype = "Binary", # Rozkład dwumianowy (0/1)
  nboot = 1000, npermut = 1000
)

rep_binary_female <- rpt(
  is_neglecting ~ session + season + (1 | ringno),
  grname = "ringno",
  data = S_female,
  datatype = "Binary",
  nboot = 1000,
  npermut = 1000
)

print(rep_binary_male)
print(rep_binary_female)



# 1. Load and Clean Data
# -------------------------------------------------------------------------
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")

E_rep <- E %>%
  # Remove missing data first
  filter(!is.na(Incubation_Duration)) %>%
  # Create a unique Nest ID regardless of ID order
  rowwise() %>%
  mutate(
    nest_id = paste(sort(c(as.character(ID1), as.character(ID2))), collapse = "_")
  ) %>%
  ungroup() %>%
  # Keep only one record per nest per season
  group_by(nest_id, season) %>%
  summarize(inc_dist = first(as.numeric(Incubation_Duration)), .groups = "drop") %>%
  # Filter for "loyal" nests appearing in at least 2 seasons
  group_by(nest_id) %>%
  filter(n() >= 2) %>% 
  ungroup()

# Quick validation
cat("Number of unique pairs with multi-season data:", length(unique(E_rep$nest_id)), "\n")

# 2. Run Repeatability Analysis
# -------------------------------------------------------------------------
# Using Gaussian because, while discrete, the data is centered around a mean (29)
rep_seasonal <- rpt(
  inc_dist ~ (1 | nest_id), 
  grname   = "nest_id", 
  data     = E_rep, 
  datatype = "Gaussian", 
  nboot    = 1000, 
  npermut  = 1000
)

# 3. Results
# -------------------------------------------------------------------------
print(rep_seasonal)
plot(rep_seasonal)
