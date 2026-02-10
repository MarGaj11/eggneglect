rm(list = ls())

library(dplyr)
library(COINr)
library(ggplot2)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# ============================
# Load combined neglect summary
# ============================
neglect_summary <- read.csv2("./neglect_summary.csv")  
# Columns: season, session, nest, group, total_neglected_seconds, number_of_neglect_events, average_neglect_duration

# ============================
# Exclude suspicious nests
# ============================
neglect_summary <- neglect_summary %>% 
  mutate(season_nest = paste(season, "_", nest, sep = ""))

cases_excl <- c("2019_K16_4", "2019_14_18K", "2020_14_2M", "2020_A20_1",  "2020_D15_6",
                "2021_113_10", "2021_12_16", "2021_14_18K", "2021_14_6K", "2021_A20_3", 
                "2021_D15_11", "2021_D17_8", "2021_W39B", 
                "2025_K20_10_22cmplx", "2025_K21_37", "2025_W4") 

`%notin%` <- Negate(`%in%`)
neglect_summary <- neglect_summary %>% filter(season_nest %notin% cases_excl)

# ============================
# Create unique codes for COINr
# ============================
neglect_summary <- neglect_summary %>%
  mutate(
    uCode = make.unique(paste0(nest, "_", session, "_", substr(season, 3, 4))),
    uName = uCode
  )

# ============================
# Prepare COINr input data
# ============================
iData <- neglect_summary[, c("uCode", "uName", "total_neglected_seconds", "average_neglect_duration")]

# Ensure numeric columns
iData$total_neglected_seconds <- as.numeric(iData$total_neglected_seconds)
iData$average_neglect_duration <- as.numeric(iData$average_neglect_duration)

# ============================
# Define indicator metadata
# ============================
iMeta <- data.frame(
  iCode = c("total_neglected_seconds", "average_neglect_duration", "NI"),
  iName = c("Total neglect duration", "Average neglect duration", "Neglect Index"),
  Direction = c(1, 1, 1),
  Level = c(1, 1, 2),
  Weight = c(1, 1, 1),  # all weights = 1
  Type = c("Indicator", "Indicator", "Aggregate"),
  Parent = c("NI", "NI", NA)
)

# ============================
# Create COINr object
# ============================
ESI <- new_coin(iData, iMeta)
plot_framework(ESI)

# ============================
# Data treatment and normalization
# ============================
ESI <- qTreat(ESI, dset = "Raw", winmax = 5)
ESI <- qNormalise(ESI, dset = "Treated")

# ============================
# Aggregate using equal weights
# ============================
ESI <- Aggregate(ESI, dset = "Normalised", f_ag = "a_amean", use_weights = TRUE)

# ============================
# Get COINr results
# ============================
ESI_df <- get_results(ESI, dset = "Aggregated")

# ============================
# Merge COINr NI with original data
# ============================
combined_NI <- neglect_summary %>%
  inner_join(ESI_df, by = "uCode")

# Keep total_neglected_seconds for future models
# Drop unnecessary columns but keep original and aggregated columns
combined_NI <- combined_NI %>%
  select(season, session, nest, group, total_neglected_seconds, number_of_neglect_events,
         average_neglect_duration, everything(), -c(uName, season_nest))

# ============================
# Save final dataset
# ============================
write.csv2(combined_NI, "COINr_NI_neglect_only.csv")

# ============================
# Summary statistics
# ============================
combined_NI_nests <- combined_NI %>%
  summarise(unique_nests = n_distinct(nest))

total_sessions <- combined_NI %>%
  group_by(season) %>%
  summarise(number_sessions = n_distinct(nest, session))

nest_history <- combined_NI %>%
  group_by(nest) %>%
  summarise(number_seasons = n_distinct(season))

nests <- data.frame(
  Criterium = c("At least 2 seasons", "At least 3 seasons", "At least 4 seasons", "At least 5 seasons"),
  Number_of_nests = c(
    sum(nest_history$number_seasons >= 2),
    sum(nest_history$number_seasons >= 3),
    sum(nest_history$number_seasons >= 4),
    sum(nest_history$number_seasons >= 5)
  )
)

print(nests)
