rm(list = ls())

library(dplyr)
library(lubridate)
library(COINr)
library(ggplot2)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# Load the data
Ind <- read.csv2("./Ind_neglect.csv")

# Create unique codes
Ind$uCode <- with(Ind, make.unique(paste0(as.character(Ind[, 3]), "_", session, "_", nest, "_", sx)))
Ind$uName <- Ind$uCode

# Prepare data for COINr
data <- Ind[, c(13,12,11,10,9,8,7)]

iData <- data[, c(
  "uCode", "uName",
  "mean_neglect_duration", 
  "Total_Complete_Neglect"
)]

iMeta <- data.frame(
  iCode = c("mean_neglect_duration", "Total_Complete_Neglect", "NI"),
  iName = c("Average neglect duration", "Total neglect", "NI"),
  Direction = c(1, 1, 1),
  Level = c(1, 1, 2),
  Weight = c(1, 1, 1),
  Type = c("Indicator", "Indicator", "Aggregate"),
  Parent = c("NI", "NI", NA)
)

# Equal weights for the two indicators
iMeta$Weight[iMeta$Level == 1] <- c(1, 1)

# Create COIN object
ESI <- new_coin(iData, iMeta)

# Data treatment
ESI <- qTreat(ESI, dset = "Raw", winmax = 5)
ESI <- qNormalise(ESI, dset = "Treated")

# Aggregate
ESI <- Aggregate(ESI, dset = "Normalised", f_ag = "a_amean", use_weights = TRUE)

# Extract results
ESI_df <- get_results(ESI, dset = "Aggregated")

# Merge back to original data
Ind_NI <- Ind %>%
  inner_join(ESI_df, by = "uCode")

# Save output
write.csv2(Ind_NI, "COINr_Individual_NI.csv", row.names = FALSE)
