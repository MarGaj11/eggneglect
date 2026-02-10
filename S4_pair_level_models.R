# 0. Setup & Libraries
rm(list = ls())
setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

library(dplyr)      # Data manipulation
library(glmmTMB)    # Tweedie GLMMs
library(DHARMa)     # Residual diagnostics
library(sjPlot)     # Table summaries
library(boot)       # Bootstrap
library(rptR)

# 1. Load Data
data <- read.csv("./Julian_Date_COINr_NI.csv", sep=";", header=TRUE)
predictors <- read.csv("./SI_PS.csv", sep=";", header=TRUE)

# 2. Merge and Clean
merged <- data %>%
  select(Season, session, Nest, total_neglected_seconds, Incubation_Duration, Days_to_Hatch, NI) %>%
  left_join(
    predictors %>% select(Season, Nest, ID1, ID2, Similarity_Index, Pair_status, 
                          Similarity_Index_Head, Similarity_Index_Tarsus),
    by = c("Nest", "Season")
  ) %>%
  # Convert comma-decimals to numeric for all relevant columns
  mutate(across(c(total_neglected_seconds, Similarity_Index, Similarity_Index_Head, 
                  Similarity_Index_Tarsus, NI, Incubation_Duration, Days_to_Hatch), 
                ~ as.numeric(gsub(",", ".", .)))) %>%
  # Create Identifiers and Factors
  mutate(
    Nest_Season = as.factor(paste(Nest, "_", Season)),
    Pair_ID = as.factor(paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_")),
    Season = as.factor(Season),
    Nest = as.factor(Nest),
    Pair_status = as.factor(Pair_status),
    # Time variable (t=0 at day 31, t=1 at day 0)
    t = (31 - Days_to_Hatch) / 31
  )

# 3. GLMM Tweedie Models using total_neglected_seconds
model_tarsus    <- glmmTMB(NI ~ t + Similarity_Index_Tarsus + (1 | Pair_ID), 
                           data = merged, family = tweedie(link = "log"))

model_head      <- glmmTMB(NI ~ t + Similarity_Index_Head + (1 | Pair_ID), 
                           data = merged, family = tweedie(link = "log"))

model_wing      <- glmmTMB(NI ~ t + Similarity_Index + (1 | Pair_ID), 
                           data = merged, family = tweedie(link = "log"))

model_pair_bond <- glmmTMB(NI ~ t + Pair_status + (1 | Pair_ID), 
                           data = merged, family = tweedie(link = "log"))

# 4. View GLMM summaries
summary(model_tarsus)
summary(model_head)
summary(model_wing)
summary(model_pair_bond)

# 5. DHARMa Residual Diagnostics
res_wing <- simulateResiduals(model_wing)
plot(res_wing)

res_head <- simulateResiduals(model_head)
plot(res_head)

res_tarsus <- simulateResiduals(model_tarsus)
plot(res_tarsus)

res_pair_bond <- simulateResiduals(model_pair_bond)
plot(res_pair_bond)

# 6. Save Cleaned Data
write.csv2(merged, "NI_cleaned_data", row.names = FALSE)

# 7. Pair-status summary
combined_NI_statystyki <- merged %>%
  group_by(Pair_status) %>%
  summarise(
    unique_nests = n_distinct(Nest)
  )

# 8. Effect of season
model_season <- glmmTMB(NI ~ t + Season + (1 | Pair_ID), 
                        family = tweedie(link = "log"), 
                        data = merged)
summary(model_season)

res_pair_bond <- simulateResiduals(model_season)
plot(res_pair_bond)


rep_neglect <- rpt(log(NI + 0.001) ~ t + (1 | Pair_ID), 
                 grname = "Pair_ID", 
                 data = merged_filt, 
                 datatype = "Gaussian", 
                 nboot = 1000, npermut = 1000)

# 6. Results
print(rep_neglect)
plot(rep_neglect) # Distribution of the repeatability R


#analysis if neglect influences incubarion durarion

# 9. Pair-season-level analysis using total_neglected_seconds
merged_filtered <- merged %>%
  filter(Season %in% c(2019, 2021, 2025))

# Only include pair-seasons with all 3 sessions
complete_pair_seasons <- merged %>%
  group_by(Pair_ID, Season) %>%
  summarise(n_sessions = n_distinct(session), .groups = "drop") %>%
  filter(n_sessions == 3)

pair_season_summary <- merged %>%
  semi_join(complete_pair_seasons, by = c("Pair_ID", "Season")) %>%
  group_by(Pair_ID, Season) %>%
  summarise(
    mean_total_neglected_seconds = mean(total_neglected_seconds, na.rm = TRUE),
    total_total_neglected_seconds = sum(total_neglected_seconds, na.rm = TRUE),
    Incubation_Duration = first(Incubation_Duration),
    Days_to_Hatch = first(Days_to_Hatch),
    .groups = "drop"
  )

# 10. Bootstrap correlation
boot_correlation <- function(data, indices) {
  d <- data[indices, ]
  cor(d$total_total_neglected_seconds, d$Incubation_Duration, use = "complete.obs")
}

set.seed(123)
boot_res_cor <- boot(data = pair_season_summary,
                     statistic = boot_correlation,
                     R = 10000)

ci_perc_cor <- if(!is.na(sd(boot_res_cor$t))) boot.ci(boot_res_cor, type="perc")$percent[4:5] else NA
ci_bca_cor  <- if(!is.na(sd(boot_res_cor$t))) boot.ci(boot_res_cor, type="bca")$bca[4:5] else NA
p_val_cor <- mean(abs(boot_res_cor$t) >= abs(boot_res_cor$t0), na.rm = TRUE)

cor_results <- list(
  observed_correlation = boot_res_cor$t0,
  percentile_CI = ci_perc_cor,
  BCa_CI = ci_bca_cor,
  p_value = p_val_cor,
  bootstrap_object = boot_res_cor
)

print(cor_results)

