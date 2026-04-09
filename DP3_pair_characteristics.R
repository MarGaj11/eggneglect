rm(list = ls())

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

library(zoo)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(purrr)

# -------------------------------------------------------------------------
# This script calculates predictors used to explain variation in neglect
# between pairs. We start with similarity indices based on wing length,
# head length (THL), and tarsus length.
# -------------------------------------------------------------------------

# --------------------------
# 1. Load main capturing data
# --------------------------
capturing_main <- read_excel("C:/Users/Martyna/Dropbox/Hornsund metabase/Hornsund/SHO_LIAK_capturing.xlsx", guess_max = 10000)

# --------------------------
# 2. Load 2025-specific data
# --------------------------
capturing_2025 <- read_excel("C:/Users/Martyna/Dropbox/Liak2025/FieldData/CapturingDt2025.xlsx")

# --------------------------
# 3. Standardize columns and combine
# --------------------------
capturing_2025_clean <- capturing_2025 %>%
  select(Nest, RingNo, Sx, Season, Wing, THL, Tarsus, Age) 

capturing <- bind_rows(capturing_main, capturing_2025_clean)
write.csv2(capturing, "capturing_allseasons.csv")

# Filter out first-year birds and select relevant columns
capturing_filtered <- capturing %>%
  filter(Age != "1") %>%
  select(Nest, RingNo, R_leg, L_leg, Sx, Season, Wing, THL) %>%
  mutate(ID = paste(Nest, Season, sep = "_"))

# Fill missing wing measurements within individuals across years
capturing_filtered <- capturing_filtered %>%
  arrange(RingNo, Season) %>%   
  group_by(RingNo) %>%
  fill(Wing, .direction = "downup") %>%  
  ungroup() %>%
  distinct(Season, Nest, RingNo, .keep_all = TRUE)

# Aggregate measurements per pair
capturing_filtered <- capturing_filtered %>%
  group_by(Nest, Season) %>%
  distinct(RingNo, .keep_all = TRUE) %>%
  summarise(
    Pair_ID = paste(sort(RingNo), collapse = "_"),
    Wings   = paste(sort(Wing), collapse = "_"),
    .groups = "drop"
  ) %>%
  separate(Pair_ID, into = c("ID1", "ID2"), sep = "_", convert = TRUE) %>%
  separate(Wings, into = c("Wing1", "Wing2"), sep = "_", convert = TRUE)

# --------------------------
# 4. Summarize pairs across multiple seasons
# --------------------------
summarise_pairs <- function(df) {
  df %>%
    group_by(ID1, ID2) %>%
    mutate(
      Seasons_together = n_distinct(Season),
      MeanWing1 = if_else(Seasons_together > 1, mean(Wing1, na.rm = TRUE), Wing1),
      MeanWing2 = if_else(Seasons_together > 1, mean(Wing2, na.rm = TRUE), Wing2)
    ) %>%
    ungroup() %>%
    select(-Seasons_together) %>%
    arrange(Nest, Season)
}

result <- summarise_pairs(capturing_filtered)
result <- result %>%
  mutate(WingDiff = abs(MeanWing1 - MeanWing2))

# --------------------------
# 5. Bootstrap population to normalize difference
# --------------------------
capturing_non_na <- capturing %>% filter(!is.na(Wing))
set.seed(123)
n_boot <- 1000
boot_deltas <- replicate(n_boot, abs(diff(sample(capturing_non_na$Wing, 2, replace = TRUE))))
mean_boot_delta <- mean(boot_deltas)

# --------------------------
# 6. Calculate Similarity Index for wings
# --------------------------
result <- result %>%
  mutate(Similarity_Index = 1 - (WingDiff / mean_boot_delta),
         Similarity_Index = pmin(pmax(Similarity_Index, 0), 1)) %>%
  filter(!is.na(Similarity_Index))

# --------------------------
# 7. Determine pair status (New vs Old)
# --------------------------
classify_pairs <- function(df) {
  df %>%
    arrange(Nest, Season) %>% 
    group_by(Nest) %>%
    mutate(
      Pair_status = map_chr(row_number(), function(i) {
        current_IDs <- c(ID1[i], ID2[i])
        if (i == 1) return("Old")
        prev_IDs <- c(ID1[1:(i-1)], ID2[1:(i-1)]) %>% na.omit()
        if (sum(current_IDs %in% prev_IDs) == 1) return("New")
        return("Old")
      })
    ) %>%
    ungroup()
}

pair_status <- classify_pairs(capturing_filtered)

# --------------------------
# 8. Merge Similarity Index and pair status
# --------------------------
merged_df <- result %>%
  select(Nest, Season, ID1, ID2, Similarity_Index) %>%
  left_join(pair_status %>% select(Nest, Season, ID1, ID2, Pair_status),
            by = c("Nest", "Season", "ID1", "ID2"))

# --------------------------
# 9. Calculate Similarity Index for Head (THL)
# --------------------------
capturing_head <- capturing %>%
  filter(Age == ">2") %>%
  select(Nest, RingNo, Sx, Season, THL) %>%
  mutate(ID = paste(Nest, Season, sep = "_")) %>%
  arrange(RingNo, Season) %>%
  group_by(RingNo) %>%
  fill(THL, .direction = "downup") %>%
  ungroup() %>%
  group_by(Nest, Season) %>%
  distinct(RingNo, .keep_all = TRUE) %>%
  summarise(
    Pair_ID = paste(sort(RingNo), collapse = "_"),
    Heads   = paste(sort(THL), collapse = "_"),
    .groups = "drop"
  ) %>%
  separate(Pair_ID, into = c("ID1", "ID2"), sep = "_", convert = TRUE) %>%
  separate(Heads, into = c("Head1", "Head2"), sep = "_", convert = TRUE)

summarise_pairs_head <- function(df) {
  df %>%
    group_by(Nest, ID1, ID2) %>%
    summarise(
      Seasons_together = n_distinct(Season),
      MeanHead1 = ifelse(Seasons_together > 1, mean(Head1, na.rm = TRUE), first(Head1)),
      MeanHead2 = ifelse(Seasons_together > 1, mean(Head2, na.rm = TRUE), first(Head2)),
      .groups = "drop"
    )
}

result_2 <- summarise_pairs_head(capturing_head)
result_2 <- result_2 %>%
  mutate(HeadDiff = abs(MeanHead1 - MeanHead2))

# Bootstrap for head
capturing_head_non_na <- capturing %>% filter(!is.na(THL))
set.seed(123)
boot_deltas_head <- replicate(n_boot, abs(diff(sample(capturing_head_non_na$THL, 2, replace = TRUE))))
mean_boot_delta_head <- mean(boot_deltas_head)

result_2 <- result_2 %>%
  mutate(Similarity_Index_Head = 1 - (HeadDiff / mean_boot_delta_head),
         Similarity_Index_Head = pmin(pmax(Similarity_Index_Head, 0), 1)) %>%
  filter(!is.na(Similarity_Index_Head))

# --------------------------
# 10. Calculate Similarity Index for Tarsus
# --------------------------
capturing_tarsus <- capturing %>%
  filter(Age == ">2") %>%
  select(Nest, RingNo, Sx, Season, Tarsus) %>%
  mutate(ID = paste(Nest, Season, sep = "_")) %>%
  arrange(RingNo, Season) %>%
  group_by(RingNo) %>%
  fill(Tarsus, .direction = "downup") %>%
  ungroup() %>%
  group_by(Nest, Season) %>%
  distinct(RingNo, .keep_all = TRUE) %>%
  summarise(
    Pair_ID = paste(sort(RingNo), collapse = "_"),
    Tarsus = paste(sort(Tarsus), collapse = "_"),
    .groups = "drop"
  ) %>%
  separate(Pair_ID, into = c("ID1", "ID2"), sep = "_", convert = TRUE) %>%
  separate(Tarsus, into = c("Tarsus1", "Tarsus2"), sep = "_", convert = TRUE)

summarise_pairs_tarsus <- function(df) {
  df %>%
    group_by(Nest, ID1, ID2) %>%
    summarise(
      Seasons_together = n_distinct(Season),
      MeanTarsus1 = ifelse(Seasons_together > 1, mean(Tarsus1, na.rm = TRUE), first(Tarsus1)),
      MeanTarsus2 = ifelse(Seasons_together > 1, mean(Tarsus2, na.rm = TRUE), first(Tarsus2)),
      .groups = "drop"
    )
}

result_3 <- summarise_pairs_tarsus(capturing_tarsus)
result_3 <- result_3 %>%
  mutate(TarsusDiff = abs(MeanTarsus1 - MeanTarsus2))

# Bootstrap for Tarsus
capturing_tarsus_non_na <- capturing %>% filter(!is.na(Tarsus))
set.seed(123)
boot_deltas_tarsus <- replicate(n_boot, abs(diff(sample(capturing_tarsus_non_na$Tarsus, 2, replace = TRUE))))
mean_boot_delta_tarsus <- mean(boot_deltas_tarsus)

result_3 <- result_3 %>%
  mutate(Similarity_Index_Tarsus = 1 - (TarsusDiff / mean_boot_delta_tarsus),
         Similarity_Index_Tarsus = pmin(pmax(Similarity_Index_Tarsus, 0), 1)) %>%
  filter(!is.na(Similarity_Index_Tarsus))

# --------------------------
# 11. Merge all Similarity Indices
# --------------------------
all_together <- merged_df %>%
  left_join(result_2 %>% select(Nest, ID1, ID2, Similarity_Index_Head),
            by = c("Nest", "ID1", "ID2")) %>%
  left_join(result_3 %>% select(Nest, ID1, ID2, Similarity_Index_Tarsus),
            by = c("Nest", "ID1", "ID2"))

write.csv2(all_together, "SI_PS.csv")
