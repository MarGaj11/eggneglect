rm(list = ls())
library(dplyr)
library(lubridate)
library(tidyr)

setwd("C:/Users/user/Dropbox/Negotiation")

data <- readRDS("./DPA_d05_status_activity data.rds")

head(data)

# 1. Helper function to merge a single bird's overlapping or back-to-back intervals
merge_bird_intervals <- function(df) {
  # 1. Sprawdzenie podstawowe
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  # 2. Oczyszczenie i upewnienie się, że mamy format czasu
  df <- df %>% 
    filter(!is.na(start), !is.na(end)) %>% 
    arrange(start)
  
  if(nrow(df) == 0) return(NULL)
  
  # 3. Logika scalania - używamy uproszczonego podejścia
  merged_df <- df[1, ]
  
  if(nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      # Pobieramy ostatni dodany wiersz do porównania
      last_idx <- nrow(merged_df)
      
      # WYMUSZAMY porównanie jako POSIXct (to zapobiega błędowi TRUE/FALSE)
      if (as.numeric(df$start[i]) <= as.numeric(merged_df$end[last_idx])) {
        merged_df$end[last_idx] <- max(merged_df$end[last_idx], df$end[i])
      } else {
        merged_df <- rbind(merged_df, df[i, ])
      }
    }
  }
  
  return(merged_df)
}

# 2. Function to compute overlaps between Male and Female absence
get_neglect_overlaps <- function(nest_data) {
  # 1. Przygotowanie danych
  nest_data <- nest_data %>%
    mutate(
      onset_dt = ymd_hms(gsub("\\.", "-", onset), quiet = TRUE),
      end_dt   = ymd_hms(gsub("\\.", "-", end_trimmed), quiet = TRUE)
    ) %>%
    filter(!is.na(onset_dt), !is.na(end_dt), end_dt > onset_dt) %>%
    select(-any_of(c("start", "end"))) %>% 
    rename(start = onset_dt, end = end_dt)
  
  # 2. Scalanie interwałów dla płci
  male_merged   <- merge_bird_intervals(nest_data %>% filter(sx == "m"))
  female_merged <- merge_bird_intervals(nest_data %>% filter(sx == "f"))
  
  # KLUCZOWA POPRAWKA: bezpieczne sprawdzenie NULL lub braku wierszy
  if (is.null(male_merged) || is.null(female_merged)) return(NULL)
  if (nrow(male_merged) == 0 || nrow(female_merged) == 0) return(NULL)
  
  overlaps <- list()
  for (i in 1:nrow(male_merged)) {
    for (j in 1:nrow(female_merged)) {
      # Wymuszamy format numeryczny do porównania (bezpieczeństwo przed NA)
      s_ov <- max(as.numeric(male_merged$start[i]), as.numeric(female_merged$start[j]))
      e_ov <- min(as.numeric(male_merged$end[i]), as.numeric(female_merged$end[j]))
      
      if (s_ov < e_ov) {
        overlaps[[length(overlaps) + 1]] <- data.frame(
          start = as_datetime(s_ov), 
          end = as_datetime(e_ov)
        )
      }
    }
  }
  
  if (length(overlaps) == 0) return(NULL)
  bind_rows(overlaps) %>% arrange(start)
}

# 3. Main processing
neglect_summary <- data %>%
  filter(
    basic_interval %in% c("absence", "end_absence", "start_absence", "col_attendance"),
    session %in% c("incubation1", "incubation2", "incubation3")
  ) %>%
  group_by(season, session, nest, group) %>%
  group_modify(~{
    ov <- get_neglect_overlaps(.x)
    
    # We use a 1-second gap threshold because video data is continuous
    # Only merge if they are literally touching
    merged_df <- if(!is.null(ov)) merge_bird_intervals(ov) else NULL
    
    if (is.null(merged_df) || nrow(merged_df) == 0) {
      res <- tibble(total_neglected_seconds = 0, number_of_neglect_events = 0, average_neglect_duration = 0)
    } else {
      durs <- as.numeric(difftime(merged_df$end, merged_df$start, units = "secs"))
      res <- tibble(
        total_neglected_seconds = sum(durs),
        number_of_neglect_events = nrow(merged_df),
        average_neglect_duration = mean(durs)
      )
    }
    res
  }) %>%
  ungroup()

write.csv2(neglect_summary, "neglect_summary_robust.csv", row.names = FALSE)
