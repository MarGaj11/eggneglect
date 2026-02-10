rm(list = ls())
library(dplyr)
library(lubridate)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# Load the data
data <- readRDS("./DPA_d05_status_activity data.rds")

# Filter relevant intervals
neglect_data <- data %>%
  filter(
    basic_interval %in% c("absence", "end_absence", "start_absence", "col_attendance"),
    session %in% c("incubation1", "incubation2", "incubation3")
  )

# Function to compute overlap intervals
get_neglect_overlaps <- function(nest_data) {
  if (!all(c("m", "f") %in% unique(nest_data$sx))) return(NULL)
  
  nest_data <- nest_data %>%
    mutate(
      onset_dt = ymd_hms(gsub("\\.", "-", onset), quiet = TRUE),
      end_dt   = ymd_hms(gsub("\\.", "-", end_trimmed), quiet = TRUE)
    ) %>%
    filter(!is.na(onset_dt), !is.na(end_dt), end_dt > onset_dt)
  
  male_intervals <- nest_data %>% filter(sx == "m")
  female_intervals <- nest_data %>% filter(sx == "f")
  
  if (nrow(male_intervals) == 0 || nrow(female_intervals) == 0) return(NULL)
  
  overlaps <- list()
  for (i in seq_len(nrow(male_intervals))) {
    for (j in seq_len(nrow(female_intervals))) {
      start_overlap <- max(male_intervals$onset_dt[i], female_intervals$onset_dt[j])
      end_overlap <- min(male_intervals$end_dt[i], female_intervals$end_dt[j])
      if (!is.na(start_overlap) && !is.na(end_overlap) && end_overlap > start_overlap) {
        overlaps[[length(overlaps) + 1]] <- data.frame(start = start_overlap, end = end_overlap)
      }
    }
  }
  
  if (length(overlaps) == 0) return(NULL)
  bind_rows(overlaps) %>% arrange(start)
}

# Function to merge overlap intervals into distinct events
merge_events <- function(overlap_df, gap_threshold = 60) {
  if (is.null(overlap_df) || nrow(overlap_df) == 0) return(NULL)
  
  merged <- list(overlap_df[1, ])
  for (i in 2:nrow(overlap_df)) {
    last_event <- merged[[length(merged)]]
    current_event <- overlap_df[i, ]
    gap <- as.numeric(difftime(current_event$start, last_event$end, units = "secs"))
    if (!is.na(gap) && gap <= gap_threshold) {
      merged[[length(merged)]]$end <- max(last_event$end, current_event$end)
    } else {
      merged[[length(merged) + 1]] <- current_event
    }
  }
  bind_rows(merged)
}

# Apply to each nest/session/group
neglect_summary <- neglect_data %>%
  group_by(season, session, nest, group) %>%
  group_modify(~{
    overlap_df <- get_neglect_overlaps(.x)
    merged_df  <- merge_events(overlap_df, gap_threshold = 60)
    
    if (is.null(merged_df) || nrow(merged_df) == 0) {
      total_sec <- 0
      n_events <- 0
      avg_sec <- 0
    } else {
      event_durations <- as.numeric(difftime(merged_df$end, merged_df$start, units = "secs"))
      event_durations <- event_durations[!is.na(event_durations) & event_durations > 0]
      total_sec <- sum(event_durations, na.rm = TRUE)
      n_events <- length(event_durations)
      avg_sec <- ifelse(n_events > 0, mean(event_durations, na.rm = TRUE), 0)
    }
    
    tibble(
      total_neglected_seconds = total_sec,
      number_of_neglect_events = n_events,
      average_neglect_duration = avg_sec
    )
  }) %>%
  ungroup()

# Save results
write.csv2(neglect_summary, "neglect_summary.csv", row.names = FALSE)
