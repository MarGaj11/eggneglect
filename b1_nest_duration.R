rm(list = ls())
library(dplyr)
library(lubridate)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# Load the data
data <- readRDS("./DPA_d05_status_activity data.rds")

# Filter only nest attendance moments for both parents
nest_data <- data %>%
  filter(
    basic_interval %in% c("nest_attendance"),
    session %in% c("incubation1", "incubation2", "incubation3")
  )

# Function to compute overlap intervals (both parents in nest)
get_nest_overlaps <- function(nest_data) {
  
  # Ensure both sexes are present
  if (!all(c("m", "f") %in% unique(nest_data$sx))) {
    return(NULL)
  }
  
  nest_data <- nest_data %>%
    mutate(
      onset_dt = ymd_hms(gsub("\\.", "-", onset), quiet = TRUE),
      end_dt   = ymd_hms(gsub("\\.", "-", end_trimmed), quiet = TRUE)
    ) %>%
    filter(!is.na(onset_dt), !is.na(end_dt), end_dt > onset_dt)
  
  male_intervals <- nest_data %>% filter(sx == "m")
  female_intervals <- nest_data %>% filter(sx == "f")
  
  if (nrow(male_intervals) == 0 || nrow(female_intervals) == 0) {
    return(NULL)
  }
  
  overlaps <- list()
  
  for (i in seq_len(nrow(male_intervals))) {
    for (j in seq_len(nrow(female_intervals))) {
      
      start_overlap <- max(male_intervals$onset_dt[i],
                           female_intervals$onset_dt[j])
      end_overlap <- min(male_intervals$end_dt[i],
                         female_intervals$end_dt[j])
      
      # Keep only positive-duration overlaps
      if (!is.na(start_overlap) && !is.na(end_overlap) && end_overlap > start_overlap) {
        overlaps[[length(overlaps) + 1]] <- data.frame(
          start = start_overlap,
          end   = end_overlap
        )
      }
    }
  }
  
  if (length(overlaps) == 0) return(NULL)
  
  overlap_df <- bind_rows(overlaps) %>%
    filter(!is.na(start), !is.na(end), end > start) %>%
    arrange(start)
  
  if (nrow(overlap_df) == 0) return(NULL)
  
  return(overlap_df)
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
      # Extend the last event
      merged[[length(merged)]]$end <- max(last_event$end, current_event$end)
    } else {
      # Start a new event
      merged[[length(merged) + 1]] <- current_event
    }
  }
  
  merged_df <- bind_rows(merged) %>%
    filter(!is.na(start), !is.na(end), end > start)
  
  if (nrow(merged_df) == 0) return(NULL)
  
  return(merged_df)
}

# Apply to each nest/session/group
nest_together_summary <- nest_data %>%
  group_by(season, session, nest, group) %>%
  group_modify(~{
    
    overlap_df <- get_nest_overlaps(.x)
    merged_df  <- merge_events(overlap_df, gap_threshold = 60)
    
    if (is.null(merged_df)) {
      return(tibble(
        total_both_nest_seconds = 0,
        number_of_both_nest_events = 0,
        average_both_nest_duration = 0
      ))
    }
    
    # Compute positive-duration events
    event_durations <- as.numeric(difftime(merged_df$end, merged_df$start, units = "secs"))
    event_durations <- event_durations[event_durations > 0]
    
    n_events <- length(event_durations)
    total_seconds <- sum(event_durations)
    avg_duration <- ifelse(n_events > 0, mean(event_durations), 0)
    
    tibble(
      total_both_nest_seconds = total_seconds,
      number_of_both_nest_events = n_events,
      average_both_nest_duration = avg_duration
    )
  }) %>%
  ungroup()

# Save results
write.csv2(nest_together_summary, "nest_both_parents_summary.csv", row.names = FALSE)
