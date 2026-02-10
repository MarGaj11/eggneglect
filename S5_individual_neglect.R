rm(list = ls())

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(purrr)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

data <- readRDS("./DPA_d05_status_activity data.rds")

# --- Fix missing sex labels safely ---
data <- data %>%
  group_by(nest, season) %>%
  mutate(
    has_m = any(sx == "m", na.rm = TRUE),
    has_f = any(sx == "f", na.rm = TRUE),
    sx = case_when(
      is.na(sx) & has_f & !has_m ~ "m",
      is.na(sx) & has_m & !has_f ~ "f",
      TRUE ~ sx
    )
  ) %>%
  select(-has_m, -has_f) %>%
  ungroup()

# --- Basic formatting ---
data <- data %>%
  mutate(
    Season_Nest = paste(season, nest, sep = "_"),
    onset = as_datetime(onset),
    end   = as_datetime(end),
    session_start = as_datetime(session_start),
    session_end_trimmed = as_datetime(session_end_trimmed),
    duration_sec = as.numeric(end - onset, units = "secs"),
    basic_interval = case_when(
      grepl("nest_attendance", basic_interval) ~ "nest_attendance",
      grepl("absence", basic_interval)         ~ "absence",
      grepl("col_attendance", basic_interval)  ~ "col_attendance",
      TRUE ~ basic_interval
    )
  ) %>%
  filter(session %in% c("incubation1", "incubation2", "incubation3")) %>%
  arrange(Season_Nest, session, onset)

# --- Exclusions ---
cases_excl <- c("2019_K16_4", "2019_14_18K", "2020_14_2M", "2020_A20_1", "2020_D15_6",
                "2021_113_10", "2021_12_16", "2021_14_18K", "2021_14_6K", "2021_A20_3",
                "2021_D15_11", "2021_D17_8", "2021_W39B",
                "2025_K20_10_22cmplx", "2025_K21_37", "2025_W40",
                "2023_12_16", "2023_K20_8")

`%notin%` <- Negate(`%in%`)
data <- data %>% filter(Season_Nest %notin% cases_excl)

data <- data %>%
  filter(!(Season_Nest == "2025_K20_8" & session == "incubation3" |
             Season_Nest == "2021_101" & session == "incubation1" |
             Season_Nest == "2020_D16_10" & session == "incubation1"))

# --- Incubation sequence for fallback logic ---
incubation_sequence <- data %>%
  filter(basic_interval == "nest_attendance") %>%
  group_by(Season_Nest, session) %>%
  arrange(onset) %>%
  mutate(
    next_expected_sex = lead(sx),
    next_expected_sex = if_else(
      is.na(next_expected_sex),
      if_else(sx == "m", "f", "m"),
      next_expected_sex
    )
  ) %>%
  ungroup()

# --- GAP DETECTION (STARTS FROM FIRST ACTIVITY, NOT SESSION START) ---
detect_true_gaps <- function(df, session_end_time) {
  
  if (nrow(df) == 0) return(tibble())
  
  # 🔴 FIX: start from first recorded activity, not session_start
  start_time <- min(df$onset, na.rm = TRUE)
  
  if (is.character(session_end_time)) {
    end_time <- as.POSIXct(session_end_time, format = "%Y.%m.%d %H:%M:%S")
  } else {
    end_time <- as.POSIXct(session_end_time)
  }
  
  if (is.na(end_time) || is.na(start_time)) return(tibble())
  
  if (end_time <= start_time) {
    warning("Session end is before or equal to start - skipping")
    return(tibble())
  }
  
  atts <- df %>%
    filter(basic_interval == "nest_attendance") %>%
    arrange(onset)
  
  if (nrow(atts) == 0) {
    return(tibble(
      gap_start = start_time,
      gap_end   = end_time
    ))
  }
  
  gaps <- list()
  
  # First gap
  time_diff <- as.numeric(difftime(atts$onset[1], start_time, units = "secs"))
  if (time_diff > 1) {
    gap_end_time <- min(atts$onset[1], end_time)
    if (gap_end_time > start_time) {
      gaps[[length(gaps) + 1]] <- c(start_time, gap_end_time)
    }
  }
  
  current_end <- atts$end[1]
  
  for (i in 2:nrow(atts)) {
    if (atts$onset[i] > current_end) {
      gap_start_time <- max(current_end, start_time)
      gap_end_time <- min(atts$onset[i], end_time)
      if (gap_end_time > gap_start_time) {
        gaps[[length(gaps) + 1]] <- c(gap_start_time, gap_end_time)
      }
      current_end <- atts$end[i]
    } else {
      current_end <- max(current_end, atts$end[i])
    }
  }
  
  # Last gap
  if (current_end < end_time) {
    time_diff_end <- as.numeric(difftime(end_time, current_end, units = "secs"))
    if (time_diff_end > 1) {
      gaps[[length(gaps) + 1]] <- c(current_end, end_time)
    }
  }
  
  if (length(gaps) == 0) return(tibble())
  
  tibble(
    gap_start = as.POSIXct(vapply(gaps, `[`, FUN.VALUE = start_time, 1)),
    gap_end   = as.POSIXct(vapply(gaps, `[`, FUN.VALUE = end_time,   2))
  )
}

# --- Workload function ---
calculate_workload_up_to <- function(df, timepoint) {
  df %>%
    filter(basic_interval == "nest_attendance", onset < timepoint) %>%
    group_by(sx) %>%
    summarise(
      total_time = sum(pmin(end, timepoint) - onset),
      .groups = "drop"
    )
}

MAX_SHIFT_SEC <- 6 * 60 * 60   # 6 hours

# --- ALL NESTS TOGETHER (NO NORMAL / SPECIAL SPLIT) ---
neglect_table <- data %>%
  group_by(Season_Nest, nest, season, session) %>%
  group_modify(~ {
    session_end_trimmed <- .x$session_end_trimmed[1]
    
    gaps <- detect_true_gaps(.x, session_end_trimmed)
    if (nrow(gaps) == 0) return(tibble())
    
    session_seq <- incubation_sequence %>%
      filter(Season_Nest == .y$Season_Nest, session == .y$session)
    if (nrow(session_seq) == 0) return(tibble())
    
    gaps %>%
      rowwise() %>%
      mutate(
        anyone_on_nest = any(.x$basic_interval == "nest_attendance" &
                               .x$end > gap_start &
                               .x$onset < gap_end),
        
        responsible_sex = {
          last_att <- .x %>% 
            filter(basic_interval == "nest_attendance", end <= gap_start) %>% 
            arrange(end) %>% 
            dplyr::slice_tail(n = 1)
          
          next_att <- .x %>%
            filter(basic_interval == "nest_attendance", onset >= gap_end) %>%
            arrange(onset) %>%
            dplyr::slice(1)
          
          workload <- calculate_workload_up_to(.x, gap_start)
          
          time_m <- if (nrow(workload %>% filter(sx == "m")) > 0) {
            as.numeric(workload %>% filter(sx == "m") %>% pull(total_time), units = "secs")
          } else 0
          
          time_f <- if (nrow(workload %>% filter(sx == "f")) > 0) {
            as.numeric(workload %>% filter(sx == "f") %>% pull(total_time), units = "secs")
          } else 0
          
          gap_duration <- as.numeric(gap_end - gap_start, units = "secs")
          res <- NA_character_
          
          if (nrow(last_att) == 1 && nrow(next_att) == 1) {
            if (last_att$sx[1] == next_att$sx[1]) {
              if (gap_duration > MAX_SHIFT_SEC && time_m > 0 && time_f > 0) {
                res <- if (time_m > time_f) "f" else "m"
              } else {
                res <- last_att$sx[1]
              }
            } else {
              last_shift_duration <- as.numeric(last_att$end - last_att$onset, units = "secs")
              
              if (last_shift_duration > MAX_SHIFT_SEC) {
                res <- next_att$sx[1]
              } else if (time_m > 0 || time_f > 0) {
                total_time <- time_m + time_f
                ratio_m <- time_m / total_time
                time_diff <- abs(time_m - time_f)
                
                if (gap_duration >= 30 * 240 && (ratio_m > 0.6 || time_diff > MAX_SHIFT_SEC)) {
                  res <- if (time_m > time_f) "f" else "m"
                } else {
                  res <- next_att$sx[1]
                }
              } else {
                res <- next_att$sx[1]
              }
            }
          } else if (nrow(last_att) == 1) {
            last_shift_duration <- as.numeric(last_att$end - last_att$onset, units = "secs")
            
            if (last_shift_duration > MAX_SHIFT_SEC) {
              res <- if_else(last_att$sx[1] == "m", "f", "m")
            } else if (time_m > 0 || time_f > 0) {
              total_time <- time_m + time_f
              ratio_m <- time_m / total_time
              time_diff <- abs(time_m - time_f)
              
              if (gap_duration >= 30 * 240 && (ratio_m > 0.6 || time_diff > MAX_SHIFT_SEC)) {
                res <- if (time_m > time_f) "f" else "m"
              } else {
                res <- if_else(last_att$sx[1] == "m", "f", "m")
              }
            } else {
              res <- if_else(last_att$sx[1] == "m", "f", "m")
            }
          } else if (nrow(next_att) == 1) {
            res <- next_att$sx[1]
          } else {
            res <- session_seq$next_expected_sex[1]
          }
          
          if (is.na(res) && nrow(last_att) == 1) res <- if_else(last_att$sx[1] == "m", "f", "m")
          res
        },
        
        ringno = {
          bird_id <- .x %>% 
            filter(sx == responsible_sex) %>% 
            pull(ringno) %>% 
            unique() %>% 
            head(1)
          if (length(bird_id) == 0) NA_character_ else bird_id
        },
        
        duration_sec = as.numeric(gap_end - gap_start, units = "secs"),
        
        workload_info = {
          workload_diag <- calculate_workload_up_to(.x, gap_start)
          tm <- if (nrow(workload_diag %>% filter(sx == "m")) > 0) {
            as.numeric(workload_diag %>% filter(sx == "m") %>% pull(total_time), units = "secs")
          } else 0
          tf <- if (nrow(workload_diag %>% filter(sx == "f")) > 0) {
            as.numeric(workload_diag %>% filter(sx == "f") %>% pull(total_time), units = "secs")
          } else 0
          paste0("M:", round(tm/3600, 1), "h; F:", round(tf/3600, 1), "h")
        }
      ) %>%
      filter(!anyone_on_nest) %>% 
      ungroup()
  }) %>%
  ungroup()

# --- INDIVIDUAL NEGLECT SUMMARY ---
individual_neglect <- neglect_table %>%
  group_by(Season_Nest, nest, season, session, sx = responsible_sex, ringno) %>%
  summarise(
    Total_Complete_Neglect = sum(duration_sec),
    n_gaps = n(),
    mean_neglect_duration = mean(duration_sec, na.rm = TRUE),
    .groups = "drop"
  )

all_combinations <- data %>%
  filter(basic_interval == "nest_attendance") %>%
  distinct(Season_Nest, nest, season, session, sx, ringno)

# --- INCUBATION STATS ---
incubation_stats <- data %>%
  filter(basic_interval == "nest_attendance") %>%
  group_by(Season_Nest, nest, season, session, sx, ringno) %>%
  summarise(
    n_incubations = n(),
    mean_incubation_duration = mean(duration_sec, na.rm = TRUE),
    .groups = "drop"
  )

individual_neglect_full <- all_combinations %>%
  left_join(
    individual_neglect,
    by = c("Season_Nest", "nest", "season", "session", "sx", "ringno")
  ) %>%
  left_join(
    incubation_stats,
    by = c("Season_Nest", "nest", "season", "session", "sx", "ringno")
  ) %>%
  mutate(
    Total_Complete_Neglect = replace_na(Total_Complete_Neglect, 0),
    n_gaps = replace_na(n_gaps, 0),
    n_incubations = replace_na(n_incubations, 0),
    mean_incubation_duration = replace_na(mean_incubation_duration, 0)
  ) %>%
  arrange(Season_Nest, session, sx)

write.csv2(individual_neglect_full, "Ind_neglect.csv", row.names = FALSE)

# --- PLOTS ---
dir.create("Weryfikacja", showWarnings = FALSE)

plot_gantt <- function(id) {
  
  df <- data %>% filter(Season_Nest == id)
  gaps <- neglect_table %>% filter(Season_Nest == id)
  
  p <- ggplot() +
    geom_segment(
      data = df %>% filter(basic_interval == "absence"),
      aes(onset, sx, xend = end, yend = sx),
      color = "grey70", linewidth = 5
    ) +
    geom_segment(
      data = gaps,
      aes(gap_start, responsible_sex, xend = gap_end, yend = responsible_sex),
      color = "red", linewidth = 5
    ) +
    scale_y_discrete(limits = c("m", "f")) +
    geom_segment(
      data = df %>% filter(basic_interval == "nest_attendance"),
      aes(onset, sx, xend = end, yend = sx),
      color = "darkgreen", linewidth = 5
    ) +
    facet_wrap(~ session, scales = "free_x", ncol = 1) +
    theme_minimal() +
    labs(title = paste(id, "- Red = neglect periods"),
         y = "Sex", x = "Time")
  
  ggsave(
    filename = paste0("Weryfikacja/WIO_", id, ".jpg"),
    plot = p,
    width = 14,
    height = 8
  )
}

walk(unique(data$Season_Nest), safely(plot_gantt))
