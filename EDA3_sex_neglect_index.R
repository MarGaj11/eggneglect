# PURPOSE: Calculate sex-separated Egg Neglection Index
# Output:
#   - EDA_sex_neglect_wide.RDS   : one row per season x session x nest x sex (+ ringno)
#   - VERIFICATION_SEX/          : Gantt plots with m/f rows, gaps coloured by responsible sex
rm(list = ls())

act_dt <- readRDS("DPA_d05_status_activity data.rds")

library(readxl)
library(tidyverse)
library(lubridate)
library(purrr)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# ── 1. LOAD & PREPARE DATA ────────────────────────────────────────────────────

# Load raw data and filter to attendance intervals only (same as Script 1)
act_dt <- readRDS("DPA_d05_status_activity data.rds")

act_dt_nest <- act_dt %>%
  mutate(loop_unit = paste0(season, "_", session, "_", nest)) %>%
  filter(
    basic_interval %in% c("start_nest_attendance", "nest_attendance", "end_nest_attendance"),
    session %in% c("incubation1", "incubation2", "incubation3")
  )

# ── 2. EXCLUSIONS (identical to Script 1) ────────────────────────────────────

act_dt_nest <- act_dt_nest %>%
  mutate(season_nest = paste(season, "_", nest, sep = ""))

cases_excl <- c(
  "2019_K16_4", "2019_14_18K", "2020_14_2M", "2020_A20_1", "2020_D15_6",
  "2021_113_10", "2021_12_16", "2021_14_18K", "2021_14_6K", "2021_A20_3",
  "2021_D15_11", "2021_D17_8", "2021_W39B", "2023_12_16",
  "2025_K20_10_22cmplx", "2025_K21_37", "2025_W4", "2025_K19_15_22cmplx",
  "2025_W40", "2025_A20_4"
)

`%notin%` <- Negate(`%in%`)
act_dt_nest <- act_dt_nest %>% filter(season_nest %notin% cases_excl)

act_dt_nest <- act_dt_nest %>%
  filter(
    !(season_nest == "2025_D15_11" & session == "incubation1"),
    !(season_nest == "2025_14_10M" & session == "incubation1"),
    !(season_nest == "2025_K20_8"  & session == "incubation3"),
    !(season_nest == "2023_12_16"  & session == "incubation2")
  )

# ── 3. FIX MISSING SEX LABELS ────────────────────────────────────────────────

act_dt_nest <- act_dt_nest %>%
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

# ── 4. HELPER CONSTANTS ───────────────────────────────────────────────────────

MAX_SHIFT_SEC <- 6 * 60 * 60  # 6 hours

# Incubation sequence — used as fallback for sex attribution
incubation_sequence <- act_dt_nest %>%
  group_by(season, session, nest) %>%
  arrange(as_datetime(onset)) %>%
  mutate(
    next_expected_sex = lead(sx),
    next_expected_sex = if_else(
      is.na(next_expected_sex),
      if_else(sx == "m", "f", "m"),
      next_expected_sex
    )
  ) %>%
  ungroup()

# ── 5. CUMULATIVE WORKLOAD HELPER ─────────────────────────────────────────────

calculate_workload_up_to <- function(df, timepoint) {
  df %>%
    filter(as_datetime(onset) < timepoint) %>%
    group_by(sx) %>%
    summarise(
      total_time = sum(
        as.numeric(
          pmin(as_datetime(end_trimmed), timepoint) - as_datetime(onset),
          units = "secs"
        )
      ),
      .groups = "drop"
    )
}

# ── 6. GAP DETECTION (identical logic to Script 1) ───────────────────────────
# Uses session_start / session_end_trimmed as window boundaries.
# Sorts and merges overlapping attendance intervals before finding gaps.

loop_units <- unique(act_dt_nest$loop_unit)
gaps_list  <- list()

for (j in seq_along(loop_units)) {
  
  act_dt_temp <- act_dt_nest %>% filter(loop_unit == loop_units[j])
  
  # Session window boundaries — same as Script 1
  wide_start <- as_datetime(min(act_dt_temp$session_start))
  wide_end   <- as_datetime(max(act_dt_temp$session_end_trimmed))
  
  # Build attendance df with datetimes, sort by start
  df <- data.frame(
    start = as_datetime(act_dt_temp$onset),
    end   = as_datetime(act_dt_temp$end_trimmed)
  ) %>% arrange(start)
  
  # Merge overlapping intervals — identical to Script 1
  merged <- df[1, ]
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      if (df$start[i] <= merged$end[nrow(merged)]) {
        merged$end[nrow(merged)] <- max(merged$end[nrow(merged)], df$end[i])
      } else {
        merged <- rbind(merged, df[i, ])
      }
    }
  }
  
  # Gaps are spaces between merged intervals within the session window
  gaps <- data.frame(
    start = c(wide_start, merged$end),
    end   = c(merged$start, wide_end)
  ) %>% filter(start < end)
  
  if (nrow(gaps) == 0) {
    gaps_list[[j]] <- data.frame(
      gap_start = NA,
      gap_end   = NA,
      season    = unique(act_dt_temp$season),
      session   = unique(act_dt_temp$session),
      nest      = unique(act_dt_temp$nest),
      dur       = 0
    )
  } else {
    gaps_list[[j]] <- gaps %>%
      rename(gap_start = start, gap_end = end) %>%
      mutate(
        season  = unique(act_dt_temp$season),
        session = unique(act_dt_temp$session),
        nest    = unique(act_dt_temp$nest),
        dur     = as.numeric(difftime(gap_end, gap_start, units = "secs"))
      )
  }
}

gaps_df <- plyr::ldply(gaps_list, data.frame) %>%
  filter(!is.na(gap_start)) %>%
  mutate(
    gap_start = as_datetime(gap_start),
    gap_end   = as_datetime(gap_end),
    season    = as.character(season)
  )

# ── 7. FINAL SEX ATTRIBUTION (STRICT INITIAL GAP & TIERED LOGIC) ────────────

# 1. Baseline Shift Durations (UNCHANGED)
nest_baselines <- act_dt_nest %>%
  mutate(
    onset_dt = as_datetime(onset),
    end_dt   = as_datetime(end_trimmed),
    shift_dur = as.numeric(difftime(end_dt, onset_dt, units = "secs"))
  ) %>%
  group_by(season, session, nest) %>%
  summarise(
    avg_shift = mean(shift_dur, na.rm = TRUE),
    dynamic_long_threshold = avg_shift * 1.5, 
    .groups = "drop"
  )

gaps_attributed <- gaps_df %>%
  left_join(nest_baselines, by = c("season", "session", "nest")) %>%
  rowwise() %>%
  mutate(
    responsible_sex = {
      
      # --- HARD-CODED OVERRIDES FOR KNOWN PROBLEM NESTS ---
      if (season == "2021" && nest == "D17_16" && dur > 18000) { res <- "f" } 
      else if (season == "2025" && nest == "A20_1" && dur > 18000) { res <- "m" }
      else if (season == "2023" && nest == "K20_8" && dur > 18000) { res <- "m" }
      else if (season == "2021" && nest == "W39" && dur > 18000) { res <- "m" }
      else if (season == "2021" && nest == "K8" && session == "incubation2" && dur > 18000) { res <- "f" }
      else if (season == "2021" && nest == "K8"&& session =="incubation3" && dur > 18000) { res <- "m" }
      else if (season == "2021" && nest == "D15_6"&& session =="incubation2" && dur > 18000) { res <- "f" }
      else if (season == "2023" && nest == "A20_4"&& session =="incubation3" && dur > 18000) { res <- "f" }
      else {
        n_data <- act_dt_nest %>%
          filter(season == .data$season, session == .data$session, nest == .data$nest) %>%
          mutate(onset_dt = as_datetime(onset), end_dt = as_datetime(end_trimmed))
        
        # Who was here immediately before and after?
        last_att <- n_data %>% filter(end_dt <= gap_start) %>% arrange(end_dt) %>% slice_tail(n = 1)
        next_att <- n_data %>% filter(onset_dt >= gap_end) %>% arrange(onset_dt) %>% slice(1)
        
        # 1. THE "BEGINNING" RULE (Session starts with a gap)
        # If there is NO bird recorded before this gap, 
        # it is the fault of the partner of the bird who finally arrives.
        if (nrow(last_att) == 0) {
          if (nrow(next_att) == 1) {
            # If Female arrives first, blame Male. If Male arrives first, blame Female.
            res <- if(next_att$sx == "m") "f" else "m"
          } else {
            res <- NA_character_
          }
        }
        
        # 2. TIER 1: SHORT GAPS (< 1 hour)
        else if (dur < 18000) {
          if (nrow(last_att) == 1 && nrow(next_att) == 1 && last_att$sx == next_att$sx) {
            res <- last_att$sx
          } else if (nrow(next_att) == 1) {
            res <- next_att$sx
          } else {
            res <- last_att$sx
          }
        } 
        
        # 3. TIER 2: LONG GAPS (> 1 hour)
        else {
          p_sex <- if(nrow(last_att) > 0) (if(last_att$sx == "m") "f" else "m") else NA_character_
          L_LIMIT <- if(!is.na(dynamic_long_threshold)) dynamic_long_threshold else 21600
          
          p_last_seen <- n_data %>% filter(sx == p_sex, end_dt <= gap_start) %>% arrange(end_dt) %>% slice_tail(n = 1)
          p_abs_dur <- if(nrow(p_last_seen) > 0) as.numeric(difftime(gap_start, p_last_seen$end_dt, units = "secs")) else 0
          
          if (!is.na(p_sex) && p_abs_dur > L_LIMIT) {
            res <- p_sex
          } else if (nrow(last_att) == 1 && nrow(next_att) == 1 && last_att$sx == next_att$sx) {
            res <- last_att$sx
          } else {
            res <- if(nrow(next_att) == 1) next_att$sx else last_att$sx
          }
        }
      }
      res
    }
  ) %>%
  ungroup()
# ── 8. WIDE SUMMARY: season × session × nest × sex ───────────────────────────

# All valid bird × session combinations (preserves zero-gap birds)
all_combinations <- act_dt_nest %>%
  distinct(season, session, nest, sx) %>%
  mutate(season = as.character(season))

# Incubation stats per bird
incubation_stats <- act_dt_nest %>%
  mutate(
    season    = as.character(season),
    shift_dur = as.numeric(difftime(
      as_datetime(end_trimmed), as_datetime(onset), units = "secs"
    ))
  ) %>%
  group_by(season, session, nest, sx) %>%
  summarise(
    n_incubations            = n(),
    mean_incubation_duration = mean(shift_dur, na.rm = TRUE),
    .groups = "drop"
  )

# Neglect summary per responsible sex
neglect_summary <- gaps_attributed %>%
  filter(!is.na(responsible_sex)) %>%
  group_by(season, session, nest, sx = responsible_sex) %>%
  summarise(
    n_gaps           = n(),
    sum_neglect_sec  = sum(dur, na.rm = TRUE),
    mean_neglect_sec = mean(dur, na.rm = TRUE),
    .groups = "drop"
  )

neglect_wide <- all_combinations %>%
  left_join(neglect_summary,  by = c("season", "session", "nest", "sx")) %>%
  left_join(incubation_stats, by = c("season", "session", "nest", "sx")) %>%
  mutate(
    n_gaps                   = replace_na(n_gaps, 0),
    sum_neglect_sec          = replace_na(sum_neglect_sec, 0),
    mean_neglect_sec         = replace_na(mean_neglect_sec, 0),
    n_incubations            = replace_na(n_incubations, 0),
    mean_incubation_duration = replace_na(mean_incubation_duration, 0)
  )

# ── 9. JOIN PHENOLOGY & PAIR BOND (identical to Script 1) ────────────────────

act_pheno <- act_dt %>%
  select(season, session, nest, session_start) %>%
  distinct() %>%
  mutate(season = as.numeric(season))

pheno_dt <- read_excel("C:/Users/Martyna/Dropbox/Hornsund metabase/Hornsund/SHO_LIAK_phenology_productivity.xlsx") %>%
  select(Season, Nest, Hatch_succ, Hatch_date, Lay_day) %>%
  mutate(season = as.numeric(Season)) %>%
  rename(nest = Nest, hatch_succ = Hatch_succ, hatch_date = Hatch_date)

act_pheno_temp <- left_join(act_pheno, pheno_dt, by = c("season", "nest")) %>%
  mutate(
    session_start   = as_datetime(session_start),
    hatch_date      = as_datetime(hatch_date),
    day_prior_hatch = -as.numeric(round(difftime(session_start, hatch_date, units = "days"), 0))
  )

pair_bond <- read.csv2("SI_PS.csv") %>%
  rename(nest = Nest, season = Season)

date_data <- readRDS("DPA_d05_status_activity data.rds") %>%
  mutate(
    session_start       = ymd_hms(session_start),
    session_end_trimmed = ymd_hms(session_end_trimmed),
    mean_date           = as.POSIXct(
      (as.numeric(session_start) + as.numeric(session_end_trimmed)) / 2,
      origin = "1970-01-01"
    )
  ) %>%
  select(season, session, nest, mean_date, session_end_trimmed) %>%
  distinct() %>%
  mutate(
    Julian_Date_point = yday(mean_date),
    season            = as.character(season)
  )

neglect_wide <- neglect_wide %>%
  mutate(season = as.numeric(season)) %>%
  left_join(act_pheno_temp, by = c("season", "session", "nest")) %>%
  left_join(pair_bond,      by = c("season", "nest")) %>%
  mutate(
    season              = as.character(season),
    J_hatchday          = yday(hatch_date),
    Incubation_Duration = J_hatchday - Lay_day
  ) %>%
  left_join(date_data, by = c("season", "session", "nest")) %>%
  arrange(season, nest, session, sx)

saveRDS(neglect_wide, "EDA_sex_neglect_wide.RDS")
message("Wide table saved: ", nrow(neglect_wide), " rows.")

# ── 10. GANTT VERIFICATION PLOTS ─────────────────────────────────────────────
# Two rows per session panel: m (top) and f (bottom)
# Green = nest attended | Grey = absence
# Blue  = gap attributed to male | Red = gap attributed to female

dir.create("VERIFICATION_SEX", showWarnings = FALSE)

# Full behaviour data for plotting (not filtered to attendance only)
gantt_data <- act_dt %>%
  mutate(
    season_nest = paste(season, "_", nest, sep = ""),
    onset       = as_datetime(onset),
    end_trimmed = as_datetime(end_trimmed),
    basic_interval = case_when(
      grepl("nest_attendance", basic_interval) ~ "nest_attendance",
      grepl("absence",         basic_interval) ~ "absence",
      TRUE ~ basic_interval
    )
  ) %>%
  filter(
    session %in% c("incubation1", "incubation2", "incubation3"),
    season_nest %notin% cases_excl
  )

plot_gantt_sex <- function(sn, ssn) {
  
  df   <- gantt_data    %>% filter(season == sn, nest == ssn)
  gaps <- gaps_attributed %>% filter(season == as.character(sn), nest == ssn)
  
  if (nrow(df) == 0) return(invisible(NULL))
  
  p <- ggplot() +
    
    # Absence (grey)
    geom_segment(
      data = df %>% filter(basic_interval == "absence"),
      aes(x = onset, xend = end_trimmed, y = sx, yend = sx),
      color = "grey75", linewidth = 5
    ) +
    
    # Nest attendance (green)
    geom_segment(
      data = df %>% filter(basic_interval == "nest_attendance"),
      aes(x = onset, xend = end_trimmed, y = sx, yend = sx),
      color = "darkgreen", linewidth = 5
    ) +
    
    # Gaps coloured by responsible sex, drawn on that sex's row
    geom_segment(
      data = gaps %>% filter(!is.na(responsible_sex)),
      aes(x = gap_start, xend = gap_end,
          y = responsible_sex, yend = responsible_sex,
          color = responsible_sex),
      linewidth = 5, alpha = 0.85
    ) +
    
    # Duration labels above each gap
    geom_text(
      data = gaps %>% filter(!is.na(responsible_sex)),
      aes(x     = gap_start + (gap_end - gap_start) / 2,
          y     = responsible_sex,
          label = paste0(round(dur / 60, 1), "m")),
      vjust = -0.8, size = 2.8, color = "black"
    ) +
    
    scale_color_manual(
      values = c("m" = "#2166ac", "f" = "#d6604d"),
      labels = c("m" = "Male responsible", "f" = "Female responsible"),
      name   = "Gap attributed to"
    ) +
    scale_y_discrete(limits = c("m", "f"), labels = c("m" = "Male", "f" = "Female")) +
    facet_wrap(~ session, scales = "free_x", ncol = 1) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    ) +
    labs(
      title = paste0(sn, "_", ssn, "  |  Green = nest attended  |  Coloured = gap (sex responsible)"),
      y = NULL, x = "Time"
    )
  
  ggsave(
    filename  = paste0("VERIFICATION_SEX/VSEX_", sn, "_", ssn, ".jpg"),
    plot      = p,
    width     = 14,
    height    = 3.5 * length(unique(df$session)),
    limitsize = FALSE
  )
}

nest_ids <- act_dt_nest %>% distinct(season, nest)
results  <- map2(nest_ids$season, nest_ids$nest, safely(plot_gantt_sex))

failed <- keep(results, ~ !is.null(.x$error))
if (length(failed) > 0) {
  message(length(failed), " plot(s) failed:")
  walk(failed, ~ message(.x$error$message))
} else {
  message("All verification plots saved to ./VERIFICATION_SEX/")
}

# Pair negotiation overlap durations (same logic as classic script)
overlaps_list <- list()

for (k in seq_along(loop_units)) {
  
  act_dt_temp <- act_dt_nest %>%
    filter(loop_unit == loop_units[k])
  
  female <- act_dt_temp %>% filter(sx == "f")
  male   <- act_dt_temp %>% filter(sx == "m")
  
  overlaps <- data.frame()
  
  if (nrow(male) > 0 && nrow(female) > 0) {
    
    male <- male %>%
      mutate(
        onset = as_datetime(onset),
        end_trimmed = as_datetime(end_trimmed)
      )
    
    female <- female %>%
      mutate(
        onset = as_datetime(onset),
        end_trimmed = as_datetime(end_trimmed)
      )
    
    for (i in 1:nrow(male)) {
      for (j in 1:nrow(female)) {
        
        start_overlap <- max(male$onset[i], female$onset[j])
        end_overlap   <- min(male$end_trimmed[i], female$end_trimmed[j])
        
        if (!is.na(start_overlap) &&
            !is.na(end_overlap) &&
            start_overlap < end_overlap) {
          
          overlaps <- rbind(
            overlaps,
            data.frame(
              overlap_start = start_overlap,
              overlap_end   = end_overlap
            )
          )
        }
      }
    }
  }
  
  if (nrow(overlaps) == 0) {
    overlaps_list[[k]] <- data.frame(
      season = act_dt_temp$season[1],
      session = act_dt_temp$session[1],
      nest = act_dt_temp$nest[1],
      dur = 0
    )
  } else {
    overlaps_list[[k]] <- overlaps %>%
      mutate(
        season = act_dt_temp$season[1],
        session = act_dt_temp$session[1],
        nest = act_dt_temp$nest[1],
        dur = as.numeric(difftime(overlap_end, overlap_start, units = "secs"))
      )
  }
}

overlaps_summary <- plyr::ldply(overlaps_list, data.frame) %>%
  group_by(season, session, nest) %>%
  summarise(
    n_overlaps    = sum(dur > 0),
    mean_overlaps = mean(dur[dur > 0], na.rm = TRUE),
    sum_overlaps  = sum(dur, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_overlaps = replace_na(mean_overlaps, 0))

# Convert sex-wide into nest-level informative table
neglect_table_final <- neglect_wide %>%
  
  group_by(season, session, nest) %>%
  summarise(
    # Classic neglect
    n_gaps            = sum(n_gaps, na.rm = TRUE),
    mean_gaps         = mean(mean_neglect_sec[sum_neglect_sec > 0], na.rm = TRUE),
    sum_gaps          = sum(sum_neglect_sec, na.rm = TRUE),
    
    # Male-specific
    n_gaps_m          = sum(n_gaps[sx == "m"], na.rm = TRUE),
    sum_gaps_m        = sum(sum_neglect_sec[sx == "m"], na.rm = TRUE),
    mean_gaps_m       = mean(mean_neglect_sec[sx == "m"], na.rm = TRUE),
    
    # Female-specific
    n_gaps_f          = sum(n_gaps[sx == "f"], na.rm = TRUE),
    sum_gaps_f        = sum(sum_neglect_sec[sx == "f"], na.rm = TRUE),
    mean_gaps_f       = mean(mean_neglect_sec[sx == "f"], na.rm = TRUE),
    
    # Incubation behaviour
    n_incubations_m   = sum(n_incubations[sx == "m"], na.rm = TRUE),
    n_incubations_f   = sum(n_incubations[sx == "f"], na.rm = TRUE),
    
    mean_incubation_m = mean(mean_incubation_duration[sx == "m"], na.rm = TRUE),
    mean_incubation_f = mean(mean_incubation_duration[sx == "f"], na.rm = TRUE),
    
    # Metadata (preserve one value per nest-session)
    session_start         = first(session_start),
    hatch_succ            = first(hatch_succ),
    hatch_date            = first(hatch_date),
    Lay_day               = first(Lay_day),
    day_prior_hatch       = first(day_prior_hatch),
    Pair_status           = first(Pair_status),
    Similarity_Index      = first(Similarity_Index),
    Similarity_Index_Head = first(Similarity_Index_Head),
    Similarity_Index_Tarsus = first(Similarity_Index_Tarsus),
    J_hatchday            = first(J_hatchday),
    Incubation_Duration   = first(Incubation_Duration),
    mean_date             = first(mean_date),
    session_end_trimmed   = first(session_end_trimmed),
    Julian_Date_point     = first(Julian_Date_point),
    .groups = "drop"
  ) %>%
  
  left_join(overlaps_summary,
            by = c("season", "session", "nest")) %>%
  
  mutate(
    mean_gaps = replace_na(mean_gaps, 0),
    mean_gaps_m = replace_na(mean_gaps_m, 0),
    mean_gaps_f = replace_na(mean_gaps_f, 0),
    
    neglect_bias_male = sum_gaps_m - sum_gaps_f,
    neglect_ratio_male = ifelse(sum_gaps == 0, 0, sum_gaps_m / sum_gaps)
  ) %>%
  
  arrange(season, nest, session)

head(neglect_table_final)
saveRDS(neglect_table_final, "EDA_egg_neglect_index_sex.RDS")

gaps_attributed_hourly <- gaps_attributed %>%
  filter(!is.na(responsible_sex), !is.na(gap_start)) %>%
  mutate(
    hour   = hour(as_datetime(gap_start)),
    season = as.character(season)
  )

saveRDS(gaps_attributed_hourly, "EDA_sex_gaps_hourly.RDS")
