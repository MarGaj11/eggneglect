# PURPOSE: to calculate Egg Neglection Index, add the Pair characteristics, and Julian point dates
rm(list = ls())

# libs
library(readxl)
library(tidyverse)
library(tibble)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# Data
# filtered/cleaned and trimmed to 48h

act_dt <- readRDS("DPA_d05_status_activity data.rds")
act_dt_nest <- act_dt %>% 
  mutate(loop_unit = paste0(season, "_", session, "_", nest)) %>% 
  filter(basic_interval %in% c("start_nest_attendance", "nest_attendance", "end_nest_attendance"),
         session %in% c("incubation1", "incubation2", "incubation3")) 

#filtering nests 
act_dt_nest <- act_dt_nest %>% 
  mutate(season_nest = paste(season, "_", nest, sep = ""))

cases_excl <- c("2019_K16_4", "2019_14_18K", "2020_14_2M", "2020_A20_1", "2020_D15_6",
                "2021_113_10", "2021_12_16", "2021_14_18K", "2021_14_6K", "2021_A20_3", 
                "2021_D15_11", "2021_D17_8", "2021_W39B", "2023_12_16", 
                "2025_K20_10_22cmplx", "2025_K21_37", "2025_W4", "2025_K19_15_22cmplx", "2025_W40", "2025_A20_4")
`%notin%` <- Negate(`%in%`)

act_dt_nest <- act_dt_nest %>% filter(season_nest %notin% cases_excl)

act_dt_nest <- act_dt_nest %>%
  filter(!(season_nest == "2025_D15_11" & session == "incubation1") &
           !(season_nest == "2025_14_10M" & session == "incubation1") &
           !(season_nest == "2025_K20_8"  & session == "incubation3") &
           !(season_nest == "2023_12_16"  & session == "incubation2"))

# Variables to calculate ----
# incubation gap - absolute egg neglect (no parent in the nest) - mean and total
# negotiation intervals - (time spent in the nest by both partners) - mean and total/number 



# incubation gaps ----

loop_units <- unique(act_dt_nest$loop_unit)
gaps_table <- list()

for(j in 1:length(loop_units)) {
  
  # 1. Get the current unit
  act_dt_temp <- act_dt_nest %>% 
    filter(loop_unit == loop_units[j]) 
  
  # 2. CONVERT session intervals to datetime immediately
  wide_start <- lubridate::as_datetime(min(act_dt_temp$session_start)) 
  wide_end   <- lubridate::as_datetime(max(act_dt_temp$session_end_trimmed))
  
  # 3. Create df with proper datetime formats
  df <- data.frame(
    start = lubridate::as_datetime(act_dt_temp$onset),
    end   = lubridate::as_datetime(act_dt_temp$end_trimmed)
  )
  
  # Sort and merge overlapping intervals
  df <- df %>% arrange(start)
  
  # Initialize merged_intervals with the first row
  merged_intervals <- df[1,]
  
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      if (df$start[i] <= merged_intervals$end[nrow(merged_intervals)]) {
        merged_intervals$end[nrow(merged_intervals)] <- max(merged_intervals$end[nrow(merged_intervals)], df$end[i])
      } else {
        merged_intervals <- rbind(merged_intervals, df[i,])
      }
    }
  }
  
  # 4. Gaps calculation (Now that all components are datetimes, filter will work)
  gaps <- data.frame(
    start = c(wide_start, merged_intervals$end),
    end   = c(merged_intervals$start, wide_end)
  ) %>%
    filter(start < end)
  
  if(nrow(gaps) == 0) {
    gaps_table[[j]] <- data.frame(
      start = NA,
      end = NA, 
      season = unique(act_dt_temp$season),
      session = unique(act_dt_temp$session),
      nest = unique(act_dt_temp$nest),
      dur = 0)
    
  } else {
    gaps_table[[j]] <- gaps %>% 
      mutate(season = unique(act_dt_temp$season),
             session = unique(act_dt_temp$session),
             nest = unique(act_dt_temp$nest),
             dur = as.numeric(difftime(end, start, units = "secs"))) # Specified units for safety
  }  
}

saveRDS(gaps_table, "EDA_egg_neglect_index_gaps_table.RDS")

gaps_table_sum <- plyr::ldply(gaps_table, data.frame) %>%
  group_by(season, session, nest) %>%
  summarise(
    n_gaps     = sum(!is.na(start)),          # count only real gaps
    mean_gaps  = as.numeric(mean(dur[!is.na(start)], na.rm = TRUE)),
    sum_gaps   = as.numeric(sum(dur, na.rm = TRUE)),
    .groups    = "drop"
  ) %>%
  mutate(mean_gaps = ifelse(n_gaps == 0, 0, mean_gaps))  # NaN → 0 for zero-gap nests

# overlaps ----

overlaps_table <- list() 

for(k in 1:length(loop_units)) {
  
  act_dt_temp <- act_dt_nest %>% 
    filter(loop_unit == loop_units[k]) 
  
  # Rozdziel płeć
  female <- act_dt_temp %>% filter(sx == "f")
  male <- act_dt_temp %>% filter(sx == "m") 
  
  overlaps <- data.frame()
  
  # SPRAWDZENIE: Czy mamy kogo porównywać?
  if (nrow(male) > 0 && nrow(female) > 0) {
    
    # Konwersja na daty
    male$onset <- lubridate::as_datetime(male$onset)
    male$end_trimmed <- lubridate::as_datetime(male$end_trimmed)
    female$onset <- lubridate::as_datetime(female$onset)
    female$end_trimmed <- lubridate::as_datetime(female$end_trimmed)
    
    for (i in 1:nrow(male)) {
      for (j in 1:nrow(female)) {
        start_overlap <- max(male$onset[i], female$onset[j])
        end_overlap   <- min(male$end_trimmed[i], female$end_trimmed[j])
        
        if (!is.na(start_overlap) && !is.na(end_overlap) && start_overlap < end_overlap) {
          overlaps <- rbind(overlaps, data.frame(
            overlap_start = start_overlap,
            overlap_end   = end_overlap
          ))
        }
      }
    }
  }
  
  # Zapisywanie wyników
  if(nrow(overlaps) == 0) {
    overlaps_table[[k]] <- data.frame(
      overlap_start = NA, 
      overlap_end = NA,
      season = act_dt_temp$season[1], # Bezpieczniejsze niż unique
      session = act_dt_temp$session[1],
      nest = act_dt_temp$nest[1],
      dur = 0)
  } else {
    overlaps_table[[k]] <- overlaps %>% 
      mutate(
        season = act_dt_temp$season[1],
        session = act_dt_temp$session[1],
        nest = act_dt_temp$nest[1],
        dur = as.numeric(difftime(
          lubridate::as_datetime(overlap_end), 
          lubridate::as_datetime(overlap_start), 
          units = "secs"
        ))
      ) 
  }
} 

saveRDS(overlaps_table, "EDA_egg_neglect_index_overlaps_table.RDS")

overlaps_table_sum <- plyr::ldply(overlaps_table, data.frame)
overlaps_table_sum <- overlaps_table_sum %>% 
  group_by(season, session, nest) %>% 
  summarise(n_overlaps = n(),
            mean_overlaps = as.numeric(mean(dur, na.rm = TRUE)),
            sum_overlaps = as.numeric(sum(dur, na.rm = TRUE)))


# neglection time ----

# merge gaps and overlaps table
neglect_table <- left_join(overlaps_table_sum, gaps_table_sum, by = c("season", "session", "nest"))

# add pheno data - n days prior hatching
act_pheno <- act_dt %>% 
  select(season, session, nest, session_start) %>% 
  distinct()

pheno_dt <- read_excel("C:/Users/Martyna/Dropbox/Hornsund metabase/Hornsund/SHO_LIAK_phenology_productivity.xlsx")
pheno_dt <- pheno_dt %>% 
  select(Season, Nest, Hatch_succ, Hatch_date, Lay_day) %>% 
  mutate(season = as.numeric(Season)) %>% 
  rename(nest = "Nest",
         hatch_succ = "Hatch_succ",
         hatch_date = "Hatch_date") 

act_pheno <- act_pheno %>%
  mutate(season = as.numeric(season))

act_pheno_temp <- left_join(act_pheno, pheno_dt, by = c("season", "nest"))
act_pheno_temp <- act_pheno_temp %>% 
  mutate(
    # 1. Konwersja na format daty/czasu
    session_start = lubridate::as_datetime(session_start),
    hatch_date    = lubridate::as_datetime(hatch_date),
    
    # 2. Obliczenie różnicy
    day_prior_hatch = as.numeric(round(difftime(session_start, hatch_date, units = "days"), 0)),
    day_prior_hatch = -day_prior_hatch 
  )

neglect_table <- neglect_table %>%
  mutate(season = as.numeric(season))

neglect_table <- left_join(neglect_table, act_pheno_temp, by = c("season", "session", "nest")) 


# add pair status data
pair_bond <- read.csv2("SI_PS.csv")

pair_bond <- pair_bond %>%
  rename(nest = Nest,
         season = Season)

# merge with neglect table
neglect_table <- left_join(neglect_table, pair_bond, by = c("season", "nest"))

neglect_table <- neglect_table %>%
  filter(!session == "hatching")

neglect_table <- neglect_table %>%
  mutate(season = as.character(season),
         J_hatchday = yday(hatch_date),
         Incubation_Duration = J_hatchday - Lay_day)


#adding Julian Date point

date_data <- readRDS("./DPA_d05_status_activity data.rds") %>%
  mutate(
    # lubridate poradzi sobie z kropkami automatycznie
    session_start = ymd_hms(session_start),
    session_end_trimmed = ymd_hms(session_end_trimmed),
    
    # Obliczanie średniej daty na obiektach czasowych
    mean_date = as.POSIXct((as.numeric(session_start) + as.numeric(session_end_trimmed)) / 2, 
                           origin = "1970-01-01")
  ) %>%
  select(season, session, nest, mean_date, session_end_trimmed) %>%
  distinct() %>%
  mutate(
    Julian_Date_point = yday(mean_date),
    season = as.character(season),
  )

neglect_table <- neglect_table %>%
  left_join(date_data, by = c("season", "session", "nest"))


saveRDS(neglect_table, "EDA_egg_neglect_index_neglect_table.RDS")

#creating a data_frame for exits frequency/per hour analysis
nest_hour_grid <- act_dt_nest %>%
  group_by(season, session, nest) %>%
  summarise(
    start = lubridate::as_datetime(min(session_start)),
    end   = lubridate::as_datetime(max(session_end_trimmed)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(datetime_hour = list(seq(
    lubridate::floor_date(start, "hour"),
    lubridate::floor_date(end,   "hour"),
    by = "1 hour"
  ))) %>%
  unnest(datetime_hour) %>%
  mutate(
    date = as.Date(datetime_hour),
    hour = lubridate::hour(datetime_hour)
  ) %>%
  select(season, session, nest, date, hour)

neglect_vars <- neglect_table %>%
  select(season, session, nest, Pair_status, day_prior_hatch, Similarity_Index, ID1, ID2) %>%
  distinct()
# Join gaps against the full grid
gaps_table_sum_hourly <- plyr::ldply(gaps_table, data.frame) %>%
  filter(!is.na(start)) %>%
  mutate(
    start = lubridate::as_datetime(start),
    date  = as.Date(start),
    hour  = lubridate::hour(start)
  ) %>%
  group_by(season, session, nest, date, hour) %>%
  summarise(n_gaps_per_hour = n(), .groups = "drop") %>%
  right_join(nest_hour_grid, by = c("season", "session", "nest", "date", "hour")) %>%
  mutate(n_gaps_per_hour = replace_na(n_gaps_per_hour, 0))

gaps_table_sum_hourly <- gaps_table_sum_hourly %>%
  left_join(neglect_vars, by = c("season", "session", "nest"))

saveRDS(gaps_table_sum_hourly, "EDA_exit_hourly.RDS")

# --- VISUAL VERIFICATION OF GAPS (nest level) ---
dir.create("VERIFICATION", showWarnings = FALSE)

# Prepare gaps data with season_nest identifier
gaps_plot_data <- plyr::ldply(gaps_table, data.frame) %>%
  filter(!is.na(start)) %>%
  mutate(season_nest = paste(season, nest, sep = "_"),
         start = lubridate::as_datetime(start),
         end   = lubridate::as_datetime(end)) %>%
  # Join session boundaries from activity data
  left_join(
    act_dt_nest %>%
      mutate(season_nest = paste(season, nest, sep = "_")) %>%
      group_by(season_nest, session) %>%
      summarise(
        session_start = lubridate::as_datetime(min(session_start)),
        session_end   = lubridate::as_datetime(max(session_end_trimmed)),
        .groups = "drop"
      ),
    by = c("season_nest", "session")
  ) %>%
  # Clip gaps to video window
  mutate(
    start = pmax(start, session_start),
    end   = pmin(end,   session_end),
    dur   = as.numeric(difftime(end, start, units = "secs"))
  ) %>%
  filter(start < end)

# Prepare activity data — merge both sexes into one "nest" track
act_plot_data <- act_dt_nest %>%
  mutate(season_nest = paste(season, nest, sep = "_"),
         onset       = lubridate::as_datetime(onset),
         end_trimmed = lubridate::as_datetime(end_trimmed)) %>%
  group_by(season_nest, session, onset, end_trimmed) %>%
  summarise(occupied = any(basic_interval %in% c("nest_attendance", "start_nest_attendance")),
            .groups = "drop")

plot_gantt <- function(id) {
  
  df   <- act_plot_data %>% filter(season_nest == id)
  gaps <- gaps_plot_data %>% filter(season_nest == id)
  
  p <- ggplot() +
    # Nest occupied (green)
    geom_segment(
      data = df %>% filter(occupied),
      aes(x = onset, xend = end_trimmed, y = "nest", yend = "nest"),
      color = "darkgreen", linewidth = 8
    ) +
    # Egg neglect periods (red shading)
    geom_rect(
      data = gaps,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      fill = "red", alpha = 0.3
    ) +
    # Duration labels above each gap
    geom_text(
      data = gaps,
      aes(x = start + (end - start) / 2,
          y = Inf,
          label = paste0(round(dur / 60, 1), "min")),
      vjust = 1.5,
      size  = 3,
      color = "red"
    ) +
    scale_y_discrete(limits = "nest") +
    facet_wrap(~ session, scales = "free_x", ncol = 1) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold"),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste(id, "— Green = nest occupied | Red = egg neglect"),
      y     = NULL,
      x     = "Time"
    )
  
  ggsave(
    filename  = paste0("VERIFICATION/VER_", id, ".jpg"),
    plot      = p,
    width     = 14,
    height    = 3 * length(unique(df$session)),
    limitsize = FALSE
  )
}

# Run for all nests
results <- map(unique(act_plot_data$season_nest), safely(plot_gantt))

# Report failures
failed <- keep(results, ~ !is.null(.x$error))
if (length(failed) > 0) {
  message(length(failed), " plot(s) failed:")
  walk(failed, ~ message(.x$error$message))
} else {
  message("All plots saved to ./VERIFICATION/")
}
