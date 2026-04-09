# PURPOSE: to check if distribution of hours in the data is equal
library(rethinking)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)


rm(list = ls())

#ngl_dt_trimmed <- readRDS("./EDA_egg_neglect_index_neglect_table_trimmed.RDS") 
ngl_dt <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS") 

#n_gaps <- counts, poisson distribution
#sum_gaps <- positive, continous, skewed(gamma-like)

#firstly we want to equal the hours - so all the videos have the same distribution of hours of the video
#Define a continuous time window that is well represented across nests
ngl_dt_hours <- ngl_dt %>%
  rowwise() %>%
  mutate(hour_seq = list(seq(
    floor_date(session_start, "hour"),
    floor_date(session_end_trimmed, "hour"),
    by = "1 hour"
  ))) %>%
  unnest(hour_seq) %>%
  ungroup() %>%
  mutate(hour = hour(hour_seq),
         vid_duration = as.numeric(difftime(session_end_trimmed, session_start, units = "mins")))

coverage <- ngl_dt_hours %>%
  group_by(season, session, nest, hour) %>%
  summarise(n_obs = n(), .groups = "drop") %>%   # how many times this hour appears per nest
  group_by(season, session, hour) %>%
  summarise(
    n_nests = n_distinct(nest),   # coverage
    total_obs = sum(n_obs),       # total counts (effort)
    mean_obs_per_nest = mean(n_obs)
  )

ggplot(coverage, aes(hour, n_nests)) +
  geom_col() +
  facet_grid(season ~ session)

ggplot(coverage, aes(hour, mean_obs_per_nest)) +
  geom_col() +
  geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
  facet_grid(season ~ session)


#ngl_dt_hours_trimmed <- ngl_dt_trimmed %>%
#  rowwise() %>%
#  mutate(hour_seq = list(seq(
#    floor_date(session_start_analysis, "hour"),
#    floor_date(session_end_trimmed, "hour"),
#    by = "1 hour"
#  ))) %>%
#  unnest(hour_seq) %>%
#  ungroup() %>%
#  mutate(hour = hour(hour_seq),
#         vid_duration = as.numeric(difftime(session_end_trimmed, session_start_analysis, units = "mins")))

#coverage_trimmed <- ngl_dt_hours_trimmed %>%
#  group_by(season, session, nest, hour) %>%
# summarise(n_obs = n(), .groups = "drop") %>%   # how many times this hour appears per nest
#  group_by(season, session, hour) %>%
#  summarise(
#    n_nests = n_distinct(nest),   # coverage
#   total_obs = sum(n_obs),       # total counts (effort)
#    mean_obs_per_nest = mean(n_obs)
#  )

#ggplot(coverage_trimmed, aes(hour, n_nests)) +
#  geom_col() +
#  facet_grid(season ~ session)

#ggplot(coverage_trimmed, aes(hour, mean_obs_per_nest)) +
#  geom_col() +
#  geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
#  facet_grid(season ~ session)

#I decided to not use trimmed data as test shows that the difference between hours is not significant
#and without trimming the duration of the video can be equal for every nest


#testing if distribution of hour is equal
hourly_data <- coverage %>%
  group_by(hour) %>%
  summarise(total_effort = sum(total_obs))

chi_test <- chisq.test(hourly_data$total_effort)

print(chi_test)

chi_results_table <- coverage %>%
  group_by(season, session) %>%
  summarise(
    chi_stat = chisq.test(total_obs)$statistic,
    p_value = chisq.test(total_obs)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    Equal = ifelse(p_value > 0.05, "yes", "no")
  )

print(chi_results_table)