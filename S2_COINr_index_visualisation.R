rm(list = ls())
library(dplyr)
library(lubridate)
library(ggplot2)
library(readxl)
library(ggthemes)
library(readr)

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")

# --------------------------
# Load COINr NI data
# --------------------------
data <- read.csv2("./COINr_NI_neglect_only.csv")  # semicolon-separated
data <- data[, -1]  # remove extra index column
colnames(data)[1:3] <- c("Season", "session", "Nest")
data$Season <- as.character(data$Season)

# --------------------------
# Load hatching info
# --------------------------
hatch_data <- read_excel("C:/Users/Martyna/Dropbox/Hornsund metabase/Hornsund/SHO_LIAK_phenology_productivity.xlsx") %>%
  select(Season, Nest, Lay_day, Hatch_day)

# Merge hatching info with NI data
merged_data <- data %>%
  left_join(hatch_data, by = c("Season", "Nest"))

# --------------------------
# Load mean session date info
# --------------------------
date_data <- readRDS("./DPA_d05_status_activity data.rds") %>%
  mutate(
    session_start = as.POSIXct(gsub("\\.", "-", session_start)),
    session_end   = as.POSIXct(gsub("\\.", "-", session_end)),
    Mean_Date = as.POSIXct((as.numeric(session_start) + as.numeric(session_end)) / 2, origin = "1970-01-01")
  ) %>%
  select(season, session, nest, Mean_Date) %>%
  distinct() %>%
  rename(Season = season, Nest = nest, mean_date = Mean_Date) %>%
  mutate(Julian_Date_point = yday(mean_date),
         Season = as.character(Season))

# Merge with NI + hatching info
df <- merged_data %>%
  left_join(date_data, by = c("Season", "session", "Nest")) %>%
  mutate(
    NI = as.numeric(gsub(",", ".", NI)),  # fix decimal commas
    Incubation_Duration = Hatch_day - Lay_day,
    Days_to_Hatch = Hatch_day - Julian_Date_point
  ) %>%
  filter(!is.na(Days_to_Hatch))

# --------------------------
# Sample sizes for labeling
# --------------------------
season_counts <- df %>%
  group_by(Season) %>%
  summarise(n = n(), .groups = "drop")  # divide by 2 if 2 sexes

season_labels <- setNames(
  paste0(season_counts$Season, " (n=", season_counts$n, ")"),
  season_counts$Season
)

# Colors for seasons
all_colors <- c(
  "2019" = "#A50026",
  "2020" = "#FDAE61",
  "2021" = "#A6D96A",
  "2023" = "#006837",
  "2025" = "#313695"
)

# --------------------------
# Plot NI vs Days to Hatch
# --------------------------
p <- ggplot(df, aes(x = Days_to_Hatch, y = NI, color = Season, fill = Season)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  scale_color_manual(values = all_colors, labels = season_labels) +
  scale_fill_manual(values = all_colors, labels = season_labels) +
  labs(
    x = "Days before hatching",
    y = "Neglect Index (NI)",
    color = "Season",
    fill = "Season"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 13, face = "plain", color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )
p

# --------------------------
# Save plot
# --------------------------
ggsave("figure_NI_plot.png", plot = p, width = 8, height = 5, units = "in", dpi = 400)

# --------------------------
# Save merged dataset
# --------------------------
write.csv2(df, "Julian_Date_COINr_NI.csv", row.names = FALSE)
