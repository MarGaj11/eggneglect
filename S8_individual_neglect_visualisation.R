rm(list = ls())

library(dplyr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(readxl)
library(lme4)

############### Individual stability and visualisation

setwd("C:/Users/Martyna/Dropbox/Negotiation/NI_files")
data <- read.csv("./Ind_neglect_Julian_Date.csv", sep=";", header=TRUE)

pair_ID <- read_excel("C:/Users/Martyna/Dropbox/Hornsund metabase/Hornsund/SHO_LIAK_capturing.xlsx", guess_max = 10000) %>%
  filter(Age != "1") %>%
  select(Nest, RingNo, R_leg, L_leg, Sx, Season, Wing, THL)

hatch <- read.csv("./Julian_Date_COINr_NI.csv", sep=";", header=TRUE) %>%
  mutate(Season_Nest = paste(Season, Nest, sep = "_"))

# Merge incubation timing
data <- data %>%
  mutate(Season_Nest = paste(season, nest, sep = "_")) %>%
  select(Season_Nest, nest, season, ringno, session, sx, NI, Partner_ringno) %>%
  left_join(
    hatch %>% select(Season_Nest, session, Days_to_Hatch),
    by = c("Season_Nest", "session")
  ) %>%
  mutate(
    Season = season,
    Nest = nest,
    Sx = sx
  )

# Prepare capture metadata
pair_ID <- pair_ID %>%
  mutate(Season_Nest = paste(Season, Nest, sep = "_"))

# Merge individual IDs
merged <- data %>%
  select(Season_Nest, Season, Nest, session, Sx, NI, Days_to_Hatch, Partner_ringno) %>%
  left_join(
    pair_ID %>% select(Season_Nest, Sx, RingNo),
    by = c("Season_Nest", "Sx")
  ) %>%
  distinct(Season_Nest, Sx, RingNo, session, .keep_all = TRUE) %>%
  mutate(
    Season = as.factor(Season),
    session = as.factor(session),
    Sx = as.factor(Sx),
    NI = as.numeric(gsub(",", ".", NI))
  )


# -------------------- PLOT 1: NI vs days to hatch by season and sex --------------------

colors <- c("2019" = "#A50026",
            "2020" = "#FDAE61",
            "2021" = "#A6D96A",
            "2023" = "#006837",
            "2025" = "#313695")

sample_sizes <- merged %>%
  group_by(Season) %>%
  summarise(n = n()/2, .groups = "drop")

season_labels <- setNames(
  paste0(sample_sizes$Season, " (n = ", sample_sizes$n, ")"),
  sample_sizes$Season
)

sex_labels <- data.frame(
  Sx = c("f", "m"),
  label = c("\u2640", "\u2642"),
  x = c(0, 0),
  y = c(max(merged$NI, na.rm = TRUE), max(merged$NI, na.rm = TRUE))
)

p2 <- ggplot(merged, aes(x = Days_to_Hatch, y = NI)) +
  geom_point(aes(color = Season), alpha = 0.5, size = 2.5) +
  geom_smooth(aes(color = Season, fill = Season),
              method = "lm", se = TRUE, alpha = 0.15) +
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  scale_color_manual(values = colors, labels = season_labels) +
  scale_fill_manual(values = colors, labels = season_labels) +
  labs(
    x = "Number of days before hatching",
    y = "Neglect Index (NI)"
  ) +
  geom_text(
    data = sex_labels,
    aes(x = x, y = y, label = label),
    hjust = 1, vjust = 1,
    fontface = "bold",
    size = 10,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~Sx) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

p2

ggsave("figure_NI_plot.png",
       plot = p2,
       width = 8,
       height = 5,
       units = "in",
       dpi = 400)

# -------------------- PLOT 2: Individual stability between seasons --------------------

individual_stability_data <- merged %>%
  select(Nest, RingNo, Sx, Season, session, NI) %>%
  group_by(Nest, RingNo, Sx, Season) %>%
  filter(n_distinct(session) > 1) %>%
  summarise(
    Mean_Neglect = mean(NI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Nest_ID = as.factor(Nest))

p_individual_stability <- ggplot(individual_stability_data,
                                 aes(x = Season, y = Mean_Neglect,
                                     group = RingNo,
                                     color = Sx)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  facet_wrap(~ Nest_ID, ncol = 4) +
  scale_color_manual(values = c("m" = "#0072B2", "f" = "#D55E00"),
                     name = "Sex") +
  labs(
    title = "Individual neglect across seasons",
    subtitle = "Birds with data from >1 season and >1 incubation session",
    y = "Mean Neglect Index",
    x = "Season"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_individual_stability
