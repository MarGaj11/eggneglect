# PURPOSE: Model nest neglect (Frequency vs. Duration) in relation to Pair Status and Similarity 
# LIBRARIES
rm(list = ls())
library(rethinking)
library(dplyr)
library(ggplot2)
library(tidyr)


# =========================================================
# 1. DATA PREPARATION & GLOBAL SETTINGS
# =========================================================
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")
Y <- readRDS("./EDA_egg_neglect_index_gaps_table.RDS")

# --- Filter short gaps and rebuild E summaries ---
Y_filtered <- bind_rows(Y) %>%          # flatten list to single dataframe
  filter(dur >= 300)                     # keep only gaps >= 5 min

Y_summary <- Y_filtered %>%
  group_by(season, session, nest) %>%
  summarise(
    n_gaps    = n(),
    sum_gaps  = sum(dur),
    mean_gaps = mean(dur),
    .groups   = "drop"
  )

E <- E %>%
  select(-n_gaps, -sum_gaps, -mean_gaps) %>%
  left_join(Y_summary, by = c("season", "session", "nest")) %>%
  mutate(
    n_gaps    = replace_na(n_gaps,    0),
    sum_gaps  = replace_na(sum_gaps,  0),
    mean_gaps = replace_na(mean_gaps, 0)
  )

nest_summary <- E %>%
  group_by(season, Pair_status) %>%
  summarise(n_distinct_nests = n_distinct(nest), .groups = 'drop')

E <- E %>%
  mutate(
    pair_01     = ifelse(Pair_status == "Old", 2, 1),
    has_neglect = ifelse(sum_gaps > 0, 1, 0)
  )

E <- E %>%
  mutate(
    session = case_when(
      day_prior_hatch >= 21 & day_prior_hatch <= 30 ~ "early",
      day_prior_hatch >= 11 & day_prior_hatch <= 20 ~ "mid",
      day_prior_hatch >= 0  & day_prior_hatch <= 10 ~ "late",
      TRUE ~ NA_character_),
    dph_c = day_prior_hatch - mean(day_prior_hatch)
  )

# Filter global datasets
D_all     <- E %>% filter(!is.na(Pair_status), !is.na(n_gaps), !is.na(session))
D_nonzero <- D_all %>% filter(sum_gaps > 0)


# =========================================================
# 2. PAIR STATUS ANALYSIS (Aggregated Sessions)
# =========================================================

m_status_freq <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a[pair_01] + b* dph_c,
    a[pair_01] ~ dnorm(0, 0.8),
    b ~ dnorm(0, 0.5)
  ), data = D_all
)

m_status_dur <- quap(
  alist(
    sum_gaps ~ dnorm(mu, sigma),
    mu <- a[pair_01] + b* dph_c,
    a[pair_01] ~ dnorm(6, 0.5),
    b ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = D_nonzero
)


# --- 2.2 Summaries ---
post_f <- extract.samples(m_status_freq)
post_d <- extract.samples(m_status_dur)

cat("\n--- PAIR STATUS RESULTS (Global Intercept) ---\n")
cat("FREQ: Rate ratio (Old/New):", round(mean(exp(post_f$b)), 3), " | 89% PI:", round(PI(exp(post_f$b)), 3), "\n")
cat("DUR:  Duration Ratio (Old/New):", round(mean(exp(post_d$b)), 3), " | 89% PI:", round(PI(exp(post_d$b)), 3), "\n")
cat("DUR:  P(Old Duration > New Duration):", round(mean(post_d$b > 0), 3), "\n")

# --- 2.3 Visualization Preparation ---

# Prepare Frequency Model Estimate Dataframe
plot_freq_df <- data.frame(
  status = factor(c("New", "Old"), levels = c("New", "Old")),
  mu   = c(mean(exp(post_f$a)), mean(exp(post_f$a + post_f$b))),
  low  = c(PI(exp(post_f$a))[1], PI(exp(post_f$a + post_f$b))[1]),
  high = c(PI(exp(post_f$a))[2], PI(exp(post_f$a + post_f$b))[2])
)

plot_dur_df <- data.frame(
  status = factor(c("New", "Old"), levels = c("New", "Old")),
  mu   = c(mean(exp(post_d$a)), mean(exp(post_d$a + post_d$b))),
  low  = c(PI(exp(post_d$a))[1], PI(exp(post_d$a + post_d$b))[1]),
  high = c(PI(exp(post_d$a))[2], PI(exp(post_d$a + post_d$b))[2])
)


# plot posteriors
f_sim <- sim(m_status_freq, data = list(pair_01 = 1:2, dph_c = c(0,0)), n = 2000)
d_sim <- sim(m_status_dur,  data = list(pair_01 = 1:2, dph_c = c(0,0)), n = 2000)

df_post_f <- data.frame(
  Pair_status = rep(c("New", "Old"), each = 2000),
  val = as.vector(f_sim)
)

df_post_d <- data.frame(
  Pair_status = rep(c("New", "Old"), each = 2000),
  val = as.vector(d_sim)
)


df_post_d_clean <- df_post_d
df_post_d_clean$val[df_post_d_clean$val <= 0] <- 1 # Zamiana na 1s, by log10 działał

# --- Wykres A: Bez zmian (Poisson nie ma problemu z log) ---
p_freq <- ggplot() +
  geom_jitter(data = D_all, aes(x = factor(pair_01, labels=c("New", "Old")), y = n_gaps),
              width = 0.2, alpha = 0.2, size = 1, color = "gray70") +
  geom_violin(data = df_post_f, aes(x = Pair_status, y = val), 
              fill = "gray95", color = "black", linewidth = 0.5, alpha = 0.5) +
  geom_errorbar(data = plot_freq_df, aes(x = status, ymin = low, ymax = high), 
                width = 0.05, linewidth = 1, color = "black") +
  geom_point(data = plot_freq_df, aes(x = status, y = mu), 
             size = 3.5, color = "black") +
  coord_cartesian(ylim = c(0, 5)) +
  labs(title = "A) Gap Frequency", x = "Pair Status", y = "Number of gaps") +
  theme_classic(base_size = 12)

# --- Wykres B: Poprawiony pod skalę logarytmiczną ---
p_dur <- ggplot() +
  geom_jitter(data = D_nonzero, aes(x = factor(pair_01, labels=c("New", "Old")), y = sum_gaps),
              width = 0.2, alpha = 0.2, size = 1, color = "gray70") +
  # Używamy wyczyszczonych danych (df_post_d_clean)
  geom_violin(data = df_post_d_clean, aes(x = Pair_status, y = val), 
              fill = "gray95", color = "black", linewidth = 0.5, alpha = 0.5) +
  geom_errorbar(data = plot_dur_df, aes(x = status, ymin = low, ymax = high), 
                width = 0.05, linewidth = 1, color = "black") +
  geom_point(data = plot_dur_df, aes(x = status, y = mu), 
             size = 3.5, color = "black") +
  # Ustawienie limitów wewnątrz skali, by uniknąć wyrzucania punktów przez jitter
  scale_y_log10(labels = scales::label_comma(), 
                breaks = c(1, 100, 1000, 10000, 100000)) +
  labs(title = "B) Neglect Duration", x = "Pair Status", y = "Duration (sec)") +
  theme_classic(base_size = 12)


# Połączenie wykresów
final <- p_freq + p_dur
final

ggsave(
  "Figure_4_pair_status_collapsed.png",
  final,
  width = 8,
  height = 6,
  dpi = 600
)

# =========================================================
# 3. PAIR SIMILARITY ANALYSIS - here I use only incubation1 as it shows the highest number of exits 
# =========================================================

# --- 3.1 Models (incubation1 only) ---
D_sim <- D_all %>% filter(!is.na(Similarity_Index), session == "early")
D_sim_nonzero <- D_sim %>% filter(sum_gaps > 0)

m_sim_freq <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a + b * Similarity_Index,
    a ~ dnorm(0, 0.8),
    b ~ dnorm(0, 0.5)
  ), data = D_sim
)

m_sim_dur <- quap(
  alist(
    sum_gaps ~ dgamma2(mu, scale),
    log(mu) <- a + b * Similarity_Index,
    a ~ dnorm(7, 1.5), 
    b ~ dnorm(0, 0.5),
    scale ~ dexp(1)
  ), 
  data = D_sim_nonzero,
  control = list(maxit = 2000)
)

# --- 3.2 Summaries ---
post_sim_f <- extract.samples(m_sim_freq)
post_sim_d <- extract.samples(m_sim_dur)

cat("\n--- PAIR SIMILARITY RESULTS ---\n")
cat("FREQ: P(Slope b > 0):", round(mean(post_sim_f$b > 0), 3), "\n")
cat("DUR:  P(Slope b > 0):", round(mean(post_sim_d$b > 0), 3), "\n")
cat("DUR:  Ratio (High/Low Similarity):", round(mean(exp(post_sim_d$b)), 3), "\n")

# 1. Prepare Prediction Data for Regressions
xseq <- seq(0, 1, length.out = 100)

# Predictions for Frequency
link_f <- link(m_sim_freq, data = list(Similarity_Index = xseq))
plot_sim_f <- data.frame(
  x = xseq,
  mu = apply(link_f, 2, mean),
  low = apply(link_f, 2, PI, prob = 0.89)[1,],
  high = apply(link_f, 2, PI, prob = 0.89)[2,]
)

# Predictions for Duration
link_d <- link(m_sim_dur, data = list(Similarity_Index = xseq))
plot_sim_d <- data.frame(
  x = xseq,
  mu = apply(link_d, 2, mean),
  low = apply(link_d, 2, PI, prob = 0.89)[1,],
  high = apply(link_d, 2, PI, prob = 0.89)[2,]
)

# --- 1. SET A SHARED SIZE LIMIT ---
# Adjust these based on your data (e.g., if your max overlap is 15, set to 15)
max_obs <- 10 

# --- 1. SET THE SHARED LEGEND RULES ---
# This ensures "n=5" looks the same in both plots
shared_size_scale <- scale_size_continuous(
  name = "n obs",          # Title of the legend
  range = c(1.5, 8),       # Physical size of bubbles
  limits = c(1, 10),       # Force 1 to 10 scale (adjust 10 to your max N)
  breaks = c(1, 2, 5, 10)  # Which numbers to show in the legend
)

# --- 2. PLOT C: FREQUENCY ---
p_sim_freq <- ggplot() +
  geom_count(data = D_sim, aes(x = Similarity_Index, y = n_gaps), 
             alpha = 0.2, color = "gray40") +
  geom_ribbon(data = plot_sim_f, aes(x = x, ymin = low, ymax = high), 
              fill = "black", alpha = 0.1) +
  geom_line(data = plot_sim_f, aes(x = x, y = mu), 
            color = "black", linewidth = 1.2) +
  
  shared_size_scale +      # Add the shared legend here
  
  labs(title = "A)", x = "Similarity Index", y = "Number of exits") +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",        # Put the legend back
    legend.background = element_blank()
  )

# --- 3. PLOT D: DURATION ---
p_sim_dur <- ggplot() +
  geom_count(data = D_sim_nonzero, aes(x = Similarity_Index, y = sum_gaps), 
             alpha = 0.2, color = "gray40") +
  geom_ribbon(data = plot_sim_d, aes(x = x, ymin = low, ymax = high), 
              fill = "black", alpha = 0.1) +
  geom_line(data = plot_sim_d, aes(x = x, y = mu), 
            color = "black", linewidth = 1.2) +
  
  scale_y_log10(labels = scales::label_number()) +
  shared_size_scale +      # Add the same shared legend here
  
  labs(title = "B)", x = "Similarity Index", y = "Duration of neglect in second (log scale)") +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",        # Put the legend back
    legend.background = element_blank()
  )

# Combine
p_final <- p_sim_freq + p_sim_dur + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "right")

print(p_final)

ggsave(
  "Figure_3_similarity_index.png",
  p_final,
  width = 8,
  height = 6,
  dpi = 600
)


