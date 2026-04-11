# PURPOSE: Model nest neglect (Frequency vs. Duration) in relation to Pair Status and Similarity 
# LIBRARIES
rm(list = ls())
library(rethinking)
library(dplyr)
library(ggplot2)


# =========================================================
# 1. DATA PREPARATION & GLOBAL SETTINGS
# =========================================================
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")

E <- E %>%
  mutate(
    pair_01     = ifelse(Pair_status == "Old", 1, 0),
    has_neglect = ifelse(sum_gaps > 0, 1, 0)
  )

# Filter global datasets
D_all     <- E %>% filter(!is.na(Pair_status), !is.na(n_gaps), !is.na(session))
D_nonzero <- D_all %>% filter(sum_gaps > 0)


# =========================================================
# 2. PAIR STATUS ANALYSIS (Aggregated Sessions)
# =========================================================

# --- 2.1 Models (session_id removed) ---
# Frequency Model (Poisson)
m_status_freq <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a + b * pair_01,
    a ~ dnorm(0, 0.8),
    b ~ dnorm(0, 0.5)
  ), data = D_all
)

# Duration Model (Gamma Hurdle)
m_status_dur <- quap(
  alist(
    sum_gaps ~ dgamma2(mu, scale), 
    log(mu) <- a + b * pair_01,
    a ~ dnorm(7, 1.5), 
    b ~ dnorm(0, 1),
    scale ~ dexp(1)
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
  mu = c(mean(exp(post_f$a)), mean(exp(post_f$a + post_f$b))),
  low = c(PI(exp(post_f$a))[1], PI(exp(post_f$a + post_f$b))[1]),
  high = c(PI(exp(post_f$a))[2], PI(exp(post_f$a + post_f$b))[2])
)

# Prepare Duration Model Estimate Dataframe
plot_dur_df <- data.frame(
  status = factor(c("New", "Old"), levels = c("New", "Old")),
  mu = c(mean(exp(post_d$a)), mean(exp(post_d$a + post_d$b))),
  low = c(PI(exp(post_d$a))[1], PI(exp(post_d$a + post_d$b))[1]),
  high = c(PI(exp(post_d$a))[2], PI(exp(post_d$a + post_d$b))[2])
)

# --- 2.4 Plotting ---

library(ggplot2)
library(patchwork)

# Plot A: Frequency
p_freq <- ggplot() +
  geom_violin(data = D_all, aes(x = Pair_status, y = n_gaps), 
              fill = "gray90", color = "gray80", alpha = 0.5) +
  geom_jitter(data = D_all, aes(x = Pair_status, y = n_gaps),
              width = 0.2, alpha = 0.3, color = "gray40", size = 1) +
  geom_errorbar(data = plot_freq_df, aes(x = status, ymin = low, ymax = high), 
                width = 0.1, linewidth = 1.2, color = "black") +
  geom_point(data = plot_freq_df, aes(x = status, y = mu), 
             size = 4, color = "black") +
  labs(title = "A)", x = "Pair Status", y = "Number of exits") +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(panel.grid = element_blank())

# Plot B: Duration
p_dur <- ggplot() +
  geom_violin(data = D_nonzero, aes(x = Pair_status, y = sum_gaps), 
              fill = "gray90", color = "gray80", alpha = 0.5) +
  geom_jitter(data = D_nonzero, aes(x = Pair_status, y = sum_gaps),
              width = 0.2, alpha = 0.3, color = "gray40", size = 1) +
  geom_errorbar(data = plot_dur_df, aes(x = status, ymin = low, ymax = high), 
                width = 0.1, linewidth = 1.2, color = "black") +
  geom_point(data = plot_dur_df, aes(x = status, y = mu), 
             size = 4, color = "black") +
  scale_y_log10(labels = scales::label_number()) +
  labs(title = "B)", x = "Pair Status", y = "Duration of neglect (s, log scale)") +
  theme_bw(base_size = 12, base_family = "sans") +
  theme(panel.grid = element_blank())

# Final Assembly
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
D_sim <- D_all %>% filter(!is.na(Similarity_Index), session == "incubation1")
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
