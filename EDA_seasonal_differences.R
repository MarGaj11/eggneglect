# PURPOSE: model seasonal differences in neglect, both n_gaps and sum_gaps
rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)
library(splines)

E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")
head(E)

# =========================================================
# LINEAR MODEL: overall n_gaps changes toward hatch + season
# =========================================================

# =========================================================
# 1. DATA
# =========================================================

E_fit_lin <- E %>%
  filter(
    !is.na(day_prior_hatch),
    !is.na(n_gaps),
    !is.na(season),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    dph_c = day_prior_hatch - mean(day_prior_hatch),
    season = factor(season)
  )

season_levels <- levels(E_fit_lin$season)
n_seasons <- length(season_levels)

d_lin <- list(
  n_gaps = E_fit_lin$n_gaps,
  dph_c = E_fit_lin$dph_c,
  season_id = as.integer(E_fit_lin$season)
)

# =========================================================
# 2. PRIOR PREDICTIVE SIMULATION
# =========================================================

set.seed(123)

day_seq <- seq(2, 28, length.out = 100)
dph_seq_c <- day_seq - mean(E_fit_lin$day_prior_hatch)

n_sim <- 1000
prior_mu <- array(NA, dim = c(n_sim, length(day_seq), n_seasons))

for(i in 1:n_sim){
  a <- rnorm(n_seasons, 1, 0.5)
  b <- rnorm(1, 0, 0.05)
  
  for(s in 1:n_seasons){
    prior_mu[i,,s] <- exp(a[s] + b * dph_seq_c)
  }
}

matplot(day_seq, t(prior_mu[1:20,,1]),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected n_gaps",
        main = paste("Prior predictive:", season_levels[1]))

# =========================================================
# 3. MODEL
# =========================================================

m_lin <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    
    log(lambda) <- a[season_id] + b * dph_c,
    
    a[season_id] ~ dnorm(1, 0.5),
    b ~ dnorm(0, 0.05)
  ),
  data = d_lin,
  start = list(
    a = rep(1, n_seasons),
    b = 0
  )
)

precis(m_lin, depth = 2)

# =========================================================
# 4. POSTERIOR PREDICTIONS
# =========================================================

post <- extract.samples(m_lin)

post_pred <- list()
mu_list <- list()
PI_list <- list()

for(s in 1:n_seasons){
  
  post_pred[[s]] <- link(
    m_lin,
    data = list(
      season_id = rep(s, length(day_seq)),
      dph_c = dph_seq_c
    )
  )
  
  mu_list[[s]] <- apply(post_pred[[s]], 2, mean)
  PI_list[[s]] <- apply(post_pred[[s]], 2, PI, prob = 0.89)
}

# =========================================================
# 5. POSTERIOR PLOT
# =========================================================

# nice colorblind-friendly palette
cols <- c(
  "#D55E00",  # vermillion
  "#0072B2",  # blue
  "#009E73",  # green
  "#CC79A7",  # purple
  "#E69F00"   # orange
)

plot(NULL,
     xlim = c(28, 2),
     ylim = c(0, max(E_fit_lin$n_gaps) * 1.2),
     xlab = "Days prior hatch",
     ylab = "Number of nest exits",
     main = "Season differences in nest exit rates")

# -----------------------------
# raw data colored by season
# -----------------------------
for(s in 1:n_seasons){
  
  idx <- E_fit_lin$season == season_levels[s]
  
  points(
    jitter(E_fit_lin$day_prior_hatch[idx], amount = 0.25),
    jitter(E_fit_lin$n_gaps[idx], amount = 0.1),
    pch = 16,
    col = col.alpha(cols[s], 0.45),
    cex = 0.8
  )
}

# -----------------------------
# posterior fits by season
# -----------------------------
for(s in 1:n_seasons){
  shade(PI_list[[s]], day_seq, col = col.alpha(cols[s], 0.12))
  lines(day_seq, mu_list[[s]], lwd = 3.5, col = cols[s])
}

legend("topleft",
       legend = season_levels,
       col = cols[1:n_seasons],
       pch = 16,
       lwd = 3,
       pt.cex = 1,
       bty = "n")

# =========================================================
# 6. INTERPRETATION OF SEASON EFFECTS
# =========================================================

cat("\n=== SEASON BASELINES ===\n")

for(s in 1:n_seasons){
  cat(
    season_levels[s], ":",
    round(mean(post$a[,s]), 3),
    "PI:", round(PI(post$a[,s]), 3), "\n"
  )
}

cat("\n=== COMMON SLOPE ===\n")
cat(
  round(mean(post$b), 4),
  "PI:", round(PI(post$b), 4),
  "P(<0):", round(mean(post$b < 0), 3), "\n"
)

# =========================================================
# 7. SEASON DIFFERENCES ON ORIGINAL SCALE
# =========================================================

cat("\n=== EXPECTED EXIT RATES BY SEASON ===\n")

season_rate <- exp(colMeans(post$a))

for(s in 1:n_seasons){
  cat(
    season_levels[s], ":",
    round(season_rate[s], 2),
    "expected gaps at mean incubation day\n"
  )
}


# Extract posterior samples
post <- extract.samples(m_lin)

# Function to compare two seasons
compare_seasons <- function(s1_idx, s2_idx, names) {
  diff <- post$a[,s1_idx] - post$a[,s2_idx]
  prob_gt <- mean(diff > 0)
  ci <- PI(diff, prob = 0.89)
  
  cat(names[s1_idx], "vs", names[s2_idx], ":\n")
  cat("  Mean Diff (log):", round(mean(diff), 3), "\n")
  cat("  89% PI:", round(ci, 3), "\n")
  cat("  Prob(S1 > S2):", round(prob_gt, 3), "\n\n")
}

# Run comparisons for all pairs
cat("=== PAIRWISE SEASONAL COMPARISONS ===\n\n")
for(i in 1:(n_seasons-1)) {
  for(j in (i+1):n_seasons) {
    compare_seasons(i, j, season_levels)
  }
}

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
p <- ggplot(E, aes(x = day_prior_hatch, y = n_gaps, color = Season, fill = Season)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 0.8) +
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  scale_color_manual(values = all_colors) +
  scale_fill_manual(values = all_colors) +
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


#seasonal differences for both sex sum_gaps 
library(rethinking)
library(tidyverse)
G <- readRDS("EDA_sex_gaps_hourly.RDS")
S <- readRDS("EDA_sex_neglect_wide.RDS")
dens(S$sum_neglect_sec)

D_neg <- S %>%
  filter(
    !is.na(sum_neglect_sec),
    !is.na(day_prior_hatch),
    !is.na(sx),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    sex_id = ifelse(sx == "f", 1, 2),
    dph_c = day_prior_hatch - mean(day_prior_hatch),
    
    # add tiny constant in case zeros exist
    neglect_log = log(sum_neglect_sec + 1)
  )

d_neg <- list(
  y = D_neg$neglect_log,
  sex_id = D_neg$sex_id,
  dph_c = D_neg$dph_c
)
# =========================================================
# SEASONAL DIFFERENCES IN SUM_GAPS
# =========================================================

D_season <- D_neg %>%
  mutate(
    season = as.character(season),
    sex_season = paste(sx, season, sep = "_"),
    sex_season_id = as.integer(factor(sex_season))
  )

group_lookup <- D_season %>%
  distinct(sex_season_id, sx, season) %>%
  arrange(sx, season)

n_groups <- length(unique(D_season$sex_season_id))

d_season <- list(
  y = D_season$neglect_log,
  sex_season_id = D_season$sex_season_id,
  dph_c = D_season$dph_c
)

m_neg_season <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    
    mu <- a[sex_season_id] + b[sex_season_id] * dph_c,
    
    a[sex_season_id] ~ dnorm(6, 0.5),
    b[sex_season_id] ~ dnorm(0, 0.08),
    sigma ~ dexp(1)
  ),
  data = d_season,
  start = list(
    a = rep(6, n_groups),
    b = rep(0, n_groups),
    sigma = 1
  )
)

precis(m_neg_season, depth = 2)

post_season <- extract.samples(m_neg_season)

par(mfrow = c(1,2))

season_levels <- sort(unique(D_season$season))
season_cols <- rainbow(length(season_levels))

for(sex_now in c("f", "m")) {
  
  # collect all posterior intercepts for axis range
  vals_now <- c()
  
  for(i in seq_along(season_levels)) {
    id_now <- group_lookup %>%
      filter(sx == sex_now, season == season_levels[i]) %>%
      pull(sex_season_id)
    
    vals_now <- c(vals_now, post_season$a[, id_now])
  }
  
  plot(NULL,
       xlim = range(vals_now),
       ylim = c(0, 1.2),
       xlab = "Log total neglect (sec)",
       ylab = "Posterior density",
       main = ifelse(sex_now == "f", "Females", "Males"))
  
  for(i in seq_along(season_levels)) {
    
    id_now <- group_lookup %>%
      filter(sx == sex_now, season == season_levels[i]) %>%
      pull(sex_season_id)
    
    dens(post_season$a[, id_now],
         col = season_cols[i],
         lwd = 3,
         add = TRUE)
  }
  
  legend("topright",
         legend = season_levels,
         col = season_cols,
         lwd = 3,
         cex = 0.8,
         bty = "n")
}

par(mfrow = c(1,1))

# =========================================================
# GGPLOT2 VERSION: Posterior Densities by Sex and Season
# =========================================================
library(ggplot2)

# 1. Prepare data for plotting
# We convert the list of samples into a clean long data frame
post_df <- data.frame()

for(i in 1:n_groups) {
  # Get group info from your lookup
  info <- group_lookup[group_lookup$sex_season_id == i, ]
  
  # Extract samples for the intercept 'a' for this group
  samples <- data.frame(
    neglect_log = post_season$a[, i],
    Sex = ifelse(info$sx == "f", "Females", "Males"),
    Season = as.character(info$season)
  )
  
  post_df <- rbind(post_df, samples)
}

# 2. Create the Plot
ggplot(post_df, aes(x = neglect_log, fill = Season, color = Season)) +
  # Use facet_wrap to separate Females and Males into different grids
  facet_wrap(~Sex, scales = "free_y") +
  
  # Add density curves with transparency (alpha)
  geom_density(alpha = 0.3, size = 1) +
  
  # Style and Labels
  theme_minimal() +
  labs(
    title = "Posterior Distribution of Log Total Neglect",
    subtitle = "Comparing seasonal intercepts between sexes",
    x = "Posterior Intercept (Log Seconds)",
    y = "Density"
  ) +
  
  # Use a nice color palette
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  
  # Final touches to make it look "pretty"
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.spacing = unit(2, "lines")
  )
# =========================================================
# KEY RESULTS: SEASON DIFFERENCES WITHIN SEX
# =========================================================

cat("\n=== SEASON DIFFERENCES IN TOTAL NEGLECT ===\n")

season_pairs <- combn(season_levels, 2, simplify = FALSE)

for(sex_now in c("f", "m")) {
  
  cat("\n---", ifelse(sex_now == "f", "FEMALES", "MALES"), "---\n")
  
  for(pair in season_pairs) {
    
    id1 <- group_lookup %>%
      filter(sx == sex_now, season == pair[1]) %>%
      pull(sex_season_id)
    
    id2 <- group_lookup %>%
      filter(sx == sex_now, season == pair[2]) %>%
      pull(sex_season_id)
    
    diff_post <- post_season$a[, id1] - post_season$a[, id2]
    
    cat(
      pair[1], "-", pair[2],
      "| Mean diff =", round(mean(diff_post), 3),
      "| PI =", round(PI(diff_post), 3),
      "| P(>0) =", round(mean(diff_post > 0), 3),
      "\n"
    )
  }
}



























# PURPOSE: model total gap duration toward hatching with season differences
# HURDLE APPROACH:
#   1) probability of any gap
#   2) duration conditional on gap > 0

rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)

E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")

# =========================================================
# 1. DATA
# =========================================================

E_fit <- E %>%
  filter(
    !is.na(day_prior_hatch),
    !is.na(sum_gaps),
    !is.na(season),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    dph_c = day_prior_hatch - mean(day_prior_hatch),
    has_gap = ifelse(sum_gaps > 0, 1, 0),
    season = factor(season)
  )

E_pos <- E_fit %>%
  filter(sum_gaps > 0)

season_levels <- levels(E_fit$season)
n_seasons <- length(season_levels)

# data for occurrence model
d_occ <- list(
  has_gap = E_fit$has_gap,
  dph_c = E_fit$dph_c,
  season_id = as.integer(E_fit$season)
)

# data for positive-duration model
d_dur <- list(
  sum_gaps = E_pos$sum_gaps,
  dph_c = E_pos$dph_c,
  season_id = as.integer(E_pos$season)
)

# prediction sequence
day_seq <- seq(2, 28, length.out = 100)
dph_seq_c <- day_seq - mean(E_fit$day_prior_hatch)

# =========================================================
# 2. PRIOR PREDICTIVE CHECKS
# =========================================================

set.seed(123)

n_sim <- 1000

# occurrence
prior_p <- array(NA, dim = c(n_sim, length(day_seq), n_seasons))

for(i in 1:n_sim){
  a <- rnorm(n_seasons, 0, 1)
  b <- rnorm(1, 0, 0.2)
  
  for(s in 1:n_seasons){
    logit_p <- a[s] + b * dph_seq_c
    prior_p[i,,s] <- inv_logit(logit_p)
  }
}

matplot(day_seq, t(prior_p[1:20,,1]),
        type = "l", lty = 1,
        xlim = c(28, 2),
        ylim = c(0, 1),
        xlab = "Days prior hatch",
        ylab = "P(any gap)",
        main = paste("Prior predictive occurrence:", season_levels[1]))

# duration
prior_mu <- array(NA, dim = c(n_sim, length(day_seq), n_seasons))

for(i in 1:n_sim){
  a <- rnorm(n_seasons, 8, 1)
  b <- rnorm(1, 0, 0.05)
  
  for(s in 1:n_seasons){
    prior_mu[i,,s] <- exp(a[s] + b * dph_seq_c)
  }
}

matplot(day_seq, t(prior_mu[1:20,,1]),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected duration if gap > 0",
        main = paste("Prior predictive duration:", season_levels[1]))

# =========================================================
# 3. MODEL 1: OCCURRENCE
# =========================================================

m_gap_occurrence <- quap(
  alist(
    has_gap ~ dbinom(1, p),
    
    logit(p) <- a[season_id] + b * dph_c,
    
    a[season_id] ~ dnorm(0, 1),
    b ~ dnorm(0, 0.2)
  ),
  data = d_occ,
  start = list(
    a = rep(0, n_seasons),
    b = 0
  )
)

precis(m_gap_occurrence, depth = 2)

# =========================================================
# 4. MODEL 2: POSITIVE DURATION
# =========================================================

m_gap_duration <- quap(
  alist(
    sum_gaps ~ dgamma2(mu, scale),
    
    log(mu) <- a[season_id] + b * dph_c,
    
    a[season_id] ~ dnorm(8, 1),
    b ~ dnorm(0, 0.05),
    scale ~ dexp(1)
  ),
  data = d_dur,
  start = list(
    a = rep(8, n_seasons),
    b = 0,
    scale = 1000
  )
)

precis(m_gap_duration, depth = 2)

# =========================================================
# 5. POSTERIOR PREDICTIONS
# =========================================================

# ---------- occurrence ----------
p_mu <- list()
p_PI <- list()

for(s in 1:n_seasons){
  
  link_s <- link(
    m_gap_occurrence,
    data = list(
      season_id = rep(s, length(day_seq)),
      dph_c = dph_seq_c
    )
  )
  
  p_mu[[s]] <- apply(link_s, 2, mean)
  p_PI[[s]] <- apply(link_s, 2, PI, prob = 0.89)
}

# ---------- duration ----------
dur_mu <- list()
dur_PI <- list()

for(s in 1:n_seasons){
  
  link_s <- link(
    m_gap_duration,
    data = list(
      season_id = rep(s, length(day_seq)),
      dph_c = dph_seq_c
    )
  )
  
  dur_mu[[s]] <- apply(link_s, 2, mean)
  dur_PI[[s]] <- apply(link_s, 2, PI, prob = 0.89)
}

# ---------- combined expectation ----------
comb_mu <- list()

for(s in 1:n_seasons){
  comb_mu[[s]] <- p_mu[[s]] * dur_mu[[s]]
}

# =========================================================
# 6. PLOTTING
# =========================================================

cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")

# -------------------------
# A) probability of any gap
# -------------------------
plot(NULL,
     xlim = c(28, 2),
     ylim = c(0, 1),
     xlab = "Days prior hatch",
     ylab = "Probability of any gap",
     main = "Season differences in gap occurrence")

for(s in 1:n_seasons){
  shade(p_PI[[s]], day_seq, col = col.alpha(cols[s], 0.12))
  lines(day_seq, p_mu[[s]], lwd = 3.5, col = cols[s])
}

legend("topleft",
       legend = season_levels,
       col = cols[1:n_seasons],
       lwd = 3,
       bty = "n")

# -------------------------
# B) duration if gap > 0
# -------------------------
plot(NULL,
     xlim = c(28, 2),
     ylim = c(0, max(E_pos$sum_gaps) * 1.1),
     xlab = "Days prior hatch",
     ylab = "Duration given gap",
     main = "Season differences in positive gap duration")

for(s in 1:n_seasons){
  
  idx <- E_pos$season == season_levels[s]
  
  points(
    jitter(E_pos$day_prior_hatch[idx], amount = 0.25),
    E_pos$sum_gaps[idx],
    pch = 16,
    col = col.alpha(cols[s], 0.45),
    cex = 0.8
  )
  
  shade(dur_PI[[s]], day_seq, col = col.alpha(cols[s], 0.12))
  lines(day_seq, dur_mu[[s]], lwd = 3.5, col = cols[s])
}

legend("topleft",
       legend = season_levels,
       col = cols[1:n_seasons],
       pch = 16,
       lwd = 3,
       bty = "n")

# -------------------------
# C) overall expected total
# -------------------------
plot(NULL,
     xlim = c(28, 2),
     ylim = c(0, max(E_fit$sum_gaps) * 1.1),
     xlab = "Days prior hatch",
     ylab = "Expected total gap duration",
     main = "Overall expected gap duration")

for(s in 1:n_seasons){
  lines(day_seq, comb_mu[[s]], lwd = 3.5, col = cols[s])
}

legend("topleft",
       legend = season_levels,
       col = cols[1:n_seasons],
       lwd = 3,
       bty = "n")

# =========================================================
# 7. INTERPRETATION
# =========================================================

cat("\n=== OCCURRENCE MODEL ===\n")
post_occ <- extract.samples(m_gap_occurrence)

cat("Slope:\n")
cat(
  round(mean(post_occ$b), 4),
  "PI:", round(PI(post_occ$b), 4),
  "P(<0):", round(mean(post_occ$b < 0), 3), "\n"
)

cat("\n=== DURATION MODEL ===\n")
post_dur <- extract.samples(m_gap_duration)

cat("Slope:\n")
cat(
  round(mean(post_dur$b), 4),
  "PI:", round(PI(post_dur$b), 4),
  "P(<0):", round(mean(post_dur$b < 0), 3), "\n"
)

# =========================================================
# 8. PAIRWISE SEASON COMPARISONS
# =========================================================

compare_seasons <- function(post_a, s1_idx, s2_idx, names, label) {
  diff <- post_a[,s1_idx] - post_a[,s2_idx]
  prob_gt <- mean(diff > 0)
  ci <- PI(diff, prob = 0.89)
  
  cat(label, "\n")
  cat(names[s1_idx], "vs", names[s2_idx], "\n")
  cat("  Mean Diff:", round(mean(diff), 3), "\n")
  cat("  89% PI:", round(ci, 3), "\n")
  cat("  Prob(S1 > S2):", round(prob_gt, 3), "\n\n")
}

cat("\n=== PAIRWISE OCCURRENCE COMPARISONS ===\n\n")
for(i in 1:(n_seasons-1)) {
  for(j in (i+1):n_seasons) {
    compare_seasons(post_occ$a, i, j, season_levels, "Occurrence")
  }
}

cat("\n=== PAIRWISE DURATION COMPARISONS ===\n\n")
for(i in 1:(n_seasons-1)) {
  for(j in (i+1):n_seasons) {
    compare_seasons(post_dur$a, i, j, season_levels, "Duration")
  }
}
