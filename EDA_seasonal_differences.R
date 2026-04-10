# ==============================================================================
# PROJECT: Seasonal Differences in Nest Neglect (2019 - 2025)
# PURPOSE: Bayesian analysis of n_gaps (Poisson) and Hurdle Model (sum_gaps)
# ==============================================================================

rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)

# Load Data
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")
all_colors <- c("2019"="#A50026", "2020"="#FDAE61", "2021"="#A6D96A", "2023"="#006837", "2025"="#313695")

# ==============================================================================
# MODULE 1: FREQUENCY MODEL (n_gaps)
# ==============================================================================

# 1.1 Data Prep ----
E_fit_lin <- E %>%
  filter(!is.na(day_prior_hatch), !is.na(n_gaps), !is.na(season),
         day_prior_hatch >= 2, day_prior_hatch <= 28) %>%
  mutate(dph_c = day_prior_hatch - mean(day_prior_hatch),
         season = factor(season))

season_levels <- levels(E_fit_lin$season)
n_seasons <- length(season_levels)
d_lin <- list(n_gaps = E_fit_lin$n_gaps, dph_c = E_fit_lin$dph_c, season_id = as.integer(E_fit_lin$season))

# 1.2 Prior Predictive Simulation (n_gaps) ----
set.seed(123)
day_seq <- seq(2, 28, length.out = 100)
dph_seq_c <- day_seq - mean(E_fit_lin$day_prior_hatch)
n_sim <- 1000
prior_mu <- array(NA, dim = c(n_sim, length(day_seq), n_seasons))

for(i in 1:n_sim){
  a <- rnorm(n_seasons, 1, 0.5); b <- rnorm(1, 0, 0.05)
  for(s in 1:n_seasons) prior_mu[i,,s] <- exp(a[s] + b * dph_seq_c)
}
matplot(day_seq, t(prior_mu[1:20,,1]), type = "l", lty = 1, xlim = c(28, 2),
        xlab = "Days prior hatch", ylab = "Expected n_gaps", 
        main = paste("Prior predictive n_gaps:", season_levels[1]))

# 1.3 Model Fitting (Poisson) ----
m_lin <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a[season_id] + b * dph_c,
    a[season_id] ~ dnorm(1, 0.5),
    b ~ dnorm(0, 0.05)
  ), data = d_lin, start = list(a = rep(1, n_seasons), b = 0)
)

# 1.4 Posterior Plot (n_gaps) ----
post <- extract.samples(m_lin)
mu_list <- list(); PI_list <- list()
for(s in 1:n_seasons){
  p_link <- link(m_lin, data = list(season_id = rep(s, 100), dph_c = dph_seq_c))
  mu_list[[s]] <- apply(p_link, 2, mean)
  PI_list[[s]] <- apply(p_link, 2, PI, prob = 0.89)
}

cols_br <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00")
plot(NULL, xlim = c(28, 2), ylim = c(0, max(E_fit_lin$n_gaps)*1.2), 
     xlab = "Days prior hatch", ylab = "Number of nest exits", main = "Season differences in n_gaps")
for(s in 1:n_seasons){
  idx <- E_fit_lin$season == season_levels[s]
  points(jitter(E_fit_lin$day_prior_hatch[idx], amount = 0.25), 
         jitter(E_fit_lin$n_gaps[idx], amount = 0.1), pch = 16, col = col.alpha(cols_br[s], 0.45), cex = 0.8)
  shade(PI_list[[s]], day_seq, col = col.alpha(cols_br[s], 0.12))
  lines(day_seq, mu_list[[s]], lwd = 3.5, col = cols_br[s])
}
legend("topright", legend = season_levels, col = cols_br[1:n_seasons], pch = 16, lwd = 3, bty = "n")

# ==============================================================================
# MODULE 2: TWO-PART HURDLE MODEL (sum_gaps)
# ==============================================================================

# 2.1 Data Prep ----
E_fit_h <- E %>%
  filter(!is.na(day_prior_hatch), !is.na(sum_gaps), !is.na(season),
         day_prior_hatch >= 2, day_prior_hatch <= 28) %>%
  mutate(dph_c = day_prior_hatch - mean(day_prior_hatch),
         has_gap = ifelse(sum_gaps > 0, 1, 0),
         season = factor(season))
E_pos <- E_fit_h %>% filter(sum_gaps > 0)

# 2.2 Model 1: Occurrence ----
m_gap_occurrence <- quap(
  alist(has_gap ~ dbinom(1, p), logit(p) <- a[season_id] + b * dph_c,
        a[season_id] ~ dnorm(0, 1), b ~ dnorm(0, 0.2)),
  data = list(has_gap=E_fit_h$has_gap, dph_c=E_fit_h$dph_c, season_id=as.integer(E_fit_h$season))
)

# 2.3 Model 2: Duration ----
m_gap_duration <- quap(
  alist(sum_gaps ~ dgamma2(mu, scale), log(mu) <- a[season_id] + b * dph_c,
        a[season_id] ~ dnorm(8, 1), b ~ dnorm(0, 0.05), scale ~ dexp(1)),
  data = list(sum_gaps=E_pos$sum_gaps, dph_c=E_pos$dph_c, season_id=as.integer(E_pos$season)),
  start = list(a = rep(8, n_seasons), b = 0, scale = 1000)
)

# 2.4 Posterior Plots (Occurrence & Duration) ----
post_occ <- extract.samples(m_gap_occurrence)
post_dur <- extract.samples(m_gap_duration)

# Plot Occurrence
p_mu <- list(); p_PI <- list()
for(s in 1:n_seasons){
  l <- link(m_gap_occurrence, data = list(season_id = rep(s, 100), dph_c = dph_seq_c))
  p_mu[[s]] <- apply(l, 2, mean); p_PI[[s]] <- apply(l, 2, PI, prob = 0.89)
}
plot(NULL, xlim = c(28, 2), ylim = c(0, 1), xlab = "Days prior hatch", ylab = "P(any gap)", main = "Gap Occurrence")
for(s in 1:n_seasons){
  shade(p_PI[[s]], day_seq, col = col.alpha(cols_br[s], 0.12))
  lines(day_seq, p_mu[[s]], lwd = 3.5, col = cols_br[s])
}

# Plot Duration
dur_mu <- list(); dur_PI <- list()
for(s in 1:n_seasons){
  l <- link(m_gap_duration, data = list(season_id = rep(s, 100), dph_c = dph_seq_c))
  dur_mu[[s]] <- apply(l, 2, mean); dur_PI[[s]] <- apply(l, 2, PI, prob = 0.89)
}
plot(NULL, xlim = c(28, 2), ylim = c(0, 10000), xlab = "Days prior hatch", ylab = "Duration given gap", main = "Gap Duration")
for(s in 1:n_seasons){
  idx <- E_pos$season == season_levels[s]
  points(jitter(E_pos$day_prior_hatch[idx], 0.25), E_pos$sum_gaps[idx], pch=16, col=col.alpha(cols_br[s], 0.45), cex=0.8)
  shade(dur_PI[[s]], day_seq, col = col.alpha(cols_br[s], 0.12))
  lines(day_seq, dur_mu[[s]], lwd = 3.5, col = cols_br[s])
}

# ==============================================================================
# MODULE 3: COMPARISONS & SUMMARIES
# ==============================================================================

cat("\n=== n_gaps COMPARISONS ===\n")
for(i in 1:(n_seasons-1)) {
  for(j in (i+1):n_seasons) {
    diff <- post$a[,i] - post$a[,j]
    cat(season_levels[i], "vs", season_levels[j], ": LogDiff", round(mean(diff), 3), "| Prob(S1>S2):", round(mean(diff > 0), 3), "\n")
  }
}

cat("\n=== HURDLE COMPARISONS (Duration) ===\n")
for(i in 1:(n_seasons-1)) {
  for(j in (i+1):n_seasons) {
    diff <- post_dur$a[,i] - post_dur$a[,j]
    cat(season_levels[i], "vs", season_levels[j], ": Diff", round(mean(diff), 3), "| Prob(S1>S2):", round(mean(diff > 0), 3), "\n")
  }
}
