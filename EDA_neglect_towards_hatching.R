# PURPOSE: model n_gaps and sum_gaps towards hatching
rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)
library(splines)

S <- readRDS("EDA_sex_neglect_wide.RDS")
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")

# Global sequence for predictions (0 is hatching, 28 is early incubation)
DAY_SEQ <- seq(2, 28, length.out = 100)

# =========================================================
# 1. MODEL: OVERALL FREQUENCY (n_gaps)
# =========================================================

# --- 1.1 Data Preparation ---
E_fit_overall <- E %>%
  filter(!is.na(day_prior_hatch), !is.na(n_gaps),
         day_prior_hatch >= 2, day_prior_hatch <= 28) %>%
  mutate(dph_c = day_prior_hatch - mean(day_prior_hatch))

d_overall <- list(n_gaps = E_fit_overall$n_gaps, dph_c = E_fit_overall$dph_c)

# --- 1.2 Prior Predictive Simulation ---
set.seed(123)
dph_seq_c_ov <- DAY_SEQ - mean(E_fit_overall$day_prior_hatch)
n_sim <- 100
prior_mu_ov <- matrix(NA, nrow = n_sim, ncol = length(DAY_SEQ))

for(i in 1:n_sim){
  a <- rnorm(1, 1, 0.5)
  b <- rnorm(1, 0, 0.05)
  prior_mu_ov[i,] <- exp(a + b * dph_seq_c_ov)
}

matplot(DAY_SEQ, t(prior_mu_ov), type = "l", lty = 1, xlim = c(28, 2),
        xlab = "Days prior hatch", ylab = "Expected n_gaps",
        main = "Prior predictive: Overall n_gaps")

# --- 1.3 Fit Model ---
m_overall <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a + b * dph_c,
    a ~ dnorm(1, 0.5),
    b ~ dnorm(0, 0.05)
  ), data = d_overall, start = list(a = 0, b = 0)
)

# --- 1.4 Posterior Plotting ---
post_mu_ov <- link(m_overall, data = list(dph_c = dph_seq_c_ov))
mu_ov   <- apply(post_mu_ov, 2, mean)
PI_mu_ov <- apply(post_mu_ov, 2, PI, prob = 0.89)

plot(NULL, xlim = c(28, 2), ylim = c(0, max(E_fit_overall$n_gaps) * 1.1),
     xlab = "Days prior hatch", ylab = "n_gaps", main = "Overall n_gaps Trend")
points(jitter(E_fit_overall$day_prior_hatch, 0.25), jitter(E_fit_overall$n_gaps, 0.1), 
       pch = 16, col = col.alpha("black", 0.3))
shade(PI_mu_ov, DAY_SEQ, col = col.alpha("black", 0.2))
lines(DAY_SEQ, mu_ov, lwd = 3, col = "darkgrey")

# --- 1.5 Interpretation ---
post_ov <- extract.samples(m_overall)
cat("\n=== OVERALL SLOPE ===\n")
cat("Slope:", round(mean(post_ov$b), 4), "PI:", round(PI(post_ov$b), 4),
    "P(<0):", round(mean(post_ov$b < 0), 3), "\n")
cat("Multiplicative change per day: exp(b) =", round(exp(mean(post_ov$b)), 3), "\n")


# =========================================================
# 2. MODEL: SEX-SPECIFIC FREQUENCY (n_gaps)
# =========================================================

# --- 2.1 Data Preparation ---
E_fit_sex <- S %>%
  filter(!is.na(day_prior_hatch), !is.na(sx),
         day_prior_hatch >= 2, day_prior_hatch <= 28) %>%
  mutate(sex_id = ifelse(sx == "f", 1, 2),
         dph_c = day_prior_hatch - mean(day_prior_hatch))

d_sex <- list(n_gaps = E_fit_sex$n_gaps, sex_id = E_fit_sex$sex_id, dph_c = E_fit_sex$dph_c)

# --- 2.2 Prior Predictive Simulation ---
prior_f <- matrix(NA, nrow = n_sim, ncol = length(DAY_SEQ))
prior_m <- matrix(NA, nrow = n_sim, ncol = length(DAY_SEQ))

for(i in 1:n_sim){
  a <- rnorm(2, 1, 0.5); b <- rnorm(2, 0, 0.1)
  prior_f[i,] <- exp(a[1] + b[1] * dph_seq_c_ov)
  prior_m[i,] <- exp(a[2] + b[2] * dph_seq_c_ov)
}

par(mfrow = c(1,2))
matplot(DAY_SEQ, t(prior_f), type = "l", lty = 1, xlim = c(28, 2), main = "Prior: Female")
matplot(DAY_SEQ, t(prior_m), type = "l", lty = 1, xlim = c(28, 2), main = "Prior: Male")
par(mfrow = c(1,1))

# --- 2.3 Fit Model ---
m_sex <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a[sex_id] + b[sex_id] * dph_c,
    a[sex_id] ~ dnorm(1, 0.5),
    b[sex_id] ~ dnorm(0, 0.1)
  ), data = d_sex, start = list(a = c(0,0), b = c(0,0))
)

# --- 2.4 Posterior Plotting ---
post_f <- link(m_sex, data = list(sex_id = rep(1, 100), dph_c = dph_seq_c_ov))
post_m <- link(m_sex, data = list(sex_id = rep(2, 100), dph_c = dph_seq_c_ov))

mu_f <- apply(post_f, 2, mean); mu_m <- apply(post_m, 2, mean)
PI_f <- apply(post_f, 2, PI);   PI_m <- apply(post_m, 2, PI)

# 1. Prepare Prediction Data for n_gaps
plot_gaps_df <- data.frame(
  day = DAY_SEQ,
  f_mu = mu_f, f_low = PI_f[1,], f_high = PI_f[2,],
  m_mu = mu_m, m_low = PI_m[1,], m_high = PI_m[2,]
)

# 2. Colors (Matching your adjusted current plot)
col_female <- "darkorange" 
col_male   <- "navyblue" 

# 3. Create the Plot
p1 <- ggplot() +
  # Ribbons
  geom_ribbon(data = plot_gaps_df, aes(x = day, ymin = f_low, ymax = f_high),  
              fill = col_female, alpha = 0.15) +
  geom_ribbon(data = plot_gaps_df, aes(x = day, ymin = m_low, ymax = m_high),  
              fill = col_male, alpha = 0.15) +
  labs(title = "A)")+
  
  # Mean Lines
  geom_line(data = plot_gaps_df, aes(x = day, y = f_mu), color = col_female, linewidth = 1.2) +
  geom_line(data = plot_gaps_df, aes(x = day, y = m_mu), color = col_male, linewidth = 1.2) +
  
  # Bubble Points (Raw n_gaps data)
  geom_count(data = E_fit_sex, aes(x = day_prior_hatch, y = n_gaps, 
                                   color = as.factor(sex_id)), alpha = 0.4) +
  
  # Manual Color Scales
  scale_color_manual(values = c("1" = col_female, "2" = col_male), guide = "none") +
  
  # Sex Labels (Using same multipliers as your neglect plot for consistency)
  annotate("text", x = 27.8, y = max(plot_gaps_df$m_mu) * 1.25, 
           label = "\u2642", size = 9, color = col_male, fontface = "bold", family = "sans") +
  annotate("text", x = 27.8, y = min(plot_gaps_df$f_mu) * 3, 
           label = "\u2640", size = 9, color = col_female, fontface = "bold", family = "sans") +
  
  # Scaling (Linear scale for n_gaps, reversed X for time)
  scale_x_reverse(breaks = seq(2, 28, by = 5)) +
  ylim(0, 15) +
  scale_size_continuous(name = "n obs", range = c(1, 7)) +
  labs(x = "Days prior hatch", y = "Number of neglect gaps") +
  
  # --- THE "PRETTY" ADJUSTMENTS (Identical to your neglect plot) ---
  theme_bw(base_size = 12, base_family = "sans") + 
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.90, 0.85),
    legend.background = element_blank(),
    
    # Matching fonts and sizes exactly
    axis.title.x = element_text(size = 13, face = "plain", family = "sans", color = "black"),
    axis.title.y = element_text(size = 13, face = "plain", family = "sans", color = "black"),
    axis.text.x = element_text(size = 11, family = "sans", color = "black"),
    axis.text.y = element_text(size = 11, family = "sans", color = "black"),
    
    legend.title = element_text(size = 11, face = "plain", family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    plot.margin = margin(10, 10, 10, 10)
  )
p1
# --- 2.5 Interpretation & Contrasts ---
post_s_samp <- extract.samples(m_sex)
diff_a <- post_s_samp$a[,1] - post_s_samp$a[,2]
gap_diff <- exp(post_s_samp$a[,1]) - exp(post_s_samp$a[,2])

cat("\n=== SEX FREQUENCY RESULTS ===\n")
cat("Female slope:", round(mean(post_s_samp$b[,1]), 4), "P(<0):", mean(post_s_samp$b[,1]<0), "\n")
cat("Male slope:  ", round(mean(post_s_samp$b[,2]), 4), "P(<0):", mean(post_s_samp$b[,2]<0), "\n")
cat("Mean Gap Diff (F-M):", round(mean(gap_diff), 2), "PI:", round(PI(gap_diff), 2), "\n")
cat("Prob Female > Male (Baseline):", mean(diff_a > 0), "\n")


# =========================================================
# 3. HURDLE NEGLECT MODEL
#    A) probability of any neglect
#    B) duration conditional on neglect > 0
# =========================================================
col_female <- "darkorange" 
col_male   <- "navyblue" 
# ---------------------------------------------------------
# 3A. PROBABILITY OF NEGLECT (includes zeros)
# ---------------------------------------------------------

D_occ <- S %>%
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
    any_neglect = ifelse(sum_neglect_sec > 0, 1, 0)
  )

d_occ <- list(
  y = D_occ$any_neglect,
  sex_id = D_occ$sex_id,
  dph_c = D_occ$dph_c
)

m_occ <- quap(
  alist(
    y ~ dbinom(1, p),
    logit(p) <- a[sex_id] + b[sex_id] * dph_c,
    a[sex_id] ~ dnorm(0, 1),
    b[sex_id] ~ dnorm(0, 0.2)
  ),
  data = d_occ,
  start = list(a = c(0,0), b = c(0,0))
)

post_f_occ <- link(m_occ, data = list(
  sex_id = rep(1, 100),
  dph_c = dph_seq_c_ov
))

post_m_occ <- link(m_occ, data = list(
  sex_id = rep(2, 100),
  dph_c = dph_seq_c_ov
))

mu_f_occ <- apply(post_f_occ, 2, mean)
mu_m_occ <- apply(post_m_occ, 2, mean)
PI_f_occ <- apply(post_f_occ, 2, PI)
PI_m_occ <- apply(post_m_occ, 2, PI)

plot_occ <- data.frame(
  day = DAY_SEQ,
  f_mu = mu_f_occ, f_low = PI_f_occ[1,], f_high = PI_f_occ[2,],
  m_mu = mu_m_occ, m_low = PI_m_occ[1,], m_high = PI_m_occ[2,]
)

p2 <- ggplot() +
  geom_ribbon(data = plot_occ,
              aes(x = day, ymin = f_low, ymax = f_high),
              fill = col_female, alpha = 0.15) +
  geom_ribbon(data = plot_occ,
              aes(x = day, ymin = m_low, ymax = m_high),
              fill = col_male, alpha = 0.15) +
  labs(title = "C)") +
  
  geom_line(data = plot_occ,
            aes(x = day, y = f_mu),
            color = col_female, linewidth = 1.2) +
  geom_line(data = plot_occ,
            aes(x = day, y = m_mu),
            color = col_male, linewidth = 1.2) +
  
  geom_count(data = D_occ,
             aes(x = day_prior_hatch,
                 y = any_neglect,
                 color = as.factor(sex_id)),
             alpha = 0.4) +
  
  scale_color_manual(values = c("1" = col_female, "2" = col_male),
                     guide = "none") +
  
  annotate("text", x = 28.8, y = 0.88,
           label = "\u2642", size = 9,
           color = col_male, fontface = "bold") +
  annotate("text", x = 28.8, y = 0.75,
           label = "\u2640", size = 9,
           color = col_female, fontface = "bold") +
  
  scale_x_reverse(breaks = seq(2, 28, by = 5)) +
  ylim(0, 1) +
  scale_size_continuous(name = "n obs", range = c(1, 6)) +
  labs(x = "Days prior hatch", y = "Probability of neglect") +
  
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.90, 0.81),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
p2

# --- 3A.5 Interpretation: Probability of Any Neglect ---
post_occ_samp <- extract.samples(m_occ)

# female - male baseline contrast on logit scale
diff_a_occ <- post_occ_samp$a[,1] - post_occ_samp$a[,2]

# convert baseline to probability scale at centered day
prob_f <- inv_logit(post_occ_samp$a[,1])
prob_m <- inv_logit(post_occ_samp$a[,2])
prob_diff <- prob_f - prob_m

cat("\n=== PROBABILITY OF NEGLECT RESULTS ===\n")
cat(
  "Female slope:",
  round(mean(post_occ_samp$b[,1]), 4),
  "PI:", round(PI(post_occ_samp$b[,1]), 4),
  "P(<0):", round(mean(post_occ_samp$b[,1] < 0), 3),
  "\n"
)

cat(
  "Male slope:  ",
  round(mean(post_occ_samp$b[,2]), 4),
  "PI:", round(PI(post_occ_samp$b[,2]), 4),
  "P(<0):", round(mean(post_occ_samp$b[,2] < 0), 3),
  "\n"
)

cat(
  "Mean Probability Diff (F-M):",
  round(mean(prob_diff), 3),
  "PI:", round(PI(prob_diff), 3),
  "\n"
)

cat(
  "Prob Female > Male (Baseline):",
  round(mean(diff_a_occ > 0), 3),
  "\n"
)

# ---------------------------------------------------------
# 3B. DURATION GIVEN NEGLECT > 0
# ---------------------------------------------------------

D_dur <- S %>%
  filter(
    !is.na(sum_neglect_sec),
    sum_neglect_sec > 0,
    !is.na(day_prior_hatch),
    !is.na(sx),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    sex_id = ifelse(sx == "f", 1, 2),
    dph_c = day_prior_hatch - mean(day_prior_hatch),
    neglect_log = log(sum_neglect_sec)
  )

d_dur <- list(
  y = D_dur$neglect_log,
  sex_id = D_dur$sex_id,
  dph_c = D_dur$dph_c
)

m_dur <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a[sex_id] + b[sex_id] * dph_c,
    a[sex_id] ~ dnorm(6, 0.5),
    b[sex_id] ~ dnorm(0, 0.1),
    sigma ~ dexp(1)
  ),
  data = d_dur,
  start = list(a = c(6,6), b = c(0,0), sigma = 1)
)

post_f_dur <- exp(link(m_dur, data = list(
  sex_id = rep(1, 100),
  dph_c = dph_seq_c_ov
)))

post_m_dur <- exp(link(m_dur, data = list(
  sex_id = rep(2, 100),
  dph_c = dph_seq_c_ov
)))

mu_f_dur <- apply(post_f_dur, 2, mean)
mu_m_dur <- apply(post_m_dur, 2, mean)
PI_f_dur <- apply(post_f_dur, 2, PI)
PI_m_dur <- apply(post_m_dur, 2, PI)

plot_df <- data.frame(
  day = DAY_SEQ,
  f_mu = mu_f_dur, f_low = PI_f_dur[1,], f_high = PI_f_dur[2,],
  m_mu = mu_m_dur, m_low = PI_m_dur[1,], m_high = PI_m_dur[2,]
)

p3 <- ggplot() +
  # 1. Uncertainty Ribbons
  geom_ribbon(data = plot_df,
              aes(x = day, ymin = f_low, ymax = f_high),
              fill = col_female, alpha = 0.15) +
  geom_ribbon(data = plot_df,
              aes(x = day, ymin = m_low, ymax = m_high),
              fill = col_male, alpha = 0.15) +
  
  # 2. Regression Lines
  geom_line(data = plot_df,
            aes(x = day, y = f_mu),
            color = col_female, linewidth = 1.2) +
  geom_line(data = plot_df,
            aes(x = day, y = m_mu),
            color = col_male, linewidth = 1.2) +
  
  # 3. Observation Bubbles 
  # Using round() ensures that close values overlap to create "n obs" > 1
  geom_point(data = D_dur,
             aes(x = day_prior_hatch,
                 y = sum_neglect_sec,
                 color = as.factor(sex_id)),
             alpha = 0.4, 
             size = 2) + #
  
  # 4. Sex Annotations (Placed at the start of the lines)
  annotate("text", x = 27.8, y = plot_df$m_mu[100] * 3, 
           label = "\u2642", size = 10, color = col_male, fontface = "bold", family = "sans") +
  annotate("text", x = 27.8, y = plot_df$f_mu[100] * 0.3, 
           label = "\u2640", size = 10, color = col_female, fontface = "bold", family = "sans") +
  
  # 5. Scales and Labels
  scale_color_manual(values = c("1" = col_female, "2" = col_male), guide = "none") +
  scale_y_log10(labels = scales::label_number(),
                breaks = scales::trans_breaks("log10", function(x) 10^x)) +
  scale_x_reverse(breaks = seq(2, 28, by = 5)) +
  
  # 6. Size Legend (Fixed to whole numbers 1 and 2 only)
  scale_size_continuous(name = "n obs", 
                        range = c(1.5, 4), 
                        breaks = c(1, 2),
                        labels = c("1", "2")) +
  
  labs(title = "B)", 
       x = "Days prior hatch",
       y = "Duration of neglect (sec, log scale)") +
  
  # 7. Theme and Font Adjustments
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.90, 0.85),      # Nudged down slightly to avoid border
    legend.background = element_blank(), 
    legend.key = element_blank(),        # Transparent legend keys
    
    # Matching all text to Arial/Sans
    axis.title.x = element_text(size = 13, family = "sans"),
    axis.title.y = element_text(size = 13, family = "sans"),
    axis.text.x = element_text(size = 11, color = "black", family = "sans"),
    axis.text.y = element_text(size = 11, color = "black", family = "sans"),
    
    legend.title = element_text(size = 11, family = "sans"),
    legend.text = element_text(size = 10, family = "sans"),
    plot.margin = margin(10, 10, 10, 10)
  )

# View plot
p3


# --- 3B.5 Interpretation: Duration Given Neglect ---
post_dur_samp <- extract.samples(m_dur)

# baseline contrast on log scale
diff_a_dur <- post_dur_samp$a[,1] - post_dur_samp$a[,2]

# convert to seconds scale
dur_f <- exp(post_dur_samp$a[,1])
dur_m <- exp(post_dur_samp$a[,2])
dur_diff <- dur_f - dur_m

cat("\n=== DURATION GIVEN NEGLECT RESULTS ===\n")
cat(
  "Female slope:",
  round(mean(post_dur_samp$b[,1]), 4),
  "PI:", round(PI(post_dur_samp$b[,1]), 4),
  "P(<0):", round(mean(post_dur_samp$b[,1] < 0), 3),
  "\n"
)

cat(
  "Male slope:  ",
  round(mean(post_dur_samp$b[,2]), 4),
  "PI:", round(PI(post_dur_samp$b[,2]), 4),
  "P(<0):", round(mean(post_dur_samp$b[,2] < 0), 3),
  "\n"
)

cat(
  "Mean Duration Diff (F-M sec):",
  round(mean(dur_diff), 1),
  "PI:", round(PI(dur_diff), 1),
  "\n"
)

cat(
  "Prob Female > Male (Baseline):",
  round(mean(diff_a_dur > 0), 3),
  "\n"
)




#save plots
library(patchwork)

#combined_plot <- "
#AB
#C#
#"

#final_fig <- (p1 + p3 + p2) +
#  plot_layout(design = combined_plot)

#final_fig
# 5. Save as one high-res file
#ggsave("Figure_1.png", final_fig, width = 14, height = 12, dpi = 600)

ggsave(
  "Figure_1A_gap_frequency.png",
  p1,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  "Figure_1B_neglect_duration.png",
  p2,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  "Figure_1C_probability_neglect.png",
  p3,
  width = 7,
  height = 6,
  dpi = 600
)
