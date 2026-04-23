# PURPOSE: model n_gaps and sum_gaps towards hatching
rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)
library(splines)

S <- readRDS("EDA_sex_neglect_wide.RDS")
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")
G <- readRDS("./EDA_sex_gaps_table.RDS")

#plot to see the duration of absences
G$duration_category <- cut(
  G$dur, 
  breaks = c(0, 300, 600, 900, 1800, 2700, 3600, 7200, Inf), 
  labels = c("0-5 min", "5-10 min", "10-15 min", "15-30 min", "30-45 min", "45-60 min", "1-2 hours", "2+ hours"),
  right = FALSE
)

library(ggplot2)

p <- ggplot(G, aes(x = duration_category)) +
  geom_bar(fill = "orange") +
  facet_wrap(~responsible_sex)+
  labs(title = "Duration of absences",
       x = "Duration", 
       y = "No. of absences") +
  theme_classic()
p


# 1. Combine unique nest-season combinations from both S and E
# This ensures we don't miss a nest that might be in one dataset but not the other
all_nest_records <- bind_rows(
  S %>% select(nest, season),
  #E %>% select(nest, season)
) %>% 
  distinct()

# 2. Count total unique nests
total_unique_nests <- n_distinct(all_nest_records$nest)

# 3. Identify nests followed for more than one season
multi_season_summary <- all_nest_records %>%
  group_by(nest) %>%
  summarise(n_seasons = n_distinct(season)) %>%
  filter(n_seasons > 1)

session_counts_per_year <- E %>%
  group_by(season) %>%
  summarise(
    total_sampling_events = n()
  )

# 4. Final counts
n_multi_season <- nrow(multi_season_summary)

# Print results
cat("Total unique nests across all seasons:", total_unique_nests, "\n")
cat("Number of nests followed for >1 season:", n_multi_season, "\n")

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

G_summary <- G %>%
  filter(dur >= 300) %>%
  group_by(season, session, nest, responsible_sex) %>%
  summarise(
    n_gaps = n(),
    sum_neglect_sec = sum(dur),
    mean_neglect_sec = mean(dur),
    .groups = "drop"
  ) %>%
  rename(sx = responsible_sex)

# Replace the neglect columns in S with the filtered versions
S_filt <- S %>%
  select(-n_gaps, -sum_neglect_sec, -mean_neglect_sec) %>%
  left_join(G_summary, by = c("season", "session", "nest", "sx")) %>%
  mutate(
    sum_neglect_sec = replace_na(sum_neglect_sec, 0),
    mean_neglect_sec = replace_na(mean_neglect_sec, 0),
    n_gaps = replace_na(n_gaps, 0L)
  )


# --- 2.1 Data Preparation ---
E_fit_sex <- S_filt %>%
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


# Symulujemy obserwacje (n_gaps) dla każdego punktu w oryginalnych danych
# m_sex to Twój dopasowany model quap
post_samples <- sim(m_sex, data = d_sex)
# Bierzemy losową próbkę symulacji (np. pierwszą), aby pokazać realistyczny rozrzut
E_fit_sex$post_n_gaps <- post_samples[1, ]

# 2. Tworzenie wykresu
p1 <- ggplot() +
  # --- 1. Wstęgi Niepewności ---
  geom_ribbon(data = plot_gaps_df, aes(x = day, ymin = f_low, ymax = f_high),  
              fill = col_female, alpha = 0.15) +
  geom_ribbon(data = plot_gaps_df, aes(x = day, ymin = m_low, ymax = m_high),  
              fill = col_male, alpha = 0.15) +
  
  # --- 2. Linie Trendu (Średnie) ---
  geom_line(data = plot_gaps_df, aes(x = day, y = f_mu), 
            color = col_female, linewidth = 1.2) +
  geom_line(data = plot_gaps_df, aes(x = day, y = m_mu), 
            color = col_male, linewidth = 1.2) +
  
  # --- 3. Bąbelki Posterior Predictive ---
  geom_count(data = E_fit_sex, 
             aes(x = day_prior_hatch, y = post_n_gaps, color = as.factor(sex_id)), 
             alpha = 0.3,
             position = position_jitter(width = 0.15, height = 0.05)) +
  
  # --- 4. Przywrócone Adnotacje Płci ---
  annotate("text", x = 27.8, y = max(plot_gaps_df$m_mu) * 1.4, 
           label = "\u2642", size = 10, color = col_male, fontface = "bold", family = "sans") +
  annotate("text", x = 27.8, y = 0.7, 
           label = "\u2640", size = 10, color = col_female, fontface = "bold", family = "sans") +
  
  # --- 5. Skalowanie i Etykiety ---
  scale_color_manual(values = c("1" = col_female, "2" = col_male), guide = "none") +
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(
    # Zoom na istotny zakres (model Poissona rzadko wychodzi poza 4-5)
    limits = c(-0.3, 2.2), 
    breaks = seq(0, 5, by = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_size_continuous(name = "n obs", range = c(1, 7)) +
  labs(title = "A)", x = "Days prior hatch", y = "Number of neglect gaps") +
  
  # --- 6. Motyw Graficzny ---
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.90, 0.85),
    legend.background = element_blank(),
    axis.title = element_text(size = 13, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
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
cat("Mean Exit Diff (F-M):", round(mean(gap_diff), 2), "PI:", round(PI(gap_diff), 2), "\n")
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

head(G)
head(S)

D_occ <- S_filt %>%
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

post_samples_occ <- sim(m_occ, data = d_occ)

# Wybieramy jedną próbkę symulacji (0 lub 1) dla każdego wiersza danych
D_occ$post_any_neglect <- post_samples_occ[1, ]

# --- 2. Tworzenie wykresu p2 ---
p2 <- ggplot() +
  # Uncertainty Ribbons (95% PI dla średniej)
  geom_ribbon(data = plot_occ,
              aes(x = day, ymin = f_low, ymax = f_high),
              fill = col_female, alpha = 0.15) +
  geom_ribbon(data = plot_occ,
              aes(x = day, ymin = m_low, ymax = m_high),
              fill = col_male, alpha = 0.15) +
  
  # Mean Lines
  geom_line(data = plot_occ,
            aes(x = day, y = f_mu),
            color = col_female, linewidth = 1.2) +
  geom_line(data = plot_occ,
            aes(x = day, y = m_mu),
            color = col_male, linewidth = 1.2) +
  
  # Bubble Points - TERAZ Z POSTERIOR PREDICTIVE (0 lub 1)
  geom_count(data = D_occ,
             aes(x = day_prior_hatch,
                 y = post_any_neglect, # Zmienione na dane symulowane
                 color = as.factor(sex_id)),
             alpha = 0.3,
             # Dodajemy jitter, żeby bąbelki na 0 i 1 nie były "płaskimi" plackami
             position = position_jitter(width = 0.15, height = 0.03)) +
  
  # Manual Color Scales
  scale_color_manual(values = c("1" = col_female, "2" = col_male),
                     guide = "none") +
  
  # Sex Labels
  annotate("text", x = 27.8, y = 0.7,
           label = "\u2642", size = 10,
           color = col_male, fontface = "bold") +
  annotate("text", x = 27.8, y = 0.20,
           label = "\u2640", size = 10,
           color = col_female, fontface = "bold") +
  
  # Scaling and Labels
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  # Dodajemy mały margines w coord_cartesian, żeby jitter nie "uciekał" poza oś
  coord_cartesian(ylim = c(-0.1, 1.1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_size_continuous(name = "n obs", range = c(1, 6)) +
  labs(title = "C)", x = "Days prior hatch", y = "Probability of neglect") +
  
  # Theme
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
post_samples_log_dur <- sim(m_dur, data = d_dur)


D_dur$post_sum_neglect_sec <- exp(post_samples_log_dur[1, ])

p3 <- ggplot() +
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
  
  # 3. Bubble Points - TERAZ Z POSTERIOR PREDICTIVE
  # geom_count automatycznie pogrupuje podobne wartości symulowane
  geom_count(data = D_dur,
             aes(x = day_prior_hatch,
                 y = post_sum_neglect_sec, 
                 color = as.factor(sex_id)),
             alpha = 0.3) + 
  
  # 4. Sex Annotations
  annotate("text", x = 27.8, y = plot_df$m_mu[100] * 2.4, 
           label = "\u2642", size = 10, color = col_male, fontface = "bold", family = "sans") +
  annotate("text", x = 27.8, y = plot_df$f_mu[100] * 0.4, 
           label = "\u2640", size = 10, color = col_female, fontface = "bold", family = "sans") +

# 5. Scales and Labels
scale_color_manual(values = c("1" = col_female, "2" = col_male), guide = "none") +
  
  scale_size_continuous(range = c(0.5, 2), guide = "none") + 
  
  scale_y_log10(
    breaks = c(100, 1000, 10000, 100000), 
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  scale_x_reverse(breaks = seq(0, 30, by = 5)) +
  labs(title = "B)", 
       x = "Days prior hatch",
       y = "Duration of neglect events (sec, log scale)") +
  
  # 7. Theme
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none", # Najpewniejszy sposób na usunięcie wszystkich legend
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11, color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

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




# Each letter on a new line stacks them vertically
combined_plot <- "
A
B
C
"

library(patchwork)
# Assign your plots to the design
# Note: Ensure p1 is your 'A', p3 is 'B', and p2 is 'C' based on your layout
final_fig <- (p1 / p3 / p2) + 
  plot_layout(design = combined_plot)

# Alternatively, the shorthand for vertical stacking is simply:
# final_fig <- p1 / p3 / p2

# Save the figure - notice I swapped width and height 
# to better suit a vertical orientation
p1
ggsave(
  "Figure_1A_n_gaps.png",
  p1,
  width = 7,
  height = 6,
  dpi = 600
)

ggsave(
  "Figure_1C_probability_neglect.png",
  p2,
  width = 7,
  height = 6,
  dpi = 600
)
p3
ggsave(
  "Figure_1B_neglect_duration.png",
  p3,
  width = 7,
  height = 6,
  dpi = 600
)






# 1. Function to strip existing text annotations from a plot
remove_annotations <- function(p) {
  p$layers <- Filter(function(x) !inherits(x$geom, "GeomText"), p$layers)
  return(p)
}

# 2. Clean the plots (this removes the 'ghost' labels)
p1_clean <- remove_annotations(p1)
p2_clean <- remove_annotations(p2)
p3_clean <- remove_annotations(p3)

# 3. Re-define the style (No grid, tiny legend)
final_style <- theme_bw(base_size = 9) + 
  theme(
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.2, "cm")
  )

# 4. Re-assemble with SINGLE set of labels
pA_final <- p1_clean + 
  annotate("text", x = 27.8, y = 1.8, label = "\u2642", size = 5, color = col_male) +
  annotate("text", x = 27.8, y = 0.7, label = "\u2640", size = 5, color = col_female) +
  scale_size_continuous(name = "n obs", range = c(0.5, 2.5)) +
  final_style + theme(legend.position = c(0.88, 0.85), axis.title.x = element_blank())

pB_final <- p3_clean + 
  annotate("text", x = 27.8, y = plot_df$m_mu[100] * 3, label = "\u2642", size = 5, color = col_male) +
  annotate("text", x = 27.8, y = plot_df$f_mu[100] * 0.3, label = "\u2640", size = 5, color = col_female) +
  scale_size_continuous(range = c(1, 1)) +
  final_style + theme(legend.position = "none", axis.title.x = element_blank())

pC_final <- p2_clean + 
  annotate("text", x = 27.8, y = 0.7, label = "\u2642", size = 5, color = col_male) +
  annotate("text", x = 27.8, y = 0.20, label = "\u2640", size = 5, color = col_female) +
  scale_size_continuous(name = "n obs", range = c(0.5, 2.5)) +
  final_style + theme(legend.position = c(0.88, 0.72)) # Top right as requested

# 5. Combine and Save
final_fig <- pA_final / pB_final / pC_final
ggsave("Figure_1_patched.png", final_fig, width = 3.5, height = 10, units = "in", dpi = 600)
