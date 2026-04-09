# ==============================================================================
# PURPOSE: Model diurnal exit patterns for nest level and both sexes
# ORGANIZED BY: Data Prep -> Model -> Priors -> Stats -> Plots
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. INITIAL SETUP & LIBRARIES
# ------------------------------------------------------------------------------
library(rethinking)
library(tidyverse)
library(lubridate)
library(splines)
library(dplyr)

# Global Visual Settings
stage_names  <- c("early", "mid", "late")
stage_colors <- c("#3B6D11", "#185FA5", "#993556")
sex_labels   <- c("Female", "Male")
sex_col      <- c("darkorange", "navyblue")

# ------------------------------------------------------------------------------
# 2. NEST-LEVEL DATA PREPARATION
# ------------------------------------------------------------------------------
E <- readRDS("./EDA_exit_hourly.RDS")

E <- E %>%
  mutate(
    incubation_stage = case_when(
      day_prior_hatch >= 21 & day_prior_hatch <= 30 ~ "early",
      day_prior_hatch >= 11 & day_prior_hatch <= 20 ~ "mid",
      day_prior_hatch >= 0  & day_prior_hatch <= 10 ~ "late",
      TRUE ~ NA_character_
    ),
    stage_id = as.integer(factor(incubation_stage, levels = stage_names)),
    pair_id  = as.integer(factor(Pair_status, levels = c("New", "Old"))),
    stage_pair_id = as.integer(factor(
      paste(incubation_stage, Pair_status, sep = "_"),
      levels = c("early_New", "early_Old", "mid_New", "mid_Old", "late_New", "late_Old")
    ))
  )

# Define Spline Basis (Hour sequence 0-23)
hour_seq <- 0:23
basis    <- bs(hour_seq, df = 6, degree = 3, intercept = FALSE)
basis_cols <- as.data.frame(basis)
colnames(basis_cols) <- paste0("B", 1:6)
basis_list <- as.list(basis_cols)

E_basis <- bind_cols(E, basis_cols[match(E$hour, hour_seq), ])

# ------------------------------------------------------------------------------
# 3. MODEL: OVERALL DAILY EXIT SPLINE
# ------------------------------------------------------------------------------
d_m2 <- c(list(n_gaps_per_hour = E_basis$n_gaps_per_hour), as.list(E_basis[, paste0("B",1:6)]))

m2 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a + w1*B1 + w2*B2 + w3*B3 + w4*B4 + w5*B5 + w6*B6,
    a ~ dnorm(0, 0.8),
    w1 ~ dnorm(0, 0.5), w2 ~ dnorm(0, 0.5), w3 ~ dnorm(0, 0.5),
    w4 ~ dnorm(0, 0.5), w5 ~ dnorm(0, 0.5), w6 ~ dnorm(0, 0.5)
  ),
  data = d_m2
)

# --- PRIOR PREDICTIVE SIMULATION (M2) ---
set.seed(10)
plot(NULL, xlim=c(0,23), ylim=c(0,5), xlab="Hour", ylab="Expected exits", xaxt="n", main="M2 Prior Predictive Check")
axis(1, at=0:23)
for(i in 1:50){
  a_sim <- rnorm(1, 0, 0.8)
  ws_sim <- rnorm(6, 0, 0.5)
  mu_sim <- exp(a_sim + as.matrix(basis_cols) %*% ws_sim)
  lines(hour_seq, mu_sim, col=col.alpha("black", 0.25))
}

# --- POSTERIOR FIT PLOT (M2) ---
pred <- link(m2, data = basis_list)
mu   <- apply(pred, 2, mean)
PI_mu <- apply(pred, 2, PI, prob=0.89)

Eb <- E %>% group_by(hour) %>% summarise(n_gaps_per_hour = mean(n_gaps_per_hour), .groups="drop")
plot(Eb$hour, Eb$n_gaps_per_hour, pch=16, xlab="Hour", ylab="Mean exits", main="M2 Posterior Fit")
shade(PI_mu, hour_seq)
lines(hour_seq, mu, lwd=2)

# ------------------------------------------------------------------------------
# 4. MODEL: STAGE-SPECIFIC DIURNAL SPLINES
# ------------------------------------------------------------------------------
E_fit3 <- E_basis %>% filter(!is.na(stage_id))
d_m3 <- c(list(n_gaps_per_hour = E_fit3$n_gaps_per_hour, stage_id = E_fit3$stage_id),
          as.list(E_fit3[, paste0("B",1:6)]))

m3 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[stage_id] + w1[stage_id]*B1 + w2[stage_id]*B2 + w3[stage_id]*B3 +
      w4[stage_id]*B4 + w5[stage_id]*B5 + w6[stage_id]*B6,
    a[stage_id] ~ dnorm(0, 0.8),
    w1[stage_id] ~ dnorm(0, 0.5), w2[stage_id] ~ dnorm(0, 0.5), w3[stage_id] ~ dnorm(0, 0.5),
    w4[stage_id] ~ dnorm(0, 0.5), w5[stage_id] ~ dnorm(0, 0.5), w6[stage_id] ~ dnorm(0, 0.5)
  ),
  data = d_m3,
  start = list(a=rep(0,3), w1=rep(0,3), w2=rep(0,3), w3=rep(0,3), w4=rep(0,3), w5=rep(0,3), w6=rep(0,3))
)

# --- POSTERIOR PLOTS & STATS (M3) ---
Eg <- E_fit3 %>% group_by(stage_id, hour) %>% summarise(mean_gaps=mean(n_gaps_per_hour), .groups="drop")
plot(NULL, xlim=c(0,23), ylim=c(0,0.5), xlab="Hour", ylab="Expected exits", xaxt="n", main="M3: Stage-Specific Splines")
axis(1, at=0:23)

mu_list_m3 <- list()
for(s in 1:3){
  pred_s <- link(m3, data = c(list(stage_id = rep(s,24)), basis_list))
  mu_s <- apply(pred_s, 2, mean)
  PI_s <- apply(pred_s, 2, PI, prob=0.89)
  mu_list_m3[[stage_names[s]]] <- apply(pred_s, 1, mean)
  shade(PI_s, hour_seq, col=col.alpha(stage_colors[s], 0.15))
  lines(hour_seq, mu_s, lwd=2, col=stage_colors[s])
  pts <- Eg %>% filter(stage_id == s)
  points(pts$hour, pts$mean_gaps, pch=16, col=stage_colors[s], cex=0.8)
}
legend("topright", legend=stage_names, col=stage_colors, lwd=2, pch=16)

# --- POSTERIOR DENSITIES & CONTRASTS (M3) ---
xrange <- range(unlist(mu_list_m3))
plot(NULL, xlim=xrange, ylim=c(0,120), xlab="Expected exits/hour", ylab="Density", main="Posterior distributions by stage")
for(i in 1:3) dens(mu_list_m3[[i]], col=stage_colors[i], lwd=2, add=TRUE)

cat("\n--- M3 STAGE CONTRASTS ---\n")
diff_early_mid <- mu_list_m3[["early"]] - mu_list_m3[["mid"]]
diff_early_late <- mu_list_m3[["early"]] - mu_list_m3[["late"]]
diff_mid_late <- mu_list_m3[["mid"]] - mu_list_m3[["late"]]
cat("Early vs Mid: ", round(mean(diff_early_mid),3), PI(diff_early_mid), round(mean(diff_early_mid>0),2), "\n")
cat("Early vs Late:", round(mean(diff_early_late),3), PI(diff_early_late), round(mean(diff_early_late>0),2), "\n")
cat("Mid vs Late:  ", round(mean(diff_mid_late),3), PI(diff_mid_late), round(mean(diff_mid_late>0),2), "\n")

# --- PEAK HOUR SUMMARIES (M3) ---
for(s in 1:3){
  pred_s <- link(m3, data = c(list(stage_id=rep(s,24)), basis_list))
  peak_hour <- apply(pred_s, 1, which.max) - 1
  cat(stage_names[s], "peak:", round(mean(peak_hour),2), PI(peak_hour), "\n")
}

# ------------------------------------------------------------------------------
# 5. SEX-SPECIFIC DATA PREPARATION
# ------------------------------------------------------------------------------
S <- readRDS("EDA_sex_neglect_wide.RDS")
G <- readRDS("EDA_sex_gaps_hourly.RDS")

sex_hourly_raw <- G %>%
  filter(!is.na(responsible_sex)) %>%
  group_by(season, session, nest, sx = responsible_sex, hour) %>%
  summarise(n_gaps_per_hour = n(), .groups = "drop")

meta <- S %>%
  mutate(
    season = as.character(season),
    incubation_stage = case_when(
      day_prior_hatch >= 21 & day_prior_hatch <= 30 ~ "early",
      day_prior_hatch >= 11 & day_prior_hatch <= 20 ~ "mid",
      day_prior_hatch >= 0  & day_prior_hatch <= 10 ~ "late",
      TRUE ~ NA_character_
    ),
    stage_id = as.integer(factor(incubation_stage, levels = stage_names)),
    sex_id   = as.integer(factor(sx, levels = c("f", "m"))),
    stage_sex_id = as.integer(factor(paste(incubation_stage, sx, sep = "_"),
                                     levels = c("early_f", "early_m", "mid_f", "mid_m", "late_f", "late_m")))
  ) %>%
  select(season, session, nest, sx, day_prior_hatch, incubation_stage, stage_id, sex_id, stage_sex_id)

E_sex <- meta %>% crossing(hour = 0:23) %>%
  left_join(sex_hourly_raw, by = c("season", "session", "nest", "sx", "hour")) %>%
  mutate(n_gaps_per_hour = replace_na(n_gaps_per_hour, 0))

E_by_sex <- E_sex %>% group_by(sex_id, hour) %>% summarise(mean_gaps = mean(n_gaps_per_hour), .groups = "drop")
E_stage_sex <- E_sex %>% group_by(stage_sex_id, hour) %>% summarise(mean_gaps = mean(n_gaps_per_hour), .groups = "drop")

# Re-scaling Splines for Sex Data
E_fit_sex <- E_sex %>% filter(!is.na(sex_id))
B_sex <- bs(E_fit_sex$hour, df = 6, degree = 3, intercept = TRUE)
B_pred_sex <- bs(hour_seq, knots = attr(B_sex, "knots"), Boundary.knots = attr(B_sex, "Boundary.knots"), degree = 3, intercept = TRUE)
B_scaled_sex <- scale(B_sex)
B_pred_scaled_sex <- scale(B_pred_sex, center = attr(B_scaled_sex, "scaled:center"), scale = attr(B_scaled_sex, "scaled:scale"))
basis_cols_sex <- setNames(as.data.frame(B_scaled_sex), paste0("B", 1:ncol(B_sex)))
basis_list_sex <- as.list(setNames(as.data.frame(B_pred_scaled_sex), paste0("B", 1:ncol(B_sex))))

# ------------------------------------------------------------------------------
# 6. MODEL 1 (SEX DATA): SEX ONLY
# ------------------------------------------------------------------------------
d_m1_sex <- c(list(n_gaps_per_hour = E_fit_sex$n_gaps_per_hour, sex_id = E_fit_sex$sex_id), as.list(basis_cols_sex))

m1_sex <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[sex_id] + w1[sex_id]*B1 + w2[sex_id]*B2 + w3[sex_id]*B3 + 
      w4[sex_id]*B4 + w5[sex_id]*B5 + w6[sex_id]*B6,
    a[sex_id] ~ dnorm(0, 1),
    w1[sex_id] ~ dnorm(0, 0.4), w2[sex_id] ~ dnorm(0, 0.4), w3[sex_id] ~ dnorm(0, 0.4),
    w4[sex_id] ~ dnorm(0, 0.4), w5[sex_id] ~ dnorm(0, 0.4), w6[sex_id] ~ dnorm(0, 0.4)
  ),
  data = d_m1_sex, start = list(a=rep(0,2), w1=rep(0,2), w2=rep(0,2), w3=rep(0,2), w4=rep(0,2), w5=rep(0,2), w6=rep(0,2))
)

hour_seq <- 0:23

# Predictions for Females (sex_id = 1)
post_f <- link(m1_sex, data = c(list(sex_id = rep(1, 24)), basis_list_sex))
mu_f   <- apply(post_f, 2, mean)
PI_f   <- apply(post_f, 2, PI, prob = 0.89)

# Predictions for Males (sex_id = 2)
post_m <- link(m1_sex, data = c(list(sex_id = rep(2, 24)), basis_list_sex))
mu_m   <- apply(post_m, 2, mean)
PI_m   <- apply(post_m, 2, PI, prob = 0.89)

# Combine into a clean dataframe for ggplot
plot_hour_df <- data.frame(
  hour   = hour_seq,
  f_mu   = mu_f, 
  f_low  = PI_f[1,], 
  f_high = PI_f[2,],
  m_mu   = mu_m, 
  m_low  = PI_m[1,], 
  m_high = PI_m[2,]
)

# --- 2. DEFINE COLORS ---
col_female <- "darkorange"
col_male   <- "navyblue"

# --- 3. CREATE THE PLOT ---
p_spline <- ggplot(plot_hour_df, aes(x = hour)) +
  # Uncertainty Ribbons (89% PI)
  geom_ribbon(aes(ymin = f_low, ymax = f_high), 
              fill = col_female, alpha = 0.15) +
  geom_ribbon(aes(ymin = m_low, ymax = m_high), 
              fill = col_male, alpha = 0.15) +
  
  # Mean Lines (Spline trends)
  geom_line(aes(y = f_mu), color = col_female, linewidth = 1.2) +
  geom_line(aes(y = m_mu), color = col_male, linewidth = 1.2) +
  
  # Sex Labels in the expanded left gutter
  # Y-positions are set to align near the start of the lines (~0.2 range)
  annotate("text", x = 0.5, y = 0.17, 
           label = "\u2642", size = 10, color = col_male, fontface = "bold", family = "sans") +
  annotate("text", x = 0.5, y = 0.07, 
           label = "\u2640", size = 10, color = col_female, fontface = "bold", family = "sans") +
  
  # X-Axis: Expanded to -2 to make room for labels
  scale_x_continuous(breaks = seq(0, 23, by = 2), limits = c(0, 23)) +
  
  # Y-Axis: Fixed to model scale (0 to 0.3)
  coord_cartesian(ylim = c(0, 0.25), clip = "off") +
  
  labs(x = "Hour of day", 
       y = "Expected exits")+
  
  # Final Paper Theme Styling
  theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 13, family = "sans"),
    axis.title.y = element_text(size = 13, family = "sans"),
    axis.text.x  = element_text(size = 11, color = "black", family = "sans"),
    axis.text.y  = element_text(size = 11, color = "black", family = "sans"),
    plot.margin  = margin(10, 10, 10, 10) # Extra left margin for the sex icons
  )

# Display the plot
print(p_spline)

# Sex Rates & Peaks
mu_f_m1 <- link(m1_sex, data = c(list(sex_id = rep(1, 24)), basis_list_sex))
mu_m_m1 <- link(m1_sex, data = c(list(sex_id = rep(2, 24)), basis_list_sex))
avg_rate_f_m1 <- apply(mu_f_m1, 1, mean); avg_rate_m_m1 <- apply(mu_m_m1, 1, mean); diff_fm_m1 <- avg_rate_f_m1 - avg_rate_m_m1
cat("\n--- M1 SEX SUMMARY ---\n")
cat("Female mean rate:", round(mean(avg_rate_f_m1), 3), "Male mean rate:", round(mean(avg_rate_m_m1), 3), "Diff:", round(mean(diff_fm_m1), 3), "\n")
for (s in 1:2) {
  pred_s <- link(m1_sex, data = c(list(sex_id = rep(s, 24)), basis_list_sex))
  peak_h <- apply(pred_s, 1, which.max) - 1
  cat(sex_labels[s], "Peak Hour:", round(mean(peak_h), 1), "\n")
}

ggsave(
  "Figure_2_diurnal_exit_pattern.png",
  p_spline,
  width = 7,
  height = 6,
  dpi = 600
)


# ------------------------------------------------------------------------------
# 7. MODEL 2 (SEX DATA): SEX X STAGE
# ------------------------------------------------------------------------------
E_fit2_sex <- E_sex %>% filter(!is.na(stage_sex_id))
d_m2_sex <- c(list(n_gaps_per_hour = E_fit2_sex$n_gaps_per_hour, stage_sex_id = E_fit2_sex$stage_sex_id), as.list(basis_cols_sex))

m2_sex <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[stage_sex_id] + w1[stage_sex_id]*B1 + w2[stage_sex_id]*B2 + w3[stage_sex_id]*B3 +
      w4[stage_sex_id]*B4 + w5[stage_sex_id]*B5 + w6[stage_sex_id]*B6,
    a[stage_sex_id] ~ dnorm(0, 1),
    w1[stage_sex_id] ~ dnorm(0, 0.4), w2[stage_sex_id] ~ dnorm(0, 0.4), w3[stage_sex_id] ~ dnorm(0, 0.4),
    w4[stage_sex_id] ~ dnorm(0, 0.4), w5[stage_sex_id] ~ dnorm(0, 0.4), w6[stage_sex_id] ~ dnorm(0, 0.4)
  ),
  data = d_m2_sex, start = list(a=rep(0,6), w1=rep(0,6), w2=rep(0,6), w3=rep(0,6), w4=rep(0,6), w5=rep(0,6), w6=rep(0,6))
)

# --- PLOT M2: INTERACTION ---
combos <- c("Early F", "Early M", "Mid F", "Mid M", "Late F", "Late M")
stage_cols_6 <- rep(stage_colors, each = 2); sex_lty_6 <- rep(c(1, 2), 3)
plot(NULL, xlim = c(0, 23), ylim = c(0, 0.5), xlab = "Hour", ylab = "Expected exits", xaxt = "n", main = "Daily exit pattern — Sex x Stage (M2)")
axis(1, at = 0:23)
for (s in 1:6) {
  pred_s <- link(m2_sex, data = c(list(stage_sex_id = rep(s, 24)), basis_list_sex))
  mu_s <- apply(pred_s, 2, mean); PI_s <- apply(pred_s, 2, PI, prob = 0.89)
  shade(PI_s, hour_seq, col = col.alpha(stage_cols_6[s], 0.15))
  lines(hour_seq, mu_s, lwd = 2, col = stage_cols_6[s], lty = sex_lty_6[s])
  pts <- E_stage_sex %>% filter(stage_sex_id == s); points(pts$hour, pts$mean_gaps, pch = ifelse(sex_lty_6[s] == 1, 16, 1), col = stage_cols_6[s], cex = 0.7)
}
legend("topright", legend = combos, col = stage_cols_6, lty = sex_lty_6, pch = rep(c(16, 1), 3), lwd = 2, cex = 0.7, bg = "white")

# --- M2 CONTRASTS & STATS ---
mu_list_m2_sex <- lapply(1:6, function(s) {
  pred_s <- link(m2_sex, data = c(list(stage_sex_id = rep(s, 24)), basis_list_sex))
  apply(pred_s, 1, mean)
})
names(mu_list_m2_sex) <- c("early_f", "early_m", "mid_f", "mid_m", "late_f", "late_m")

cat("\n--- M2 SEX x STAGE SUMMARY ---\n")
for (i in seq_along(stage_names)) {
  f_n <- paste0(stage_names[i], "_f"); m_n <- paste0(stage_names[i], "_m")
  diff_s <- mu_list_m2_sex[[f_n]] - mu_list_m2_sex[[m_n]]
  cat(toupper(stage_names[i]), "| F-M Diff:", round(mean(diff_s), 3), "P(F > M) =", round(mean(diff_s > 0), 2), "\n")
}

# --- PEAK HOURS (M2) ---
cat("\n--- PEAK HOURS (M2) ---\n")
for (s in 1:6) {
  pred_s <- link(m2_sex, data = c(list(stage_sex_id = rep(s, 24)), basis_list_sex))
  peak_h <- apply(pred_s, 1, which.max) - 1
  cat(str_pad(combos[s], 8, "right"), ": mean peak =", round(mean(peak_h), 1), "\n")
}

# ------------------------------------------------------------------------------
# 8. M2 VISUALIZATIONS: CONTRASTS & DENSITIES
# ------------------------------------------------------------------------------

# 1. Sex Contrast Plots by Stage
par(mfrow = c(1, 3))
for (i in seq_along(stage_names)) {
  f_n <- paste0(stage_names[i], "_f"); m_n <- paste0(stage_names[i], "_m")
  diff_s <- mu_list_m2_sex[[f_n]] - mu_list_m2_sex[[m_n]]
  dens(diff_s, col = stage_colors[i], lwd = 3, xlab = "Difference (F - M)", main = paste("Contrast:", toupper(stage_names[i])))
  abline(v = 0, lty = 2); mtext(paste0("P(F > M) = ", round(mean(diff_s > 0), 2)), cex = 0.8)
}
par(mfrow = c(1, 1))

# 2. Stage Comparisons within Sexes
cat("\n--- STAGE COMPARISONS BY SEX ---\n")
stage_pairs <- list(c("early", "mid"), c("mid", "late"), c("early", "late"))
sexes <- c("f", "m")
for (sex in sexes) {
  cat("\n---", ifelse(sex == "f", "FEMALES", "MALES"), "---\n")
  for (pair in stage_pairs) {
    g1 <- paste0(pair[1], "_", sex); g2 <- paste0(pair[2], "_", sex)
    diff_post <- mu_list_m2_sex[[g1]] - mu_list_m2_sex[[g2]]
    cat(str_pad(paste(toupper(pair[1]), "-", toupper(pair[2])), 15, "right"), "| Mean diff =", round(mean(diff_post), 3), "| P(>0) =", round(mean(diff_post > 0), 3), "\n")
  }
}

# Density Plots: Stage comparisons
par(mfrow = c(2, 3))
for (sex in sexes) {
  for (pair in stage_pairs) {
    g1 <- paste0(pair[1], "_", sex); g2 <- paste0(pair[2], "_", sex)
    diff_post <- mu_list_m2_sex[[g1]] - mu_list_m2_sex[[g2]]
    dens(diff_post, lwd = 3, col = ifelse(sex == "f", "darkred", "darkblue"), main = paste(ifelse(sex == "f", "Females", "Males"), "\n", toupper(pair[1]), "-", toupper(pair[2])))
    abline(v = 0, lty = 2); mtext(paste0("P(>0) = ", round(mean(diff_post > 0), 2)), cex = 0.8)
  }
}
par(mfrow = c(1, 1))

# Final Distribution Overview
par(mfrow = c(1, 2))
for (sex in c("f", "m")) {
  plot(NULL, xlim = range(unlist(mu_list_m2_sex)), ylim = c(0, 100), xlab = "Mean exit rate", ylab = "Density", main = ifelse(sex == "f", "Females", "Males"))
  for (i in seq_along(stage_names)) {
    g <- paste0(stage_names[i], "_", sex)
    dens(mu_list_m2_sex[[g]], col = stage_colors[i], lwd = 3, add = TRUE)
  }
  legend("topright", legend = toupper(stage_names), col = stage_colors, lwd = 3, bty = "n")
}
par(mfrow = c(1, 1))

