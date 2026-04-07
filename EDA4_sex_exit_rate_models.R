# PURPOSE: model n_gaps within 24 hours for both sexes
rm(list = ls())
library(rethinking)
library(tidyverse)
library(lubridate)
library(splines)

S <- readRDS("EDA_sex_neglect_wide.RDS")
head(S)
G <- readRDS("EDA_sex_gaps_hourly.RDS")

hour_seq     <- 0:23
sex_col      <- c("#D4537E", "#185FA5")
sex_labels   <- c("Female", "Male")
stage_names  <- c("early", "mid", "late")
stage_colors <- c("#3B6D11", "#185FA5", "#993556")

# ── BUILD HOURLY TABLE ────────────────────────────────────────────────────────

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
      day_prior_hatch >=  0 & day_prior_hatch <= 10 ~ "late",
      TRUE ~ NA_character_
    ),
    stage_id = as.integer(factor(incubation_stage, levels = stage_names)),
    sex_id   = as.integer(factor(sx, levels = c("f", "m"))),
    stage_sex_id = as.integer(factor(
      paste(incubation_stage, sx, sep = "_"),
      levels = c("early_f", "early_m", "mid_f", "mid_m", "late_f", "late_m")
    ))
  ) %>%
  select(season, session, nest, sx, day_prior_hatch, incubation_stage, stage_id, sex_id, stage_sex_id)

# Łączenie i wypełnianie zerami
hour_grid <- meta %>% crossing(hour = 0:23)

E_sex <- hour_grid %>%
  left_join(sex_hourly_raw, by = c("season", "session", "nest", "sx", "hour")) %>%
  mutate(
    n_gaps_per_hour = replace_na(n_gaps_per_hour, 0)
  )

# --- BRAKUJĄCY ELEMENT 1: Agregacja do wykresów ---
E_by_sex <- E_sex %>%
  group_by(sex_id, hour) %>%
  summarise(mean_gaps = mean(n_gaps_per_hour), .groups = "drop")

E_stage_sex <- E_sex %>%
  group_by(stage_sex_id, hour) %>%
  summarise(mean_gaps = mean(n_gaps_per_hour), .groups = "drop")

# ── SPLINES ───────────────────────────────────────────────────────────────────

E_fit <- E_sex %>% filter(!is.na(sex_id))
B <- bs(E_fit$hour, df = 6, degree = 3, intercept = TRUE)

B_pred <- bs(hour_seq, knots = attr(B, "knots"), 
             Boundary.knots = attr(B, "Boundary.knots"), 
             degree = 3, intercept = TRUE)

n_basis <- ncol(B)
B_scaled      <- scale(B)
B_pred_scaled <- scale(B_pred, 
                       center = attr(B_scaled, "scaled:center"), 
                       scale  = attr(B_scaled, "scaled:scale"))

basis_cols <- setNames(as.data.frame(B_scaled), paste0("B", seq_len(n_basis)))
basis_list <- as.list(setNames(as.data.frame(B_pred_scaled), paste0("B", seq_len(n_basis))))

# ── M1: SEX ONLY ──────────────────────────────────────────────────────────────

# Naprawa d_m1: musi zawierać kolumny basis_cols (dane obserwowane), a nie pred
d_m1 <- c(
  list(n_gaps_per_hour = E_fit$n_gaps_per_hour, sex_id = E_fit$sex_id),
  as.list(basis_cols)
)

m1 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[sex_id] + 
      w1[sex_id]*B1 + w2[sex_id]*B2 + w3[sex_id]*B3 + 
      w4[sex_id]*B4 + w5[sex_id]*B5 + w6[sex_id]*B6,
    a[sex_id] ~ dnorm(0, 1),
    w1[sex_id] ~ dnorm(0, 0.4),
    w2[sex_id] ~ dnorm(0, 0.4),
    w3[sex_id] ~ dnorm(0, 0.4),
    w4[sex_id] ~ dnorm(0, 0.4),
    w5[sex_id] ~ dnorm(0, 0.4),
    w6[sex_id] ~ dnorm(0, 0.4)
  ),
  data = d_m1,
  start = list(a=rep(0,2), w1=rep(0,2), w2=rep(0,2), w3=rep(0,2), w4=rep(0,2), w5=rep(0,2), w6=rep(0,2))
)

# Wykres M1
plot(NULL, xlim = c(0, 23), ylim = c(0, max(E_by_sex$mean_gaps) * 1.5),
     xlab = "Hour", ylab = "Expected exits", main = "Daily exit pattern by sex")
axis(1, at = 0:23)

for (s in 1:2) {
  pred_s <- link(m1, data = c(list(sex_id = rep(s, 24)), basis_list))
  mu_s <- apply(pred_s, 2, mean)
  PI_s <- apply(pred_s, 2, PI, prob = 0.89)
  lines(hour_seq, mu_s, lwd = 2, col = sex_col[s])
  shade(PI_s, hour_seq, col = col.alpha(sex_col[s], 0.15))
  pts <- E_by_sex %>% filter(sex_id == s)
  points(pts$hour, pts$mean_gaps, pch = 16, col = sex_col[s])
}

mu_f_m1 <- link(m1, data = c(list(sex_id = rep(1, 24)), basis_list))
mu_m_m1 <- link(m1, data = c(list(sex_id = rep(2, 24)), basis_list))

# Średnia liczba wyjść na godzinę (całodobowo)
avg_rate_f_m1 <- apply(mu_f_m1, 1, mean)
avg_rate_m_m1 <- apply(mu_m_m1, 1, mean)
diff_fm_m1    <- avg_rate_f_m1 - avg_rate_m_m1

cat("\n=== MODEL 1: SEX SUMMARY ===\n")
cat("Female mean rate:", round(mean(avg_rate_f_m1), 3), "PI:", round(PI(avg_rate_f_m1), 3), "\n")
cat("Male mean rate:  ", round(mean(avg_rate_m_m1), 3), "PI:", round(PI(avg_rate_m_m1), 3), "\n")
cat("Difference (F-M):", round(mean(diff_fm_m1), 3), "PI:", round(PI(diff_fm_m1), 3), 
    "| P(F > M) =", round(mean(diff_fm_m1 > 0), 2), "\n")

# --- GODZINY SZCZYTU M1 ---
for (s in 1:2) {
  pred_s <- link(m1, data = c(list(sex_id = rep(s, 24)), basis_list))
  peak_hour <- apply(pred_s, 1, which.max) - 1
  cat(sex_labels[s], "Peak Hour: mean =", round(mean(peak_hour), 1), 
      "PI:", round(PI(peak_hour), 1), "\n")
}

# Przygotowanie okna graficznego (1 wiersz, 2 kolumny)
par(mfrow = c(1, 2))

# A. Nakładające się rozkłady średniej aktywności
xr <- range(c(avg_rate_f_m1, avg_rate_m_m1))
plot(NULL, xlim = xr, ylim = c(0, 150), 
     xlab = "Mean expected exits / hour", ylab = "Density",
     main = "Posterior Mean Rate (M1)")
dens(avg_rate_f_m1, col = sex_col[1], lwd = 3, add = TRUE)
dens(avg_rate_m_m1, col = sex_col[2], lwd = 3, add = TRUE)
legend("topright", legend = sex_labels, col = sex_col, lwd = 3, bty = "n")

# B. Rozkład różnicy (F - M)
dens(diff_fm_m1, col = "purple", lwd = 3,
     xlab = "Difference (exits/hour)", 
     main = "Contrast: Female - Male")
abline(v = 0, lty = 2, lwd = 1.5) # Linia zero - brak różnicy









# ── M2: SEX X STAGE ──────────────────────────────────────────────────────────

# Przygotowanie danych do M2
E_fit2 <- E_sex %>% filter(!is.na(stage_sex_id))
# Ponieważ używamy splajnów na tych samych danych 'hour', możemy użyć basis_cols
d_m2 <- c(
  list(n_gaps_per_hour = E_fit2$n_gaps_per_hour, stage_sex_id = E_fit2$stage_sex_id),
  as.list(basis_cols)
)

m2 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[stage_sex_id] + 
      w1[stage_sex_id]*B1 + w2[stage_sex_id]*B2 + w3[stage_sex_id]*B3 +
      w4[stage_sex_id]*B4 + w5[stage_sex_id]*B5 + w6[stage_sex_id]*B6,
    a[stage_sex_id] ~ dnorm(0, 1),
    w1[stage_sex_id] ~ dnorm(0, 0.4),
    w2[stage_sex_id] ~ dnorm(0, 0.4),
    w3[stage_sex_id] ~ dnorm(0, 0.4),
    w4[stage_sex_id] ~ dnorm(0, 0.4),
    w5[stage_sex_id] ~ dnorm(0, 0.4),
    w6[stage_sex_id] ~ dnorm(0, 0.4)
  ),
  data = d_m2,
  start = list(a=rep(0,6), w1=rep(0,6), w2=rep(0,6), w3=rep(0,6), w4=rep(0,6), w5=rep(0,6), w6=rep(0,6))
)

# --- BRAKUJĄCY ELEMENT 2: Definicje do wykresu M2 ---
combos <- c("Early F", "Early M", "Mid F", "Mid M", "Late F", "Late M")
stage_cols_6 <- rep(stage_colors, each = 2)
sex_lty_6    <- rep(c(1, 2), 3)
# Przygotowujemy kolory dla cieniowania (półprzezroczyste)
# Używamy stage_cols_6, aby kolor cienia odpowiadał etapowi inkubacji
fill_cols <- col.alpha(stage_cols_6, 0.15)

plot(NULL, 
     xlim = c(0, 23), 
     ylim = c(0, 0.5),
     xlab = "Hour", 
     ylab = "Expected exits", 
     xaxt = "n",
     main = "Daily exit pattern — Sex x Stage\n(points = empirical means, lines = spline fit)")

# Dodajemy oś X z pełnym zakresem godzin
axis(1, at = 0:23, labels = 0:23)

for (s in 1:6) {
  # 1. Obliczamy predykcje z modelu
  pred_s <- link(m2, data = c(list(stage_sex_id = rep(s, 24)), basis_list))
  mu_s   <- apply(pred_s, 2, mean)
  PI_s   <- apply(pred_s, 2, PI, prob = 0.89)
  
  # 2. Rysujemy cieniowanie przedziału ufności
  shade(PI_s, hour_seq, col = fill_cols[s])
  
  # 3. Rysujemy linię średniej z modelu
  # lty = 1 dla samic (nieparzyste s), lty = 2 dla samców (parzyste s)
  lines(hour_seq, mu_s, lwd = 2, col = stage_cols_6[s], lty = sex_lty_6[s])
  
  # 4. Dodajemy punkty empiryczne (średnie z danych)
  # Filtrujemy dane dla konkretnej grupy stage_sex_id
  pts <- E_stage_sex %>% filter(stage_sex_id == s)
  points(pts$hour, pts$mean_gaps, 
         pch = ifelse(sex_lty_6[s] == 1, 16, 1), # pełne kółko dla F, puste dla M
         col = stage_cols_6[s], 
         cex = 0.7)
}

# Legenda
legend("topright", 
       legend = combos, 
       col = stage_cols_6, 
       lty = sex_lty_6, 
       pch = rep(c(16, 1), 3), 
       lwd = 2, 
       cex = 0.7, 
       bg = "white")

# --- PODSUMOWANIE M2: ŚREDNIE DLA 6 GRUP ---
mu_list_m2 <- lapply(1:6, function(s) {
  pred_s <- link(m2, data = c(list(stage_sex_id = rep(s, 24)), basis_list))
  apply(pred_s, 1, mean)
})
names(mu_list_m2) <- c("early_f", "early_m", "mid_f", "mid_m", "late_f", "late_m")

cat("\n=== MODEL 2: SEX x STAGE SUMMARY ===\n")

# Kontrasty: Samica - Samiec w każdym etapie
for (i in seq_along(stage_names)) {
  f_name <- paste0(stage_names[i], "_f")
  m_name <- paste0(stage_names[i], "_m")
  
  diff_stage <- mu_list_m2[[f_name]] - mu_list_m2[[m_name]]
  
  cat(str_pad(toupper(stage_names[i]), 6, "right"), "| F-M Difference:", 
      round(mean(diff_stage), 3), "PI:", round(PI(diff_stage), 3),
      "| P(F > M) =", round(mean(diff_stage > 0), 2), "\n")
}

# --- GODZINY SZCZYTU M2 ---
cat("\n--- PEAK HOURS (M2) ---\n")
for (s in 1:6) {
  pred_s <- link(m2, data = c(list(stage_sex_id = rep(s, 24)), basis_list))
  peak_hour <- apply(pred_s, 1, which.max) - 1
  cat(str_pad(combos[s], 8, "right"), ": mean peak =", round(mean(peak_hour), 1), 
      "PI:", round(PI(peak_hour), 1), "\n")
}


# Przygotowanie okna graficznego (1 wiersz, 3 kolumny)
par(mfrow = c(1, 3))

for (i in seq_along(stage_names)) {
  f_name <- paste0(stage_names[i], "_f")
  m_name <- paste0(stage_names[i], "_m")
  
  # Obliczamy różnicę dla danego etapu
  diff_stage <- mu_list_m2[[f_name]] - mu_list_m2[[m_name]]
  
  # Wykres gęstości różnicy
  dens(diff_stage, col = stage_colors[i], lwd = 3,
       xlab = "Difference (F - M)",
       main = paste("Contrast:", toupper(stage_names[i])))
  
  # Dodajemy linię zero i informację o prawdopodobieństwie kierunku
  abline(v = 0, lty = 2)
  
  # Opcjonalnie: dodaj tekst z % pewności, że F > M
  p_f_gt_m <- mean(diff_stage > 0)
  mtext(paste0("P(F > M) = ", round(p_f_gt_m, 2)), cex = 0.8)
}

# Resetowanie ustawień graficznych
par(mfrow = c(1, 1))

#differences between session for males and females

# Definicje kontrastów etapów
stage_pairs <- list(
  c("early", "mid"),
  c("mid", "late"),
  c("early", "late")
)

sexes <- c("f", "m")

for (sex in sexes) {
  
  cat("\n---", ifelse(sex == "f", "FEMALES", "MALES"), "---\n")
  
  for (pair in stage_pairs) {
    
    g1 <- paste0(pair[1], "_", sex)
    g2 <- paste0(pair[2], "_", sex)
    
    diff_post <- mu_list_m2[[g1]] - mu_list_m2[[g2]]
    
    cat(
      str_pad(paste(toupper(pair[1]), "-", toupper(pair[2])), 15, "right"),
      "| Mean diff =", round(mean(diff_post), 3),
      "| PI =", round(PI(diff_post), 3),
      "| P(>0) =", round(mean(diff_post > 0), 3),
      "\n"
    )
  }
}

par(mfrow = c(2, 3))

for (sex in sexes) {
  for (pair in stage_pairs) {
    
    g1 <- paste0(pair[1], "_", sex)
    g2 <- paste0(pair[2], "_", sex)
    
    diff_post <- mu_list_m2[[g1]] - mu_list_m2[[g2]]
    
    dens(
      diff_post,
      lwd = 3,
      col = ifelse(sex == "f", "darkred", "darkblue"),
      main = paste(
        ifelse(sex == "f", "Females", "Males"),
        "\n",
        toupper(pair[1]), "-", toupper(pair[2])
      ),
      xlab = "Difference in mean exit rate"
    )
    
    abline(v = 0, lty = 2)
    
    mtext(
      paste0("P(>0) = ", round(mean(diff_post > 0), 2)),
      cex = 0.8
    )
  }
}

par(mfrow = c(1,1))

par(mfrow = c(1, 2))

stage_pairs <- c("early", "mid", "late")
cols <- stage_colors

for (sex in c("f", "m")) {
  
  plot(NULL,
       xlim = range(unlist(mu_list_m2)),
       ylim = c(0, 100),
       xlab = "Mean exit rate",
       ylab = "Density",
       main = ifelse(sex == "f", "Females", "Males"))
  
  for (i in seq_along(stage_pairs)) {
    
    g <- paste0(stage_pairs[i], "_", sex)
    
    dens(
      mu_list_m2[[g]],
      col = cols[i],
      lwd = 3,
      add = TRUE
    )
  }
  
  legend("topright",
         legend = toupper(stage_pairs),
         col = cols,
         lwd = 3,
         bty = "n")
}

par(mfrow = c(1,1))
# =========================================================
# LINEAR MODEL: n_gaps changes toward hatch
# =========================================================
# =========================================================
# 1. DATA
# =========================================================

E_fit_lin <- S %>%
  filter(
    !is.na(day_prior_hatch),
    !is.na(sx),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    sex_id = ifelse(sx == "f", 1, 2),
    
    # centered predictor helps interpretation + sampling
    dph_c = day_prior_hatch - mean(day_prior_hatch)
  )

d_lin <- list(
  n_gaps = E_fit_lin$n_gaps,
  sex_id = E_fit_lin$sex_id,
  dph_c = E_fit_lin$dph_c
)

# =========================================================
# 2. PRIOR PREDICTIVE SIMULATION
# =========================================================

set.seed(123)

day_seq <- seq(2, 28, length.out = 100)
dph_seq_c <- day_seq - mean(E_fit_lin$day_prior_hatch)

n_sim <- 100
prior_f <- matrix(NA, nrow = n_sim, ncol = length(day_seq))
prior_m <- matrix(NA, nrow = n_sim, ncol = length(day_seq))

for(i in 1:n_sim){
  
  a <- rnorm(2, 0, 0.3)
  b <- rnorm(2, 0, 0.05)
  
  prior_f[i,] <- exp(a[1] + b[1] * dph_seq_c)
  prior_m[i,] <- exp(a[2] + b[2] * dph_seq_c)
}

par(mfrow = c(1,2))

matplot(day_seq, t(prior_f),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected n_gaps",
        main = "Prior predictive: Female")

matplot(day_seq, t(prior_m),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected n_gaps",
        main = "Prior predictive: Male")

par(mfrow = c(1,1))

# =========================================================
# 3. MODEL
# =========================================================

m_lin <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    
    log(lambda) <- a[sex_id] + b[sex_id] * dph_c,
    
    a[sex_id] ~ dnorm(0, 0.3),
    b[sex_id] ~ dnorm(0, 0.05)
  ),
  data = d_lin,
  start = list(
    a = c(0, 0),
    b = c(0, 0)
  )
)

precis(m_lin, depth = 2)

# =========================================================
# 4. POSTERIOR PREDICTIONS
# =========================================================

post_f <- link(
  m_lin,
  data = list(
    sex_id = rep(1, length(day_seq)),
    dph_c = dph_seq_c
  )
)

post_m <- link(
  m_lin,
  data = list(
    sex_id = rep(2, length(day_seq)),
    dph_c = dph_seq_c
  )
)

mu_f <- apply(post_f, 2, mean)
mu_m <- apply(post_m, 2, mean)

PI_f <- apply(post_f, 2, PI, prob = 0.89)
PI_m <- apply(post_m, 2, PI, prob = 0.89)

# =========================================================
# 5. POSTERIOR PLOT
# =========================================================

# =========================================================
# POSTERIOR + RAW DATA
# =========================================================

plot(NULL,
     xlim = c(28, 2),
     ylim = c(0, max(E_fit_lin$n_gaps, mu_f, mu_m) * 1.2),
     xlab = "Days prior hatch",
     ylab = "n_gaps",
     main = "Linear change in n_gaps across incubation")

# -----------------------------
# raw data (jittered)
# -----------------------------
points(
  jitter(E_fit_lin$day_prior_hatch, amount = 0.25),
  jitter(E_fit_lin$n_gaps, amount = 0.1),
  pch = ifelse(E_fit_lin$sex_id == 1, 16, 1),
  col = ifelse(E_fit_lin$sex_id == 1, "darkred", "darkblue"),
  cex = 0.7
)

# -----------------------------
# posterior intervals
# -----------------------------
shade(PI_f, day_seq, col = col.alpha("darkred", 0.2))
shade(PI_m, day_seq, col = col.alpha("darkblue", 0.2))

# posterior means
lines(day_seq, mu_f, lwd = 3, col = "darkred")
lines(day_seq, mu_m, lwd = 3, col = "darkblue")

legend("topleft",
       legend = c("Female raw", "Male raw", "Female fit", "Male fit"),
       pch = c(16, 1, NA, NA),
       lty = c(NA, NA, 1, 1),
       col = c("darkred", "darkblue", "darkred", "darkblue"),
       lwd = c(NA, NA, 3, 3),
       bty = "n")

# =========================================================
# 6. INTERPRETATION OF SLOPES
# =========================================================

post <- extract.samples(m_lin)

cat("\n=== SLOPES ===\n")

cat("Female slope:",
    round(mean(post$b[,1]), 4),
    "PI:", round(PI(post$b[,1]), 4),
    "P(<0):", round(mean(post$b[,1] < 0), 3), "\n")

cat("Male slope:",
    round(mean(post$b[,2]), 4),
    "PI:", round(PI(post$b[,2]), 4),
    "P(<0):", round(mean(post$b[,2] < 0), 3), "\n")

# =========================================================
# 7. SLOPE DIFFERENCE
# =========================================================

slope_diff <- post$b[,1] - post$b[,2]

cat("\n=== SLOPE DIFFERENCE (F - M) ===\n")
cat(
  round(mean(slope_diff), 4),
  "PI:", round(PI(slope_diff), 4),
  "P(F > M):", round(mean(slope_diff > 0), 3), "\n"
)

# =========================================================
# SEX DIFFERENCES AT INCUBATION 1
# =========================================================
# =========================================================
# SEX DIFFERENCE ONLY IN INCUBATION1
# adjusted for day_prior_hatch
# =========================================================
# =========================================================
# 1. FILTER ONLY INCUBATION1
# =========================================================

E_inc1 <- S %>%
  filter(
    session == "incubation1",
    !is.na(day_prior_hatch),
    !is.na(sx),
    day_prior_hatch >= 2,
    day_prior_hatch <= 28
  ) %>%
  mutate(
    sex_id = ifelse(sx == "f", 1, 2),
    dph_c = day_prior_hatch - mean(day_prior_hatch)
  )

d_inc1 <- list(
  n_gaps = E_inc1$n_gaps,
  sex_id = E_inc1$sex_id,
  dph_c = E_inc1$dph_c
)

# =========================================================
# 2. MODEL
# same slope, different intercepts
# =========================================================

m_inc1 <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    
    log(lambda) <- a[sex_id] + b * dph_c,
    
    a[sex_id] ~ dnorm(0, 1),
    b ~ dnorm(0, 0.05)
  ),
  data = d_inc1,
  start = list(
    a = c(0, 0),
    b = 0
  )
)

precis(m_inc1, depth = 2)

post_inc1 <- extract.samples(m_inc1)

sex_diff_inc1 <- post_inc1$a[,1] - post_inc1$a[,2]

cat("\n=== INCUBATION1 SEX DIFFERENCE (F - M) ===\n")
cat(
  round(mean(sex_diff_inc1), 3),
  "PI:", round(PI(sex_diff_inc1), 3),
  "P(F > M):", round(mean(sex_diff_inc1 > 0), 3), "\n"
)

#plot 
day_seq <- seq(
  min(E_inc1$day_prior_hatch),
  max(E_inc1$day_prior_hatch),
  length.out = 100
)

dph_seq_c <- day_seq - mean(E_inc1$day_prior_hatch)

post_f <- link(
  m_inc1,
  data = list(
    sex_id = rep(1, length(day_seq)),
    dph_c = dph_seq_c
  )
)

post_m <- link(
  m_inc1,
  data = list(
    sex_id = rep(2, length(day_seq)),
    dph_c = dph_seq_c
  )
)

mu_f <- apply(post_f, 2, mean)
mu_m <- apply(post_m, 2, mean)

PI_f <- apply(post_f, 2, PI)
PI_m <- apply(post_m, 2, PI)

plot(NULL,
     xlim = c(max(day_seq), min(day_seq)),
     ylim = c(0, max(E_inc1$n_gaps, mu_f, mu_m) * 1.2),
     xlab = "Days prior hatch",
     ylab = "n_gaps",
     main = "Incubation1 only")

points(
  jitter(E_inc1$day_prior_hatch, 0.2),
  jitter(E_inc1$n_gaps, 0.1),
  pch = ifelse(E_inc1$sex_id == 1, 16, 1),
  col = ifelse(E_inc1$sex_id == 1, "darkred", "darkblue"),
  cex = 0.7
)

shade(PI_f, day_seq, col = col.alpha("darkred", 0.2))
shade(PI_m, day_seq, col = col.alpha("darkblue", 0.2))

lines(day_seq, mu_f, col = "darkred", lwd = 3)
lines(day_seq, mu_m, col = "darkblue", lwd = 3)

legend("topleft",
       legend = c("Female", "Male"),
       col = c("darkred", "darkblue"),
       pch = c(16, 1),
       lwd = 3,
       bty = "n")

# =========================================================
# 1. Posterior distributions of a[sex_id] (intercepts)
# =========================================================

par(mfrow = c(1,2))

# Female posterior
dens(post_inc1$a[,1], col = "darkred", lwd = 3,
     xlab = "Intercept (log(lambda))",
     main = "Incubation1: Female vs Male")
# Male posterior
dens(post_inc1$a[,2], col = "darkblue", lwd = 3, add = TRUE)

abline(v = 0, lty = 2) # zero reference
legend("topright",
       legend = c("Female", "Male"),
       col = c("darkred", "darkblue"),
       lwd = 3,
       bty = "n")

# =========================================================
# 2. Posterior difference: Female - Male
# =========================================================

dens(sex_diff_inc1, col = "purple", lwd = 3,
     xlab = "F - M (log(lambda))",
     main = "Posterior difference (F - M)")
abline(v = 0, lty = 2)

# Add text for probability
p_f_gt_m <- mean(sex_diff_inc1 > 0)
mtext(paste0("P(F > M) = ", round(p_f_gt_m, 3)), cex = 0.8)

par(mfrow = c(1,1))

# ── M3: SEX X SEASON ─────────────────────────────────────────────

# =========================================================
# 1. DATA PREPARATION
# =========================================================
E_sex <- E_sex %>%
  mutate(
    season = as.character(season),
    sex_season = paste(sx, season, sep = "_"),
    sex_season_id = as.integer(factor(sex_season))
  )

E_fit3 <- E_sex %>%
  filter(!is.na(season))

season_levels <- sort(unique(E_fit3$season))
n_seasons <- length(season_levels)

group_lookup <- E_fit3 %>%
  distinct(sex_season_id, sx, season) %>%
  arrange(sx, season)

n_groups <- length(unique(E_fit3$sex_season_id))

# empirical means for plotting
E_season_sex <- E_fit3 %>%
  group_by(season, sx, sex_season_id, hour) %>%
  summarise(mean_gaps = mean(n_gaps_per_hour), .groups = "drop")

# =========================================================
# 2. MODEL
# =========================================================
d_m3 <- c(
  list(
    n_gaps_per_hour = E_fit3$n_gaps_per_hour,
    sex_season_id   = E_fit3$sex_season_id
  ),
  as.list(basis_cols)
)

m3_season <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    
    log(lambda) <- a[sex_season_id] +
      w1[sex_season_id]*B1 +
      w2[sex_season_id]*B2 +
      w3[sex_season_id]*B3 +
      w4[sex_season_id]*B4 +
      w5[sex_season_id]*B5 +
      w6[sex_season_id]*B6,
    
    a[sex_season_id]  ~ dnorm(0, 1),
    w1[sex_season_id] ~ dnorm(0, 0.4),
    w2[sex_season_id] ~ dnorm(0, 0.4),
    w3[sex_season_id] ~ dnorm(0, 0.4),
    w4[sex_season_id] ~ dnorm(0, 0.4),
    w5[sex_season_id] ~ dnorm(0, 0.4),
    w6[sex_season_id] ~ dnorm(0, 0.4)
  ),
  data = d_m3,
  start = list(
    a  = rep(0, n_groups),
    w1 = rep(0, n_groups),
    w2 = rep(0, n_groups),
    w3 = rep(0, n_groups),
    w4 = rep(0, n_groups),
    w5 = rep(0, n_groups),
    w6 = rep(0, n_groups)
  )
)

# =========================================================
# 3. PRIOR PREDICTIVE CHECK
# =========================================================
set.seed(123)

plot(NULL,
     xlim = c(0, 23),
     ylim = c(0, 15),
     xlab = "Hour",
     ylab = "Expected exits",
     main = "Prior predictive simulation — Season model")

for(i in 1:20){
  
  a  <- rnorm(1, 0, 1)
  ws <- rnorm(6, 0, 0.4)
  
  log_lambda <- a +
    ws[1]*basis_list$B1 +
    ws[2]*basis_list$B2 +
    ws[3]*basis_list$B3 +
    ws[4]*basis_list$B4 +
    ws[5]*basis_list$B5 +
    ws[6]*basis_list$B6
  
  lines(hour_seq, exp(log_lambda),
        col = col.alpha("black", 0.3))
}

# =========================================================
# 4. RAW DATA + SPLINE FITS BY SEASON
# =========================================================
sex_cols <- c(f = "red", m = "blue")
sex_lty  <- c(f = 1, m = 2)

ncol_plot <- 2
nrow_plot <- ceiling(n_seasons / ncol_plot)

par(mfrow = c(nrow_plot, ncol_plot),
    mar = c(4, 4, 3, 1))

for(seas in season_levels) {
  
  season_groups <- group_lookup %>%
    filter(season == seas)
  
  pts_season <- E_season_sex %>%
    filter(season == seas)
  
  plot(NULL,
       xlim = c(0, 23),
       ylim = c(0, max(E_season_sex$mean_gaps) * 1.4),
       xlab = "Hour",
       ylab = "Expected exits",
       xaxt = "n",
       main = paste("Season:", seas))
  
  axis(1, at = 0:23)
  
  for(i in 1:nrow(season_groups)) {
    
    id_now <- season_groups$sex_season_id[i]
    sx_now <- season_groups$sx[i]
    
    pred_i <- link(
      m3_season,
      data = c(
        list(sex_season_id = rep(id_now, 24)),
        basis_list
      )
    )
    
    mu_i <- apply(pred_i, 2, mean)
    PI_i <- apply(pred_i, 2, PI, prob = 0.89)
    
    shade(PI_i, hour_seq,
          col = col.alpha(sex_cols[sx_now], 0.12))
    
    lines(hour_seq, mu_i,
          col = sex_cols[sx_now],
          lty = sex_lty[sx_now],
          lwd = 2)
    
    pts <- pts_season %>%
      filter(sex_season_id == id_now)
    
    points(pts$hour, pts$mean_gaps,
           pch = ifelse(sx_now == "f", 16, 1),
           col = sex_cols[sx_now],
           cex = 0.7)
  }
}

# =========================================================
# 5. KEY PLOT: POSTERIOR DISTRIBUTIONS OF ALL SEASONS
#    WITHIN EACH SEX (ONE COMMON SCALE)
# =========================================================
season_cols <- rainbow(n_seasons)

# collect posterior means
mu_by_group <- list()

for(i in 1:nrow(group_lookup)) {
  
  id_now <- group_lookup$sex_season_id[i]
  label  <- paste(group_lookup$sx[i], group_lookup$season[i], sep = "_")
  
  pred_i <- link(
    m3_season,
    data = c(list(sex_season_id = rep(id_now, 24)), basis_list)
  )
  
  mu_by_group[[label]] <- apply(pred_i, 1, mean)
}

x_range <- range(unlist(mu_by_group))

par(mfrow = c(1, 2),
    mar = c(4, 4, 3, 1))

for(sex_now in c("f", "m")) {
  
  plot(NULL,
       xlim = x_range,
       ylim = c(0,100),
       xlab = "Mean exit rate",
       ylab = "Posterior density",
       main = ifelse(sex_now == "f", "Females", "Males"))
  
  for(i in seq_along(season_levels)) {
    
    lab_now <- paste(sex_now, season_levels[i], sep = "_")
    
    dens(mu_by_group[[lab_now]],
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
# 6. NUMERICAL POSTERIOR SUMMARY:
#    SEASON DIFFERENCES WITHIN SEX
# =========================================================
cat("\n=== MODEL 3: SEASON DIFFERENCES WITHIN SEX ===\n")

season_pairs <- combn(season_levels, 2, simplify = FALSE)

for(sex_now in c("f", "m")) {
  
  cat("\n---", ifelse(sex_now == "f", "FEMALES", "MALES"), "---\n")
  
  for(pair in season_pairs) {
    
    g1 <- paste(sex_now, pair[1], sep = "_")
    g2 <- paste(sex_now, pair[2], sep = "_")
    
    diff_post <- mu_by_group[[g1]] - mu_by_group[[g2]]
    
    cat(
      pair[1], "-", pair[2],
      "| Mean diff =", round(mean(diff_post), 3),
      "| PI =", round(PI(diff_post), 3),
      "| P(>0) =", round(mean(diff_post > 0), 3),
      "\n"
    )
  }
}


##################### repeatability
library(rptR)
str(S)
S$ringno <- as.factor(S$ringno)
S$log_n_gaps <- log(S$n_gaps + 0.1)

model_gauss <- rpt(n_gaps ~ session + season + (1 | ringno), 
                   grname = "ringno", 
                   data = S, 
                   datatype = "Poisson", 
                   nboot = 1000, 
                   npermut = 1000)

print(model_gauss)

S_male <- S %>%
  filter(sx == "m")
S_female <- S %>%
  filter(sx == "f")

rep_male <- rpt(n_gaps ~ session + season + (1 | ringno), 
                   grname = "ringno", 
                   data = S_male, 
                   datatype = "Poisson", 
                   nboot = 1000, 
                   npermut = 1000)

rep_female <- rpt(n_gaps ~ session + season + (1 | ringno), 
                grname = "ringno", 
                data = S_female, 
                datatype = "Poisson", 
                nboot = 1000, 
                npermut = 1000)

print(rep_male)
print(rep_female)


