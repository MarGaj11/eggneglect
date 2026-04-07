# PURPOSE: model sum_gaps for both sexes
# =========================================================
# TOTAL NEGLECT ACROSS INCUBATION
# Response: sum_neglect_sec
# =========================================================

rm(list = ls())
library(rethinking)
library(tidyverse)
G <- readRDS("EDA_sex_gaps_hourly.RDS")
S <- readRDS("EDA_sex_neglect_wide.RDS")
dens(S$sum_neglect_sec)
# =========================================================
# 1. DATA
# =========================================================

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
# 2. PRIOR PREDICTIVE
# =========================================================

set.seed(123)

day_seq <- seq(2, 28, length.out = 100)
dph_seq_c <- day_seq - mean(D_neg$day_prior_hatch)

n_sim <- 100
prior_f <- matrix(NA, nrow = n_sim, ncol = length(day_seq))
prior_m <- matrix(NA, nrow = n_sim, ncol = length(day_seq))

for(i in 1:n_sim){
  
  a <- rnorm(2, 6, 0.5)
  b <- rnorm(2, 0, 0.08)
  
  prior_f[i,] <- exp(a[1] + b[1] * dph_seq_c)
  prior_m[i,] <- exp(a[2] + b[2] * dph_seq_c)
}

par(mfrow = c(1,2))

matplot(day_seq, t(prior_f),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected neglect (sec)",
        main = "Prior predictive: Female")

matplot(day_seq, t(prior_m),
        type = "l", lty = 1,
        xlim = c(28, 2),
        xlab = "Days prior hatch",
        ylab = "Expected neglect (sec)",
        main = "Prior predictive: Male")

par(mfrow = c(1,1))

# =========================================================
# 3. MODEL
# =========================================================

m_neg <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    
    mu <- a[sex_id] + b[sex_id] * dph_c,
    
    a[sex_id] ~ dnorm(6, 0.5),
    b[sex_id] ~ dnorm(0, 0.08),
    sigma ~ dexp(1)
  ),
  data = d_neg,
  start = list(
    a = c(6, 6),
    b = c(0, 0),
    sigma = 1
  )
)

precis(m_neg, depth = 2)

# =========================================================
# 4. POSTERIOR PREDICTIONS (LOG Y AXIS)
# =========================================================

post_f <- link(
  m_neg,
  data = list(
    sex_id = rep(1, length(day_seq)),
    dph_c = dph_seq_c
  )
)

post_m <- link(
  m_neg,
  data = list(
    sex_id = rep(2, length(day_seq)),
    dph_c = dph_seq_c
  )
)

mu_f <- exp(apply(post_f, 2, mean))
mu_m <- exp(apply(post_m, 2, mean))

PI_f <- apply(exp(post_f), 2, PI, prob = 0.89)
PI_m <- apply(exp(post_m), 2, PI, prob = 0.89)

# small constant avoids log(0)
eps <- 1

plot(NULL,
     xlim = c(28, 2),
     ylim = c(1, max(D_neg$sum_neglect_sec, mu_f, mu_m) * 1.2),
     log = "y",
     xlab = "Days prior hatch",
     ylab = "Total neglect (sec, log scale)",
     main = "Total neglect across incubation")

points(
  jitter(D_neg$day_prior_hatch, 0.25),
  pmax(jitter(D_neg$sum_neglect_sec, 50), eps),
  pch = ifelse(D_neg$sex_id == 1, 16, 1),
  col = ifelse(D_neg$sex_id == 1, "darkred", "darkblue"),
  cex = 0.7
)

# posterior intervals
shade(pmax(PI_f, eps), day_seq, col = col.alpha("darkred", 0.2))
shade(pmax(PI_m, eps), day_seq, col = col.alpha("darkblue", 0.2))

# posterior means
lines(day_seq, pmax(mu_f, eps), lwd = 3, col = "darkred")
lines(day_seq, pmax(mu_m, eps), lwd = 3, col = "darkblue")

legend("topleft",
       legend = c("Female", "Male"),
       col = c("darkred", "darkblue"),
       pch = c(16, 1),
       lwd = 3,
       bty = "n")

#slope interpretations 
post <- extract.samples(m_neg)

cat("\n=== TOTAL NEGLECT SLOPES ===\n")

cat("Female slope:",
    round(mean(post$b[,1]), 4),
    "PI:", round(PI(post$b[,1]), 4),
    "P(<0):", round(mean(post$b[,1] < 0), 3), "\n")

cat("Male slope:",
    round(mean(post$b[,2]), 4),
    "PI:", round(PI(post$b[,2]), 4),
    "P(<0):", round(mean(post$b[,2] < 0), 3), "\n")

slope_diff <- post$b[,1] - post$b[,2]

cat("\n=== SLOPE DIFFERENCE (F - M) ===\n")
cat(
  round(mean(slope_diff), 4),
  "PI:", round(PI(slope_diff), 4),
  "P(F > M):", round(mean(slope_diff > 0), 3), "\n"
)

# =========================================================
# SEASONAL DIFFERENCES IN TOTAL NEGLECT
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


# =========================================================
# SEX DIFFERENCE IN TOTAL NEGLECT
# adjusted for day_prior_hatch
# =========================================================

D_sex <- S %>%
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
    neglect_log = log(sum_neglect_sec + 1)
  )

d_sex <- list(
  y = D_sex$neglect_log,
  sex_id = D_sex$sex_id,
  dph_c = D_sex$dph_c
)

# =========================================================
# MODEL
# same slope, different intercepts for sexes
# =========================================================

m_sex_neg <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    
    mu <- a[sex_id] + b * dph_c,
    
    a[sex_id] ~ dnorm(6, 0.5),
    b ~ dnorm(0, 0.08),
    sigma ~ dexp(1)
  ),
  data = d_sex,
  start = list(
    a = c(6, 6),
    b = 0,
    sigma = 1
  )
)

precis(m_sex_neg, depth = 2)

post_sex <- extract.samples(m_sex_neg)

sex_diff <- post_sex$a[,2] - post_sex$a[,1]

cat("\n=== SEX DIFFERENCE IN TOTAL NEGLECT (M - F) ===\n")
cat(
  "Mean diff =", round(mean(sex_diff), 3),
  "PI =", round(PI(sex_diff), 3),
  "P(M > F) =", round(mean(sex_diff > 0), 3),
  "\n"
)


# =========================================================
# SEX × SESSION DIFFERENCES IN TOTAL NEGLECT
# =========================================================

session_levels <- c("incubation1", "incubation2", "incubation3")
session_cols   <- c("#3B6D11", "#185FA5", "#993556")

D_session <- S %>%
  filter(
    !is.na(sum_neglect_sec),
    !is.na(session),
    !is.na(sx),
    session %in% session_levels
  ) %>%
  mutate(
    session = factor(session, levels = session_levels),
    sex_session = paste(sx, session, sep = "_"),
    sex_session_id = as.integer(factor(
      sex_session,
      levels = c(
        "f_incubation1", "m_incubation1",
        "f_incubation2", "m_incubation2",
        "f_incubation3", "m_incubation3"
      )
    )),
    neglect_log = log(sum_neglect_sec + 1)
  )


d_session <- list(
  y = D_session$neglect_log,
  sex_session_id = D_session$sex_session_id
)

m_session <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    
    mu <- a[sex_session_id],
    
    a[sex_session_id] ~ dnorm(6, 1),
    sigma ~ dexp(1)
  ),
  data = d_session,
  start = list(
    a = rep(6, 6),
    sigma = 1
  )
)

precis(m_session, depth = 2)


post_session <- extract.samples(m_session)

combos <- c(
  "Inc1 F", "Inc1 M",
  "Inc2 F", "Inc2 M",
  "Inc3 F", "Inc3 M"
)

cols_6 <- rep(session_cols, each = 2)
lty_6  <- rep(c(1, 2), 3)

par(mfrow = c(1, 3))

for(i in 1:3) {
  
  id_f <- (i - 1) * 2 + 1
  id_m <- id_f + 1
  
  dens(post_session$a[, id_f],
       col = session_cols[i],
       lwd = 3,
       xlab = "Log total neglect",
       main = paste("Session", i))
  
  dens(post_session$a[, id_m],
       col = session_cols[i],
       lwd = 3,
       lty = 2,
       add = TRUE)
  
  legend("topright",
         legend = c("Female", "Male"),
         col = c(session_cols[i], session_cols[i]),
         lty = c(1, 2),
         lwd = 3,
         bty = "n")
}

par(mfrow = c(1,1))


cat("\n=== SEX DIFFERENCES WITHIN SESSION ===\n")

for(i in 1:3) {
  
  id_f <- (i - 1) * 2 + 1
  id_m <- id_f + 1
  
  diff_post <- post_session$a[, id_m] - post_session$a[, id_f]
  ratio_post <- exp(diff_post)
  
  cat(
    "\n", session_levels[i], "\n",
    "M - F =", round(mean(diff_post), 3),
    "| PI =", round(PI(diff_post), 3),
    "| P(M > F) =", round(mean(diff_post > 0), 3),
    "\nRatio =", round(mean(ratio_post), 2),
    "| PI =", round(PI(ratio_post), 2),
    "\n"
  )
}

# =========================================================
# HOURLY NEGLECT PROPORTION
# =========================================================
# =========================================================
# EXACT HOURLY SPLIT OF NEGLECT
# =========================================================

split_gap_to_hours <- function(gap_start, gap_end) {
  
  hour_starts <- seq(
    floor_date(gap_start, "hour"),
    floor_date(gap_end, "hour"),
    by = "1 hour"
  )
  
  tibble(hour_start = hour_starts) %>%
    mutate(
      overlap_start = pmax(hour_start, gap_start),
      overlap_end   = pmin(hour_start + hours(1), gap_end),
      neglect_sec   = as.numeric(difftime(overlap_end, overlap_start, units = "secs")),
      neglect_hour  = hour(hour_start)  # renamed column
    ) %>%
    filter(neglect_sec > 0)
}

G_hourly <- G %>%
  rowwise() %>%
  mutate(split = list(
    split_gap_to_hours(gap_start, gap_end)
  )) %>%
  unnest(split) %>%
  ungroup()

G_hourly <- G_hourly %>%
  mutate(date = as.Date(hour_start))

neg_hourly <- G_hourly %>%
  group_by(season, session, nest, responsible_sex, date, neglect_hour) %>%
  summarise(neglect_sec = sum(neglect_sec), .groups = "drop") %>%
  mutate(neglect_prop = pmin(neglect_sec / 3600, 1))


D_diurnal <- neg_hourly %>%
  mutate(
    sex_id        = ifelse(responsible_sex == "f", 1, 2),
    hour_circ     = neglect_hour / 24,
    neglect_logit = log((neglect_prop + 1e-4) / (1 - neglect_prop + 1e-4))
  )
# =========================================================
# DIURNAL NEGLECT PROPORTION (SPLINES, SEX-SPECIFIC)
# =========================================================
library(splines)

# ----------------------------
# DATA
# ----------------------------
neg_hourly_complete <- neg_hourly %>%
  # 1. Ensure every combination of nest/sex/date has all 24 hours
  group_by(season, session, nest, responsible_sex, date) %>%
  complete(neglect_hour = 0:23, fill = list(neglect_sec = 0)) %>%
  ungroup() %>%
  # 2. Recalculate proportions with the new zeros
  mutate(neglect_prop = pmin(neglect_sec / 3600, 1))

# Now create D_diurnal from the COMPLETE dataset
D_diurnal <- neg_hourly_complete %>%
  mutate(
    sex_id        = ifelse(responsible_sex == "f", 1, 2),
    hour_circ      = neglect_hour / 24,
    # logit(0) is -Inf, so the 1e-4 offset is CRITICAL here
    neglect_logit = log((neglect_prop + 1e-4) / (1 - neglect_prop + 1e-4))
  )

# spline basis
n_knots <- 5
knots <- seq(0, 1, length.out = n_knots)
B <- bs(D_diurnal$hour_circ, knots = knots[-c(1, n_knots)], degree = 3, intercept = FALSE)  # ← intercept=FALSE
n_basis <- ncol(B)
colnames(B) <- paste0("b", 1:n_basis)
D_diurnal <- cbind(D_diurnal, B)

# ----------------------------
# PRIOR PREDICTIVE CHECK
# ----------------------------
set.seed(123)
n_sim <- 100
day_seq <- seq(0, 23, length.out = 100)
hour_circ_seq <- day_seq / 24
B_seq <- bs(hour_circ_seq, knots = knots[-c(1, n_knots)], degree = 3, intercept = FALSE)

prior_f <- matrix(NA, nrow = n_sim, ncol = length(day_seq))
prior_m <- matrix(NA, nrow = n_sim, ncol = length(day_seq))

for (i in 1:n_sim) {
  a <- rnorm(2, -8, 2)    # ← was 0.3, now matches model prior
  b <- matrix(rnorm(n_basis * 2, 0, 0.8), ncol = 2)  # ← was 0.2, now matches model prior
  prior_f[i,] <- 1 / (1 + exp(-(a[1] + B_seq %*% b[,1])))
  prior_m[i,] <- 1 / (1 + exp(-(a[2] + B_seq %*% b[,2])))
}

# Check prior range is sensible
cat("Prior neglect proportion range - Female:",
    round(range(prior_f), 3), "\n")
cat("Prior neglect proportion range - Male:",
    round(range(prior_m), 3), "\n")

par(mfrow = c(1, 2))
matplot(day_seq, t(prior_f), type = 'l', col = col.alpha("darkred", 0.3),
        xlab = "Hour", ylab = "Neglect proportion",
        main = "Prior predictive - Female", ylim = c(0, 0.1))  # ← ylim to 1, data goes to 1
matplot(day_seq, t(prior_m), type = 'l', col = col.alpha("darkblue", 0.3),
        xlab = "Hour", ylab = "Neglect proportion",
        main = "Prior predictive - Male", ylim = c(0, 0.1))
par(mfrow = c(1, 1))

# Sprawdź wymiar! Powinien wyjść 7
n_basis <- ncol(B) 

# Przygotuj listę danych tak, by każda kolumna była osobnym wektorem
d_list_quap <- list(
  y = D_diurnal$neglect_logit,
  sex_id = as.integer(D_diurnal$sex_id)
)

# Dynamiczne dodawanie kolumn b1, b2... b7 do listy
for (i in 1:n_basis) {
  d_list_quap[[paste0("b", i)]] <- B[, i]
}

# Model z jawnym wypisaniem parametrów (dla 7 baz)
m_diurnal_final <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a[sex_id] + 
      bf1[sex_id]*b1 + bf2[sex_id]*b2 + bf3[sex_id]*b3 + 
      bf4[sex_id]*b4 + bf5[sex_id]*b5 + bf6[sex_id]*b6,
    
    # Intercepts (Baseline neglect around -8.5 logit)
    a[sex_id] ~ dnorm(-8, 2),
    
    # Spline coefficients (Effect of hour per sex)
    # I renamed them bf (basis factor) to distinguish from the data columns b1..b6
    bf1[sex_id] ~ dnorm(0, 0.8),
    bf2[sex_id] ~ dnorm(0, 0.8),
    bf3[sex_id] ~ dnorm(0, 0.8),
    bf4[sex_id] ~ dnorm(0, 0.8),
    bf5[sex_id] ~ dnorm(0, 0.8),
    bf6[sex_id] ~ dnorm(0, 0.8),
    
    sigma ~ dexp(1)
  ), 
  data = d_list_quap, 
  chains = 4, cores = 4, iter = 2000
)

# Now check again - depth=2 is required to see the [1] and [2] indices
precis(m_diurnal_final, depth = 2)


# ----------------------------
# POSTERIOR PREDICTIONS
# ----------------------------
# Important: We must use the same names (bf1, bf2...) used in the model alist
# but we map them to the B_seq (the actual basis values)
post_f <- link(m_diurnal_final, data = list(
  sex_id = rep(1, 100),
  b1 = B_seq[,1], b2 = B_seq[,2], b3 = B_seq[,3],
  b4 = B_seq[,4], b5 = B_seq[,5], b6 = B_seq[,6]
))

post_m <- link(m_diurnal_final, data = list(
  sex_id = rep(2, 100),
  b1 = B_seq[,1], b2 = B_seq[,2], b3 = B_seq[,3],
  b4 = B_seq[,4], b5 = B_seq[,5], b6 = B_seq[,6]
))

# Convert logit predictions to proportions: inv_logit = 1 / (1 + exp(-x))
mu_f <- apply(inv_logit(post_f), 2, mean)
mu_m <- apply(inv_logit(post_m), 2, mean)

PI_f <- apply(inv_logit(post_f), 2, PI, prob = 0.89)
PI_m <- apply(inv_logit(post_m), 2, PI, prob = 0.89)

# ----------------------------
# PLOT
# ----------------------------
# Setting ylim slightly above the max seen in the data if neglect is rare
# or keeping it at 1 for the full scale.
plot(NULL, xlim = c(0, 23), ylim = c(0,0.0005
                                  ),
     xlab = "Hour of Day", ylab = "Neglect Proportion",
     main = "Diurnal Neglect by Sex (Spline Fit)")

# Add uncertainty intervals
shade(PI_f, day_seq, col = col.alpha("darkred", 0.15))
shade(PI_m, day_seq, col = col.alpha("darkblue", 0.15))

# Add mean lines
lines(day_seq, mu_f, col = "darkred", lwd = 3)
lines(day_seq, mu_m, col = "darkblue", lwd = 3)

# Add raw data points (only if prop > 0 to see the events clearly, or use alpha)
points(jitter(D_diurnal$neglect_hour[D_diurnal$sex_id == 1], 0.3), 
       D_diurnal$neglect_prop[D_diurnal$sex_id == 1],
       pch = 16, col = col.alpha("darkred", 0.3), cex = 0.6)

points(jitter(D_diurnal$neglect_hour[D_diurnal$sex_id == 2], 0.3), 
       D_diurnal$neglect_prop[D_diurnal$sex_id == 2],
       pch = 1, col = col.alpha("darkblue", 0.3), cex = 0.6)

legend("topright", legend = c("Female", "Male"),
       col = c("darkred", "darkblue"), pch = c(16, 1), lwd = 3, bty = "n")

# ----------------------------
# SUMMARY & DIFFERENCES
# ----------------------------
cat("Female peak hour:", round(day_seq[which.max(mu_f)], 1), "\n")
cat("Male peak hour:",   round(day_seq[which.max(mu_m)], 1), "\n")

# Difference in proportions (Male - Female)
# We compute the difference on the proportion scale, not logit scale
diff_prop <- inv_logit(post_m) - inv_logit(post_f)
diff_mean <- apply(diff_prop, 2, mean)
diff_PI   <- apply(diff_prop, 2, PI, prob = 0.89)

cat("Max sex difference in neglect proportion (M-F):", round(max(abs(diff_mean)), 4), "\n")


##################### repeatability
library(rptR)

# ------------------------------
# Overall repeatability (Gamma)
# ------------------------------
S$is_neglecting <- ifelse(S$sum_neglect_sec > 0, 1, 0)

rep_binary <- rpt(
  is_neglecting ~ session + season + (1 | ringno),
  grname = "ringno",
  data = S,
  datatype = "Binary", # Rozkład dwumianowy (0/1)
  nboot = 1000, npermut = 1000
)

print(rep_binary)

# ------------------------------
# Repeatability by sex
# ------------------------------
S_male <- S %>% filter(sx == "m")
S_female <- S %>% filter(sx == "f")

rep_binary_male <- rpt(
  is_neglecting ~ session + season + (1 | ringno),
  grname = "ringno",
  data = S_male,
  datatype = "Binary", # Rozkład dwumianowy (0/1)
  nboot = 1000, npermut = 1000
)

rep_binary_female <- rpt(
  sum_neglect_sec ~ session + season + (1 | ringno),
  grname = "ringno",
  data = S_female,
  datatype = "Gamma",
  nboot = 1000,
  npermut = 1000
)

print(rep_binary_male)
print(rep_binary_female)


