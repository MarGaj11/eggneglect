# PURPOSE: to check if distribution of hours in the data is equal, model n_gaps within 24 hours
library(devtools)
library(rethinking)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)


rm(list = ls())

#ngl_dt_trimmed <- readRDS("./EDA_egg_neglect_index_neglect_table_trimmed.RDS") 
ngl_dt <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS") 

#n_gaps <- counts, poisson distribution
#sum_gaps <- positive, continous, skewed(gamma-like)

#firstly we want to equal the hours - so all the videos have the same distribution of hours of the video
#Define a continuous time window that is well represented across nests
ngl_dt_hours <- ngl_dt %>%
  rowwise() %>%
  mutate(hour_seq = list(seq(
    floor_date(session_start, "hour"),
    floor_date(session_end_trimmed, "hour"),
    by = "1 hour"
  ))) %>%
  unnest(hour_seq) %>%
  ungroup() %>%
  mutate(hour = hour(hour_seq),
         vid_duration = as.numeric(difftime(session_end_trimmed, session_start, units = "mins")))

coverage <- ngl_dt_hours %>%
  group_by(season, session, nest, hour) %>%
  summarise(n_obs = n(), .groups = "drop") %>%   # how many times this hour appears per nest
  group_by(season, session, hour) %>%
  summarise(
    n_nests = n_distinct(nest),   # coverage
    total_obs = sum(n_obs),       # total counts (effort)
    mean_obs_per_nest = mean(n_obs)
  )

ggplot(coverage, aes(hour, n_nests)) +
  geom_col() +
  facet_grid(season ~ session)

ggplot(coverage, aes(hour, mean_obs_per_nest)) +
  geom_col() +
  geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
  facet_grid(season ~ session)


#ngl_dt_hours_trimmed <- ngl_dt_trimmed %>%
#  rowwise() %>%
#  mutate(hour_seq = list(seq(
#    floor_date(session_start_analysis, "hour"),
#    floor_date(session_end_trimmed, "hour"),
#    by = "1 hour"
#  ))) %>%
#  unnest(hour_seq) %>%
#  ungroup() %>%
#  mutate(hour = hour(hour_seq),
#         vid_duration = as.numeric(difftime(session_end_trimmed, session_start_analysis, units = "mins")))

#coverage_trimmed <- ngl_dt_hours_trimmed %>%
#  group_by(season, session, nest, hour) %>%
 # summarise(n_obs = n(), .groups = "drop") %>%   # how many times this hour appears per nest
#  group_by(season, session, hour) %>%
#  summarise(
#    n_nests = n_distinct(nest),   # coverage
 #   total_obs = sum(n_obs),       # total counts (effort)
#    mean_obs_per_nest = mean(n_obs)
#  )

#ggplot(coverage_trimmed, aes(hour, n_nests)) +
#  geom_col() +
#  facet_grid(season ~ session)

#ggplot(coverage_trimmed, aes(hour, mean_obs_per_nest)) +
#  geom_col() +
#  geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
#  facet_grid(season ~ session)

#I decided to not use trimmed data as test shows that the difference between hours is not significant
#and without trimming the duration of the video can be equal for every nest


#testing if distribution of hour is equal
hourly_data <- coverage %>%
  group_by(hour) %>%
  summarise(total_effort = sum(total_obs))

chi_test <- chisq.test(hourly_data$total_effort)

print(chi_test)

chi_results_table <- coverage %>%
  group_by(season, session) %>%
  summarise(
    chi_stat = chisq.test(total_obs)$statistic,
    p_value = chisq.test(total_obs)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    Equal = ifelse(p_value > 0.05, "yes", "no")
  )

print(chi_results_table)

#the distribution of hours/season/sessions in the whole data is equal due to the chi-square test, so we can see the 
#distribution of exits 

library(rethinking)
library(dplyr)
library(splines)

# =========================================================
# 0. DATA PREP
# =========================================================
E <- readRDS("./EDA_exit_hourly.RDS")

E <- E %>%
  mutate(
    incubation_stage = case_when(
      day_prior_hatch >= 21 & day_prior_hatch <= 30 ~ "early",
      day_prior_hatch >= 11 & day_prior_hatch <= 20 ~ "mid",
      day_prior_hatch >= 0  & day_prior_hatch <= 10 ~ "late",
      TRUE ~ NA_character_
    ),
    stage_id = as.integer(factor(incubation_stage,
                                 levels = c("early", "mid", "late"))),
    pair_id = as.integer(factor(Pair_status,
                                levels = c("New", "Old"))),
    stage_pair_id = as.integer(factor(
      paste(incubation_stage, Pair_status, sep = "_"),
      levels = c("early_New", "early_Old",
                 "mid_New", "mid_Old",
                 "late_New", "late_Old")
    ))
  )

hour_seq <- 0:23

# spline basis
basis <- bs(hour_seq,
            df = 6,
            degree = 3,
            intercept = FALSE)

basis_cols <- as.data.frame(basis)
colnames(basis_cols) <- paste0("B", 1:6)

basis_list <- as.list(basis_cols)

E_basis <- bind_cols(E, basis_cols[match(E$hour, hour_seq), ])

stage_names  <- c("early", "mid", "late")
stage_colors <- c("steelblue", "darkgreen", "firebrick")

# =========================================================
# 1. M2: OVERALL DAILY EXIT SPLINE
# =========================================================
d_m2 <- c(list(n_gaps_per_hour = E_basis$n_gaps_per_hour), as.list(E_basis[, paste0("B",1:6)]))

m2 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a + w1*B1 + w2*B2 + w3*B3 + w4*B4 + w5*B5 + w6*B6,
    a ~ dnorm(0,0.8),
    w1 ~ dnorm(0,0.5),
    w2 ~ dnorm(0,0.5),
    w3 ~ dnorm(0,0.5),
    w4 ~ dnorm(0,0.5),
    w5 ~ dnorm(0,0.5),
    w6 ~ dnorm(0,0.5)
  ),
  data = d_m2
)

# prior predictive
set.seed(10)
plot(NULL, xlim=c(0,23), ylim=c(0,5), xlab="Hour", ylab="Expected exits", xaxt="n")
axis(1, at=0:23)
for(i in 1:50){
  a <- rnorm(1,0,0.8)
  ws <- rnorm(6,0,0.5)
  mu <- exp(a + as.matrix(basis_cols) %*% ws)
  lines(hour_seq, mu, col=col.alpha("black",0.25))
}

# posterior fit
pred <- link(m2, data = basis_list)
mu <- apply(pred, 2, mean)
PI_mu <- apply(pred, 2, PI, prob=0.89)

Eb <- E %>% group_by(hour) %>% summarise(n_gaps_per_hour = mean(n_gaps_per_hour), .groups="drop")
plot(Eb$hour, Eb$n_gaps_per_hour, pch=16, xlab="Hour", ylab="Mean exits")
shade(PI_mu, hour_seq)
lines(hour_seq, mu, lwd=2)

# =========================================================
# 2. M3: STAGE-SPECIFIC SPLINES
# =========================================================
E_fit3 <- E_basis %>% filter(!is.na(stage_id))

d_m3 <- c(
  list(n_gaps_per_hour = E_fit3$n_gaps_per_hour,
       stage_id = E_fit3$stage_id),
  as.list(E_fit3[, paste0("B",1:6)])
)

m3 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[stage_id] +
      w1[stage_id]*B1 + w2[stage_id]*B2 + w3[stage_id]*B3 +
      w4[stage_id]*B4 + w5[stage_id]*B5 + w6[stage_id]*B6,
    a[stage_id] ~ dnorm(0,0.8),
    w1[stage_id] ~ dnorm(0,0.5),
    w2[stage_id] ~ dnorm(0,0.5),
    w3[stage_id] ~ dnorm(0,0.5),
    w4[stage_id] ~ dnorm(0,0.5),
    w5[stage_id] ~ dnorm(0,0.5),
    w6[stage_id] ~ dnorm(0,0.5)
  ),
  data = d_m3,
  start = list(
    a=rep(0,3), w1=rep(0,3), w2=rep(0,3),
    w3=rep(0,3), w4=rep(0,3), w5=rep(0,3), w6=rep(0,3)
  )
)

Eg <- E_fit3 %>% group_by(stage_id, hour) %>% summarise(mean_gaps=mean(n_gaps_per_hour), .groups="drop")

plot(NULL, xlim=c(0,23), ylim=c(0,0.5), xlab="Hour", ylab="Expected exits", xaxt="n")
axis(1, at=0:23)

mu_list_m3 <- list()
for(s in 1:3){
  pred_s <- link(m3, data = c(list(stage_id = rep(s,24)), basis_list))
  mu_s <- apply(pred_s,2,mean)
  PI_s <- apply(pred_s,2,PI, prob=0.89)
  mu_list_m3[[stage_names[s]]] <- apply(pred_s,1,mean)
  shade(PI_s, hour_seq, col=col.alpha(stage_colors[s],0.15))
  lines(hour_seq, mu_s, lwd=2, col=stage_colors[s])
  pts <- Eg %>% filter(stage_id==s)
  points(pts$hour, pts$mean_gaps, pch=16, col=stage_colors[s], cex=0.8)
}
legend("topright", legend=stage_names, col=stage_colors, lwd=2, pch=16)

# posterior distributions
xrange <- range(unlist(mu_list_m3))
plot(NULL, xlim=xrange, ylim=c(0,120), xlab="Expected exits/hour", ylab="Density",
     main="Posterior distributions by incubation stage")
for(i in 1:3){
  dens(mu_list_m3[[i]], col=stage_colors[i], lwd=2, add=TRUE)
}
legend("topright", legend=stage_names, col=stage_colors, lwd=2)

# contrasts
par(mfrow=c(1,3))
diff_early_mid <- mu_list_m3[["early"]] - mu_list_m3[["mid"]]
diff_early_late <- mu_list_m3[["early"]] - mu_list_m3[["late"]]
diff_mid_late <- mu_list_m3[["mid"]] - mu_list_m3[["late"]]

dens(diff_early_mid, col="purple", lwd=2, xlab="Difference", main="Early - Mid")
abline(v=0,lty=2)
dens(diff_early_late, col="orange", lwd=2, xlab="Difference", main="Early - Late")
abline(v=0,lty=2)
dens(diff_mid_late, col="brown", lwd=2, xlab="Difference", main="Mid - Late")
abline(v=0,lty=2)
par(mfrow=c(1,1))

cat("Early vs Mid: ", round(mean(diff_early_mid),3), PI(diff_early_mid), round(mean(diff_early_mid>0),2), "\n")
cat("Early vs Late:", round(mean(diff_early_late),3), PI(diff_early_late), round(mean(diff_early_late>0),2), "\n")
cat("Mid vs Late:  ", round(mean(diff_mid_late),3), PI(diff_mid_late), round(mean(diff_mid_late>0),2), "\n")

# peak hour summaries
for(s in 1:3){
  pred_s <- link(m3, data = c(list(stage_id=rep(s,24)), basis_list))
  peak_hour <- apply(pred_s,1,which.max)-1
  cat(stage_names[s], "peak:", round(mean(peak_hour),2), PI(peak_hour), "\n")
}

# =========================================================
# 3. M4: STAGE × PAIR STATUS SPLINES (Old vs new)
# =========================================================
E_fit4 <- E_basis %>% filter(!is.na(stage_pair_id))

d_m4 <- c(
  list(n_gaps_per_hour = E_fit4$n_gaps_per_hour,
       stage_pair_id = E_fit4$stage_pair_id),
  as.list(E_fit4[, paste0("B",1:6)])
)

m4 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    log(lambda) <- a[stage_pair_id] +
      w1[stage_pair_id]*B1 + w2[stage_pair_id]*B2 + w3[stage_pair_id]*B3 +
      w4[stage_pair_id]*B4 + w5[stage_pair_id]*B5 + w6[stage_pair_id]*B6,
    a[stage_pair_id] ~ dnorm(0,1),
    w1[stage_pair_id] ~ dnorm(0,0.4),
    w2[stage_pair_id] ~ dnorm(0,0.4),
    w3[stage_pair_id] ~ dnorm(0,0.4),
    w4[stage_pair_id] ~ dnorm(0,0.4),
    w5[stage_pair_id] ~ dnorm(0,0.4),
    w6[stage_pair_id] ~ dnorm(0,0.4)
  ),
  data = d_m4,
  start = list(
    a=rep(0,6), w1=rep(0,6), w2=rep(0,6),
    w3=rep(0,6), w4=rep(0,6), w5=rep(0,6), w6=rep(0,6)
  )
)

combos <- c("early_New", "early_Old", "mid_New", "mid_Old", "late_New", "late_Old")
stage_cols <- c("steelblue","steelblue","darkgreen","darkgreen","firebrick","firebrick")
stage_lty <- c(1,2,1,2,1,2)

plot(NULL, xlim=c(0,23), ylim=c(0,0.5), xlab="Hour", ylab="Expected exits", xaxt="n",
     main="Exit rate by incubation stage and pair status")
axis(1, at=0:23)

for(s in 1:6){
  pred_s <- link(m4, data = c(list(stage_pair_id = rep(s,24)), basis_list))
  mu_s <- apply(pred_s,2,mean)
  PI_s <- apply(pred_s,2,PI, prob=0.89)
  shade(PI_s, hour_seq, col=col.alpha(stage_cols[s],0.1))
  lines(hour_seq, mu_s, lwd=2, col=stage_cols[s], lty=stage_lty[s])
}
legend("topright", legend=combos, col=stage_cols, lty=stage_lty, lwd=2, cex=0.8)


# =========================================================
# 3. M5: SEASONAL DIFFERENCES IN HOURLY EXIT RATE
# =========================================================
# Goal:
# Test whether exit rate and diel timing differ among seasons
# using season-specific Poisson spline curves.
# =========================================================


# =========================================================
# 3.1 PREPARE SEASONAL SPLINE DATA
# =========================================================

season_levels <- c("2019", "2020", "2021", "2023", "2025")
season_cols   <- c("steelblue", "darkgreen", "firebrick",
                   "goldenrod", "purple")

E_basis <- E %>%
  mutate(
    season_id = as.integer(factor(season, levels = season_levels))
  ) %>%
  bind_cols(basis_cols[match(E$hour, hour_seq), ])

n_seasons <- length(season_levels)

d_m5 <- list(
  n_gaps_per_hour = E_basis$n_gaps_per_hour,
  season_id = E_basis$season_id,
  B1 = E_basis$B1,
  B2 = E_basis$B2,
  B3 = E_basis$B3,
  B4 = E_basis$B4,
  B5 = E_basis$B5,
  B6 = E_basis$B6
)


# =========================================================
# 3.2 FIT SEASON-SPECIFIC SPLINE MODEL
# =========================================================

m5 <- quap(
  alist(
    n_gaps_per_hour ~ dpois(lambda),
    
    log(lambda) <- a[season_id] +
      w1[season_id]*B1 +
      w2[season_id]*B2 +
      w3[season_id]*B3 +
      w4[season_id]*B4 +
      w5[season_id]*B5 +
      w6[season_id]*B6,
    
    a[season_id]  ~ dnorm(0, 0.8),
    w1[season_id] ~ dnorm(0, 0.5),
    w2[season_id] ~ dnorm(0, 0.5),
    w3[season_id] ~ dnorm(0, 0.5),
    w4[season_id] ~ dnorm(0, 0.5),
    w5[season_id] ~ dnorm(0, 0.5),
    w6[season_id] ~ dnorm(0, 0.5)
  ),
  data = d_m5
)

precis(m5, depth = 2)


# =========================================================
# 3.3 PLOT FITTED DIEL CURVES + RAW DATA
# =========================================================

E_raw <- E_basis %>%
  group_by(season, hour) %>%
  summarise(mean_exit = mean(n_gaps_per_hour), .groups = "drop")

plot(NULL,
     xlim = c(0, 23),
     ylim = c(0, 0.4),
     xlab = "Hour",
     ylab = "Expected exits",
     xaxt = "n",
     main = "Seasonal diel exit curves")

axis(1, at = 0:23)

for(s in 1:n_seasons){
  
  pred <- link(
    m5,
    data = list(
      season_id = rep(s, 24),
      B1 = basis_cols$B1,
      B2 = basis_cols$B2,
      B3 = basis_cols$B3,
      B4 = basis_cols$B4,
      B5 = basis_cols$B5,
      B6 = basis_cols$B6
    )
  )
  
  mu <- apply(pred, 2, mean)
  
  # fitted spline
  lines(hour_seq, mu,
        lwd = 3,
        col = season_cols[s])
  
  # raw observed means
  raw_s <- E_raw %>%
    filter(season == season_levels[s])
  
  points(raw_s$hour,
         raw_s$mean_exit,
         pch = 16,
         cex = 0.8,
         col = season_cols[s])
}

legend("topright",
       legend = season_levels,
       col = season_cols,
       lwd = 3,
       pch = 16,
       bty = "n")


# =========================================================
# 3.4 POSTERIOR DISTRIBUTIONS OF EXIT RATE BY SEASON
# =========================================================

mu_list_m5 <- list()

for(s in 1:n_seasons){
  
  pred <- link(
    m5,
    data = list(
      season_id = rep(s, 24),
      B1 = basis_cols$B1,
      B2 = basis_cols$B2,
      B3 = basis_cols$B3,
      B4 = basis_cols$B4,
      B5 = basis_cols$B5,
      B6 = basis_cols$B6
    )
  )
  
  mu_list_m5[[ season_levels[s] ]] <- rowMeans(pred)
}

xrange <- range(unlist(mu_list_m5))

plot(NULL,
     xlim = xrange,
     ylim = c(0, 120),
     xlab = "Expected exits/hour",
     ylab = "Density",
     main = "Posterior distributions by season")

for(i in 1:n_seasons){
  dens(mu_list_m5[[i]],
       col = season_cols[i],
       lwd = 2,
       add = TRUE)
}

legend("topright",
       legend = season_levels,
       col = season_cols,
       lwd = 2,
       bty = "n")


# =========================================================
# 3.5 PAIRWISE SEASONAL CONTRASTS
# =========================================================

season_pairs <- combn(season_levels, 2, simplify = FALSE)

par(mfrow = c(2, 5))

for(pair in season_pairs){
  
  diff <- mu_list_m5[[ pair[1] ]] -
    mu_list_m5[[ pair[2] ]]
  
  dens(diff,
       lwd = 2,
       xlab = "Difference",
       main = paste(pair[1], "-", pair[2]))
  
  abline(v = 0, lty = 2)
  
  cat(pair[1], "vs", pair[2], ": ",
      round(mean(diff), 3),
      PI(diff),
      round(mean(diff > 0), 2), "\n")
}

par(mfrow = c(1, 1))


# =========================================================
# 3.6 PEAK EXIT HOUR BY SEASON
# =========================================================

for(s in 1:n_seasons){
  
  pred_s <- link(
    m5,
    data = list(
      season_id = rep(s, 24),
      B1 = basis_cols$B1,
      B2 = basis_cols$B2,
      B3 = basis_cols$B3,
      B4 = basis_cols$B4,
      B5 = basis_cols$B5,
      B6 = basis_cols$B6
    )
  )
  
  peak_hour <- apply(pred_s, 1, which.max) - 1
  
  cat(season_levels[s],
      "peak:",
      round(mean(peak_hour), 2),
      PI(peak_hour), "\n")
}


# =========================================================
# 4. M6: PAIR SIMILARITY AND TOTAL NUMBER OF EXITS
# =========================================================
# Goal:
# Test whether total number of nest exits (n_gaps)
# changes linearly with pair similarity.
# =========================================================

ext_dt <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")



# keep complete cases only
D_sim <- ext_dt %>%
  filter(
    !is.na(Similarity_Index),
    !is.na(n_gaps),
    session == "incubation1"
  )

# =========================================================
# 4.1 FIT POISSON LINEAR MODEL
# =========================================================

m6 <- quap(
  alist(
    n_gaps ~ dpois(lambda),
    log(lambda) <- a + b*Similarity_Index,
    
    a ~ dnorm(0, 0.8),
    b ~ dnorm(0, 0.5)
  ),
  data = D_sim
)

precis(m6)

# =========================================================
# 4.2 PLOT RELATIONSHIP
# =========================================================

plot(D_sim$Similarity_Index,
     D_sim$n_gaps,
     pch = 16,
     xlab = "Pair similarity index",
     ylab = "Total number of exits")

xseq <- seq(0, 1, length.out = 100)

pred <- link(
  m6,
  data = list(Similarity_Index = xseq)
)

mu <- apply(pred, 2, mean)
PI_mu <- apply(pred, 2, PI, prob = 0.89)

shade(PI_mu, xseq)
lines(xseq, mu, lwd = 3)

#difference between 0 and 1
post <- extract.samples(m6)

lambda_0 <- exp(post$a + post$b*0)
lambda_1 <- exp(post$a + post$b*1)

diff_sim <- lambda_1 - lambda_0

mean(diff_sim)
PI(diff_sim)
mean(diff_sim > 0)
