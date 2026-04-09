
rm(list = ls())
library(dplyr)
library(rethinking)

# 1. Data Prep
E <- readRDS("./EDA_egg_neglect_index_neglect_table.RDS")

nest_gaps <- E %>%
  group_by(season, nest) %>%
  filter(n() == 3) %>%
  summarize(
    total_sum_gaps = sum(sum_gaps, na.rm = TRUE),
    incubation_duration = first(Incubation_Duration),
    .groups = "drop"
  ) %>%
  mutate(
    sum_gaps_k = total_sum_gaps / 1000,
    sum_gaps_c = sum_gaps_k - mean(sum_gaps_k)
  ) %>%
  filter(!is.na(incubation_duration))

d_quap <- list(
  incubation_duration = nest_gaps$incubation_duration,
  sum_gaps_c          = nest_gaps$sum_gaps_c
)

# 2. Your "Perfect" Prior Simulation
set.seed(123)
sum_gaps_seq <- seq(min(nest_gaps$sum_gaps_c), max(nest_gaps$sum_gaps_c), length.out = 50)
n_sim <- 50

plot(NULL, xlim = range(sum_gaps_seq), ylim = c(28, 32),
     xlab = "Total sum_gaps (centered, thousands)",
     ylab = "Incubation Duration",
     main = "Prior predictive: Incubation Duration vs Neglect")

for(i in 1:n_sim){
  a <- rnorm(1, 29, 0.5)   # Your original intercept
  b <- rnorm(1, 0.1, 0.05)  # Your original slope
  lines(sum_gaps_seq, a + b * sum_gaps_seq, col = col.alpha("grey", 0.6))
}
abline(h = c(28, 32), col = "red", lty = 2)

# 3. Model Fitting with your exact Priors
m_quap <- quap(
  alist(
    incubation_duration ~ dnorm(mu, sigma),
    mu <- a + b * sum_gaps_c,
    a ~ dnorm(29, 0.5),    # Matches your simulation
    b ~ dnorm(0.1, 0.05),  # Matches your simulation
    sigma ~ dexp(1)        # Switched to Exponential to avoid 'vmmin' non-finite errors
  ),
  data = d_quap,
  start = list(a = 29, b = 0.1, sigma = 1)
)

# 4. Summaries
cat("--- Parameter Estimates ---\n")
print(precis(m_quap))
# 1. Sample from the posterior to get uncertainty in the mean (mu)
post_mu <- link(m_quap, data = list(sum_gaps_c = sum_gaps_seq))
mu_mean <- apply(post_mu, 2, mean)
mu_PI   <- apply(post_mu, 2, PI, prob = 0.89)

# 2. Sample from the posterior to get uncertainty in the actual predictions (sigma included)
post_sim <- sim(m_quap, data = list(sum_gaps_c = sum_gaps_seq))
pred_PI  <- apply(post_sim, 2, PI, prob = 0.89)

# 3. Plot everything
plot(incubation_duration ~ sum_gaps_c, data = d_quap, 
     col = col.alpha(rangi2, 0.7), pch = 16,
     xlab = "Neglect (centered, thousands)", ylab = "Incubation Duration (days)")

# Add the MAP line
lines(sum_gaps_seq, mu_mean, lwd = 2)

# Add the uncertainty of the mean (the "hazy" region around the line)
shade(mu_PI, sum_gaps_seq, col = col.alpha("black", 0.15))

# Add the uncertainty of the predictions (where 89% of actual nests should fall)
shade(pred_PI, sum_gaps_seq, col = col.alpha("black", 0.05))

