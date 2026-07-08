# =============================================================================
# DEMO: Mann-Kendall Trend Analysis + Bayesian Markov Chain
#       for Extreme Heat Wave Frequency Time Series
#
# IMPORTANT: This script uses SIMULATED data for demonstration purposes.
# It illustrates the two-layer analytical framework (trend + persistence)
# applied in the Benue Basin study, without disclosing proprietary results.
#
# Method 1 : Mann-Kendall non-parametric trend test + Sen's slope
# Method 2 : Bayesian Markov chain (annual states) with posterior inference
#
# Author  : KOUPNA II Higelin Saint-Clair
# Context : PhD Thesis — Climatology-Epidemiology
# AI note : Coding framework developed with AI as partner.
#           All scientific validation performed by the author.
# =============================================================================

# ── 0. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("Kendall", "trend", "MCMCpack", "ggplot2",
"patchwork", "tidyverse", "lubridate")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

set.seed(2026)

# ── 1. Simulate a 44-year HWF series with structural breakpoint ───────────────
# Mimics a Sahelian context: low frequency pre-2003, accelerating post-2003
years   <- 1981:2024
n       <- length(years)

# Pre-break period (1981-2002): low and stable frequency
hwf_pre  <- rpois(22, lambda = 0.4)

# Post-break period (2003-2024): increasing frequency trend
lambda_post <- seq(0.8, 2.5, length.out = 22)
hwf_post    <- rpois(22, lambda = lambda_post)

hwf <- c(hwf_pre, hwf_post)

df <- tibble(year = years, hwf = hwf,
             period = ifelse(year <= 2002, "1981-2002", "2003-2024"))

cat("=== Simulated HWF series ===\n")
cat("Pre-break  mean:", round(mean(hwf_pre),  2), "waves/yr\n")
cat("Post-break mean:", round(mean(hwf_post), 2), "waves/yr\n")
cat("Change         :", round((mean(hwf_post)/mean(hwf_pre)-1)*100, 1), "%\n\n")

# ── 2. Mann-Kendall test + Sen's slope ───────────────────────────────────────
mk_result  <- MannKendall(hwf)
sen_result <- sens.slope(hwf)

cat("=== Mann-Kendall Test ===\n")
cat("tau  =", round(mk_result$tau, 4),  "\n")
cat("S    =", mk_result$S,              "\n")
cat("p    =", round(mk_result$sl, 4),   "\n")
cat("Sen slope =", round(as.numeric(sen_result$estimates), 4),
"waves/year\n\n")

# ── 3. Pettitt change point test ──────────────────────────────────────────────
pett <- pettitt.test(hwf)
bp_year <- years[pett$estimate]

cat("=== Pettitt Change Point ===\n")
cat("Breakpoint year:", bp_year, "\n")
cat("p-value        :", round(pett$p.value, 4), "\n\n")

# ── 4. Bayesian Markov Chain — annual states ──────────────────────────────────
# State definition (adapt thresholds to your actual classification)
# S0 = no wave | S1 = 1 wave | S2 = 2+ waves
state_vec <- case_when(
  hwf == 0 ~ 0L,
  hwf == 1 ~ 1L,
  TRUE     ~ 2L
)

# Build observed transition counts
n_states <- 3
trans_counts <- matrix(0, n_states, n_states)
for (i in 1:(n - 1)) {
  from <- state_vec[i]  + 1
  to   <- state_vec[i+1] + 1
  trans_counts[from, to] <- trans_counts[from, to] + 1
}

cat("=== Observed Transition Counts ===\n")
rownames(trans_counts) <- colnames(trans_counts) <- c("S0","S1","S2")
print(trans_counts)

# Bayesian estimation: Dirichlet posterior for each row
# Prior: Dirichlet(1,1,1) — non-informative
n_mcmc   <- 20000
post_mats <- vector("list", n_states)

for (i in 1:n_states) {
  alpha_post  <- trans_counts[i, ] + 1   # add 1 = Dirichlet(1,1,1) prior
  post_mats[[i]] <- MCMCpack::rdirichlet(n_mcmc, alpha_post)
}

# Posterior summary: transition matrix medians
trans_post_median <- sapply(1:n_states, function(i)
  apply(post_mats[[i]], 2, median))
trans_post_median <- t(trans_post_median)

cat("\n=== Posterior Transition Matrix (medians) ===\n")
rownames(trans_post_median) <- colnames(trans_post_median) <- c("S0","S1","S2")
print(round(trans_post_median, 3))

# Stationary distribution
stat_dist_samples <- t(sapply(1:n_mcmc, function(s) {
  P <- rbind(post_mats[[1]][s,], post_mats[[2]][s,], post_mats[[3]][s,])
  eigen_decomp <- eigen(t(P))
  stat <- Re(eigen_decomp$vectors[, 1])
  stat / sum(stat)
}))

cat("\n=== Stationary Distribution (posterior median + 95% CI) ===\n")
for (j in 1:n_states) {
  cat(sprintf("S%d: %.3f [%.3f, %.3f]\n", j-1,
    median(stat_dist_samples[, j]),
    quantile(stat_dist_samples[, j], 0.025),
    quantile(stat_dist_samples[, j], 0.975)))
}

# ── 5. Plots ──────────────────────────────────────────────────────────────────
BLUE <- "#2E5496"; RED <- "#C0392B"; ORG <- "#E8A045"

# 5a — HWF chronology with breakpoint
df_period_means <- df %>%
  group_by(period) %>%
  summarise(mean_hwf = mean(hwf), .groups = "drop")

# Sen slope trend line
sen_val <- as.numeric(sen_result$estimates)
mid_yr  <- mean(years)
trend_line <- tibble(
  year = years,
  trend = median(hwf) + sen_val * (years - mid_yr)
)

p1 <- ggplot(df, aes(x = year, y = hwf, fill = period)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_line(data = trend_line, aes(x = year, y = trend),
            inherit.aes = FALSE,
            color = RED, linewidth = 1.1, linetype = "dashed") +
  geom_vline(xintercept = bp_year - 0.5,
             linetype = "dotted", color = "grey30", linewidth = 0.9) +
  annotate("text", x = bp_year - 0.5, y = max(hwf) * 0.92,
           label = paste0("Breakpoint\n", bp_year),
           hjust = 1.1, size = 3, color = "grey30", fontface = "italic") +
  annotate("text", x = 2018, y = max(hwf) * 0.88,
           label = sprintf("MK tau = %.3f***\nSen = %.3f waves/yr",
                           mk_result$tau, sen_val),
           size = 3.2, color = RED, hjust = 0, fontface = "bold") +
  scale_fill_manual(values = c("1981-2002" = "#7FABD4",
"2003-2024" = BLUE),
                    name = "Period") +
  scale_x_continuous(breaks = seq(1981, 2024, 5)) +
  labs(title    = "Annual heat wave frequency — simulated 44-year series",
       subtitle = "Mann-Kendall trend test | Sen's slope | Pettitt breakpoint",
       x = NULL, y = "Heat waves per year",
       caption  = "Simulated data for demonstration — not actual study results") +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", color = BLUE),
        plot.caption = element_text(color = "grey60", size = 8),
        legend.position = "top")

# 5b — Posterior transition probabilities (heatmap)
trans_df <- expand.grid(from = c("S0","S1","S2"),
                        to   = c("S0","S1","S2")) %>%
  mutate(prob = as.vector(t(trans_post_median)),
         label = sprintf("%.2f", prob))

p2 <- ggplot(trans_df, aes(x = to, y = from, fill = prob)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 5, fontface = "bold",
            color = ifelse(trans_df$prob > 0.5, "white", "grey20")) +
  scale_fill_gradient(low = "#EAF0FA", high = BLUE,
                      name = "Transition\nprobability",
                      limits = c(0, 1)) +
  scale_x_discrete(position = "top",
                   labels = c("S0\n(no wave)", "S1\n(1 wave)", "S2\n(≥2 waves)")) +
  scale_y_discrete(limits = rev,
                   labels = c("S2\n(≥2 waves)", "S1\n(1 wave)", "S0\n(no wave)")) +
  labs(title    = "Posterior Bayesian Markov transition matrix",
       subtitle = "Dirichlet posterior | 20,000 MCMC samples | Non-informative prior",
       x = "State at year t+1", y = "State at year t",
       caption  = "Simulated data for demonstration — not actual study results") +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", color = BLUE),
        plot.caption = element_text(color = "grey60", size = 8),
        axis.text    = element_text(size = 9),
        panel.grid   = element_blank())

# 5c — Stationary distribution with credible intervals
stat_df <- tibble(
  state  = c("S0 — No wave", "S1 — 1 wave", "S2 — ≥2 waves"),
  median = apply(stat_dist_samples, 2, median),
  lo     = apply(stat_dist_samples, 2, quantile, 0.025),
  hi     = apply(stat_dist_samples, 2, quantile, 0.975)
)

p3 <- ggplot(stat_df, aes(x = reorder(state, -median),
                           y = median, fill = state)) +
  geom_col(alpha = 0.85, width = 0.55) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.15, linewidth = 0.9, color = "grey30") +
  geom_text(aes(label = sprintf("%.1f%%", median * 100)),
            vjust = -1.8, fontface = "bold", size = 4, color = "grey20") +
  scale_fill_manual(values = c("#7FABD4", BLUE, RED), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  labs(title    = "Stationary distribution of annual heat wave states",
       subtitle = "Posterior median + 95% credible intervals",
       x = NULL, y = "Long-run probability",
       caption  = "Simulated data for demonstration — not actual study results") +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", color = BLUE),
        plot.caption = element_text(color = "grey60", size = 8))

# ── 6. Combine and export ─────────────────────────────────────────────────────
combined <- (p1) / (p2 | p3) +
  plot_annotation(
    title    = "Mann-Kendall + Bayesian Markov Chain — Analytical Framework Demo",
    subtitle = paste0(
"Framework: PhD thesis — Heat waves in the Benue River Basin, Cameroon\n",
"KOUPNA II Higelin Saint-Clair | 2026"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, color = BLUE),
      plot.subtitle = element_text(size = 9, color = "grey50"))
  )

ggsave("demo_mk_markov_output.png",
       combined, width = 14, height = 10, dpi = 300, bg = "white")

cat("\n✅ Figure exported: demo_mk_markov_output.png\n")
cat("   Mann-Kendall tau  :", round(mk_result$tau, 4), "\n")
cat("   Sen slope         :", round(sen_val, 4), "waves/yr\n")
cat("   Pettitt breakpoint:", bp_year, "\n")
