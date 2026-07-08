# =============================================================================
# DEMO: Bayesian GEV Analysis for Annual Temperature Maxima
# Method: Ratio-of-Uniforms with MDI prior (revdbayes package)
#
# IMPORTANT: This script uses SIMULATED data for demonstration purposes.
# It illustrates the analytical framework applied in the Benue Basin study
# without disclosing proprietary station records or classified results.
#
# Author  : KOUPNA II Higelin Saint-Clair
# Context : PhD Thesis — Climatology-Epidemiology
# AI note : Coding framework developed with AI as partner.
#           All scientific validation performed by the author.
# =============================================================================

# 0. Packages
pkgs <- c("revdbayes", "ggplot2", "patchwork")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

set.seed(2026)

# 1. Simulate annual maxima (Sahel-like climate, n = 44 years)
# NOTE: parameters below are ILLUSTRATIVE, not the actual study values
n_years <- 44
mu_sim  <- 42.0   # location (°C) — illustrative
sig_sim <-  1.2   # scale    (°C)
xi_sim  <-  0.15  # shape    — Fréchet family (xi > 0)

# GEV random variates via quantile function
u        <- runif(n_years)
tmax_sim <- mu_sim + sig_sim * (((-log(u))^(-xi_sim) - 1) / xi_sim)

cat("Simulated maxima summary:\n")
print(round(summary(tmax_sim), 2))

# 2. Bayesian GEV estimation (Ratio-of-Uniforms, MDI prior)
n_post    <- 50000
gev_bayes <- revdbayes::rpost_rcpp(
  n     = n_post,
  model = "gev",
  prior = revdbayes::set_prior(prior = "mdi", model = "gev"),
  data  = tmax_sim,
  nrep  = 100
)

post <- as.data.frame(gev_bayes$sim_vals)
names(post) <- c("mu", "sigma", "xi")

cat("\nPosterior quantiles:\n")
print(round(apply(post, 2, quantile, c(0.025, 0.5, 0.975)), 3))
cat("P(xi > 0) =", round(mean(post$xi > 0) * 100, 1), "% — Fréchet family\n")

# 3. Return level function
rl_fn <- function(p, mu, sigma, xi) {
  yp <- -log(1 - p)
  if (abs(xi) < 1e-6) mu - sigma * log(yp)
  else mu + sigma * (yp^(-xi) - 1) / xi
}

# 4. Return level curve
periods  <- 2:200
rl_curve <- lapply(periods, function(T) {
  p   <- 1 - 1/T
  rl  <- apply(post, 1, function(r) rl_fn(p, r[1], r[2], r[3]))
  data.frame(period = T,
             median = median(rl),
             lo     = quantile(rl, 0.025),
             hi     = quantile(rl, 0.975))
})
rl_curve <- do.call(rbind, rl_curve)

# 5. Plots
BLUE <- "#2E5496"; RED <- "#C0392B"

# 5a — Posterior densities
post_long <- data.frame(
  value = c(post$mu, post$sigma, post$xi),
  param = rep(c("mu — location (°C)",
"sigma — scale (°C)",
"xi — shape [Fréchet if > 0]"),
              each = n_post)
)

p1 <- ggplot(post_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = BLUE, alpha = 0.7, color = "white") +
  geom_density(color = RED, linewidth = 0.9) +
  geom_vline(
    data = subset(post_long, param == "xi — shape [Fréchet if > 0]"),
    aes(xintercept = 0),
    linetype = "dashed", color = "orange", linewidth = 0.8) +
  facet_wrap(~param, scales = "free", nrow = 1) +
  labs(
    title    = "Posterior distributions of GEV parameters",
    subtitle = paste0("Bayesian | Ratio-of-Uniforms | MDI prior | N = ",
                      format(n_post, big.mark = ","), " samples"),
    x = NULL, y = "Density",
    caption  = "Simulated data — not actual study results"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text   = element_text(face = "bold"),
        plot.title   = element_text(face = "bold", color = BLUE),
        plot.caption = element_text(color = "grey60", size = 8))

# 5b — Return level plot
obs_df <- data.frame(
  x = n_years / (n_years + 1 - rank(tmax_sim)),
  y = sort(tmax_sim)
)

p2 <- ggplot(rl_curve, aes(x = period)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = BLUE, alpha = 0.15) +
  geom_line(aes(y = median), color = BLUE, linewidth = 1.2) +
  geom_point(data = obs_df, aes(x = x, y = y),
             color = RED, size = 2.5, alpha = 0.8,
             shape = 21, fill = "white", stroke = 1.2) +
  geom_vline(xintercept = c(10, 50),
             linetype = "dotted", color = "grey40", linewidth = 0.7) +
  annotate("text", x = c(10, 50), y = min(rl_curve$lo) + 0.3,
           label = c("T10", "T50"), color = "grey30",
           size = 3.2, hjust = -0.2, fontface = "italic") +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200)) +
  labs(
    title    = "Return level plot — Bayesian GEV (Fréchet family)",
    subtitle = "Blue band = 95% credible interval | Dots = empirical maxima",
    x = "Return period (years, log scale)",
    y = "Return level (°C)",
    caption  = "Simulated data — not actual study results"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title   = element_text(face = "bold", color = BLUE),
        plot.caption = element_text(color = "grey60", size = 8))

# 6. Export
combined <- p1 / p2 +
  patchwork::plot_annotation(
    title    = "Bayesian GEV Analysis — Analytical Framework Demo",
    subtitle = paste0(
"PhD thesis: Heat waves in the Benue River Basin, Cameroon\n",
"KOUPNA II Higelin Saint-Clair | 2026"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, color = BLUE),
      plot.subtitle = element_text(size = 9, color = "grey50"))
  )

ggsave("demo_bayesian_gev_output.png",
       combined, width = 13, height = 9, dpi = 300, bg = "white")

cat("\n✅ Figure saved: demo_bayesian_gev_output.png\n")

