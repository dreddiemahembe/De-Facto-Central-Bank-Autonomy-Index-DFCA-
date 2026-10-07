## Simulated country-year panel with a known latent "capture" factor.
## Used by test_simulate.R and demo_figures.R. Nothing here describes a real
## country: C01 to C06 are a simulated high-capture group, C07 onwards a
## low-capture group, and C01 carries an injected outlier run that stands in for
## a hyperinflation episode.
##
## Needs build_index.R to be sourced first (it defines the indicator list RAW).

simulate_panel <- function(seed = 20261002, n_countries = 20, years = 1985:2023) {
  set.seed(seed)
  N_C <- n_countries; YRS <- years
  iso  <- rep(sprintf("C%02d", 1:N_C), each = length(YRS))
  year <- rep(YRS, times = N_C)

  ## latent capture: higher for the high-capture group, drifting over time
  dominant <- rep(c(rep(1, 6), rep(0, N_C - 6)), each = length(YRS))
  tenure   <- ave(year, iso, FUN = function(v) seq_along(v))
  latent   <- 0.35 * dominant + 0.012 * tenure * dominant +
              rnorm(length(iso), 0, 0.15)
  latent   <- (latent - min(latent)) / (max(latent) - min(latent))

  df <- data.frame(
    iso, year, dominant, tenure, latent,
    ncg_gdp    = 8 * latent + rnorm(length(iso), 0, 1.2),
    dncg_gdp   = 2 * latent + rnorm(length(iso), 0, 0.8),
    cos_gdp    = 5 * latent + rnorm(length(iso), 0, 1.5),
    cb_contrib = 0.6 * latent + rnorm(length(iso), 0, 0.2)
  )
  df$breach <- rbinom(nrow(df), 1, pmin(0.95, pmax(0.02, latent)))

  ## an outlier run in C01, as a hyperinflation would produce
  hyp <- df$iso == "C01" & df$year %in% 2003:2008
  df$ncg_gdp[hyp]    <- df$ncg_gdp[hyp] * 12
  df$cb_contrib[hyp] <- df$cb_contrib[hyp] * 8

  ## 8% of cells missing at random, as in real balance-sheet data
  for (v in RAW) df[[v]][sample(nrow(df), round(0.08 * nrow(df)))] <- NA
  df
}
