## Test harness. Simulates a panel with a known latent "capture" factor, runs
## every component of build_index.R, prints the results and STOPS with an error
## if a check fails. Run from this folder:  Rscript test_simulate.R
##
## The data are simulated, so the recovery and validation checks show that the
## code works. They are not evidence about any real central bank.

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(f)) setwd(dirname(normalizePath(f)))

source("build_index.R")

failures <- character(0)
check <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok  " else "FAIL", label))
  if (!isTRUE(ok)) failures <<- c(failures, label)
}
errors <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

## ---- simulate the panel -----------------------------------------------------

source("simulate_panel.R")
df <- simulate_panel()
N_C <- length(unique(df$iso)); YRS <- sort(unique(df$year))

cat("=== panel:", nrow(df), "country-years,", N_C, "countries ===\n")

## ---- 1. input checks --------------------------------------------------------
cat("\n--- input checks ---\n")
check("a missing column is rejected", errors(build_index(df[, setdiff(names(df), "cos_gdp")])))
check("duplicate country-years are rejected", errors(build_index(rbind(df, df[1, ]))))

## ---- 2. the aggregation rules -----------------------------------------------
cat("\n--- aggregation rules ---\n")
M <- rbind(c(1, 1, 1, 0), c(0.5, 0.5, 0.5, 0.5), c(0, 0, 0, 0), c(1, 1, 1, 1))
colnames(M) <- names(DFCA_DIMS)
cap_ar <- aggregate_capture(M, method = "arithmetic")
cap_ge <- aggregate_capture(M, method = "geometric")
toy <- data.frame(case = c("three dimensions fully captured, one clean", "all 0.5",
                           "fully autonomous", "fully captured"),
                  DFCA_arithmetic = round(1 - cap_ar, 3), DFCA_geometric = round(1 - cap_ge, 3))
print(toy, row.names = FALSE)
check("geometric is never more autonomous than arithmetic", all(1 - cap_ge <= 1 - cap_ar + 1e-12))
check("a bank captured on 3 of 4 dimensions is not rated highly autonomous",
      (1 - cap_ge)[1] < 0.6)
check("ends of the scale are kept (0 stays 0, 1 stays 1)",
      isTRUE(all.equal((1 - cap_ge)[3:4], c(1, 0))))
check("pure geometric mean (floor = 0) matches the textbook value",
      isTRUE(all.equal(aggregate_scores(matrix(c(0.2, 0.8), 1), method = "geometric", floor = 0),
                       sqrt(0.2 * 0.8))))

## ---- 3. baseline index ------------------------------------------------------
cat("\n--- baseline index ---\n")
idx <- build_index(df)
cat("rows:", nrow(idx), "  non-missing DFCA:", sum(!is.na(idx$dfca)), "\n")
rec <- cor(idx$capture, df$latent, use = "complete.obs")
cat("correlation of capture score with the latent factor:", round(rec, 3), "\n")
check("DFCA lies in [0, 1]", all(idx$dfca >= 0 & idx$dfca <= 1, na.rm = TRUE))
check("DFCA = 1 - capture", isTRUE(all.equal(idx$dfca, 1 - idx$capture)))
check("orientation: capture rises with the latent factor", rec > 0.5)

## recovery under alternative choices, for the README
alt <- function(label, ...) {
  i <- build_index(df, ...)
  cat(sprintf("  %-40s cor with latent %.3f   median DFCA %.3f\n", label,
              cor(i$capture, df$latent, use = "complete.obs"), median(i$dfca, na.rm = TRUE)))
}
cat("alternatives:\n")
alt("arithmetic across dimensions", agg_across = "arithmetic")
for (fl in c(0, 0.05, 0.1, 0.25)) alt(sprintf("geometric, floor %.2f", fl), geo_floor = fl)

## ---- 4. diagnostics ---------------------------------------------------------
cat("\n--- diagnostics ---\n")
dg <- diagnostics(idx)
cat("dimension correlations:\n"); print(dg$dimension_correlations)
cat("PCA variance share:", dg$pca_variance_share[1:3], "\n")
cat("first component loadings:", dg$first_component_loadings, "\n")
cat("within-country share of variance:\n"); print(dg$within_country_share)
cat("complete country-years:", dg$complete_country_years, "of", nrow(idx), "\n")
off <- dg$dimension_correlations[upper.tri(dg$dimension_correlations)]
check("dimensions are related but not redundant (correlations below 0.8)", all(off < 0.8))
check("within-country shares lie in [0, 1]", all(dg$within_country_share >= 0 & dg$within_country_share <= 1))

## ---- 5. outliers ------------------------------------------------------------
cat("\n--- winsorising ---\n")
raw_max <- max(df$ncg_gdp, na.rm = TRUE)
win_max <- max(winsorise(df$ncg_gdp), na.rm = TRUE)
cat("ncg_gdp max before:", round(raw_max, 1), " after:", round(win_max, 1), "\n")
check("winsorising tames the hyperinflation run", win_max < raw_max / 5)

## ---- 6. coverage, reduced core, interpolation --------------------------------
cat("\n--- coverage and missing data ---\n")
au <- coverage_audit(df)
cat("countries:", nrow(au), "  median longest complete run:", median(au$longest_run), "years\n")
cm <- coverage_matrix(df)
check("coverage matrix has one row per country and one column per year",
      all(dim(cm) == c(N_C, length(YRS))))
check("coverage matrix counts are between 0 and 5", all(cm >= 0 & cm <= 5))

core <- build_index(df, dims = DFCA_DIMS_CORE)
both <- !is.na(core$dfca) & !is.na(idx$dfca)
cat("non-missing DFCA, full index:", sum(!is.na(idx$dfca)),
    "  reduced core:", sum(!is.na(core$dfca)), "\n")
cat("correlation of core with full index where both exist:",
    round(cor(core$dfca[both], idx$dfca[both]), 3), "\n")
check("the reduced core has wider coverage", sum(!is.na(core$dfca)) > sum(!is.na(idx$dfca)))
check("the reduced core tracks the full index", cor(core$dfca[both], idx$dfca[both]) > 0.5)

tp <- data.frame(iso = "T1", year = 2000:2009,
                 ncg_gdp = c(1, NA, 3, NA, NA, NA, 7, 8, NA, 10),
                 dncg_gdp = 1, cos_gdp = 1, cb_contrib = 1, breach = 0)
im <- impute_within_country(tp, vars = "ncg_gdp", max_gap = 2,
                            breaks = data.frame(iso = "T1", year = 2008))
check("a one-year gap is filled by linear interpolation", isTRUE(all.equal(im$ncg_gdp[2], 2)))
check("a gap longer than max_gap is left alone", all(is.na(im$ncg_gdp[4:6])))
check("a gap that spans a definitional break is left alone", is.na(im$ncg_gdp[9]))
check("imputed cells are flagged", identical(which(im$imp_ncg_gdp), 2L))

## ---- 7. sensitivity and rank bands ------------------------------------------
cat("\nrunning sensitivity analysis (12 specifications x 40 weight draws)...\n")
sens <- sensitivity(df, reps_weights = 40)
sens2 <- sensitivity(df, reps_weights = 40)
check("the sensitivity run is reproducible", isTRUE(all.equal(sens$dfca, sens2$dfca)))
cat("draws:", length(unique(sens$draw)), " rows:", nrow(sens), "\n\n")

rb <- rank_bands(sens, year_sel = 2015)
n_2015 <- sum(!is.na(idx$dfca[idx$year == 2015]))
cat("--- rank uncertainty, 2015 (rank 1 = most autonomous) ---\n")
print(head(rb[, c("iso", "rk", "rank_p05", "rank_p95", "dfca")], 10), row.names = FALSE)
cat("\ncountries ranked in 2015:", rb$n_ranked[1], "of", N_C, "\n")
cat("median width of the 90% rank interval:",
    round(median(rb$rank_p95 - rb$rank_p05), 2), "positions out of", rb$n_ranked[1], "\n")
check("rank bands cover exactly the countries with a 2015 value", rb$n_ranked[1] == n_2015)
check("each rank interval contains the median rank", all(rb$rank_p05 <= rb$rk & rb$rk <= rb$rank_p95))

lb <- level_bands(sens, "C01")
check("level bands bracket the median", all(lb$p05 <= lb$median & lb$median <= lb$p95))

## ---- 8. validation against a simulated external measure ----------------------
cat("\n--- validation (simulated external measures: circular by construction) ---\n")
ext <- data.frame(iso = df$iso, year = df$year,
                  cbie_index = pmax(0, pmin(1, 0.75 - 0.3 * df$latent +
                                              rnorm(nrow(df), 0, 0.12))),
                  inflation  = 3 + 25 * df$latent + rnorm(nrow(df), 0, 4))
v <- validate(idx[, c("iso", "year", "dfca")], ext)
cat("correlations:\n"); print(v$correlations)
cat("\ninflation on de jure index and DFCA:\n"); print(round(v$incremental, 4))
check("DFCA enters the inflation regression with a negative sign",
      v$incremental["dfca", "Estimate"] < 0)

## ---- result -----------------------------------------------------------------
cat("\n")
if (length(failures)) {
  cat("FAILED:", length(failures), "check(s)\n")
  quit(status = 1)
}
cat("All checks passed.\n")
