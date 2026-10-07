## ---------------------------------------------------------------------------
## De facto central bank autonomy index (DFCA): reference implementation
## Version 1.1, October 2026
##
## Input:  a country-year panel (data frame) with columns iso, year and the raw
##         indicator columns named in DFCA_DIMS
## Output: the index, its dimension scores, sensitivity bands and diagnostics
##
## Status: tested on simulated data only (see test_simulate.R). No result for a
## real country has been produced with this code.
##
## Author: Dr Eddie Mahembe, Underhill Corporate Solutions
##
## Design choices are explained in README.md; CHANGES.md lists what differs from
## version 1.0 of the technical note.
## ---------------------------------------------------------------------------

## ---- 0. Settings -----------------------------------------------------------

## Every raw indicator is oriented so that HIGHER = MORE fiscal capture.
## The published index is reported as autonomy: DFCA = 1 - capture.

DFCA_DIMS <- list(direct       = c("ncg_gdp", "dncg_gdp"),  # net claims on government, stock and flow
                  quasifiscal  = "cos_gdp",                 # claims on other sectors incl. public corporations
                  monetisation = "cb_contrib",              # CB credit to government in base-money growth
                  compliance   = "breach")                  # net lending positive where the statute prohibits it

## Reduced core (technical note, section 6): two indicators, much wider coverage
DFCA_DIMS_CORE <- list(direct       = "ncg_gdp",
                       monetisation = "cb_contrib")

## Names used by the technical note and by version 1.0
DIMS <- DFCA_DIMS
RAW  <- unlist(DFCA_DIMS, use.names = FALSE)

## ---- 1. Utilities ----------------------------------------------------------

check_panel <- function(df, dims = DFCA_DIMS) {
  vars <- unlist(dims, use.names = FALSE)
  need <- c("iso", "year", vars)
  if (!nrow(df)) stop("The panel has no rows.")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  if (anyDuplicated(df[, c("iso", "year")]))
    stop("Duplicate iso-year rows. The panel needs one row per country and year.")
  bad <- vars[!vapply(df[vars], is.numeric, logical(1))]
  if (length(bad)) stop("Non-numeric indicator columns: ", paste(bad, collapse = ", "))
  invisible(TRUE)
}

winsorise <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

## Normalisation: three options, all mapping to [0, 1] with higher = more capture
normalise <- function(x, method = c("minmax", "rank", "zscore_cdf"), group = NULL) {
  method <- match.arg(method)
  f <- function(v) {
    if (all(is.na(v))) return(v)
    if (isTRUE(sd(v, na.rm = TRUE) == 0) || sum(!is.na(v)) < 2) {
      warning("an indicator has no variation in the normalisation sample; returning NA")
      return(rep(NA_real_, length(v)))
    }
    switch(method,
      minmax     = (v - min(v, na.rm = TRUE)) /
                   (max(v, na.rm = TRUE) - min(v, na.rm = TRUE)),
      rank       = (rank(v, na.last = "keep", ties.method = "average") - 1) /
                   (sum(!is.na(v)) - 1),
      zscore_cdf = pnorm((v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
    )
  }
  if (is.null(group)) f(x) else ave(x, group, FUN = f)
}

## Weighted mean of the columns of M. Scores must lie in [0, 1].
##   arithmetic: full compensation between columns.
##   geometric : penalises low columns. Scores are shifted into [floor, 1] so a
##               zero does not force the result to zero, then shifted back, so
##               0 stays 0 and 1 stays 1. floor = 0 gives the pure geometric mean.
aggregate_scores <- function(M, w = NULL, method = c("arithmetic", "geometric"),
                             floor = 0.1) {
  method <- match.arg(method)
  if (is.null(w)) w <- rep(1 / ncol(M), ncol(M))
  if (length(w) != ncol(M)) stop("weights must have one value per column")
  w <- w / sum(w)
  if (method == "arithmetic") return(as.numeric(M %*% w))
  S <- floor + (1 - floor) * M
  S <- pmax(S, 1e-12)
  g <- exp(as.numeric(log(S) %*% w))
  (g - floor) / (1 - floor)
}

## Aggregate CAPTURE scores. The geometric rule is applied to autonomy
## (1 - capture): one dimension of full capture pulls autonomy down hard, and a
## clean record elsewhere cannot offset it. Arithmetic gives the same answer in
## either direction.
aggregate_capture <- function(M, w = NULL, method = c("arithmetic", "geometric"),
                              floor = 0.1) {
  method <- match.arg(method)
  if (method == "arithmetic") aggregate_scores(M, w, "arithmetic")
  else 1 - aggregate_scores(1 - M, w, "geometric", floor)
}

## ---- 2. Build one version of the index -------------------------------------

build_index <- function(df,
                        dims = DFCA_DIMS,
                        norm = "minmax",
                        agg_within = "arithmetic",
                        agg_across = "geometric",
                        dim_weights = NULL,
                        winsor = 0.01,
                        within_country = FALSE,
                        geo_floor = 0.1) {

  check_panel(df, dims)
  vars <- unlist(dims, use.names = FALSE)
  d <- df
  for (v in vars) {
    x <- winsorise(d[[v]], winsor)
    d[[paste0("n_", v)]] <- normalise(x, norm,
                                      group = if (within_country) d$iso else NULL)
  }

  dim_scores <- sapply(names(dims), function(k) {
    cols <- paste0("n_", dims[[k]])
    M <- as.matrix(d[, cols, drop = FALSE])
    if (ncol(M) == 1) as.numeric(M)
    else aggregate_capture(M, method = agg_within, floor = geo_floor)
  })
  dim_scores <- as.data.frame(dim_scores)

  capture <- aggregate_capture(as.matrix(dim_scores), w = dim_weights,
                               method = agg_across, floor = geo_floor)

  cbind(d[, c("iso", "year")], dim_scores, capture = capture, dfca = 1 - capture)
}

## ---- 3. Sensitivity: the index is a set of choices, so report the set -------

rdirichlet_one <- function(n, alpha = 2) {
  g <- rgamma(n, shape = alpha, rate = 1)
  g / sum(g)
}

## 12 specifications (3 normalisations x 2 within x 2 across) times reps_weights
## random weight vectors on the simplex. Returns one row per country-year-draw.
sensitivity <- function(df, dims = DFCA_DIMS, reps_weights = 200,
                        seed = 20261002, geo_floor = 0.1) {
  set.seed(seed)
  grid <- expand.grid(norm = c("minmax", "rank", "zscore_cdf"),
                      agg_within = c("arithmetic", "geometric"),
                      agg_across = c("arithmetic", "geometric"),
                      stringsAsFactors = FALSE)
  res <- list(); k <- 0
  for (i in seq_len(nrow(grid))) {
    for (r in seq_len(reps_weights)) {
      w <- rdirichlet_one(length(dims))
      idx <- build_index(df, dims = dims, norm = grid$norm[i],
                         agg_within = grid$agg_within[i],
                         agg_across = grid$agg_across[i],
                         dim_weights = w, geo_floor = geo_floor)
      k <- k + 1
      res[[k]] <- data.frame(iso = idx$iso, year = idx$year, dfca = idx$dfca,
                             draw = k)
    }
  }
  do.call(rbind, res)
}

## Country-level rank uncertainty for one year (rank 1 = most autonomous).
## Only countries with a value in that year can be ranked; the number ranked is
## returned in the n_ranked column, because the intervals are positions out of
## that number, not out of all countries.
rank_bands <- function(sens, year_sel) {
  s <- sens[sens$year == year_sel & !is.na(sens$dfca), ]
  if (!nrow(s)) stop("no non-missing values in year ", year_sel)
  s$rk <- ave(s$dfca, s$draw, FUN = function(v) rank(-v, ties.method = "average"))
  med <- aggregate(cbind(rk, dfca) ~ iso, data = s, FUN = median)
  lo  <- aggregate(rk ~ iso, data = s, FUN = function(v) quantile(v, 0.05))
  hi  <- aggregate(rk ~ iso, data = s, FUN = function(v) quantile(v, 0.95))
  names(lo)[2] <- "rank_p05"; names(hi)[2] <- "rank_p95"
  out <- merge(merge(med, lo, by = "iso"), hi, by = "iso")
  out$n_ranked <- nrow(out)
  out[order(out$rk), ]
}

## Index level with its uncertainty band for one country over time
level_bands <- function(sens, iso_sel) {
  s <- sens[sens$iso == iso_sel & !is.na(sens$dfca), ]
  agg <- aggregate(dfca ~ year, data = s, FUN = function(v)
    c(median = median(v), p05 = quantile(v, 0.05), p95 = quantile(v, 0.95)))
  data.frame(year = agg$year, median = agg$dfca[, 1],
             p05 = agg$dfca[, 2], p95 = agg$dfca[, 3])
}

## ---- 4. Coverage audit (technical note, section 4.3: do this first) ---------

## Rows with a missing country-year are added so every country has every year.
complete_panel <- function(df, years = NULL) {
  if (is.null(years)) years <- seq(min(df$year), max(df$year))
  full <- expand.grid(iso = unique(df$iso), year = years, stringsAsFactors = FALSE)
  merge(full, df, by = c("iso", "year"), all.x = TRUE)
}

## Per country: first and last year, the longest unbroken run of years in which
## EVERY indicator is observed, and the share of years observed per indicator.
coverage_audit <- function(df, dims = DFCA_DIMS) {
  vars <- unlist(dims, use.names = FALSE)
  df <- complete_panel(df)
  df <- df[order(df$iso, df$year), ]
  one <- function(g) {
    ok <- complete.cases(g[vars])
    r <- rle(ok)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    best <- which(r$values)
    if (length(best)) {
      b <- best[which.max(r$lengths[best])]
      run_len <- r$lengths[b]; run_from <- g$year[starts[b]]; run_to <- g$year[ends[b]]
    } else { run_len <- 0; run_from <- NA; run_to <- NA }
    c(setNames(vapply(vars, function(v) mean(!is.na(g[[v]])), numeric(1)), paste0("share_", vars)),
      complete_years = sum(ok), longest_run = run_len, run_from = run_from, run_to = run_to)
  }
  parts <- lapply(split(df, df$iso), one)
  data.frame(iso = names(parts), do.call(rbind, parts), row.names = NULL)
}

## Country-by-year count of observed indicators (0 to number of indicators)
coverage_matrix <- function(df, dims = DFCA_DIMS) {
  vars <- unlist(dims, use.names = FALSE)
  df <- complete_panel(df)
  df$n_obs <- rowSums(!is.na(df[vars]))
  xtabs(n_obs ~ iso + year, data = df)
}

## ---- 5. Missing data (technical note, section 6) ------------------------------

## Linear interpolation inside a country series, for gaps of at most max_gap
## years with an observed value on both sides. Never across a definitional
## break: breaks is a data frame (iso, year) giving the first year measured under
## the new definition. Every interpolated cell is flagged in imp_<variable>.
impute_within_country <- function(df, vars = RAW, max_gap = 2, breaks = NULL) {
  df <- complete_panel(df)
  df <- df[order(df$iso, df$year), ]
  for (v in vars) df[[paste0("imp_", v)]] <- FALSE
  for (iso_i in unique(df$iso)) {
    rows <- which(df$iso == iso_i)
    yrs <- df$year[rows]
    bk <- if (is.null(breaks)) numeric(0) else breaks$year[breaks$iso == iso_i]
    for (v in vars) {
      x <- df[[v]][rows]
      obs <- which(!is.na(x))
      if (length(obs) < 2) next
      for (j in seq_len(length(obs) - 1)) {
        a <- obs[j]; b <- obs[j + 1]
        gap <- b - a - 1
        if (gap < 1 || gap > max_gap) next
        if (any(bk > yrs[a] & bk <= yrs[b])) next
        fill <- (a + 1):(b - 1)
        x[fill] <- x[a] + (x[b] - x[a]) * (fill - a) / (b - a)
        df[[paste0("imp_", v)]][rows[fill]] <- TRUE
      }
      df[[v]][rows] <- x
    }
  }
  df
}

## ---- 6. Validation ---------------------------------------------------------

## external needs columns iso, year, cbie_index (de jure) and inflation.
validate <- function(idx, external) {
  m <- merge(idx, external, by = c("iso", "year"))
  num <- m[, sapply(m, is.numeric)]
  list(
    n = nrow(m),
    correlations = round(cor(num, use = "pairwise.complete.obs"), 3),
    ## Does the de facto index add anything beyond the de jure index?
    incremental = summary(lm(inflation ~ cbie_index + dfca, data = m))$coefficients
  )
}

## ---- 7. Diagnostics on dimension structure ---------------------------------

## Share of the total variance of x that lies within countries.
within_share <- function(x, g) {
  ok <- !is.na(x)
  x <- x[ok]; g <- g[ok]
  tot <- sum((x - mean(x))^2)
  if (tot == 0) return(NA_real_)
  sum((x - ave(x, g))^2) / tot
}

diagnostics <- function(idx, dims = DFCA_DIMS) {
  dm <- idx[, names(dims), drop = FALSE]
  out <- list(
    coverage = colSums(!is.na(dm)),
    complete_country_years = sum(complete.cases(dm)),
    within_country_share = round(c(vapply(dm, within_share, numeric(1), g = idx$iso),
                                   dfca = within_share(idx$dfca, idx$iso)), 3)
  )
  if (ncol(dm) >= 2) {
    out$dimension_correlations <- round(cor(dm, use = "pairwise.complete.obs"), 3)
    pc <- prcomp(na.omit(dm), scale. = TRUE)
    out$pca_variance_share <- round(summary(pc)$importance[2, ], 3)
    out$first_component_loadings <- round(pc$rotation[, 1], 3)
  }
  out
}
