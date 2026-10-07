## Two figures from the SIMULATED panel, to show what the output looks like.
## Run from this folder:  Rscript demo_figures.R
## Writes figure_rank_bands_simulated.png and figure_levels_simulated.png.

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(f)) setwd(dirname(normalizePath(f)))

source("build_index.R")
source("simulate_panel.R")

df0  <- simulate_panel()
df   <- impute_within_country(df0, max_gap = 2)   # gaps of up to two years, flagged in imp_*
sens <- sensitivity(df, reps_weights = 100)

## colours: categorical slots 1 and 2 of the reference palette, text in neutral ink
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
grp <- unique(df[, c("iso", "dominant")])
col_of <- function(iso) ifelse(grp$dominant[match(iso, grp$iso)] == 1, ORANGE, BLUE)

heading <- function(title, ..., top = 4.4) {
  sub <- c(...)
  mtext(title, side = 3, line = top, adj = 0, font = 2, cex = 1.15, col = INK)
  for (i in seq_along(sub))
    mtext(sub[i], side = 3, line = top - 1.25 * i, adj = 0, cex = 0.82, col = INK2)
}

## ---- Figure 1: rank uncertainty, 2015 ------------------------------------------
rb <- rank_bands(sens, 2015)
rb <- rb[order(rb$rk), ]
n <- nrow(rb)

png("figure_rank_bands_simulated.png", width = 1600, height = 1100, res = 200, bg = SURFACE)
par(mar = c(4.2, 4.6, 6.2, 1.5), family = "sans", col.axis = INK2, col.lab = INK2, fg = INK2)
plot(NA, xlim = c(1, n), ylim = c(n + 0.5, 0.5), xlab = "", ylab = "", axes = FALSE)
abline(v = seq(2, n, by = 2), col = GRID, lwd = 0.5)
segments(rb$rank_p05, seq_len(n), rb$rank_p95, seq_len(n), col = col_of(rb$iso), lwd = 1.3, lend = 1)
points(rb$rk, seq_len(n), pch = 21, cex = 1.5, bg = col_of(rb$iso), col = SURFACE, lwd = 2)
axis(1, at = seq(1, n, by = 2), tick = FALSE, line = -0.4, col.axis = INK2)
axis(2, at = seq_len(n), labels = rb$iso, las = 1, tick = FALSE, line = -0.4, col.axis = INK2, cex.axis = 0.85)
mtext("Rank, 1 = most autonomous", side = 1, line = 2.4, cex = 0.85, col = INK2)
heading("Rankings carry wide uncertainty bands",
        "SIMULATED DATA, not real countries. Dot = median rank, line = 5th to 95th percentile",
        sprintf("across %s versions of the index. %d countries have a 2015 value.",
                format(length(unique(sens$draw)), big.mark = ","), n))
legend("topright", legend = c("simulated high-capture group", "simulated low-capture group"),
       pch = 21, pt.bg = c(ORANGE, BLUE), col = SURFACE, pt.lwd = 2, pt.cex = 1.4,
       bty = "n", text.col = INK2, cex = 0.85)
invisible(dev.off())

## ---- Figure 2: levels over time with bands -------------------------------------
shade <- function(x, ylo, yhi, col) {
  ## one polygon per unbroken run of years
  runs <- cumsum(c(1, diff(x) != 1))
  for (r in unique(runs)) {
    k <- runs == r
    if (sum(k) < 2) next
    polygon(c(x[k], rev(x[k])), c(ylo[k], rev(yhi[k])), border = NA,
            col = adjustcolor(col, alpha.f = 0.12))
  }
}
trace <- function(x, y, col) {
  runs <- cumsum(c(1, diff(x) != 1))
  for (r in unique(runs)) {
    k <- runs == r
    lines(x[k], y[k], col = col, lwd = 1.3, lend = 1, lty = 1)
    if (sum(k) == 1) points(x[k], y[k], pch = 21, bg = col, col = SURFACE, cex = 1.1, lwd = 2)
  }
}
a <- level_bands(sens, "C01")   # high-capture group, with the injected outlier run
b <- level_bands(sens, "C15")   # low-capture group

png("figure_levels_simulated.png", width = 1600, height = 1000, res = 200, bg = SURFACE)
par(mar = c(4.2, 4.2, 7.0, 6.2), family = "sans", col.axis = INK2, col.lab = INK2, fg = INK2)
plot(NA, xlim = range(df$year), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)
abline(h = seq(0, 1, by = 0.25), col = GRID, lwd = 0.5)
rect(2003 - 0.5, 0, 2008 + 0.5, 1, col = adjustcolor("#52514e", alpha.f = 0.06), border = NA)
text(2005.5, 0.985, "injected outlier run", cex = 0.72, col = INK2, adj = c(0.5, 1))
shade(a$year, a$p05, a$p95, ORANGE); shade(b$year, b$p05, b$p95, BLUE)
trace(a$year, a$median, ORANGE);     trace(b$year, b$median, BLUE)
axis(1, at = seq(1990, 2020, by = 10), tick = FALSE, line = -0.4, col.axis = INK2)
axis(2, at = seq(0, 1, by = 0.25), las = 1, tick = FALSE, line = -0.4, col.axis = INK2)
mtext("DFCA, 1 = fully autonomous", side = 2, line = 2.6, cex = 0.85, col = INK2)
## direct labels at the right-hand end, in ink, beside a coloured key
ends <- function(d, lab, col) {
  y <- d$median[nrow(d)]
  segments(max(df$year) + 0.6, y, max(df$year) + 1.8, y, col = col, lwd = 1.3, xpd = NA)
  text(max(df$year) + 2.2, y, lab, adj = 0, cex = 0.8, col = INK, xpd = NA)
}
ends(a, "C01", ORANGE); ends(b, "C15", BLUE)
heading("An index level should carry its band",
        "SIMULATED DATA, not real countries. Line = median DFCA, band = 5th to 95th percentile",
        sprintf("across %s versions of the index. Gaps of up to two years are interpolated; longer gaps stay blank.",
                format(length(unique(sens$draw)), big.mark = ",")),
        top = 5.4)
legend(x = min(df$year), y = 1.17, xpd = NA, horiz = TRUE,
       legend = c("C01, simulated high-capture group", "C15, simulated low-capture group"),
       lty = 1, lwd = 1.3, col = c(ORANGE, BLUE), bty = "n", text.col = INK2, cex = 0.8)
invisible(dev.off())

cat("wrote figure_rank_bands_simulated.png and figure_levels_simulated.png\n")
