# Stratifies the CSD result by winter for both networks, converting each dataset's time stamps to local time.
# Inputs: 02_MART/intl_ews.rds, icetarget_ews.rds.
# Output: 02_MART/ews_winter_strat.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")

auc <- function(score, pos) {
  o <- is.finite(score); score <- score[o]; pos <- pos[o]
  if (length(unique(pos)) != 2L) return(NA_real_)
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

add_winter <- function(ws, tzi) {
  setDT(ws)
  if (!"t_end" %in% names(ws)) stop("Column t_end missing — rerun scripts 09/50")
  ws[, tm := as.POSIXct(t_end * 60, origin = "1970-01-01", tz = tzi)]
  ws[, `:=`(yr = as.integer(format(tm, "%Y")), mo = as.integer(format(tm, "%m")))]
  ws[, winter := ifelse(mo >= 9L, yr, yr - 1L)]
  ws[]
}

strat <- function(ws, lab) {
  cat(sprintf("\n=== [%s] by winter ===\n", lab))
  s <- ws[, {
    on <- type == "onset"
    p <- if (sum(on & is.finite(var_tau)) >= 5L)
           suppressWarnings(wilcox.test(var_tau[on & is.finite(var_tau)],
                                        mu = 0, alternative = "greater")$p.value)
         else NA_real_
    .(n_stations = uniqueN(stn), onset = sum(on), control = sum(!on),
      var_median = round(median(var_tau[on], na.rm = TRUE), 3),
      var_AUC = round(auc(var_tau, on), 3),
      ar1_AUC = round(auc(ar1_tau, on), 3),
      `τ>0_p` = signif(p, 3))
  }, by = winter][order(winter)]
  print(s)
  s[, set := lab][]
}

A <- add_winter(readRDS(file.path(MART, "intl_ews.rds")),      "UTC")
B <- add_winter(readRDS(file.path(MART, "icetarget_ews.rds")), "Asia/Seoul")
cat("Primary 121 stations, windows", nrow(A), "· independent 127 stations, windows", nrow(B), "\n")

RA <- strat(A, "primary 121 stations")
RB <- strat(B, "independent 127 stations")

cat("\n=== Pooled (all winters) ===\n")
for (x in list(list(A, "primary 121 stations"), list(B, "independent 127 stations"))) {
  w <- x[[1]]; on <- w$type == "onset"
  cat(sprintf("  %-16s var_AUC %.3f · ar1_AUC %.3f\n", x[[2]],
              auc(w$var_tau, on), auc(w$ar1_tau, on)))
}

saveRDS(rbindlist(list(RA, RB), fill = TRUE), file.path(MART, "ews_winter_strat.rds"))
cat("\nSaved: 02_MART/ews_winter_strat.rds\n")
