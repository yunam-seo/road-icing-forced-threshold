# Applies the CSD procedure of script 09 to synthetic linear stochastic series whose recovery rate declines towards the end of the window (fold-type scaling), including a power curve and block gaps.
# Input: settings from 09_intl_ews.R.
# Output: 02_MART/ews_positive_control.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
sys.source(code_file("09_intl_ews.R"), envir = environment(), toplevel.env = environment())
set.seed(2026)

W <- cfg$W
LAM0 <- 0.05

sim_series <- function(lam1, miss = 0) {
  tt <- seq(0, 1, length.out = W + 1L)
  lam <- lam1 + (LAM0 - lam1) * sqrt(pmax(0, 1 - tt))
  x <- numeric(W + 1L)
  x[1] <- rnorm(1, 0, sqrt(1 / (2 * LAM0)))
  e <- rnorm(W)
  for (i in 2:(W + 1L)) x[i] <- x[i - 1L] * (1 - lam[i]) + e[i - 1L]
  if (miss > 0) {
    nb <- floor((W + 1L) / 30L)
    drop <- sample.int(nb, max(1L, round(nb * miss)))
    for (b in drop) x[((b - 1L) * 30L + 1L):min(b * 30L, W + 1L)] <- NA_real_
  }
  x
}

auc <- function(s, p) {
  o <- is.finite(s); s <- s[o]; p <- p[o]
  np <- as.numeric(sum(p)); nn <- as.numeric(sum(!p))
  if (!np || !nn) return(NA_real_)
  r <- rank(s); (sum(r[p]) - np * (np + 1) / 2) / (np * nn)
}

run_setting <- function(lam1, n, miss = 0) {
  on <- t(vapply(seq_len(n),      function(i) ews_window(sim_series(lam1, miss)), numeric(4)))
  ct <- t(vapply(seq_len(2L * n), function(i) ews_window(sim_series(LAM0, miss)), numeric(4)))
  d <- rbind(data.table(type = "onset", var_tau = on[, 1], ar1_tau = on[, 2]),
             data.table(type = "control", var_tau = ct[, 1], ar1_tau = ct[, 2]))
  p <- d$type == "onset"
  list(auc_var = auc(d$var_tau, p), auc_ar1 = auc(d$ar1_tau, p),
       pos_var = mean(d$var_tau[p] > 0, na.rm = TRUE),
       pos_ar1 = mean(d$ar1_tau[p] > 0, na.rm = TRUE),
       n_valid = sum(is.finite(d$var_tau[p])))
}

log_path <- file.path(ROOT, "05_REPORT/LOG/71_ews_poscontrol.log")
con <- file(log_path, open = "wt", encoding = "UTF-8"); sink(con, split = TRUE)

cat("=== Positive control: script-09 procedure applied unchanged to a linear stochastic process with a declining recovery rate ===\n")
cat(sprintf("Settings (script 09): W=%d rw=%d step=%d bw=%g cov_thr=%.2f · 2 controls per event\n",
            cfg$W, cfg$rw, cfg$step, cfg$bw, cfg$cov_thr))
cat(sprintf("Control recovery rate λ0=%.3f (lag-1 AR1 %.3f)\n\n", LAM0, 1 - LAM0))

N_MAIN <- 2427L
cat("[1] Same size as the observations (onset", N_MAIN, "· control", 2L * N_MAIN, ")\n")
m <- run_setting(lam1 = 0.005, n = N_MAIN)
cat(sprintf("  λ1=0.005 (window-end AR1 0.995): discrimination AUC rolling variance %.3f · autocorrelation %.3f | τ>0 share %.3f / %.3f\n\n",
            m$auc_var, m$auc_ar1, m$pos_var, m$pos_ar1))

cat("[2] Power curve (onset 600 · control 1200)\n")
res <- list(); k <- 0L
for (l1 in c(0.005, 0.010, 0.020, 0.030, 0.040, 0.050)) {
  r <- run_setting(lam1 = l1, n = 600L); k <- k + 1L
  res[[k]] <- data.table(setting = "no gaps", lam1 = l1, window_end_ar1 = round(1 - l1, 3),
                         ar1_rise = round((1 - l1) - (1 - LAM0), 3),
                         auc_rolling_variance = round(r$auc_var, 3), auc_lag1_autocorrelation = round(r$auc_ar1, 3),
                         frac_positive_tau_variance = round(r$pos_var, 3), n_valid_windows = r$n_valid)
  cat(sprintf("  λ1=%.3f (AR1 %.3f, rise %+.3f): AUC variance %.3f · autocorrelation %.3f\n",
              l1, 1 - l1, (1 - l1) - (1 - LAM0), r$auc_var, r$auc_ar1))
}
cat("\n[3] Robustness to block gaps (onset 600 · control 1200)\n")
cat("    Note: at 20 % missing, window coverage is 0.792, below the script-09 threshold (0.8), so all windows are excluded.\n")
cat("        This reflects the strict no-interpolation rule; results are shown for 5 % and 10 %, which pass the threshold.\n")
for (lm in list(c(0.005, 0.05), c(0.005, 0.10), c(0.020, 0.10))) {
  l1 <- lm[1]; ms <- lm[2]
  r <- run_setting(lam1 = l1, n = 600L, miss = ms); k <- k + 1L
  res[[k]] <- data.table(setting = sprintf("block gaps %.0f%%", ms * 100), lam1 = l1,
                         window_end_ar1 = round(1 - l1, 3),
                         ar1_rise = round((1 - l1) - (1 - LAM0), 3),
                         auc_rolling_variance = round(r$auc_var, 3), auc_lag1_autocorrelation = round(r$auc_ar1, 3),
                         frac_positive_tau_variance = round(r$pos_var, 3), n_valid_windows = r$n_valid)
  cat(sprintf("  λ1=%.3f · missing %.0f %%: AUC variance %.3f · autocorrelation %.3f (valid windows %d/600)\n",
              l1, ms * 100, r$auc_var, r$auc_ar1, r$n_valid))
}

S <- rbindlist(res)
fwrite(S, file.path(MART, "ews_positive_control.csv"), bom = TRUE)
cat("\n=== Summary table ===\n"); print(S)
cat(sprintf("[Power at the size of this experiment] AUC %.3f / %.3f at λ1=0.005\n", m$auc_var, m$auc_ar1))
cat("Saved: 02_MART/ews_positive_control.csv\n")
sink(); close(con)
cat("Log:", log_path, "\n")
