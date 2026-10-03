# Station-cluster bootstrap (B = 5,000) 95% confidence intervals of the CSD discrimination AUC (Table 2, Fig. 4d).
# Inputs: 02_MART/intl_ews.rds, 02_MART/icetarget_ews.rds.
# Output: 02_MART/table2_auc_ci.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
E <- readRDS(file.path(MART, "intl_ews.rds")); setDT(E)
I <- readRDS(file.path(MART, "icetarget_ews.rds")); setDT(I)
E[, stn := as.integer(sub("_.*", "", eid))]

auc <- function(a, b) {
  r <- rank(c(a, b))
  (sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2) / (length(a) * length(b))
}

boot_ci <- function(D, nm, B, seed = 2026) {
  set.seed(seed)
  s <- unique(D$stn)
  v <- replicate(B, {
    ss <- sample(s, length(s), replace = TRUE)
    d <- D[stn %in% ss]
    o <- d[type == "onset"][[nm]]; c_ <- d[type == "control"][[nm]]
    o <- o[is.finite(o)]; c_ <- c_[is.finite(c_)]
    if (!length(o) || !length(c_)) NA_real_ else auc(o, c_)
  })
  quantile(v, c(0.025, 0.975), na.rm = TRUE)
}

SETS <- list(list("primary", E), list("independent", I))
IND <- list(list("rolling variance tau", "var_tau"), list("lag-1 autocorrelation tau", "ar1_tau"))
cat("== Station-cluster bootstrap, B = 5,000 ==\n")
rows <- list()
for (S in SETS) for (K in IND) {
  D <- S[[2]]; nm <- K[[2]]
  o <- D[type == "onset"][[nm]]; c_ <- D[type == "control"][[nm]]
  o <- o[is.finite(o)]; c_ <- c_[is.finite(c_)]
  a <- auc(o, c_); ci <- boot_ci(D, nm, 5000L)
  cat(sprintf("  %-10s %-16s AUC %.3f  CI %.3f–%.3f\n", S[[1]], K[[1]], a, ci[1], ci[2]))
  rows[[length(rows) + 1L]] <- data.table(
    network = S[[1]], indicator = K[[1]], n_onset = length(o), n_control = length(c_),
    AUC = round(a, 3), ci_low = round(ci[1], 3), ci_high = round(ci[2], 3))
}
out <- rbindlist(rows)
fwrite(out, file.path(MART, "table2_auc_ci.csv"), bom = TRUE)
cat("\nSaved: 02_MART/table2_auc_ci.csv\n")
cat("All CIs within 0.45–0.55: ")
cat(if (all(out$ci_low >= 0.45 & out$ci_high <= 0.55)) "yes\n" else "no\n")
