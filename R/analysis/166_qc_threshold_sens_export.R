# Exports the sensitivity of the independent-network CSD result to the |Ts-Ta| QC tolerance (none, 30, 15, 10 degC) and counts onsets without tau.
# Inputs: 02_MART/icetarget_ews_thr_sens.rds, icetarget_ews.rds.
# Output: 02_MART/qc_threshold_sens.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
auc <- function(score, pos) {
  o <- is.finite(score); score <- score[o]; pos <- pos[o]
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
S <- readRDS(file.path(MART, "icetarget_ews_thr_sens.rds"))
thr <- c(A = NA, B = 30, C = 15, D = 10)
O <- rbindlist(lapply(names(S), function(k) { w <- as.data.table(S[[k]])
  data.table(variant = k, tolerance_c = thr[[k]], onset = sum(w$type == "onset"),
             n_onset_tau_computed = sum(w$type == "onset" & is.finite(w$var_tau)),
             var_AUC = auc(w$var_tau, w$type == "onset"), ar1_AUC = auc(w$ar1_tau, w$type == "onset")) }))
print(O)
M <- as.data.table(readRDS(file.path(MART, "icetarget_ews.rds")))
stopifnot(nrow(M) == nrow(S$A), sum(M$type == "onset" & is.finite(M$var_tau)) == O[variant == "A", n_onset_tau_computed])
fwrite(O, file.path(MART, "qc_threshold_sens.csv"), bom = TRUE)
cat("Saved: 02_MART/qc_threshold_sens.csv\n")
