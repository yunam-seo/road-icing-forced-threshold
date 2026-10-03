# Runs one-sided Wilcoxon signed-rank tests of onset Kendall tau > 0 for both networks.
# Inputs: 02_MART/intl_ews.rds, icetarget_ews.rds.
# Output: 02_MART/csd_wilcoxon_onesided.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
A <- as.data.table(readRDS(file.path(MART, "intl_ews.rds")))
B <- as.data.table(readRDS(file.path(MART, "icetarget_ews.rds")))
cat("Primary network type:", paste(names(table(A$type)), table(A$type), collapse = " "), "\n")
cat("Independent network type:", paste(names(table(B$type)), table(B$type), collapse = " "), "\n")
one <- function(d, net) rbindlist(lapply(c("var_tau", "ar1_tau"), function(v) {
  x <- d[type == "onset"][[v]]; x <- x[is.finite(x)]
  w <- suppressWarnings(wilcox.test(x, mu = 0, alternative = "greater"))
  data.table(network = net, indicator = v, n = length(x), V = unname(w$statistic), p = w$p.value)
}))
O <- rbind(one(A, "primary"), one(B, "independent"))
print(O)
fwrite(O, file.path(MART, "csd_wilcoxon_onesided.csv"), bom = TRUE)
cat("Saved: 02_MART/csd_wilcoxon_onesided.csv\n")
