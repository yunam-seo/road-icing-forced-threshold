# Recomputes the trigger rows of the nowcast bootstrap table with the shared threshold rule in _trigger_rule.R; AUC rows are unchanged.
# Inputs: 02_MART/nowcast_oof_store.rds, intl_ew_features.rds, nowcast_bootstrap_ci.csv, trigger_unified.csv.
# Output: 02_MART/nowcast_bootstrap_ci.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
source(code_file("_trigger_rule.R"))
B <- 1000L; STEP <- 30L; ppd <- 24 * 60 / STEP

feats <- c("cold_dose", "dep_dose", "ft", "ts", "rh", "dts_1h")
D <- readRDS(file.path(MART, "intl_ew_features.rds"))
D <- D[complete.cases(D[, ..feats])]; setorder(D, stn, tm_min)
store <- readRDS(file.path(MART, "nowcast_oof_store.rds"))
stns <- unique(D$stn)
cat(sprintf("Evaluation points %s · stations %d · replicates B=%d\n\n", format(nrow(D), big.mark = ","),
            length(stns), B))

EV <- D[, .(stn, tm = tm_min, tno)]
tg <- make_trigger(EV, STEP)

fire_of <- function(pr, H, k, thr) {
  dt <- copy(EV)[, p := pr]
  dt[, fire := {
    a <- p >= thr
    if (k == 1L) a else {
      ok <- a
      for (j in 1:(k - 1L)) {
        prev <- shift(a, j) & (tm - shift(tm, j) == j * STEP)
        ok <- ok & (prev %in% TRUE)
      }
      ok }
  }, by = stn]
  dt
}
per_station <- function(pr, H, k, thr) {
  dt <- fire_of(pr, H, k, thr)
  on <- dt[!is.na(tno) & tno > 0 & tno <= H]
  on[, oid := paste0(stn, "_", tm + tno)]
  ev <- on[, .(det = any(fire), lead = if (any(fire)) as.numeric(max(tno[fire])) else NA_real_),
           by = .(stn, oid)]
  st <- ev[, .(n_ons = .N, n_det = sum(det)), by = stn]
  fa <- dt[, .(n_fa = sum(fire & (is.na(tno) | tno > H)), n_pt = .N), by = stn]
  list(st = merge(st, fa, by = "stn", all = TRUE), lead = ev[!is.na(lead)])
}

draw <- replicate(B, sample(stns, length(stns), replace = TRUE), simplify = FALSE)
CFG <- list(list(360L, 1L, 0.80, "T+6 POD 0.80"), list(360L, 1L, 0.70, "T+6 POD 0.70"),
            list(180L, 1L, 0.80, "T+3 POD 0.80"), list(360L, 2L, 0.55, "T+6 persistence k=2"))
rows <- list()
for (a in CFG) {
  H <- a[[1]]; k <- a[[2]]; target <- a[[3]]; lab <- a[[4]]
  pr <- store[[as.character(H)]]$pr
  thr <- tg$thr_for_pod(pr, H, k, target)
  ref <- tg$trig(pr, H, k, thr)
  r <- per_station(pr, H, k, thr)
  st <- r$st
  for (cn in c("n_ons", "n_det", "n_fa", "n_pt")) set(st, which(is.na(st[[cn]])), cn, 0)
  agg <- function(ss) {
    x <- st[match(ss, st$stn)]
    c(sum(x$n_det) / sum(x$n_ons), sum(x$n_fa) / (sum(x$n_pt) / ppd))
  }
  bs <- vapply(draw, function(ss) agg(as.integer(ss)), numeric(2))
  lb <- vapply(draw, function(ss) {
    v <- r$lead[stn %in% as.integer(ss)]$lead
    if (!length(v)) NA_real_ else median(v)}, numeric(1))
  ci <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
  c1 <- ci(bs[1, ]); c2 <- ci(bs[2, ]); c3 <- ci(lb)
  cat(sprintf("  %-16s POD %.3f (%.3f-%.3f) · false alarms %.2f (%.2f-%.2f) · lead %.0f (%.0f-%.0f)\n",
              lab, ref$pod, c1[1], c1[2], ref$far, c2[1], c2[2],
              median(ref$lead), c3[1], c3[2]))
  rows <- c(rows, list(
    data.table(indicator = paste(lab, "POD"), value = round(ref$pod, 3),
               ci_low = round(c1[1], 3), ci_high = round(c1[2], 3)),
    data.table(indicator = paste(lab, "false alarms per station-day"), value = round(ref$far, 2),
               ci_low = round(c2[1], 2), ci_high = round(c2[2], 2)),
    data.table(indicator = paste(lab, "first-trigger lead (min)"), value = round(median(ref$lead)),
               ci_low = round(c3[1]), ci_high = round(c3[2]))))
}

f <- file.path(MART, "nowcast_bootstrap_ci.csv")
old <- fread(f)
keep <- old[grepl("AUC$", indicator)]
S <- rbind(keep, rbindlist(rows))
fwrite(S, f, bom = TRUE)

cat("\n== Comparison with trigger_unified.csv ==\n")
TU <- fread(file.path(MART, "trigger_unified.csv"))
MAP <- c("T+6 POD 0.80" = "T+6 POD 0.80", "T+6 POD 0.70" = "T+6 POD 0.70",
         "T+3 POD 0.80" = "T+3 POD 0.80", "T+6 persistence k=2" = "T+6 persistence k=2 POD 0.55")
bad <- 0L
for (lab in names(MAP)) {
  u <- TU[configuration == MAP[[lab]]]
  fa <- S[indicator == paste(lab, "false alarms per station-day")]$value
  ld <- S[indicator == paste(lab, "first-trigger lead (min)")]$value
  okf <- abs(fa - u$false_alarms_per_station_day) < 0.005; okl <- abs(ld - u$first_trigger_lead_median_min) < 0.5
  if (!okf || !okl) bad <- bad + 1L
  cat(sprintf("  %-16s false alarms %.2f vs %.2f %s · lead %.0f vs %.0f %s\n",
              lab, fa, u$false_alarms_per_station_day, if (okf) "OK" else "MISMATCH", ld, u$first_trigger_lead_median_min,
              if (okl) "OK" else "MISMATCH"))
}
cat(sprintf("\nMismatches: %d\n", bad))
cat("Saved: 02_MART/nowcast_bootstrap_ci.csv\n")
