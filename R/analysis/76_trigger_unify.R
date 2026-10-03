# Computes every operational trigger metric reported in the paper with one threshold rule (_trigger_rule.R).
# Inputs: 02_MART/intl_ew_features.rds, nowcast_oof_store.rds, intl_ew_pattern_store.rds.
# Output: 02_MART/trigger_unified.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table); library(lightgbm)})
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
source(code_file("_trigger_rule.R"))

feats <- c("cold_dose", "dep_dose", "ft", "ts", "rh", "dts_1h")
step <- 30L

D <- readRDS(file.path(MART, "intl_ew_features.rds"))
D <- D[complete.cases(D[, ..feats])]; setorder(D, stn, tm_min)
cat(sprintf("Evaluation points %s · stations %d\n", format(nrow(D), big.mark = ","), uniqueN(D$stn)))
TG <- make_trigger(data.table(stn = D$stn, tm = D$tm_min, tno = D$tno), step)

S72 <- readRDS(file.path(MART, "nowcast_oof_store.rds"))
S41 <- readRDS(file.path(MART, "intl_ew_pattern_store.rds"))

auc <- function(s, p) { o <- is.finite(s); s <- s[o]; p <- p[o]
  np <- as.numeric(sum(p)); nn <- as.numeric(sum(!p))
  if (!np || !nn) return(NA_real_); r <- rank(s); (sum(r[p]) - np * (np + 1) / 2) / (np * nn) }

rows <- list()
report <- function(label, tbl, pr, H, k, target) {
  th <- TG$thr_for_pod(pr, H, k, target)
  if (is.na(th)) { cat(sprintf("  %-28s target POD %.2f not reachable\n", label, target)); return(invisible()) }
  r <- TG$trig(pr, H, k, th)
  cat(sprintf("  %-28s POD %.3f · false alarms %6.2f per station-day · median lead %3.0f min (IQR %.0f-%.0f) · onsets %d\n",
              label, r$pod, r$far, median(r$lead), quantile(r$lead, .25), quantile(r$lead, .75), r$n_onset))
  rows[[length(rows) + 1L]] <<- data.table(
    group = tbl, configuration = label, lead_h = H / 60, persistence_k = k, target_pod = target,
    achieved_pod = round(r$pod, 3), false_alarms_per_station_day = round(r$far, 2),
    first_trigger_lead_median_min = round(median(r$lead)), first_trigger_lead_q1_min = round(quantile(r$lead, .25)),
    first_trigger_lead_q3_min = round(quantile(r$lead, .75)), onset = r$n_onset, probability_threshold = th)
}

cat("\n=== [0] RWIS-only OOF agreement between the two stores (checked on every run) ===\n")
for (H in c("180", "360")) {
  a <- S72[[H]]$pr; b <- S41[[H]]$b$pr; y <- S72[[H]]$y
  stopifnot(identical(y, S41[[H]]$b$y))
  cat(sprintf("  T+%-2s h  AUC %.4f vs %.4f | rank correlation %.5f · max probability difference %.2e\n",
              as.integer(H) / 60, auc(a, y == 1L), auc(b, y == 1L),
              cor(a, b, method = "spearman"), max(abs(a - b))))
}

cat("\n=== [1] Table 4 — RWIS-only operational trigger ===\n")
report("T+6 POD 0.80",          "Table 4", S72[["360"]]$pr, 360L, 1L, 0.80)
report("T+6 POD 0.70",          "Table 4", S72[["360"]]$pr, 360L, 1L, 0.70)
report("T+3 POD 0.80",          "Table 4", S72[["180"]]$pr, 180L, 1L, 0.80)
report("T+6 persistence k=2 POD 0.55", "Table 4", S72[["360"]]$pr, 360L, 2L, 0.55)

cat("\n=== [2] Fig. 6b — RWIS vs RWIS + surrounding fields ===\n")
report("T+6 RWIS",           "Fig. 6b (with and without surrounding fields)", S41[["360"]]$b$pr, 360L, 1L, 0.80)
report("T+6 RWIS + surrounding fields",  "Fig. 6b (with and without surrounding fields)", S41[["360"]]$p$pr, 360L, 1L, 0.80)
report("T+3 RWIS",           "Fig. 6b (with and without surrounding fields)", S41[["180"]]$b$pr, 180L, 1L, 0.80)
report("T+3 RWIS + surrounding fields",  "Fig. 6b (with and without surrounding fields)", S41[["180"]]$p$pr, 180L, 1L, 0.80)

cat("\n=== [3] Model ladder (H=3h, POD 0.80) ===\n")
H3 <- 180L
y3 <- as.integer(!is.na(D$tno) & D$tno <= H3)
stns <- unique(D$stn)
fold <- setNames(sample(rep(1:5, length.out = length(stns))), stns)
D[, fold := fold[as.character(stn)]]
params <- list(objective = "binary", learning_rate = 0.05, num_leaves = 31,
               min_data_in_leaf = 100, feature_fraction = 0.8,
               bagging_fraction = 0.8, bagging_freq = 1, verbosity = -1)
oof <- function(fs) {
  pr <- numeric(nrow(D))
  for (kf in 1:5) {
    tr <- which(D$fold != kf); te <- which(D$fold == kf)
    m <- lgb.train(params, lgb.Dataset(as.matrix(D[tr, ..fs]), label = y3[tr]), nrounds = 300)
    p <- predict(m, as.matrix(D[te, ..fs])); if (is.matrix(p)) p <- p[, 1]; pr[te] <- p
  }
  pr
}
sc0 <- -D$ts
cat(sprintf("  [Check] Model 0 AUC %.3f\n", auc(sc0, y3 == 1L)))
sc1 <- oof(c("ts", "rh"));  cat(sprintf("  [Check] Model 1 AUC %.3f\n", auc(sc1, y3 == 1L)))
sc2 <- S72[["180"]]$pr;     cat(sprintf("  [Check] Model 2 AUC %.3f (reusing the script-72 store)\n", auc(sc2, y3 == 1L)))
report("Model 0 temperature threshold",       "Fig. 5 (information ladder)", sc0, H3, 1L, 0.80)
report("Model 1 current state",       "Fig. 5 (information ladder)", sc1, H3, 1L, 0.80)
report("Model 2 (current + recent + cumulative)", "Fig. 5 (information ladder)", sc2, H3, 1L, 0.80)

out <- rbindlist(rows)
fwrite(out, file.path(MART, "trigger_unified.csv"))
cat("\n=== Summary table ===\n"); print(out[, .(group, configuration, achieved_pod, false_alarms_per_station_day, first_trigger_lead_median_min, onset)])
cat("\nSaved: 02_MART/trigger_unified.csv\n")
