# Applies the nowcast trained on the primary network, without retraining or recalibration, to the independent 127-station network.
# Inputs: 02_MART/intl_ew_features.rds, icetarget_event_table.rds; 1-min records of the independent network (RICE_VALID_FILE).
# Outputs: 02_MART/nowcast_independent127.csv, nowcast_ind127_features.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table); library(lightgbm)})
MART <- file.path(ROOT, "02_MART")
source(code_file("_trigger_rule.R"))
ICE <- valid_file()
E7 <- new.env(); sys.source(code_file("07_intl_cumstress.R"), envir = E7)
E3 <- new.env(); sys.source(code_file("03_intl_events_cleaned.R"), envir = E3)
cfg7 <- E7$cfg; cfg3 <- E3$cfg
set.seed(2026)
feats <- c("cold_dose", "dep_dose", "ft", "ts", "rh", "dts_1h")
Hs <- c(180L, 360L); step_eval <- 30L; ppd <- 24 * 60 / step_eval

log_path <- file.path(ROOT, "05_REPORT/LOG/73_nowcast_ind127.log")
con <- file(log_path, open = "wt", encoding = "UTF-8"); sink(con, split = TRUE)
cat("=== Nowcast applied to the independent 127 stations (no retraining or recalibration) ===\n")
cat(sprintf("[Reused] tau=%d gap_max=%d · evaluation grid %d min · features %s\n\n",
            cfg7$tau, cfg7$gap_max, step_eval, paste(feats, collapse = ",")))

D <- readRDS(file.path(MART, "intl_ew_features.rds"))
D <- D[complete.cases(D[, ..feats])]
cat(sprintf("[1] Training data (121 stations): %s points · %d stations\n", format(nrow(D), big.mark = ","), uniqueN(D$stn)))
params <- list(objective = "binary", learning_rate = 0.05, num_leaves = 31,
               min_data_in_leaf = 100, feature_fraction = 0.8,
               bagging_fraction = 0.8, bagging_freq = 1, verbosity = -1)
models <- list()
for (H in Hs) {
  y <- as.integer(!is.na(D$tno) & D$tno <= H)
  models[[as.character(H)]] <- lgb.train(
    params, lgb.Dataset(as.matrix(D[, ..feats]), label = y), nrounds = 300)
  cat(sprintf("    H=%d min model trained (positive rate %.4f)\n", H, mean(y)))
}

cat("\n[2] Building the feature table for the independent 127 stations ...\n")
d <- readRDS(ICE); setDT(d)
d <- d[, .(STN_ID, tnum, year, mon, TR1, TA, HM, SFS)]
d <- d[mon %in% cfg3$freeze_mon]
source(code_file("00_analysis_window.R"))
d <- filter_winter(d, ym = list(y = d$year, m = d$mon))
d[, STN_ID := as.integer(STN_ID)]; d[, tm_min := as.integer(tnum %/% 60)]; d[, tnum := NULL]
setorder(d, STN_ID, tm_min); d <- unique(d, by = c("STN_ID", "tm_min"))
ev127 <- readRDS(file.path(MART, "icetarget_event_table.rds")); setDT(ev127)
ecol <- if ("stn" %in% names(ev127)) "stn" else "STN_ID"
ocol <- if ("onset_min" %in% names(ev127)) "onset_min" else names(ev127)[grep("onset", names(ev127))][1]
cat(sprintf("    Independent event table: %s events · %d stations (columns %s/%s)\n",
            format(nrow(ev127), big.mark = ","), uniqueN(ev127[[ecol]]), ecol, ocol))

vcache <- file.path(MART, "nowcast_ind127_features.rds")
if (file.exists(vcache)) { V <- readRDS(vcache); cat("    (using feature-table cache)\n") } else {
out <- list(); k <- 0L
for (sid in sort(unique(d$STN_ID))) {
  oo <- sort(unique(ev127[get(ecol) == sid][[ocol]]))
  if (!length(oo)) next
  s <- d[STN_ID == sid]
  s[(TR1 < cfg7$ts_range[1]) | (TR1 > cfg7$ts_range[2]), TR1 := NA_real_]
  s[(TA  < cfg7$ta_range[1]) | (TA  > cfg7$ta_range[2]), TA  := NA_real_]
  s[(HM  < cfg7$hm_range[1]) | (HM  > cfg7$hm_range[2]), HM  := NA_real_]
  p <- s[, .(tm_min, ts = TR1, ta = TA, rh = HM, cls = E7$map_sfs_class(as.integer(SFS)))]
  if (nrow(p) < 120L) next
  L <- E7$leaky(p, cfg7$tau, cfg7$gap_max)
  M <- merge(p[, .(tm_min, ts, ta, rh)], L[, .(tm_min, cold_dose, dep_dose, ft, cls)], by = "tm_min")
  setkey(M, tm_min)
  M[, ts_1h := M[data.table(tm_min = M$tm_min - 60L), on = "tm_min", roll = "nearest"]$ts]
  M[, dts_1h := ts - ts_1h]
  E <- M[cls %in% c("DRY", "WET", "SNOW") & (tm_min %% step_eval == 0L) & is.finite(cold_dose)]
  if (!nrow(E)) next
  j <- findInterval(E$tm_min, oo)
  nxt <- ifelse(j < length(oo), oo[pmin(j + 1L, length(oo))], NA_integer_)
  E[, tno := nxt - tm_min]; E[is.na(tno) | tno < 0, tno := NA_integer_]
  k <- k + 1L
  out[[k]] <- E[, .(stn = sid, tm_min, cold_dose, dep_dose, ft, ts, rh, dts_1h, tno)]
}
V <- rbindlist(out, fill = TRUE)
V <- V[complete.cases(V[, ..feats])]
saveRDS(V, vcache)
}
cat(sprintf("    Evaluation points %s · stations %d\n", format(nrow(V), big.mark = ","), uniqueN(V$stn)))

auc <- function(s, p) {
  o <- is.finite(s); s <- s[o]; p <- p[o]
  np <- as.numeric(sum(p)); nn <- as.numeric(sum(!p))
  if (!np || !nn) return(NA_real_)
  r <- rank(s); (sum(r[p]) - np * (np + 1) / 2) / (np * nn)
}
TGV <- make_trigger(data.table(stn = V$stn, tm = V$tm_min, tno = V$tno), step_eval)
cat("\n[3] External application (127 stations never used in training)\n")
rows <- list(); kk <- 0L
for (H in Hs) {
  pr <- predict(models[[as.character(H)]], as.matrix(V[, ..feats]))
  if (is.matrix(pr)) pr <- pr[, 1]
  y <- as.integer(!is.na(V$tno) & V$tno <= H)
  a <- auc(pr, y == 1L)
  dt <- data.table(stn = V$stn, tm = V$tm_min, tno = V$tno, p = pr, y = y)
  setorder(dt, stn, tm)
  dt[, onset_id := ifelse(!is.na(tno) & tno <= H, paste0(stn, "_", tm + tno), NA_character_)]
  ons <- unique(dt[!is.na(onset_id)]$onset_id)
  pick <- TGV$thr_for_pod(pr, H, 1L, 0.80)
  if (is.na(pick)) {
    pick <- min(dt$p)
    cat("    No threshold reaches POD 0.80 - reporting the operating point at the lowest threshold\n")
  }
  hit <- dt$p >= pick
  pod <- uniqueN(dt[hit & !is.na(onset_id)]$onset_id) / length(ons)
  fa <- mean(hit[is.na(dt$onset_id)]) * ppd
  lead <- dt[hit & !is.na(onset_id), .(l = max(tno)), by = onset_id]$l
  lead <- lead[is.finite(lead)]
  cat(sprintf("  T+%d h: AUC %.3f (positive rate %.4f · onsets %d) | at POD %.2f false alarms %.2f per station-day · median lead %.0f min\n",
              H / 60, a, mean(y), length(ons), pod, fa, median(lead)))
  kk <- kk + 1L
  rows[[kk]] <- data.table(lead = sprintf("T+%d h", H / 60), station_set = "independent 127 stations (external application)",
                           AUC = round(a, 3), positive_rate = round(mean(y), 4),
                           POD = round(pod, 2), false_alarms_per_station_day = round(fa, 2),
                           first_trigger_lead_median_min = round(median(lead)))
}
S <- rbindlist(rows)
fwrite(S, file.path(MART, "nowcast_independent127.csv"), bom = TRUE)
ref <- fread(file.path(MART, "nowcast_bootstrap_ci.csv"), encoding = "UTF-8")
g <- function(k) as.numeric(ref[indicator == k]$value)
cat(sprintf("\n[Reference] Station-held-out CV on the 121-station set: T+3 h AUC %.3f, T+6 h %.3f\n",
            g("T+3 AUC"), g("T+6 AUC")))
cat("=== Summary table ===\n"); print(S)
cat("Saved: 02_MART/nowcast_independent127.csv\n")
sink(); close(con)
cat("Log:", log_path, "\n")
