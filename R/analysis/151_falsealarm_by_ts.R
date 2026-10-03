# Compares false alarms per station-day by surface-temperature bin between the RWIS-only and surrounding-field models at the same operating threshold.
# Inputs: 02_MART/intl_ew_features.rds, intl_ew_pattern_store.rds, trigger_unified.csv.
# Output: 02_MART/falsealarm_by_ts.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
feats <- c("cold_dose", "dep_dose", "ft", "ts", "rh", "dts_1h")
D <- readRDS(file.path(MART, "intl_ew_features.rds"))
D <- D[complete.cases(D[, ..feats])]; setorder(D, stn, tm_min)
S41 <- readRDS(file.path(MART, "intl_ew_pattern_store.rds"))
TR <- fread(file.path(MART, "trigger_unified.csv"), encoding = "UTF-8")[group == "Fig. 6b (with and without surrounding fields)"]
D[, ts_bin := cut(ts, c(-Inf, -2, 0, 2, Inf), labels = c("Ts<=-2", "-2<Ts<=0", "0<Ts<=2", "Ts>2"), right = TRUE)]
OUT <- list()
for (H in c(180L, 360L)) {
  st <- S41[[as.character(H)]]; stopifnot(length(st$b$pr) == nrow(D), identical(st$b$y, st$p$y))
  nonev <- is.na(D$tno) | D$tno > H
  for (m in c("b", "p")) {
    cfg <- sprintf("T+%d %s", H / 60, if (m == "b") "RWIS" else "RWIS + surrounding fields")
    thr <- TR[configuration == cfg]$probability_threshold; stopifnot(length(thr) == 1L)
    fire <- st[[m]]$pr >= thr
    far_all <- mean(fire[nonev]) * 48
    stopifnot(abs(far_all - TR[configuration == cfg]$false_alarms_per_station_day) < 0.006)
    OUT[[length(OUT) + 1L]] <- data.table(H_h = H / 60, model = if (m == "b") "RWIS" else "RWIS + surrounding fields",
      D[nonev, .(n_evaluation_points = .N), by = ts_bin][, false_alarms_per_station_day := sapply(ts_bin, function(b)
        mean(fire[nonev & D$ts_bin == b]) * 48)])
    OUT[[length(OUT)]] <- rbind(OUT[[length(OUT)]], data.table(H_h = H / 60, model = OUT[[length(OUT)]]$model[1],
      ts_bin = "full", n_evaluation_points = sum(nonev), false_alarms_per_station_day = far_all))
  }
}
O <- rbindlist(OUT)
O[, ts_bin := factor(ts_bin, c("Ts<=-2", "-2<Ts<=0", "0<Ts<=2", "Ts>2", "full"))]
W <- dcast(O, H_h + ts_bin + n_evaluation_points ~ model, value.var = "false_alarms_per_station_day")
W[, reduction_pct := round(100 * (RWIS - `RWIS + surrounding fields`) / RWIS, 1)]
setorder(W, H_h, ts_bin); print(W)
fwrite(W, file.path(MART, "falsealarm_by_ts.csv"), bom = TRUE)
cat("Saved: 02_MART/falsealarm_by_ts.csv\n")
