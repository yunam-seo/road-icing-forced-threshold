# Computes values quoted in the text: station metadata summary, Fig. 2 statistics, clustering and trajectory-classification statistics.
# Inputs: station metadata of the independent network (RICE_VALID_META); meta/stn_exclusion_intl.csv; 02_MART/intl_event_table.rds, intl_lagged_table.rds, intl_physical_diag.csv, intl_continuum.csv.
# Outputs: 02_MART/station_meta_summary.csv, fig2_stats.csv, cluster_stats.csv, trajectory_classification.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table); library(cluster); library(diptest); library(lightgbm)})
MART <- file.path(ROOT, "02_MART")

DDIR <- rwis_dir()
EX <- fread(meta_file("stn_exclusion_intl.csv"))
sid <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", list.files(DDIR, pattern = "[.]csv$")))
sid <- setdiff(sid, EX$STN_ID)
M <- fread(valid_meta(), encoding = "UTF-8")
M <- M[STN_ID %in% sid]; stopifnot(nrow(M) == length(sid))
# KMA station metadata column names (Korean): elevation, expressway name, station type, sensor type, year built, latitude, longitude
setnames(M, c("\uace0\ub3c4", "\uace0\uc18d\ub3c4\ub85c\uba85", "\uad00\uce21\uc18c\uc885\ub958", "\uc13c\uc11c\uc885\ub958", "\uad6c\ucd95\uc5f0\ub3c4", "\uc704\ub3c4", "\uacbd\ub3c4"),
         c("elev_txt", "route", "stn_type", "sensor_type", "built_year", "lat", "lon"))
M[, elev := as.numeric(sub("m$", "", elev_txt))]
SM <- rbindlist(list(
  data.table(item = "analysed stations", value = nrow(M)),
  data.table(item = "number of expressway routes", value = uniqueN(M$route)),
  M[, .(value = .N), by = stn_type][, .(item = paste0("station type ", stn_type), value)],
  M[, .(value = .N), by = sensor_type][, .(item = paste0("sensor type ", sensor_type), value)],
  M[, .(value = .N), by = built_year][order(built_year), .(item = paste0("construction year ", built_year), value)],
  data.table(item = c("elevation min", "elevation max", "elevation median", "latitude min", "latitude max", "longitude min", "longitude max"),
             value = c(min(M$elev), max(M$elev), median(M$elev), min(M$lat), max(M$lat), min(M$lon), max(M$lon)))))
print(SM); fwrite(SM, file.path(MART, "station_meta_summary.csv"), bom = TRUE)

ev <- readRDS(file.path(MART, "intl_event_table.rds")); ev[, eid := paste0(stn, "_", onset_min)]
ph <- fread(file.path(MART, "intl_physical_diag.csv"))
D <- merge(ph, ev[, .(eid, pre_state, ts0, rh0, m_liq)], by = "eid")
eta2 <- function(yv, g) { g <- factor(g); gm <- tapply(yv, g, mean, na.rm = TRUE)
  ss_b <- sum(table(g) * (gm - mean(yv, na.rm = TRUE))^2, na.rm = TRUE)
  ss_t <- sum((yv - mean(yv, na.rm = TRUE))^2, na.rm = TRUE); ss_b / ss_t }
Pc <- fread(file.path(MART, "intl_continuum.csv"))
F2 <- data.table(item = c("eta2 ts0", "eta2 rh0", "eta2 m_liq", "eta2 ts_amp", "dip p PC1"),
                 value = c(sapply(c("ts0", "rh0", "m_liq", "ts_amp"), function(v) eta2(D[[v]], D$pre_state)), dip.test(Pc$PC1)$p.value))
print(F2); fwrite(F2, file.path(MART, "fig2_stats.csv"), bom = TRUE)

set.seed(2026)
lg <- readRDS(file.path(MART, "intl_lagged_table.rds")); setorder(lg, eid, lag)
fill_event <- function(v, x) { if (all(is.na(v))) return(v); approx(x[!is.na(v)], v[!is.na(v)], xout = x, rule = 2)$y }
cov <- lg[, .(cov = mean(!is.na(ts))), by = eid]; keep <- cov[cov >= 0.7, eid]
d <- lg[eid %in% keep]; for (v in c("ts", "m_dep", "rh")) d[, (v) := fill_event(get(v), lag), by = eid]
w <- dcast(melt(d, id.vars = c("eid", "lag"), measure.vars = c("ts", "m_dep", "rh")), eid ~ variable + lag, value.var = "value")
feat <- as.matrix(w[, -1]); rownames(feat) <- w$eid; feat <- feat[complete.cases(feat), , drop = FALSE]; Z <- scale(feat)
sil <- sapply(2:6, function(k) { km <- kmeans(Z, centers = k, nstart = 25, iter.max = 50); mean(silhouette(km$cluster, dist(Z))[, 3]) })
kbest <- (2:6)[which.max(sil)]
km <- kmeans(Z, centers = kbest, nstart = 50, iter.max = 100)
asg <- merge(data.table(eid = rownames(Z), cluster = km$cluster), ev[, .(eid, pre_state)], by = "eid")
ari <- function(a, b) { t <- table(a, b); n <- sum(t); a_ <- rowSums(t); b_ <- colSums(t)
  idx <- sum(choose(t, 2)); e <- sum(choose(a_, 2)) * sum(choose(b_, 2)) / choose(n, 2)
  m <- (sum(choose(a_, 2)) + sum(choose(b_, 2))) / 2; (idx - e) / (m - e) }
CL <- data.table(item = c(paste0("silhouette k", 2:6), "best k", "ARI cluster vs prior state", "clustering input events", "clustering input features",
                          "nstart for k selection", "nstart final"),
                 value = c(sil, kbest, ari(asg$cluster, asg$pre_state), nrow(Z), ncol(Z), 25, 50))
print(CL); fwrite(CL, file.path(MART, "cluster_stats.csv"), bom = TRUE)

set.seed(2026)
w2 <- dcast(melt(lg, id.vars = c("eid", "lag"), measure.vars = c("ts", "m_dep", "rh")), eid ~ variable + lag, value.var = "value")
w2 <- merge(w2, ev[, .(eid, pre_state)], by = "eid"); w2 <- w2[complete.cases(w2)]
w2[, y := as.integer(factor(pre_state, levels = c("DRY", "WET", "SNOW"))) - 1L]
Xc <- as.matrix(w2[, setdiff(names(w2), c("eid", "pre_state", "y")), with = FALSE])
tr <- sample(nrow(w2), round(.7 * nrow(w2)))
m <- lgb.train(list(objective = "multiclass", num_class = 3, learning_rate = 0.05, num_leaves = 15, min_data_in_leaf = 20, verbosity = -1),
               lgb.Dataset(Xc[tr, ], label = w2$y[tr]), nrounds = 200)
pp <- predict(m, Xc[-tr, ]); if (is.null(dim(pp))) pp <- matrix(pp, ncol = 3, byrow = TRUE)
pred <- max.col(pp) - 1L; yte <- w2$y[-tr]
TC <- data.table(item = c("accuracy", "majority-class baseline", "balanced accuracy", "events (complete cases)", "test fraction", "number of features"),
                 value = c(mean(pred == yte), max(table(yte)) / length(yte),
                        mean(sapply(0:2, function(k) mean(pred[yte == k] == k))), nrow(w2), 0.3, ncol(Xc)))
print(TC); fwrite(TC, file.path(MART, "trajectory_classification.csv"), bom = TRUE)

cat("Saved: station_meta_summary.csv · fig2_stats.csv · cluster_stats.csv · trajectory_classification.csv\n")
