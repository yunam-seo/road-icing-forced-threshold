# Computes station-cluster bootstrap CIs of nowcast AUC by lead time and of trigger metrics.
# Input: 02_MART/intl_ew_features.rds.
# Outputs: 02_MART/nowcast_bootstrap_ci.csv, nowcast_oof_store.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table); library(lightgbm)})
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
feats <- c("cold_dose", "dep_dose", "ft", "ts", "rh", "dts_1h")
Hs <- c(60L, 120L, 180L, 360L, 720L); step <- 30L; ppd <- 24 * 60 / step
B <- 1000L

D <- readRDS(file.path(MART, "intl_ew_features.rds"))
D <- D[complete.cases(D[, ..feats])]; setorder(D, stn, tm_min)
stns <- unique(D$stn)
fold <- setNames(sample(rep(1:5, length.out = length(stns))), stns)
D[, fold := fold[as.character(stn)]]
params <- list(objective = "binary", learning_rate = 0.05, num_leaves = 31,
               min_data_in_leaf = 100, feature_fraction = 0.8,
               bagging_fraction = 0.8, bagging_freq = 1, verbosity = -1)

oof_prob <- function(H) {
  y <- as.integer(!is.na(D$tno) & D$tno <= H); pr <- numeric(nrow(D))
  for (kf in 1:5) {
    tr <- which(D$fold != kf); te <- which(D$fold == kf)
    m <- lgb.train(params, lgb.Dataset(as.matrix(D[tr, ..feats]), label = y[tr]), nrounds = 300)
    p <- predict(m, as.matrix(D[te, ..feats])); if (is.matrix(p)) p <- p[, 1]; pr[te] <- p
  }
  list(y = y, pr = pr)
}

auc_w <- function(ord_y, w) {
  neg <- cumsum(w * (1 - ord_y))
  negb <- neg - w * (1 - ord_y)
  P <- sum(w * ord_y); N <- sum(w) - P
  if (P <= 0 || N <= 0) return(NA_real_)
  sum(w * ord_y * negb) / (P * N)
}

log_path <- file.path(ROOT, "05_REPORT/LOG/72_nowcast_ci.log")
con <- file(log_path, open = "wt", encoding = "UTF-8"); sink(con, split = TRUE)
cat("=== Station-cluster bootstrap 95 % CI of nowcast performance ===\n")
cat(sprintf("Evaluation points %s · stations %d · replicates B=%d · folds and hyperparameters as in script 19\n\n",
            format(nrow(D), big.mark = ","), length(stns), B))

sid <- as.character(D$stn)
draw <- replicate(B, sample(stns, length(stns), replace = TRUE), simplify = FALSE)

cache <- file.path(MART, "nowcast_oof_store.rds")
store <- if (file.exists(cache)) readRDS(cache) else list()
rows <- list(); k <- 0L
cat("[A] Station-held-out CV AUC by lead time", if (length(store)) "(using OOF cache)" else "", "\n")
for (H in Hs) {
  o <- if (!is.null(store[[as.character(H)]])) store[[as.character(H)]] else oof_prob(H)
  store[[as.character(H)]] <- o
  ord <- order(o$pr)
  oy <- o$y[ord]; os <- sid[ord]
  pt <- auc_w(oy, rep(1, length(oy)))
  bs <- vapply(draw, function(ss) {
    tab <- table(ss); w <- as.numeric(tab[os]); w[is.na(w)] <- 0
    auc_w(oy, w)
  }, numeric(1))
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  k <- k + 1L
  rows[[k]] <- data.table(indicator = sprintf("T+%d AUC", H / 60), value = round(pt, 3),
                          ci_low = round(ci[1], 3), ci_high = round(ci[2], 3))
  cat(sprintf("  T+%-2d h: AUC %.3f (95%% CI %.3f–%.3f)  positive rate %.4f\n",
              H / 60, pt, ci[1], ci[2], mean(o$y)))
}

saveRDS(store, cache)

cat("\n[B] Operational trigger metrics (threshold fixed on the full data; only stations resampled)\n")
trig_station <- function(H, k_sust, pod_target) {
  y <- store[[as.character(H)]]$y; p <- store[[as.character(H)]]$pr
  dt <- data.table(stn = D$stn, tm = D$tm_min, tno = D$tno, p = p)
  setorder(dt, stn, tm)
  dt[, onset_id := ifelse(!is.na(tno), paste0(stn, "_", tm + tno), NA_character_)]
  ons_all <- unique(dt[!is.na(onset_id) & tno <= H]$onset_id)
  thrs <- quantile(dt$p, seq(0.90, 0.9995, length.out = 120), names = FALSE)
  pick <- NA_real_
  for (th in rev(thrs)) {
    dt[, hit := p >= th]
    if (k_sust > 1L) dt[, hit := hit & shift(hit, 1L, fill = FALSE), by = stn]
    det <- uniqueN(dt[hit == TRUE & !is.na(onset_id) & tno <= H]$onset_id)
    if (det / length(ons_all) >= pod_target) { pick <- th; break }
  }
  dt[, hit := p >= pick]
  if (k_sust > 1L) dt[, hit := hit & shift(hit, 1L, fill = FALSE), by = stn]
  ev <- dt[!is.na(onset_id) & tno <= H, .(det = any(hit), lead = if (any(hit)) as.numeric(max(tno[hit])) else NA_real_),
           by = .(stn, onset_id)]
  st <- ev[, .(n_ons = .N, n_det = sum(det)), by = stn]
  fa <- dt[, .(n_fa = sum(hit & (is.na(tno) | tno > H)), n_pt = .N), by = stn]
  list(st = merge(st, fa, by = "stn", all = TRUE), lead = ev[!is.na(lead)], thr = pick)
}
boot_trig <- function(H, k_sust, pod_target, label) {
  r <- trig_station(H, k_sust, pod_target)
  st <- r$st
  for (cn in c("n_ons", "n_det", "n_fa", "n_pt")) set(st, which(is.na(st[[cn]])), cn, 0)
  agg <- function(ss) {
    x <- st[match(ss, st$stn)]
    pod <- sum(x$n_det) / sum(x$n_ons)
    fa <- sum(x$n_fa) / (sum(x$n_pt) / ppd)
    c(pod, fa)
  }
  pt <- agg(as.character(st$stn))
  bs <- vapply(draw, function(ss) agg(as.character(ss)), numeric(2))
  ldv <- r$lead$lead
  lb <- vapply(draw, function(ss) {
    v <- r$lead[stn %in% as.integer(ss)]$lead
    if (!length(v)) NA_real_ else median(v)
  }, numeric(1))
  ci <- function(v) quantile(v, c(0.025, 0.975), na.rm = TRUE)
  c1 <- ci(bs[1, ]); c2 <- ci(bs[2, ]); c3 <- ci(lb)
  cat(sprintf("  %-18s POD %.2f (%.2f–%.2f) · false alarms %.2f/day (%.2f–%.2f) · median lead %.0f min (%.0f–%.0f)\n",
              label, pt[1], c1[1], c1[2], pt[2], c2[1], c2[2],
              median(ldv), c3[1], c3[2]))
  list(data.table(indicator = paste(label, "POD"), value = round(pt[1], 3),
                  ci_low = round(c1[1], 3), ci_high = round(c1[2], 3)),
       data.table(indicator = paste(label, "false alarms per station-day"), value = round(pt[2], 2),
                  ci_low = round(c2[1], 2), ci_high = round(c2[2], 2)),
       data.table(indicator = paste(label, "first-trigger lead (min)"), value = round(median(ldv)),
                  ci_low = round(c3[1]), ci_high = round(c3[2])))
}
for (a in list(list(360L, 1L, 0.80, "T+6 POD 0.80"), list(360L, 1L, 0.70, "T+6 POD 0.70"),
               list(180L, 1L, 0.80, "T+3 POD 0.80"), list(360L, 2L, 0.55, "T+6 persistence k=2"))) {
  rr <- boot_trig(a[[1]], a[[2]], a[[3]], a[[4]])
  for (x in rr) { k <- k + 1L; rows[[k]] <- x }
}

S <- rbindlist(rows)
fwrite(S, file.path(MART, "nowcast_bootstrap_ci.csv"), bom = TRUE)
cat("\n=== Summary table ===\n"); print(S)
cat("Saved: 02_MART/nowcast_bootstrap_ci.csv\n")
sink(); close(con)
cat("Log:", log_path, "\n")
