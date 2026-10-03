# Computes the decay-time (tau_d) sweep of the cold dose at onset and 30 and 60 min before, with station-cluster bootstrap CIs and hour-stratified AUC, and summarises onset hour by prior state.
# Inputs: 02_MART/intl_cumstress.rds, intl_event_table.rds, intl_lagged_table.rds; 1-min RWIS files.
# Outputs: 02_MART/tau_sweep_lagged.csv, tau_sweep_hour.csv, tau_sweep_fig2c.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
source(code_file("07_intl_cumstress.R"))
TAUS <- c(1, 2, 5, 10, 15, 20, 30, 60, 180, 1440, 4320)
LAGS <- c(0L, 30L, 60L)
CS <- as.data.table(readRDS(file.path(MART, "intl_cumstress.rds")))[, .(stn, tm_min, type)]
files <- list.files(DDIR, pattern = "[.]csv$", full.names = TRUE)
sof <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(files)))
accum <- function(tm, x, tau, gap) {
  n <- length(tm); dt <- c(0, diff(tm)); u <- ifelse(!is.na(dt) & dt <= gap & !is.na(x), x * dt, 0)
  E <- cumsum(dt / tau); blk <- floor(E / 500); S <- numeric(n); carry <- 0; Eend <- 0
  for (b in unique(blk)) { ii <- which(blk == b); E0 <- E[ii[1]]
    acc <- carry * exp(-(E0 - Eend)) + cumsum(u[ii] * exp(E[ii] - E0)); S[ii] <- acc * exp(-(E[ii] - E0))
    carry <- S[ii[length(ii)]]; Eend <- E[ii[length(ii)]] }
  S }

R <- list()
for (i in seq_along(files)) { sid <- sof[i]; q <- CS[stn == sid]; if (!nrow(q)) next
  d <- read_station(files[i]); if (is.null(d) || nrow(d) < 10L) next
  xc <- pmax(0, -d$ts); out <- data.table(tm_min = d$tm_min, ts = d$ts)
  for (tau in TAUS) out[[paste0("t", tau)]] <- accum(d$tm_min, xc, tau, cfg$gap_max)
  out[, obs_tm := tm_min]; setkey(out, tm_min)
  for (L in LAGS) {
    j <- out[.(q$tm_min - L), roll = "nearest", on = "tm_min"]
    ok <- abs(j$obs_tm - (q$tm_min - L)) <= 2L
    blk <- cbind(q, lag = L, j[, -c("tm_min", "obs_tm")]); blk[!ok, (c("ts", paste0("t", TAUS))) := NA_real_]
    R[[length(R) + 1L]] <- blk } }
R <- rbindlist(R)
p <- as.POSIXct(R$tm_min * 60, origin = "1970-01-01", tz = "UTC")
auc <- function(s, pos) { o <- is.finite(s); s <- s[o]; pos <- pos[o]; r <- rank(s)
  np <- as.numeric(sum(pos)); nn <- as.numeric(sum(!pos)); if (!np || !nn) return(NA_real_)
  (sum(r[pos]) - np * (np + 1) / 2) / (np * nn) }

R[, hour_raw := as.integer(format(p, "%H"))]

B <- 500L; stns <- unique(R$stn)
res <- list()
for (L in LAGS) { X <- R[lag == L]; pos <- X$type == "onset"
  base <- sapply(TAUS, function(t) auc(X[[paste0("t", t)]], pos)); ib <- which.max(base)
  idx <- split(seq_len(nrow(X)), X$stn)
  BS <- matrix(NA_real_, B, length(TAUS))
  for (b in seq_len(B)) { s <- sample(names(idx), length(idx), replace = TRUE); ii <- unlist(idx[s], use.names = FALSE)
    BS[b, ] <- sapply(TAUS, function(t) auc(X[[paste0("t", t)]][ii], pos[ii])) }
  D <- BS - BS[, ib]
  res[[length(res) + 1L]] <- data.table(evaluation_time = sprintf("t0-%d", L), tau_min = TAUS, AUC = base,
    ci_low = apply(BS, 2, quantile, .025), ci_high = apply(BS, 2, quantile, .975),
    diff_from_best = base - base[ib], diff_ci_low = apply(D, 2, quantile, .025), diff_ci_high = apply(D, 2, quantile, .975),
    auc_instant_ts = auc(-X$ts, pos), auc_instant_subzero = auc(pmax(0, -X$ts), pos), n = sum(is.finite(X$t60)))
}
S1 <- rbindlist(res); print(S1, digits = 4)
fwrite(S1, file.path(MART, "tau_sweep_lagged.csv"), bom = TRUE)

X0 <- R[lag == 0L]; X0[, hb := (hour_raw %/% 3L) * 3L]
cat("\nHour distribution (raw time stamps, 3-h bins, %)\n")
print(dcast(X0[, .N, by = .(type, hb)][, pct := round(100 * N / sum(N), 1), by = type], hb ~ type, value.var = "pct"))
strat <- function(col) { g <- X0[, .(a = auc(get(col), type == "onset"), w = sum(type == "onset")), by = hb][is.finite(a)]
  sum(g$a * g$w) / sum(g$w) }
S3 <- data.table(tau_min = TAUS, auc_overall = sapply(TAUS, function(t) auc(X0[[paste0("t", t)]], X0$type == "onset")),
                 auc_hour_stratified = sapply(TAUS, function(t) strat(paste0("t", t))))
S3 <- rbind(S3, data.table(tau_min = NA, auc_overall = auc(-X0$ts, X0$type == "onset"), auc_hour_stratified = strat("ts") * -1 + 1))
cat("\nHour-stratified AUC (last row = instant −Ts)\n"); print(S3, digits = 4)
fwrite(S3, file.path(MART, "tau_sweep_hour.csv"), bom = TRUE)

lg <- readRDS(file.path(MART, "intl_lagged_table.rds")); ev <- readRDS(file.path(MART, "intl_event_table.rds"))
N4 <- lg[lag <= 1440, .(n = sum(!is.na(ts)), md = median(ts, na.rm = TRUE)), by = .(pre_state, lag)]
W4 <- dcast(N4, lag ~ pre_state, value.var = c("md", "n"))
W4[, snow_lowest := md_SNOW <= pmin(md_DRY, md_WET)]
print(W4[order(-lag)])
eh <- as.integer(format(as.POSIXct(ev$onset_min * 60, origin = "1970-01-01", tz = "UTC"), "%H"))
cat("\nOnset hour by prior state (raw time stamps): median and quartiles\n")
print(data.table(pre = ev$pre_state, h = as.numeric(eh))[, .(n = .N, q1 = as.numeric(quantile(h, .25)), med = as.numeric(median(h)), q3 = as.numeric(quantile(h, .75)),
  night_18_06_pct = round(100 * mean(h >= 18 | h < 6), 1)), by = pre])
fwrite(W4, file.path(MART, "tau_sweep_fig2c.csv"), bom = TRUE)
cat("Saved: tau_sweep_lagged.csv · tau_sweep_hour.csv · tau_sweep_fig2c.csv\n")
