# Tests whether the relaxation time depends on the detrending bandwidth (30, 60, 120, 240 min).
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Outputs: 02_MART/lambda_bandwidth_sens.csv, lambda_bandwidth_paired.csv, variance_bandwidth_paired.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
DDIR <- rwis_dir()
ARGS<-commandArgs(trailingOnly=TRUE)
NSTN <- if(length(ARGS)>=1L) as.integer(ARGS[1]) else NA_integer_
W <- 720L; RW <- 360L; STEP <- 30L
BWS <- c(30, 60, 120, 240)
cov_thr <- 0.8; ctrl_per <- 2L; ctrl_buf <- 1440L
fmon <- c(11L, 12L, 1L, 2L, 3L)

map_cls <- function(s) {o <- rep(NA_character_, length(s)); o[s == 1L] <- "DRY"
  o[s %in% c(2L, 11:14)] <- "WET"; o[s == 3L] <- "SNOW"; o[s %in% c(4L, 15L)] <- "ICE"; o}
read_ts <- function(f) {
  d <- fread(f, select = c("datetime", "TR1", "SFS", "mon"),
             na.strings = c("-999", "-999.0", ""))
  d <- d[mon %in% fmon]; d <- filter_winter(d)
  if (!nrow(d)) return(NULL)
  d[(TR1 < -40) | (TR1 > 40), TR1 := NA_real_]
  d[, tm_min := as.integer(as.POSIXct(strptime(sprintf("%012.0f", as.numeric(datetime)),
     "%Y%m%d%H%M", tz = "UTC"))) %/% 60L]
  d <- d[!is.na(tm_min)]; setorder(d, tm_min)
  d[, cls := map_cls(SFS)]; setkey(d, tm_min); d[, .(tm_min, ts = TR1, cls)]}
gtrend <- function(v, bw) {ok <- which(!is.na(v))
  if (length(ok) < 3) return(rep(NA_real_, length(v)))
  sapply(seq_along(v), function(i) {w <- exp(-0.5 * ((ok - i) / bw)^2)
    sum(w * v[ok]) / sum(w)})}

end_stats <- function(v, bw) {
  v <- tail(v, W + 1L)
  if (mean(!is.na(v)) < cov_thr) return(c(NA_real_, NA_real_))
  r <- v - gtrend(v, bw)
  st <- seq(1L, length(v) - RW, by = STEP)
  aa <- vv <- rep(NA_real_, length(st))
  for (k in seq_along(st)) {
    seg <- r[st[k]:(st[k] + RW - 1L)]
    vv[k] <- var(seg, na.rm = TRUE)
    a <- seg[-length(seg)]; b <- seg[-1]; g <- is.finite(a) & is.finite(b)
    aa[k] <- if (sum(g) > 30) suppressWarnings(cor(a[g], b[g])) else NA_real_}
  af <- aa[is.finite(aa)]; vf <- vv[is.finite(vv)]
  c(if (length(af)) tail(af, 1) else NA_real_,
    if (length(vf)) tail(vf, 1) else NA_real_)}

ev <- readRDS(file.path(MART, "intl_event_table.rds")); setkey(ev, stn)
files <- list.files(DDIR, pattern = "\\.csv$", full.names = TRUE)
sof <- as.integer(sub(".*rwis_(\\d+)_.*", "\\1", basename(files)))
sel <- which(sof %in% ev$stn)
if(!is.na(NSTN)) sel <- sel[1:min(NSTN,length(sel))]

acc <- list()
for (bw in BWS) acc[[as.character(bw)]] <- list(on = numeric(0), ct = numeric(0),
                                                von = numeric(0), vct = numeric(0),
                                                son = integer(0), sct = integer(0))
n_on <- n_ct <- 0L

for (ii in sel) {
  d <- try(read_ts(files[ii]), silent = TRUE)
  if (inherits(d, "try-error") || is.null(d) || nrow(d) < 5000) next
  sid <- sof[ii]; evs <- ev[.(sid)][!is.na(onset_min)]
  if (!nrow(evs)) next
  grab <- function(t_end) {g <- data.table(tm_min = (t_end - W):t_end, key = "tm_min")
    d[g, on = "tm_min"]$ts}
  ons <- lapply(evs$onset_min - 1L, grab)
  oo <- sort(unique(evs$onset_min))
  cand <- d[cls %in% c("DRY", "WET", "SNOW") & tm_min > d$tm_min[1] + W]
  cts <- list()
  if (nrow(cand)) {
    jj <- findInterval(cand$tm_min, oo)
    dl <- ifelse(jj >= 1, abs(cand$tm_min - oo[pmax(jj, 1L)]), Inf)
    dh <- ifelse(jj < length(oo), abs(cand$tm_min - oo[pmin(jj + 1L, length(oo))]), Inf)
    far <- cand[pmin(dl, dh) >= ctrl_buf]
    nc <- min(nrow(far), ctrl_per * nrow(evs))
    if (nc > 0) cts <- lapply(far[sample(.N, nc)]$tm_min, grab)}
  n_on <- n_on + length(ons); n_ct <- n_ct + length(cts)
  for (bw in BWS) {
    k <- as.character(bw)
    if (length(ons)) {E <- vapply(ons, end_stats, c(0, 0), bw = bw)
      acc[[k]]$on <- c(acc[[k]]$on, E[1, ]); acc[[k]]$von <- c(acc[[k]]$von, E[2, ])
      acc[[k]]$son <- c(acc[[k]]$son, rep(sid, length(ons)))}
    if (length(cts)) {E <- vapply(cts, end_stats, c(0, 0), bw = bw)
      acc[[k]]$ct <- c(acc[[k]]$ct, E[1, ]); acc[[k]]$vct <- c(acc[[k]]$vct, E[2, ])
      acc[[k]]$sct <- c(acc[[k]]$sct, rep(sid, length(cts)))}}
  cat(sprintf("  station %d  onsets %d  controls %d\n", sid, length(ons), length(cts)))
}

lam <- function(a) {a <- a[is.finite(a) & a > 0 & a < 1]; -1 / log(a)}
lam_keep <- function(a) ifelse(is.finite(a) & a > 0 & a < 1, -1 / log(a), NA_real_)
out <- rbindlist(lapply(BWS, function(bw) {
  k <- as.character(bw)
  lo <- lam(acc[[k]]$on); lc <- lam(acc[[k]]$ct)
  ao <- acc[[k]]$on[is.finite(acc[[k]]$on)]; ac <- acc[[k]]$ct[is.finite(acc[[k]]$ct)]
  p <- if (length(lo) > 5 && length(lc) > 5)
         suppressWarnings(wilcox.test(lo, lc)$p.value) else NA_real_
  data.table(bandwidth_min = bw, n_onset = length(lo), n_control = length(lc),
             ar1_onset = round(median(ao), 4), ar1_control = round(median(ac), 4),
             inv_lambda_onset = round(median(lo), 1), inv_lambda_control = round(median(lc), 1),
             p_value = signif(p, 3))}))

cat("\n== Dependence of relaxation time λ⁻¹ on the detrending bandwidth ==\n"); print(out)
r60 <- out[bandwidth_min == 60]
cat(sprintf("\nbw=60 (main text) λ⁻¹ onset %.1f min\n", r60$inv_lambda_onset))
cat(sprintf("Onset λ⁻¹ for bw 30→240: %s min\n",
            paste(out$inv_lambda_onset, collapse = " → ")))
cat(sprintf("λ⁻¹/bw ratio: %s\n",
            paste(round(out$inv_lambda_onset / out$bandwidth_min, 2), collapse = " · ")))
cat("
== Station-paired (median per station → paired signed-rank test) ==
")
pair <- rbindlist(lapply(BWS, function(bw) {
  k <- as.character(bw)
  stopifnot(length(acc[[k]]$son) == length(acc[[k]]$on),
            length(acc[[k]]$sct) == length(acc[[k]]$ct))
  O <- data.table(stn = acc[[k]]$son, v = lam_keep(acc[[k]]$on))
  C <- data.table(stn = acc[[k]]$sct, v = lam_keep(acc[[k]]$ct))
  O <- O[is.finite(v), .(on = median(v)), by = stn]
  C <- C[is.finite(v), .(ct = median(v)), by = stn]
  S <- merge(O, C, by = "stn")
  p <- if (nrow(S) > 5) suppressWarnings(
         wilcox.test(S$on, S$ct, paired = TRUE)$p.value) else NA_real_
  data.table(bandwidth_min = bw, n_stations = nrow(S),
             onset = round(median(S$on), 1), control = round(median(S$ct), 1),
             diff_median = round(median(S$on - S$ct), 2),
             n_onset_greater = sum(S$on > S$ct), p_value = signif(p, 3))}))
print(pair)
cat(sprintf("
Direction consistency: sign of diff_median = %s
",
            paste(ifelse(pair$diff_median > 0, "+", "-"), collapse = " ")))
cat("
== Stationary variance (window-end variance), same approach ==
")
pv <- rbindlist(lapply(BWS, function(bw) {
  k <- as.character(bw)
  O <- data.table(stn = acc[[k]]$son, v = acc[[k]]$von)[is.finite(v), .(on = median(v)), by = stn]
  C <- data.table(stn = acc[[k]]$sct, v = acc[[k]]$vct)[is.finite(v), .(ct = median(v)), by = stn]
  S <- merge(O, C, by = "stn")
  p <- if (nrow(S) > 5) suppressWarnings(
         wilcox.test(S$on, S$ct, paired = TRUE)$p.value) else NA_real_
  data.table(bandwidth_min = bw, n_stations = nrow(S),
             onset = round(median(S$on), 4), control = round(median(S$ct), 4),
             diff_median = round(median(S$on - S$ct), 4),
             n_onset_greater = sum(S$on > S$ct), p_value = signif(p, 3))}))
print(pv)
cat(sprintf("
Stationary variance, direction signs = %s
",
            paste(ifelse(pv$diff_median > 0, "+", "-"), collapse = " ")))
pv[, n_stations := length(sel)]
pair[, n_stations := length(sel)]
out[, n_stations := length(sel)]
fwrite(pv, file.path(MART, "variance_bandwidth_paired.csv"))
fwrite(pair, file.path(MART, "lambda_bandwidth_paired.csv"))
fwrite(out, file.path(MART, "lambda_bandwidth_sens.csv"))
cat("Saved: 02_MART/lambda_bandwidth_sens.csv\n")
