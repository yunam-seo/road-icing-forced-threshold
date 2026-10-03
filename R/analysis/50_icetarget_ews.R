# Applies the event definition and CSD analysis unchanged to the independent 127-station network, with |Ts-Ta| QC variants.
# Input: 1-min records of the independent network (RICE_VALID_FILE).
# Outputs: 02_MART/icetarget_event_table.rds, icetarget_ews.rds, icetarget_ews_thr_sens.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))

MART <- file.path(ROOT, "02_MART")
ICE  <- valid_file()

E3 <- new.env(); sys.source(code_file("03_intl_events_cleaned.R"), envir = E3)
E9 <- new.env(); sys.source(code_file("09_intl_ews.R"),            envir = E9)
cfg3 <- E3$cfg; cfg9 <- E9$cfg
cat("[Reused] event definition p_min=", cfg3$p_min, " g_min=", cfg3$g_min,
    " / EWS W=", cfg9$W, " rw=", cfg9$rw, " bw=", cfg9$bw, "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
test <- length(args) >= 1L && args[1] == "test"

cat("\n[1] Loading ICE_TARGET ...\n")
d <- readRDS(ICE); setDT(d)
keep <- c("STN_ID", "tnum", "year", "mon", "TR1", "TA", "HM", "SFS", "IFT")
d <- d[, ..keep]; invisible(gc(FALSE))
n_all <- nrow(d)
d <- d[mon %in% cfg3$freeze_mon]
source(code_file("00_analysis_window.R"))
n_pre <- nrow(d)
d <- filter_winter(d, ym = list(y = d$year, m = d$mon))
cat(sprintf("  Winter window (>=%d) applied: %s -> %s rows\n", WINTER_MIN,
            format(n_pre, big.mark = ","), format(nrow(d), big.mark = ",")))
d[, STN_ID := as.integer(STN_ID)]
d[, tm_min := as.integer(tnum %/% 60)]
d[, tnum := NULL]
setorder(d, STN_ID, tm_min)
n_dup <- sum(duplicated(d, by = c("STN_ID", "tm_min")))
if (n_dup > 0L) d <- unique(d, by = c("STN_ID", "tm_min"))
stns <- sort(unique(d$STN_ID)); if (test) stns <- head(stns, 5L)
cat(sprintf("  All %s rows -> freezing season %s rows · %d stations (%d duplicates removed)\n",
            format(n_all, big.mark = ","), format(nrow(d), big.mark = ","),
            length(unique(d$STN_ID)), n_dup))

prep_stn <- function(s, thr) {
  s <- copy(s)
  s[(TR1 < cfg3$ts_range[1]) | (TR1 > cfg3$ts_range[2]), TR1 := NA_real_]
  s[(TA  < cfg3$ta_range[1]) | (TA  > cfg3$ta_range[2]), TA  := NA_real_]
  s[(HM  < cfg3$hm_range[1]) | (HM  > cfg3$hm_range[2]), HM  := NA_real_]
  if (is.finite(thr))
    s[!is.na(TR1) & !is.na(TA) & abs(TR1 - TA) > thr, TR1 := NA_real_]
  s[, dprev := TR1 - shift(TR1)]; s[, dnext := shift(TR1, type = "lead") - TR1]
  s[!is.na(dprev) & !is.na(dnext) & abs(dprev) > cfg3$spike_max &
    abs(dnext) > cfg3$spike_max & sign(dprev) != sign(dnext), TR1 := NA_real_]
  s[, cls := E3$map_sfs_class(as.integer(SFS))]
  s[, .(stn = STN_ID, tm_min, ts = TR1, ta = TA, hm = HM, ift = IFT, cls)]
}

run_variant <- function(thr, lab) {
  cat(sprintf("\n[Variant %s] |Ts-Ta| threshold %s\n", lab,
              if (is.finite(thr)) paste0(thr, "℃") else "not applied"))
  set.seed(2026)
  EV <- list(); WS <- list(); ke <- kw <- 0L
  for (i in seq_along(stns)) {
    sid <- stns[i]
    p <- prep_stn(d[STN_ID == sid], thr)
    ev <- E3$extract_events_station(p)
    if (is.null(ev) || !nrow(ev)) next
    ke <- ke + 1L; EV[[ke]] <- ev

    dd <- p[, .(tm_min, ts, cls)]; setkey(dd, tm_min)
    if (nrow(dd) < cfg9$W) next
    tmin <- dd$tm_min[1]
    grab <- function(t_end) {
      g <- data.table(tm_min = (t_end - cfg9$W):t_end, key = "tm_min")
      dd[g, on = "tm_min"]$ts
    }
    for (j in seq_len(nrow(ev))) {
      e <- E9$ews_window(grab(ev$onset_min[j] - 1L))
      kw <- kw + 1L
      WS[[kw]] <- data.table(type = "onset", stn = sid,
                             eid = paste0(sid, "_", ev$onset_min[j]),
                             t_end = ev$onset_min[j] - 1L, t(e))
    }
    oo <- sort(unique(ev$onset_min))
    cand <- dd[cls %in% c("DRY", "WET", "SNOW") & tm_min > tmin + cfg9$W]
    if (nrow(cand)) {
      jj <- findInterval(cand$tm_min, oo)
      dl <- ifelse(jj >= 1, abs(cand$tm_min - oo[pmax(jj, 1L)]), Inf)
      dh <- ifelse(jj < length(oo), abs(cand$tm_min - oo[pmin(jj + 1L, length(oo))]), Inf)
      far <- cand[pmin(dl, dh) >= cfg9$ctrl_buffer]
      nc <- min(nrow(far), cfg9$ctrl_per_event * nrow(ev))
      if (nc > 0) for (t_end in far[sample(.N, nc)]$tm_min) {
        e <- E9$ews_window(grab(t_end))
        kw <- kw + 1L
        WS[[kw]] <- data.table(type = "control", stn = sid, eid = NA_character_,
                               t_end = t_end, t(e))
      }
    }
    if (i %% 20L == 0L) message("  ", i, "/", length(stns))
  }
  list(ev = rbindlist(EV, fill = TRUE), ws = rbindlist(WS, fill = TRUE))
}

auc <- function(score, pos) {
  o <- is.finite(score); score <- score[o]; pos <- pos[o]
  if (!length(unique(pos)) == 2L) return(NA_real_)
  r <- rank(score); n1 <- sum(pos); n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
summarize <- function(ws, lab) {
  cat(sprintf("\n=== [%s] EWS summary ===\n", lab))
  print(ws[, .(n = .N,
               var_tau = round(median(var_tau, na.rm = TRUE), 3),
               pos_var = round(mean(var_tau > 0, na.rm = TRUE), 2),
               ar1_tau = round(median(ar1_tau, na.rm = TRUE), 3),
               pos_ar1 = round(mean(ar1_tau > 0, na.rm = TRUE), 2)), by = type][order(type)])
  for (v in c("var_tau", "ar1_tau")) {
    x <- ws[type == "onset"][[v]]; y <- ws[type == "control"][[v]]
    p1 <- suppressWarnings(wilcox.test(x, y)$p.value)
    xf <- x[is.finite(x)]
    p2 <- suppressWarnings(wilcox.test(xf, mu = 0, alternative = "greater")$p.value)
    cat(sprintf("  %s: onset vs control Wilcoxon p=%.3g · onset τ>0 sign test p=%.3g · AUC %.3f\n",
                v, p1, p2, auc(ws[[v]], ws$type == "onset")))
  }
}

VAR <- list(A = list(thr = Inf, lab = "A direct transfer (not applied)"),
            B = list(thr = 30,  lab = "B main threshold 30 °C"),
            C = list(thr = 15,  lab = "C stricter 15℃"),
            D = list(thr = 10,  lab = "D stricter 10℃"))
OUT <- list()
for (k in names(VAR)) {
  r <- run_variant(VAR[[k]]$thr, VAR[[k]]$lab)
  cat(sprintf("\n[%s] %d events · %d stations\n", k, nrow(r$ev), uniqueN(r$ev$stn)))
  if (k == "A") print(table(r$ev$pre_state))
  summarize(r$ws, VAR[[k]]$lab)
  OUT[[k]] <- r
}

cat("\n=== Threshold sensitivity summary ===\n")
print(rbindlist(lapply(names(VAR), function(k) {
  w <- OUT[[k]]$ws
  data.table(variant = VAR[[k]]$lab, n_events = nrow(OUT[[k]]$ev),
             var_AUC = round(auc(w$var_tau, w$type == "onset"), 3),
             ar1_AUC = round(auc(w$ar1_tau, w$type == "onset"), 3))
})))

if (!test) {
  saveRDS(OUT$A$ev, file.path(MART, "icetarget_event_table.rds"))
  saveRDS(OUT$A$ws, file.path(MART, "icetarget_ews.rds"))
  saveRDS(lapply(OUT, `[[`, "ws"), file.path(MART, "icetarget_ews_thr_sens.rds"))
  cat("\nSaved: 02_MART/icetarget_{event_table,ews,ews_thr_sens}.rds\n")
}
