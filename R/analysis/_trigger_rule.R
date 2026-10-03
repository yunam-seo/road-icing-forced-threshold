# Defines the operational trigger rule: the largest probability threshold that reaches a target POD, found by bisection.
# Usage: tg <- make_trigger(dt_eval, step_min); th <- tg$thr_for_pod(pr, H, k, pod); tg$trig(pr, H, k, th).

suppressPackageStartupMessages(library(data.table))

make_trigger <- function(dt_eval, step_min = 30L) {
  D <- as.data.table(dt_eval)[, .(stn, tm, tno)]
  setorder(D, stn, tm)
  ppd <- 24 * 60 / step_min

  trig <- function(pr, H, k, thr) {
    dt <- copy(D)[, p := pr]
    dt[, fire := {
      a <- p >= thr
      if (k == 1L) a else {
        ok <- a
        for (j in 1:(k - 1L)) {
          prev <- shift(a, j) & (tm - shift(tm, j) == j * step_min)
          ok <- ok & (prev %in% TRUE)
        }
        ok }
    }, by = stn]
    on <- dt[!is.na(tno) & tno > 0 & tno <= H]
    on[, oid := paste0(stn, "_", tm + tno)]
    n_onset <- uniqueN(on$oid)
    fired <- on[fire == TRUE]
    lead <- if (nrow(fired)) fired[, .(lead = max(tno)), by = oid]$lead else numeric(0)
    list(pod = length(lead) / n_onset,
         far = mean(dt[(is.na(tno) | tno > H)]$fire) * ppd,
         lead = lead, n_onset = n_onset)
  }

  thr_for_pod <- function(pr, H, k, target, iter = 40L) {
    lo <- min(pr, na.rm = TRUE); hi <- max(pr, na.rm = TRUE)
    if (trig(pr, H, k, lo)$pod < target) return(NA_real_)
    if (trig(pr, H, k, hi)$pod >= target) return(hi)
    for (i in seq_len(iter)) {
      mid <- (lo + hi) / 2
      if (trig(pr, H, k, mid)$pod >= target) lo <- mid else hi <- mid
    }
    lo
  }

  best_reachable <- function(pr, H, k) {
    r <- trig(pr, H, k, min(pr, na.rm = TRUE)); r$thr <- min(pr, na.rm = TRUE); r
  }

  list(trig = trig, thr_for_pod = thr_for_pod, best_reachable = best_reachable, ppd = ppd)
}
