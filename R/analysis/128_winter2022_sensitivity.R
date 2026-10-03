# Repeats the CSD discrimination analysis with the incomplete 2022/23 winter included and compares it with the main analysis window.
# Inputs: 1-min RWIS files; meta/stn_exclusion_intl.csv.
# Output: 02_MART/winter2022_sensitivity.csv.
# Reads the full 1-min archive twice (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
EXF <- meta_file("stn_exclusion_intl.csv")

source(code_file("03_intl_events_cleaned.R")); CFG_EV <- cfg
source(code_file("09_intl_ews.R"));            CFG_EWS <- cfg

files0 <- list.files(CFG_EV$data_dir, pattern = "[.]csv$", full.names = TRUE)
EX <- fread(EXF)
sof0 <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(files0)))
files0 <- files0[!(sof0 %in% EX$STN_ID)]
cat(sprintf("Cleaned stations %d (excluded %d)\n", length(files0), nrow(EX)))

auc <- function(s, p) {
  o <- is.finite(s); s <- s[o]; p <- p[o]; r <- rank(s)
  (sum(r[p]) - sum(p) * (sum(p) + 1) / 2) / (sum(p) * sum(!p))
}

build_ews <- function(ev, files) {
  setkey(ev, stn)
  sof <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(files)))
  out <- list(); k <- 0L
  for (i in seq_along(files)) {
    sid <- sof[i]; evs <- ev[stn == sid]; if (!nrow(evs)) next
    d <- read_ts(files[i]); if (is.null(d) || nrow(d) < cfg$W) next
    setkey(d, tm_min); tmin <- d$tm_min[1]
    grab <- function(t_end) {
      g <- data.table(tm_min = (t_end - cfg$W):t_end, key = "tm_min")
      d[g, on = "tm_min"]$ts
    }
    for (j in seq_len(nrow(evs))) {
      e <- ews_window(grab(evs$onset_min[j] - 1L))
      k <- k + 1L
      out[[k]] <- data.table(type = "onset", stn = sid, t_end = evs$onset_min[j] - 1L, t(e))
    }
    oo <- sort(unique(evs$onset_min))
    cand <- d[cls %in% c("DRY", "WET", "SNOW") & tm_min > tmin + cfg$W]
    if (nrow(cand)) {
      jj <- findInterval(cand$tm_min, oo)
      dl <- ifelse(jj >= 1, abs(cand$tm_min - oo[pmax(jj, 1L)]), Inf)
      dh <- ifelse(jj < length(oo), abs(cand$tm_min - oo[pmin(jj + 1L, length(oo))]), Inf)
      far <- cand[pmin(dl, dh) >= cfg$ctrl_buffer]
      nc <- min(nrow(far), cfg$ctrl_per_event * nrow(evs))
      if (nc > 0) for (t_end in far[sample(.N, nc)]$tm_min) {
        e <- ews_window(grab(t_end))
        k <- k + 1L
        out[[k]] <- data.table(type = "control", stn = sid, t_end = t_end, t(e))
      }
    }
    if (i %% 20L == 0L) message("  ", i, "/", length(files))
  }
  rbindlist(out, fill = TRUE)
}

run_one <- function(wmin, label) {
  assign("WINTER_MIN", as.integer(wmin), envir = globalenv())
  cat(sprintf("\n===== %s (WINTER_MIN=%d) =====\n", label, wmin))
  .qc$consist <- 0L; .qc$spike <- 0L; .qc$ice_onset <- 0L
  set.seed(2026)
  assign("cfg", CFG_EV, envir = globalenv())
  ev <- build_intl(files0)
  cat(sprintf("Events %d · stations %d\n", nrow(ev), length(unique(ev$stn))))
  set.seed(2026)
  assign("cfg", CFG_EWS, envir = globalenv())
  res <- build_ews(ev, files0)
  pp <- res$type == "onset"
  data.table(analysis = label, winter_window = wmin, n_events = nrow(ev),
             var_AUC = auc(res$var_tau, pp), ar1_AUC = auc(res$ar1_tau, pp))
}

t0 <- Sys.time()
a <- run_one(2023L, "main analysis window (2023/24–2025/26)")
b <- run_one(2022L, "including the 2022/23 winter")

S <- rbind(a, b)
S[, `:=`(var_AUC = round(var_AUC, 4), ar1_AUC = round(ar1_AUC, 4))]
S <- rbind(S, data.table(analysis = "absolute change", winter_window = NA_integer_,
                         n_events = b$n_events - a$n_events,
                         var_AUC = round(b$var_AUC - a$var_AUC, 4),
                         ar1_AUC = round(b$ar1_AUC - a$ar1_AUC, 4)))
print(S)
fwrite(S, file.path(MART, "winter2022_sensitivity.csv"), bom = TRUE)
cat(sprintf("\nElapsed %.1f min\nSaved: 02_MART/winter2022_sensitivity.csv\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
