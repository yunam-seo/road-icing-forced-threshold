# Measures how much the CSD discrimination AUC changes when only the random seed for control times is changed.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files; meta/stn_exclusion_intl.csv.
# Output: 02_MART/control_seed_spread.csv.
# Reads the full 1-min archive several times (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
EXF <- meta_file("stn_exclusion_intl.csv")
SEEDS <- c(2026L, 2027L, 2028L)

source(code_file("09_intl_ews.R"))

files <- list.files(DDIR, pattern = "[.]csv$", full.names = TRUE)
EX <- fread(EXF)
sof <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(files)))
files <- files[!(sof %in% EX$STN_ID)]
sof <- sof[!(sof %in% EX$STN_ID)]
cat(sprintf("Stations %d (excluded %d) · seeds %s\n", length(files), nrow(EX),
            paste(SEEDS, collapse = ", ")))

auc <- function(s, p) {
  o <- is.finite(s); s <- s[o]; p <- p[o]; r <- rank(s)
  (sum(r[p]) - sum(p) * (sum(p) + 1) / 2) / (sum(p) * sum(!p))
}

ev <- readRDS(file.path(MART, "intl_event_table.rds"))
setDT(ev); setkey(ev, stn)

run_seed <- function(seed) {
  set.seed(seed)
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
      k <- k + 1L; out[[k]] <- data.table(type = "onset", t(e))
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
        k <- k + 1L; out[[k]] <- data.table(type = "control", t(e))
      }
    }
    if (i %% 20L == 0L) message("  seed ", seed, ": ", i, "/", length(files))
  }
  res <- rbindlist(out, fill = TRUE); pp <- res$type == "onset"
  data.table(seed = seed, n_control_times = sum(!pp),
             var_AUC = auc(res$var_tau, pp), ar1_AUC = auc(res$ar1_tau, pp))
}

t0 <- Sys.time()
S <- rbindlist(lapply(SEEDS, run_seed))
S[, `:=`(var_AUC = round(var_AUC, 4), ar1_AUC = round(ar1_AUC, 4))]
print(S)
cat(sprintf("\nSpread from changing the seed only — variance %.4f · AR1 %.4f\n",
            diff(range(S$var_AUC)), diff(range(S$ar1_AUC))))
fwrite(S, file.path(MART, "control_seed_spread.csv"), bom = TRUE)
cat(sprintf("Elapsed %.1f min\nSaved: 02_MART/control_seed_spread.csv\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
