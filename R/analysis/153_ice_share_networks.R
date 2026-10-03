# Computes the share of 1-min ICE observations in the primary and independent networks under three definitions.
# Inputs: 1-min RWIS files; meta/stn_exclusion_intl.csv; 1-min records of the independent network (RICE_VALID_FILE).
# Output: 02_MART/ice_share_networks.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
source(code_file("00_analysis_window.R"))
DD <- rwis_dir()
EX <- fread(meta_file("stn_exclusion_intl.csv"))
ICE <- valid_file()
cls_of <- function(s) { o <- rep(NA_character_, length(s))
  o[s == 1L] <- "DRY"; o[s %in% c(2L, 11:14)] <- "WET"; o[s == 3L] <- "SNOW"; o[s %in% c(4L, 15L)] <- "ICE"; o }
summ <- function(stn, sfs, ts) {
  cl <- cls_of(as.integer(sfs)); okc <- !is.na(cl); okt <- okc & is.finite(ts) & ts >= -40 & ts <= 40
  data.table(stn = stn[1], n_cls = sum(okc), ice_cls = sum(cl[okc] == "ICE"),
             n_ts = sum(okt), ice_ts = sum(cl[okt] == "ICE")) }

fs <- list.files(DD, pattern = "[.]csv$", full.names = TRUE)
sid <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(fs)))
fs <- fs[!(sid %in% EX$STN_ID)]; sid <- sid[!(sid %in% EX$STN_ID)]
A <- rbindlist(lapply(seq_along(fs), function(i) {
  d <- fread(fs[i], select = c("datetime", "mon", "TR1", "SFS"), na.strings = c("-999", "-999.0", ""))
  d <- d[mon %in% c(11L, 12L, 1L, 2L, 3L)]; d <- filter_winter(d)
  if (!nrow(d)) return(NULL)
  summ(rep(sid[i], nrow(d)), d$SFS, d$TR1) }))
cat("Primary network stations", nrow(A), "\n")

d <- readRDS(ICE); setDT(d); d <- d[, .(STN_ID, year, mon, TR1, SFS)]; invisible(gc(FALSE))
d <- d[mon %in% c(11L, 12L, 1L, 2L, 3L)]
d <- filter_winter(d, ym = list(y = d$year, m = d$mon))
B <- d[, summ(STN_ID, SFS, TR1), by = STN_ID][, stn := NULL]
setnames(B, "STN_ID", "stn"); rm(d); invisible(gc(FALSE))
cat("Independent network stations", nrow(B), "\n")

res <- function(X, net) data.table(network = net, n_stations = nrow(X),
  ice_share_classified_pct = 100 * sum(X$ice_cls) / sum(X$n_cls),
  ice_share_station_mean_pct = 100 * mean(X$ice_cls / X$n_cls, na.rm = TRUE),
  ice_share_valid_ts_pct = 100 * sum(X$ice_ts) / sum(X$n_ts),
  n_classified_minutes = sum(X$n_cls), n_ice_minutes = sum(X$ice_cls))
O <- rbind(res(A, "primary"), res(B, "independent"))
print(O)
fwrite(O, file.path(MART, "ice_share_networks.csv"), bom = TRUE)
cat("Saved: 02_MART/ice_share_networks.csv\n")
