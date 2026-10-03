# Builds the pre-onset trajectory table: surface and air conditions at fixed lags before each onset, without interpolation.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Output: 02_MART/intl_lagged_table.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))

cfg <- list(
  data_dir = rwis_dir(),
  out_dir  = file.path(ROOT, "02_MART"),
  ts_range = c(-40,40), ta_range = c(-40,40), hm_range = c(0,100),
  freeze_mon = c(11L,12L,1L,2L,3L),
  lags = c(0L,15L,30L,60L,120L,180L,240L,360L,480L,600L,720L,960L,1200L,1440L),
  tol  = 5L
)

map_sfs_class <- function(sfs) {
  out <- rep(NA_character_, length(sfs))
  out[sfs==1L] <- "DRY"; out[sfs %in% c(2L,11L,12L,13L,14L)] <- "WET"
  out[sfs==3L] <- "SNOW"; out[sfs %in% c(4L,15L)] <- "ICE"; out
}

frost_point <- function(Ta, RH) {
  aw<-17.62; bw<-243.12; ai<-22.46; bi<-272.62
  e <- (pmax(RH,1e-3)/100) * 6.112*exp(aw*Ta/(bw+Ta))
  l <- log(e/6.112)
  bi*l/(ai-l)
}

read_station <- function(fpath) {
  d <- fread(fpath, select=c("datetime","STN_ID","TR1","TA","HM","RN05_YN","SFS","mon"),
             na.strings=c("-999","-999.0",""))
  d <- d[mon %in% cfg$freeze_mon]; d <- filter_winter(d)
  if (!nrow(d)) return(NULL)
  d[(TR1<cfg$ts_range[1])|(TR1>cfg$ts_range[2]), TR1:=NA_real_]
  d[(TA <cfg$ta_range[1])|(TA >cfg$ta_range[2]), TA :=NA_real_]
  d[(HM <cfg$hm_range[1])|(HM >cfg$hm_range[2]), HM :=NA_real_]
  d[, tm_min := as.integer(as.POSIXct(strptime(
        sprintf("%012.0f", as.numeric(datetime)), "%Y%m%d%H%M", tz="UTC"))) %/% 60L]
  d <- d[!is.na(tm_min)]
  res <- d[, .(tm_min, ts=TR1, ta=TA, rh=HM, rain=as.integer(RN05_YN>0))]
  setkey(res, tm_min); res
}

build_lagged <- function(nfile=NULL) {
  ev <- readRDS(file.path(cfg$out_dir, "intl_event_table.rds"))
  setkey(ev, stn)
  files <- list.files(cfg$data_dir, pattern="\\.csv$", full.names=TRUE)
  if (!is.null(nfile)) files <- head(files, nfile)
  stn_of <- as.integer(sub(".*rwis_(\\d+)_.*", "\\1", basename(files)))
  out <- list(); k <- 0L
  for (i in seq_along(files)) {
    sid <- stn_of[i]; evs <- ev[stn==sid]
    if (!nrow(evs)) next
    d <- read_station(files[i]); if (is.null(d)) next
    for (j in seq_len(nrow(evs))) {
      o <- evs$onset_min[j]
      targ <- o - cfg$lags
      q <- data.table(tm_min=targ, key="tm_min")
      m <- d[q, .(ttar=i.tm_min, tobs=x.tm_min, ts, ta, rh, rain),
             roll="nearest", on="tm_min"]
      m[abs(tobs - ttar) > cfg$tol, c("ts","ta","rh","rain") := NA]
      k <- k + 1L
      out[[k]] <- data.table(
        eid = paste0(sid,"_",o), stn=sid, onset_min=o,
        pre_state=evs$pre_state[j], lag=o - m$ttar,
        ts=m$ts, ta=m$ta, rh=m$rh, rain=m$rain)
    }
    if (i %% 20L == 0L) message("  ", i, "/", length(files))
  }
  lg <- rbindlist(out)
  lg[, m_liq := ts - 0]
  lg[, m_dep := ts - frost_point(ta, rh)]
  lg[]
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  test <- length(args)>=1L && args[1]=="test"
  lg <- build_lagged(if (test) 8L else NULL)
  cat("lag rows:", nrow(lg), " events:", uniqueN(lg$eid), " stations:", uniqueN(lg$stn), "\n")
  cat("Missing rate by lag (ts):\n")
  print(lg[, .(miss=round(mean(is.na(ts)),3)), by=lag][order(lag)])
  if (!test) {
    saveRDS(lg, file.path(cfg$out_dir, "intl_lagged_table.rds"))
    cat("Saved: 02_MART/intl_lagged_table.rds\n")
  }
}
