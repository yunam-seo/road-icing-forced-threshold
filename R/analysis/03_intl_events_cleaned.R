# Extracts 1-min icing onsets ({DRY, WET, SNOW} -> ICE) for the primary network after QC and station exclusion.
# Inputs: 1-min RWIS files; meta/stn_exclusion_intl.csv.
# Output: 02_MART/intl_event_table.{rds,csv}.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))

cfg <- list(
  data_dir = rwis_dir(),
  out_dir  = file.path(ROOT, "02_MART"),
  p_min = 10L, g_min = 60L,
  ts_range = c(-40, 40), ta_range = c(-40, 40), hm_range = c(0, 100),
  tsta_max = 30,
  spike_max = 10,
  ice_ts_max = 5,
  freeze_mon = c(11L, 12L, 1L, 2L, 3L)
)

.qc <- new.env(); .qc$consist <- 0L; .qc$spike <- 0L; .qc$ice_onset <- 0L

map_sfs_class <- function(sfs) {
  out <- rep(NA_character_, length(sfs))
  out[sfs == 1L] <- "DRY"
  out[sfs %in% c(2L, 11L, 12L, 13L, 14L)] <- "WET"
  out[sfs == 3L] <- "SNOW"
  out[sfs %in% c(4L, 15L)] <- "ICE"
  return(out)
}

read_cleaned_station <- function(fpath) {
  d <- fread(fpath, select = c("datetime","STN_ID","TR1","TA","HM","SFS","IFT",
                               "RN05_YN","mon"),
             na.strings = c("-999","-999.0",""))
  d <- d[mon %in% cfg$freeze_mon]
  d <- filter_winter(d)
  if (nrow(d) == 0L) return(NULL)
  d[, tm_min := as.integer(as.POSIXct(strptime(
        sprintf("%012.0f", as.numeric(datetime)), "%Y%m%d%H%M", tz="UTC"))) %/% 60L]
  d <- d[!is.na(tm_min)]
  setorder(d, tm_min)
  d[(TR1 < cfg$ts_range[1]) | (TR1 > cfg$ts_range[2]), TR1 := NA_real_]
  d[(TA  < cfg$ta_range[1]) | (TA  > cfg$ta_range[2]), TA  := NA_real_]
  d[(HM  < cfg$hm_range[1]) | (HM  > cfg$hm_range[2]), HM  := NA_real_]
  nc <- d[!is.na(TR1) & !is.na(TA) & abs(TR1 - TA) > cfg$tsta_max, .N]
  d[!is.na(TR1) & !is.na(TA) & abs(TR1 - TA) > cfg$tsta_max, TR1 := NA_real_]
  d[, dprev := TR1 - shift(TR1)]; d[, dnext := shift(TR1, type="lead") - TR1]
  ns <- d[!is.na(dprev) & !is.na(dnext) & abs(dprev) > cfg$spike_max &
          abs(dnext) > cfg$spike_max & sign(dprev) != sign(dnext), .N]
  d[!is.na(dprev) & !is.na(dnext) & abs(dprev) > cfg$spike_max &
    abs(dnext) > cfg$spike_max & sign(dprev) != sign(dnext), TR1 := NA_real_]
  d[, c("dprev","dnext") := NULL]
  .qc$consist <- .qc$consist + nc; .qc$spike <- .qc$spike + ns
  d[, cls := map_sfs_class(SFS)]
  return(d[, .(stn=STN_ID, tm_min, ts=TR1, ta=TA, hm=HM, ift=IFT, cls)])
}

extract_events_station <- function(d) {
  if (is.null(d) || nrow(d) < 2L) return(NULL)
  v <- d[!is.na(cls)]; if (nrow(v) < 2L) return(NULL)
  setorder(v, tm_min)
  seg <- cumsum(c(TRUE, diff(v$tm_min) > 1L))
  out <- list()
  for (g in unique(seg)) {
    idx <- which(seg == g); s <- v[idx]
    n <- nrow(s); if (n < 2L) next
    is_ice <- s$cls == "ICE"
    if (!any(is_ice)) next
    r <- rle(is_ice); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
    ir <- which(r$values)
    runs <- data.table(st = starts[ir], en = ends[ir])
    merged <- list(); cur <- runs[1]
    if (nrow(runs) >= 2L) for (i in 2:nrow(runs)) {
      if (runs$st[i] - cur$en - 1L < cfg$g_min) cur$en <- runs$en[i]
      else { merged[[length(merged)+1L]] <- cur; cur <- runs[i] }
    }
    merged[[length(merged)+1L]] <- cur; merged <- rbindlist(merged)
    for (i in seq_len(nrow(merged))) {
      st <- merged$st[i]; en <- merged$en[i]
      if (st <= 1L) next
      pre <- s$cls[st-1L]
      if (!(pre %in% c("DRY","WET","SNOW"))) next
      if (!is.na(s$ts[st]) && s$ts[st] > cfg$ice_ts_max) {
        .qc$ice_onset <- .qc$ice_onset + 1L; next }
      dur <- sum(s$cls[st:en] == "ICE")
      if (dur < cfg$p_min) next
      out[[length(out)+1L]] <- data.table(
        stn = s$stn[1], onset_min = s$tm_min[st], pre_state = pre,
        ts0 = s$ts[st], ta0 = s$ta[st], rh0 = s$hm[st],
        ice_dur = as.integer(dur),
        ift_pos = isTRUE(is.finite(s$ift[st]) && s$ift[st] > 0),
        m_liq = s$ts[st] - 0)
    }
  }
  if (!length(out)) return(NULL)
  return(rbindlist(out))
}

build_intl <- function(files) {
  res <- vector("list", length(files))
  for (i in seq_along(files)) {
    res[[i]] <- extract_events_station(read_cleaned_station(files[i]))
    if (i %% 10L == 0L) message("  ", i, "/", length(files))
  }
  ev <- rbindlist(res, fill = TRUE)
  if (nrow(ev) > 0L) {
    ev[, onset_tm := as.POSIXct(onset_min*60, origin="1970-01-01", tz="UTC")]
    ev[, year := data.table::year(onset_tm)]
    setorder(ev, stn, onset_tm)
  }
  return(ev)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  full <- length(args) >= 1L && args[1] == "full"
  files <- list.files(cfg$data_dir, pattern = "\\.csv$", full.names = TRUE)

  EXF <- meta_file("stn_exclusion_intl.csv")
  n0 <- length(files)
  if (file.exists(EXF)) {
    EX <- fread(EXF)
    sof <- as.integer(sub(".*rwis_(\\d+)_.*", "\\1", basename(files)))
    files <- files[!(sof %in% EX$STN_ID)]
    cat(sprintf("[Excluded] %d stations (%s) — %d -> %d\n", n0 - length(files),
                paste(EX$STN_ID, collapse = ", "), n0, length(files)))
  } else {
    stop("Station exclusion list not found: ", EXF)
  }

  use <- if (full) files else head(files, 3L)
  cat("Cleaned stations:", length(files), " / used:", length(use),
      if (full) "(full)" else "(test)", "\n")
  ev <- build_intl(use)
  cat(sprintf("\n[QC] internal consistency (|Ts-Ta|>%g) %d min · spikes %d min · ICE onsets (Ts>%g) excluded: %d\n",
              cfg$tsta_max, .qc$consist, .qc$spike, cfg$ice_ts_max, .qc$ice_onset))
  cat("Events:", nrow(ev), " stations:", length(unique(ev$stn)), "\n")
  if (nrow(ev) > 0L) {
    cat("pre_state:\n"); print(table(ev$pre_state))
    cat("Years:\n"); print(table(ev$year))
    cat("IFT>0:", round(mean(ev$ift_pos), 3), "\n")
    if (full) {
      dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(ev, file.path(cfg$out_dir, "intl_event_table.rds"))
      fwrite(ev, file.path(cfg$out_dir, "intl_event_table.csv"))
      cat("Saved: 02_MART/intl_event_table.{rds,csv}\n")
    }
  }
}
