# Computes missing-data rates per variable and QC removal counts within the analysis window.
# Inputs: 1-min RWIS files; meta/stn_exclusion_intl.csv.
# Output: 02_MART/data_availability.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
EXF <- meta_file("stn_exclusion_intl.csv")
source(code_file("00_analysis_window.R"))

DDIR <- rwis_dir()
FREEZE <- c(11L, 12L, 1L, 2L, 3L)
ROAD <- c(surface_temperature = "TR1", surface_state = "SFS", ice_film_thickness = "IFT")
MET <- c(air_temperature = "TA", relative_humidity = "HM", precipitation_occurrence = "RN05_YN", wind_speed = "WS", wind_direction = "WD", pressure = "PA")

files <- list.files(DDIR, pattern = "[.]csv$", full.names = TRUE)
EX <- fread(EXF)
sof <- as.integer(sub(".*rwis_([0-9]+)_.*", "\\1", basename(files)))
files <- files[!(sof %in% EX$STN_ID)]
cat(sprintf("Stations %d (excluded %d)\n", length(files), nrow(EX)))

cols <- c("datetime", "mon", unname(ROAD), unname(MET))
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- fread(files[i], select = cols, na.strings = c("-999", "-999.0", ""))
  d <- d[mon %in% FREEZE]
  d <- filter_winter(d)
  if (!nrow(d)) next
  tm <- as.POSIXct(strptime(sprintf("%012.0f", as.numeric(d$datetime)),
                            "%Y%m%d%H%M", tz = "UTC"))
  d <- d[!is.na(tm)]
  if (!nrow(d)) next
  d[, tm_min := as.integer(tm) %/% 60L]
  setorder(d, tm_min)
  r <- as.list(sapply(c(ROAD, MET), function(cc) mean(is.na(d[[cc]]))))
  ts <- d$TR1; ta <- d$TA
  dmax_raw <- suppressWarnings(max(abs(ts - ta), na.rm = TRUE))
  ts[ts < -40 | ts > 40] <- NA_real_; ta[ta < -40 | ta > 40] <- NA_real_
  n_consist <- sum(!is.na(ts) & !is.na(ta) & abs(ts - ta) > 30)
  ts[!is.na(ts) & !is.na(ta) & abs(ts - ta) > 30] <- NA_real_
  dp <- ts - shift(ts); dn <- shift(ts, type = "lead") - ts
  n_spike <- sum(!is.na(dp) & !is.na(dn) & abs(dp) > 10 & abs(dn) > 10 &
                 sign(dp) != sign(dn))
  rows[[i]] <- c(list(stn = sof[!(sof %in% EX$STN_ID)][i], n_obs_minutes = nrow(d), max_abs_ts_minus_ta_raw = dmax_raw,
                      n_qc_consistency = n_consist, n_qc_spike = n_spike), r)
  if (i %% 20L == 0L) message("  ", i, "/", length(files))
}
D <- rbindlist(rows[!vapply(rows, is.null, TRUE)], fill = TRUE)

tot <- sum(D$n_obs_minutes)
out <- data.table(item = "total observed minutes", value = tot, unit = "min")
out <- rbind(out, data.table(item = "observed minutes, median per station",
                             value = median(D$n_obs_minutes), unit = "min"))
out <- rbind(out, data.table(item = "QC removed: surface–air temperature difference",
                             value = sum(D$n_qc_consistency), unit = "min"))
out <- rbind(out, data.table(item = "QC removed: spikes",
                             value = sum(D$n_qc_spike), unit = "min"))
out <- rbind(out, data.table(item = "max |Ts-Ta| before QC", value = round(max(D$max_abs_ts_minus_ta_raw[is.finite(D$max_abs_ts_minus_ta_raw)]), 1), unit = "℃"))
out <- rbind(out, data.table(item = "QC removed share (%)", value = signif(100 * (sum(D$n_qc_consistency) + sum(D$n_qc_spike)) / tot, 3), unit = "%"))
out <- rbind(out, data.table(item = "QC removed, total",
                             value = sum(D$n_qc_consistency) + sum(D$n_qc_spike),
                             unit = "min"))
for (nm in names(ROAD)) {
  out <- rbind(out, data.table(item = sprintf("missing rate %s", nm),
                               value = round(100 * mean(D[[nm]]), 2), unit = "%"))
}
met_mean <- mean(vapply(names(MET), function(cc) mean(D[[cc]]), 0))
out <- rbind(out, data.table(item = sprintf("missing rate, mean of %d meteorological variables", length(MET)),
                             value = round(100 * met_mean, 2), unit = "%"))
for (nm in names(MET)) {
  out <- rbind(out, data.table(item = sprintf("missing rate %s", nm),
                               value = round(100 * mean(D[[nm]]), 2), unit = "%"))
}
print(out)
fwrite(out, file.path(MART, "data_availability.csv"), bom = TRUE)
cat(sprintf("\nStations %d · total observed minutes %s\nSaved: 02_MART/data_availability.csv\n",
            nrow(D), format(tot, big.mark = ",")))
