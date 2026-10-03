# Summarises events per network: counts, stations, prior states, IFT agreement and events per winter.
# Inputs: 02_MART/intl_event_table.rds, icetarget_event_table.rds.
# Outputs: 02_MART/table1_event_summary.csv, table2_network_summary.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")
source(code_file("00_analysis_window.R"))
ev <- as.data.table(readRDS(file.path(MART, "intl_event_table.rds")))
p <- as.POSIXct(ev$onset_min * 60, origin = "1970-01-01", tz = "UTC")
ev[, winter := winter_label(as.integer(format(p, "%Y")), as.integer(format(p, "%m")))]
ev[, winter_lab := sprintf("%d/%02d", winter, (winter + 1L) %% 100L)]
O <- rbindlist(list(
  data.table(item = "n_events", value = nrow(ev)),
  data.table(item = "stations with events", value = uniqueN(ev$stn)),
  data.table(item = paste0("prior state ", c("WET", "SNOW", "DRY")),
             value = as.numeric(sapply(c("WET", "SNOW", "DRY"), function(s) sum(ev$pre_state == s)))),
  data.table(item = "IFT>0 share (%)", value = round(100 * mean(ev$ift_pos > 0, na.rm = TRUE), 1)),
  ev[, .(value = as.numeric(.N)), by = .(winter, winter_lab)][order(winter), .(item = paste0("winter ", winter_lab), value)]))
print(O)
cat("Calendar years (reference):", paste(names(table(format(p, "%Y"))), table(format(p, "%Y")), collapse = " "), "\n")
stopifnot(sum(O[grepl("^winter", item)]$value) == nrow(ev))
fwrite(O, file.path(MART, "table1_event_summary.csv"), bom = TRUE)
ind <- as.data.table(readRDS(file.path(MART, "icetarget_event_table.rds")))
T2 <- rbindlist(lapply(list(primary = ev, independent = ind), function(x) data.table(n_events = nrow(x), n_stations = uniqueN(x$stn),
                WET_pct = round(100 * mean(x$pre_state == "WET"), 1))), idcol = "network")
print(T2); fwrite(T2, file.path(MART, "table2_network_summary.csv"), bom = TRUE)
cat("Saved: 02_MART/table2_network_summary.csv
")
cat("Saved: 02_MART/table1_event_summary.csv\n")
