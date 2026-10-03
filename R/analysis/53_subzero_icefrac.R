# Computes how often the road surface is at or below 0 degC and how much of that time is icy.
# Inputs: 1-min RWIS files; meta/stn_exclusion_intl.csv.
# Output: 02_MART/subzero_icefrac.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
DD <- rwis_dir()
EX <- fread(meta_file("stn_exclusion_intl.csv"))

map_sfs <- function(s) { o <- rep(NA_character_, length(s))
  o[s == 1L] <- "DRY"; o[s %in% c(2L, 11:14)] <- "WET"
  o[s == 3L] <- "SNOW"; o[s %in% c(4L, 15L)] <- "ICE"; o }

fs <- list.files(DD, pattern = "\\.csv$", full.names = TRUE)
sid <- as.integer(sub(".*rwis_(\\d+)_.*", "\\1", basename(fs)))
fs <- fs[!(sid %in% EX$STN_ID)]
cat("n_stations", length(fs), "\n")

acc <- vector("list", length(fs))
for (i in seq_along(fs)) {
  d <- fread(fs[i], select = c("datetime", "mon", "TR1", "SFS"),
             na.strings = c("-999", "-999.0", ""))
  d <- d[mon %in% c(11L, 12L, 1L, 2L, 3L)]; d <- filter_winter(d)
  d <- d[is.finite(TR1) & TR1 >= -40 & TR1 <= 40]
  d[, cls := map_sfs(as.integer(SFS))]
  d <- d[!is.na(cls)]
  if (!nrow(d)) next
  acc[[i]] <- data.table(n = nrow(d), n_sub = sum(d$TR1 <= 0),
                         n_ice = sum(d$cls == "ICE"),
                         n_sub_ice = sum(d$TR1 <= 0 & d$cls == "ICE"))
  if (i %% 20L == 0L) message("  ", i, "/", length(fs))
}
A <- rbindlist(acc)
tot <- A[, lapply(.SD, sum)]
cat(sprintf("\nValid observed minutes  %s\n", format(tot$n, big.mark = ",")))
cat(sprintf("Surface temp ≤0℃        %s (%.2f%%)\n", format(tot$n_sub, big.mark = ","),
            100 * tot$n_sub / tot$n))
cat(sprintf("Icing (ICE)             %s (%.3f%%)\n", format(tot$n_ice, big.mark = ","),
            100 * tot$n_ice / tot$n))
cat(sprintf("Icing among ≤0℃         %s (%.2f%%)  -> %.1f%% of sub-zero time is ice-free\n",
            format(tot$n_sub_ice, big.mark = ","), 100 * tot$n_sub_ice / tot$n_sub,
            100 * (1 - tot$n_sub_ice / tot$n_sub)))
fwrite(tot, file.path(ROOT, "02_MART/subzero_icefrac.csv"), bom = TRUE)
cat("\nSaved: 02_MART/subzero_icefrac.csv\n")
