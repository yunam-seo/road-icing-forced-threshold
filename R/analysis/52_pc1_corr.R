# Correlates the first principal component of onset conditions with dose, precipitation, wind and pressure tendency.
# Inputs: 02_MART/intl_continuum.csv, intl_physical_diag.csv, intl_event_table.rds.
# Output: 02_MART/intl_pc1_corr.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")

C <- fread(file.path(MART, "intl_continuum.csv"))
P <- fread(file.path(MART, "intl_physical_diag.csv"))
J <- merge(C, P[, .(eid, ws_mean, ws_onset, pa_tend6, noct_ts_ta)], by = "eid")
ev <- readRDS(file.path(MART, "intl_event_table.rds")); ev[, eid := paste0(stn, "_", onset_min)]
J <- merge(J, ev[, .(eid, rh0, m_liq)], by = "eid")
cat(sprintf("Continuum %d rows · physical diagnostics %d rows -> joined %d rows\n", nrow(C), nrow(P), nrow(J)))

VARS <- c(cold_dose = "cumulative cold dose", precip_frac = "24 h share of minutes with precipitation",
          ws_mean = "wind speed (mean)", ws_onset = "wind speed (at onset)",
          pa_tend6 = "6h pressure tendency", ts_amp = "surface temperature amplitude",
          onset_hour = "onset hour", noct_ts_ta = "night Ts-Ta",
          rh0 = "relative humidity at onset", m_liq = "liquid freezing margin (surface temperature at onset)")

res <- rbindlist(lapply(names(VARS), function(v) {
  if (!v %in% names(J)) return(NULL)
  x <- J[[v]]; ok <- is.finite(J$PC1) & is.finite(x)
  data.table(variable = v, description = VARS[[v]], n = sum(ok),
             r = round(cor(J$PC1[ok], x[ok]), 3),
             p = signif(suppressWarnings(cor.test(J$PC1[ok], x[ok])$p.value), 3))
}))
res <- res[order(-abs(r))]
print(res)

fwrite(res, file.path(MART, "intl_pc1_corr.csv"), bom = TRUE)
cat("\nSaved: 02_MART/intl_pc1_corr.csv\n")
