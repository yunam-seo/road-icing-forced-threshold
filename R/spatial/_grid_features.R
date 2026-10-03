# Defines the nine surrounding-field features (30-min change, 10 km anomaly, upwind dew-point depression and temperature difference) from extracted grid components.
# Input columns: station-cell, previous-time, 10 km disc-mean and upwind cell values, plus the wind-observation time gap.

FEAT_SPATIAL_FINAL <- c("d_dd","d_td","d_ta","ta_anom10","td_anom10",
                        "dd_up10","dd_up5","dd_up20","dta_up10")
WIND_GAP_MAX_MIN <- 30L

build_final_features <- function(G, keys) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  need <- c("ta_vc","td_vc","ta_vc_prev","td_vc_prev","ta_vreg","td_vreg",
            paste0(rep(c("ta_up","td_up"), each = 3), c(5, 10, 20)), "wind_gap_min")
  miss <- setdiff(need, names(G)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  O <- G[, ..keys]
  O[, d_ta := G$ta_vc - G$ta_vc_prev]
  O[, d_td := G$td_vc - G$td_vc_prev]
  O[, d_dd := (G$ta_vc - G$td_vc) - (G$ta_vc_prev - G$td_vc_prev)]
  O[, ta_anom10 := G$ta_vc - G$ta_vreg]
  O[, td_anom10 := G$td_vc - G$td_vreg]
  masked <- is.na(G$wind_gap_min) | G$wind_gap_min > WIND_GAP_MAX_MIN
  for (k in c(5, 10, 20)) {
    v <- G[[paste0("ta_up", k)]] - G[[paste0("td_up", k)]]; v[masked] <- NA_real_
    O[, (paste0("dd_up", k)) := v] }
  v <- G$ta_up10 - G$ta_vc; v[masked] <- NA_real_
  O[, dta_up10 := v]
  O[, wind_masked := masked]
  O[]
}
