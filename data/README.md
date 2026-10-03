# Road icing and early warning signals — derived data

Aggregated data behind every figure and table of the manuscript, and the critical-slowing-down
(CSD) indicator values of all icing-onset events and control times in both station networks.
The analysis code is in this repository (see the top-level README); `MANIFEST.csv` maps every
file here to the display item it serves and to the script that produced it.

## What is included

| Folder | Contents |
|---|---|
| `display_items/` | The values plotted in Figs. 1–6 and printed in Tables 1–4 and Supplementary Tables 1–7. Distribution panels are given as aggregates (bin counts, quartiles, decile means, box-plot statistics). |
| `event_level/` | CSD indicator values (Kendall's τ of rolling variance and lag-1 autocorrelation) for every icing onset and control time, primary and independent networks. These are statistics computed from the observations, not observations. |
| `MANIFEST.csv` | One row per file: rows, columns, display items, source file in the analysis tree, producing script, contents |

Every table cell in the manuscript (Tables 1–4, Supplementary Tables 1–7; 405 cells) was checked
against these files with no mismatch. Supplementary Table 7 is `display_items/trigger_operating_points.csv`.

## How these files relate to the released code

Each file is exported from one output of the analysis tree; `MANIFEST.csv` gives that output
(`source_file`) and the script that produces it (`produced_by`). The released code reproduces
those outputs from the RWIS records. The export itself — selecting the rows and columns used in the
article and renaming the columns to English — is a separate step that is not part of the released
code; the mapping is fully described by the manifest and the column dictionary below.

## What is not included, and why

The 1-min road weather information system (RWIS) observations are produced by the Korea
Meteorological Administration (KMA). Records from 15 December 2025 onwards are publicly available
through the KMA API Hub (https://apihub.kma.go.kr); earlier records were provided for this study by
the National Institute of Meteorological Sciences of KMA and are available from KMA on request.
They are not redistributed here, and neither are event-level values that reproduce them
(road surface temperature, air temperature and humidity at or before onset, cumulative doses,
station coordinates). These event-level data are available from the corresponding author on
reasonable request, subject to the KMA conditions of use.

The station map (Fig. 1a) and the analysis diagram (Fig. 1c) have no separate data file.
Fig. 2a is shown as smoothed densities; the data file gives the counts of the same values in
0.5 °C bins.

## Conventions

- **Time.** `window_end_min` and the onset minute inside `event_id` count minutes since
  1970-01-01 00:00 on the station clock, Korea Standard Time (KST, UTC+9), in both event-level
  files. The analysis stored the independent network's minutes as UTC epoch minutes; they were
  shifted by +540 min at export so that the two files share this convention. To recover a
  timestamp, treat the value as minutes since 1970-01-01 00:00 UTC and read the result as KST
  (for example, in R: `as.POSIXct(m * 60, origin = "1970-01-01", tz = "UTC")`).
- **Event identifier.** `event_id` = `<station_id>_<onset minute>`; control times have no event identifier.
- **Winter.** November–March; `2023/24` denotes November 2023 to March 2024.
- **Operating points.** Every probability of detection (POD), false-alarm count and trigger lead
  comes from one rule: the largest probability threshold with POD ≥ target.

## Column dictionary

Columns are listed file by file. AUC is the area under the receiver operating characteristic
curve for discriminating icing onsets from control times (or, for prediction, events from
non-events); τ is Kendall's rank correlation of an indicator with time within the analysis
window; `ci_low`, `ci_high` bound the 95% interval; `n…` columns are counts.

### event_level/csd_indicators_primary.csv
event_id (`<station_id>_<onset minute>`; empty for control times); type (`onset`: window ending 1 min before onset, or `control`); station_id (KMA RWIS station identifier); window_end_min (end of the analysis window, minutes since 1970-01-01 00:00 KST); variance_tau, ar1_tau (τ of rolling variance and lag-1 autocorrelation, AR1); variance_window_end, ar1_window_end (rolling variance and AR1 at the end of the window)

### event_level/csd_indicators_independent.csv
Same columns as csd_indicators_primary.csv. The analysis stored this network's minutes as UTC epoch minutes; they were shifted by +540 min here so that window_end_min and the onset minute inside event_id follow the same KST-clock convention as the primary file

### Figures
| File | Columns |
|---|---|
| fig1b_events_per_winter | network (`primary`, `independent`); winter; n_events |
| fig2a_onset_ts_histogram | prior_state (DRY, WET, SNOW); bin_low_c, bin_high_c (road surface temperature bin, °C); n_events |
| fig2b_eta2 | variable; eta2 (share of variation explained by prior surface state, η²) |
| fig2c_trajectories | prior_state; minutes_before_onset; median_c, q1_c, q3_c (median and quartiles of road surface temperature, °C); n (events) |
| fig2d_pc1_histogram | bin_low, bin_high (first principal component, PC1); n_events |
| fig2e_pc1_deciles | pc_decile; pc1_median; precipitation_fraction_24h (share of minutes with precipitation in the preceding 24 h); rh_onset_pct (relative humidity at onset, %); ts_onset_c (road surface temperature at onset, °C); log1p_cold_dose (ln(1 + cold dose)); n_events. The figure rescales each variable to 0–1 across the decile means |
| fig3abc_box_statistics | indicator (cold dose, log10(1 + deposition dose), freeze–thaw accumulation); group (`control`, `at_onset`: the indicator evaluated at the onset minute); whisker_low, q1, median, q3, whisker_high (box-plot statistics, whiskers at 1.5 IQR); n |
| fig3_tau_sweep | evaluation_time (`t0-L`: cold dose evaluated L min before onset); tau_min (decay time constant τd, min); AUC, ci_low, ci_high; diff_from_best, diff_ci_low, diff_ci_high (AUC difference from the best τd, paired station bootstrap, B = 500); auc_instant_ts, auc_instant_subzero (AUC of the instantaneous road surface temperature and of max(0, −Ts) at the same time); n |
| fig3_ladder_auc | stage; description; AUC (station-held-out CV, T+3 h) |
| fig4_csd_sensitivity | W_h (analysis window, h); bw (detrending bandwidth, min); var_med, ar1_med (median τ of rolling variance and AR1 in onset windows); var_auc, ar1_auc; n_stations |
| fig5_ladder_metrics | model; ROC_AUC; PR_AUC (area under the precision–recall curve); BSS (Brier skill score against the climatological base rate); POD; FA_per_day (false alarms per station-day); FAR_point, CSI_point (point-wise false-alarm ratio and critical success index); lead_med_min (median first-trigger lead, min) |
| fig5_base_rate | item; value (T+3 h positive rate); n_evaluation_points |
| fig6_holdout_auc | holdout (`stn`, `day`, `winter` held out); model (RWIS only, or with concurrent surrounding fields or fields 1 h or 3 h earlier); AUC; n (evaluation points); n_events; n_folds; gain_vs_rwis_only |
| fig6_auc_by_lead | lead_min, lead_h (forecast horizon); auc_rwis_only; auc_with_surrounding_fields; difference |

### Tables
| File | Columns |
|---|---|
| table1_event_summary | item; value |
| table2_network_summary | network; n_events; n_stations (stations with events); WET_pct (share of events with prior state WET, %) |
| table2_csd_auc_ci | network; indicator; n_onset; n_control; AUC; ci_low, ci_high (station-cluster bootstrap, B = 5,000) |
| table2_wilcoxon | network; indicator (`var_tau`, `ar1_tau`); n; V, p (one-sided Wilcoxon signed-rank test of τ > 0) |
| table3_auc_and_trigger_ci | indicator (AUC by forecast horizon, or an alert metric); value; ci_low, ci_high (station-cluster bootstrap, B = 1,000) |
| trigger_operating_points | group (display item); configuration; lead_h (forecast horizon, h); persistence_k (consecutive evaluation times required); target_pod; achieved_pod; false_alarms_per_station_day; first_trigger_lead_median_min, first_trigger_lead_q1_min, first_trigger_lead_q3_min; n_onsets; probability_threshold (an alert is issued when the score is at or above it; for the Model 0 row the score is the negative road surface temperature, so −1.2 means an alert when Ts ≤ 1.2 °C) |

### Supplementary tables
| File | Columns |
|---|---|
| supplementary_table1_csd_indicator_suite | control_design (`random`, `hour_matched`); indicator; code; onset_median, control_median (median τ); AUC; frac_positive_tau_onset (share of onset windows with τ > 0); n_onset; n_control |
| supplementary_table2_positive_control | setting (gaps imposed); lam1 (recovery rate at the end of the synthetic series); window_end_ar1; ar1_rise (AR1 increase over the window); auc_rolling_variance, auc_lag1_autocorrelation; frac_positive_tau_variance; n_valid_windows. The positive control is a linear stochastic process whose recovery rate declines, not a simulated fold transition |
| supplementary_table3_state_variables | state_variable; indicator; AUC; frac_positive_tau_onset; n_valid_windows_onset; n_valid_windows_control |
| supplementary_table4_lead_dependence | variable; indicator; lead_min (window end before onset, min); AUC; n_valid_windows_onset |
| supplementary_table5_winter2022_sensitivity | analysis; winter_window (first winter included, by the year of its November; empty in the difference row); n_events; var_AUC, ar1_AUC |
| supplementary_table6_independent_application | lead; station_set; AUC; positive_rate; POD; false_alarms_per_station_day; first_trigger_lead_median_min |

## Licence

Released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0).
Please cite the accompanying article.
