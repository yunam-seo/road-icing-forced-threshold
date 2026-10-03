# Manifest — derived data, display items and the scripts that produce them

Each row is one file in `data/` (the derived dataset `road-icing-derived-data`): which display
item it backs, which script here produces it, and what it contains. The same information is in
`data/MANIFEST.csv`.

| file | display item | produced by | contents |
|---|---|---|---|
| `event_level/csd_indicators_primary.csv` | Fig. 4a,b; Table 2 | `09_intl_ews.R` | Kendall tau of rolling variance and lag-1 autocorrelation before onsets and at control times, primary network |
| `event_level/csd_indicators_independent.csv` | Table 2 | `50_icetarget_ews.R` | Same indicators, independent network (minutes shifted by +540 to the KST-clock convention of the primary file) |
| `display_items/fig1b_events_per_winter.csv` | Fig. 1b | `64_fig1_overview.R` | Icing-onset events per winter and network |
| `display_items/fig2a_onset_ts_histogram.csv` | Fig. 2a | `65_fig2_conditions.R` | Counts of road surface temperature at onset in 0.5 degC bins by prior surface state |
| `display_items/fig2b_eta2.csv` | Fig. 2b | `65_fig2_conditions.R` | Proportion of variation in conditions at icing onset explained by prior surface state |
| `display_items/fig2c_trajectories.csv` | Fig. 2c | `65_fig2_conditions.R` | Median and interquartile range of road surface temperature before onset by prior surface state |
| `display_items/fig2d_pc1_histogram.csv` | Fig. 2d | `65_fig2_conditions.R` | Histogram of PC1 of the onset conditions |
| `display_items/fig2e_pc1_deciles.csv` | Fig. 2e | `65_fig2_conditions.R` | Mean physical characteristics in PC1 deciles (the figure shows them rescaled to 0-1) |
| `display_items/fig3abc_box_statistics.csv` | Fig. 3a-c | `66_fig3_cumulative.R` | Box-plot statistics (whiskers 1.5 IQR, outliers not shown) of the cumulative antecedent-condition indicators |
| `display_items/table1_event_summary.csv` | Table 1 | `150_event_summary_table1.R` | Event counts, prior states, IFT agreement, events per winter |
| `display_items/table2_network_summary.csv` | Table 2 | `150_event_summary_table1.R` | Events, stations and WET share per network |
| `display_items/table2_csd_auc_ci.csv` | Table 2; Fig. 4d | `119_table2_ci.R` | CSD discrimination AUC with station-cluster bootstrap CI (B = 5,000) |
| `display_items/table2_wilcoxon.csv` | Table 2 | `148_csd_wilcoxon_export.R` | One-sided Wilcoxon signed-rank tests of tau > 0 |
| `display_items/table3_auc_and_trigger_ci.csv` | Table 3 | `72_nowcast_bootstrap_ci.R; 105_bootstrap_trigger_rows.R` | Station-held-out CV AUC by lead and trigger metrics with bootstrap CI (B = 1,000) |
| `display_items/trigger_operating_points.csv` | Table 4; Supplementary Table 7; Fig. 5d; Fig. 6b | `76_trigger_unify.R` | Every operating point in the paper, one threshold rule |
| `display_items/supplementary_table1_csd_indicator_suite.csv` | Supplementary Table 1 | `56_ews_indicator_suite.R` | Six indicators x two control designs |
| `display_items/supplementary_table2_positive_control.csv` | Supplementary Table 2 | `71_ews_positive_control.R` | Positive control: linear stochastic process with a declining recovery rate |
| `display_items/supplementary_table3_state_variables.csv` | Supplementary Table 3 | `57_ews_state_variables.R` | Alternative state variables |
| `display_items/supplementary_table4_lead_dependence.csv` | Supplementary Table 4 | `58_frc_variance_lead.R` | Indicator trends vs window-end lead |
| `display_items/supplementary_table5_winter2022_sensitivity.csv` | Supplementary Table 5 | `128_winter2022_sensitivity.R` | Including the incomplete 2022/23 winter |
| `display_items/supplementary_table6_independent_application.csv` | Supplementary Table 6 | `73_nowcast_independent127.R` | Nowcast applied unchanged to the independent network |
| `display_items/fig3_tau_sweep.csv` | Fig. 3d | `155_tau_sweep.R` | Cold-dose AUC across decay time constants with paired bootstrap differences (B = 500) |
| `display_items/fig3_ladder_auc.csv` | Fig. 3 | `22_intl_ew_baseline.R` | AUC of threshold, dose-only, current, dynamics and full models |
| `display_items/fig4_csd_sensitivity.csv` | Fig. 4c | `11_intl_ews_sensitivity.R` | CSD result across window length x detrending bandwidth |
| `display_items/fig5_ladder_metrics.csv` | Fig. 5 | `55_model_ladder_metrics.R` | ROC-AUC, PR-AUC, BSS and trigger metrics at equal POD |
| `display_items/fig5_base_rate.csv` | Fig. 5b | `63_fig5_ladder.R` | Climatological base rate of the T+3 h target |
| `display_items/fig6_holdout_auc.csv` | Fig. 6c | `145_grid_holdout_final.R` | Event-control AUC with station/date/winter held out, RWIS only vs surrounding fields |
| `display_items/fig6_auc_by_lead.csv` | Fig. 6a | `129_grid_auc_export.R` | AUC by lead with and without surrounding fields |

## Display item to script

| display item | script |
|---|---|
| Fig. 1 | `64_fig1_overview.R` |
| Fig. 2 | `65_fig2_conditions.R` |
| Fig. 3 | `155_tau_sweep.R`, `22_intl_ew_baseline.R`, `66_fig3_cumulative.R` |
| Fig. 4 | `09_intl_ews.R`, `119_table2_ci.R`, `11_intl_ews_sensitivity.R`, `62_fig4_csd.R` |
| Fig. 5 | `55_model_ladder_metrics.R`, `63_fig5_ladder.R`, `76_trigger_unify.R` |
| Fig. 6 | `129_grid_auc_export.R`, `145_grid_holdout_final.R`, `67_fig6_falsealarm.R`, `76_trigger_unify.R` |
| Table 1 | `150_event_summary_table1.R` |
| Table 2 | `09_intl_ews.R`, `119_table2_ci.R`, `148_csd_wilcoxon_export.R`, `150_event_summary_table1.R`, `50_icetarget_ews.R` |
| Table 3 | `105_bootstrap_trigger_rows.R`, `72_nowcast_bootstrap_ci.R` |
| Table 4 | `76_trigger_unify.R` |
| Supplementary Table 1 | `56_ews_indicator_suite.R` |
| Supplementary Table 2 | `71_ews_positive_control.R` |
| Supplementary Table 3 | `57_ews_state_variables.R` |
| Supplementary Table 4 | `58_frc_variance_lead.R` |
| Supplementary Table 5 | `128_winter2022_sensitivity.R` |
| Supplementary Table 6 | `73_nowcast_independent127.R` |
| Supplementary Table 7 | `76_trigger_unify.R` |
