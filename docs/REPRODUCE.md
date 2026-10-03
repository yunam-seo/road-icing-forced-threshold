# Reproduction

Run order, what each stage produces and which stages are heavy. The stages are listed in the
order in which they must be run. Times are measured on one desktop machine reading the
observations from a network drive; treat them as relative.

Set the paths first (README section 4). `RICE_ROOT` (output tree, created if absent) and
`RICE_DATA` (the 1-min RWIS records) run the primary-network stages; the station metadata the
scripts read is shipped in `meta/`. The independent-network stages (`50_`, `73_`, `153_`, `166_`)
also need `RICE_VALID_FILE` and `RICE_VALID_META`, and the spatial extension needs
`RICE_GRID_FIELDS` and `RICE_GRID_LATLON`.

| stage | what it does | scripts | measured |
|---|---|---|---|
| Event definition | Icing onsets and the tables every later stage reads. | `03_intl_events_cleaned.R`, `04_intl_lagged_build.R`, `05_intl_cluster.R`, `13_intl_physical_diag.R` | 16 min |
| Cumulative stress and CSD | Doses, freeze-thaw accumulation, CSD indicators and their sensitivity. | `07_intl_cumstress.R`, `09_intl_ews.R`, `11_intl_ews_sensitivity.R` | 61 min |
| Nowcast | Evaluation table, operational metrics, calibration and the information ladder. | `12_intl_earlywarning.R`, `19_intl_ew_operational.R`, `20_intl_ew_calibration.R`, `22_intl_ew_baseline.R`, `23_intl_onset_continuity.R` | 20 min |
| Independent network | The same event definition and CSD procedure on the 127-station set. | `50_icetarget_ews.R`, `51_ews_winter_strat.R`, `52_pc1_corr.R`, `53_subzero_icefrac.R` | 19 min |
| Ladder and controls | Indicator suite, alternative state variables, lead dependence, positive control. | `55_model_ladder_metrics.R`, `56_ews_indicator_suite.R`, `57_ews_state_variables.R`, `58_frc_variance_lead.R`, `71_ews_positive_control.R` | 83 min |
| Bootstrap and external application | Station-cluster confidence intervals; the nowcast applied without retraining. | `72_nowcast_bootstrap_ci.R`, `73_nowcast_independent127.R` | 17 min |
| Spatial extension | Needs RICE_GRID_FIELDS and RICE_GRID_LATLON. Reads the 500 m fields; the heaviest stage. | `30_road_grid500_match.R`, `36_grid_full_points.R`, `37_grid_full_extract.R`, `38_grid_full_model.R`, `47_grid_upwind_extract.R`, `45_spatial_corr_scale.R`, `143_grid_ew_scale_extract.R`, `41_intl_ew_pattern_compare.R`, `144_grid_past_scale_extract.R`, `145_grid_holdout_final.R` | 158 min |
| Operational triggers | One trigger rule produces every operating point reported in the article. | `76_trigger_unify.R`, `105_bootstrap_trigger_rows.R` | 1 min |
| Robustness and exported values | Model adequacy, winter and seed sensitivity, QC thresholds, values quoted in the text. | `90_lambda_bandwidth_sens.R`, `98_ar1_model_adequacy.R`, `119_table2_ci.R`, `128_winter2022_sensitivity.R`, `129_grid_auc_export.R`, `130_data_availability.R`, `131_control_seed_spread.R`, `146_leaky_approx_check.R`, `148_csd_wilcoxon_export.R`, `150_event_summary_table1.R`, `151_falsealarm_by_ts.R`, `153_ice_share_networks.R`, `155_tau_sweep.R`, `159_text_values.R`, `166_qc_threshold_sens_export.R`, `169_ladder_tau_compare.R`, `162_derived_numbers.py` | 153 min |
| Figures | Read derived outputs only; cheap. | `64_fig1_overview.R`, `65_fig2_conditions.R`, `66_fig3_cumulative.R`, `62_fig4_csd.R`, `63_fig5_ladder.R`, `67_fig6_falsealarm.R` | 1 min |

## Notes

**Start here.** `R/analysis/03_intl_events_cleaned.R full` builds the event table; nothing else
runs without it. The argument `full` processes every station (omit it for a short test on a few
stations).

**The spatial stage is optional, with a known cost.** The 500 m objectively analysed fields are
a separate, much larger archive. Without them, skip `R/spatial/`; the CSD analysis, the cumulative
indicators, the information ladder (Figs. 2–4 and 5a–c, Tables 1–2, Supplementary Tables 1–5) reproduce
from the RWIS records alone. The operating-point results share one trigger store that includes the
surrounding-field model, so Fig. 5d, Fig. 6, Tables 3–4 and Supplementary Tables 6–7 need the
spatial stage.

**Order within a stage matters in two places.** `76_trigger_unify.R` must run after the nowcast
probabilities exist, because every operating point in the article comes from it; and
`41_intl_ew_pattern_compare.R` must run before `129_grid_auc_export.R` and `76_trigger_unify.R`,
which read the models it stores. Running them in the wrong order produces numbers that look
plausible but mix two versions of the surrounding-field model.

**Missing values are not interpolated** and **operating points come from one rule** — see README
section 5 before changing anything in the chain.

## Full order

 1. `03_intl_events_cleaned.R`
 2. `04_intl_lagged_build.R`
 3. `05_intl_cluster.R`
 4. `13_intl_physical_diag.R`
 5. `07_intl_cumstress.R`
 6. `09_intl_ews.R`
 7. `11_intl_ews_sensitivity.R`
 8. `12_intl_earlywarning.R`
 9. `19_intl_ew_operational.R`
10. `20_intl_ew_calibration.R`
11. `22_intl_ew_baseline.R`
12. `23_intl_onset_continuity.R`
13. `50_icetarget_ews.R`
14. `51_ews_winter_strat.R`
15. `52_pc1_corr.R`
16. `53_subzero_icefrac.R`
17. `55_model_ladder_metrics.R`
18. `56_ews_indicator_suite.R`
19. `57_ews_state_variables.R`
20. `58_frc_variance_lead.R`
21. `71_ews_positive_control.R`
22. `72_nowcast_bootstrap_ci.R`
23. `73_nowcast_independent127.R`
24. `30_road_grid500_match.R`
25. `36_grid_full_points.R`
26. `37_grid_full_extract.R`
27. `38_grid_full_model.R`
28. `47_grid_upwind_extract.R`
29. `45_spatial_corr_scale.R`
30. `143_grid_ew_scale_extract.R`
31. `41_intl_ew_pattern_compare.R`
32. `144_grid_past_scale_extract.R`
33. `145_grid_holdout_final.R`
34. `76_trigger_unify.R`
35. `90_lambda_bandwidth_sens.R`
36. `98_ar1_model_adequacy.R`
37. `105_bootstrap_trigger_rows.R`
38. `119_table2_ci.R`
39. `128_winter2022_sensitivity.R`
40. `129_grid_auc_export.R`
41. `130_data_availability.R`
42. `131_control_seed_spread.R`
43. `146_leaky_approx_check.R`
44. `148_csd_wilcoxon_export.R`
45. `150_event_summary_table1.R`
46. `151_falsealarm_by_ts.R`
47. `153_ice_share_networks.R`
48. `155_tau_sweep.R`
49. `159_text_values.R`
50. `166_qc_threshold_sens_export.R`
51. `169_ladder_tau_compare.R`
52. `162_derived_numbers.py`
53. `64_fig1_overview.R`
54. `65_fig2_conditions.R`
55. `66_fig3_cumulative.R`
56. `62_fig4_csd.R`
57. `63_fig5_ladder.R`
58. `67_fig6_falsealarm.R`
