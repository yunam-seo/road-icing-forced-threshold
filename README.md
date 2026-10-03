# Early warning of road icing from current, antecedent and surrounding meteorological conditions — analysis code

Code accompanying the article. It produces every number, figure and table reported there.
Assembling the manuscript files themselves is not part of it.

**Observations.** One-minute road weather information system (RWIS) records from 121 stations on
expressways in the Republic of Korea, three freezing seasons (2023–2026): 48,117,793 station-minutes
and 2,427 icing onsets. A second, independent set of 127 stations with **no station in common**
(3,091 onsets) is used for external replication.

---

## 1. What is here, and what is not

This repository **is** the analysis. It **is not** a redistribution of the observations: the RWIS
records and the 500 m objectively analysed fields are produced by the Korea Meteorological
Administration (KMA), and the scripts read them from paths you configure (section 4).

The derived datasets behind every figure and table — the CSD indicator values of every icing onset
and control time in both station sets, and the aggregated values plotted or tabulated — are
included in `data/` of this repository (`road-icing-derived-data`), with a column dictionary
(`data/README.md`) and a manifest (`data/MANIFEST.csv`) mapping each file to the script here that
produced it. Event-level observations (temperatures, humidity, doses, coordinates) are not part of
that dataset; see its README.

## 2. Layout

```
_config.R, _config.py   path configuration, read once by every script
R/analysis/             event definition, cumulative stress, CSD, nowcast   (44 files)
R/spatial/              spatial extension; needs the objectively analysed fields   (11 files)
R/figures/              the six article figures   (7 files)
docs/REPRODUCE.md       execution order and which stages are heavy
docs/MANIFEST.md        every figure and table mapped to the script that produces it
```

**Reading the file names.** The leading number is the order in which a script was written, not the
order in which it runs — `docs/REPRODUCE.md` gives the run order. Gaps in the numbering are scripts
that were not needed for this article and are therefore not published. Two tags appear in names and
in intermediate file names: `intl_` marks the primary 121-station network, and `icetarget_` the
independent 127-station network. Scripts beginning with `_` are shared definitions, not entry
points.

## 3. Environment

R 4.5.1 with

| package | version | used for |
|---|---|---|
| cluster | 2.1.8.1 | silhouette widths in the k-means check |
| data.table | 1.18.4 | all data handling |
| diptest | 0.77.1 | unimodality (dip) test behind the continuum result |
| lightgbm | 4.6.0 | nowcast and ladder models |
| maps | 3.4.3 | station map in Fig. 1 |
| terra | 1.8.54 | reading the 500 m objectively analysed NetCDF fields (`R/spatial` only) |

Python 3 with no third-party packages is needed only for `R/analysis/162_derived_numbers.py`, which
collects values quoted in the text.

## 4. Configuration

Paths come from environment variables; `_config.R` reads them once and every script sources it.
Only the first two are needed for the core analysis.

| variable | meaning |
|---|---|
| `RICE_ROOT` | output tree (02_MART, 05_REPORT); created if absent |
| `RICE_DATA` | input tree (01_ORG/RWIS_1MIN by default) |
| `RICE_META` | station metadata folder (default: RICE_DATA/00_META if present, else meta/ in this repository) |
| `RICE_RWIS` | folder of the 1-min RWIS files, if it is elsewhere |
| `RICE_VALID_FILE` | 1-min records of the independent network (.rds) |
| `RICE_VALID_META` | station metadata of the independent network (.csv) |
| `RICE_GRID_FIELDS` | folder of the 500 m objectively analysed fields (R/spatial only) |
| `RICE_GRID_LATLON` | grid coordinate table (.csv) |

```bash
export RICE_ROOT=/path/to/output-tree      # 02_MART and 05_REPORT are created here
export RICE_DATA=/path/to/input-tree       # holds 01_ORG/RWIS_1MIN (station metadata is shipped in meta/)
Rscript R/analysis/03_intl_events_cleaned.R full
```

On Windows PowerShell: `$env:RICE_ROOT = "D:\output-tree"`.

If a path is missing, the script stops and names the variable to set; it does not fall back to a
guess.

Two optional variables: `FIGURE_VECTOR=1` makes the figure scripts write vector PDFs (in
`05_REPORT/FIGURES/VECTOR/`) instead of PNGs, and `RICE_CKPT` sets the checkpoint folder of the two
longest spatial extractions (default `02_MART/_ckpt`), from which they resume.

## 5. Two conventions that matter

**Missing values are not interpolated in the main analyses.** Interpolation depresses rolling
variance and lag-1 autocorrelation, which would contaminate the critical-slowing-down (CSD)
diagnosis the article's central negative result rests on. The two analyses treat gaps differently,
as the article states: for event extraction, a break in the 1-min continuity splits the record into
separate continuous segments; in the CSD analysis, missing minutes are left missing inside fixed
analysis windows, and a window is used only if at least 80% of its minutes are observed. The
exception, stated in the article, is the supplementary k-means check of 24-h trajectories
(`05_intl_cluster.R`, `159_text_values.R`), which fills gaps inside an event by linear
interpolation.

**Operating points come from one rule.** Every probability of detection, false-alarm rate and
trigger lead time in the article is produced by `R/analysis/_trigger_rule.R`, which finds the
largest threshold satisfying `POD >= target` by bisection rather than by scanning a quantile grid;
grid selection makes the reported false-alarm rate depend on the grid spacing.
`R/analysis/76_trigger_unify.R` computes all of these numbers together so that they cannot drift
apart.

## 6. Reproduction

`docs/REPRODUCE.md` gives the order and flags the heavy stages. In short: the event table comes
first; the cumulative-stress, CSD and nowcast chains branch from it; the spatial extension needs
the objectively analysed fields; the figures read only derived outputs and are cheap.

## 7. Which script makes which display item

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

`docs/MANIFEST.md` repeats this with the derived data files.

## 8. Data availability

The 1-min RWIS observations and the 500 m objectively analysed fields are produced by the Korea
Meteorological Administration and are available through its open data service; records outside the
standard distribution window are available from KMA on request and cannot be redistributed by the
authors. The derived datasets needed to reproduce the findings are in `data/`.

## 9. Licence and citation

The code is released under the MIT licence (`LICENSE`): use, modify and redistribute freely,
including commercially, provided the copyright notice is kept, and without warranty. The MIT
licence does not require citation; if this code contributes to your work, please cite the article
(`CITATION.cff`).

The derived datasets in `data/` are released under CC BY 4.0. The observations themselves remain
subject to KMA's terms.
