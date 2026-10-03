# -*- coding: utf-8 -*-
# Computes derived values quoted in the text (ranges, ratios, differences) from 02_MART tables using fixed rules.
# Inputs: CSV outputs in 02_MART.
# Output: 02_MART/derived_numbers.csv.

import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../../_config.py"), encoding="utf-8").read())
import csv
import io
import math

MART = os.path.join(ROOT, "02_MART/")
META = os.path.join(META_DIR, "")

def rows(n, base=MART):
    return list(csv.DictReader(io.open(base + n, encoding="utf-8-sig")))

def kv(n, k="item", v="value"):
    return {r[k]: r[v] for r in rows(n)}

def ceil3(x):
    return math.ceil(x * 1000 - 1e-9) / 1000

OUT = []
def put(name, value, rule):
    OUT.append({"item": name, "value": value, "rule": rule})

z = rows("subzero_icefrac.csv")[0]
put("sub-zero surface temperature share pct", 100 * float(z["n_sub"]) / float(z["n"]), "subzero_icefrac n_sub/n")
put("ICE share among sub-zero pct", 100 * float(z["n_sub_ice"]) / float(z["n_sub"]), "n_sub_ice/n_sub")
put("non-ICE share among sub-zero pct", 100 - 100 * float(z["n_sub_ice"]) / float(z["n_sub"]), "100 − above")

LT = {(r["tau_min"], r["configuration"]): float(r["AUC"]) for r in rows("ladder_tau_compare.csv")}
L = {r["stage"]: float(r["AUC"]) for r in rows("model_ladder_full.csv")}
put("ladder gain from cumulative antecedent conditions", LT[("1440", "full")] - LT[("", "plus_dynamics")],
    "ladder_tau_compare, unrounded: full(tau_d 1440) - plus_dynamics")
put("ladder total gain over current", LT[("1440", "full")] - L["current"],
    "full from ladder_tau_compare (unrounded) - current from model_ladder_full")

S1 = [float(r["AUC"]) for r in rows("ews_indicator_suite.csv")]
put("6 indicators AUC min", min(S1), "ews_indicator_suite AUC min"); put("6 indicators AUC max", max(S1), "max")
C9 = rows("ews_sensitivity_9combo.csv")
v9 = [float(r[c]) for r in C9 for c in ("var_auc", "ar1_auc")]
put("9 combinations AUC min", min(v9), "var_auc·ar1_auc min"); put("9 combinations AUC max", max(v9), "max")
WS = rows("ews_winter_strat.csv")
for col, nm in (("auc_variance_trend", "variance"), ("auc_ar1_trend", "AR1")):
    v = [float(r[col]) for r in WS]
    put("winter-stratified %s AUC min" % nm, min(v), "two networks × 3 winters min"); put("winter-stratified %s AUC max" % nm, max(v), "max")

W22 = rows("winter2022_sensitivity.csv")
put("including 2022/23 AUC change upper bound", ceil3(max(abs(float(W22[2]["var_AUC"])), abs(float(W22[2]["ar1_AUC"])))), "absolute change max, rounded up at the third decimal")
SS = rows("control_seed_spread.csv")
spread = max(max(float(r[c]) for r in SS) - min(float(r[c]) for r in SS) for c in ("var_AUC", "ar1_AUC"))
put("control-time resampling AUC spread upper bound", ceil3(spread), "max of max−min across seeds, rounded up at the third decimal")

T = {r["configuration"]: float(r["false_alarms_per_station_day"]) for r in rows("trigger_unified.csv") if r["group"] == "Fig. 6b (with and without surrounding fields)"}
r6 = round(100 * (T["T+6 RWIS"] - T["T+6 RWIS + surrounding fields"]) / T["T+6 RWIS"])
r3 = round(100 * (T["T+3 RWIS"] - T["T+3 RWIS + surrounding fields"]) / T["T+3 RWIS"])
put("false alarm reduction T+6 pct", r6, "Fig. 6b (RWIS − RWIS + surrounding fields)/RWIS, rounded"); put("false alarm reduction T+3 pct", r3, "")
FB = rows("falsealarm_by_ts.csv")
for h in ("3", "6"):
    v = [float(r["reduction_pct"]) for r in FB if r["H_h"] in (h, h + ".0") and r["ts_bin"] in ("Ts<=-2", "-2<Ts<=0")]
    put("cold-surface false alarm reduction T+%s min" % h, round(min(v)), "Ts≤0 two bins min"); put("cold-surface false alarm reduction T+%s max" % h, round(max(v)), "max")

HO = {(r["holdout"], r["model"]): r for r in rows("grid_holdout_auc.csv")}
put("event-control number of control times", int(HO[("stn", "RWIS")]["n"]) - int(HO[("stn", "RWIS")]["n_events"]), "n − n_events")
dk = [round(float(HO[(s, "initial_3km")]["AUC"]), 3) - round(float(HO[(s, "concurrent surrounding fields")]["AUC"]), 3) for s in ("stn", "day", "winter")]
put("AUC difference vs initial 3 km min", round(min(dk), 3), "difference of values rounded at the third decimal"); put("AUC difference vs initial 3 km max", round(max(dk), 3), "")

EX = rows("stn_exclusion_intl.csv", META)
dist = lambda r: float(r["distance_to_route_m"] if "distance_to_route_m" in r else r["dist_m"])
for i, r in enumerate(sorted(EX, key=lambda r: -dist(r)), 1):
    put("excluded station distance %d" % i, dist(r), "stn_exclusion_intl distance_to_route_m (m)")
    put("excluded station distance km %d" % i, dist(r) / 1000, "distance_to_route_m / 1000")
IS = {r["network"]: float(r["ice_share_classified_pct"]) for r in rows("ice_share_networks.csv")}
put("ICE share ratio independent/primary", round(IS["independent"], 2) / round(IS["primary"], 2), "ratio of rounded values")
nc = [int(r["n_control"]) for r in rows("ews_indicator_suite.csv")]
put("6 indicators number of controls min", min(nc), ""); put("6 indicators number of controls max", max(nc), "")

for s in ("stn", "day", "winter"):
    for m in ("concurrent surrounding fields", "surrounding fields 1 h earlier", "surrounding fields 3 h earlier"):
        put("surrounding-field gain %s %s" % (s, m),
            round(round(float(HO[(s, m)]["AUC"]), 3) - round(float(HO[(s, "RWIS")]["AUC"]), 3), 3), "difference of values rounded at the third decimal")

QS = [r for r in rows("qc_threshold_sens.csv") if r["tolerance_c"] not in ("", "NA")]
for col, nm in (("var_AUC", "variance"), ("ar1_AUC", "AR1")):
    v = [round(float(r[col]), 3) for r in QS]
    put("QC sensitivity %s AUC min" % nm, min(v), "30·15·10 ℃ rounded min"); put("QC sensitivity %s AUC max" % nm, max(v), "max")
B30 = [r for r in QS if float(r["tolerance_c"]) == 30][0]
put("independent network onsets without τ", int(B30["onset"]) - int(B30["n_onset_tau_computed"]), "30 ℃: onsets − onsets with τK (analysis-window coverage < 0.8)")

LT = rows("ladder_tau_compare.csv")
for part, nm in (("dose_only", "cumulative indicators only"), ("full", "full")):
    v = [round(float(r["AUC"]), 3) for r in LT if r["configuration"] == part and r["tau_min"] in ("20", "60")]
    put("short τd %s AUC min" % nm, min(v), "τd 20·60 min rounded min"); put("short τd %s AUC max" % nm, max(v), "max")
    v2 = [round(float(r["AUC"]), 3) for r in LT if r["configuration"] == part]
    put("τd full range %s AUC min" % nm, min(v2), "τd 20·60·1440 min rounded min"); put("τd full range %s AUC max" % nm, max(v2), "max")

SC = rows("spatial_corr_scale.csv")
put("spatial correlation min r floor3", math.floor(min(float(r["r"]) for r in SC) * 1000 + 1e-9) / 1000, "spatial_corr_scale r min, rounded down at the third decimal")
put("spatial correlation max distance km", max(float(r["km"]) for r in SC), "spatial_corr_scale km max")

with io.open(MART + "derived_numbers.csv", "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["item", "value", "rule"]); w.writeheader(); w.writerows(OUT)
for r in OUT:
    print("%-28s %s" % (r["item"], r["value"]))
print("Saved: 02_MART/derived_numbers.csv (%d items)" % len(OUT))

