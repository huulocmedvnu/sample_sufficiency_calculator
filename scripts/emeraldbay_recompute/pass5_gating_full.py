#!/usr/bin/env python
"""Full-atlas (52-line) downsample-and-measure gating verification for EmeraldBay.

The regime-gating CLASSIFICATION (OVER vs UNDER) was already run over all 52 lines from the
per-condition moments (3,971 groups, 217 predicted OVER). The empirical downsample-and-measure CHECK
of those predictions, however, was previously restricted to the five lines shared with Tahoe-100M
(only the shared-line per-cell coordinates had been retained). Now that pass2b re-projected the FULL
atlas (all_cells.npz, 52 lines, 1,831,648 cells), this verifies every predicted-OVER group across the
entire atlas: classify each (sample x line) group from its PCA-diagonal plug-in, and for each one
predicted OVER, downsample to n* and check the realized RMS angle meets 1.05*theta_gate.

Writes the full-atlas result back into fixtures/emeraldbay_calibration.json (gating_over_accuracy +
gating_accuracy_scope). Run: python scripts/emeraldbay_recompute/pass5_gating_full.py
"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute   # OVER classification + n* via the single engine (two-arm, real per-line pool)
OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
# EB_CALIB_FIX must point at the SAME fixture pass3 wrote (this pass updates it in place with the
# full-atlas gating result). Default: the committed fixtures/emeraldbay_calibration.json.
FIX = os.environ.get("EB_CALIB_FIX",
                     os.path.join(os.path.dirname(__file__), "..", "..", "fixtures", "emeraldbay_calibration.json"))
THETA_GATE = 0.20
REPS = 200
rng = np.random.default_rng(0)

d = int(np.load(os.path.join(OUT, "basis.npz"))["n_comps"])
z = np.load(os.path.join(OUT, "all_cells.npz"), allow_pickle=True)
coords = z["coords"].astype(np.float64); samp = z["sample"].astype(str); line = z["line"].astype(str)
print(f"[gating-full] {len(coords):,} cells, {len(np.unique(line))} lines", flush=True)

line_base = {l: coords[line == l].mean(0) for l in np.unique(line)}
line_n = {l: int((line == l).sum()) for l in np.unique(line)}        # per-line control-pool size = n_c
grp = np.char.add(np.char.add(samp, "|"), line)                      # (sample x line) unit, as pseudobulk
uniq = np.unique(grp)

n_groups = n_over = met = 0
for g in uniq:
    mask = grp == g; N = int(mask.sum())
    if N < 100:
        continue
    n_groups += 1
    l = g.rsplit("|", 1)[1]
    C = coords[mask]; base = line_base[l]
    v = C.mean(0) - base; m = float(np.linalg.norm(v))
    if m <= 0:
        continue
    u = v / m
    ell = C.var(0)                                                   # PCA-diagonal plug-in (per-PC var)
    res = compute(C.mean(0)[None, :], base[None, :], np.array([N], float),
                  np.array([float(line_n[l])], float), np.diag(ell), theta=THETA_GATE)
    if res["regime"][0] != "OVER":                                  # engine: two-arm quota, real n_c
        continue
    nstar = float(res["n_star"][0])
    n_over += 1
    n = int(min(max(np.ceil(nstar), 5), N - 1))
    th2 = np.empty(REPS)
    for j in range(REPS):
        idx = rng.choice(N, n, replace=False)
        vn = C[idx].mean(0) - base; un = vn / np.linalg.norm(vn)
        th2[j] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
    if np.sqrt(th2.mean()) <= THETA_GATE * 1.05:
        met += 1

acc = round(met / n_over, 4) if n_over else 1.0
print(f"[gating-full] all 52 lines: groups>=100 = {n_groups}; predicted OVER = {n_over}; "
      f"met tolerance = {met} -> accuracy {acc}", flush=True)

fix = json.load(open(FIX))
fix["gating_n_groups"] = int(n_groups)
fix["gating_n_over"] = int(n_over)
fix["gating_pct_over"] = round(100 * n_over / n_groups, 1)
fix["gating_over_accuracy"] = acc
fix["gating_accuracy_scope"] = (f"downsample-and-measure verified over the FULL atlas (all 52 lines): "
                                f"{n_over} predicted OVER, {met} met tolerance")
json.dump(fix, open(FIX, "w"), indent=1)
print(f"[gating-full] updated {FIX}", flush=True)
