#!/usr/bin/env python
"""Full recalibration of the Tahoe headline from the MARGINAL sigma^2=2.406 to the theory-preferred
WITHIN-condition sigma^2=0.9567. Magnitudes m and depths N0 are sigma^2-independent, so every n* just
rescales; the regime tallies, per-drug/per-dose tables, constants, and screen budget are recomputed.
EmeraldBay tables (held-out curves, invariance) use EmeraldBay's own Sigma and are unaffected.

Regenerates fixtures/tahoe_quota_per_condition.csv, tahoe_per_drug.csv, tahoe_per_dose.csv,
tahoe_constants.json and prints every number needed to propagate into the manuscript/docs.
"""
import os, json, numpy as np, pandas as pd
FIX = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
SIG_NEW = 0.9567; SIG_OLD = 2.406; d = 50; th = 0.1
const_new = 2 * (d - 1) * SIG_NEW / th ** 2                     # 9375.66
GHOST = 50000.0

qc = pd.read_csv(os.path.join(FIX, "tahoe_quota_per_condition.csv"))
qc["n_star"] = const_new / qc["m"] ** 2
def regime(n, N0): return np.where(n < N0, "OVER", np.where(n > GHOST, "Ghost", "UNDER"))
qc["regime"] = regime(qc["n_star"].values, qc["N0"].values)
qc.to_csv(os.path.join(FIX, "tahoe_quota_per_condition.csv"), index=False)

ns = qc["n_star"].values; N0 = qc["N0"].values; m = qc["m"].values
pO = 100 * (ns < N0).mean(); pG = 100 * (ns > GHOST).mean(); pU = 100 - pO - pG
med = np.median(ns); iqr = (np.percentile(ns, 25), np.percentile(ns, 75))
p90 = float(np.percentile(ns, 90))
boundary_m = float(np.sqrt(const_new / np.median(N0)))
depth_coeff = float(np.sqrt(const_new * th ** 2 / np.median(N0)))   # theta(N0)=coeff/m

# ---- per-drug ----
pd_rows = []
for drug, g in qc.groupby("drug"):
    o = int((g.n_star < g.N0).sum()); gh = int((g.n_star > GHOST).sum()); u = len(g) - o - gh
    pd_rows.append(dict(drug=drug, median_m=round(g.m.median(), 2), median_nstar=int(round(g.n_star.median())),
                        n_conditions=len(g), OVER=o, UNDER=u, Ghost=gh))
pdd = pd.DataFrame(pd_rows).sort_values("OVER", ascending=False)
pdd.to_csv(os.path.join(FIX, "tahoe_per_drug.csv"), index=False)
drugs_over_0 = int((pdd.OVER == 0).sum())

# ---- per-dose ----
pdose = []
for dose, g in qc.groupby("dose_uM"):
    o = 100 * (g.n_star < g.N0).mean(); gh = 100 * (g.n_star > GHOST).mean(); u = 100 - o - gh
    pdose.append(dict(dose_uM=dose, n=len(g), pct_OVER=round(o, 1), pct_UNDER=round(u, 1),
                      pct_Ghost=round(gh, 1), pct_under_or_ghost=round(u + gh, 1),
                      median_m=round(g.m.median(), 2), median_nstar=int(round(g.n_star.median()))))
pd.DataFrame(pdose).to_csv(os.path.join(FIX, "tahoe_per_dose.csv"), index=False)

# ---- budget (100x3x2=600 conditions drawn from empirical distribution) ----
NC = 600
flat_M = p90 * NC / 1e6
adapt10 = float(np.minimum(ns, 10000).mean()) * NC / 1e6
adapt20 = float(np.minimum(ns, 20000).mean()) * NC / 1e6
ns02 = ns * (0.1 / 0.2) ** 2
theta02_10 = float(np.minimum(ns02, 10000).mean()) * NC / 1e6

const = json.load(open(os.path.join(FIX, "tahoe_constants.json")))
const.update(sigma2=SIG_NEW, nstar_const=round(const_new, 1), nstar_formula=f"n* = {round(const_new):,}/m^2",
             calibration="within-condition residual variance (marginal 2.406 documented as conservative proxy)",
             over_under_boundary_m=round(boundary_m, 3), ghost_boundary_m=round(float(np.sqrt(const_new / GHOST)), 3),
             depth_fixed_coeff=round(depth_coeff, 3),
             depth_fixed_resolution=f"theta(N0) = {round(depth_coeff,3)}/m rad at median N0={int(np.median(N0))}",
             median_nstar=int(round(med)), median_nstar_iqr=[int(round(iqr[0])), int(round(iqr[1]))],
             pct_OVER=round(pO, 1), pct_UNDER=round(pU, 1), pct_Ghost=round(pG, 1),
             pct_under_or_ghost=round(pU + pG, 1), drugs_over_in_0_conditions=drugs_over_0,
             screen_budget_100x3x2=dict(flat_p90=p90, flat_M=round(flat_M, 1),
                                        adaptive_cap10k_M=round(adapt10, 1), adaptive_cap20k_M=round(adapt20, 1),
                                        theta02_cap10k_M=round(theta02_10, 1)))
json.dump(const, open(os.path.join(FIX, "tahoe_constants.json"), "w"), indent=1)

# ---- Discussion dual-use exemplar: homoharringtonine 5uM x NCI-H460 (m=15.53, N0=6060) ----
mm = 15.53; N0h = 6060; nh = const_new / mm ** 2
print("=== HEADLINE (within-condition sigma^2=0.9567) ===")
print(f"n* = {round(const_new):,}/m^2   (was 23,577)")
print(f"spectrum: {pO:.1f}% OVER / {pU:.1f}% UNDER / {pG:.1f}% Ghost  (97.5->{pU+pG:.1f}% under-or-ghost)")
print(f"median n* = {med:,.0f} (IQR {iqr[0]:,.0f}-{iqr[1]:,.0f}); boundary m = {boundary_m:.2f}; theta(N0)={depth_coeff:.3f}/m")
print(f"drugs OVER in 0 conditions = {drugs_over_0} (was 132)")
print(f"budget: flat p90={p90:,.0f} -> {flat_M:.1f}M; adaptive cap10k={adapt10:.1f}M; cap20k={adapt20:.1f}M; theta0.2 cap10k={theta02_10:.1f}M")
print(f"dual-use homoharringtonine 5uM NCI-H460: n*={nh:.0f} (was 98), multiplex={N0h/nh:.0f}x (was 60x), "
      f"dry-save lin={100*(1-nh/N0h):.0f}% quad={100*(1-(nh/N0h)**2):.1f}%")
print("\n=== TABLE 1 (per-drug) ===")
for name in ["Panobinostat","Homoharringtonine","Harringtonine","Idarubicin (hydrochloride)","Trametinib","palbociclib","4EGI-1","(S)-Crizotinib"]:
    r = pdd[pdd.drug == name]
    if len(r): r = r.iloc[0]; print(f"  {name:26s} med_m={r.median_m:.2f} med_n*={r.median_nstar:,} OVER={r.OVER} UNDER={r.UNDER} Ghost={r.Ghost}")
print("\n=== TABLE 2 (per-dose) ==="); print(pd.DataFrame(pdose).to_string(index=False))
print("\nregenerated: tahoe_quota_per_condition.csv, tahoe_per_drug.csv, tahoe_per_dose.csv, tahoe_constants.json")
