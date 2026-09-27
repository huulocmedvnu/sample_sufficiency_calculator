"""
calibrate.py -- verifiable demonstration of the cell-quota calculator on the Tahoe-100M calibration.

Uses the committed fixtures (per-PC within-condition variance, example per-drug perturbation vectors in
NCI-H460). Every quota comes from src/calculator.py, which calls the single engine (src/engine.py); this
script assembles no formula of its own. The recommended quota is the two-arm form
n_t* = 1/(m^2 theta^2/tr(P Sigma P) - 1/n_c): the control-pool size n_c is REQUIRED, because against
Tahoe's real shared DMSO pool (median n_c ~ 3,113) most conditions are control-pool-limited (n* = inf) -- the
error this tool exists to prevent. The large-pool limit (d-1)sigma^2/theta^2 and the equal-arm 2x it are
shown as the special cases they are.

Run:  python src/calibrate.py
"""
import os
import json
import math
import numpy as np
from calculator import (cell_quota, cell_quota_large_pool, cell_quota_equal_arm,
                        cell_quota_report, resource_allocation)

FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
CAL = json.load(open(os.path.join(FIX, "tahoe_calibration.json")))
CONST = json.load(open(os.path.join(FIX, "tahoe_constants.json")))

SIGMA2 = float(CAL["sigma2_mean_within_condition"])          # 0.9567 within-condition residual (headline)
SIGMA2_MARGINAL = float(CAL["sigma2_mean"])                  # 2.406 marginal (conservative bound only)
NUM_DIMS = int(CAL["num_dimensions"])                        # 50
BASELINE = float(CONST["median_N0"])                         # 1296 median cells / (drug x dose x line)
SPEC = json.load(open(os.path.join(FIX, "unified_spectrum.json")))
TAHOE = next(s for s in SPEC["headline"] if s["name"] == "Tahoe-100M")
N_C_VEHICLE = int(TAHOE["median_n_c"])                        # Tahoe's real shared DMSO pool (median, all wells pooled)
ELL = np.array(CAL["per_component_variance"]) * (SIGMA2 / SIGMA2_MARGINAL)   # within-cond Sigma diagonal
VECS = {k: np.array(v) for k, v in CAL["example_perturbation_vectors"].items()}   # drug vectors (NCI-H460)


def main():
    mags = dict(sorted({k: float(np.linalg.norm(v)) for k, v in VECS.items()}.items(), key=lambda kv: -kv[1]))
    print(f"sigma^2 (within-condition per-cell PCA variance) = {SIGMA2}   d = {NUM_DIMS}   "
          f"large-pool bound (d-1)sigma^2/theta^2 = {CONST['nstar_const']:,} / m^2   (equal-arm = 2x; "
          f"the operational quota is the two-arm form with the REAL n_c)   line: {CAL['example_line']}")
    print("magnitudes m = ||v||:", {k: round(v, 2) for k, v in mags.items()})

    print("\nLARGE-POOL quota n* (cells) by tolerance:")
    print(f"{'tol(rad)':>9s} {'tol(deg)':>9s} " + " ".join(f"{k[:16]:>16s}" for k in mags))
    for theta in (0.05, 0.10, 0.20, 0.30):
        row = [cell_quota_large_pool(VECS[k], ELL, theta) for k in mags]
        print(f"{theta:>9.2f} {np.degrees(theta):>9.1f} " + " ".join(f"{n:>16,.0f}" for n in row))

    # m^2 law: large-pool quota ratio of two drugs == (m_a/m_b)^2
    a, b = list(mags)[0], list(mags)[-1]
    na = cell_quota_large_pool(VECS[a], ELL, 0.1); nb = cell_quota_large_pool(VECS[b], ELL, 0.1)
    print(f"\n[verify] {b}/{a} quota ratio = {nb/na:.2f}  ==  (m_{a[:4]}/m_{b[:4]})^2 = {(mags[a]/mags[b])**2:.2f}")

    # --- the two-arm quota against the REAL shared DMSO vehicle: control-pool limitation ---
    print(f"\n=== two-arm quota vs Tahoe's real shared DMSO vehicle (n_c = {N_C_VEHICLE}) @ 0.1 rad ===")
    print(f"{'drug':22s} {'m':>6s} {'m_min':>7s} {'large-pool n*':>14s} {'shared-vehicle n*':>18s}")
    for name in mags:
        v = VECS[name]; rep = cell_quota_report(v, ELL, 0.1, control_pool_size=N_C_VEHICLE)
        lp = cell_quota_large_pool(v, ELL, 0.1)
        sv = rep["required_cells_treated"]
        sv_str = "inf (POOL-LIMITED)" if not math.isfinite(sv) else f"{sv:,.0f}"
        print(f"{name[:22]:22s} {mags[name]:>6.2f} {rep['m_min']:>7.2f} {lp:>14,.0f} {sv_str:>18s}")
    print(f"Across the {TAHOE['n_conditions']:,} (drug x dose x line) conditions at 0.1 rad against the real vehicle: "
          f"{TAHOE['pct_over']}% over, {TAHOE['pct_pool_limited']}% control-pool-limited, "
          f"{TAHOE['pct_under']}% treated-depth-limited (fixtures/unified_spectrum.json).")

    # --- DUAL-SIDED resource allocation (large per-line pool: every quota finite) ---
    print(f"\n=== dual-sided allocation @ 0.1 rad, baseline N0={BASELINE:.0f}, large per-line pool ===")
    print(f"{'signature':30s} {'m':>6s} {'n*':>9s} {'regime':>6s} {'dry(lin)':>9s} {'dry(quad)':>10s} {'multiplex':>10s}")
    for name in mags:
        v = VECS[name]
        lin = resource_allocation(v, ELL, 0.1, math.inf, BASELINE, "linear")
        quad = resource_allocation(v, ELL, 0.1, math.inf, BASELINE, "quadratic")
        reg = "OVER" if lin["wet_lab_multiplex_gain"] > 1 else "UNDER"
        mult = f"{lin['wet_lab_multiplex_gain']:.1f}x" if reg == "OVER" else "-"
        print(f"{name[:30]:30s} {mags[name]:>6.2f} {lin['required_cells_per_well']:>9,.0f} {reg:>6s} "
              f"{100*lin['dry_lab_compute_reduction_ratio']:>8.0f}% {100*quad['dry_lab_compute_reduction_ratio']:>9.0f}% {mult:>10s}")

    # --- ANISOTROPIC vs isotropic (real Sigma), and the tail-controlled quota ---
    print(f"\n=== anisotropic vs isotropic large-pool quota, real Sigma "
          f"(trace={ELL.sum():.0f}, range {ELL.min():.2f}..{ELL.max():.2f}) @ 0.1 rad ===")
    print(f"{'drug':22s} {'m':>6s} {'n*_iso':>9s} {'n*_aniso':>10s} {'d_eff':>6s}")
    iso = ELL.mean() * np.ones(NUM_DIMS)
    for name in mags:
        v = VECS[name]; m = mags[name]
        n_iso = cell_quota_large_pool(v, iso, 0.1)
        rep = cell_quota_report(v, ELL, 0.1, control_pool_size=math.inf)
        print(f"{name[:22]:22s} {m:>6.2f} {n_iso:>9,.0f} {cell_quota_large_pool(v, ELL, 0.1):>10,.0f} "
              f"{rep['effective_dimensions']:>6.1f}")
    strong = max(VECS, key=lambda k: np.linalg.norm(VECS[k]))
    rc = cell_quota_report(VECS[strong], ELL, 0.1, control_pool_size=math.inf, confidence=0.05)
    mean_n = rc["required_cells_treated"]; conf_n = rc["required_cells_treated_confident"]
    print(f"[tail] {strong}: large-pool n*={mean_n:,.0f}; 95%-confident (P(theta>0.1)<=0.05) n*={conf_n:,.0f} "
          f"(x{conf_n/mean_n:.2f}, d_eff={rc['effective_dimensions']:.1f})")


if __name__ == "__main__":
    main()
