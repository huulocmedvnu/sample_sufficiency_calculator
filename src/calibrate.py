"""
calibrate.py — empirical calibration / verifiable demonstration of the cell-quota calculator.

Uses the FRESH Tahoe-100M recompute (from the raw 100.6M-cell dataset via the theislab vevo_100m
recipe: normalize_total(1e4) -> log1p -> HVG(2000) -> PCA(50); see scripts/tahoe_recompute/). The
committed fixtures carry the per-cell PCA variance sigma^2, the per-component variances (anisotropic
Sigma diagonal), and example per-drug perturbation vectors (in the NCI-H460 line).

Calibrated constants (fixtures/chemical_within_cov.json): sigma^2 = 0.9567 (within-condition residual,
real diagonal Sigma; reproduces from the cached per-cell PCA coords). The quota is the two-arm form
n_t* = 1 / (m^2 theta^2 / tr(P Sigma P) - 1/n_c); with a large control pool it approaches the large-pool
bound (d-1) sigma^2 / theta^2 = 4,688 / m^2. Against Tahoe's REAL shared DMSO vehicle (n_c ~ 1,514 cells
shared by ~94 conditions), the magnitude floor m_min = 1.75 exceeds the median magnitude 1.27, so of the
56,827 (drug x dose x line) conditions at 0.1 rad 69% are control-pool-limited (unresolvable at any
treated depth), 10% over-sampled, 14% treated-depth-limited; referencing the per-line mean instead gives
22% over-sampled and is reported as a sensitivity. The single computation lives in src/engine.py.

Run:  python src/calibrate.py
"""
import os
import json
import numpy as np
from calculator import (calculate_experimental_cell_quota, rms_angular_error,
                        calculate_optimal_resource_allocation,
                        calculate_cell_quota_anisotropic)

FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
CAL = json.load(open(os.path.join(FIX, "tahoe_calibration.json")))
CONST = json.load(open(os.path.join(FIX, "tahoe_constants.json")))

SIGMA2_MARGINAL = float(CAL["sigma2_mean"])   # 2.406 marginal (conservative bound only)
SIGMA2 = float(CAL["sigma2_mean_within_condition"])  # 0.9567 within-condition residual (headline)
NUM_DIMS = int(CAL["num_dimensions"])         # 50
BASELINE = float(CONST["median_N0"])          # 1296 median cells per (drug x dose x line) condition (post-filter)
# anisotropic Sigma diagonal: rescale the marginal per-PC diagonal to the within-condition scale
# (the switch preserves Sigma's SHAPE and only its trace changes, so aniso/iso ratios are unaffected)
ELL = np.array(CAL["per_component_variance"]) * (SIGMA2 / SIGMA2_MARGINAL)
VECS = {k: np.array(v) for k, v in CAL["example_perturbation_vectors"].items()}  # drug vectors (NCI-H460)


def main():
    mags = {k: float(np.linalg.norm(v)) for k, v in VECS.items()}
    mags = dict(sorted(mags.items(), key=lambda kv: -kv[1]))
    print(f"sigma^2 (per-cell PCA variance) = {SIGMA2}   d = {NUM_DIMS} dims   "
          f"->   n* = {CONST['nstar_const']:.0f} / m^2   (example line: {CAL['example_line']})")
    print("magnitudes m = ||v|| (PCA space):", {k: round(v, 2) for k, v in mags.items()})

    print("\nCELL QUOTA n* (cells PER ARM) by tolerance:")
    print(f"{'tolerance(rad)':>14s} {'tolerance(deg)':>14s} " + " ".join(f"{k[:16]:>16s}" for k in mags))
    for theta in (0.05, 0.10, 0.20, 0.30):
        row = [calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, m, theta) for m in mags.values()]
        print(f"{theta:>14.2f} {np.degrees(theta):>14.1f} " + " ".join(f"{n:>16,.0f}" for n in row))

    # m^2 law: quota ratio of two drugs == (m_a/m_b)^2
    names = list(mags)
    a, b = names[0], names[-1]
    na = calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, mags[a], 0.1)
    nb = calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, mags[b], 0.1)
    print(f"\n[verify] {b}/{a} quota ratio = {nb/na:.2f}  ==  (m_{a[:4]}/m_{b[:4]})^2 = "
          f"{(mags[a]/mags[b])**2:.2f}")

    # ---- DUAL-SIDED resource allocation vs the per-condition baseline N0 ----
    print(f"\n=== DUAL-SIDED resource allocation @ tolerance=0.1 rad, baseline N0={BASELINE:.0f} "
          f"cells/condition ===")
    examples = dict(mags); examples[f"median condition (m={CONST['median_m']:.2f})"] = CONST["median_m"]
    print(f"{'signature':30s} {'m':>6s} {'n*/cond':>9s} {'regime':>7s} "
          f"{'dry(lin)':>9s} {'dry(quad)':>10s} {'multiplex':>10s}")
    for lab, m in examples.items():
        lin = calculate_optimal_resource_allocation(SIGMA2, NUM_DIMS, m, 0.1, BASELINE, "linear")
        quad = calculate_optimal_resource_allocation(SIGMA2, NUM_DIMS, m, 0.1, BASELINE, "quadratic")
        reg = "OVER" if lin["wet_lab_multiplex_gain"] > 1 else "UNDER"
        mult = f"{lin['wet_lab_multiplex_gain']:.1f}x" if reg == "OVER" else "-"
        print(f"{lab[:30]:30s} {m:>6.2f} {lin['required_cells_per_well']:>9,.0f} {reg:>7s} "
              f"{100*lin['dry_lab_compute_reduction_ratio']:>8.0f}% "
              f"{100*quad['dry_lab_compute_reduction_ratio']:>9.0f}% {mult:>10s}")
    print(f"\nHonest read: across the full 379-drug x 3-dose x 50-line panel, at 0.1 rad only "
          f"{CONST['pct_OVER']:.1f}% of conditions are OVER-sampled; "
          f"{CONST['pct_under_or_ghost']:.1f}% are UNDER-sampled or Ghost (n*>50k). "
          f"{CONST['drugs_over_in_0_conditions']}/{CONST['n_drugs']} drugs are resolvable in ZERO conditions at "
          f"this depth. Downsampling savings are polynomial, never exponential.")
    anisotropic_demo(mags)


def anisotropic_demo(mags):
    """ANISOTROPIC quota (docs/THEORY.md) on the real PCA-space covariance, vs the isotropic form."""
    d = NUM_DIMS; theta = 0.1
    print(f"\n=== ANISOTROPIC quota (docs/THEORY.md) vs isotropic, real Sigma "
          f"(trace={ELL.sum():.0f}, range {ELL.min():.2f}..{ELL.max():.2f}) ===")
    v0 = np.zeros(d); v0[0] = 5.10
    chk = calculate_cell_quota_anisotropic(ELL.mean() * np.ones(d), v0, theta)
    print(f"[consistency] Sigma=sigma^2 I: anisotropic n*={chk['required_cells_per_arm']:.0f} "
          f"== isotropic n*={calculate_experimental_cell_quota(ELL.mean(), d, 5.10, theta):.0f}  "
          f"(ratio {chk['isotropic_ratio']:.3f})")
    print(f"\n{'drug':22s} {'m':>6s} {'n*_iso':>9s} {'n*_aniso':>10s} {'aniso/iso':>10s} {'d_eff':>6s}")
    for name, v in VECS.items():
        m = float(np.linalg.norm(v))
        r = calculate_cell_quota_anisotropic(ELL, v, theta)
        n_iso = calculate_experimental_cell_quota(ELL.mean(), d, m, theta)
        print(f"{name[:22]:22s} {m:>6.2f} {n_iso:>9,.0f} {r['required_cells_per_arm']:>10,.0f} "
              f"{r['isotropic_ratio']:>10.3f} {r['effective_dimensions']:>6.1f}")
    # tail / confidence quota on the strongest example drug
    strong = max(VECS, key=lambda k: np.linalg.norm(VECS[k]))
    rc = calculate_cell_quota_anisotropic(ELL, VECS[strong], theta, confidence=0.05)
    print(f"\n[tail] {strong}: mean n*={rc['required_cells_per_arm']:,.0f}; "
          f"95%-confident (P(theta>0.1)<=0.05) n*={rc['required_cells_per_arm_confident']:,.0f} "
          f"(x{rc['required_cells_per_arm_confident']/rc['required_cells_per_arm']:.2f}, "
          f"d_eff={rc['effective_dimensions']:.1f})")


if __name__ == "__main__":
    main()
