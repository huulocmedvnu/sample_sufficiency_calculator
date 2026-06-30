"""
calibrate.py — empirical calibration / verifiable demonstration of the cell-quota calculator.

Pulls real perturbation magnitudes from the cached Tahoe-100M drug-similarity array
(`sig_excl3_corrected.npz`, the batch-clean 292-drug matrix, RESEARCH_LOG §27/§30) and combines them
with the per-cell PCA variance (sigma^2) calibrated from the plate-6 checkpoint subsample, to print
the concrete cell quota n* for a STRONG signature (Resveratrol, a validated mTOR-inhibitor hit, §30)
versus a WEAK signature.

Honest note: these are computed from the data, not assumed. With the real per-cell PCA variance
(sigma^2 ~ 7.66) the quotas land in the THOUSANDS of cells/arm at realistic tolerances (0.05-0.2 rad),
not the hundreds — single-cell embeddings are noisy. (An earlier external draft quoting n*=361/2661
is not reproducible from this data and is not used here.)

Run:
    python src/calibrate.py
    DATA=/path/to/sig_excl3_corrected.npz python src/calibrate.py
"""
import os
import json
import numpy as np
from calculator import (calculate_experimental_cell_quota, rms_angular_error,
                        calculate_optimal_resource_allocation,
                        calculate_cell_quota_anisotropic)

FIXTURE = os.path.join(os.path.dirname(__file__), "..", "fixtures", "tahoe_calibration.json")

# sigma^2 = mean per-cell variance per PCA dim, calibrated once from the plate-6 checkpoint subsample
# (60k cells projected through the excl3 PCA, pca_excl3.npz): mean over 50 dims = 7.66 (range 0.03-31.7).
SIGMA2 = 7.66
NUM_DIMS = 50
# baseline cells per (drug,line) well in the excl3 atlas (15,200 groups, 28.7M cells): MEDIAN = 1394.
BASELINE_CELLS_PER_WELL = 1394

# Default path to the cached array. We use the RAW perturbation array (sig_excl3.npz), NOT the
# consensus-corrected one: sigma^2 (=7.66) was measured in the RAW PCA space, so the magnitude m must
# come from the same RAW space for a consistent power calculation. (The consensus correction is a
# batch-removal step for the similarity matrix; it compresses magnitudes into a different scale.)
DEFAULT_DATA = os.environ.get(
    "DATA", "/mnt/hdd2/loc-tran/obgyn/outputs/drugsim_cache/sig_excl3.npz")

# fallback magnitudes (raw perturbation ||v|| in PCA space) if the array is unavailable
FALLBACK = {"Resveratrol": 2.969, "weak(p25)": 1.328}


def load_magnitudes(path):
    if not os.path.exists(path):
        print(f"[calibrate] {path} not found -> using documented fallback magnitudes")
        return FALLBACK
    z = np.load(path, allow_pickle=True)
    drugs = list(z["drugs"]); mag = z["mag"]
    di = {d: i for i, d in enumerate(drugs)}
    out = {"Resveratrol": float(mag[di["Resveratrol"]])} if "Resveratrol" in di else {}
    out["weak(p25)"] = float(np.nanpercentile(mag, 25))
    return out


def main():
    mags = load_magnitudes(DEFAULT_DATA)
    print(f"sigma^2 (per-cell PCA variance) = {SIGMA2}   d = {NUM_DIMS} dims")
    print("magnitudes m = ||v|| (PCA space):", {k: round(v, 3) for k, v in mags.items()})
    print("\nCELL QUOTA n* (cells PER ARM) by tolerance:")
    print(f"{'tolerance(rad)':>14s} {'tolerance(deg)':>14s} " +
          " ".join(f"{k:>16s}" for k in mags))
    for theta in (0.05, 0.10, 0.20, 0.30):
        row = [calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, m, theta) for m in mags.values()]
        print(f"{theta:>14.2f} {np.degrees(theta):>14.1f} " +
              " ".join(f"{n:>16,.0f}" for n in row))

    # verification: the ratio weak/strong should equal (m_strong/m_weak)^2 (the m^2 law)
    if "Resveratrol" in mags:
        r = (mags["Resveratrol"] / mags["weak(p25)"]) ** 2
        nr = calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, mags["Resveratrol"], 0.1)
        nw = calculate_experimental_cell_quota(SIGMA2, NUM_DIMS, mags["weak(p25)"], 0.1)
        print(f"\n[verify] weak/strong quota ratio = {nw/nr:.2f}  ==  (m_strong/m_weak)^2 = {r:.2f}")
        print(f"[verify] strong-signature quota at 0.1 rad = {nr:,.0f} cells/arm "
              f"(RMS angle achieved by 1000 cells = {np.degrees(rms_angular_error(SIGMA2,NUM_DIMS,mags['Resveratrol'],1000)):.1f} deg)")

    # ---- DUAL-SIDED resource allocation vs the real atlas baseline (N0 = 1394 cells/well) ----
    print(f"\n=== DUAL-SIDED resource allocation @ tolerance=0.1 rad, baseline N0={BASELINE_CELLS_PER_WELL} cells/well ===")
    examples = dict(mags); examples["strong-cytotoxic(max)"] = 14.35
    print(f"{'signature':24s} {'m':>6s} {'n* /well':>9s} {'regime':>12s} "
          f"{'dry save(lin)':>13s} {'dry save(quad)':>14s} {'wet multiplex':>13s}")
    for lab, m in examples.items():
        lin = calculate_optimal_resource_allocation(SIGMA2, NUM_DIMS, m, 0.1,
                                                    BASELINE_CELLS_PER_WELL, "linear")
        quad = calculate_optimal_resource_allocation(SIGMA2, NUM_DIMS, m, 0.1,
                                                     BASELINE_CELLS_PER_WELL, "quadratic")
        reg = "OVER" if lin["wet_lab_multiplex_gain"] > 1 else "UNDER"
        mult = f"{lin['wet_lab_multiplex_gain']:.1f}x" if reg == "OVER" else "-"
        print(f"{lab:24s} {m:>6.2f} {lin['required_cells_per_well']:>9,.0f} {reg:>12s} "
              f"{100*lin['dry_lab_compute_reduction_ratio']:>12.0f}% "
              f"{100*quad['dry_lab_compute_reduction_ratio']:>13.0f}% {mult:>13s}")
    print("\nHonest read: at a tight 0.1-rad tolerance the Tahoe atlas is UNDER-sampled for moderate/weak\n"
          "signatures (n* > N0 -> 0% savings, cells are NOT redundant); only very strong perturbers are\n"
          "over-sampled and safely downsamplable. Downsampling here applies to centroid/pseudobulk\n"
          "analysis only, and the speedup is polynomial (linear-to-quadratic), never exponential.")
    anisotropic_demo()


def anisotropic_demo():
    """ANISOTROPIC quota (docs/THEORY.md) on the real PCA-space covariance, vs the isotropic form."""
    if not os.path.exists(FIXTURE):
        print("\n[aniso] fixture missing; skipping anisotropic demo"); return
    fx = json.load(open(FIXTURE))
    ell = np.array(fx["per_component_variance"]); d = fx["num_dimensions"]
    theta = 0.1
    print(f"\n=== ANISOTROPIC quota (docs/THEORY.md) vs isotropic, real Sigma (trace={ell.sum():.0f}, "
          f"range {ell.min():.2f}..{ell.max():.2f}) ===")
    # consistency: Sigma = sigma^2 I must reproduce the isotropic formula
    v0 = np.zeros(d); v0[0] = 2.97
    chk = calculate_cell_quota_anisotropic(ell.mean() * np.ones(d), v0, theta)
    print(f"[consistency] Sigma=sigma^2 I: anisotropic n*={chk['required_cells_per_arm']:.0f} "
          f"== isotropic n*={calculate_experimental_cell_quota(ell.mean(),d,2.97,theta):.0f}  "
          f"(ratio {chk['isotropic_ratio']:.3f})")
    print(f"\n{'drug':18s} {'m':>6s} {'n*_iso':>9s} {'n*_aniso':>10s} {'aniso/iso':>10s} {'d_eff':>6s}")
    for name, vec in fx["example_perturbation_vectors"].items():
        v = np.array(vec); m = float(np.linalg.norm(v))
        r = calculate_cell_quota_anisotropic(ell, v, theta)
        n_iso = calculate_experimental_cell_quota(ell.mean(), d, m, theta)
        print(f"{name[:18]:18s} {m:>6.2f} {n_iso:>9,.0f} {r['required_cells_per_arm']:>10,.0f} "
              f"{r['isotropic_ratio']:>10.3f} {r['effective_dimensions']:>6.1f}")
    # tail / confidence quota
    v = np.array(fx["example_perturbation_vectors"]["Resveratrol"])
    rc = calculate_cell_quota_anisotropic(ell, v, theta, confidence=0.05)
    print(f"\n[tail] Resveratrol: mean n*={rc['required_cells_per_arm']:,.0f}; "
          f"95%-confident (P(theta>0.1)<=0.05) n*={rc['required_cells_per_arm_confident']:,.0f} "
          f"(x{rc['required_cells_per_arm_confident']/rc['required_cells_per_arm']:.2f}, d_eff={rc['effective_dimensions']:.1f})")
    print("Honest finding: on THIS atlas the anisotropic MEAN-quota correction is small (~2-3%) because\n"
          "drug directions carry little variance along themselves (not aligned with the top noise PCs);\n"
          "the anisotropic value-add is d_eff<<d-1 and the rigorous tail quota. The correction would be\n"
          "LARGE on data where perturbations align with dominant (cell-cycle/lineage) axes.")


if __name__ == "__main__":
    main()
