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
import numpy as np
from calculator import calculate_experimental_cell_quota, rms_angular_error

# sigma^2 = mean per-cell variance per PCA dim, calibrated once from the plate-6 checkpoint subsample
# (60k cells projected through the excl3 PCA): mean over the 50 dims = 7.66 (per-dim range 0.03-31.7).
SIGMA2 = 7.66
NUM_DIMS = 50

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


if __name__ == "__main__":
    main()
