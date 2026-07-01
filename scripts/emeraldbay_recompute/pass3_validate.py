#!/usr/bin/env python
"""Pass 3 of the EmeraldBay recompute: the INDEPENDENT validation of docs/THEORY.md.

Using the Pass-2 outputs (per-condition moments + retained per-cell PCA coordinates for the five
shared cell lines), this computes:
  * the WITHIN-condition per-PC variance (sigma^2), the platform-specific plug-in;
  * held-out ANGULAR-ERROR CURVES: for high-cell-count (condition x line) groups, subsample n cells,
    measure the realized RMS angle of the estimated direction (relative to the per-line baseline)
    against the closed form  theta^2(n) = tr(P Sigma P)/m^2 * (1/n - 1/N)  (finite-population corrected);
  * regime GATING at theta = 0.20 rad: classify OVER/UNDER and verify predicted-OVER groups meet the
    tolerance when downsampled to n*.

Writes fixtures/emeraldbay_calibration.json (schema consumed by tests/test_emeraldbay_integration.py).
"""
import os, json, numpy as np
import pyarrow.parquet as pq

OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
META = "/mnt/hdd2/loc-tran/eb_work/meta/metadata"
REPO_FIX = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures", "emeraldbay_calibration.json")
THETA_GATE = 0.20
REPS = 300
rng = np.random.default_rng(0)

b = np.load(os.path.join(OUT, "basis.npz")); d = int(b["n_comps"])
pb = np.load(os.path.join(OUT, "pseudobulk.npz"), allow_pickle=True)
sums, sumsq, counts = pb["sums"].astype(float), pb["sumsq"].astype(float), pb["counts"].astype(float)
sh = np.load(os.path.join(OUT, "shared_cells.npz"), allow_pickle=True)
coords = sh["coords"].astype(float); samp = sh["sample"].astype(str); cline = sh["line"].astype(str)
s2cond = json.load(open(os.path.join(OUT, "sample2cond.json")))
clm = pq.read_table(f"{META}/cell_line_metadata.parquet").to_pandas()
cvcl2name = dict(clm.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur")["cell_name"])

# --- within-condition sigma^2 (pooled residual variance per PC over all conditions) ---
ss_res = (sumsq - sums**2 / counts[:, None]).sum(0)         # per-PC pooled within-condition SS
dof = counts.sum() - len(counts)
per_pc = ss_res / dof
sigma2_within = float(per_pc.mean())

# --- per-cell condition label + per-line baseline (mean over all that line's cells) ---
def cond_name(s):
    import ast
    try:
        return str(ast.literal_eval(s2cond[s])[0][0])
    except Exception:
        return s2cond.get(s, "?")
cond_cell = np.array([cond_name(s) for s in samp])
baseline = {l: coords[cline == l].mean(0) for l in np.unique(cline)}


def curve(mask, line):
    C = coords[mask]; N = len(C); base = baseline[line]
    v = C.mean(0) - base; m = float(np.linalg.norm(v)); u = v / m
    P = np.eye(d) - np.outer(u, u)
    Sig = np.cov(C.T)
    trPSP = float(np.trace(P @ Sig @ P))
    ns = np.unique(np.geomspace(20, max(N // 2, 40), 9).round().astype(int))
    rows, xs, ys = [], [], []
    for n in ns:
        if n >= N:
            continue
        th2 = np.empty(REPS)
        for j in range(REPS):
            idx = rng.choice(N, n, replace=False)
            vn = C[idx].mean(0) - base; un = vn / np.linalg.norm(vn)
            th2[j] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
        realized = float(np.sqrt(th2.mean()))
        predicted = float(np.sqrt(max(trPSP / m**2 * (1.0 / n - 1.0 / N), 0)))
        rows.append([int(n), realized, predicted]); xs.append(1.0 / n - 1.0 / N); ys.append(realized**2)
    xs, ys = np.array(xs), np.array(ys)
    slope = float((xs @ ys) / (xs @ xs))                    # through-origin fit
    ss_tot = float(((ys - ys.mean())**2).sum()); ss_res2 = float(((ys - slope * xs)**2).sum())
    r2 = 1.0 - ss_res2 / ss_tot if ss_tot > 0 else 1.0
    rel = [abs(r - p) / p for _, r, p in rows if p > 0]
    return dict(N=int(N), m=round(m, 4), slope=round(slope, 5), expected_slope=round(trPSP / m**2, 5),
                r2=round(r2, 5), intercept=0.0, mean_rel_err=round(float(np.mean(rel)), 5),
                trPSP=trPSP, rows=rows)


# --- pick held-out groups: highest-m (condition x line) with enough cells, across the 5 shared lines ---
cand = []
for l in np.unique(cline):
    for cd in np.unique(cond_cell[cline == l]):
        mask = (cline == l) & (cond_cell == cd)
        N = int(mask.sum())
        if N < 400:
            continue
        v = coords[mask].mean(0) - baseline[l]
        cand.append((float(np.linalg.norm(v)), cd, l, N))
cand.sort(reverse=True)
# take the top-3 by magnitude (typically DMSO_T0 time-zero groups) PLUS the strongest DRUG group
# (non-DMSO_T0) so the validation is not restricted to time-zero populations.
picks = cand[:3]
for c in cand:
    if c[1] != "DMSO_T0":
        picks.append(c); break
heldout = {}
for m, cd, l, N in picks:
    key = f"{cd}|{cvcl2name.get(l, l)}"
    heldout[key] = curve((cline == l) & (cond_cell == cd), l)
    print(f"[curve] {key}: N={N} m={m:.2f} r2={heldout[key]['r2']:.4f} "
          f"mean_rel_err={heldout[key]['mean_rel_err']:.4f} slope {heldout[key]['slope']} vs {heldout[key]['expected_slope']}")

# --- gating at theta=0.20 over all shared-line (condition x line) groups with >=100 cells ---
n_over = met = 0; groups = 0
for l in np.unique(cline):
    for cd in np.unique(cond_cell[cline == l]):
        mask = (cline == l) & (cond_cell == cd); N = int(mask.sum())
        if N < 100:
            continue
        groups += 1
        C = coords[mask]; base = baseline[l]; v = C.mean(0) - base; m = np.linalg.norm(v)
        if m <= 0:
            continue
        u = v / m; P = np.eye(d) - np.outer(u, u); trPSP = np.trace(P @ np.cov(C.T) @ P)
        nstar = 2.0 * trPSP / (m**2 * THETA_GATE**2)
        if N >= nstar:                                       # predicted OVER
            n_over += 1
            n = int(min(max(np.ceil(nstar), 5), N - 1))
            th2 = np.empty(200)
            for j in range(200):
                idx = rng.choice(N, n, replace=False)
                vn = C[idx].mean(0) - base; un = vn / np.linalg.norm(vn)
                th2[j] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
            if np.sqrt(th2.mean()) <= THETA_GATE * 1.05:
                met += 1
gate_acc = round(met / n_over, 4) if n_over else 1.0
print(f"[gating] groups>=100 cells: {groups}; predicted OVER: {n_over}; met tol: {met} -> acc {gate_acc}")

fix = dict(
    description="EmeraldBay external calibration (from-raw recompute; scripts/emeraldbay_recompute/): "
                "within-condition per-PC variance + held-out angular-error curves (independent validation "
                "of THEORY.md). Baseline = per-line mean; magnitudes are 5-day survivor norms.",
    num_dimensions=d, per_component_variance=[round(float(x), 5) for x in per_pc],
    sigma2_mean=round(sigma2_within, 4), heldout_curves=heldout,
    gating_tolerance=THETA_GATE, gating_over_accuracy=gate_acc)
json.dump(fix, open(REPO_FIX, "w"), indent=1)
print(f"[pass3-eb] sigma^2(within)={sigma2_within:.4f}; wrote {REPO_FIX}")
