#!/usr/bin/env python
"""Full-atlas (50-line) downsample-and-measure gating verification for Tahoe-100M.

The EmeraldBay gating check (pass5_gating_full.py) was extended to all 52 lines; this is the Tahoe
counterpart, run over the whole atlas from the per-cell memmap written by pass4a_stream.py. For every
(sample x line) condition classified OVER (N >= n*, per-line-mean baseline, theta_star = 0.1 rad), we
downsample to n* and check that the realized RMS angle meets 1.05*theta_star -- the same protocol used
for EmeraldBay, now applied across all 50 Tahoe cell lines rather than any subset.

Classification (n*) is read from pass4b's per-condition records (cheap); only the OVER conditions are
gathered from the 20 GB memmap and downsampled. Writes the result into fixtures/tahoe_direct_curves_full.json.
Env: THETA(0.1), MIN_CELLS(400), REPS(200).
"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute   # quota + OVER classification: two-arm engine, real per-line pool n_c
BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"
_dose = "/mnt/hdd2/loc-tran/tahoe_work/out_dose"
BASIS = _dose if os.path.exists(os.path.join(_dose, "basis.npz")) else "/mnt/hdd2/loc-tran/tahoe_work/out"
HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "..", "..", "fixtures", "tahoe_direct_curves_full.json")
THETA = float(os.environ.get("THETA", "0.1")); MIN_CELLS = int(os.environ.get("MIN_CELLS", "400"))
REPS = int(os.environ.get("REPS", "200"))
rng = np.random.default_rng(0)

NC = int(np.load(os.path.join(BASIS, "basis.npz"))["n_comps"])
z = np.load(os.path.join(BASE, "ckpt.npz"), allow_pickle=True)
cursor = int(z["cursor"]); cond_keys = list(z["cond_keys"])
line_sum = dict(z["line_sum"].item()); line_cnt = dict(z["line_cnt"].item())
baseline = {l: line_sum[l] / line_cnt[l] for l in line_sum}

# candidates = every condition with enough cells; the ENGINE (below) makes the OVER/UNDER call.
# (No quota math here -- classification is engine.compute's job, per the single-source-of-truth rule.)
recs = json.load(open(os.path.join(BASE, "curves_progress.json")))["records"]
over = [r for r in recs if r["N"] >= MIN_CELLS]
lines_over = sorted(set(r["line"] for r in over))
print(f"[gating-tahoe] {len(recs)} conditions; predicted OVER at theta={THETA}: {len(over)} "
      f"across {len(lines_over)} lines", flush=True)

coords = np.memmap(os.path.join(BASE, "coords.f32"), dtype=np.float32, mode="r", shape=(cursor, NC))
condid = np.memmap(os.path.join(BASE, "condid.i32"), dtype=np.int32, mode="r", shape=(cursor,))
print("[gating-tahoe] grouping rows by condition id ...", flush=True)
order = np.argsort(condid, kind="stable")
cid_sorted = condid[order]
bnd = np.searchsorted(cid_sorted, np.arange(len(cond_keys) + 1))

met = 0; checked = 0
for r in over:
    j = r["cid"]; a, b = bnd[j], bnd[j + 1]
    line = cond_keys[j].rsplit("|", 1)[1]
    if line not in baseline:
        continue
    C = np.asarray(coords[order[a:b]], dtype=np.float64); N = len(C)
    base = baseline[line]; v = C.mean(0) - base; m = float(np.linalg.norm(v))
    if m <= 0:
        continue
    u = v / m; ell = C.var(0)
    # exact quota + regime via the single engine (two-arm form, real per-line control pool n_c)
    res = compute(C.mean(0)[None, :], base[None, :], np.array([N], float),
                  np.array([float(line_cnt[line])], float), np.diag(ell), theta=THETA)
    nstar = float(res["n_star"][0])
    if res["regime"][0] != "OVER":         # engine says N < n* (or not detectable) -> not predicted-over
        continue
    checked += 1
    n = int(min(max(np.ceil(nstar), 5), N - 1))
    th2 = np.empty(REPS)
    for k in range(REPS):
        idx = rng.choice(N, n, replace=False)
        vn = C[idx].mean(0) - base; un = vn / np.linalg.norm(vn)
        th2[k] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
    if np.sqrt(th2.mean()) <= THETA * 1.05:
        met += 1
    if checked % 1000 == 0:
        print(f"  verified {checked}/{len(over)} (met {met})", flush=True)

acc = round(met / checked, 4) if checked else 1.0
print(f"[gating-tahoe] FULL atlas ({len(lines_over)} lines): predicted OVER {checked}; "
      f"met tolerance {met} -> accuracy {acc}", flush=True)

fix = json.load(open(FIX))
fix["gating"] = dict(theta_star=THETA, baseline="per-line mean", scope=f"full atlas, all {len(lines_over)} lines",
                     min_cells=MIN_CELLS, n_predicted_over=checked, n_met_tolerance=met, accuracy=acc,
                     criterion="realized RMS angle <= 1.05*theta_star when downsampled to n*")
json.dump(fix, open(FIX, "w"), indent=1)
print(f"[gating-tahoe] updated {FIX}", flush=True)
