#!/usr/bin/env python
"""OVER downsample-and-measure gate for Tahoe-100M against its REAL shared DMSO_TF vehicle (headline estimand).

Regenerates fixtures/verification_over.json (previously produced by an uncommitted one-off; this is the
committed generator). Protocol (unchanged from pass4c, which gates the per-line-mean sensitivity):
for every (drug x dose x line) condition the ENGINE classifies OVER (n* < N0) on the shared-vehicle
control, gather its treated cells from the per-cell PCA memmap (pass4a_stream.py), draw REPS random
subsets of n* cells, and measure the RMS angle between the subset direction and the full-depth direction,
both referenced to the same pooled DMSO_TF centroid. Gate: RMS angle <= 1.05*theta_star.

Note (audit 2026-09): because truth and subsample share the control centroid, this gate exercises only
the treated-arm term of the two-arm quota; it does not test the control-pool (1/n_c) term.

Env: THETA(0.1), REPS(200). Needs /mnt/hdd2/loc-tran/tahoe_work/pass4/{coords.f32,condid.i32,ckpt.npz}.
"""
import os, sys, ast, json, time, numpy as np, pyarrow.parquet as pq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from engine import compute
from pipelines import chemical

BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"
FIX = os.path.join(HERE, "..", "..", "fixtures", "verification_over.json")
THETA = float(os.environ.get("THETA", "0.1")); REPS = int(os.environ.get("REPS", "200"))
rng = np.random.default_rng(0)

# --- spectrum classification: the engine on the shared-vehicle pipeline (single source of truth) ---
d = chemical.load("Tahoe-100M", control="shared-dmso")
res = compute(d["mu_t"], d["mu_c"], d["n_t"], d["n_c"], d["Sigma"], theta=THETA)
over_idx = np.where(res["regime"] == "OVER")[0]
print(f"[gate-vehicle] {len(d['cond_id'])} conditions; predicted OVER at theta={THETA}: {len(over_idx)}; "
      f"median n_c {np.median(d['n_c']):.0f}", flush=True)

# --- memmap layout: cond_keys are 'sample|line'; map each (drug, dose, line) condition to its rows ---
z = np.load(os.path.join(BASE, "ckpt.npz"), allow_pickle=True)
cursor = int(z["cursor"]); cond_keys = [str(k) for k in z["cond_keys"]]
NC = int(np.load(f"{chemical.TAHOE}/basis.npz")["n_comps"])
sm = pq.read_table(f"{chemical.TAHOE_META}/sample_metadata.parquet").to_pandas()
sm["dose"] = sm.drugname_drugconc.map(lambda s: (lambda t: float(t[0][1]))(ast.literal_eval(s)) if s else np.nan)
S2 = {r.sample: (r.drug, r.dose) for r in sm.itertuples()}
groups = {}
for j, k in enumerate(cond_keys):
    s, line = k.rsplit("|", 1); mt = S2.get(s)
    if mt is None or mt[0] == chemical.CTRL: continue
    groups.setdefault(f"{mt[0]}|{mt[1]}|{line}", []).append(j)

coords = np.memmap(os.path.join(BASE, "coords.f32"), dtype=np.float32, mode="r", shape=(cursor, NC))
print("[gate-vehicle] sorting condition ids ...", flush=True)
condid = np.array(np.memmap(os.path.join(BASE, "condid.i32"), dtype=np.int32, mode="r", shape=(cursor,)))
order = np.argsort(condid, kind="stable"); cid_sorted = condid[order]; del condid
bnd = np.searchsorted(cid_sorted, np.arange(len(cond_keys) + 1)); del cid_sorted

def gather(js):
    return np.concatenate([np.asarray(coords[order[bnd[j]:bnd[j + 1]]], dtype=np.float64) for j in js])

met = checked = 0; rms_all = []; exceptions = []; t0 = time.time()
for i in over_idx:
    key = str(d["cond_id"][i]); js = groups.get(key)
    if not js: continue
    C = gather(js); N = len(C)
    if abs(N - d["n_t"][i]) > 0 or N < 6: continue
    mu_c = d["mu_c"][i]
    if checked < 5:                        # basis consistency: memmap centroid must equal the pipeline centroid
        assert np.allclose(C.mean(0), d["mu_t"][i], atol=1e-3), key
    v = C.mean(0) - mu_c; u = v / np.linalg.norm(v)
    n = int(min(max(np.ceil(res["n_star"][i]), 5), N - 1))
    th2 = np.empty(REPS)
    for k in range(REPS):
        vn = C[rng.choice(N, n, replace=False)].mean(0) - mu_c; un = vn / np.linalg.norm(vn)
        th2[k] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
    rms = float(np.sqrt(th2.mean())); rms_all.append(rms); checked += 1
    if rms <= THETA * 1.05: met += 1
    else: exceptions.append(dict(cond=key, n_star=round(float(res["n_star"][i]), 1), N0=int(N), m=round(float(res["m_corr"][i]), 2), rms_deg=round(np.degrees(rms), 2)))
    if checked % 500 == 0:
        print(f"  verified {checked}/{len(over_idx)} (met {met}) {time.time() - t0:.0f}s", flush=True)

rms_all = np.array(rms_all)
out = {"note": "OVER/UNDER downsample-and-measure verification on the predicted-OVER set (Tahoe, shared DMSO "
               "vehicle with ALL DMSO_TF wells per plate pooled, real diagonal Sigma, theta=0.1). Downsample "
               "treated cells to n*, RMS angle vs full-depth direction (same control centroid for both).",
       "generator": "scripts/tahoe_recompute/pass4d_gating_vehicle.py", "reps": REPS,
       f"tahoe_shared_dmso_theta{THETA}": dict(
           predicted_over=int(checked), verified=int(met), pass_fraction=round(met / checked, 4) if checked else None,
           gate="RMS angle <= 1.05*theta", median_rms=round(float(np.median(rms_all)), 4),
           p95_rms=round(float(np.percentile(rms_all, 95)), 4), max_rms=round(float(rms_all.max()), 4),
           exceptions=len(exceptions), exception_list=exceptions,
           median_n_c=int(np.median(d["n_c"])))}
json.dump(out, open(FIX, "w"), indent=1)
print(f"[gate-vehicle] predicted OVER {checked}; met {met} ({100 * met / max(checked, 1):.2f}%); "
      f"median RMS {np.median(rms_all):.4f}; max {rms_all.max():.4f}; wrote {FIX}", flush=True)
