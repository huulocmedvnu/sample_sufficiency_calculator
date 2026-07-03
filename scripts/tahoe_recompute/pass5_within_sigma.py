#!/usr/bin/env python
"""Compute the Tahoe-100M WITHIN-condition per-cell variance (the plug-in the theory calls for),
and re-tally the sample-sufficiency spectrum with the rescaled quota. Acts on the manuscript TODO:
the headline used the marginal sigma^2=2.406; this reports the within-condition value and the
resulting spectrum. Reads the per-cell memmap from pass4a_stream.py. Writes fixtures/tahoe_within_sigma.json.
"""
import numpy as np, csv, json, os
BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"
z = np.load(BASE + "/ckpt.npz", allow_pickle=True)
cursor = int(z["cursor"]); NCOND = len(z["cond_keys"]); NC = 50
coords = np.memmap(BASE + "/coords.f32", dtype=np.float32, mode="r", shape=(cursor, NC))
condid = np.memmap(BASE + "/condid.i32", dtype=np.int32, mode="r", shape=(cursor,))

csum = np.zeros((NCOND, NC)); csq = np.zeros((NCOND, NC)); cnt = np.zeros(NCOND)
CH = 4_000_000
for s in range(0, cursor, CH):
    cc = np.asarray(condid[s:s + CH]); X = np.asarray(coords[s:s + CH], dtype=np.float64)
    cnt += np.bincount(cc, minlength=NCOND)
    np.add.at(csum, cc, X)
    np.add.at(csq, cc, X * X)
    print(f"  ..{min(s+CH,cursor):,}/{cursor:,}", flush=True)

ok = cnt >= 2
ss = (csq[ok] - csum[ok] ** 2 / cnt[ok, None]).sum(0)
dof = (cnt[ok] - 1).sum()
per_pc = ss / dof
sig2_within = float(per_pc.mean())

d = 50; th = 0.1
const_new = 2 * (d - 1) * sig2_within / th ** 2
rows = list(csv.DictReader(open("fixtures/tahoe_quota_per_condition.csv")))
m = np.array([float(r["m"]) for r in rows]); N0 = np.array([float(r["N0"]) for r in rows])
nstar = const_new / m ** 2
over = float((nstar < N0).mean() * 100); ghost = float((nstar > 50000).mean() * 100)
under = 100 - over - ghost
out = dict(sigma2_within=round(sig2_within, 4), sigma2_marginal=2.406,
           scale_factor=round(sig2_within / 2.406, 4), n_cells=int(cursor), n_conditions_used=int(ok.sum()),
           quota_const_new=round(const_new, 1), quota_const_old=23577,
           per_pc_min=round(float(per_pc.min()), 3), per_pc_max=round(float(per_pc.max()), 3),
           pct_over=round(over, 1), pct_under=round(under, 1), pct_ghost=round(ghost, 1),
           median_nstar=int(np.median(nstar)), boundary_m=round(float(np.sqrt(const_new / np.median(N0))), 2),
           note="within-condition per-cell variance over all cells (sample x line groups); "
                "m and N0 unchanged, n* rescaled by sigma2_within/2.406")
json.dump(out, open("fixtures/tahoe_within_sigma.json", "w"), indent=1)
print("\n=== RESULT ===")
print(f"within-condition sigma^2 = {sig2_within:.4f}  (marginal 2.406; scale {sig2_within/2.406:.3f})")
print(f"new quota constant = {const_new:,.0f} / m^2   (was 23,577)")
print(f"spectrum: {over:.1f}% OVER / {under:.1f}% UNDER / {ghost:.1f}% Ghost")
print(f"median n* = {np.median(nstar):,.0f}  (was 14,570); boundary m = {np.sqrt(const_new/np.median(N0)):.2f} (was 4.27)")
print("wrote fixtures/tahoe_within_sigma.json")
