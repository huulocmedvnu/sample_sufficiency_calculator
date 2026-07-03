#!/usr/bin/env python
"""APPLICATIONS B & C: reliability audit of the similarity graph, and resource/cost translation.

B. RELIABILITY AUDIT. At the acquired depth, a condition's direction is resolved to tolerance theta*
   only if N0 >= n*(theta*). Since n* scales as 1/theta*^2, we report, across the full 56,827-condition
   Tahoe panel, the fraction of conditions RESOLVED at theta* in {0.1, 0.2, 0.3} rad, and -- on the
   drug-similarity k-NN graph -- the fraction of neighbour edges whose BOTH endpoints are resolved.
   This quantifies how much of a naively-built MoA/similarity graph rests on under-powered directions.

C. RESOURCE / COST TRANSLATION. For a 100-drug x 3-dose x 2-line = 600-condition screen drawing
   magnitudes from the empirical Tahoe distribution, we contrast uniform loading at a "safe" depth
   (90th-percentile n*) against quota-guided allocation (cap 10k), and translate the cell delta into
   sequencing cost at a stated per-cell price. We also report the multiplexing dividend for the
   over-sampled conditions.
"""
import os, json, ast, numpy as np, pandas as pd
import pyarrow.parquet as pq

TW = "/mnt/hdd2/loc-tran/tahoe_work"; OUT = os.path.join(TW, "out_dose"); META = os.path.join(TW, "meta/metadata")
FIXDIR = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
qc = pd.read_csv(os.path.join(FIXDIR, "tahoe_quota_per_condition.csv"))   # drug,dose_uM,cell_line,m,n_star,N0,regime
COST_PER_CELL = float(os.environ.get("COST_PER_CELL", "0.30"))            # all-in 10x + sequencing, USD (assumption)

# ---------- B. reliability audit ----------
ns = qc["n_star"].values; N0 = qc["N0"].values
audit = {}
for th in (0.1, 0.2, 0.3):
    nstar_th = ns * (0.1 / th) ** 2                      # n* scales as 1/theta^2
    resolved = N0 >= nstar_th
    audit[th] = dict(pct_resolved=round(100 * resolved.mean(), 2))
# edge-level reliability on the 5 uM direction graph (both endpoints resolved), k=10
b = np.load(os.path.join(OUT, "basis.npz")); comps = b["components"].astype(float); pca_mean = b["pca_mean"].astype(float)
pb = np.load(os.path.join(OUT, "pseudobulk.npz"), allow_pickle=True)
keys = [str(k) for k in pb["cond_keys"]]; sums = pb["sums"].astype(float); counts = pb["counts"].astype(float)
sm = pq.read_table(f"{META}/sample_metadata.parquet").to_pandas()
sm["dose"] = sm.drugname_drugconc.map(lambda s: ast.literal_eval(s)[0][1]); s2 = {r.sample: (r.drug, r.dose) for r in sm.itertuples()}
cs, cc, ls, lc = {}, {}, {}, {}
for i, k in enumerate(keys):
    s, line = k.rsplit("|", 1); meta = s2.get(s)
    if meta is None: continue
    drug, dose = meta
    ls[line] = ls.get(line, 0) + sums[i]; lc[line] = lc.get(line, 0) + counts[i]
    if drug == "DMSO_TF": continue
    cs[(drug, dose, line)] = cs.get((drug, dose, line), 0) + sums[i]; cc[(drug, dose, line)] = cc.get((drug, dose, line), 0) + counts[i]
bp = {l: (ls[l] / lc[l] - pca_mean) @ comps.T for l in ls}
# Headline is the within-condition residual variance (0.9567); basis.npz b["sigma2"] is the marginal
# 2.406, retained only as a conservative bound. Override with SIGMA2 env to reproduce either.
sig = float(os.environ.get("SIGMA2", 0.9567)); d = int(b["n_comps"]); const = 2 * (d - 1) * sig / 0.01
rows = []
for (drug, dose, line), s_ in cs.items():
    if abs(dose - 5.0) > 1e-9: continue
    v = (s_ / cc[(drug, dose, line)] - pca_mean) @ comps.T - bp[line]; m = np.linalg.norm(v)
    if m > 0: rows.append((drug, v / m, cc[(drug, dose, line)], const / m**2))
U = np.stack([r[1] for r in rows]); dr = np.array([r[0] for r in rows]); N0g = np.array([r[2] for r in rows]); nsg = np.array([r[3] for r in rows])
S = U @ U.T; same = dr[:, None] == dr[None, :]; S[same] = -2; np.fill_diagonal(S, -2)
nn = np.argpartition(-S, 10, axis=1)[:, :10]
for th in (0.1, 0.2, 0.3):
    res = N0g >= nsg * (0.1 / th) ** 2
    both = res[:, None] & res[nn]                        # edge reliable iff both endpoints resolved
    audit[th]["pct_graph_edges_both_resolved"] = round(100 * both.mean(), 2)
    audit[th]["pct_conditions_resolved_5uM"] = round(100 * res.mean(), 2)
print("[B] reliability audit (Tahoe, full panel + 5 uM graph):")
for th, a in audit.items():
    print(f"  theta={th}: {a['pct_resolved']}% of conditions resolved; "
          f"{a['pct_graph_edges_both_resolved']}% of k-NN edges have both endpoints resolved")

# ---------- C. cost / throughput ----------
NC = 600
p90 = float(np.percentile(ns, 90))
flat_cells = p90 * NC
adapt_cells = float(np.minimum(ns, 10000).mean()) * NC
frac_over = float((ns < N0).mean())                      # over-sampled fraction (theta=0.1)
mult = float(np.median((N0 / ns)[ns < N0])) if (ns < N0).any() else 0.0
cost = dict(cost_per_cell_usd=COST_PER_CELL, screen="100 drugs x 3 doses x 2 lines = 600 conditions",
            flat_p90_cells=int(flat_cells), flat_p90_cost_usd=round(flat_cells * COST_PER_CELL),
            adaptive_cap10k_cells=int(adapt_cells), adaptive_cap10k_cost_usd=round(adapt_cells * COST_PER_CELL),
            cell_savings_factor=round(flat_cells / adapt_cells, 2),
            dollars_saved=round((flat_cells - adapt_cells) * COST_PER_CELL),
            pct_conditions_over_sampled=round(100 * frac_over, 2),
            median_multiplex_gain_over=round(mult, 1))
print(f"\n[C] cost (@ ${COST_PER_CELL}/cell): uniform-safe {flat_cells/1e6:.1f}M cells (${cost['flat_p90_cost_usd']:,}) "
      f"vs quota-guided {adapt_cells/1e6:.1f}M (${cost['adaptive_cap10k_cost_usd']:,}) -> "
      f"{cost['cell_savings_factor']}x, ${cost['dollars_saved']:,} saved")
print(f"[C] {cost['pct_conditions_over_sampled']}% of conditions over-sampled; median multiplexing gain there = {cost['median_multiplex_gain_over']}x")

json.dump(dict(reliability_audit=audit, cost=cost), open(os.path.join(FIXDIR, "tahoe_applications.json"), "w"), indent=1)
print("\n[wrote] fixtures/tahoe_applications.json")
