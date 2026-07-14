#!/usr/bin/env python
"""APPLICATIONS B & C: reliability audit of the similarity graph, and resource/cost translation.

All sufficiency numbers come from the SINGLE engine (src/engine.py) via the chemical pipeline; this
script assembles no quota math of its own and needs no atlas cache (the pipeline reads the committed
per-condition sufficient statistics). It regenerates fixtures/tahoe_applications.json.

B. RELIABILITY AUDIT. At the acquired depth, a condition's direction is resolved to tolerance theta*
   only if N0 >= n*(theta*) and the shared vehicle itself resolves it (not control-pool-limited); since
   n* scales as 1/theta*^2 we report, across the full 56,827-condition Tahoe panel (real shared DMSO
   vehicle), the fraction RESOLVED (= OVER) at theta* in {0.1, 0.2, 0.3} rad, and -- on the drug-
   similarity k-NN graph -- the fraction of neighbour edges whose BOTH endpoints are resolved.

C. RESOURCE / COST TRANSLATION. For a 100-drug x 3-dose x 2-line = 600-condition screen drawing
   magnitudes from the Tahoe distribution, we contrast uniform loading at a "safe" depth (90th-
   percentile n*) against quota-guided allocation (cap 10k) under the per-line-mean reference (a large
   pool, so every quota is finite and treated allocation is the lever), and translate the cell delta
   into sequencing cost at a stated per-cell price.
"""
import os, sys, json, numpy as np
from sklearn.neighbors import NearestNeighbors
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute
from pipelines import chemical

FIXDIR = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
COST_PER_CELL = 0.3
NC = 600
CAP = 10000

# --- B: reliability audit against the real shared DMSO vehicle (headline) ---
sp = chemical.load("Tahoe-100M", control="shared-dmso")
v = sp["mu_t"] - sp["mu_c"]; u = v / np.linalg.norm(v, axis=1, keepdims=True)
_, idx = NearestNeighbors(n_neighbors=11).fit(u).kneighbors(u)
N = len(u); src = np.repeat(np.arange(N), 10); dst = idx[:, 1:].ravel()
rel = {}
for th in (0.1, 0.2, 0.3):
    res = compute(sp["mu_t"], sp["mu_c"], sp["n_t"], sp["n_c"], sp["Sigma"], theta=th)["regime"] == "OVER"
    rel[f"{th}"] = dict(pct_resolved=round(100 * float(res.mean()), 1),
                        pct_graph_edges_both_resolved=round(100 * float((res[src] & res[dst]).mean()), 1))

# --- C: budget under the per-line-mean reference (large pool -> finite quotas) ---
spL = chemical.load("Tahoe-100M", control="per-line-mean")
nsL = compute(spL["mu_t"], spL["mu_c"], spL["n_t"], spL["n_c"], spL["Sigma"])["n_star"]
fin = nsL[np.isfinite(nsL)]
p90 = float(np.percentile(fin, 90))
flat = p90 * NC; adapt = float(np.minimum(fin, CAP).mean()) * NC
cost = dict(cost_per_cell_usd=COST_PER_CELL, screen="100 drugs x 3 doses x 2 lines = 600 conditions",
            control="per-line mean (large pool), within-condition sigma^2",
            flat_p90_quota=round(p90), flat_p90_cells=int(flat), flat_p90_cost_usd=round(flat * COST_PER_CELL),
            adaptive_cap10k_cells=int(adapt), adaptive_cap10k_cost_usd=round(adapt * COST_PER_CELL),
            cell_savings_factor=round(flat / adapt, 2), flat_M=round(flat / 1e6, 1), adaptive_M=round(adapt / 1e6, 1))

out = dict(reliability_audit=rel, cost=cost)
json.dump(out, open(os.path.join(FIXDIR, "tahoe_applications.json"), "w"), indent=1)

for th in ("0.1", "0.2", "0.3"):
    print(f"[B] theta={th}: resolved {rel[th]['pct_resolved']}%  graph-edges both-resolved {rel[th]['pct_graph_edges_both_resolved']}%")
print(f"[C] 90th-pct quota {cost['flat_p90_quota']:,}; uniform {cost['flat_M']}M vs quota-guided {cost['adaptive_M']}M "
      f"-> {cost['cell_savings_factor']}x reduction")
print("wrote fixtures/tahoe_applications.json")
