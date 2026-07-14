#!/usr/bin/env python
"""APPLICATION A: budget-matched recovery of the drug-similarity graph and per-condition direction.

Does allocating cells by the quota n* buy a better scientific output at fixed cost? The output a
direction-based screen produces is (i) each condition's perturbation DIRECTION and (ii) the drug-
similarity graph built from those directions. We measure how faithfully both are recovered under two
sampling policies at a FIXED total budget:

  * FLAT     : every condition receives B / K cells.
  * ADAPTIVE : cap-based water-filling by the quota -- condition c receives min(n*_c, C), with C set so
               the total equals B (cheap/resolvable conditions get exactly what they need; the rest are
               capped). Cells are spent where the direction is actually resolvable.

Two metrics vs the full-depth ground truth:
  - within-tolerance fraction: share of conditions whose estimated direction lies within theta* = 0.1
    rad of its true direction (the quota's own guarantee);
  - similarity-graph fidelity: mean Jaccard overlap of each condition's k-nearest-neighbour set (by
    direction cosine) with its full-depth neighbour set.

Estimated directions use the shipped sampling model u_hat = normalize(v + e), e ~ N(0, 2 Sigma / n),
Sigma = diag(per-PC variance) (the committed global sigma^2). Data: the from-raw Tahoe-100M recompute
(dose = 5 uM, strongest signal). (We also report full-depth k-NN accuracy for the annotated moa-fine
label as an honest aside: it is near chance, because most single-condition directions are weak.)
"""
import os, sys, json, ast, numpy as np, pandas as pd
import pyarrow.parquet as pq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute   # the per-condition quota n* comes from the single engine (two-arm, real n_c)

TW = "/mnt/hdd2/loc-tran/tahoe_work"
OUT = os.path.join(TW, "out_dose"); META = os.path.join(TW, "meta/metadata")
FIXDIR = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
DOSE = 5.0; K_NN = 10; N_DRAWS = 10; MAXNODES = 4000; GHOSTCAP = 50000; THETA = 0.1
rng = np.random.default_rng(0)

b = np.load(os.path.join(OUT, "basis.npz"))
comps = b["components"].astype(np.float64); pca_mean = b["pca_mean"].astype(np.float64)
ppv = b["per_pc_var"].astype(np.float64); d = int(b["n_comps"]); sig = float(b["sigma2"])
Sigma = np.diag(ppv)                                      # global within-condition per-PC variance (diagonal)
pb = np.load(os.path.join(OUT, "pseudobulk.npz"), allow_pickle=True)
keys = [str(k) for k in pb["cond_keys"]]; sums = pb["sums"].astype(np.float64); counts = pb["counts"].astype(np.float64)
sm = pq.read_table(f"{META}/sample_metadata.parquet").to_pandas()
sm["dose"] = sm.drugname_drugconc.map(lambda s: ast.literal_eval(s)[0][1])
s2 = {r.sample: (r.drug, r.dose) for r in sm.itertuples()}
dm = pq.read_table(f"{META}/drug_metadata.parquet").to_pandas().set_index("drug")
moa = dm["moa-fine"].fillna("unclear").to_dict()

cond_sum, cond_cnt = {}, {}; line_sum, line_cnt = {}, {}
for i, k in enumerate(keys):
    s, line = k.rsplit("|", 1); meta = s2.get(s)
    if meta is None:
        continue
    drug, dose = meta
    line_sum[line] = line_sum.get(line, 0) + sums[i]; line_cnt[line] = line_cnt.get(line, 0) + counts[i]
    if drug == "DMSO_TF":
        continue
    key = (drug, dose, line)
    cond_sum[key] = cond_sum.get(key, 0) + sums[i]; cond_cnt[key] = cond_cnt.get(key, 0) + counts[i]

base_pca = {l: (line_sum[l] / line_cnt[l] - pca_mean) @ comps.T for l in line_sum}
rows = []
for (drug, dose, line), sm_ in cond_sum.items():
    if abs(dose - DOSE) > 1e-9:
        continue
    N = cond_cnt[(drug, dose, line)]
    cen = (sm_ / N - pca_mean) @ comps.T
    v = cen - base_pca[line]; m = float(np.linalg.norm(v))
    if m <= 0:
        continue
    # per-condition quota n* from the engine (two-arm; n_c = the per-line control pool); inf -> cap by min() below
    res = compute(cen[None, :], base_pca[line][None, :], np.array([N], float),
                  np.array([float(line_cnt[line])], float), Sigma, theta=THETA)
    rows.append((drug, line, str(moa.get(drug, "unclear")), v, m, float(res["n_star"][0])))
drugs = np.array([r[0] for r in rows]); labels = np.array([r[2] for r in rows])
V = np.stack([r[3] for r in rows]); mags = np.array([r[4] for r in rows]); nstar = np.array([r[5] for r in rows])
if len(rows) > MAXNODES:
    sel = rng.choice(len(rows), MAXNODES, replace=False)
    drugs, labels, V, mags, nstar = drugs[sel], labels[sel], V[sel], mags[sel], nstar[sel]
K = len(V)
U = V / np.linalg.norm(V, axis=1, keepdims=True)
sameinstance = drugs[:, None] == drugs[None, :]
print(f"[app-A] {K} conditions (drug x line at {DOSE} uM); median m = {np.median(mags):.2f}")

# full-depth similarity-graph neighbour sets (exclude same-drug neighbours)
def neighbours(Un):
    S = Un @ Un.T; S[sameinstance] = -2.0; np.fill_diagonal(S, -2.0)
    return np.argpartition(-S, K_NN, axis=1)[:, :K_NN]
NN_full = neighbours(U)
full_sets = [set(r) for r in NN_full]

def moa_acc(Un):
    usable = ~pd.Series(labels).str.lower().isin(["unclear", "none", ""]).values
    S = Un @ Un.T; S[sameinstance] = -2.0; np.fill_diagonal(S, -2.0)
    nn = np.argpartition(-S, 5, axis=1)[:, :5]
    return float((labels[nn][usable] == labels[usable, None]).mean())
moa_ceiling = moa_acc(U)
maj = float(pd.Series(labels[~pd.Series(labels).str.lower().isin(['unclear','none',''])]).value_counts(normalize=True).iloc[0])

sd = np.sqrt(2 * ppv)                                     # per-PC noise sd, before /sqrt(n)

def metrics_at(n_per_cond):
    wt, jac = [], []
    for _ in range(N_DRAWS):
        e = rng.standard_normal((K, d)) * (sd[None, :] / np.sqrt(n_per_cond)[:, None])
        Uh = V + e; Uh /= np.linalg.norm(Uh, axis=1, keepdims=True)
        ang = np.arccos(np.clip((Uh * U).sum(1), -1, 1))
        wt.append(float((ang <= THETA).mean()))
        nn = neighbours(Uh)
        jac.append(float(np.mean([len(full_sets[i] & set(nn[i])) / len(full_sets[i] | set(nn[i])) for i in range(K)])))
    return (float(np.mean(wt)), float(np.std(wt)), float(np.mean(jac)), float(np.std(jac)))

print(f"[app-A] full-depth moa-fine k-NN accuracy = {moa_ceiling:.3f} (majority baseline {maj:.3f}) "
      f"-> single-condition directions do NOT recover fine MoA; the recoverable output is the graph itself")
curve = []
for C in [50, 100, 200, 400, 800, 1600, 3200, 6400, 12800]:
    adap = np.minimum(nstar, C); B = adap.sum(); flat = np.full(K, B / K)
    fwt, fws, fj, fjs = metrics_at(flat); awt, aws, aj, ajs = metrics_at(adap)
    curve.append(dict(cap=C, budget_cells=int(B), cells_per_cond_flat=round(B / K, 1),
                      flat_within_tol=round(fwt, 4), adaptive_within_tol=round(awt, 4),
                      flat_graph_jaccard=round(fj, 4), adaptive_graph_jaccard=round(aj, 4)))
    print(f"  B={B/1e6:5.2f}M ({B/K:6.0f}/cond)  within-tol flat={fwt:.3f} adap={awt:.3f} | "
          f"graph-Jaccard flat={fj:.3f} adap={aj:.3f}")

def savings(metric, target):
    xs = [c["budget_cells"] for c in curve]
    bf = next((xs[i] for i in range(len(curve)) if curve[i]["flat_" + metric] >= target), None)
    ba = next((xs[i] for i in range(len(curve)) if curve[i]["adaptive_" + metric] >= target), None)
    return bf, ba, (round(bf / ba, 2) if bf and ba else None)
bf, ba, sav = savings("graph_jaccard", 0.5)
bf2, ba2, sav2 = savings("within_tol", 0.5)
summary = dict(analysis="Budget-matched recovery of similarity graph and per-condition direction",
               dose_uM=DOSE, n_conditions=K, k_nn=K_NN, theta=THETA,
               moa_fine_full_depth_acc=round(moa_ceiling, 3), moa_majority_baseline=round(maj, 3),
               graph_jaccard_target=0.5, flat_cells_to_graph_target=bf, adaptive_cells_to_graph_target=ba,
               graph_cell_savings_factor=sav,
               within_tol_target=0.5, flat_cells_to_tol_target=bf2, adaptive_cells_to_tol_target=ba2,
               within_tol_cell_savings_factor=sav2, curve=curve)
json.dump(summary, open(os.path.join(FIXDIR, "tahoe_moa_recovery.json"), "w"), indent=1)
print(f"\n[app-A] to recover 50% of the full-depth similarity graph: flat {bf:,} vs adaptive {ba:,} cells "
      f"-> {sav}x fewer with adaptive; wrote fixtures/tahoe_moa_recovery.json")
