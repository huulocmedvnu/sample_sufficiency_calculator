#!/usr/bin/env python
"""Pass 6 (EmeraldBay): downsampling-invariance exemplars (manuscript Table 4 + INVARIANCE.md).

Regenerates the "downsampling stability tracks the quota diagnosis" numbers on the CURRENT committed
EmeraldBay basis (now the frozen-HVG basis). Closes a reproducibility gap: the original generator for
outputs/emeraldbay_invariance.json was not committed.

Takes one OVER-sampled and one UNDER-sampled (condition x line) group at theta=0.1 rad, downsamples both
by R random replicates to the OVER group's quota n* = 2 tr(P Sigma P)/(m^2 theta^2), and reports, versus
full depth:
  * Euclidean distance of the condition centroid from the per-line baseline (relative drift);
  * cosine similarity of its direction to its 3 nearest neighbours in the drug-similarity graph;
  * replicate SD of both.
Drug-similarity graph = cosine between condition direction vectors (centroid - per-line mean), built over
ALL 52-line (condition x line) groups with >= MIN_CELLS cells (per-cell coords from all_cells.npz), so the
nearest neighbours are the true nearest among the whole atlas, not just the 5 Tahoe-shared lines. Set
EB_COORDS to another .npz (e.g. shared_cells.npz) to restrict the graph.

Writes outputs/emeraldbay_invariance.json. Reads eb_work/out (basis, all_cells, sample2cond).
"""
import os, json, ast, numpy as np

OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
META = "/mnt/hdd2/loc-tran/eb_work/meta/metadata"
COORDS_NPZ = os.environ.get("EB_COORDS", os.path.join(OUT, "all_cells.npz"))  # full 52-line atlas
HERE = os.path.dirname(__file__)
FIGDIR = os.path.join(HERE, "..", "..", "outputs"); os.makedirs(FIGDIR, exist_ok=True)
OUTJSON = os.path.join(FIGDIR, "emeraldbay_invariance.json")
THETA = 0.1
REPS = 300
MIN_CELLS = 100
rng = np.random.default_rng(0)

import pyarrow.parquet as pq
b = np.load(os.path.join(OUT, "basis.npz")); d = int(b["n_comps"])
sh = np.load(COORDS_NPZ, allow_pickle=True)
coords = sh["coords"].astype(float); samp = sh["sample"].astype(str); cline = sh["line"].astype(str)
s2cond = json.load(open(os.path.join(OUT, "sample2cond.json")))
clm = pq.read_table(f"{META}/cell_line_metadata.parquet").to_pandas()
cvcl2name = dict(clm.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur")["cell_name"])


def cond_name(s):
    try:
        return str(ast.literal_eval(s2cond[s])[0][0])
    except Exception:
        return s2cond.get(s, "?")


cond_cell = np.array([cond_name(s) for s in samp])
baseline = {l: coords[cline == l].mean(0) for l in np.unique(cline)}

# --- build the drug-similarity graph over shared-line (condition x line) groups ---
groups = {}                                    # key -> (mask, direction unit vector, m)
for l in np.unique(cline):
    for cd in np.unique(cond_cell[cline == l]):
        mask = (cline == l) & (cond_cell == cd); N = int(mask.sum())
        if N < MIN_CELLS:
            continue
        v = coords[mask].mean(0) - baseline[l]; m = float(np.linalg.norm(v))
        if m <= 0:
            continue
        key = f"{cd}|{cvcl2name.get(l, l)}"
        groups[key] = dict(mask=mask, u=v / m, m=m, line=l, N=N)
keys = list(groups)
U = np.stack([groups[k]["u"] for k in keys])   # (G, d) unit directions
S = U @ U.T                                     # cosine similarity graph


def nn3(key):
    i = keys.index(key)
    order = np.argsort(-S[i]); order = order[order != i]
    return [keys[j] for j in order[:3]]


def quota(key):
    g = groups[key]; C = coords[g["mask"]]
    Sig = np.cov(C.T); u = g["u"]
    P = np.eye(d) - np.outer(u, u)
    trPSP = float(np.trace(P @ Sig @ P))
    return 2.0 * trPSP / (g["m"] ** 2 * THETA ** 2), trPSP


def downsample_metrics(key, n_d, neighbours):
    g = groups[key]; C = coords[g["mask"]]; base = baseline[g["line"]]; N = len(C)
    full_v = C.mean(0) - base; full_dist = float(np.linalg.norm(full_v)); full_u = full_v / full_dist
    full_cos = {nb: float(full_u @ groups[nb]["u"]) for nb in neighbours}
    dists, cos_by_nb = [], {nb: [] for nb in neighbours}
    for _ in range(REPS):
        idx = rng.choice(N, min(n_d, N - 1), replace=False)
        v = C[idx].mean(0) - base; dd = float(np.linalg.norm(v)); un = v / dd
        dists.append(dd)
        for nb in neighbours:
            cos_by_nb[nb].append(float(un @ groups[nb]["u"]))
    dists = np.array(dists)
    return dict(
        N0=N, n_d=int(n_d), full_dist=full_dist,
        dist_rel_drift=float(abs(dists.mean() - full_dist) / full_dist),
        dist_sd=float(dists.std()),
        neighbours=neighbours,
        cos_full={nb: round(full_cos[nb], 4) for nb in neighbours},
        cos_drift={nb: round(float(abs(np.mean(cos_by_nb[nb]) - full_cos[nb])), 4) for nb in neighbours},
        cos_sd=round(float(np.mean([np.std(cos_by_nb[nb]) for nb in neighbours])), 4),
    )


OVER = os.environ.get("EB_OVER", "DMSO_T0|HS-578T")
UNDER = os.environ.get("EB_UNDER", "Galunisertib|AN3-CA")
for req in (OVER, UNDER):
    if req not in groups:
        raise SystemExit(f"group '{req}' not found among {len(keys)} shared-line groups; "
                         f"available drug examples e.g. {[k for k in keys if 'DMSO' not in k][:8]}")

nstar_over, trPSP_over = quota(OVER); nstar_under, _ = quota(UNDER)
n_d = int(round(nstar_over))                            # common depth = the OVER group's own quota
nn_over, nn_under = nn3(OVER), nn3(UNDER)               # each group vs ITS OWN 3 nearest neighbours
res_over = downsample_metrics(OVER, n_d, nn_over)
res_under = downsample_metrics(UNDER, n_d, nn_under)

summary = dict(
    theta=THETA, reps=REPS, basis="frozen-HVG (committed)",
    graph_scope=dict(source=os.path.basename(COORDS_NPZ), n_lines=int(np.unique(cline).size),
                     n_groups=len(keys), n_cells=int(len(coords))),
    over=dict(key=OVER, m=round(groups[OVER]["m"], 3), nstar=int(round(nstar_over)), **res_over),
    under=dict(key=UNDER, m=round(groups[UNDER]["m"], 3), nstar=int(round(nstar_under)), **res_under),
)
json.dump(summary, open(OUTJSON, "w"), indent=1)

print(f"[pass6-inv] graph: {len(keys)} groups over {np.unique(cline).size} lines "
      f"({len(coords):,} cells, {os.path.basename(COORDS_NPZ)})")
print(f"[pass6-inv] OVER  {OVER}: N0={res_over['N0']} m={groups[OVER]['m']:.2f} n*={int(round(nstar_over))}")
print(f"           dist rel drift={res_over['dist_rel_drift']*100:.2f}%  cos drift<= {max(res_over['cos_drift'].values()):.4f}  "
      f"cos SD={res_over['cos_sd']}  dist SD={res_over['dist_sd']:.3f}")
print(f"[pass6-inv] UNDER {UNDER}: N0={res_under['N0']} m={groups[UNDER]['m']:.2f} n*={int(round(nstar_under))} "
      f"({nstar_under/res_under['N0']:.0f}x short)")
print(f"           dist rel drift={res_under['dist_rel_drift']*100:.1f}%  cos drift<= {max(res_under['cos_drift'].values()):.3f}  "
      f"cos SD={res_under['cos_sd']}  dist SD={res_under['dist_sd']:.3f}")
print(f"[pass6-inv] neighbours of OVER: {nn_over}")
print(f"[pass6-inv] wrote {OUTJSON}")
