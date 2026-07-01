#!/usr/bin/env python
"""Pass 3 of the Tahoe-100M recompute: turn the per-(drug x line) pseudobulk centroids (Pass 2) into
perturbation magnitudes, cell quotas n*, and the per-cell-line UNDER / Ghost / OVER breakdown.

  centroid(cond)  = sum(cond) / count(cond)                       # HVG log-CP10K space (2000-d)
  coord(cond)     = (centroid - pca_mean) @ PCsᵀ                  # 50-d PCA embedding (Pass 1 basis)
  v(drug,line)    = coord(drug,line) - coord(DMSO_TF, line)       # perturbation vector vs matched vehicle
  m               = ||v||
  n*(drug,line)   = 2 (d-1) sigma^2 / (m^2 theta^2)               # isotropic cell quota, theta=0.1 rad

Classification per condition (baseline N0 = actual cells acquired for that condition):
  OVER   : n* <  N0                 (already saturated; downsampling safe)
  UNDER  : N0 <= n* <= GHOST        (needs more cells)
  Ghost  : n* >  GHOST (=50,000)    (too faint to orient at routine depth)

Aggregated per cell line over all real drugs -> the "cell-line-specific under-sampling" table.
Outputs $OUT/quota_per_condition.csv, $OUT/per_cell_line.csv, $OUT/summary.json (+ printed tables).
"""
import os, json, numpy as np, pandas as pd
import pyarrow.parquet as pq

OUT = os.environ.get("OUT", "/mnt/hdd2/loc-tran/tahoe_work/out")
META = "/mnt/hdd2/loc-tran/tahoe_work/meta/metadata"
THETA = float(os.environ.get("THETA", "0.1"))
GHOST = float(os.environ.get("GHOST", "50000"))
CTRL = "DMSO_TF"

b = np.load(os.path.join(OUT, "basis.npz"))
comps = b["components"].astype(np.float64)          # (50, 2000)
pca_mean = b["pca_mean"].astype(np.float64)         # (2000,)
sigma2 = float(b["sigma2"]); d = int(b["n_comps"])
pb = np.load(os.path.join(OUT, "pseudobulk.npz"), allow_pickle=True)
keys = [str(k) for k in pb["cond_keys"]]
sums = pb["sums"].astype(np.float64); counts = pb["counts"].astype(np.float64)

centroids = sums / counts[:, None]
coords = (centroids - pca_mean) @ comps.T           # (n_cond, 50)
info = {k: (coords[i], counts[i]) for i, k in enumerate(keys)}

# CVCL -> cell_name
clm = pq.read_table(f"{META}/cell_line_metadata.parquet").to_pandas()
cvcl2name = dict(clm.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur")["cell_name"])

# authoritative N0 = post-filter cells per (drug, cell_line) from obs_metadata (Pass 0).
# Prefer the committed fixture (the locked reference vector); fall back to the run output.
_REPO = os.path.join(os.path.dirname(__file__), "..", "..")
_cc_fixture = os.path.join(_REPO, "fixtures", "tahoe_condition_counts.csv")
_cc_path = _cc_fixture if os.path.exists(_cc_fixture) else os.path.join(OUT, "condition_counts.csv")
cc = pd.read_csv(_cc_path)
N0map = {(str(r.drug), str(r.cell_line)): float(r.cells_post) for r in cc.itertuples()}
print(f"[pass3] N0 reference: {_cc_path} ({len(N0map)} conditions)")

const = 2.0 * (d - 1) * sigma2 / THETA**2           # n* = const / m^2
rows = []
for k in keys:
    drug, line = k.rsplit("|", 1)
    if drug == CTRL:
        continue
    ck = f"{CTRL}|{line}"
    if ck not in info:
        continue
    v = info[k][0] - info[ck][0]
    m = float(np.linalg.norm(v))
    if m <= 0:
        continue
    nstar = const / m**2
    N0 = N0map.get((drug, line), float(info[k][1]))   # authoritative obs_metadata count; fallback = stream tally
    reg = "OVER" if nstar < N0 else ("Ghost" if nstar > GHOST else "UNDER")
    rows.append(dict(drug=drug, cvcl=line, cell_line=cvcl2name.get(line, line),
                     m=m, n_star=nstar, N0=N0, regime=reg))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "quota_per_condition.csv"), index=False)

# ---- per cell line ----
def frac(s, r): return 100.0 * (s == r).mean()
pl = []
for cl, g in df.groupby("cell_line"):
    n = len(g)
    pl.append(dict(cell_line=cl, n_drugs=n,
                   pct_OVER=frac(g.regime, "OVER"), pct_UNDER=frac(g.regime, "UNDER"),
                   pct_Ghost=frac(g.regime, "Ghost"),
                   pct_under_or_ghost=100.0 * g.regime.isin(["UNDER", "Ghost"]).mean(),
                   median_m=g.m.median(), median_nstar=g.n_star.median()))
pl = pd.DataFrame(pl).sort_values("pct_under_or_ghost", ascending=False)
pl.to_csv(os.path.join(OUT, "per_cell_line.csv"), index=False)

# ---- overall + case studies ----
overall = dict(
    n_conditions=int(len(df)), n_drugs=int(df.drug.nunique()), n_lines=int(df.cell_line.nunique()),
    sigma2=sigma2, d=d, theta=THETA, nstar_const=const,
    pct_OVER=float(frac(df.regime, "OVER")), pct_UNDER=float(frac(df.regime, "UNDER")),
    pct_Ghost=float(frac(df.regime, "Ghost")),
    median_nstar=float(df.n_star.median()), median_m=float(df.m.median()),
    median_N0=float(df.N0.median()),
)
# per-drug across lines (transpose view): lines OVER-sampled out of the lines present
case = []
for drug, g in df.groupby("drug"):
    case.append(dict(drug=drug, median_m=g.m.median(), median_nstar=g.n_star.median(),
                     n_lines=len(g), lines_OVER=int((g.regime == "OVER").sum())))
case = pd.DataFrame(case).sort_values("median_m", ascending=False)
case.to_csv(os.path.join(OUT, "per_drug.csv"), index=False)
json.dump(overall, open(os.path.join(OUT, "summary.json"), "w"), indent=2)

pd.set_option("display.width", 200); pd.set_option("display.max_rows", 60)
print("=== OVERALL (theta=%.2f rad, d=%d, sigma^2=%.3f, n*=%.0f/m^2) ===" % (THETA, d, sigma2, const))
print(json.dumps(overall, indent=2))
print("\n=== PER CELL LINE (%% of %d drugs) — sorted by under+ghost ===" % df.drug.nunique())
print(pl.round(1).to_string(index=False))
print("\n=== CASE-STUDY DRUGS (by magnitude) ===")
print(case.assign(median_m=case.median_m.round(2), median_nstar=case.median_nstar.round(0)).head(12).to_string(index=False))
print("\n[wrote] quota_per_condition.csv, per_cell_line.csv, per_drug.csv, summary.json in", OUT)
