#!/usr/bin/env python
"""Pass 3 (dose-resolved) — turn per-(sample x line) pseudobulk centroids into per-(drug x dose x line)
cell quotas and the OVER / UNDER / Ghost spectrum.

The perturbation CONDITION is (drug x dose x cell_line). `sample` identifies drug x dose x plate; we
pool plate-replicates of the same drug+dose and subtract the PLATE-MATCHED DMSO_TF vehicle (per plate
x line), so the effect vector is batch-referenced within plate.

  coord(sample,line) = (centroid - pca_mean) @ PCsᵀ                       # 50-d PCA (Pass 1 basis)
  v(drug,dose,line)  = Σ_p c_p·coord(sample_p,line)/Σc_p  −  Σ_p c_p·coord(DMSO@plate_p,line)/Σc_p
  m = ||v||,  n* = 2(d-1)σ²/(m²θ²) = 23,577/m²  (θ=0.1),  N0 = Σ_p c_p (cells for that condition)
  OVER n*<N0 ; UNDER N0<=n*<=50k ; Ghost n*>50k

Outputs $OUT/{quota_per_condition,per_cell_line,per_drug,per_dose}.csv, summary.json.
"""
import os, json, ast, numpy as np, pandas as pd, pyarrow.parquet as pq

OUT = os.environ.get("OUT", "/mnt/hdd2/loc-tran/tahoe_work/out_dose")
META = "/mnt/hdd2/loc-tran/tahoe_work/meta/metadata"
THETA = float(os.environ.get("THETA", "0.1")); GHOST = float(os.environ.get("GHOST", "50000"))
CTRL = "DMSO_TF"

b = np.load(os.path.join(OUT, "basis.npz"))
comps = b["components"].astype(np.float64); pca_mean = b["pca_mean"].astype(np.float64)
sig = float(b["sigma2"]); d = int(b["n_comps"]); const = 2 * (d - 1) * sig / THETA**2

pb = np.load(os.path.join(OUT, "pseudobulk.npz"), allow_pickle=True)
keys = [str(k) for k in pb["cond_keys"]]; sums = pb["sums"].astype(np.float64); counts = pb["counts"].astype(np.float64)
coords = (sums / counts[:, None] - pca_mean) @ comps.T          # (n_keys, 50)

# sample -> (drug, dose, plate) and CVCL -> cell_name
sm = pq.read_table(f"{META}/sample_metadata.parquet").to_pandas()
def parse_dose(s):
    try: return float(ast.literal_eval(s)[0][1])
    except Exception: return np.nan
sm["dose"] = sm.drugname_drugconc.map(parse_dose)
S2meta = {r.sample: (r.drug, r.dose, r.plate) for r in sm.itertuples()}
clm = pq.read_table(f"{META}/cell_line_metadata.parquet").to_pandas()
cvcl2name = dict(clm.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur")["cell_name"])

# index coords by (sample,line); collect DMSO per (plate,line) and treatment per (drug,dose,line)
dmso = {}                                     # (plate, line) -> coord
cond = {}                                     # (drug, dose, line) -> list[(coord, count, plate)]
for i, k in enumerate(keys):
    samp, line = k.rsplit("|", 1)
    meta = S2meta.get(samp)
    if meta is None:
        continue
    drug, dose, plate = meta
    if drug == CTRL:
        dmso[(plate, line)] = coords[i]
    else:
        cond.setdefault((drug, dose, line), []).append((coords[i], counts[i], plate))

rows = []
for (drug, dose, line), items in cond.items():
    tot = sum(c for _, c, _ in items)
    if tot <= 0:
        continue
    treat = sum(c * co for co, c, _ in items) / tot
    parts = [c * dmso[(p, line)] for co, c, p in items if (p, line) in dmso]
    if len(parts) != len(items):                # require plate-matched control for every part
        continue
    ctrl = sum(parts) / tot
    m = float(np.linalg.norm(treat - ctrl))
    if m <= 0:
        continue
    nstar = const / m**2
    reg = "OVER" if nstar < tot else ("Ghost" if nstar > GHOST else "UNDER")
    rows.append(dict(drug=drug, dose_uM=dose, cvcl=line, cell_line=cvcl2name.get(line, line),
                     m=m, n_star=nstar, N0=tot, regime=reg))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "quota_per_condition.csv"), index=False)


def frac(s, r): return 100.0 * (s == r).mean()
def spectrum(g):
    return dict(n=len(g), pct_OVER=round(frac(g.regime, "OVER"), 1), pct_UNDER=round(frac(g.regime, "UNDER"), 1),
                pct_Ghost=round(frac(g.regime, "Ghost"), 1),
                pct_under_or_ghost=round(100 * g.regime.isin(["UNDER", "Ghost"]).mean(), 1),
                median_m=round(g.m.median(), 2), median_nstar=int(g.n_star.median()))

per_line = pd.DataFrame([{**{"cell_line": cl}, **spectrum(g)} for cl, g in df.groupby("cell_line")]
                        ).sort_values("pct_under_or_ghost", ascending=False)
per_line.to_csv(os.path.join(OUT, "per_cell_line.csv"), index=False)
per_dose = pd.DataFrame([{**{"dose_uM": dz}, **spectrum(g)} for dz, g in df.groupby("dose_uM")]).sort_values("dose_uM")
per_dose.to_csv(os.path.join(OUT, "per_dose.csv"), index=False)
# per drug: pooled over doses+lines
pd_rows = []
for drug, g in df.groupby("drug"):
    pd_rows.append(dict(drug=drug, median_m=round(g.m.median(), 2), median_nstar=int(g.n_star.median()),
                        n_conditions=len(g), OVER=int((g.regime == "OVER").sum()),
                        UNDER=int((g.regime == "UNDER").sum()), Ghost=int((g.regime == "Ghost").sum())))
per_drug = pd.DataFrame(pd_rows).sort_values("OVER", ascending=False)
per_drug.to_csv(os.path.join(OUT, "per_drug.csv"), index=False)

overall = dict(n_conditions=int(len(df)), n_drugs=int(df.drug.nunique()),
               n_doses=int(df.dose_uM.nunique()), n_lines=int(df.cell_line.nunique()),
               sigma2=sig, d=d, theta=THETA, nstar_const=round(const, 1),
               **{k: v for k, v in spectrum(df).items()}, median_N0=int(df.N0.median()))
json.dump(overall, open(os.path.join(OUT, "summary.json"), "w"), indent=2)

pd.set_option("display.width", 200)
print("=== OVERALL (drug x dose x line; theta=%.2f, sigma^2=%.3f, n*=%.0f/m^2) ===" % (THETA, sig, const))
print(json.dumps(overall, indent=2))
print("\n=== PER DOSE ===\n", per_dose.to_string(index=False))
print("\n=== PER CELL LINE (worst 8 / best 4) ===")
print(pd.concat([per_line.head(8), per_line.tail(4)]).to_string(index=False))
print("\n=== PER DRUG (strongest 6 / weakest 4) ===")
print(pd.concat([per_drug.head(6), per_drug.tail(4)]).to_string(index=False))
print(f"\ndrugs OVER in 0 conditions: {(per_drug.OVER==0).sum()}/{len(per_drug)}")
print("[wrote]", OUT)
