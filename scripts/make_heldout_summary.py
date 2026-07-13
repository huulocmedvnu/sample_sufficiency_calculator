"""
Unified per-dataset held-out falsification summary, each dataset INDEPENDENTLY, on ALL its eligible
groups with N >= 400 cells (no shared-cell-line restriction). One row per dataset:

  n(N>=400) | median slope ratio (all) | in-regime n (rho^2>=3) | median ratio (in-regime) | median R^2
  (in-regime) | Spearman[1-ratio vs 1/rho^2]  (the predicted low-SNR breakdown: deficit grows with 1/rho^2)

The held-out test itself (downsample -> realized theta^2(n) vs (1/n - 1/N), through-origin slope vs the
a-priori tr(PSP)/m^2) is run per group by each atlas's own recompute pass; here we only read the committed
per-group fixtures and aggregate at the N>=400 cut. Sources:
  Tahoe      fixtures/tahoe_direct_curves_full.json   (per_group; min_cells already 400)
  EmeraldBay fixtures/emeraldbay_falsification_full.json (per_group_slope)
  TRADE      fixtures/trade_{jurkat,hepg2}_falsification.json (rows)
  Orion      fixtures/orion_{HCT116,HEK293T}_heldout.json (per_group; from the pass4 re-stream)

Emits fixtures/heldout_summary.json + fixtures/heldout_summary.md.
"""
import os, json, numpy as np
from scipy.stats import spearmanr
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
NMIN = 400


def agg(rows, label):
    """rows: list of dicts with N, rho2, ratio, r2. Filter N>=NMIN, summarize."""
    rows = [r for r in rows if r.get("N", 0) >= NMIN]
    if not rows:
        return dict(dataset=label, n=0)
    ratio = np.array([r["ratio"] for r in rows], float)
    rho2 = np.array([r["rho2"] for r in rows], float)
    r2 = np.array([r["r2"] for r in rows], float)
    inr = rho2 >= 3
    deficit = 1.0 - ratio
    # Spearman[1-ratio, 1/rho^2] is the scale-free universal statistic: rho^2 scales differ hugely
    # across datasets (Tahoe ~1-2, TRADE up to ~28,000), so no fixed rho^2 threshold is "high SNR"
    # uniformly. The rank correlation confirms the monotone low-SNR breakdown on every dataset.
    sp = float(spearmanr(1.0 / np.maximum(rho2, 1e-9), deficit).correlation) if len(rows) > 2 else None
    return dict(dataset=label, n=len(rows),
                median_ratio_all=round(float(np.median(ratio)), 3),
                n_in_regime=int(inr.sum()),
                median_ratio_in_regime=round(float(np.median(ratio[inr])), 3) if inr.any() else None,
                median_r2_in_regime=round(float(np.median(r2[inr])), 4) if inr.any() else None,
                spearman_deficit_vs_inv_rho2=round(sp, 3) if sp is not None else None)


def load(name):
    p = os.path.join(FX, name)
    return json.load(open(p)) if os.path.exists(p) else None


ROWS = []
# Tahoe
t = load("tahoe_direct_curves_full.json")
if t:
    ROWS.append(agg(t["per_group"], "Tahoe-100M"))
# EmeraldBay
e = load("emeraldbay_falsification_full.json")
if e:
    ROWS.append(agg(e["per_group_slope"], "EmeraldBay"))
# Orion (present only after the pass4 re-stream)
for line in ("HCT116", "HEK293T"):
    o = load(f"orion_{line}_heldout.json")
    ROWS.append(agg(o["per_group"], f"Orion {line}") if o else dict(dataset=f"Orion {line}", n=0, pending=True))
# TRADE
for line, disp in (("jurkat", "TRADE Jurkat"), ("hepg2", "TRADE HepG2")):
    tr = load(f"trade_{line}_falsification.json")
    if tr:
        ROWS.append(agg(tr["rows"], disp))

hdr = ("| dataset | groups (N≥400) | median ratio (all) | n(ρ²≥3) | median ratio (ρ²≥3) | median R² (ρ²≥3) | "
       "Spearman deficit vs 1/ρ² |")
sep = "|" + "|".join(["---"] * 7) + "|"
out = [hdr, sep]
for r in ROWS:
    if r.get("pending"):
        out.append(f"| {r['dataset']} | *(re-stream running)* | | | | | |")
    elif r["n"] == 0:
        out.append(f"| {r['dataset']} | 0 (none eligible) | | | | | |")
    else:
        out.append(f"| {r['dataset']} | {r['n']:,} | {r['median_ratio_all']} | {r['n_in_regime']} | "
                   f"{r['median_ratio_in_regime']} | {r['median_r2_in_regime']} | "
                   f"{r['spearman_deficit_vs_inv_rho2']} |")
md = "\n".join(out)
print(md)
open(os.path.join(FX, "heldout_summary.md"), "w").write(md + "\n")
json.dump(dict(min_cells=NMIN, rows=ROWS), open(os.path.join(FX, "heldout_summary.json"), "w"), indent=1)
print("\nwrote fixtures/heldout_summary.{md,json}")
