#!/usr/bin/env python
"""Dataset Inventory & Scale Audit generator (platform-agnostic).

Builds the README "## II. Dataset Inventory & Scale Audit" tables for the two public reference
perturbation atlases used for empirical validation (Tahoe-100M and EmeraldBay). GLOBAL figures are
dataset-level (published/metadata); CAPTURED per-well distributions are computed from the cached arrays
produced by the companion calibration pipeline (default $DATA_DIR). Cell-line identities are anonymized
to "Representative Line 1-5" in the public tables; real identifiers are retained in the committed
calibration fixtures for reproducibility.

  python scripts/dataset_inventory.py            # uses $DATA_DIR (default below)
"""
import io, glob, os, sys
import numpy as np, pandas as pd
import pyarrow.parquet as pq

DATA_DIR = os.environ.get("DATA_DIR", "data")   # companion-pipeline cache
README = os.path.join(os.path.dirname(__file__), "..", "README.md")
# five established cancer lines common to both atlases -> anonymized labels (identities kept in fixtures)
REF_BY_CVCL = {"CVCL_0332": "Representative Line 1", "CVCL_0028": "Representative Line 2",
               "CVCL_0293": "Representative Line 3", "CVCL_0179": "Representative Line 4",
               "CVCL_1094": "Representative Line 5"}
REF_BY_NAME = {"HS-578T": "Representative Line 1", "AN3-CA": "Representative Line 2",
               "HEC-1-A": "Representative Line 3", "BT-474": "Representative Line 4",
               "C-33 A": "Representative Line 5"}


def dist(c):
    c = np.asarray(c)
    return dict(n=len(c), total=int(c.sum()), min=int(c.min()), median=int(np.median(c)),
                mean=int(round(c.mean())), max=int(c.max()))


# ---- Atlas A (Tahoe-100M): captured = the five dose-matched plates 3/6/9/12/13 ----
grp = {}
for ck in glob.glob(f"{DATA_DIR}/drugsim_cache/plate*.npz"):
    z = np.load(ck, allow_pickle=True); p = int(z["plate"])
    for d, l, c in zip(z["drugs"], z["lines"], z["counts"]):
        grp[(str(d), str(l), p)] = int(c)
A = dist(list(grp.values()))
# per-condition replication R (distinct plates per drug-line); pooling does not change cells when R=1
from collections import Counter as _C
_perdl = _C((str_d, str_l) for (str_d, str_l, _p) in grp if str_d != "DMSO_TF")
import numpy as _np
_R = _np.array(list(_perdl.values()))
R1_pct = 100.0 * (_R == 1).mean(); Rmax = int(_R.max())
A_x = []
for cv, lab in REF_BY_CVCL.items():
    veh = [v for k, v in grp.items() if k[1] == cv and k[0] == "DMSO_TF"]
    act = [v for k, v in grp.items() if k[1] == cv and k[0] != "DMSO_TF"]
    A_x.append((lab, len(veh), sum(veh), len(act), sum(act)))

# ---- Atlas B (EmeraldBay): global from metadata; captured = the five representative lines ----
from huggingface_hub import HfFileSystem
fs = HfFileSystem(); REPO = "datasets/tahoebio/EmeraldBay"
dmeta = pq.read_table(io.BytesIO(fs.cat_file(f"{REPO}/metadata/drug_metadata.parquet"))).to_pandas()
ss = pq.read_table(io.BytesIO(fs.cat_file(f"{REPO}/metadata/summary_statistics.parquet"))).to_pandas()
B_lines, B_wells, B_cond = ss.cell_line.nunique(), len(ss), ss.condition.nunique()
z = np.load(f"{DATA_DIR}/emeraldbay_gyn.npz", allow_pickle=True)
ebd = pd.DataFrame({"drug": z["drugs"], "line": z["lines"], "dose": z["doses"]})
ebd["well"] = ebd.drug.astype(str) + "|" + ebd.line.astype(str) + "|" + ebd.dose.astype(str)
B = dist(ebd.groupby("well").size().values)
B_x = []
for name, lab in REF_BY_NAME.items():
    sub = ebd[ebd.line == name]; veh = sub[sub.drug.str.contains("DMSO")]; act = sub[~sub.drug.str.contains("DMSO")]
    B_x.append((lab, veh.well.nunique(), len(veh), act.well.nunique(), len(act)))
A_x.sort(); B_x.sort()

M = ["## II. Dataset Inventory & Scale Audit\n",
     "*Scale ledger of two large-scale, independent, multi-line reference perturbation atlases used "
     "purely for empirical validation. **GLOBAL** rows are dataset-level (published/metadata); "
     "**CAPTURED** rows are computed from the cached arrays of the companion calibration pipeline. "
     "Per-well cell-count distributions are reported for CAPTURED data only, since global per-well "
     "counts are not in either atlas's metadata. Cell-line identities are anonymized to Representative "
     "Lines 1-5 (five established cancer lines common to both atlases; identities retained in the "
     "committed fixtures).*\n",
     "### A. Global inventory\n",
     "| Reference atlas | Total cells (global) | Unique perturbagens | Unique cell lines | Wells (line × condition) |",
     "|---|---:|---:|---:|---:|",
     "| **Atlas A (Tahoe-100M)** | ~100,000,000 | 379 | 50 | ~56,850 |",
     f"| **Atlas B (EmeraldBay)** | ~1,831,756 | {len(dmeta)} molecules ({B_cond} conditions) | {B_lines} | {B_wells:,} |",
     "",
     "### B. Per-well cell-count distribution (CAPTURED data only)\n",
     "| Atlas (captured scope) | Wells | Cells | Min | Median | Mean | Max |",
     "|---|---:|---:|---:|---:|---:|---:|",
     f"| Atlas A — 5 dose-matched plates | {A['n']:,} | {A['total']:,} | {A['min']} | {A['median']:,} | {A['mean']:,} | {A['max']:,} |",
     f"| Atlas B — 5 representative lines | {B['n']:,} | {B['total']:,} | {B['min']} | {B['median']} | {B['mean']} | {B['max']:,} |",
     "",
     f"*Atlas A captured median **{A['median']:,}** cells/well sits far below the median quota "
     "n★≈23,934 required at θ=0.1 rad — i.e. most wells are under-sampled at tight tolerance "
     "(see `docs/SCALE_AUDIT.md`). Atlas A conditions are essentially unreplicated (R=1 for {:.0f}% of drug-line conditions, max R={}), so pooling replicate wells does not change this gating — only ~0.08% of conditions cross UNDER->OVER when pooled (docs/SCALE_AUDIT.md sec 5).*\n".format(R1_pct, Rmax),
     "### C. Cell-line cross-tabulation — vehicle vs active perturbations (CAPTURED)\n",
     "| Cell line | Atlas | Vehicle wells | Vehicle cells | Active wells | Active cells | Total cells |",
     "|---|---|---:|---:|---:|---:|---:|"]
for lab, vw, vc, aw, ac in A_x:
    M.append(f"| {lab} | Atlas A | {vw} | {vc:,} | {aw} | {ac:,} | {vc+ac:,} |")
for lab, vw, vc, aw, ac in B_x:
    M.append(f"| {lab} | Atlas B | {vw} | {vc:,} | {aw} | {ac:,} | {vc+ac:,} |")
M.append("")
md = "\n".join(M)
print(md)

txt = open(README).read()
mk = "## II. Dataset Inventory & Scale Audit"
if mk in txt:
    txt = txt[:txt.index(mk)].rstrip() + "\n\n"
open(README, "w").write(txt.rstrip() + "\n\n" + md + "\n")
print(f"\n[done] regenerated README section II (generic labels)", file=sys.stderr)
import os as _o; _o._exit(0)
