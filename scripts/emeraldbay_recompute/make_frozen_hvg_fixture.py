#!/usr/bin/env python
"""Materialize the EmeraldBay FROZEN HVG set as a committed fixture.

The frozen HVG(2000) is the seurat HVG selection computed on a 20-shard consensus pool (the union of
the committed 10-shard fit set and a disjoint 10-shard set, ~315k cells -- twice the cells of the
committed per-fit selection, so the ranking near the 2000-gene cutoff is far more stable). This is the
set that raises the EmeraldBay basis reproducibility from mean cos^2 = 0.79 (per-fit HVG) to ~0.92
(shared HVG); see scripts/confirm_frozen_hvg.py which computed it.

Source : eb_work/out/frozen_hvg_confirm.json  (key "frozen_hvg_ids")
Output : fixtures/emeraldbay_frozen_hvg.json   (JSON list of 2000 gene token ids, sorted)
"""
import os, json, numpy as np

SRC = "/mnt/hdd2/loc-tran/eb_work/out/frozen_hvg_confirm.json"
HERE = os.path.dirname(__file__)
DST = os.path.join(HERE, "..", "..", "fixtures", "emeraldbay_frozen_hvg.json")

s = json.load(open(SRC))
ids = np.unique(np.asarray(s["frozen_hvg_ids"], dtype=np.int64))
assert ids.size == 2000, f"expected 2000 HVGs, got {ids.size}"
json.dump(ids.tolist(), open(DST, "w"))
print(f"[frozen-hvg] wrote {ids.size} token ids -> {DST}")
print(f"[frozen-hvg] provenance: 20-shard consensus (A={s['A_shards']}, B={s['B_shards']}), "
      f"{s['A_cells']+s['B_cells']:,} cells; baseline cos^2={s['baseline_cos2']:.3f} -> "
      f"frozen cos^2={s['frozen_cos2']:.3f}")
