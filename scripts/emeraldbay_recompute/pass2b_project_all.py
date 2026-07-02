#!/usr/bin/env python
"""Full-atlas EmeraldBay projection: re-stream ALL 116 shards and retain per-cell PCA(50) coordinates
for ALL 52 cell lines (not just the 5 shared with Tahoe), so the falsification test (pass4) can run on
the entire 1.83M-cell atlas rather than the 141,720-cell shared-line slice.

Writes one file per shard to $ALLDIR/shard_XXX.npz (coords, sample, line) -> naturally RESUMABLE
(a shard whose file already exists is skipped, so an interrupted 57.7 GB stream resumes for free).
After all shards, concatenates into $OUT_EB/all_cells.npz and regenerates sample2cond.json.

Run: python scripts/emeraldbay_recompute/pass2b_project_all.py
"""
import os, json, time, numpy as np, scipy.sparse as sp
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

REPO = "tahoebio/EmeraldBay"
CACHE = "/mnt/hdd2/loc-tran/eb_work/dl2"
OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
ALLDIR = os.path.join(OUT, "allcells"); os.makedirs(ALLDIR, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
N_TOTAL = 116

b = np.load(os.path.join(OUT, "basis.npz"))
HVG = b["hvg_token_ids"].astype(np.int64); GW = int(b["gene_width"]); TS = float(b["target_sum"])
NH = HVG.size; comps = b["components"].astype(np.float64); pca_mean = b["pca_mean"].astype(np.float64)
c0 = pca_mean @ comps.T; hvg_map = np.full(GW, -1, dtype=np.int64); hvg_map[HVG] = np.arange(NH)
s2cond = {}


def project_shard(path):
    t = pq.read_table(path, columns=["genes", "expressions", "sample", "cell_line", "drugname_drugconc"])
    for s, c in zip(t.column("sample").to_pylist(), t.column("drugname_drugconc").to_pylist()):
        s2cond.setdefault(str(s), str(c))
    g = t.column("genes").combine_chunks(); e = t.column("expressions").combine_chunks()
    tok = np.asarray(g.values, dtype=np.int64); val = np.asarray(e.values, dtype=np.float64)
    off = np.asarray(g.offsets, dtype=np.int64); ncell = len(off) - 1
    cell_id = np.repeat(np.arange(ncell), np.diff(off))
    real = tok >= 3; tok, val, cid = tok[real], val[real], cell_id[real]
    tot = np.bincount(cid, weights=val, minlength=ncell); tot[tot == 0] = 1.0
    hj = hvg_map[tok]; keep = hj >= 0; cH, hH = cid[keep], hj[keep]
    normv = np.log1p(val[keep] * (TS / tot[cH]))
    M = sp.csr_matrix((normv, (cH, hH)), shape=(ncell, NH))
    coords = (M @ comps.T - c0).astype(np.float32)
    samp = np.asarray(t.column("sample").to_pylist()).astype(str)
    line = np.asarray(t.column("cell_line").to_pylist()).astype(str)
    return coords, samp, line


t0 = time.time()
for i in range(N_TOTAL):
    outf = os.path.join(ALLDIR, f"shard_{i:03d}.npz")
    if os.path.exists(outf):
        continue
    f = f"expression_data/train-{i:05d}-of-00116.parquet"; p = None
    for a in range(5):
        try:
            p = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE); break
        except Exception as ex:
            print(f"  [warn] dl shard {i} try {a}: {str(ex)[:80]}", flush=True); time.sleep(5)
    if p is None:
        print(f"[skip] shard {i} undownloadable", flush=True); continue
    try:
        coords, samp, line = project_shard(p)
        np.savez(outf, coords=coords, sample=samp, line=line)
    except Exception as ex:
        print(f"  [warn] proc shard {i}: {str(ex)[:120]}", flush=True)
    finally:
        try: os.remove(p)
        except OSError: pass
    r = (i + 1) / max(time.time() - t0, 1)
    print(f"[eb-all] shard {i+1}/{N_TOTAL} done | {r*60:.1f} shards/min", flush=True)

# --- concatenate ---
files = sorted(f for f in os.listdir(ALLDIR) if f.startswith("shard_"))
print(f"[eb-all] concatenating {len(files)} shard files...", flush=True)
Cs, Ss, Ls = [], [], []
for fn in files:
    z = np.load(os.path.join(ALLDIR, fn), allow_pickle=True)
    Cs.append(z["coords"]); Ss.append(z["sample"]); Ls.append(z["line"])
coords = np.concatenate(Cs); sample = np.concatenate(Ss); line = np.concatenate(Ls)
np.savez(os.path.join(OUT, "all_cells.npz"), coords=coords, sample=sample, line=line)
if s2cond:
    json.dump(s2cond, open(os.path.join(OUT, "sample2cond.json"), "w"))
print(f"[eb-all] DONE: {len(coords)} cells, {len(np.unique(line))} lines -> {OUT}/all_cells.npz", flush=True)
