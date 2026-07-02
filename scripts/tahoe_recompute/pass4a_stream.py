#!/usr/bin/env python
"""Full Tahoe-100M re-stream + projection (stage A of the direct-curve falsification).

Streams ALL 3,388 expression shards from HuggingFace (one at a time, deleted after use), projects
every cell into the SAME calibrated PCA(50) basis (Pass 1), and writes per-cell coordinates to a disk
MEMMAP (coords.f32 ~20 GB, cond_id.i32) plus a per-line baseline accumulator. Condition unit =
(sample x cell_line_id) = one drug treatment at one dose in one line.

RESUMABLE: a checkpoint (cursor, done-shard bitmap, condition key table, per-line sums/counts) is saved
every CKPT_EVERY shards; re-running resumes without re-downloading finished shards. Stage B
(pass4b_curves.py) then runs the subsample-and-measure curves per condition from the memmap.

Env: N_WORKERS(8), CKPT_EVERY(25), MAXN(101_000_000).
"""
import os, json, time, threading, queue, numpy as np, scipy.sparse as sp
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

REPO = "tahoebio/Tahoe-100M"
CACHE = "/mnt/hdd2/loc-tran/tahoe_work/dl2"
BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"; os.makedirs(BASE, exist_ok=True); os.makedirs(CACHE, exist_ok=True)
_dose = "/mnt/hdd2/loc-tran/tahoe_work/out_dose"
BASIS = _dose if os.path.exists(os.path.join(_dose, "basis.npz")) else "/mnt/hdd2/loc-tran/tahoe_work/out"
N_TOTAL = 3388
N_WORKERS = int(os.environ.get("N_WORKERS", "8"))
CKPT_EVERY = int(os.environ.get("CKPT_EVERY", "25"))
MAXN = int(os.environ.get("MAXN", "101000000"))
CKPT = os.path.join(BASE, "ckpt.npz")

b = np.load(os.path.join(BASIS, "basis.npz"))
HVG = b["hvg_token_ids"].astype(np.int64); GW = int(b["gene_width"]); TS = float(b["target_sum"])
NH = HVG.size; NC = int(b["n_comps"]); comps = b["components"].astype(np.float64)
pca_mean = b["pca_mean"].astype(np.float64); c0 = pca_mean @ comps.T
hvg_map = np.full(GW, -1, dtype=np.int64); hvg_map[HVG] = np.arange(NH)

coords_mm = np.memmap(os.path.join(BASE, "coords.f32"), dtype=np.float32, mode="r+" if os.path.exists(os.path.join(BASE, "coords.f32")) else "w+", shape=(MAXN, NC))
condid_mm = np.memmap(os.path.join(BASE, "condid.i32"), dtype=np.int32, mode="r+" if os.path.exists(os.path.join(BASE, "condid.i32")) else "w+", shape=(MAXN,))

# --- restore checkpoint ---
if os.path.exists(CKPT):
    z = np.load(CKPT, allow_pickle=True)
    cursor = int(z["cursor"]); done = set(z["done"].tolist())
    cond_keys = list(z["cond_keys"]); cond_id = {k: i for i, k in enumerate(cond_keys)}
    cond_drug = dict(z["cond_drug"].item()) if "cond_drug" in z else {}
    line_sum = dict(z["line_sum"].item()); line_cnt = dict(z["line_cnt"].item())
    print(f"[resume] cursor={cursor:,} done_shards={len(done)} conds={len(cond_keys)}", flush=True)
else:
    cursor = 0; done = set(); cond_keys = []; cond_id = {}; cond_drug = {}
    line_sum = {}; line_cnt = {}


def save_ckpt():
    np.savez(CKPT, cursor=cursor, done=np.array(sorted(done)), cond_keys=np.array(cond_keys),
             cond_drug=np.array(cond_drug, dtype=object), line_sum=np.array(line_sum, dtype=object),
             line_cnt=np.array(line_cnt, dtype=object))
    coords_mm.flush(); condid_mm.flush()


def project_shard(path):
    t = pq.read_table(path, columns=["genes", "expressions", "sample", "cell_line_id", "drug"])
    g = t.column("genes").combine_chunks(); e = t.column("expressions").combine_chunks()
    tok = np.asarray(g.values, dtype=np.int64); val = np.asarray(e.values, dtype=np.float64)
    off = np.asarray(g.offsets, dtype=np.int64); ncell = len(off) - 1
    cid = np.repeat(np.arange(ncell), np.diff(off))
    real = tok >= 3; tok, val, cc = tok[real], val[real], cid[real]
    tot = np.bincount(cc, weights=val, minlength=ncell); tot[tot == 0] = 1.0
    hj = hvg_map[tok]; keep = hj >= 0
    M = sp.csr_matrix((np.log1p(val[keep] * (TS / tot[cc[keep]])), (cc[keep], hj[keep])), shape=(ncell, NH))
    coords = (M @ comps.T - c0).astype(np.float32)
    samp = np.asarray(t.column("sample").to_pylist()).astype(str)
    line = np.asarray(t.column("cell_line_id").to_pylist()).astype(str)
    drug = np.asarray(t.column("drug").to_pylist()).astype(str)
    return coords, samp, line, drug


# --- threaded downloader ---
lock = threading.Lock(); todo = iter([i for i in range(N_TOTAL) if i not in done])
ready = queue.Queue(maxsize=16); n_to_do = N_TOTAL - len(done)


def worker():
    while True:
        with lock:
            try: i = next(todo)
            except StopIteration: return
        f = f"data/train-{i:05d}-of-03388.parquet"; p = None
        for a in range(5):
            try:
                p = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE); break
            except Exception as ex:
                time.sleep(5)
        ready.put((i, p))


for _ in range(N_WORKERS):
    threading.Thread(target=worker, daemon=True).start()

t0 = time.time(); processed = 0
for k in range(n_to_do):
    i, p = ready.get()
    if p is None:
        print(f"[warn] shard {i} undownloadable", flush=True); done.add(i); continue
    try:
        coords, samp, line, drug = project_shard(p)
        n = len(coords)
        if cursor + n > MAXN:
            print(f"[FATAL] MAXN exceeded at shard {i} (cursor {cursor:,}+{n:,}); raise MAXN", flush=True)
            save_ckpt(); raise SystemExit(1)
        # map conditions
        cond = np.char.add(np.char.add(samp, "|"), line)
        ids = np.empty(n, dtype=np.int32)
        for u in np.unique(cond):
            if u not in cond_id:
                cond_id[u] = len(cond_keys); cond_keys.append(u)
                cond_drug[u] = str(drug[cond == u][0])
            ids[cond == u] = cond_id[u]
        coords_mm[cursor:cursor + n] = coords
        condid_mm[cursor:cursor + n] = ids
        cursor += n
        for l in np.unique(line):
            ml = line == l; s = coords[ml].sum(0)
            line_sum[l] = line_sum.get(l, np.zeros(NC)) + s
            line_cnt[l] = line_cnt.get(l, 0) + int(ml.sum())
    except Exception as ex:
        print(f"[warn] proc shard {i}: {str(ex)[:120]}", flush=True)
    finally:
        try: os.remove(p)
        except OSError: pass
    done.add(i); processed += 1
    if processed % CKPT_EVERY == 0:
        save_ckpt()
        r = processed / max(time.time() - t0, 1)
        eta = (n_to_do - processed) / max(r, 1e-9) / 60
        print(f"[tahoe-stream] {len(done)}/{N_TOTAL} shards | cursor={cursor:,} | conds={len(cond_keys)} "
              f"| {r*60:.1f} shards/min | ETA {eta:.0f} min", flush=True)

save_ckpt()
json.dump({"n_cells": int(cursor), "n_conditions": len(cond_keys), "n_shards": len(done)},
          open(os.path.join(BASE, "stream_summary.json"), "w"), indent=1)
print(f"[tahoe-stream] DONE: {cursor:,} cells, {len(cond_keys)} conditions, {len(done)} shards", flush=True)
