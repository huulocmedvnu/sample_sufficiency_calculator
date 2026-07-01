#!/usr/bin/env python
"""Pass 2 of the Tahoe-100M recompute: stream ALL 3388 expression shards and accumulate the
per-(drug x cell_line) pseudobulk centroid in the HVG log-CP10K space defined by Pass 1.

For each cell (theislab recipe, applied streaming):
  total = sum of raw counts over real genes (token_id >= 3)
  factor = TARGET_SUM / total                      # normalize_total(1e4)
  value_g = log1p(count_g * factor)                # log1p, restricted to the 2000 HVG genes
and we accumulate sum(value) and cell-count per (drug, cell_line_id).

Resumable: checkpoints every CKPT_EVERY shards to $OUT/pass2_ckpt.npz (done-shard list + running
sums/counts). Re-running resumes from the last checkpoint. Peak disk stays small (one shard at a
time, deleted after use).

Output -> $OUT/pseudobulk.npz : cond_keys ("drug|line"), sums (n_cond, 2000), counts (n_cond,).
"""
import os, glob, time, threading, queue, numpy as np
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

N_WORKERS = int(os.environ.get("N_WORKERS", "8"))     # concurrent shard downloads (bandwidth ~caps here)

REPO = "tahoebio/Tahoe-100M"
CACHE = "/mnt/hdd2/loc-tran/tahoe_work/dl2"
OUT = os.environ.get("OUT", "/mnt/hdd2/loc-tran/tahoe_work/out")
N_TOTAL = 3388
CKPT_EVERY = int(os.environ.get("CKPT_EVERY", "50"))
CKPT = os.path.join(OUT, "pass2_ckpt.npz")
os.makedirs(CACHE, exist_ok=True)

b = np.load(os.path.join(OUT, "basis.npz"))
HVG = b["hvg_token_ids"].astype(np.int64)                 # (2000,) token ids
GENE_WIDTH = int(b["gene_width"]); TARGET_SUM = float(b["target_sum"])
NH = HVG.size
hvg_map = np.full(GENE_WIDTH, -1, dtype=np.int64)         # token_id -> hvg index (0..1999) or -1
hvg_map[HVG] = np.arange(NH)

# running accumulators as dicts keyed by "drug|line"
sums = {}; counts = {}
done = set()
if os.path.exists(CKPT):
    z = np.load(CKPT, allow_pickle=True)
    for k, s, c in zip(z["cond_keys"], z["sums"], z["counts"]):
        sums[str(k)] = s.astype(np.float64); counts[str(k)] = float(c)
    done = set(int(i) for i in z["done"])
    print(f"[pass2] resumed: {len(done)} shards done, {len(sums)} conditions so far")


def process_shard(path):
    # key on (sample x cell_line): `sample` uniquely identifies drug x dose x plate, so this
    # resolves the true (drug x dose x plate x line) unit (pooled to drug x dose x line in Pass 3).
    t = pq.read_table(path, columns=["genes", "expressions", "sample", "cell_line_id"])
    g = t.column("genes"); e = t.column("expressions")
    # flatten list-arrays: values + offsets
    gv = g.combine_chunks(); ev = e.combine_chunks()
    tokens = np.asarray(gv.values, dtype=np.int64)
    vals = np.asarray(ev.values, dtype=np.float64)
    off = np.asarray(gv.offsets, dtype=np.int64)              # (ncell+1,)
    ncell = len(off) - 1
    per_cell = np.diff(off)
    cell_id = np.repeat(np.arange(ncell), per_cell)
    real = tokens >= 3
    tok = tokens[real]; val = vals[real]; cid = cell_id[real]
    tot = np.bincount(cid, weights=val, minlength=ncell)      # per-cell total real counts
    tot[tot == 0] = 1.0
    factor = TARGET_SUM / tot
    hidx = hvg_map[tok]                                       # -1 if not HVG
    keep = hidx >= 0
    cidH = cid[keep]; hjdx = hidx[keep]
    normv = np.log1p(val[keep] * factor[cidH])
    # per-cell condition labels: sample|cell_line  (sample = drug x dose x plate)
    samp = np.asarray(t.column("sample").to_pylist()); line = np.asarray(t.column("cell_line_id").to_pylist())
    cond = np.char.add(np.char.add(samp.astype(str), "|"), line.astype(str))
    uconds, cinv = np.unique(cond, return_inverse=True)       # cinv: per-cell local cond id
    k = len(uconds)
    cond_of_entry = cinv[cidH]
    flat = cond_of_entry * NH + hjdx
    local_sum = np.bincount(flat, weights=normv, minlength=k * NH).reshape(k, NH)
    local_cnt = np.bincount(cinv, minlength=k).astype(np.float64)
    for j, ck in enumerate(uconds):
        ck = str(ck)
        if ck in sums:
            sums[ck] += local_sum[j]; counts[ck] += local_cnt[j]
        else:
            sums[ck] = local_sum[j].copy(); counts[ck] = float(local_cnt[j])


def save_ckpt(done_set, final=False):
    keys = np.array(sorted(sums.keys()))
    S = np.stack([sums[k] for k in keys]).astype(np.float32)
    C = np.array([counts[k] for k in keys], dtype=np.float64)
    tmp = CKPT + ".tmp.npz"
    np.savez(tmp, cond_keys=keys, sums=S, counts=C, done=np.array(sorted(done_set)))
    os.replace(tmp, CKPT)
    if final:
        np.savez(os.path.join(OUT, "pseudobulk.npz"), cond_keys=keys, sums=S, counts=C)


def main():
    """Prefetch shards concurrently (N_WORKERS threads) while the main thread processes them in
    completion order. Accumulators are touched only by the main thread. Bounded queue keeps at most
    ~24 shards on disk at once."""
    todo = [i for i in range(N_TOTAL) if i not in done]
    print(f"[pass2] {len(todo)} shards to stream ({N_WORKERS} download workers)", flush=True)
    idx_lock = threading.Lock(); idx_iter = iter(todo)
    ready = queue.Queue(maxsize=16)                       # (i, path|None)

    def worker():
        while True:
            with idx_lock:
                try:
                    i = next(idx_iter)
                except StopIteration:
                    return
            f = f"data/train-{i:05d}-of-03388.parquet"
            path = None
            for attempt in range(5):
                try:
                    path = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE); break
                except Exception as ex:
                    print(f"  [warn] dl shard {i} attempt {attempt}: {ex}", flush=True); time.sleep(5)
            ready.put((i, path))

    workers = [threading.Thread(target=worker, daemon=True) for _ in range(N_WORKERS)]
    for w in workers:
        w.start()

    t0 = time.time(); processed = 0; total = len(todo)
    for _ in range(total):
        i, p = ready.get()
        if p is None:
            print(f"  [ERROR] shard {i} failed after retries; skipping", flush=True)
        else:
            try:
                process_shard(p)
            except Exception as ex:
                print(f"  [warn] process shard {i}: {ex}", flush=True)
            finally:
                try: os.remove(p)
                except OSError: pass
        done.add(i); processed += 1
        if processed % CKPT_EVERY == 0 or processed == total:
            save_ckpt(done)
            el = time.time() - t0; rate = processed / max(el, 1)
            eta = (total - processed) / max(rate, 1e-9) / 3600
            print(f"[pass2] {len(done)}/{N_TOTAL} shards | {len(sums)} conds | "
                  f"{rate*3600:.0f} shards/hr | ETA {eta:.1f} h | ckpt saved", flush=True)
    save_ckpt(done, final=True)
    print(f"[pass2] DONE: {len(sums)} conditions -> {OUT}/pseudobulk.npz", flush=True)


if __name__ == "__main__":
    main()
