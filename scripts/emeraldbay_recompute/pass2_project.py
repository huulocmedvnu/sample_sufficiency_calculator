#!/usr/bin/env python
"""Pass 2 of the EmeraldBay recompute: stream all 116 shards, project every cell into the Pass-1
PCA(50) space, and accumulate per-(sample × cell_line) centroid moments (sum, sum-of-squares, count)
for magnitudes and the WITHIN-condition sigma^2. Additionally RETAIN per-cell PCA coordinates for the
five cell lines shared with Tahoe-100M (HS-578T, AN3-CA, HEC-1-A, BT-474, C-33 A) so Pass 3 can run the
held-out angular-error subsampling validation on real independent cells.

Outputs $OUT_EB/pseudobulk.npz (cond_keys "sample|line", sums (n,50), sumsq (n,50), counts (n,))
        $OUT_EB/shared_cells.npz (coords (m,50), sample (m,), line (m,))  -- 5 shared lines only.
"""
import os, time, threading, queue, numpy as np, scipy.sparse as sp
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download

N_WORKERS = int(os.environ.get("N_WORKERS", "8"))

REPO = "tahoebio/EmeraldBay"
CACHE = "/mnt/hdd2/loc-tran/eb_work/dl2"
OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
N_TOTAL = 116
SHARED = {"CVCL_0332", "CVCL_0028", "CVCL_0293", "CVCL_0179", "CVCL_1094"}  # HS-578T,AN3-CA,HEC-1-A,BT-474,C-33 A
os.makedirs(CACHE, exist_ok=True)

b = np.load(os.path.join(OUT, "basis.npz"))
HVG = b["hvg_token_ids"].astype(np.int64); GW = int(b["gene_width"]); TS = float(b["target_sum"])
NH = HVG.size; NC = int(b["n_comps"])
comps = b["components"].astype(np.float64)            # (50, 2000)
pca_mean = b["pca_mean"].astype(np.float64)           # (2000,)
c0 = pca_mean @ comps.T                                # (50,) projection offset
hvg_map = np.full(GW, -1, dtype=np.int64); hvg_map[HVG] = np.arange(NH)

sums, sumsq, counts = {}, {}, {}
shared_coords, shared_samp, shared_line = [], [], []
s2cond = {}                                            # sample -> drugname_drugconc (drug x dose)


def project_shard(path):
    t = pq.read_table(path, columns=["genes", "expressions", "sample", "cell_line",
                                     "drug", "drugname_drugconc"])
    for s, c in zip(t.column("sample").to_pylist(), t.column("drugname_drugconc").to_pylist()):
        s2cond.setdefault(str(s), str(c))
    g = t.column("genes").combine_chunks(); e = t.column("expressions").combine_chunks()
    tok = np.asarray(g.values, dtype=np.int64); val = np.asarray(e.values, dtype=np.float64)
    off = np.asarray(g.offsets, dtype=np.int64); ncell = len(off) - 1
    cell_id = np.repeat(np.arange(ncell), np.diff(off))
    real = tok >= 3
    tok, val, cid = tok[real], val[real], cell_id[real]
    tot = np.bincount(cid, weights=val, minlength=ncell); tot[tot == 0] = 1.0
    hj = hvg_map[tok]; keep = hj >= 0
    cH, hH = cid[keep], hj[keep]
    normv = np.log1p(val[keep] * (TS / tot[cH]))
    M = sp.csr_matrix((normv, (cH, hH)), shape=(ncell, NH))     # log-norm HVG (ncell x 2000)
    coords = M @ comps.T - c0                                    # (ncell x 50)
    samp = np.asarray(t.column("sample").to_pylist()); line = np.asarray(t.column("cell_line").to_pylist())
    cond = np.char.add(np.char.add(samp.astype(str), "|"), line.astype(str))
    for u in np.unique(cond):
        m = cond == u; C = coords[m]
        if u in sums:
            sums[u] += C.sum(0); sumsq[u] += (C * C).sum(0); counts[u] += C.shape[0]
        else:
            sums[u] = C.sum(0); sumsq[u] = (C * C).sum(0); counts[u] = C.shape[0]
    sh = np.isin(line, list(SHARED))
    if sh.any():
        shared_coords.append(coords[sh]); shared_samp.append(samp[sh]); shared_line.append(line[sh])


def main():
    t0 = time.time()
    lock = threading.Lock(); it = iter(range(N_TOTAL)); ready = queue.Queue(maxsize=12)

    def worker():
        while True:
            with lock:
                try:
                    i = next(it)
                except StopIteration:
                    return
            f = f"expression_data/train-{i:05d}-of-00116.parquet"; path = None
            for a in range(5):
                try:
                    path = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE); break
                except Exception as ex:
                    print(f"  [warn] dl shard {i} try {a}: {ex}", flush=True); time.sleep(5)
            ready.put((i, path))

    ws = [threading.Thread(target=worker, daemon=True) for _ in range(N_WORKERS)]
    for w in ws:
        w.start()
    for k in range(N_TOTAL):
        i, p = ready.get()
        if p is not None:
            try:
                project_shard(p)
            except Exception as ex:
                print(f"  [warn] proc shard {i}: {ex}", flush=True)
            finally:
                try: os.remove(p)
                except OSError: pass
        if (k + 1) % 10 == 0 or k == N_TOTAL - 1:
            r = (k + 1) / max(time.time() - t0, 1)
            print(f"[pass2-eb] {k+1}/{N_TOTAL} | {len(sums)} conds | {r*60:.1f}/min", flush=True)
    keys = np.array(sorted(sums.keys()))
    np.savez(os.path.join(OUT, "pseudobulk.npz"), cond_keys=keys,
             sums=np.stack([sums[k] for k in keys]), sumsq=np.stack([sumsq[k] for k in keys]),
             counts=np.array([counts[k] for k in keys], dtype=np.float64))
    np.savez(os.path.join(OUT, "shared_cells.npz"),
             coords=np.concatenate(shared_coords), sample=np.concatenate(shared_samp),
             line=np.concatenate(shared_line))
    import json
    json.dump(s2cond, open(os.path.join(OUT, "sample2cond.json"), "w"))
    print(f"[pass2-eb] DONE: {len(keys)} conditions; {sum(c.shape[0] for c in shared_coords)} shared-line cells retained")


if __name__ == "__main__":
    main()
