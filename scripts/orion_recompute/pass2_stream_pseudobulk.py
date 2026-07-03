"""
Pass 2 (FULL SCALE, resumable): stream the entire X-Atlas/Orion expression table for one cell line and
accumulate, per gene_target, the pseudobulk sufficient statistics in PCA(50) space:
    count, sum(coords), sum(coords^2).
From these fall out everything the quota needs: each perturbation centroid, the Non-Targeting centroid,
the pooled WITHIN-condition per-cell variance sigma^2 (and per-PC ell_k for the anisotropic trace), and
the global marginal sigma^2.

Projection (validated to 1e-4 vs scanpy): coord = log1p(count/total_counts*1e4)[HVG] projected through
the pre-fit basis:  coords = X_hvg_lognorm @ components.T - (pca_mean @ components.T).
normalize_total denominator = cells.total_counts (verified == raw row sum exactly).

Resumable: expression is sorted by cell; each chunk we consume only COMPLETE cells (rows with
cell < last_cell_in_chunk) and advance the offset by exactly those rows, so a checkpoint is just the
row offset + accumulators. Kill/restart continues from the checkpoint.

Outputs -> $OUT/<LINE>/pseudobulk.npz : gt_names, cnt, csum(G,50), csq(G,50), n_cells_done, offset_done.
"""
import os, argparse, time, numpy as np, scipy.sparse as sp, pandas as pd
import lance

HF = "hf://datasets/slaf-project/X-Atlas-Orion"
NTC = "Non-Targeting"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", required=True)
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/orion_work/full")
    ap.add_argument("--basis", required=True, help="path to basis.npz (components, pca_mean, hvg_ids)")
    ap.add_argument("--chunk", type=int, default=16_000_000)
    ap.add_argument("--ckpt-every", type=int, default=40)
    args = ap.parse_args()
    d = os.path.join(args.dir, args.line); os.makedirs(d, exist_ok=True)
    ckpt = os.path.join(d, "pseudobulk.npz")

    b = np.load(args.basis)
    comps = b["components"].astype(np.float64)              # (50,2000)
    hvg = b["hvg_ids"]; D = comps.shape[0]
    mean_proj = (b["pca_mean"].astype(np.float64) @ comps.T)  # (50,)
    n_genes_full = 38606
    hvg_col = np.full(n_genes_full, -1, np.int64); hvg_col[hvg] = np.arange(hvg.size)

    # cells table: total_counts + gene_target, indexed by cell_integer_id (0..n-1 contiguous)
    cds = lance.dataset(f"{HF}/data/{args.line}/cells.lance")
    cm = cds.to_table(columns=["cell_integer_id", "total_counts", "gene_target"]).to_pandas()
    cm = cm.sort_values("cell_integer_id").reset_index(drop=True)
    assert (cm["cell_integer_id"].values == np.arange(len(cm))).all(), "cell ids not 0..n-1 contiguous"
    tot = cm["total_counts"].values.astype(np.float64)
    gt_names, gt_code = np.unique(cm["gene_target"].values, return_inverse=True)
    gt_code = gt_code.astype(np.int32); G = gt_names.size
    print(f"[pass2:{args.line}] {len(cm):,} cells, {G:,} gene_targets (NTC idx "
          f"{int(np.where(gt_names==NTC)[0][0])})", flush=True)

    ex = lance.dataset(f"{HF}/data/{args.line}/expression.lance")
    total_rows = ex.count_rows()

    if os.path.exists(ckpt):
        z = np.load(ckpt, allow_pickle=True)
        csum = z["csum"].astype(np.float64); csq = z["csq"].astype(np.float64)
        cnt = z["cnt"].astype(np.float64); offset = int(z["offset_done"])
        print(f"[pass2:{args.line}] RESUME from offset {offset:,}/{total_rows:,} "
              f"({100*offset/total_rows:.1f}%)", flush=True)
    else:
        csum = np.zeros((G, D)); csq = np.zeros((G, D)); cnt = np.zeros(G); offset = 0

    def save(off):
        np.savez(ckpt, gt_names=gt_names, cnt=cnt, csum=csum, csq=csq,
                 offset_done=off, n_cells_done=int(cnt.sum()), n_genes=G, n_comps=D)

    def read_chunk(off):
        last_err = None
        for attempt in range(6):
            try:
                return ex.scanner(columns=["cell_integer_id", "gene_integer_id", "value"],
                                  limit=args.chunk, offset=off).to_table()
            except Exception as e:                       # transient network / object-store blip
                last_err = e
                print(f"[pass2:{args.line}] read retry {attempt+1} at offset {off:,}: "
                      f"{type(e).__name__}", flush=True)
                time.sleep(min(60, 5 * 2**attempt))
        raise last_err

    t0 = time.time(); ch = 0
    while offset < total_rows:
        tab = read_chunk(offset)
        n = tab.num_rows
        if n == 0:
            break
        cid = tab.column("cell_integer_id").to_numpy()
        last = cid[-1]
        end_of_table = (n < args.chunk) or (offset + n >= total_rows)
        if end_of_table:
            complete = np.ones(n, dtype=bool); consumed = n
        else:
            complete = cid < last                 # drop the trailing partial cell
            consumed = int(complete.sum())
        gid = tab.column("gene_integer_id").to_numpy()
        val = tab.column("value").to_numpy().astype(np.float64)
        cid = cid[complete]; gid = gid[complete]; val = val[complete]
        keep = hvg_col[gid] >= 0
        cid = cid[keep]; col = hvg_col[gid[keep]]; val = val[keep]
        # local dense cell indexing
        uc, inv = np.unique(cid, return_inverse=True)
        ln = np.log1p(val / tot[cid] * 1e4)
        S = sp.csr_matrix((ln, (inv, col)), shape=(uc.size, hvg.size))
        coords = S @ comps.T - mean_proj          # (uc,50)
        gc = gt_code[uc]
        np.add.at(csum, gc, coords)
        np.add.at(csq, gc, coords * coords)
        np.add.at(cnt, gc, 1.0)
        offset += consumed
        ch += 1
        if ch % args.ckpt_every == 0 or end_of_table:
            save(offset)
            el = time.time() - t0
            print(f"[pass2:{args.line}] offset={offset:>14,}/{total_rows:,} "
                  f"({100*offset/total_rows:5.1f}%)  cells={int(cnt.sum()):>9,}  "
                  f"{el/60:6.1f}min  ~{offset/el/1e6:.2f}M rows/s", flush=True)
        if end_of_table:
            break
    save(offset)
    print(f"[pass2:{args.line}] DONE  {int(cnt.sum()):,} cells  offset {offset:,}  "
          f"{(time.time()-t0)/60:.1f}min", flush=True)

if __name__ == "__main__":
    main()
