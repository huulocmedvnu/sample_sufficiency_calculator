"""
Pass 0 (pilot): download a representative HEAD-BLOCK of X-Atlas/Orion (SLAF/Lance) to local disk.

Cells in the SLAF cells table are ordered by batch; within each batch the perturbations
(and the Non-Targeting controls) are well-mixed and representative of the global composition.
So the first K cells form a valid representative subsample with co-located controls -- enough to
fit the HVG/PCA basis, estimate the within-condition per-cell variance sigma^2, and compute the
n* quota for whatever perturbations are well-sampled in the block. The full 18k-gene spectrum
needs the whole expression table (~17.4B nnz, ~6h stream); this pilot is the bounded proof.

Expression is COO (cell_integer_id uint32, gene_integer_id uint16, value uint16), sorted by cell.
We stream contiguous row chunks and stop once cell_integer_id crosses the head threshold.

Outputs (in $OUT, default /mnt/hdd2/loc-tran/orion_work/pilot/<LINE>):
  cells_meta.parquet   -- per-cell metadata for the block (incl. gene_target, sample, cell_integer_id)
  genes.parquet        -- full gene table (gene_integer_id -> gene_id, mean_counts, etc.)
  expr_coo.npz         -- COO arrays (cell_idx, gene_idx, val) for the block, cell_idx re-based to 0..n-1
"""
import os, sys, time, argparse
import numpy as np, pandas as pd
import lance

HF = "hf://datasets/slaf-project/X-Atlas-Orion"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HCT116", choices=["HCT116", "HEK293T"])
    ap.add_argument("--max-cells", type=int, default=250_000, help="head block size (cells)")
    ap.add_argument("--chunk", type=int, default=8_000_000, help="expression rows per read")
    ap.add_argument("--out", default=os.environ.get("OUT", "/mnt/hdd2/loc-tran/orion_work/pilot"))
    args = ap.parse_args()

    outdir = os.path.join(args.out, args.line)
    os.makedirs(outdir, exist_ok=True)
    K = args.max_cells

    print(f"[pass0] line={args.line}  head={K:,} cells  -> {outdir}", flush=True)

    # --- cells + genes metadata (small) ---
    cds = lance.dataset(f"{HF}/data/{args.line}/cells.lance")
    cols = ["cell_integer_id", "gene_target", "guide_target", "sample",
            "total_counts", "pct_counts_mt", "n_genes_by_counts", "pass_guide_filter"]
    cm = cds.to_table(columns=cols).to_pandas().sort_values("cell_integer_id").reset_index(drop=True)
    cm = cm[cm["cell_integer_id"] < K].reset_index(drop=True)
    cm.to_parquet(os.path.join(outdir, "cells_meta.parquet"))
    print(f"[pass0] cells in block: {len(cm):,}  (gene_targets={cm.gene_target.nunique()}, "
          f"NTC={int((cm.gene_target=='Non-Targeting').sum())}, batches={cm['sample'].nunique()})", flush=True)

    gds = lance.dataset(f"{HF}/data/{args.line}/genes.lance")
    gm = gds.to_table().to_pandas().sort_values("gene_integer_id").reset_index(drop=True)
    gm.to_parquet(os.path.join(outdir, "genes.parquet"))
    print(f"[pass0] genes: {len(gm):,}", flush=True)

    # --- expression block: contiguous stream, stop at threshold ---
    ex = lance.dataset(f"{HF}/data/{args.line}/expression.lance")
    total = ex.count_rows()
    ci_parts, gi_parts, v_parts = [], [], []
    got = 0
    offset = 0
    t0 = time.time()
    while offset < total:
        tab = ex.scanner(columns=["cell_integer_id", "gene_integer_id", "value"],
                         limit=args.chunk, offset=offset).to_table()
        n = tab.num_rows
        if n == 0:
            break
        ci = tab.column("cell_integer_id").to_numpy()
        keep = ci < K
        if keep.any():
            gi = tab.column("gene_integer_id").to_numpy()
            v = tab.column("value").to_numpy()
            ci_parts.append(ci[keep].astype(np.uint32))
            gi_parts.append(gi[keep].astype(np.uint16))
            v_parts.append(v[keep].astype(np.uint16))
            got += int(keep.sum())
        offset += n
        last = int(ci[-1])
        el = time.time() - t0
        print(f"[pass0]   offset={offset:>14,}  last_cell={last:>9,}  kept_nnz={got:>13,}  {el:6.1f}s", flush=True)
        if last >= K:
            break

    ci = np.concatenate(ci_parts); gi = np.concatenate(gi_parts); v = np.concatenate(v_parts)
    # re-base cell ids to 0..n-1 in block order
    uniq = np.unique(ci)
    remap = np.full(uniq.max() + 1, -1, dtype=np.int64)
    remap[uniq] = np.arange(uniq.size)
    cell_idx = remap[ci].astype(np.int32)
    np.savez_compressed(os.path.join(outdir, "expr_coo.npz"),
                        cell_idx=cell_idx, gene_idx=gi.astype(np.int32), val=v.astype(np.float32),
                        cell_integer_ids=uniq.astype(np.int64))
    print(f"[pass0] DONE  nnz={got:,}  cells_with_expr={uniq.size:,}  "
          f"({time.time()-t0:.0f}s, {got/1e6/(time.time()-t0):.2f}M nnz/s)", flush=True)

if __name__ == "__main__":
    main()
