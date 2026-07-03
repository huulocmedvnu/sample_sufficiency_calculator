"""
Pass 1 (pilot): fit the HVG(2000)+PCA(50) embedding on the X-Atlas/Orion head-block and estimate
both per-cell variances the quota can use:
  * sigma^2 MARGINAL   -- mean per-PC variance over all block cells (mixes perturbations)   [conservative bound]
  * sigma^2 WITHIN     -- pooled residual variance after removing each gene_target's own mean [theory-preferred]
    (this is the noise that actually enters a knockdown centroid; Non-Targeting controls included)
Also stores the per-PC WITHIN-condition variances ell_k (diagonal of Sigma in the PCA basis) for the
anisotropic trace tr(P Sigma P) = sum_k (1 - u_k^2) ell_k.

Same fixed recipe as the Tahoe/EmeraldBay recomputes: normalize_total(1e4) -> log1p ->
highly_variable_genes(2000, seurat) -> PCA(50). sigma^2 is a platform/pipeline-specific plug-in and
MUST be re-estimated here (that is the whole point of running it on a new modality).

Outputs -> $OUT/<LINE>/basis.npz : components, pca_mean, hvg_ids, sigma2_marginal, sigma2_within,
           ell_within (50,), evr, coords (n,50) [per-cell PCA], cell_integer_ids (n,), n_cells_fit.
"""
import os, argparse, numpy as np, scipy.sparse as sp
import pandas as pd, scanpy as sc, anndata as ad

TARGET_SUM = 1e4
N_COMPS = 50
MIN_CELLS_WITHIN = 20   # gene_targets with >= this many block cells contribute to the within-condition pool

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HCT116")
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/orion_work/pilot")
    args = ap.parse_args()
    d = os.path.join(args.dir, args.line)

    meta = pd.read_parquet(os.path.join(d, "cells_meta.parquet"))
    genes = pd.read_parquet(os.path.join(d, "genes.parquet"))
    z = np.load(os.path.join(d, "expr_coo.npz"))
    cell_idx, gene_idx, val = z["cell_idx"], z["gene_idx"], z["val"]
    cell_ids = z["cell_integer_ids"]                       # block-local cell order -> cell_integer_id
    n_cells = cell_ids.size
    n_genes = int(genes["gene_integer_id"].max()) + 1
    print(f"[pass1] {args.line}: {n_cells:,} cells x {n_genes:,} genes, {val.size:,} nnz", flush=True)

    # align metadata rows to the block-local cell order (by cell_integer_id)
    meta = meta.set_index("cell_integer_id").loc[cell_ids].reset_index()
    assert (meta["cell_integer_id"].values == cell_ids).all()

    X = sp.csr_matrix((val.astype(np.float32), (cell_idx, gene_idx)),
                      shape=(n_cells, n_genes))
    A = ad.AnnData(X=X)
    A.obs["gene_target"] = meta["gene_target"].values
    A.obs["sample"] = meta["sample"].values

    sc.pp.normalize_total(A, target_sum=TARGET_SUM)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat")
    hvg_mask = A.var["highly_variable"].values
    hvg_ids = np.where(hvg_mask)[0].astype(np.int64)
    A = A[:, hvg_mask].copy()
    print(f"[pass1] HVG selected: {A.n_vars}", flush=True)
    sc.pp.pca(A, n_comps=N_COMPS, zero_center=True, svd_solver="arpack")
    comps = A.varm["PCs"].T.astype(np.float64)
    pca_mean = np.asarray(A.X.mean(axis=0)).ravel().astype(np.float64)
    coords = A.obsm["X_pca"].astype(np.float64)
    evr = A.uns["pca"]["variance_ratio"].astype(np.float64)

    # marginal sigma^2
    per_pc_marg = coords.var(axis=0)
    sigma2_marg = float(per_pc_marg.mean())

    # within-condition: pooled residual after removing each gene_target's own mean
    gt = meta["gene_target"].values
    order = np.argsort(gt, kind="stable")
    gts, starts = np.unique(gt[order], return_index=True)
    ends = np.r_[starts[1:], len(gt)]
    ss = np.zeros(N_COMPS); dof = 0; used = 0
    for g0, s, e in zip(gts, starts, ends):
        idx = order[s:e]
        if idx.size < MIN_CELLS_WITHIN:
            continue
        C = coords[idx]
        ss += ((C - C.mean(0)) ** 2).sum(0)
        dof += idx.size - 1
        used += 1
    ell_within = ss / dof                          # per-PC within-condition variance
    sigma2_within = float(ell_within.mean())
    print(f"[pass1] sigma^2 marginal = {sigma2_marg:.4f}  |  within-condition = {sigma2_within:.4f} "
          f"(pooled over {used} gene_targets, dof {dof:,})", flush=True)
    print(f"[pass1] top-5 evr = {evr[:5].round(4).tolist()}", flush=True)

    np.savez(os.path.join(d, "basis.npz"),
             components=comps, pca_mean=pca_mean, hvg_ids=hvg_ids,
             sigma2_marginal=sigma2_marg, sigma2_within=sigma2_within,
             ell_within=ell_within, per_pc_marginal=per_pc_marg, evr=evr,
             coords=coords, cell_integer_ids=cell_ids, n_cells_fit=n_cells,
             target_sum=TARGET_SUM, n_comps=N_COMPS)
    meta.to_parquet(os.path.join(d, "cells_meta_aligned.parquet"))
    print(f"[pass1] wrote {d}/basis.npz + cells_meta_aligned.parquet", flush=True)

if __name__ == "__main__":
    main()
