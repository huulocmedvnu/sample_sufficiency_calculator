#!/usr/bin/env python
"""Pass 1 of the Tahoe-100M recompute: fit the HVG + PCA(50) embedding basis and sigma^2
on a multi-plate subsample, following the theislab vevo_100m recipe
(https://theislab.github.io/vevo_Tahoe_100m_analysis/vevo_100m_pca.html):

  filter (pass_filter=="full" is already applied in this release) -> normalize_total ->
  log1p -> highly_variable_genes(n_top_genes=2000) -> PCA.

Deviations (documented, for a streaming pseudobulk analysis rather than cell-level viz):
 * normalize_total target_sum = 1e4 (CP10K), a fixed target so the same transform can be applied
   shard-by-shard in Pass 2 without a global median pass.
 * PCA n_comps = 50 (the sample-sufficiency framework's embedding dimension; the tutorial uses 300
   for UMAP visualisation). sigma^2 is the mean over the 50 PCs of the per-cell coordinate variance.
 * HVG(2000) computed on a multi-plate subsample (the tutorial does per-plate HVG then consensus).

Outputs -> $OUT/basis.npz : hvg_token_ids (2000,), components (50,2000), pca_mean (2000,),
           sigma2 (float, mean per-PC per-cell variance), evr (50,), n_cells_fit.
"""
import os, glob, numpy as np, scipy.sparse as sp
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download
import scanpy as sc, anndata as ad

REPO = "tahoebio/Tahoe-100M"
CACHE = "/mnt/hdd2/loc-tran/tahoe_work/dl"
OUT = os.environ.get("OUT", "/mnt/hdd2/loc-tran/tahoe_work/out")
N_SHARDS_FIT = int(os.environ.get("N_SHARDS_FIT", "14"))   # spread across the 3388 (~14 plates)
GENE_WIDTH = 62713          # max token_id (62712) + 1; special tokens 0,1,2 dropped
TARGET_SUM = 1e4
N_COMPS = 50
os.makedirs(CACHE, exist_ok=True); os.makedirs(OUT, exist_ok=True)


def shard_to_csr(path):
    """Load one expression shard -> (CSR raw counts over GENE_WIDTH cols, drug, line arrays)."""
    t = pq.read_table(path, columns=["genes", "expressions", "drug", "cell_line_id"])
    genes = t.column("genes").to_pylist()
    expr = t.column("expressions").to_pylist()
    drug = np.array(t.column("drug").to_pylist())
    line = np.array(t.column("cell_line_id").to_pylist())
    indptr = np.zeros(len(genes) + 1, dtype=np.int64)
    idx_parts, val_parts = [], []
    for i, (g, e) in enumerate(zip(genes, expr)):
        g = np.asarray(g, dtype=np.int64); e = np.asarray(e, dtype=np.float32)
        keep = g >= 3                       # drop special tokens 0,1,2
        g = g[keep]; e = e[keep]
        idx_parts.append(g); val_parts.append(e)
        indptr[i + 1] = indptr[i] + g.size
    indices = np.concatenate(idx_parts); data = np.concatenate(val_parts)
    X = sp.csr_matrix((data, indices, indptr), shape=(len(genes), GENE_WIDTH))
    return X, drug, line


def main():
    # spread shard indices across the full range to span plates
    idxs = np.linspace(0, 3387, N_SHARDS_FIT).round().astype(int)
    idxs = sorted(set(int(i) for i in idxs))
    print(f"[pass1] fitting basis on {len(idxs)} shards: {idxs}")
    Xs, drugs, lines = [], [], []
    for i in idxs:
        f = f"data/train-{i:05d}-of-03388.parquet"
        p = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE)
        X, d, l = shard_to_csr(p)
        Xs.append(X); drugs.append(d); lines.append(l)
        os.remove(p)
        print(f"  shard {i}: {X.shape[0]} cells")
    X = sp.vstack(Xs).tocsr()
    A = ad.AnnData(X=X)
    A.obs["drug"] = np.concatenate(drugs); A.obs["cell_line_id"] = np.concatenate(lines)
    print(f"[pass1] subsample AnnData: {A.n_obs} cells x {A.n_vars} genes")

    # --- theislab recipe ---
    sc.pp.normalize_total(A, target_sum=TARGET_SUM)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat")
    hvg_mask = A.var["highly_variable"].values
    hvg_token_ids = np.where(hvg_mask)[0].astype(np.int64)          # column index == token_id
    A = A[:, hvg_mask].copy()
    print(f"[pass1] HVG selected: {A.n_vars} genes")
    sc.pp.pca(A, n_comps=N_COMPS, zero_center=True, svd_solver="arpack")
    comps = A.varm["PCs"].T.astype(np.float64)                      # (50, 2000)
    pca_mean = np.asarray(A.X.mean(axis=0)).ravel().astype(np.float64)  # zero-center mean (2000,)
    coords = A.obsm["X_pca"].astype(np.float64)                     # (n_cells, 50)
    per_pc_var = coords.var(axis=0)                                 # variance per PC across cells
    sigma2 = float(per_pc_var.mean())
    evr = A.uns["pca"]["variance_ratio"].astype(np.float64)
    print(f"[pass1] sigma^2 (mean per-PC per-cell var) = {sigma2:.4f}")
    print(f"[pass1] top-5 evr = {evr[:5].round(4).tolist()}")
    np.savez(os.path.join(OUT, "basis.npz"),
             hvg_token_ids=hvg_token_ids, components=comps, pca_mean=pca_mean,
             sigma2=sigma2, per_pc_var=per_pc_var, evr=evr, n_cells_fit=A.n_obs,
             target_sum=TARGET_SUM, n_comps=N_COMPS, gene_width=GENE_WIDTH)
    print(f"[pass1] wrote {OUT}/basis.npz")


if __name__ == "__main__":
    main()
