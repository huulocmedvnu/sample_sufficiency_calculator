#!/usr/bin/env python
"""Pass 1 of the EmeraldBay recompute: fit EmeraldBay's OWN HVG(2000)+PCA(50) embedding on a
multi-shard subsample (same theislab-style recipe as the Tahoe recompute). EmeraldBay is the
independent validation atlas (tahoebio/EmeraldBay, 1.83M cells, 116 shards, 57.7 GB).

Outputs $OUT_EB/basis.npz : hvg_token_ids (2000,), components (50,2000), pca_mean (2000,),
        sigma2_marginal, per_pc_var, evr, n_cells_fit.
(The WITHIN-condition sigma^2 -- the one the theory uses -- is computed in Pass 2 from per-condition
residuals.)
"""
import os, numpy as np, scipy.sparse as sp
import pyarrow.parquet as pq
from huggingface_hub import hf_hub_download
import scanpy as sc, anndata as ad

REPO = "tahoebio/EmeraldBay"
CACHE = "/mnt/hdd2/loc-tran/eb_work/dl"
OUT = os.environ.get("OUT_EB", "/mnt/hdd2/loc-tran/eb_work/out")
N_TOTAL = 116
N_SHARDS_FIT = int(os.environ.get("N_SHARDS_FIT", "10"))
GENE_WIDTH = 63287          # max token_id 63286 + 1; special tokens 0,1,2 dropped
TARGET_SUM = 1e4
N_COMPS = 50
os.makedirs(CACHE, exist_ok=True); os.makedirs(OUT, exist_ok=True)


def shard_to_csr(path):
    t = pq.read_table(path, columns=["genes", "expressions", "drug", "cell_line"])
    genes = t.column("genes").to_pylist(); expr = t.column("expressions").to_pylist()
    drug = np.array(t.column("drug").to_pylist()); line = np.array(t.column("cell_line").to_pylist())
    indptr = np.zeros(len(genes) + 1, dtype=np.int64); idx, val = [], []
    for i, (g, e) in enumerate(zip(genes, expr)):
        g = np.asarray(g, dtype=np.int64); e = np.asarray(e, dtype=np.float32)
        keep = g >= 3; g = g[keep]; e = e[keep]
        idx.append(g); val.append(e); indptr[i + 1] = indptr[i] + g.size
    X = sp.csr_matrix((np.concatenate(val), np.concatenate(idx), indptr), shape=(len(genes), GENE_WIDTH))
    return X, drug, line


def main():
    idxs = sorted(set(int(i) for i in np.linspace(0, N_TOTAL - 1, N_SHARDS_FIT).round()))
    print(f"[pass1-eb] fitting basis on {len(idxs)} shards: {idxs}")
    Xs, drugs, lines = [], [], []
    for i in idxs:
        f = f"expression_data/train-{i:05d}-of-00116.parquet"
        p = hf_hub_download(REPO, f, repo_type="dataset", local_dir=CACHE)
        X, d, l = shard_to_csr(p); Xs.append(X); drugs.append(d); lines.append(l); os.remove(p)
        print(f"  shard {i}: {X.shape[0]} cells")
    X = sp.vstack(Xs).tocsr()
    A = ad.AnnData(X=X); A.obs["drug"] = np.concatenate(drugs); A.obs["cell_line"] = np.concatenate(lines)
    print(f"[pass1-eb] subsample: {A.n_obs} cells x {A.n_vars} genes")
    sc.pp.normalize_total(A, target_sum=TARGET_SUM); sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat")
    hvg = np.where(A.var["highly_variable"].values)[0].astype(np.int64)
    A = A[:, A.var["highly_variable"].values].copy()
    sc.pp.pca(A, n_comps=N_COMPS, zero_center=True, svd_solver="arpack")
    comps = A.varm["PCs"].T.astype(np.float64)
    pca_mean = np.asarray(A.X.mean(axis=0)).ravel().astype(np.float64)
    coords = A.obsm["X_pca"].astype(np.float64)
    sigma2_marg = float(coords.var(axis=0).mean())
    np.savez(os.path.join(OUT, "basis.npz"), hvg_token_ids=hvg, components=comps, pca_mean=pca_mean,
             sigma2_marginal=sigma2_marg, per_pc_var=coords.var(axis=0), evr=A.uns["pca"]["variance_ratio"],
             n_cells_fit=A.n_obs, target_sum=TARGET_SUM, n_comps=N_COMPS, gene_width=GENE_WIDTH)
    print(f"[pass1-eb] HVG={A.n_vars}, sigma^2(marginal)={sigma2_marg:.4f}, top5 evr={A.uns['pca']['variance_ratio'][:5].round(4).tolist()}")
    print(f"[pass1-eb] wrote {OUT}/basis.npz")


if __name__ == "__main__":
    main()
