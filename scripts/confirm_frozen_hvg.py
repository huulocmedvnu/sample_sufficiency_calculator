#!/usr/bin/env python
"""Confirm (with real cells) that FREEZING the HVG set fixes the EmeraldBay basis instability.

Diagnostic (this session) showed EmeraldBay's 50-dim subspace only reproduces at mean cos^2 = 0.79
between two disjoint-shard fits, but when restricted to the SHARED HVGs it reproduces at ~0.92 --
i.e. the instability is HVG-selection churn, not subspace rotation. This script proves the fix
end to end: it refits two disjoint-shard bases while forcing them to use ONE frozen HVG set, and
measures the subspace agreement against the per-fit-HVG baseline.

Procedure
  1. download the committed A-shards and a disjoint set of B-shards (one at a time, deleted after)
  2. normalize_total(1e4) -> log1p each; keep A-pool and B-pool full-gene matrices
  3. BASELINE  : HVG(2000) on A-pool and on B-pool independently -> PCA(50) each -> cos^2  (~0.79)
  4. FROZEN    : HVG(2000) on the A+B pool (one frozen set) -> restrict both -> PCA(50) each -> cos^2
  5. report both, plus sigma^2 / evr reproduction, and write a JSON summary.

Env/flags: --dataset {emeraldbay,tahoe}  --n-shards N  --out JSON
"""
import os, json, argparse, numpy as np
import scipy.sparse as sp

DATASETS = {
    "emeraldbay": dict(repo="tahoebio/EmeraldBay", cache="/mnt/hdd2/loc-tran/eb_work/dl_frozen",
                       fmt="expression_data/train-{i:05d}-of-00116.parquet", n_total=116,
                       gene_width=63287, n_default=10),
    "tahoe": dict(repo="tahoebio/Tahoe-100M", cache="/mnt/hdd2/loc-tran/tahoe_work/dl_frozen",
                  fmt="data/train-{i:05d}-of-03388.parquet", n_total=3388,
                  gene_width=62713, n_default=14),
}


def committed_shards(n_total, n):
    return sorted(set(int(i) for i in np.linspace(0, n_total - 1, n).round().astype(int)))


def disjoint_shards(n_total, n, avoid):
    step = (n_total - 1) / max(n - 1, 1)
    cand = np.clip((np.linspace(0, n_total - 1, n) + step / 2).round().astype(int), 0, n_total - 1)
    out = []
    for c in map(int, cand):
        if c not in avoid and c not in out:
            out.append(c)
    j = 0
    while len(out) < n and j < n_total:
        if j not in avoid and j not in out:
            out.append(j)
        j += 1
    return sorted(out)[:n]


def load_pool(shards, cfg):
    """Download shards, normalize_total(1e4)+log1p, return one CSR pool (full-gene)."""
    import scanpy as sc, anndata as ad
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq
    os.makedirs(cfg["cache"], exist_ok=True)
    Xs = []
    for i in shards:
        p = hf_hub_download(cfg["repo"], cfg["fmt"].format(i=i), repo_type="dataset", local_dir=cfg["cache"])
        t = pq.read_table(p, columns=["genes", "expressions"])
        genes = t.column("genes").to_pylist(); expr = t.column("expressions").to_pylist()
        indptr = np.zeros(len(genes) + 1, dtype=np.int64); idx, val = [], []
        for k, (g, e) in enumerate(zip(genes, expr)):
            g = np.asarray(g, dtype=np.int64); e = np.asarray(e, dtype=np.float32)
            keep = g >= 3; g = g[keep]; e = e[keep]
            idx.append(g); val.append(e); indptr[k + 1] = indptr[k] + g.size
        Xs.append(sp.csr_matrix((np.concatenate(val), np.concatenate(idx), indptr),
                                shape=(len(genes), cfg["gene_width"])))
        try: os.remove(p)
        except OSError: pass
        print(f"    shard {i}: {Xs[-1].shape[0]} cells", flush=True)
    A = ad.AnnData(X=sp.vstack(Xs).tocsr())
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    return A


def hvg_ids(A, n=2000):
    import scanpy as sc
    B = A.copy()
    sc.pp.highly_variable_genes(B, n_top_genes=n, flavor="seurat")
    return np.where(B.var["highly_variable"].values)[0].astype(np.int64)


def pca_on(A, hvg, k=50):
    import scanpy as sc
    S = A[:, hvg].copy()
    sc.pp.pca(S, n_comps=k, zero_center=True, svd_solver="arpack")
    coords = S.obsm["X_pca"].astype(np.float64)
    return (S.varm["PCs"].T.astype(np.float64),                       # (k, len(hvg)) comps
            float(coords.var(axis=0).mean()),                        # sigma^2 marginal
            S.uns["pca"]["variance_ratio"].astype(np.float64))       # evr


def cos2(compsA, hvgA, compsB, hvgB, gene_width):
    """mean cos^2 between the two 50-dim subspaces, lifted into the shared full-gene ambient space."""
    D = gene_width
    k = min(compsA.shape[0], compsB.shape[0])
    QA = np.zeros((D, k)); QA[hvgA] = compsA[:k].T
    QB = np.zeros((D, k)); QB[hvgB] = compsB[:k].T
    s = np.clip(np.linalg.svd(QA.T @ QB, compute_uv=False), 0, 1)
    ang = np.degrees(np.arccos(s))
    return float((s ** 2).mean()), float(ang.max()), s ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="emeraldbay", choices=list(DATASETS))
    ap.add_argument("--n-shards", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    cfg = DATASETS[a.dataset]
    n = a.n_shards or cfg["n_default"]
    shA = committed_shards(cfg["n_total"], n)
    shB = disjoint_shards(cfg["n_total"], n, set(shA))
    print(f"[frozen] dataset={a.dataset}  A-shards={shA}\n[frozen] B-shards={shB}", flush=True)

    print("[frozen] loading A-pool ...", flush=True); Apool = load_pool(shA, cfg)
    print("[frozen] loading B-pool ...", flush=True); Bpool = load_pool(shB, cfg)
    print(f"[frozen] A={Apool.n_obs} cells  B={Bpool.n_obs} cells", flush=True)

    # BASELINE: independent per-fit HVG (reproduces the ~0.79 result)
    hA0, hB0 = hvg_ids(Apool), hvg_ids(Bpool)
    cA0, sA0, eA0 = pca_on(Apool, hA0); cB0, sB0, eB0 = pca_on(Bpool, hB0)
    base_c2, base_max, _ = cos2(cA0, hA0, cB0, hB0, cfg["gene_width"])
    base_overlap = int(np.intersect1d(hA0, hB0).size)

    # FROZEN: one HVG set from the A+B pool, both fits reuse it
    import anndata as ad
    pool = ad.concat([Apool, Bpool], join="outer")
    hF = hvg_ids(pool)
    cAf, sAf, eAf = pca_on(Apool, hF); cBf, sBf, eBf = pca_on(Bpool, hF)
    frz_c2, frz_max, frz_per = cos2(cAf, hF, cBf, hF, cfg["gene_width"])

    print("\n" + "=" * 64)
    print(f"{a.dataset}: n_shards/side={n}  (A={Apool.n_obs:,} cells, B={Bpool.n_obs:,} cells)")
    print("-" * 64)
    print(f"BASELINE (per-fit HVG, {base_overlap}/2000 overlap):")
    print(f"   mean cos^2 = {base_c2:.4f}   max angle = {base_max:.1f} deg")
    print(f"   sigma^2 A/B = {sA0:.4f}/{sB0:.4f}")
    print(f"FROZEN (shared HVG, 2000/2000 overlap):")
    print(f"   mean cos^2 = {frz_c2:.4f}   max angle = {frz_max:.1f} deg")
    print(f"   sigma^2 A/B = {sAf:.4f}/{sBf:.4f}")
    print(f"   directions <10deg: {(np.degrees(np.arccos(np.sqrt(frz_per))) < 10).sum()}/50")
    print(f"IMPROVEMENT: cos^2 {base_c2:.3f} -> {frz_c2:.3f}  (+{frz_c2-base_c2:.3f})")
    print("=" * 64, flush=True)

    summary = dict(dataset=a.dataset, n_shards_side=n, A_cells=int(Apool.n_obs), B_cells=int(Bpool.n_obs),
                   A_shards=shA, B_shards=shB,
                   baseline_cos2=base_c2, baseline_max_angle=base_max, baseline_hvg_overlap=base_overlap,
                   baseline_sigma2=[sA0, sB0], frozen_cos2=frz_c2, frozen_max_angle=frz_max,
                   frozen_sigma2=[sAf, sBf], frozen_evr_top5=[eAf[:5].tolist(), eBf[:5].tolist()],
                   frozen_hvg_ids=hF.tolist())
    out = a.out or f"/mnt/hdd2/loc-tran/{'eb_work' if a.dataset=='emeraldbay' else 'tahoe_work'}/out/frozen_hvg_confirm.json"
    json.dump(summary, open(out, "w"), indent=1)
    print(f"[frozen] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
