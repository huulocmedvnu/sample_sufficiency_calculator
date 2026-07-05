#!/usr/bin/env python
"""Principal-angle stability check for the subsample-fit PCA basis.

The recompute pipelines fit HVG(2000)+PCA(50) on a small, plate-spread subsample of shards and then
project *all* cells through that fixed basis (out-of-sample projection). That is only a good embedding
if the 50-dim subspace is STABLE under resampling -- i.e. a basis fit on a different, disjoint set of
shards spans (nearly) the same subspace. This script quantifies that by the principal angles between
the two 50-dim subspaces.

Two bases select their OWN 2000 HVGs, so the comparison is done in the shared full-gene ambient space:
each (50 x 2000) loading matrix is scattered into a (D x 50) column-orthonormal matrix Q at its own
gene-token columns (genes chosen by only one basis simply contribute zeros on the other side). Then for
Q_A, Q_B with orthonormal columns, the singular values of Q_A^T Q_B are the cosines of the principal
angles theta_1 <= ... <= theta_50.

Reported:
  * HVG overlap |hvg_A cap hvg_B| / 2000
  * principal angles (deg), the largest angle, and how many directions align within 10/25/45 deg
  * mean cos^2 = tr(P_A P_B)/k in [0,1]  (1 = identical subspace; the single headline number)
  * Grassmann geodesic + chordal distances
  * top-5 evr and sigma^2 side by side

Modes
-----
compare : two already-fit basis.npz files (no download). e.g. Orion pilot vs full for the same line.
    python scripts/check_basis_stability.py compare \
        --basis-a /mnt/hdd2/loc-tran/orion_work/pilot/HCT116/basis.npz \
        --basis-b /mnt/hdd2/loc-tran/orion_work/full/HCT116/basis.npz

fit     : load the committed basis A, fit a fresh basis B on a DISJOINT plate-spread set of shards
          (Tahoe / EmeraldBay), then compare. Downloads N_SHARDS shards one at a time (deleted after).
    python scripts/check_basis_stability.py fit --dataset tahoe --n-shards 14
    python scripts/check_basis_stability.py fit --dataset emeraldbay --n-shards 10
"""
import os, argparse, numpy as np


# ---------------------------------------------------------------- basis IO / geometry

def load_basis(path):
    """Return a normalized dict: comps (k,h), hvg (h,), evr, sigma2, gene_width, n_cells, src."""
    b = np.load(path, allow_pickle=True)
    keys = set(b.keys())
    hvg = b["hvg_token_ids"] if "hvg_token_ids" in keys else b["hvg_ids"]
    sig = ("sigma2" if "sigma2" in keys else
           "sigma2_within" if "sigma2_within" in keys else
           "sigma2_marginal" if "sigma2_marginal" in keys else None)
    comps = b["components"].astype(np.float64)                     # (k, h)
    hvg = np.asarray(hvg).astype(np.int64)
    gw = int(b["gene_width"]) if "gene_width" in keys else int(hvg.max()) + 1
    return {
        "comps": comps, "hvg": hvg, "gene_width": gw,
        "evr": np.asarray(b["evr"]).astype(np.float64) if "evr" in keys else None,
        "sigma2": float(b[sig]) if sig else float("nan"),
        "sigma2_key": sig,
        "n_cells": int(b["n_cells_fit"]) if "n_cells_fit" in keys else -1,
        "src": path,
    }


def lift(comps, hvg, D):
    """Scatter (k,h) loadings into a (D,k) matrix with columns orthonormal within this basis."""
    Q = np.zeros((D, comps.shape[0]), dtype=np.float64)
    Q[hvg, :] = comps.T                                           # comps.T is (h,k)
    return Q


def principal_angles(QA, QB):
    s = np.linalg.svd(QA.T @ QB, compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    return s, np.degrees(np.arccos(s))                            # cosines, angles(deg) ascending-angle


def report(A, B):
    kA, kB = A["comps"].shape[0], B["comps"].shape[0]
    if kA != kB:
        print(f"[warn] different n_comps: A={kA} B={kB}; comparing min")
    k = min(kA, kB)
    D = max(A["gene_width"], B["gene_width"], int(A["hvg"].max()) + 1, int(B["hvg"].max()) + 1)
    QA = lift(A["comps"][:k], A["hvg"], D)
    QB = lift(B["comps"][:k], B["hvg"], D)

    # orthonormality sanity (scatter should preserve it)
    oaA = np.abs(QA.T @ QA - np.eye(k)).max()
    oaB = np.abs(QB.T @ QB - np.eye(k)).max()

    cos, ang = principal_angles(QA, QB)
    cos2 = cos ** 2
    mean_cos2 = float(cos2.mean())                               # tr(P_A P_B)/k
    grassmann = float(np.sqrt((np.radians(ang) ** 2).sum()))     # geodesic distance
    chordal = float(np.sqrt((np.sin(np.radians(ang)) ** 2).sum()))
    inter = np.intersect1d(A["hvg"], B["hvg"]).size

    print("=" * 72)
    print(f"A: {A['src']}")
    print(f"   n_cells_fit={A['n_cells']:,}  sigma2({A['sigma2_key']})={A['sigma2']:.4f}")
    print(f"B: {B['src']}")
    print(f"   n_cells_fit={B['n_cells']:,}  sigma2({B['sigma2_key']})={B['sigma2']:.4f}")
    print("-" * 72)
    print(f"ambient gene dim D                 : {D:,}")
    print(f"HVG overlap                        : {inter}/{A['comps'].shape[1]} "
          f"({100*inter/A['comps'].shape[1]:.1f}%)")
    print(f"orthonormality residual (A,B)      : {oaA:.2e}, {oaB:.2e}")
    print("-" * 72)
    print(f"MEAN cos^2 (subspace overlap, 0..1): {mean_cos2:.4f}   <-- headline")
    print(f"largest principal angle            : {ang.max():.2f} deg  (cos={cos.min():.4f})")
    print(f"median principal angle             : {np.median(ang):.2f} deg")
    print(f"directions aligned  < 10 deg       : {(ang < 10).sum()}/{k}")
    print(f"directions aligned  < 25 deg       : {(ang < 25).sum()}/{k}")
    print(f"directions aligned  < 45 deg       : {(ang < 45).sum()}/{k}")
    print(f"Grassmann geodesic dist (rad)      : {grassmann:.3f}")
    print(f"chordal dist (sqrt sum sin^2)      : {chordal:.3f}")
    if A["evr"] is not None and B["evr"] is not None:
        print("-" * 72)
        print(f"top-5 evr A: {np.round(A['evr'][:5], 4).tolist()}")
        print(f"top-5 evr B: {np.round(B['evr'][:5], 4).tolist()}")
    print("-" * 72)
    print("principal angles (deg), ascending:")
    print("  " + "  ".join(f"{a:5.1f}" for a in ang))
    print("=" * 72)
    print("interpretation: mean cos^2 > ~0.9 and largest angle small => the subsample basis is")
    print("stable under resampling and the fixed out-of-sample projection is well justified.")
    print("A large first-angle jump with high overlap means a specific direction is unstable;")
    print("low HVG overlap warns the instability is in feature selection, not just rotation.")
    return {"mean_cos2": mean_cos2, "max_angle": float(ang.max()), "hvg_overlap": int(inter)}


# ---------------------------------------------------------------- fit mode (Tahoe / EmeraldBay)

DATASETS = {
    "tahoe": dict(repo="tahoebio/Tahoe-100M", cache="/mnt/hdd2/loc-tran/tahoe_work/dl_stab",
                  fmt="data/train-{i:05d}-of-03388.parquet", n_total=3388, gene_width=62713,
                  basis="/mnt/hdd2/loc-tran/tahoe_work/out/basis.npz", n_default=14),
    "emeraldbay": dict(repo="tahoebio/EmeraldBay", cache="/mnt/hdd2/loc-tran/eb_work/dl_stab",
                       fmt="expression_data/train-{i:05d}-of-00116.parquet", n_total=116, gene_width=63287,
                       basis="/mnt/hdd2/loc-tran/eb_work/out/basis.npz", n_default=10),
}


def committed_fit_shards(n_total, n):
    return sorted(set(int(i) for i in np.linspace(0, n_total - 1, n).round().astype(int)))


def disjoint_shards(n_total, n, avoid):
    """A plate-spread set of n shards that avoids `avoid` (the committed fit shards)."""
    step = (n_total - 1) / max(n - 1, 1)
    cand = np.clip((np.linspace(0, n_total - 1, n) + step / 2).round().astype(int), 0, n_total - 1)
    out = []
    for c in map(int, cand):
        if c not in avoid and c not in out:
            out.append(c)
    j = 0
    while len(out) < n and j < n_total:                          # top up if edge collisions dropped some
        if j not in avoid and j not in out:
            out.append(j)
        j += 1
    return sorted(out)[:n]


def shard_to_csr(path, gene_width):
    import scipy.sparse as sp, pyarrow.parquet as pq
    t = pq.read_table(path, columns=["genes", "expressions"])
    genes = t.column("genes").to_pylist(); expr = t.column("expressions").to_pylist()
    indptr = np.zeros(len(genes) + 1, dtype=np.int64); idx, val = [], []
    for i, (g, e) in enumerate(zip(genes, expr)):
        g = np.asarray(g, dtype=np.int64); e = np.asarray(e, dtype=np.float32)
        keep = g >= 3; g = g[keep]; e = e[keep]
        idx.append(g); val.append(e); indptr[i + 1] = indptr[i] + g.size
    return sp.csr_matrix((np.concatenate(val), np.concatenate(idx), indptr),
                         shape=(len(genes), gene_width))


def fit_basis(shards, cfg):
    """Replicate the pass1 recipe exactly on `shards`; return a basis dict compatible with load_basis."""
    import scipy.sparse as sp, scanpy as sc, anndata as ad
    from huggingface_hub import hf_hub_download
    os.makedirs(cfg["cache"], exist_ok=True)
    Xs = []
    for i in shards:
        p = hf_hub_download(cfg["repo"], cfg["fmt"].format(i=i), repo_type="dataset", local_dir=cfg["cache"])
        X = shard_to_csr(p, cfg["gene_width"]); Xs.append(X)
        try: os.remove(p)
        except OSError: pass
        print(f"  shard {i}: {X.shape[0]} cells", flush=True)
    A = ad.AnnData(X=sp.vstack(Xs).tocsr())
    print(f"[fit] subsample: {A.n_obs} cells x {A.n_vars} genes", flush=True)
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat")
    m = A.var["highly_variable"].values
    hvg = np.where(m)[0].astype(np.int64); A = A[:, m].copy()
    sc.pp.pca(A, n_comps=50, zero_center=True, svd_solver="arpack")
    coords = A.obsm["X_pca"].astype(np.float64)
    return {"comps": A.varm["PCs"].T.astype(np.float64), "hvg": hvg, "gene_width": cfg["gene_width"],
            "evr": A.uns["pca"]["variance_ratio"].astype(np.float64),
            "sigma2": float(coords.var(axis=0).mean()), "sigma2_key": "marginal",
            "n_cells": A.n_obs, "src": f"fresh-fit shards={shards}"}


# ---------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description="Principal-angle stability check for a subsample-fit PCA basis")
    sub = ap.add_subparsers(dest="mode", required=True)

    c = sub.add_parser("compare", help="compare two existing basis.npz files")
    c.add_argument("--basis-a", required=True); c.add_argument("--basis-b", required=True)

    f = sub.add_parser("fit", help="fit a fresh basis on disjoint shards and compare to the committed one")
    f.add_argument("--dataset", required=True, choices=list(DATASETS))
    f.add_argument("--basis-a", default=None, help="committed basis (default: dataset's out/basis.npz)")
    f.add_argument("--n-shards", type=int, default=None)
    f.add_argument("--shards-b", default=None, help="comma-sep shard indices to override the disjoint set")
    f.add_argument("--save-b", default=None, help="optional path to np.savez the fresh basis B")

    a = ap.parse_args()
    if a.mode == "compare":
        report(load_basis(a.basis_a), load_basis(a.basis_b))
        return

    cfg = DATASETS[a.dataset]
    n = a.n_shards or cfg["n_default"]
    basis_a = a.basis_a or cfg["basis"]
    A = load_basis(basis_a)
    committed = committed_fit_shards(cfg["n_total"], n)
    if a.shards_b:
        shards_b = sorted(int(x) for x in a.shards_b.split(","))
    else:
        shards_b = disjoint_shards(cfg["n_total"], n, set(committed))
    overlap = set(committed) & set(shards_b)
    print(f"[fit] committed A shards: {committed}")
    print(f"[fit] fresh    B shards: {shards_b}")
    print(f"[fit] shard overlap A/B : {sorted(overlap) if overlap else 'none (disjoint)'}", flush=True)
    B = fit_basis(shards_b, cfg)
    if a.save_b:
        np.savez(a.save_b, components=B["comps"], hvg_token_ids=B["hvg"], evr=B["evr"],
                 sigma2=B["sigma2"], n_cells_fit=B["n_cells"], gene_width=cfg["gene_width"],
                 n_comps=50, target_sum=1e4, shards=np.array(shards_b))
        print(f"[fit] saved B -> {a.save_b}")
    report(A, B)


if __name__ == "__main__":
    main()
