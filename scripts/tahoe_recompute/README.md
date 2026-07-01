# Tahoe-100M recompute (fresh, from the raw dataset)

Recomputes the sample-sufficiency inputs **from the raw Tahoe-100M single-cell counts**
(`tahoebio/Tahoe-100M` on Hugging Face, 100M+ cells, 379 drugs × 50 cell lines, 3,388 expression
shards), following the theislab reference workflow
(<https://theislab.github.io/vevo_Tahoe_100m_analysis/vevo_100m_pca.html>). It replaces the earlier
cached/derived results.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| Pass 1 | `pass1_basis.py` | Fit the embedding on a multi-plate subsample: `normalize_total(1e4)` → `log1p` → `highly_variable_genes(n_top_genes=2000)` → `PCA(50)`. Emits `basis.npz` (HVG token ids, PCA components + mean, **σ²** = mean per-PC per-cell variance). |
| Pass 2 | `pass2_pseudobulk.py` | Stream all 3,388 shards; apply the same per-cell transform; accumulate the per-(drug × cell_line) pseudobulk centroid in HVG log-CP10K space. Resumable (checkpoints every 50 shards); one shard on disk at a time. Emits `pseudobulk.npz`. |
| Pass 3 | `pass3_quota.py` | Project centroids to the PCA basis; `v = coord(drug,line) − coord(DMSO_TF,line)`, `m = ‖v‖`, `n* = 2(d−1)σ²/(m²θ²)` (θ=0.1). Classify each condition OVER / UNDER / Ghost against its actual cell count, and aggregate the **per-cell-line** UNDER/Ghost breakdown. Emits `quota_per_condition.csv`, `per_cell_line.csv`, `per_drug.csv`, `summary.json`. |

Run (from repo root):
```bash
export OUT=/mnt/hdd2/loc-tran/tahoe_work/out
python scripts/tahoe_recompute/pass1_basis.py        # minutes
python scripts/tahoe_recompute/pass2_pseudobulk.py   # hours (streams ~337 GB); resumable
python scripts/tahoe_recompute/pass3_quota.py         # seconds
```

## Faithful-recipe notes / documented deviations

* **normalize_total target_sum = 1e4** (fixed CP10K) so the identical transform applies shard-by-shard
  without a global median pass (the tutorial's default normalizes to the per-run median library size).
* **PCA dimension d = 50** — the sample-sufficiency framework's embedding dimension; the tutorial uses
  `n_comps=300` for UMAP visualisation. σ² and n* are reported in this 50-d space.
* **HVG(2000)** computed on a multi-plate subsample (the tutorial does per-plate HVG then keeps genes
  appearing in >2 plates); special gene tokens (id < 3) are dropped; centroids use the linear
  PCA projection of the per-condition mean (exact, since PCA is affine).
* Baseline `N0` for OVER/UNDER is each condition's **actual acquired cell count**; "Ghost" = n* > 50,000.

Outputs are written outside the repo (`$OUT`, default `/mnt/hdd2/loc-tran/tahoe_work/out`); only the
distilled numbers/fixtures are committed.
