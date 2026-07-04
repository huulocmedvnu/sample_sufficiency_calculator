# Extending the directional sample-sufficiency framework to the TRADE essential-gene atlases

**Target datasets.** *Jurkat-Essential* and *HepG2-Essential* — the two new CRISPRi Perturb-seq screens of
DepMap common-essential genes introduced by the TRADE paper (Nadig, Replogle, Pogson, … Weissman,
*Nat. Genet.* 2025; bioRxiv `2024.07.03.601903`). Raw single-cell UMI counts are deposited at
**GEO GSE264667** as per-line AnnData `.h5ad`:
`GSE264667_jurkat_raw_singlecell_01.h5ad` (~8.7 GB) and `GSE264667_hepg2_raw_singlecell_01.h5ad` (~5.2 GB),
plus `GSE264667_RAW.tar` (MTX/TSV mirror).

## 0. How the two papers relate

**TRADE** estimates, for each perturbation, the *transcriptome-wide distribution of true differential-
expression effects* — correcting noisy per-gene log-fold-changes for sampling error to recover a variance
component it summarizes as the **transcriptome-wide impact** (roughly, the total true expression change a
perturbation induces). It is a statement about the **magnitude/among-gene spread** of a perturbation's effect.

**Our framework** is orthogonal and logically downstream: given that a perturbation *has* an effect, it asks
how many cells are needed to resolve the **direction** of that effect vector to an angular tolerance
`θ★` — `n★ = 2·tr(PΣP)/(m²θ★²)`. TRADE's per-perturbation impact is a magnitude scalar; our `m` is the norm
of the same displacement in the shared PCA embedding. They should correlate strongly, giving a built-in
cross-check (a perturbation TRADE calls high-impact should carry large `m` and a small `n★`). The extension
therefore adds the *directional-resolution* axis to TRADE's *magnitude* axis on the very datasets TRADE
introduced.

## 1. Confirmed data structure (from the GSE264667 h5ad headers)

Legacy-format AnnData; obs columns (verified by remote HDF5 inspection on HepG2 = **145,473 cells × 9,624
genes**; Jurkat is larger, ~240k cells):

| obs field | meaning | used for |
|---|---|---|
| `UMI_count` (float) | total UMIs per cell | low-UMI QC |
| `mitopercent` (float) | mitochondrial fraction/percent | high-mito QC (→ mito UMIs) |
| `z_gemgroup_UMI` (float) | UMI z-scored within 10x lane | alternative UMI-QC axis |
| `gem_group` (int) | 10x lane / batch | batch bookkeeping |
| `gene` (categorical) | assigned target gene (incl. non-targeting) | **condition key** |
| `gene_id`, `gene_transcript` | target IDs | — |
| `sgID_AB` (categorical) | dual-sgRNA identity | guide-assignment QC |
| `cell_barcode` | 10x barcode | — |

`var` already carries `mean/std/fano/cv/in_matrix` on a pre-filtered ~9.6k-gene set; the shared embedding is
rebuilt from raw counts regardless. **Non-targeting controls** appear as a category in `gene` (matched by
`non.?target|control|safe|scramble|negative`) and form the control arm `μ_c` — the genetic analogue of DMSO.

Because each line is ~10⁵ cells (not 10⁸ like X-Atlas/Orion), the entire matrix loads in RAM: no streaming
projection is required, and — critically — **per-cell PCA coordinates are retained**, enabling the
downsample-and-measure falsification (§4) that the Orion streaming run could not perform.

## 2. Load-bearing pipeline topology (implemented exactly)

The order of operations determines whether the anisotropic noise functional `tr(PΣP)` even exists; aggregating
to pseudobulk *before* dimensionality reduction collapses Σ to zero. `pass1_qc_basis_quota.py` enforces:

1. **Global embedding first.** After QC, on the *pooled single-cell* matrix:
   `normalize_total(1e4)` → `log(1+x)` → `highly_variable_genes(2000, seurat)` → `PCA(50)`. One shared
   `d = 50` basis for the whole line; every cell gets coordinates `x_i ∈ ℝ⁵⁰`.
2. **Cell-level covariance Σ ∈ ℝ⁵⁰ˣ⁵⁰.** From *individual embedded cells*: center each condition on its own
   sample mean, pool the residuals, `Σ = (Σ_c R_cᵀR_c)/Σ_c(n_c−1)`. We keep the **full 50×50 Σ** (not only its
   diagonal `ℓ_k`) so `tr(PΣP) = tr Σ − uᵀΣu` uses the exact within-condition covariance — a strict upgrade
   over the PCA-diagonal approximation used elsewhere in the study.
3. **Late pseudobulk.** Only now collapse the cell dimension: `μ_g` (knockdown centroids) and `μ_c` (NTC
   centroid) as coordinate averages *within* the pre-computed embedding. `m = ‖μ_g − μ_c‖`, `u = (μ_g−μ_c)/m`.

## 3. QC alignment (exact TRADE filters)

`pass1` applies per line (constants in `LINE_QC`):

| filter | Jurkat | HepG2 | implementation |
|---|---|---|---|
| low-UMI | `UMI_count < 0.14 × median` | `< 0.18 × median` | drop cells below the fraction-of-median UMI floor |
| high-mito | mito UMIs `> 1750` | `> 3000` | mito UMIs `= mitopercent(as fraction) × UMI_count` (units auto-detected: ÷100 if the column is a percent); drop above the cap |
| guide assignment | single sgRNA **or** two sgRNAs of the **same gene** | same | parse `sgID_AB` into guide tokens; drop cells whose tokens map to >1 distinct target gene |

Notes on faithful interpretation (flagged for the record):
- The "< 14% / < 18%" low-UMI rule is read as *fraction of the median UMI count* (`UMI < frac·median`). The
  `z_gemgroup_UMI` column is retained as an alternative axis if the intended rule is a within-lane z-score
  cutoff; switching is a one-line change.
- The paper's mito threshold is expressed in **UMIs** but the h5ad stores `mitopercent`, so we reconstruct
  absolute mito UMIs as `mitopercent × UMI_count`. Unit detection (percent vs fraction) is automatic.
- Guide parsing is best-effort against the `sgID_AB` encoding; since `gene` already carries a single assigned
  target per cell, the dominant effect of this filter is to remove ambiguous/mixed assignments.

## 4. Quota spectrum and the falsification (genetic Phase C)

**Quota** (`pass1`, per knockdown vs NTC): `m` bias-corrected for the finite-sample floor
`tr(S)=trΣ(1/n_g+1/n_NTC)`; `n★_iso = 2(d−1)σ²/(m²θ★²)`, `n★_aniso = 2·tr(PΣP)/(m²θ★²)`
(equal-arm factor 2, matching the Tahoe/Orion headline; the large shared-NTC-pool variant halves it and is
reported alongside). Each knockdown is classed OVER / UNDER / Ghost against its acquired cell count; the
detectable subset (signal > 1.5× the sampling floor) is reported separately so the weak-majority estimates
are not over-read.

**Falsification** (`pass2`, the genetic counterpart of manuscript Phase C — newly *possible* here): for each
well-sampled knockdown, take the ground-truth direction from all `N` cells vs NTC, subsample `n ≤ N/2` cells
without replacement (`R` reps), and regress realized `θ²(n)` on `(1/n − 1/N)`. The slope is the **a-priori,
parameter-free** `tr(PΣP)/m²`; we compare fitted-vs-predicted and stratify by the theory's validity criterion
`ρ² = m²/(uᵀΣu) ≥ 3`, reproducing Phase C's in-regime-match / low-SNR-breakdown structure on genetic data.

## 5. Scalability

- **In-memory, no streaming.** ~10⁵ cells × 2,000 HVGs fits comfortably in RAM; PCA(50) via `arpack`. Full
  50×50 Σ and all per-condition sufficient statistics are trivial. Wall-clock is minutes per line, dominated
  by the one-time h5ad load. (Contrast Orion: 10⁸ cells forced a resumable streaming projection.)
- **Reproducible + parameterized.** `run_trade.sh` runs both lines end to end; QC constants, `θ★`, and the
  min-cell thresholds are CLI/`LINE_QC` parameters. `pass1` introspects obs column names, so it survives minor
  schema differences between the two files.
- **Deterministic.** Fixed HVG/PCA recipe and seeded subsampling; outputs are committed fixtures.

## 6. Deliverables & how they slot into the study

- `scripts/trade_recompute/` — `pass0_download.py`, `pass1_qc_basis_quota.py`, `pass2_falsification.py`,
  `run_trade.sh`.
- `fixtures/trade_{jurkat,hepg2}_quota.csv` + `summary.json`; `fixtures/trade_{jurkat,hepg2}_falsification.json`.
- Comparison against the chemical (Tahoe/EmeraldBay) and CRISPRi (X-Atlas/Orion) results: a third laboratory,
  a fourth/fifth cell line, and — with the full-Σ anisotropic quota plus a genetic Phase-C falsification — the
  first end-to-end *validated* gene-perturbation application. If σ² again lands near ~1.0 and the parameter-free
  slope matches in-regime, this closes the cross-modality generalization with a proper falsification, not just
  an application. Expected (to be confirmed by the run): essential-gene knockdowns are strong perturbations, so
  the magnitude distribution should sit **above** Orion's (median `m` larger), making a larger over-sampled
  fraction plausible than in the genome-wide Orion screens.

## 7. Risks / caveats

1. **QC-threshold interpretation** (low-UMI fraction-of-median; mito-percent→UMI reconstruction) is made
   explicit above and is one-line switchable if the paper's exact rule differs.
2. **Guide encoding** in `sgID_AB` may not split cleanly on the assumed delimiters; the loader logs the guide
   pass-rate so an anomaly is visible, and `gene`'s single assigned target is the fallback.
3. **Baseline choice.** The quota references the NTC centroid (large pool ≈ near-noiseless); the falsification
   uses the same NTC baseline, so its `m`-values are NTC-referenced (stated, matching Phase C's baseline note).
4. **Factor-2 vs large-pool.** With a big NTC pool the true control-arm noise is small; both the equal-arm and
   large-pool quota constants are reported so the convention is explicit.
5. **Scope unchanged.** `n★` governs centroid-direction precision only; magnitudes are transcriptional
   displacements of surviving cells; compute savings are polynomial.
