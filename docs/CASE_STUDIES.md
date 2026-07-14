> ⚠️ **SUPERSEDED.** Numbers here predate the audit (findings #1-#9, gaps 1-3). The authoritative source is `docs/MANUSCRIPT_DEEPSEEK.md` + `src/engine.py`. See `docs/AUDIT_DENOMINATORS.md` and `docs/ARCHITECTURE.md`.

# Empirical Case Studies — Sample-Sufficiency Spectrum (frozen source-of-truth)

**Configuration (Tahoe-100M, from-raw recompute; `scripts/tahoe_recompute/`).** σ² = 0.9567
(within-condition residual variance — the theory-preferred headline), d = 50,
θ★ = 0.1 rad (5.73°), equal arms. The governing law is **n★ = 2(d−1)σ²/(m²θ★²) = 9,376 / m²**. At the
median depth N₀ = 1,296 the over/under boundary falls at **m = 2.69**, and the depth-fixed angular
resolution is **θ(N₀) = 0.269 / m rad**. (The marginal per-cell variance σ² = 2.406 — pooling scatter
across conditions — yields a ~2.5× larger, more conservative quota, retained only as a cross-check.) The unit of analysis is the **(drug × dose × cell-line)
condition**: each magnitude m is the length ‖v̂‖ of the perturbation vector v̂ = μ̂_t − μ̂_c, the
difference of the treated and plate-matched-DMSO **centroids** (pseudobulk mean *vectors*, each a
coordinate average of that arm's cells) in the shared PCA(50) space where σ² was estimated. The two arms
usually have different cell counts, so this is a subtraction of two d-vectors, never of the raw
cell-by-gene matrices. All values are machine-derived; full tables are in `fixtures/tahoe_*.csv`.

## Per-drug spectrum (median across 150 conditions = 3 doses × 50 lines)

| Drug (MoA / target) | median m | median n★ (cells/arm) | OVER | UNDER | Ghost |
|---|---:|---:|---:|---:|---:|
| Panobinostat (HDAC) | 4.24 | 523 | 120 | 30 | 0 |
| Homoharringtonine (protein synthesis) | 5.88 | 271 | 100 | 50 | 0 |
| Harringtonine (protein synthesis) | 3.67 | 697 | 78 | 72 | 0 |
| Idarubicin (anthracycline; TOP2A) | 2.93 | 1,091 | 70 | 80 | 0 |
| Trametinib (MEK) | 2.29 | 1,789 | 86 | 64 | 0 |
| Palbociclib (CDK4/6) | 1.41 | 4,684 | 4 | 146 | 0 |
| 4EGI-1 (eIF4E) | 1.34 | 5,250 | 4 | 146 | 0 |
| Crizotinib (ALK/MET) | 0.99 | 9,632 | 1 | 146 | 3 |

The strongest compounds are now over-sampled in a majority of their 150 conditions — panobinostat in 120
and homoharringtonine in 100 — while at the other extreme only **28 of 379 drugs (7%) are over-sampled in
no condition at all**. Full per-drug counts are in `fixtures/tahoe_per_drug.csv`.

## Two levers beyond drug identity — worked conditions

Magnitude, and therefore the quota, is set not by the compound alone but jointly by dose and cell-line
background. Two worked examples make this concrete.

**Dose** (homoharringtonine in the responsive NCI-H460 line): raising the concentration enlarges m and
lowers the quota.

| condition | m | n★ | N₀ | regime | dry save (lin / quad) | multiplex |
|---|---:|---:|---:|:--:|:--:|:--:|
| Homoharringtonine 0.05 µM × NCI-H460 | 10.3 | 88 | 2,904 | OVER | 97% / 100% | 33.0× |
| Homoharringtonine 0.5 µM × NCI-H460 | 14.5 | 44 | 2,114 | OVER | 98% / 100% | 48.0× |
| Homoharringtonine 5 µM × NCI-H460 | 15.5 | 39 | 6,060 | OVER | 99% / 100% | 155.4× |

**Cell line** (same drug and dose, different background): homoharringtonine at 5 µM reaches m = 15.5
(n★ = 39, OVER) in NCI-H460 but only m = 3.1 (n★ = 951, UNDER) in NCI-H661 — the same treatment flips
regime with the background.

## Analysis

The quota is governed by squared transcriptional potency, n★ = 9,376/m², so a condition's perturbation
magnitude alone fixes its sampling regime, with the over/under boundary at m ≈ 2.69 (at the median depth).
The consequences run in two directions.

At the strong end, high-magnitude cytotoxics displace cells far off baseline — the protein-synthesis
inhibitors (homoharringtonine, harringtonine) and the HDAC inhibitor panobinostat. These saturate their
directional estimate in a few hundred cells and are over-sampled in the majority of their conditions,
which leaves those wells over-provisioned and licenses large downstream-compute reductions or multiplexing
via cell hashing.

At the weak end, n★ rises quadratically as potency falls — halving m quadruples the requirement — and at
fixed depth the achievable resolution degrades linearly as θ = 0.269/m. Weak targeted agents such as
palbociclib (CDK4/6, median m = 1.41) and crizotinib (ALK/MET, median m = 0.99) carry faint
transcriptional footprints: palbociclib is under-sampled in 146 of its 150 conditions, and crizotinib is a
"ghost" (n★ > 50,000) in only 3.

The practical lesson is that magnitude is not predicted by nominal mechanism: targeted kinase inhibitors
sit at the weak end while broad cytotoxics dominate the strong end. Sample sufficiency must therefore be
calibrated per condition from the empirical perturbation magnitude — and jointly by dose and cell line —
never inferred from drug-class labels.

## Population context (do not over-read the strong tail)

The rows above are **illustrative tiers, not a typical sample**. Across the full 56,827-condition panel
(379 drugs × 3 doses × 50 lines) at the same configuration, the **median n★ is 5,794 cells** (IQR
2,581–10,998); **10.6%** of conditions are over-sampled (n★ < N₀), 28% require 10 k–50 k, and 0.4% are
ghosts (> 50 k) — **89.4% under-sampled or worse**. At a tight ~5.7° tolerance the *typical* condition is
therefore **under-sampled**; the over-sampled cases are the high-magnitude / high-dose / responsive-line
minority. See [`SCALE_AUDIT.md`](SCALE_AUDIT.md) for the full distribution and the project budget.

## Provenance / caveats

- All m and n★ are machine-computed from the raw Tahoe-100M counts through the streaming pipeline
  (`scripts/tahoe_recompute/`) and the shipped calculator; no placeholders.
- Pipeline order (important): the PCA(50) embedding is fitted once on individual cells pooled across the
  atlas; σ² is the within-condition cell-to-cell scatter measured from single cells *in that space*
  (residual variance after centering each condition on its own centroid; the marginal cross-condition
  scatter, 2.406, is the conservative bound); the centroids μ̂ are
  then formed **downstream** as a coordinate average of embedded cells — a within-embedding average, not
  a gene-count aggregation *before* PCA (which would leave one point per condition, zero out σ², and
  destroy the anisotropic trace tr(PΣP) the quota is built on).
- Magnitudes are 24 h survivor transcriptional norms (not viability), specific to this
  platform/timepoint/embedding; the plate-matched DMSO_TF vehicle is the control.
- The condition unit is (drug × dose × cell line); pooling doses inflates N₀ and averages magnitudes, and
  is avoided here.
