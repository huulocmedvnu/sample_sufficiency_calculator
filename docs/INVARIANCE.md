# Downstream Functional Invariance (empirical demonstration)

**Claim.** If the calculator diagnoses a well as OVER-sampled (n★ ≪ N₀), downsampling it from N₀ to its
quota n★ leaves **centroid/pseudobulk** downstream metrics — cosine correlation with other drugs, and
Euclidean centroid distance from the vehicle — virtually invariant. Performing the *same* downsampling
on an UNDER-sampled well (whose signal has not saturated) makes those metrics fluctuate wildly.

**Substrate & scope.** Per-cell data from **tahoebio/EmeraldBay** (the Tahoe fixtures store only
centroids and cannot be downsampled). Vehicle baseline = per-line mean (no matched-timepoint DMSO).
Anisotropic n★ at θ★ = 0.1 rad, large-control-pool form. **n★ governs centroid metrics only** (cosine,
Euclidean centroid distance); it does **not** govern distributional **E-distance** (energy distance),
which depends on the cell-cloud shape the centroid discards and obeys a different sufficiency criterion.

## Result (R = 300 downsampling replicates; both wells downsampled to the same depth n_d = 67)

| metric | **OVER: DMSO_T0 × HS-578T** (m = 12.41) | **UNDER: GDC-6036 × HEC-1-A** (m = 0.40) |
|---|---|---|
| N₀ cells / well | 1,067 | 1,518 |
| n★ (cells/arm) | **67** (OVER, 16× spare) | **65,066** (UNDER, 43× short) |
| downsample N₀ → | 67 ( = n★) | 67 ( ≪ n★ ) |
| Euclidean centroid distance from vehicle: full → down | 12.410 → 12.524 ± 0.567 (**Δ = 0.11, 0.9%**) | 0.398 → 1.188 ± 0.194 (**Δ = 0.79, +199%**) |
| cosine vs nearest neighbour #1: full → down | Cetuximab +0.170 → +0.170 ± 0.021 (**Δ = 0.000**) | Sotorasib +0.493 → +0.167 ± 0.185 (**Δ = 0.326**) |
| cosine vs neighbour #2 | Panitumumab +0.051 → +0.052 ± 0.018 (Δ = 0.001) | Cetuximab +0.415 → +0.127 ± 0.184 (Δ = 0.288) |
| cosine vs neighbour #3 | Adagrasib +0.023 → +0.025 ± 0.018 (Δ = 0.002) | Panitumumab +0.333 → +0.099 ± 0.251 (Δ = 0.234) |
| replicate std (stability) | ≈ **0.02** (cosine), 0.57 (distance) | ≈ **0.20** (cosine), 0.19 (distance) |

## Analysis

The OVER-sampled centroid is functionally saturated: a **16-fold** reduction (1,067 → 67 cells, to its
own quota) shifts the vehicle distance by **0.9%** and the cosine correlation with each of its three
nearest neighbours by **≤ 0.002**, with replicate scatter of only ≈0.02 — the downstream similarity graph
is, for practical purposes, unchanged. Subjecting the UNDER-sampled drug to the *identical* 67-cell
downsampling collapses it: the vehicle distance more than triples (Δ = +199%) and the top-neighbour
cosine drops from +0.49 to +0.17 ± 0.19, i.e. the correlation is both biased and unstable, with replicate
scatter (≈0.20) an order of magnitude larger than the OVER case. This is the empirical confirmation of
the dry-lab ledger: downsampling to n★ is a *provably safe* compression for over-acquired wells but a
*destructive* one for under-acquired wells — and the calculator distinguishes the two a priori from
magnitude alone.

## Auditor notes (do not omit from the manuscript)

- **The OVER exemplar is a reference population (DMSO_T0), not a drug.** In EmeraldBay's representative-line subset *no
  real drug is cleanly over-sampled* — Paclitaxel itself is under-sampled there (N₀ ≈ 107 cells/representative-line,
  m = 3.3, n★ = 942). The invariance is a geometric property of any high-magnitude, well-sampled centroid.
- **The UNDER exemplar (GDC-6036, a KRAS-G12C inhibitor) is a "ghost" with biological cause:** HEC-1-A
  carries **KRAS G12D**, not G12C, so on-target activity is not expected — m = 0.40, n★ = 65,066. A
  targeted agent against the wrong genotype is exactly the case the calculator flags as un-resolvable.
- **Resveratrol is absent from EmeraldBay** (its 27-drug panel is colorectal/RAS/HER2); the UNDER drug is
  GDC-6036, the closest available real under-sampled agent.
- Metrics are 24 h-independent centroid quantities in EmeraldBay's own HVG/PCA space; E-distance is out of
  scope (see Scope above). Source: the companion calibration pipeline →
  `outputs/emeraldbay_invariance.json`.
