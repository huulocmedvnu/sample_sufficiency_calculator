# Downstream Functional Invariance (empirical demonstration)

**Claim.** If the calculator diagnoses a well as OVER-sampled (n★ ≪ N₀), downsampling it from N₀ to its
quota n★ leaves **centroid/pseudobulk** downstream metrics — cosine correlation with other drugs, and
Euclidean centroid distance from the vehicle — virtually invariant. Performing the *same* downsampling
on an UNDER-sampled well (whose signal has not saturated) makes those metrics fluctuate wildly.

**Substrate & scope.** Per-cell data from **tahoebio/EmeraldBay** (the Tahoe fixtures store only
centroids and cannot be downsampled), in EmeraldBay's own **frozen** HVG(2000)+PCA(50) space (the
reproducibility fix; see the manuscript Methods). Vehicle baseline = per-line mean (no matched-timepoint
DMSO). Anisotropic n★ at θ★ = 0.1 rad, n★ = 2 tr(PΣP)/(m²θ★²). **n★ governs centroid metrics only**
(cosine, Euclidean centroid distance); it does **not** govern distributional **E-distance** (energy
distance), which depends on the cell-cloud shape the centroid discards and obeys a different sufficiency
criterion. Regenerated on the committed (frozen-HVG) basis by `scripts/emeraldbay_recompute/pass6_invariance.py`.

## Result (R = 300 downsampling replicates; both wells downsampled to the same depth n_d = 184)

| metric | **OVER: DMSO_T0 × HS-578T** (m = 5.89) | **UNDER: Galunisertib × AN3-CA** (m = 0.57) |
|---|---|---|
| N₀ cells / well | 1,067 | 305 |
| n★ (cells/arm) | **184** (OVER, 5.8× spare) | **22,475** (UNDER, 74× short) |
| downsample N₀ → | 184 ( = n★) | 184 ( ≪ n★ ) |
| Euclidean centroid distance from vehicle: full → down | 5.890 → 5.894 ± 0.127 (**Δ = 0.004, 0.07%**) | 0.570 → 0.632 ± 0.052 (**Δ = 0.06, +11%**) |
| cosine vs nearest neighbour #1: full → down | DMSO_T0 × HEC-1-A +0.722 → +0.723 ± 0.010 (**Δ = 0.001**) | Lapatinib × AN3-CA +0.680 → +0.603 ± 0.061 (**Δ = 0.077**) |
| cosine vs neighbour #2 | DMSO_T0 × AN3-CA +0.657 → +0.658 ± 0.010 (Δ = 0.001) | Galunisertib × HEC-1-A +0.656 → +0.585 ± 0.061 (Δ = 0.071) |
| cosine vs neighbour #3 | DMSO_T0 × C-33 A +0.583 → +0.584 ± 0.010 (Δ = 0.001) | Fluorouracil × HEC-1-A +0.643 → +0.571 ± 0.061 (Δ = 0.072) |
| replicate std (stability) | ≈ **0.010** (cosine), 0.13 (distance) | ≈ **0.061** (cosine), 0.05 (distance) |

## Analysis

The OVER-sampled centroid is functionally saturated: a **5.8-fold** reduction (1,067 → 184 cells, to its
own quota) shifts the vehicle distance by **0.07%** and the cosine correlation with each of its three
nearest neighbours by **≤ 0.001**, with replicate scatter of only ≈0.01 — the downstream similarity graph
is, for practical purposes, unchanged. Subjecting the UNDER-sampled drug to the *identical* 184-cell
downsampling degrades it: the vehicle distance inflates by **+11%** (noise adds to the tiny signal) and the
neighbour cosines fall by **0.07–0.08**, i.e. the correlation is both biased and unstable, with replicate
scatter (≈0.06) several-fold larger than the OVER case. This is the empirical confirmation of the dry-lab
ledger: downsampling to n★ is a *provably safe* compression for over-acquired wells but a *destructive* one
for under-acquired wells — and the calculator distinguishes the two a priori from magnitude alone.

## Auditor notes (do not omit from the manuscript)

- **The OVER exemplar is a reference population (DMSO_T0), not a drug.** In EmeraldBay's shared-line
  subset *no real drug is cleanly over-sampled* at θ★ = 0.1 rad — drug magnitudes against the per-line-mean
  baseline are small. The invariance is a geometric property of any high-magnitude, well-sampled centroid.
- **The UNDER exemplar (Galunisertib, a TGF-βR1 inhibitor) is deeply under-sampled** — m = 0.57, n★ = 22,475
  against N₀ = 305 (74× short), exactly the case the calculator flags as un-resolvable.
- The under-group magnitude is basis-dependent (m = 0.57 in the frozen embedding); its degradation under
  downsampling is milder than for a still-weaker drug, but remains an order of magnitude worse than the
  over-sampled exemplar.
- Metrics are centroid quantities in EmeraldBay's own frozen HVG/PCA space; E-distance is out of scope (see
  Scope above). Source: `scripts/emeraldbay_recompute/pass6_invariance.py` → `outputs/emeraldbay_invariance.json`.
