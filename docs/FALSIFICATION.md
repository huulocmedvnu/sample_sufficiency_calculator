# Full-atlas falsification of the angular-error law

A "valid theory makes a prediction that can be checked against data." The cell-quota theorem itself is
proved (Lean 4 core + Monte-Carlo), so data cannot refute the *algebra*; what data can refute is the
**modelling assumption** — that real single-cell centroid noise behaves closely enough to the
anisotropic-Gaussian model that the closed form holds *quantitatively* on cells nobody chose to fit. This
document reports that test, run at full scale on both reference atlases.

The sharp, **parameter-free** prediction is that when \(n\) cells are subsampled (without replacement)
from a condition's acquired pool of \(N\) and their centroid direction is measured against the full-\(N\)
reference direction,

\[
\mathbb{E}[\theta^2(n)] \;=\; \frac{\operatorname{tr}(P\Sigma P)}{m^2}\Big(\tfrac1n-\tfrac1N\Big),
\]

a straight line through the origin in \((1/n-1/N)\) whose **slope \(\operatorname{tr}(P\Sigma P)/m^2\) is
computed from \(\Sigma\) and \(m\) independently of the subsampling curve** (no fitted parameter). Protocol
is `scripts/*/pass4*.py` and the manuscript Methods (all-\(N\) reference, per-line-mean baseline,
without-replacement resampling, geometric grid capped at \(N/2\)). Because the law is a *first-order*
(small-angle) result valid for \(\rho^2 = m^2/(u^\top\Sigma u)\gg 1\), we stratify by the theory's own
validity criterion (\(\rho^2\ge 3\)) rather than lumping in- and out-of-regime conditions together.

## Coverage (what actually ran)

| atlas | cells projected | conditions | tested (≥ cells) | m range |
|:--|--:|--:|--:|:--|
| **Tahoe-100M** (all 3,388 shards, re-streamed from raw) | 95,624,334 | 67,018 (sample × line) | 56,195 (≥ 400) | 0.25 – 13.0 |
| **EmeraldBay** (all 116 shards, 52 lines) | 1,831,648 | 4,912 (condition × line) | 1,064 (≥ 300) | 0.17 – 9.8 |

Every cell was projected through the *same* calibrated PCA(50) basis used for the manuscript's constants;
no basis was refit. Baseline is the per-line mean (near-noiseless large pool), chosen so the test is
baseline-agnostic and identical across atlases — distinct from the plate-matched DMSO reference used for
the magnitude spectrum in the main text.

## Result 1 — the σ² calibration holds across every cell

The manuscript's noise constants were fit on small subsamples (Tahoe: a 14-shard, ~0.4% subset). Projecting
*all* cells through the same basis confirms they are representative. The headline quota uses the Tahoe
**within-condition** residual variance σ² = 0.9567; the **marginal** per-cell variance σ² = 2.406 is
retained only as a conservative bound. The full-atlas re-stream re-measures the *marginal* quantity, so it
validates that conservative bound directly:

| atlas | quantity | full-atlas value | calibration | Δ |
|:--|:--|--:|--:|--:|
| Tahoe-100M | marginal σ² — conservative bound (mean over 50 PCs, all 95.6 M cells) | **2.4158** | 2.406 | +0.4 % |
| EmeraldBay | within-condition σ² (all 1.83 M cells, 52 lines) | **0.9745** | 0.963 | +1.2 % |

Both land within rounding of the calibrated values; the calibration fit on a fraction of a percent of the
data reproduces on the whole atlas. Note that EmeraldBay's full-atlas within-condition σ² ≈ 0.9745 is a
like-for-like match to the Tahoe within-condition headline σ² = 0.9567 — the two platforms agree on the
quantity that actually sets the quota.

## Result 2 — the parameter-free law is confirmed in-regime, on both atlases

Regressing the realized \(\theta^2\) on \((1/n-1/N)\) and comparing the fitted slope to the a-priori
\(\operatorname{tr}(P\Sigma P)/m^2\), for every qualifying condition:

| atlas | in-regime conditions (ρ²≥3) | median realized/predicted slope | median R² | whole-population median ratio | log–log Pearson |
|:--|--:|--:|--:|--:|--:|
| **Tahoe-100M** | 1,790 | **0.942** | **0.9987** | 0.679 | 0.974 |
| **EmeraldBay** | 33 | **0.981** | **0.9994** | 0.368 | 0.944 |

In the regime where the first-order law claims to hold, the *predicted* slope matches the *realized* slope
to within ~6 % (Tahoe, across 1,790 independent conditions) and ~2 % (EmeraldBay, 33 groups), each with
R² ≈ 0.999 — with **no free parameter**. This is the first time the direct subsample-and-measure test has
been run on Tahoe-100M itself (the primary atlas that anchors every headline constant); previously only
the *derived* quota spectrum had been applied to it.

The much lower whole-population ratios (0.68, 0.37) are **not** a failure: they are dominated by the many
real conditions that sit *below* the SNR regime, where the first-order law is known to break down (next).

## Result 3 — the predicted low-SNR breakdown is reproduced quantitatively

The theory predicts not only where the law holds but where it fails: the relative bias grows as
\(\rho^{-2}\). Across the full magnitude range the realized/predicted slope ratio falls monotonically as
the along-signal SNR drops:

| atlas | Spearman(1 − ratio, 1/ρ²) | median ratio, high SNR | median ratio, low SNR |
|:--|--:|--:|--:|
| **Tahoe-100M** | 0.567 | — | — |
| **EmeraldBay** | 0.789 | 0.532 | 0.241 |

Both correlations are positive: the deviation from the parameter-free prediction is governed by \(\rho^2\),
exactly as the second-order Delta expansion predicts (THEORY §5). The theory therefore forecasts its own
failure mode, and the data obey the forecast.

## Result 4 — the variance structure (two independent halves)

A cleaner check that avoids the "truth from the same cells" overlap: split each EmeraldBay group in two,
estimate the direction from two *disjoint* \(n\)-cell subsamples, and compare the angle between them to the
angle against the full-\(N\) reference. Two independent estimates should differ by \(\sqrt2\) times as much
as one estimate differs from a near-noiseless truth. Over 642 groups the median ratio is **1.352**
(IQR 1.29–1.41) against the target \(\sqrt2 = 1.414\) — the variance-doubling holds, with the mild
shortfall again the low-SNR saturation of the bulk of the panel.

## Honest caveats

- **Baseline choice.** The reference is the per-line mean, not a plate-matched DMSO control; this validates
  the *sampling law for whatever direction is defined*, and keeps the two atlases comparable, but the
  m-values here are not the DMSO-referenced magnitudes of the manuscript's spectrum.
- **Sampling scheme.** Both atlases use without-replacement subsampling, matching the finite-population
  \((1/n-1/N)\) prediction (an earlier Tahoe draft used with-replacement and was corrected before this run).
- **In-regime population is small on EmeraldBay** (33 of 1,064 groups) and modest on Tahoe (1,790 of
  56,195): against a per-line-mean baseline, most real perturbation directions are weak, so the *majority*
  of conditions are out-of-regime — itself consistent with the manuscript's finding that most conditions
  are under-sampled or "ghosts" at a stringent tolerance.
- **What is and isn't tested.** This corroborates the CLT/Gaussian-centroid *assumptions* on real,
  out-of-distribution cells; it does not (and cannot) re-prove the theorem.

## Provenance

Scripts: `scripts/emeraldbay_recompute/pass2b_project_all.py` (full-atlas EB projection),
`scripts/emeraldbay_recompute/pass4_falsification.py` (EB population test + σ² check),
`scripts/tahoe_recompute/pass4a_stream.py` (full Tahoe re-stream + projection to memmap),
`scripts/tahoe_recompute/pass4b_curves.py` (per-condition curves + σ² check). Fixtures:
`fixtures/emeraldbay_falsification.json` (5-line slice), `fixtures/emeraldbay_falsification_full.json`
(full 52-line atlas), `fixtures/tahoe_direct_curves_full.json`. Figures:
`outputs/emeraldbay_falsification.png`, `outputs/tahoe_direct_curves_full.png`.
