# Empirical Case Studies — Sample-Sufficiency Spectrum (frozen source-of-truth)

**Configuration (Tahoe-100M, from-raw recompute; `scripts/tahoe_recompute/`).** σ² = 2.406, d = 50,
θ★ = 0.1 rad (5.73°), equal arms. The governing law is **n★ = 2(d−1)σ²/(m²θ★²) = 23,577 / m²**. At the
median depth N₀ = 1,296 the over/under boundary falls at **m = 4.27**, and the depth-fixed angular
resolution is **θ(N₀) = 0.427 / m rad**. The unit of analysis is the **(drug × dose × cell-line)
condition**: each magnitude m is the plate-matched (treated − DMSO) perturbation-vector norm ‖v‖ in the
shared PCA(50) space where σ² was estimated. All values are machine-derived; full tables are in
`fixtures/tahoe_*.csv`.

## Per-drug spectrum (median across 150 conditions = 3 doses × 50 lines)

| Drug (MoA / target) | median m | median n★ (cells/arm) | OVER | UNDER | Ghost |
|---|---:|---:|---:|---:|---:|
| Panobinostat (HDAC) | 4.24 | 1,314 | 73 | 77 | 0 |
| Homoharringtonine (protein synthesis) | 5.88 | 682 | 72 | 78 | 0 |
| Harringtonine (protein synthesis) | 3.67 | 1,753 | 64 | 86 | 0 |
| Idarubicin (anthracycline; TOP2A) | 2.93 | 2,745 | 45 | 105 | 0 |
| Trametinib (MEK) | 2.29 | 4,498 | 45 | 104 | 1 |
| Palbociclib (CDK4/6) | 1.41 | 11,778 | 0 | 150 | 0 |
| 4EGI-1 (eIF4E) | 1.34 | 13,202 | 0 | 149 | 1 |
| Crizotinib (ALK/MET) | 0.99 | 24,222 | 0 | 122 | 28 |

Even the strongest compounds are over-sampled in only about half of their 150 conditions, and no compound
is over-sampled in a majority of them. At the other extreme, **132 of 379 drugs (35%) are over-sampled in
no condition at all**. Full per-drug counts are in `fixtures/tahoe_per_drug.csv`.

## Two levers beyond drug identity — worked conditions

Magnitude, and therefore the quota, is set not by the compound alone but jointly by dose and cell-line
background. Two worked examples make this concrete.

**Dose** (homoharringtonine in the responsive NCI-H460 line): raising the concentration enlarges m and
lowers the quota.

| condition | m | n★ | N₀ | regime | dry save (lin / quad) | multiplex |
|---|---:|---:|---:|:--:|:--:|:--:|
| Homoharringtonine 0.05 µM × NCI-H460 | 10.3 | 222 | 2,904 | OVER | 92% / 99% | 13.1× |
| Homoharringtonine 0.5 µM × NCI-H460 | 14.5 | 112 | 2,114 | OVER | 95% / 100% | 18.9× |
| Homoharringtonine 5 µM × NCI-H460 | 15.5 | 98 | 6,060 | OVER | 98% / 100% | 61.8× |

**Cell line** (same drug and dose, different background): homoharringtonine at 5 µM reaches m = 15.5
(n★ = 98, OVER) in NCI-H460 but only m = 3.1 (n★ = 2,390, UNDER) in NCI-H661 — the same treatment flips
regime with the background.

## Analysis

The quota is governed by squared transcriptional potency, n★ = 23,577/m², so a condition's perturbation
magnitude alone fixes its sampling regime, with the over/under boundary at m ≈ 4.27 (at the median depth).
The consequences run in two directions.

At the strong end, high-magnitude cytotoxics displace cells far off baseline — the protein-synthesis
inhibitors (homoharringtonine, harringtonine) and the HDAC inhibitor panobinostat. These saturate their
directional estimate in a few hundred cells and are over-sampled in roughly half of their conditions,
which leaves those wells over-provisioned and licenses large downstream-compute reductions or multiplexing
via cell hashing.

At the weak end, n★ rises quadratically as potency falls — halving m quadruples the requirement — and at
fixed depth the achievable resolution degrades linearly as θ = 0.427/m. Weak targeted agents such as
palbociclib (CDK4/6, median m = 1.41) and crizotinib (ALK/MET, median m = 0.99) carry faint
transcriptional footprints: palbociclib is under-sampled in all 150 of its conditions, and crizotinib is a
"ghost" (n★ > 50,000) in 28.

The practical lesson is that magnitude is not predicted by nominal mechanism: targeted kinase inhibitors
sit at the weak end while broad cytotoxics dominate the strong end. Sample sufficiency must therefore be
calibrated per condition from the empirical perturbation magnitude — and jointly by dose and cell line —
never inferred from drug-class labels.

## Population context (do not over-read the strong tail)

The rows above are **illustrative tiers, not a typical sample**. Across the full 56,827-condition panel
(379 drugs × 3 doses × 50 lines) at the same configuration, the **median n★ is 14,570 cells** (IQR
6,489–27,656); only **2.5%** of conditions are over-sampled (n★ < N₀), 55% require 10 k–50 k, and 8.2% are
ghosts (> 50 k) — **97.5% under-sampled or worse**. At a tight ~5.7° tolerance the *typical* condition is
therefore **under-sampled**; the over-sampled cases are the high-magnitude / high-dose / responsive-line
minority. See [`SCALE_AUDIT.md`](SCALE_AUDIT.md) for the full distribution and the project budget.

## Provenance / caveats

- All m and n★ are machine-computed from the raw Tahoe-100M counts through the streaming pipeline
  (`scripts/tahoe_recompute/`) and the shipped calculator; no placeholders.
- Magnitudes are 24 h survivor transcriptional norms (not viability), specific to this
  platform/timepoint/embedding; the plate-matched DMSO_TF vehicle is the control.
- The condition unit is (drug × dose × cell line); pooling doses inflates N₀ and averages magnitudes, and
  is avoided here.
