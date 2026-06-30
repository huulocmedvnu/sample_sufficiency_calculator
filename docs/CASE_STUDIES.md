# Empirical Case Studies — Sample-Sufficiency Spectrum (frozen source-of-truth)

**Configuration (Tahoe-100M calibration):** σ² = 7.66, d = 50, θ★ = 0.1 rad (5.73°), N₀ = 1,394
cells/well, equal arms. Governing law **n★ = 2(d−1)σ²/(m²θ★²) = 75,068 / m²**; over/under boundary at
**m = 7.34**; depth-fixed angular resolution **θ(N₀) = 0.734 / m rad**. Magnitudes m are raw
perturbation-vector norms (mean per-line ‖v‖) in the plate3-excluded PCA space — the space σ²=7.66 was
estimated in. All values machine-derived from the reference perturbation array via
`src/calculator.py::calculate_optimal_resource_allocation`.

| # | Drug (MoA / target) | m | n★ (cells/arm) | Regime | Wet-lab | Dry-lab (lin / quad) | Resolution @ N₀ | Action |
|---|---|---:|---:|---|---|---|---:|---|
| 1 | Homoharringtonine (protein-synthesis inhibitor) | 14.35 | 365 | OVER | 3.8× multiplex | 74% / 93% | 2.9° | Downsample ~3.8×; reallocate reads |
| 2 | Idarubicin (anthracycline; TOP2A) | 10.63 | 664 | OVER | 2.1× multiplex | 52% / 77% | 4.0° | Downsample ~2.1× |
| 3 | Dinaciclib (pan-CDK; CDK1/2/5/9) | 6.53 | 1,760 | UNDER (1.3×) | — | 0% / 0% | 6.4° | Tips in at θ=6.4° or +26% depth |
| 4 | Resveratrol (moderate-signal modulator) | 2.97 | 8,518 | UNDER (6.1×) | — | 0% / 0% | 14.2° | Downsampling forbidden; ~6× deeper required |
| 5 | Ribociclib (CDK4/6 inhibitor) | 0.84 | 107,547 | UNDER (77×) | — | 0% / 0% | 50.3° | Blind: MoA unresolvable at standard depth |

*Resveratrol is used purely as a moderate-magnitude exemplar (m ≈ 3); its metadata `moa-fine` is*
*"unclear" and no mechanistic claim is made.*

**Tolerance that tips an UNDER well into sufficiency at N₀=1,394:** Dinaciclib 6.4°, Resveratrol 14.2°,
Ribociclib 50.3° (θ_suff = θ★·√(n★/N₀)).

## Analysis

Across the spectrum the quota is governed entirely by squared transcriptional potency, n★ = 75,068/m²,
so a compound's perturbation magnitude alone fixes its sampling regime, with the over/under boundary
falling sharply at m ≈ 7.34. High-magnitude cytotoxics that displace cells far off baseline — the
protein-synthesis inhibitor Homoharringtonine (m = 14.4) and the anthracycline Idarubicin (m = 10.6) —
saturate their directional estimate in only 365–664 cells, leaving the standard 1,394-cell well 2–4×
over-provisioned and licensing 52–93% downstream compute reduction or 2–4× multiplexing via cell hashing.
As potency falls, n★ rises quadratically (a halving of m quadruples the requirement), and at fixed depth
the achievable angular resolution degrades linearly as θ = 0.734/m: the moderate modulator Resveratrol
(m = 3.0) is resolvable only to 14° at standard depth and strictly requires ~6× deeper sequencing —
downsampling is mathematically prohibited — while the CDK4/6 inhibitor Ribociclib (m = 0.84), despite
being a clinically pivotal targeted agent, carries a faint transcriptional footprint that renders it a
"ghost signature": its 1,394-cell direction estimate bears a ~50° error indistinguishable from a random
axis, and clearing it to 0.1 rad would demand ~108,000 cells per arm. Critically, magnitude is not
predicted by nominal mechanism — the pan-CDK/CDK9 inhibitor Dinaciclib (m = 6.5) and the CDK4/6 inhibitor
Ribociclib (m = 0.8) share a CDK-inhibitor annotation yet occupy opposite extremes of the spectrum —
establishing that sample sufficiency must be calibrated per-compound from its empirical perturbation
magnitude, never inferred from drug-class labels.

## Population context (do not over-read the strong tail)

These five profiles are **illustrative tiers spanning the spectrum, not a typical sample**. Across the
full 292-drug Tahoe panel at the same configuration, the **median n★ is 23,934 cells** (IQR
12,087–42,571); only **2%** of drugs are over-sampled (n★ < N₀ = 1,394), 64% require 10 k–50 k, and 17%
are ghosts (> 50 k). At a tight ~5.7° tolerance the *typical* perturbation is therefore **under-sampled**;
the over-sampled cases (Homoharringtonine, Idarubicin) are the m > 7.34 minority. See
[`SCALE_AUDIT.md`](SCALE_AUDIT.md) for the distribution and the project-budget ledger.

## Provenance / caveats (auditor)

- All m and n★ are machine-pulled from the reference perturbation array through the shipped calculator; no placeholders.
- **Paclitaxel is unavailable** in this array (plate3-excluded batch-clean build); the extreme-cytotoxic
  exemplar is **Homoharringtonine**, not Paclitaxel.
- m are **24 h survivor transcriptional** magnitudes (not viability), specific to this platform/timepoint/
  embedding.
- Resveratrol is a moderate-magnitude exemplar only (m ≈ 3); no mechanistic/MoA claim is attached to it.
  finding. Keep the distinction in the manuscript.
