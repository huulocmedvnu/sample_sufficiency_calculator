> ⚠️ **SUPERSEDED.** Numbers here predate the audit (findings #1-#9, gaps 1-3). The authoritative source is `docs/MANUSCRIPT_DEEPSEEK.md` + `src/engine.py`. See `docs/AUDIT_DENOMINATORS.md` and `docs/ARCHITECTURE.md`.

# Technical Blueprint & Data Supplement
### Anisotropic Sample-Sufficiency for Single-Cell Perturbation-Direction Screens
**Status:** frozen source-of-truth for manuscript drafting, reconciled to the dose-resolved from-raw
recompute (the Constants of Record in `docs/SUPPLEMENT.md`).
**Provenance:** repository `sample_sufficiency_calculator` (proof `docs/THEORY.md`, code `src/calculator.py`,
tests `tests/`), external calibrator `scripts/emeraldbay_recompute/`. All numbers below are
machine-generated and traced to `fixtures/`; do not paraphrase the figures.

---

## 0. MATHEMATICAL AUDITOR'S CORRECTIONS (read before drafting)

Six points in the original drafting brief were inaccurate relative to the executed code and results, and
have been corrected throughout this document. **Use the corrected values.**

1. **Regime-gating tolerance is θ★ = 0.20 rad, not 0.1 rad.** The gating experiment
   (`scripts/emeraldbay_recompute/`) was run at θ★ = 0.20 rad, and compliance was defined as a realized
   RMS angle at n★ of ≤ 1.05·θ★. The "0.1 rad" in the brief is incorrect.
2. **Three of the four held-out validation groups are `DMSO_T0` reference populations, not drug effects.**
   They are the time-zero (DMSO_T0) population's displacement from the per-line centroid baseline, for
   HS-578T (m = 5.89), HEC-1-A (m = 3.58), and BT-474 (m = 3.35). Only `Irinotecan × HS-578T` (m = 1.96)
   is an actual drug perturbation. The held-out test validates the *angular-error geometry*; it is
   agnostic to biological interpretation. Do **not** describe the DMSO_T0 groups as "strong drug effects".
3. **The σ² values are a like-for-like within-condition comparison.** EmeraldBay's
   σ² ≈ 0.938 is a *within-condition residual* variance (full atlas, 52 lines; 0.848 on the 5-shared-line
   slice); the Tahoe headline σ² = 0.9567 is the matching *within-condition residual* per-cell variance.
   The two agree to ~2% — a like-for-like cross-platform match. The Tahoe *marginal* per-cell variance
   σ² = 2.406 (which additionally includes between-condition and between-line structure) is retained only
   as a **conservative upper bound / cross-check**, never as the headline. All support the claim that σ² is
   a platform/pipeline-specific plug-in requiring re-estimation.
4. **The implemented tail bound is Laurent–Massart (2000), not generic Hanson–Wright.** Laurent–Massart is
   the specialization to a weighted sum of χ² (our exact case) and yields explicit constants (2, 2);
   Hanson–Wright is the general quadratic-form parent inequality. Cite Laurent–Massart for n★_δ.
5. **The validation prediction used the large-control-pool form with a finite-population correction:**
   θ²(n) = tr(PΣP)/m² · (1/n − 1/N), factor 1 (the per-line-mean baseline is a near-noiseless large
   pool)---**not** the equal-arms factor-2 planning formula. Both are correct; they describe different
   experimental designs. Keep them distinct.
6. **The resource-allocation figures (~156× multiplex; 99%/~100% compute reduction) come from the Tahoe
   calibration** (N₀ = 1,296 cells/well, θ = 0.1 rad), not EmeraldBay. Attribute them to Tahoe.

---

## 1. INTRODUCTION & PROBLEM FORMULATION

**The screening matrix.** A single-cell perturbation screen yields a sparse count matrix
**X ∈ ℝ^{N×G}** (N cells, G ≈ 2×10⁴ genes). After normalization (`normalize_total`→`log1p`), selection
of highly variable genes, and z-scaling, the data are projected into a shared **d-dimensional
principal-component space** (here d = 50). Each experimental condition (a *well*: a drug × cell-line ×
dose) is summarized by its **pseudobulk centroid** μ = (1/n)∑ᵢ xᵢ ∈ ℝ^d, the sample mean of its n cells
in PCA space.

**The perturbation vector.** For a treated condition and its reference, define
**v = μ_treated − μ_control ∈ ℝ^d**, with **magnitude m = ‖v‖** (a potency/effect-size scalar) and
**unit direction u = v/m ∈ 𝕊^{d−1}** (the mechanistic axis). In MoA and drug-similarity analyses the
scientific quantity is the *direction* u (cosine geometry); the magnitude m is a nuisance.

**The angular error as the resolution metric.** The true v is never observed---only the finite-sample
estimate v̂ = μ̂_treated − μ̂_control. The **angular error θ = arccos(û·u)** between the estimated and true
direction is the fundamental limit on declaring two perturbations co-directional (same MoA). "How many
cells are enough?" is therefore the question: *how large must n be so that the RMS angular error of v̂
falls below a tolerance θ★?* This document derives, verifies, and externally validates the closed-form
answer.

---

## 2. FORMAL MATHEMATICAL PROOF (the skeleton)

**Sampling law.** Treated and control cells are i.i.d. with per-cell covariances Σ_t, Σ_c. Because each
centroid is an average, its scatter shrinks with the cell count; by the multivariate CLT and independence,

> v̂ ~ 𝒩(v, S),  S = Σ_t/n_t + Σ_c/n_c.   (S = O(1/n); write e = v̂ − v ~ 𝒩(0, S).)

**Tangent-space projector.** Let **P = I − uuᵀ** (symmetric, idempotent, rank d−1), the orthogonal
projector onto the tangent space of the unit sphere at u. By construction **P u = 0**: it keeps only the
part of a vector perpendicular to u.

**First-order Delta method on g(x) = x/‖x‖.** The geometric idea is that noise along the signal only
lengthens v̂ and cannot rotate it. Formally, the Jacobian is Dg(x) = ‖x‖⁻¹(I − x̂x̂ᵀ); at x = v this is
**Dg(v) = (1/m) P**. The radial direction u lies in ker Dg(v), so to first order û − u ≈ (1/m) P e, and

> Cov(û − u) = (1/m²) P S P.

**Why P u = 0 makes variance along the signal inert.** Because the Jacobian annihilates the radial
component (Dg(v)·u = 0), only the perpendicular noise P e reaches the direction estimate. The
mean-squared angular error is therefore the trace of the tangent-space covariance:

> **E[θ²] = tr(P S P) / m² = (tr S − uᵀ S u) / m² + O(ρ⁻⁴),  ρ² = m²/(uᵀ S u).**

The term uᵀSu is the noise variance *along* the signal; it is subtracted out. Sample sufficiency depends
**only on the perpendicular noise** tr(PSP). (Isserlis' theorem gives the relative second-order
correction; the result holds in the high-SNR regime ρ ≫ 1.)

**Closed-form anisotropic cell quota.** Setting E[θ²] = θ★² and S = 2Σ/n for equal arms (Σ_t ≈ Σ_c ≈ Σ):

> **n★ = 2 · tr(P Σ P) / (m² · θ★²) = 2 · Σ_k (1 − u_k²) ℓ_k / (m² θ★²)**  (PCA basis, ℓ_k = per-PC variance).

Two special cases anchor the formula: a large or shared control pool (n_c → ∞) replaces the factor 2 with
1, and isotropic noise Σ = σ²I recovers n★ = 2σ²(d−1)/(m²θ★²). The effective noise dimension is
d_eff = (tr PSP)²/tr((PSP)²) ≤ d−1.

**Tail-controlled (risk-managed) quota — Laurent–Massart (2000).** The mean quota controls the average
error; to bound the *chance* of a large error, note that ‖Pe‖² = Σ_i ν_i z_i² is a weighted χ² (ν_i the
eigenvalues of PSP, z_i iid 𝒩(0,1)). Laurent–Massart's bound with **exact constants**,

> Pr( ‖Pe‖² ≥ tr(PSP) + 2‖PSP‖_F √L + 2‖PSP‖_op L ) ≤ e^{−L},  L = log(1/δ),

inverts (every term ∝ 1/n) to a confidence quota guaranteeing Pr(θ > θ★) ≤ δ:

> **n★_δ = (2/(m²θ★²)) [ tr(PΣP) + 2‖PΣP‖_F √(log 1/δ) + 2‖PΣP‖_op log(1/δ) ]
>        = n★ · (1 + O(√(log(1/δ)/d_eff)) ).**

This is a multiplicative safety margin, not a different scaling law. *(Hanson–Wright is the general
quadratic-form inequality from which the weighted-χ² specialization descends; we implement the tighter
Laurent–Massart form.)* One caveat: the bound rigorously controls ‖Pe‖², and reaches the angle only
through the first-order step θ ≈ ‖Pe‖/m; it is therefore an approximate (1−δ) guarantee on θ, exact as
ρ → ∞ (see `docs/THEORY.md` §4). In practice it is conservative---empirical coverage is 0.34% against a
nominal 10%.

---

## 3. INTERNAL & STANDALONE VERIFICATION SUMMARY (`tests/verify_theory.py`)

Three asserted layers, run in CI on Python 3.10/3.11/3.12.

**Layer 1 — Symbolic (SymPy).** g(x)=x/‖x‖ is auto-differentiated and compared to the analytical
(1/m)(I − uuᵀ). The symbolic residual **simplifies to the zero matrix**, and a numeric cross-check at a
concrete vector v = (3, −1, 2) gives **max |autodiff − analytical| = 2.8×10⁻¹⁷** (machine precision).
Conclusion: the Jacobian Dg(v) = (1/m)P is exact.

**Layer 2 — Empirical (Monte Carlo).** d = 50; n = 50,000 cells/arm; K = 50,000 i.i.d. noise draws
e ~ 𝒩(0, 2Σ/n); Σ anisotropic (eigenvalues ∝ 0.85^k in a random orthonormal basis, tr Σ = 50);
deterministic seed = 0. The exact geometric angle arccos(⟨v,v̂⟩/(‖v‖‖v̂‖)) is measured. Analytical
E[θ²] = 2.196×10⁻⁴ against Monte-Carlo 2.193×10⁻⁴ gives a **relative error of 0.13%** (Monte-Carlo
std-error 0.18%). Validity diagnostics: along-signal SNR ρ² = 3.82×10⁵; effective noise dof d_eff = 12.2;
RMS angle ≈ 0.85°.

**Layer 3 — Implementation cross-check.** For tolerance θ★ = √(E[θ²]), the shipped
`calculate_cell_quota_anisotropic(Σ, v, θ★)` returns n★ = 50,000 exactly (0.0000% error), closing the
loop proof ↔ code ↔ simulation.

*(A pytest suite of 12 unit/property tests additionally pins the scaling laws n★ ∝ 1/m², ∝ σ², ∝ 1/θ★²,
the anisotropic→isotropic reduction, the perpendicular-only property, the dual-allocation regimes,
input validation, and Laurent–Massart tail coverage — empirical Pr(θ>θ★) = 0.34% ≤ δ = 10%.)*

---

## 4. EXTERNAL VALIDATION ON AN INDEPENDENT ATLAS (tahoebio/EmeraldBay)

**Rationale.** The quota is a theorem; an independent atlas cannot re-prove it, but it can test whether
the CLT/Gaussian-centroid *assumptions* hold on real, out-of-distribution single cells and whether the
plug-in calibration transfers. EmeraldBay (1.8×10⁶ cells; **5-day** culture vs Tahoe's 24 h; pooled
MOSAIC) shares five representative cell lines and the token/expression schema with Tahoe-100M.

**Procedure (`scripts/emeraldbay_recompute/`).** All 116 shards (≈ 57.7 GB; 1,831,648 cells across 52
lines) were streamed. EmeraldBay's *own* HVG(2000)+PCA(50) embedding and its **within-condition** residual
per-PC variance were fit over the **full atlas**, giving **σ² ≈ 0.938** (mean over 50 PCs). Per-cell
coordinates were then retained for the **141,720 cells of the 5 shared lines** (5-line within-condition
σ² = 0.848), on which the held-out validation and downsample checks are run. EmeraldBay now fixes this HVG
selection to a **frozen 2000-HVG set** (`fixtures/emeraldbay_frozen_hvg.json`), raising basis
reproducibility cos² from 0.79 to 0.96 (Tahoe was already stable at 0.93→0.995 and is left per-fit). The baseline is the
**per-line mean** (no DMSO at a matched timepoint; a `DMSO_T0` population exists and forms the
high-magnitude reference groups). *Note (item 3): σ² ≈ 0.938 is a within-condition value and matches the
Tahoe within-condition headline σ² = 0.9567 to ~2% — a like-for-like cross-platform agreement (the
marginal Tahoe σ² = 2.406 is a conservative bound only); both establish σ² as a platform/pipeline-specific
plug-in requiring re-estimation.*

**Held-out angular-error subsampling curves.** For a target group of N cells, the truth direction is taken
from all N; for a grid of n, cells are subsampled without replacement (R = 300 reps) and the realized RMS
angle to truth is measured. The prediction is the **finite-population-corrected, large-pool** form
θ²(n) = tr(PΣP)/m² · (1/n − 1/N) (item 5). A linear regression of realized θ² on (1/n − 1/N) tests
linearity-through-origin, with slope equal to tr(PΣP)/m².

| Validation group (actual identity) | N | m | mean rel. err | fitted slope | expected slope tr(PΣP)/m² | intercept | R² |
|---|---:|---:|---:|---:|---:|---:|---:|
| DMSO_T0 × HS-578T (reference population) | 1067 | 5.89 | 0.8% | 0.9060 | 0.9180 | 0.0 | 0.9999 |
| DMSO_T0 × HEC-1-A (reference population) | 435 | 3.58 | 0.6% | 1.8310 | 1.8450 | 0.0 | 0.9999 |
| DMSO_T0 × BT-474 (reference population) | 461 | 3.35 | 0.8% | 1.7090 | 1.7800 | 0.0 | 0.9977 |
| Irinotecan × HS-578T (drug; low SNR) | 563 | 1.96 | 5.3% | 10.0843 | 12.8079 | 0.0 | 0.9786 |

The realized curve matches the closed form to **0.6–5.3% mean relative error**, with
realized θ² ∝ (1/n − 1/N) **linear-through-origin (R² > 0.978, intercept ≈ 0)** in all four cases. For the
strongest group (HS-578T, m = 5.89) the fitted slope matches the anisotropic tr(PΣP)/m² to **1.4%**, and
the two moderate DMSO_T0 groups to within 1–4%. The genuine low-magnitude drug group (Irinotecan,
m = 1.96) shows a **~21% slope deficit concentrated at small n**---the **expected second-order Taylor /
finite-SNR breakdown** of the first-order Delta method (THEORY §5), since validity requires
ρ² = m²/(uᵀSu) ≫ 1.

**Regime-gating accuracy.** With θ★ = **0.20 rad** (item 1), each of **3,971 (condition × line) groups
across all 52 cell lines** (≥ 100 cells) was classified OVER- vs UNDER-sampled (OVER ⇔ N₀ ≥ n★): **347
(8.7%) are predicted OVER-sampled**. The downsample-and-measure check requires per-cell coordinates,
retained for the 5 shared lines: of the **10 predicted-OVER groups in that subset, 100%** achieved a
realized RMS angle ≤ 1.05·θ★ when downsampled to n★. The gating prediction is therefore empirically
corroborated on independent data.

**Distillation / CI.** Results are stored in `fixtures/emeraldbay_calibration.json` and asserted by
`tests/test_emeraldbay_integration.py` (network-free); figure `outputs/emeraldbay_heldout_validation.png`.

**Experimental unit & replication structure (pooled-condition re-gating).** n★ is a quota *per unique
perturbation condition* (the treated arm---the total cells aggregated into the condition's centroid),
**not** per physical well (formal treatment, including the batch-variance decomposition, in
`docs/SCALE_AUDIT.md` §4). A condition could in principle reach n★ by pooling cells across replicate
wells, which raises the question of whether single-well under-sampling is rescued once replicates are
pooled. We tested this directly on both reference atlases, and **in these data it is not**, because the
conditions are essentially **unreplicated**: in Tahoe-100M a given (drug × dose) occupies **a single plate
for 86.3%** of combinations (two for 13.4%, ≥ 3 for 0.4%), and EmeraldBay contributes a single pooled
sample per (perturbagen × dose) condition. Same-dose replicate plates are largely absent. Over the 56,827
(drug × dose × line) conditions (θ★ = 0.1 rad) only **10.6% are over-sampled** (89.4% under-sampled or
ghost). The median condition requires n★ ≈ 5,794 cells but holds ≈ 1,296 (an ~4.5× deficit that the sparse
replicate plates cannot close), and the line-resolved profiles are unchanged (high-magnitude agents
over-sampled in most lines; moderate and weak agents under-sampled in nearly all). **Conclusion:** the
pooled-condition unit is the mathematically correct one, but in these single-sample-per-condition atlases
it *coincides* with the single-condition unit, so the project-level conclusion stands---**≈ 10.6% of
conditions are safely saturated at a 5.7° tolerance**. Dose is a real lever (the over-sampled fraction
rises 8.8% → 14.2% across 0.05 → 5 µM) but does not move the bulk of the panel out of under-sampling. A
screen *designed* with replicate batches and within-batch (co-plated vehicle) referencing would benefit
from pooling---n★ being per-condition---but a single-sample-per-condition atlas does not, and cross-plate
replicate pooling without referencing is further capped by the batch-variance floor of `SCALE_AUDIT.md` §4.

---

## 5. TRANSLATIONAL & RESOURCE-ALLOCATION IMPACT

A single saturation threshold n★ governs two budgets.

**Wet-lab (sequencing economics).** Past n★, angular error improves only as 1/√n (diminishing returns),
so n★ is the knee of the curve. With a fixed per-lane read budget, a well acquired at N₀ > n★ wastes
(N₀ − n★) reads, and the surplus permits multiplexing N₀/n★ additional conditions (e.g. via cell hashing).
**On the Tahoe calibration (θ = 0.1 rad)** homoharringtonine at 5 µM in the responsive NCI-H460 line
(N₀ = 6,060) reaches n★ ≈ 39---a **~156× multiplexing gain**. Weak and moderate signatures (m ≲ 2) are
instead *under*-sampled (n★ > N₀); their cells are not redundant, and the calculator reports the regime
per drug rather than assuming redundancy.

**Dry-lab (compute economics).** For **centroid/pseudobulk** analyses, downsampling an over-acquired well
to n★ preserves its perturbation direction (and hence the drug-drug similarity graph) to θ★ by
construction. The compute saved is **polynomial in N**, not exponential: linear tasks (RAM, storage,
PCA-fit) shrink by 1 − n★/N₀; quadratic tasks (cell-cell pairwise, kernels, k-NN) by 1 − (n★/N₀)²;
neighbour graphs scale as ~ N log N. On the over-sampled cytotoxic example above this is
**99% (linear) / ~100% (quadratic)** reduction.

**Population-scale framing (correction of a common misconception).** The small n★ values of strong
cytotoxics are *not* representative. Across the 56,827-condition Tahoe panel (379 drugs × 3 doses × 50
lines) at θ★ = 0.1 rad the **median n★ is 5,794 cells/arm** (only 10.6% over-sampled; 28% need 10 k–50 k;
0.4% > 50 k). At a tight MoA tolerance the *typical* well is **under-sampled**---atlases are large by
aggregating many conditions, not because any single condition is cheap. n★ is strictly
**per-arm/per-condition**: the global atlas size enters no term of the variance S = Σ_t/n_t + Σ_c/n_c;
only the per-condition sum Σ_c n★_c sets the project budget (see `SCALE_AUDIT.md`). A 100-drug × 3-dose ×
2-line screen costs ≈ 10.9 M cells under flat 90th-percentile loading versus ≈ 3.5 M under
magnitude-adaptive allocation with a 10 k/well cap---a 3.1× reduction, with the ghost conditions flagged
rather than silently under-powered.

**Scope limitations (must be stated in the manuscript).** n★ certifies **centroid/direction** information
only. It does **not** guarantee preservation of **cell-resolution** structure---per-cell UMAP local
neighbourhoods, **rare-population / cell-type detection**, or trajectory branch points---which are
governed by local density, not by this angular threshold, and can be distorted by downsampling. Global PCA
eigenvalue spectra are preserved only *approximately* (random-matrix sampling error). Compute speed-ups
are polynomial; no operation here is exponential in N. σ² must be re-estimated per platform and pipeline
(the within-condition residual variance, not the marginal).

---

## Appendix — Provenance & reproducibility (for traceability, not prose)

- Theory: `docs/THEORY.md`. Estimator: `src/calculator.py`
  (`calculate_cell_quota_anisotropic`, `calculate_experimental_cell_quota`,
  `calculate_optimal_resource_allocation`). Verification: `tests/verify_theory.py`,
  `tests/test_calculator.py`, `tests/test_emeraldbay_integration.py`; CI `.github/workflows/ci.yml`.
- External calibrator: `scripts/emeraldbay_recompute/` → `fixtures/emeraldbay_calibration.json`.
- Constants of record: d = 50; Tahoe within-condition σ² = 0.9567 (headline; marginal 2.406 retained as a
  conservative bound), N₀ = 1,296; EmeraldBay within-condition
  σ² ≈ 0.938 (full atlas; 0.848 on the 5-shared-line slice), 141,720 cells streamed for the 5 shared
  lines, 3,971 gated (condition × line) groups across 52 lines, θ★ = 0.20 rad; Monte-Carlo K = 50,000,
  seed = 0, relative error 0.13%; symbolic residual 2.8×10⁻¹⁷.
- **Honesty ledger for reviewers:** (i) the headline compares like-for-like within-condition σ² across
  atlases (EmeraldBay 0.938 vs Tahoe 0.9567); the marginal Tahoe σ² = 2.406 is retained only as a
  conservative bound (item 3); (ii) three of the four validation groups are DMSO_T0 reference populations, not drug effects,
  and only Irinotecan × HS-578T is a drug (item 2); (iii) the low-SNR group exhibits the predicted
  first-order breakdown (not hidden); (iv) resource-allocation figures are Tahoe-derived (item 6); (v)
  cell-level topology preservation is explicitly out of scope.
