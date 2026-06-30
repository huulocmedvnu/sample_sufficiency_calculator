# Technical Blueprint & Data Supplement
### Anisotropic Sample-Sufficiency for Single-Cell Perturbation-Direction Screens
**Status:** frozen source-of-truth for manuscript drafting (target: *Nature Methods* / *Bioinformatics*).
**Provenance:** repository `sample_sufficiency_calculator` (proof `docs/THEORY.md`, code `src/calculator.py`,
tests `tests/`), external calibrator `obgyn/scripts/calibrate_emeraldbay.py`. All numbers below are
machine-generated; do not paraphrase the figures.

---

## 0. MATHEMATICAL AUDITOR'S CORRECTIONS (read before drafting)

The following six items in the drafting brief were inaccurate relative to the executed code/results and
have been corrected throughout this document. **Use the corrected values.**

1. **Regime-gating tolerance is θ★ = 0.20 rad, not 0.1 rad.** The gating experiment
   (`calibrate_emeraldbay.py`) was run at θ★ = 0.20 rad; compliance was defined as realized RMS angle at
   n★ ≤ 1.15·θ★. The "0.1 rad" in the brief is incorrect.
2. **The two high-magnitude validation groups are `DMSO_T0 × HS-578T` and `DMSO_T0 × AN3-CA`, not "drug
   effects".** These are the time-zero (DMSO_T0) reference population's displacement from the per-line
   centroid baseline. Only `Encorafenib × HEC-1-A` is an actual drug perturbation. The held-out test
   validates the *angular-error geometry*; it is agnostic to biological interpretation. Do **not**
   describe HS-578T/AN3-CA as "strong drug effects".
3. **σ² = 2.1 (EmeraldBay) and σ² = 7.66 (Tahoe) are different variance definitions and are NOT a clean
   head-to-head.** EmeraldBay's 2.1 is the *within-condition residual* variance; the Tahoe calculator
   default 7.66 is a *marginal* per-cell variance (includes between-condition/between-line structure).
   Both support the claim "σ² is a platform/pipeline-specific plug-in requiring re-estimation," but a
   like-for-like within-condition Tahoe σ² was **not** computed. State this caveat; do not assert the
   gap is purely platform-driven.
4. **The tail bound implemented is Laurent–Massart (2000), not generic Hanson–Wright.** Laurent–Massart
   is the specialization to a weighted sum of χ² (our exact case) and yields explicit constants (2, 2);
   Hanson–Wright is the general quadratic-form parent inequality. Cite Laurent–Massart for n★_δ.
5. **The validation prediction used the large-control-pool form with a finite-population correction:**
   θ²(n) = tr(PΣP)/m² · (1/n − 1/N), factor 1 (the per-line-mean baseline is a near-noiseless large
   pool), **not** the equal-arms factor-2 planning formula. Both are correct; they describe different
   experimental designs. Keep them distinct.
6. **The resource-allocation figures (3.8× multiplex; 74%/93% compute reduction) are from the Tahoe
   calibration (N₀ = 1394 cells/well, θ = 0.1 rad), not EmeraldBay.** Attribute them to Tahoe.

---

## 1. INTRODUCTION & PROBLEM FORMULATION

**The screening matrix.** A single-cell perturbation screen yields a sparse count matrix
**X ∈ ℝ^{N×G}** (N cells, G ≈ 2×10⁴ genes). After normalization (`normalize_total`→`log1p`), selection
of highly variable genes, and z-scaling, the data are projected into a shared **d-dimensional principal-
component space** (here d = 50). Each experimental condition (a *well*: a drug × cell-line × dose) is
summarized by the **pseudobulk centroid** μ = (1/n)∑ᵢ xᵢ ∈ ℝ^d, the sample mean of its n cells in PCA
space.

**The perturbation vector.** For a treated condition and its reference, define
**v = μ_treated − μ_control ∈ ℝ^d**, with **magnitude m = ‖v‖** (a potency/effect-size scalar) and
**unit direction u = v/m ∈ 𝕊^{d−1}** (the mechanistic axis). In MoA / drug-similarity analyses the
scientific quantity is the *direction* u (cosine geometry); m is a nuisance.

**The angular error as the resolution metric.** v is never observed; only the finite-sample estimate
v̂ = μ̂_treated − μ̂_control. The **angular error θ = arccos(û·u)** between the estimated and true
direction is the fundamental resolution limit for declaring two perturbations co-directional (same MoA).
"How many cells are enough?" is therefore the question: *how large must n be so that the RMS angular
error of v̂ falls below a tolerance θ★?* This document derives, verifies, and externally validates the
closed-form answer.

---

## 2. FORMAL MATHEMATICAL PROOF (the skeleton)

**Sampling law.** Treated/control cells are i.i.d. with per-cell covariances Σ_t, Σ_c; by the
multivariate CLT and independence,

> v̂ ~ 𝒩(v, S),  S = Σ_t/n_t + Σ_c/n_c.   (S = O(1/n); write e = v̂ − v ~ 𝒩(0, S).)

**Tangent-space projector.** Let **P = I − uuᵀ** (symmetric, idempotent, rank d−1), the orthogonal
projector onto the tangent space of the unit sphere at u. By construction **P u = 0**.

**First-order Delta method on g(x) = x/‖x‖.** The Jacobian is
Dg(x) = ‖x‖⁻¹(I − x̂x̂ᵀ); evaluated at x = v this is **Dg(v) = (1/m) P**. The radial direction u lies in
ker Dg(v): perturbations of v̂ *along* u do not rotate the direction. Hence to first order
û − u ≈ (1/m) P e, and

> Cov(û − u) = (1/m²) P S P.

**Why P u = 0 makes variance along the signal inert.** Because the Jacobian annihilates the radial
component (Dg(v)·u = 0), only the perpendicular noise P e contributes to the direction error. The
mean-squared angular error is the trace of the tangent-space covariance:

> **E[θ²] = tr(P S P) / m² = (tr S − uᵀ S u) / m² + O(ρ⁻⁴),  ρ² = m²/(uᵀ S u).**

The term uᵀSu is the noise variance *along* the signal; it is subtracted out. Sample sufficiency depends
**only on the perpendicular noise** tr(PSP). (Isserlis' theorem gives the relative second-order
correction O(ρ⁻²); the result holds in the high-SNR regime ρ ≫ 1.)

**Closed-form anisotropic cell quota.** Setting E[θ²] = θ★² and S = 2Σ/n for equal arms (Σ_t ≈ Σ_c ≈ Σ):

> **n★ = 2 · tr(P Σ P) / (m² · θ★²) = 2 · Σ_k (1 − u_k²) ℓ_k / (m² θ★²)**  (PCA basis, ℓ_k = per-PC variance).

Special cases: a large/shared control pool (n_c → ∞) replaces the factor 2 with 1; Σ = σ²I recovers the
isotropic n★ = 2σ²(d−1)/(m²θ★²). Effective noise dimension d_eff = (tr PSP)²/tr((PSP)²) ≤ d−1.

**Tail-controlled (risk-managed) quota — Laurent–Massart (2000).** ‖Pe‖² = Σ_i ν_i z_i² is a weighted
χ² (ν_i = eigenvalues of PSP, z_i iid 𝒩(0,1)). Laurent–Massart's bound with **exact constants**,

> Pr( ‖Pe‖² ≥ tr(PSP) + 2‖PSP‖_F √L + 2‖PSP‖_op L ) ≤ e^{−L},  L = log(1/δ),

inverts (every term ∝ 1/n) to a confidence quota guaranteeing Pr(θ > θ★) ≤ δ:

> **n★_δ = (2/(m²θ★²)) [ tr(PΣP) + 2‖PΣP‖_F √(log 1/δ) + 2‖PΣP‖_op log(1/δ) ]
>        = n★ · (1 + O(√(log(1/δ)/d_eff)) ).**

This is a multiplicative safety margin, not a different scaling law. *(Hanson–Wright is the general
quadratic-form inequality from which the weighted-χ² specialization descends; we implement the tighter
Laurent–Massart form.)*

---

## 3. INTERNAL & STANDALONE VERIFICATION SUMMARY (`tests/verify_theory.py`)

Three asserted layers; run in CI on Python 3.10/3.11/3.12.

**Layer 1 — Symbolic (SymPy).** g(x)=x/‖x‖ is auto-differentiated and compared to the analytical
(1/m)(I − uuᵀ). The symbolic residual **simplifies to the zero matrix**; numeric cross-check at a
concrete vector v = (3, −1, 2) gives **max |autodiff − analytical| = 2.8×10⁻¹⁷** (machine precision).
Conclusion: the Jacobian Dg(v) = (1/m)P is exact.

**Layer 2 — Empirical (Monte Carlo).** d = 50; n = 50,000 cells/arm; K = 50,000 i.i.d. noise draws
e ~ 𝒩(0, 2Σ/n); Σ anisotropic (eigenvalues ∝ 0.85^k in a random orthonormal basis, tr Σ = 50);
deterministic seed = 0. The exact geometric angle arccos(⟨v,v̂⟩/(‖v‖‖v̂‖)) is measured.
Analytical E[θ²] = 2.196×10⁻⁴ vs Monte-Carlo 2.193×10⁻⁴ → **relative error 0.13%** (Monte-Carlo
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

**Rationale.** The quota is a theorem; an independent atlas cannot re-prove it but can test whether the
CLT/Gaussian-centroid *assumptions* hold on real, out-of-distribution single cells and whether the
plug-in calibration transfers. EmeraldBay (1.8×10⁶ cells; **5-day** culture vs Tahoe's 24 h; pooled
MOSAIC) shares the 5 gyn/breast lines and the token/expression schema.

**Procedure (`calibrate_emeraldbay.py`).** All 116 shards (≈58 GB) were streamed and discarded, retaining
**142,883 cells** of the 5 shared lines; EmeraldBay's *own* HVG(2000)+PCA(50) embedding was fit. The
baseline is the **per-line mean** (no DMSO at matched timepoint; a `DMSO_T0` population exists and forms
a high-magnitude reference group). The **within-condition** residual per-PC variance was **σ² ≈ 2.12**
(mean over 50 PCs). *Auditor caveat (item 3): this within-condition value is not directly comparable to
the Tahoe calculator default 7.66, which is marginal; both nonetheless establish σ² as a
platform/pipeline-specific plug-in requiring re-estimation.*

**Held-out angular-error subsampling curves.** For a target group of N cells, the truth direction is
taken from all N; for a grid of n, cells are subsampled without replacement (R = 300 reps) and the
realized RMS angle to truth is measured. The prediction is the **finite-population-corrected, large-pool**
form θ²(n) = tr(PΣP)/m² · (1/n − 1/N) (item 5). A linear regression of realized θ² on (1/n − 1/N) tests
linearity-through-origin (slope should equal tr(PΣP)/m²).

| Validation group (actual identity) | N | m | mean rel. err | fitted slope | expected slope tr(PΣP)/m² | intercept | R² |
|---|---:|---:|---:|---:|---:|---:|---:|
| DMSO_T0 × HS-578T (high magnitude) | 1067 | 12.41 | 1.4% | 0.7692 | 0.7628 | −4.8×10⁻⁵ | 0.9979 |
| DMSO_T0 × AN3-CA (moderate-high) | 752 | 6.39 | 0.8% | 1.3791 | 1.3753 | −2.6×10⁻⁴ | 0.9996 |
| Encorafenib × HEC-1-A (drug; low SNR) | 1117 | 3.02 | 2.4% | 6.6596 | 7.8217 | 6.8×10⁻³ | 0.9954 |

The realized curve matches the closed form to **0.8–2.4% mean relative error**, with
realized² ∝ (1/n − 1/N) **linear-through-origin (R² > 0.99, intercept ≈ 0)** in all three cases. For the
two high-magnitude groups the fitted slope matches the anisotropic tr(PΣP)/m² to **< 1%**; the genuine
low-magnitude drug group (Encorafenib, m = 3.0) shows a **~15% slope deficit concentrated at small n**
(8.3% point error at n = 20), the **expected second-order Taylor / finite-SNR breakdown** of the
first-order Delta method (THEORY §5), since validity requires ρ² = m²/(uᵀSu) ≫ 1.

**Regime-gating accuracy.** With θ★ = **0.20 rad** (item 1), each of 101 wells was classified OVER- vs
UNDER-sampled (OVER ⇔ N₀ ≥ n★). Of the **15 wells predicted OVER-sampled, 100%** achieved realized RMS
angle ≤ 1.15·θ★ when downsampled to n★ (downsample-and-measure). The gating prediction is therefore
empirically corroborated on independent data.

**Distillation / CI.** Results are stored in `fixtures/emeraldbay_calibration.json` and asserted by
`tests/test_emeraldbay_integration.py` (network-free); figure `outputs/emeraldbay_heldout_validation.png`.

---

## 5. TRANSLATIONAL & RESOURCE-ALLOCATION IMPACT

A single saturation threshold n★ governs two budgets.

**Wet-lab (sequencing economics).** Past n★, angular error improves only as 1/√n (diminishing returns);
n★ is the knee. With a fixed per-lane read budget, a well acquired at N₀ > n★ wastes (N₀ − n★) reads;
the surplus permits multiplexing N₀/n★ additional conditions (e.g., via cell hashing). **On the Tahoe
calibration (N₀ = 1394 cells/well, θ = 0.1 rad)** a maximal-magnitude cytotoxic (m ≈ 14.3) reaches
n★ ≈ 365, a **3.8× multiplexing gain**; weak/moderate signatures (m ≈ 1.3–3) are *under*-sampled
(n★ > N₀), i.e. their cells are not redundant — the calculator reports the regime per drug rather than
assuming redundancy.

**Dry-lab (compute economics).** For **centroid/pseudobulk** analyses, downsampling an over-acquired well
to n★ preserves its perturbation direction (and hence the drug-drug similarity graph) to θ★ by
construction. The compute saved is **polynomial in N**, not exponential:
linear (RAM, storage, PCA-fit): 1 − n★/N₀; quadratic (cell-cell pairwise / kernels / k-NN):
1 − (n★/N₀)²; near-linear (neighbour graphs ~ N log N). On the Tahoe over-sampled cytotoxic example this
is **74% (linear) / 93% (quadratic)** reduction.

**Population-scale framing (correction of a common misconception).** The small n★ values of strong
cytotoxics are *not* representative: across the 292-drug Tahoe panel at θ★ = 0.1 rad the **median n★ is
23,934 cells/arm** (only 2% over-sampled; 64% need 10 k–50 k; 17% > 50 k). At a tight MoA tolerance the
*typical* well is **under-sampled** — atlases are large by aggregating many conditions, not because any
single condition is cheap. n★ is strictly **per-arm/per-condition**: the global atlas size enters no term
of the variance S = Σ_t/n_t + Σ_c/n_c, only the per-well sum Σ_w n★_w sets project budget (see
`SCALE_AUDIT.md`). A 100-drug × 3-dose × 2-line screen costs ≈ 39 M cells under flat loading vs ≈ 5.7 M
under magnitude-adaptive allocation with a 10 k/well cap (85/100 drugs then capped below 0.1 rad).

**Scope limitations (must be stated in the manuscript).** n★ certifies **centroid/direction** information
only. It does **not** guarantee preservation of **cell-resolution** structure — per-cell UMAP local
neighbourhoods, **rare-population / cell-type detection**, or trajectory branch points — which are
governed by local density, not by this angular threshold, and can be distorted by downsampling. Global
PCA eigenvalue spectra are preserved only *approximately* (random-matrix sampling error). Compute
speed-ups are polynomial; no operation here is exponential in N. σ² must be re-estimated per platform/
pipeline (the within-condition residual variance, not the marginal).

---

## Appendix — Provenance & reproducibility (for traceability, not prose)

- Theory: `docs/THEORY.md`. Estimator: `src/calculator.py`
  (`calculate_cell_quota_anisotropic`, `calculate_experimental_cell_quota`,
  `calculate_optimal_resource_allocation`). Verification: `tests/verify_theory.py`,
  `tests/test_calculator.py`, `tests/test_emeraldbay_integration.py`; CI `.github/workflows/ci.yml`.
- External calibrator: `obgyn/scripts/calibrate_emeraldbay.py` → `fixtures/emeraldbay_calibration.json`.
- Constants of record: d = 50; Tahoe marginal σ² = 7.66, N₀ = 1394; EmeraldBay within-condition
  σ² ≈ 2.12, 142,883 gyn cells, 101 gated wells, θ★ = 0.20 rad; Monte-Carlo K = 50,000, seed = 0,
  relative error 0.13%; symbolic residual 2.8×10⁻¹⁷.
- **Honesty ledger for reviewers:** (i) σ² values use different variance definitions across atlases
  (item 3); (ii) the two high-m validation groups are DMSO_T0 reference populations, not drug effects
  (item 2); (iii) the low-SNR group exhibits the predicted first-order breakdown (not hidden); (iv)
  resource-allocation figures are Tahoe-derived (item 6); (v) cell-level topology preservation is
  explicitly out of scope.
