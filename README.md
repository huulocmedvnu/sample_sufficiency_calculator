> # ⚠️ Single source of truth
> **Every number in the study is produced by `src/engine.py`** (via `scripts/run_unified_spectrum.py` →
> `fixtures/unified_spectrum.json`). Any number not produced by it is stale. The authoritative narrative
> is `docs/MANUSCRIPT_DEEPSEEK.md`; the current constants of record are `docs/SUPPLEMENT.md`.
>
> **A retired thesis was removed by the control-pool audit.** The numbers `n★ = 9,376/m²`,
> `10.6% over`, `89% under`, `0.4% ghost`, the equal-arm `2·tr(PΣP)/(m²θ²)` used as a headline, and
> the marginal-σ² spectrum (`2.5/89.3/8.2`, `23,577`) all belong to that retired thesis and appear now
> **only** inside `docs/AUDIT_LOG.md` / `docs/AUDIT_DENOMINATORS.md` / `docs/ARCHITECTURE.md`. If you
> see them anywhere else, treat them as stale — do not quote them.
>
> **This README's body predates the audit** and is kept for orientation only; trust the manuscript,
> `src/engine.py`, and `docs/SUPPLEMENT.md` over anything below. Retired-artifact list: `docs/ARCHITECTURE.md`.

# Dual-Sided Sample-Sufficiency Calculator for Single-Cell Perturbation Screens

[![CI](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml/badge.svg)](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml)

A small, dependency-light calculator built on **one information-saturation threshold `n*`** (the cells
needed to resolve a drug's perturbation *direction* to an angular tolerance), read across **two
ledgers**:

> **Wet-lab:** how many cells per arm must I sequence to pin a drug's direction to within tolerance —
> and if I'm over that, how much can I multiplex instead?
>
> **Dry-lab:** when is it provably safe to downsample over-acquired wells before heavy matrix
> operations, and how much RAM/CPU does that actually save?

Strong perturbations need quadratically fewer cells than weak ones (`n* ∝ 1/m²`), so `n*` is the lever
for both sequencing budget and compute budget. The dry-lab guarantee is scoped to **centroid /
pseudobulk** analyses (the drug-similarity graph) and the speedups are **polynomial, not exponential**
— see the honest scope notes below.

---

## The math (Delta method on the angular error)

A drug's perturbation vector in a `d`-dimensional embedding is

```
v = mu_treated - mu_control
```

where each centroid `mu` is a **mean over single cells**. By the CLT a centroid from `n` cells has
covariance `Sigma_cell / n`, so for `n` cells per arm

```
Cov(v_hat) ~ Sigma_t/n + Sigma_c/n ~ (2/n) * Sigma_cell .
```

We only care about the **direction** `u = v / ||v||`. Decompose the estimation noise into components
parallel and perpendicular to `v`; only the perpendicular part rotates the direction. With the
tangent-space projector `P = I - u u^T` (rank `d-1`),

```
E[ || P (v_hat - v) ||^2 ] = tr( P Cov(v_hat) P )  ~  (2 sigma^2 / n) (d - 1)
```

under the isotropic approximation `Sigma_cell ~ sigma^2 I` (`sigma^2` = mean per-cell variance per
dim). For small angles, `theta ~ ||perp noise|| / ||v||`, giving the **RMS angular error**

```
E[theta^2]  ~  2 (d - 1) sigma^2 / ( n * m^2 ) ,      m = ||v|| .
```

### The angular error ball

Geometrically, the estimated direction lands inside a `(d-1)`-dimensional spherical cap — the
**angular error ball** — centred on the true direction with RMS radius

```
theta_RMS(n) = sqrt( 2 (d-1) sigma^2 / ( n m^2 ) )   [radians].
```

More cells shrink the ball as `1/sqrt(n)`; a bigger perturbation `m` shrinks it as `1/m`.

### Solving for the quota

Set `theta_RMS = tolerance` and invert:

```
                 2 (d - 1) sigma^2
   n*  =  ---------------------------------       [cells PER ARM]
              m^2  *  tolerance^2
```

This is `calculate_experimental_cell_quota(single_cell_variance, num_dimensions,
perturbation_magnitude, tolerance)` in `src/calculator.py`. The `1/m^2` factor is the key lever:
**halving the effect size quadruples the cells needed.**

**Approximation caveats (read before quoting a number):**
1. *Isotropic variance.* The exact form uses `tr(P Sigma P)` for the specific `v`; when `v` aligns
   with high-variance directions the true `n*` departs from the isotropic estimate (here per-PC
   variance ranges 0.91–14.0, marginal mean 2.406).
2. `tolerance` is an **RMS angular SD in radians**. Sub-degree tolerances on noisy single-cell data
   demand very large `n` — e.g. `tolerance=0.01` (~0.57°) needs ~10^5–10^6 cells/arm. Realistic
   values are typically **0.05–0.2 rad (≈3–11°)**.
3. `sigma^2` is tied to a fixed pipeline (HVG selection, normalization, PCA `d`). Re-derive it if you
   change the embedding.

---

## Calibration on Tahoe-100M (from the raw 100.6M-cell dataset)

`src/calibrate.py` reads the **fresh recompute** (`scripts/tahoe_recompute/`, following the theislab
[vevo_100m recipe](https://theislab.github.io/vevo_Tahoe_100m_analysis/vevo_100m_pca.html):
`normalize_total(1e4)` → `log1p` → HVG(2000) → PCA(50), streamed over all 100,648,790 cells). The
theory-preferred **within-condition** per-cell PCA variance — the residual noise *within* a fixed
condition, which is what actually limits direction estimation — is **`sigma^2 = 0.9567`** (mean over
50 dims), giving the headline quota law **`n* = 9,376 / m²`** at θ=0.1 rad. (The **marginal** per-cell
variance **`sigma^2 = 2.406`**, which mixes all conditions together, is retained only as a
**conservative bound**: `n* = 23,577 / m²`.) A strong cytotoxic — homoharringtonine at 5 µM in the
responsive NCI-H460 line (m=15.5) — resolves its direction in **~39 cells/arm**, whereas a
median-magnitude condition (m≈1.27) needs **~5,794**, and the requirement scales as `1/θ²` with the
tolerance. The **quota ratio between any two drugs equals `(m_a/m_b)²`** exactly — the `m²` law, which
`calibrate.py` verifies.

**The panel (379 drugs × 3 doses × 50 lines = 56,827 conditions, θ=0.1 rad).** The condition is
`(drug × dose × cell line)`; median `N0 = 1,296` cells/condition. At the within-condition headline,
**10.6%** of conditions are over-sampled, **89.0%** are under-sampled, and **0.4%** are "ghosts"
(`n* > 50,000`) — **89.4% under-sampled or worse**; **28 of 379 drugs are resolvable in no condition**
at standard depth. Full
per-condition / per-drug / per-dose / per-line tables are committed under `fixtures/tahoe_*.csv`.

> **Provenance.** All numbers are recomputed from the raw Tahoe-100M counts (`scripts/tahoe_recompute/`,
> a resumable streaming pipeline over the 337 GB expression data); the study-design layout
> (1,344 wells = 14 plates × 96; 5.0% QC-filtered) is derived from `obs_metadata`. See
> `docs/SUPPLEMENT.md` for the consolidated constants of record.

---

## Anisotropic form (peer-review standard) — `docs/THEORY.md`

The isotropic `σ²I` above is a placeholder; single-cell noise is strongly anisotropic. The rigorous
result (full derivation, distribution, finite-sample tail bound, and validity regime in
[`docs/THEORY.md`](docs/THEORY.md)) replaces `(d−1)σ²` with the **signal-orthogonal noise trace**:

```
            2 · tr(P Σ P)            2 · Σ_k (1 − u_k²) ℓ_k
n*_aniso = ────────────────  =  ──────────────────────────   (P = I − uuᵀ, u = v/‖v‖)
              m² · θ²                     m² · θ²                ℓ_k = per-PC variance
```

Variance **along** the signal is subtracted out (`u_k²ℓ_k`); only perpendicular noise rotates the
direction. This **contains the isotropic formula** as the case `Σ = σ²I` (verified in `calibrate.py`).
`calculate_cell_quota_anisotropic(covariance, perturbation_vector, tolerance, ...)` also returns the
**effective noise dimension** `d_eff = tr(PΣP)²/tr((PΣP)²)` and, with `confidence=δ`, a Hanson–Wright
**tail-controlled quota** guaranteeing `P(θ>θ*) ≤ δ` (not just the mean).

**Honest empirical finding (Tahoe-100M within-condition covariance):** on this atlas the anisotropic *mean-quota* correction is
**small — aniso/iso ratio ≈ 0.93–0.96** (drug directions carry little variance along themselves, so they
are not aligned with the dominant cell-cycle/lineage PCs, and with d=50 no single axis can move the trace
much); the full-atlas falsification (below) puts the in-regime realized/predicted slope ratio in the same
band (**0.94–0.98**). The anisotropic machinery's real value here is (i) the effective noise dimension
**`d_eff ≈ 29 ≪ 49`** and (ii) the rigorous **tail-controlled quota** — 95%-confident homoharringtonine
`n*=83` vs mean `37` cells/arm, a **~2.3× safety factor** (all from `python src/calibrate.py`). On data
where perturbations *do* align with high-variance axes the correction is large — the ratio formula in §8
of `THEORY.md` says exactly when.

## External validation on an independent atlas (EmeraldBay)

The formula is a theorem (verified symbolically + by Monte Carlo), so a second dataset can't make it
*more* true — but it **can** test whether the CLT/Gaussian-centroid *assumptions* hold on real,
independent single cells, and whether the calibration transfers. We did this on **tahoebio/EmeraldBay**
(a separate 1.8 M-cell, 5-day atlas sharing the 5 representative cell lines), re-estimating σ² and m from
EmeraldBay's *own* cells: the embedding (a **frozen** 2000-HVG PCA(50) set for reproducibility) and within-condition σ² ≈ 0.938 are computed over the **full atlas** (52 lines, 1.83M cells) — a like-for-like match to Tahoe's own **within-condition σ² ≈ 0.9567** (to ~2%) (Tahoe's marginal 2.406 is only a conservative bound), a like-for-like cross-platform match; re-estimation is still
necessary in general, confirming σ² is the platform-specific input).

Per-cell coordinates were retained for the **full 52-line atlas**, so the gating verification and the full-population falsification (below) span the whole atlas; the four illustrative **held-out angular-error curves** here use the 5 lines shared with Tahoe (141,720 cells; 5-line σ²=0.85) for direct cross-comparison. They compare the realized RMS angle from n subsampled cells to the closed form
`θ²(n)=tr(PΣP)/m²·(1/n−1/N)` (finite-population-corrected):

| group (N, m) | mean rel. err | fitted slope vs `tr(PΣP)/m²` | R² |
|---|---|---|---|
| DMSO_T0 × HS-578T (1067, m=5.89) | **0.8%** | 0.906 vs 0.918 | **0.9999** |
| DMSO_T0 × HEC-1-A (435, m=3.58) | 0.6% | 1.831 vs 1.845 | 0.9999 |
| DMSO_T0 × BT-474 (461, m=3.35) | 0.8% | 1.709 vs 1.780 | 0.9977 |
| Irinotecan × HS-578T (563, m=1.96, drug) | 5.3% | 10.08 vs 12.81 (low-SNR) | 0.9786 |

`realized² ∝ (1/n − 1/N)` is linear-through-origin with **R² > 0.997 (reference groups; 0.979 for the low-magnitude drug group)**, and the fitted slope matches the
anisotropic `tr(PΣP)/m²` to **within ~1–4% for the reference groups** (the low-magnitude drug group deviates ~21% at small
n — the expected breakdown of the first-order law at lower SNR, §5 of `THEORY.md`). **Regime gating
(full atlas):** of 3,971 (condition × line) groups across all 52 lines, 347 (8.7%) are predicted
OVER-sampled, and the downsample-and-measure check now covers **every one — 347/347 (100%)** meet the
tolerance at n* when downsampled (`pass5_gating_full.py`).

This is distilled into `fixtures/emeraldbay_calibration.json` and asserted by
`tests/test_emeraldbay_integration.py` (runs in CI, no network). The streaming calibrator that produced
it is `scripts/emeraldbay_recompute/` (one-time streaming job over 57.7 GB / 116 shards).

### Full-atlas falsification — both atlases, from raw (`docs/FALSIFICATION.md`)

Beyond the four illustrative curves, the **parameter-free** slope test was run across the *entire*
population of both atlases. All **95,624,334 Tahoe cells** were re-streamed from raw and **56,195**
conditions tested; the full 52-line EmeraldBay atlas gave **1,064** groups. In the theory's validity
regime (`ρ²≥3`) the a-priori slope `tr(PΣP)/m²` matches the realized fitted slope to a median ratio of
**0.94** (1,790 Tahoe conditions) and **0.98** (32 EmeraldBay groups), each at **R²≈0.999** with no
fitted parameter; outside it the deficit grows monotonically with `1/ρ²`, exactly as the second-order
theory predicts. The same re-projection reconfirms the noise constants on *every* cell — Tahoe marginal
(conservative-bound) **σ²=2.4158** (all 95.6M cells) vs the 2.406 marginal calibration (0.4%), and EmeraldBay within-condition **σ²=0.938**
computed over all 1.83M cells, matching its frozen-basis calibration. And the OVER/UNDER gating decision
holds on the primary atlas too: **5,503/5,503** predicted-OVER Tahoe conditions met the tolerance at n*
(per-line-mean baseline, θ=0.1; `pass4c_gating.py`).

## Cross-modality scale-up: gene perturbation (X-Atlas/Orion) — `docs/ORION_GENE_PERTURBATION.md`

The quota is a theorem about centroid-**direction** estimation, so it is indifferent to *what* moves the
centroid. We tested that on the largest public **genome-wide CRISPRi Perturb-seq** atlas, **X-Atlas/Orion**
(~8M cells, 18,903 gene knockdowns × 2 lines; SLAF/Lance format), streaming **all 46.5 billion expression
entries** from raw and mapping the framework one-to-one: condition = (target gene × line), control arm = the
pooled **Non-Targeting** centroid, `m = ‖μ_knockdown − μ_NTC‖`. Findings (`scripts/orion_recompute/`,
`fixtures/orion_*`):

- **σ² transfers across modality.** Within-condition `σ² = 0.913` (HCT116) / `1.033` (HEK293T) — same ~1.0 band
  as the chemical atlases (Tahoe 0.957, EmeraldBay 0.938), so the calibrated law is `n* ≈ 9,000/m²` in both
  modalities. On Orion the *marginal* σ² ≈ the within-condition value (unlike Tahoe's 2.5× gap) — a direct
  fingerprint that most knockdowns barely move the transcriptome.
- **No knockdown is over-sampled.** Genome-wide at θ=0.1 rad: HCT116 **0.0% OVER / 28.3% UNDER / 71.7% Ghost**
  (median deficit ~2,260×); HEK293T **0.0% / 49.4% / 50.6%** (~248×). Gene knockdowns move the transcriptome
  ~4–10× less than drugs (median m 0.14–0.35 vs 1.27), so `n* ∝ 1/m²` explodes. Robust to tolerance (≤2% OVER
  even at θ=0.3 rad).
- **Self-validation.** The strongest-magnitude knockdowns are core essential genes (ribosomal proteins,
  ribosome biogenesis, Pol II, splicing, nucleoporins, translation initiation) — known biology recovered at the
  top of the spectrum. And **HEK293T is where the anisotropic form earns its keep** (aniso/iso ratio 0.66 vs
  ~0.93 elsewhere): its noise is single-axis-dominated, the regime `docs/THEORY.md` §8 predicts.

This is an *application* of the proven framework with a transferred calibration (Phase E of the manuscript), not
an independent re-proof; the direct falsification is supplied on the essential-gene screens below.

## Validated essential-gene extension (TRADE) — `docs/TRADE_EXTENSION_PLAN.md`

The Orion streaming run could apply the quota but not re-run the held-out falsification (only summary stats are
kept). We close that on the **TRADE** essential-gene CRISPRi screens (Nadig et al., *Nat Genet* 2025; GEO
**GSE264667**) — **Jurkat-Essential** and **HepG2-Essential**, 2,393 DepMap common-essential genes/line. At
~2.5e5 cells/line these load in RAM, so per-cell coordinates are kept and the parameter-free downsample-and-measure
test runs directly (`scripts/trade_recompute/`, `fixtures/trade_*`):

- **First *validated* gene-perturbation result.** For strong knockdowns the a-priori slope `tr(PΣP)/m²` matches the
  realized downsampling slope at **1.02 (Jurkat) / 0.96 (HepG2), R²=0.996**, no fitted parameter — as tight as the
  chemical Phase C — with the predicted low-SNR breakdown (Spearman 0.92 / 0.79).
- **A new, inverted regime.** Essential-gene knockdowns are *strong* (median m 1.8–2.2, above Tahoe's 1.27) but
  *shallowly sampled* (median 48–85 cells/gene). So 61–71% are detectable, yet still only 0.2–0.6% are over-sampled
  and 14–18% ghost: the binding constraint flips from magnitude (Orion) to **depth**.
- **Biology per line.** Strongest knockdowns are the RNA exosome (EXOSC2–9)/splicing/Mediator in Jurkat, and
  translation/proteasome/ribosome (EIF2S1, PSMB5, RPL\*) in HepG2. σ² here is higher (1.5–1.9; the large-NTC-pool
  quota constant 7.3k–9.1k/m² returns to the ~9,000 band) — noted as a looser transfer than the other datasets.

## The dual-sided framework: one threshold, two ledgers

`calculate_optimal_resource_allocation(single_cell_variance, num_dimensions, perturbation_magnitude,
tolerance=0.01, baseline_cells_per_well=1296, complexity="linear")` reads the **same** saturation
threshold `n*` two ways, and returns `required_cells_per_well` and `dry_lab_compute_reduction_ratio`.

### Wet-lab ledger — budget gating & multiplexing
`n*` is the cells/arm needed to resolve a drug's direction to `tolerance`. If your current depth `N0`
exceeds `n*`, the surplus `N0 − n*` reads are wasted on a saturated estimate. The freed budget can
instead **multiplex more conditions per lane**: `wet_lab_multiplex_gain = N0 / n*` extra wells fit in
the same per-lane read budget. (If `n* > N0` you are *under*-sampled — sequence deeper, don't multiplex.)

### Dry-lab ledger — provably-safe downsampling, with honest scope
**Claim (provable).** For **centroid / pseudobulk** analyses — the drug–drug similarity graph, drug-
level PCA/clustering — downsampling a well from `N0` to `n*` cells keeps each perturbation centroid
within angular tolerance `tolerance` *by construction* (that is exactly what `n*` solves). So the
similarity matrix entries, and the drug-level embedding built from them, are preserved to that
tolerance. This is a direct corollary of the Delta-method bound above, not a separate assumption.

**Cost saved is POLYNOMIAL, not exponential.** Reducing per-well cells from `N0` to `n*` reduces:

```
RAM / storage / streaming / PCA-fit   (cost ~ N)     :  saved = 1 - (n*/N0)
cell-cell pairwise / kernels / k-NN   (cost ~ N^2)   :  saved = 1 - (n*/N0)^2
neighbour graphs / UMAP               (cost ~ N logN):  saved = 1 - (n* log n*)/(N0 log N0)
```

No standard single-cell operation is exponential in N, so no amount of downsampling yields an
exponential speedup — the gains are linear-to-quadratic. `dry_lab_compute_reduction_ratio` reports the
`max(0, …)` of the chosen `complexity` model (0 when the well is under-sampled — nothing to cut).

**Scope limit (do not over-claim).** `n*` governs **direction/centroid** information only. It does
**not** certify preservation of **cell-resolution** structure — per-cell UMAP local neighbourhoods,
rare-population or cell-type detection, trajectory branch points — which are set by *local density*,
not by this angular threshold, and **can be distorted** by downsampling. Use `n*` to thin
**over-sampled wells before pseudobulk/graph analysis**, not to subsample an atlas before cell-level
embedding. Global PCA eigenvalue spectra converge with sampling error (random-matrix theory), so they
are preserved *approximately*, not exactly.

### Calibration reality check (Tahoe-100M, θ=0.1 rad, N0=1,296 cells/condition — reproduces `python src/calibrate.py`)

| signature (5 µM, NCI-H460) | m | n*/arm | regime | dry save (lin / quad) | wet multiplex |
|---|---:|---:|:--:|:--:|:--:|
| Homoharringtonine | 15.53 | 39 | OVER | **97% / 99.9%** | **33.4×** |
| Panobinostat | 9.97 | 94 | OVER | 93% / 99% | 13.7× |
| Trametinib | 4.39 | 485 | OVER | 63% / 86% | 2.7× |
| median condition (panel-wide) | 1.27 | 5,795 | **UNDER** | 0% / 0% | — |

**Honest finding:** at a tight tolerance the atlas is *under-sampled* for **89.4%** of (drug×dose×line)
conditions (those cells are **not** redundant); only strong perturbers in responsive lines at high dose
are over-sampled enough to downsample. The calculator tells you exactly which regime each condition is
in (`regime` field). Full per-condition data → `fixtures/tahoe_quota_per_condition.csv`.

## How a wet-lab uses it to cut sequencing cost

1. **Estimate `sigma^2` once** for your platform/pipeline: run a pilot (any condition), embed the
   same way you will analyze (HVG → PCA `d`), take the mean per-dim variance of the per-cell PCA
   coordinates *within* a single condition. (For the Tahoe HVG/PCA(50) pipeline this within-condition
   value is ~0.9567; the marginal per-cell σ²≈2.406, mixing conditions, is a conservative bound.)
2. **Estimate the smallest effect size `m` you must resolve** — the perturbation magnitude (PCA-space
   `||treated − control centroid||`) of your *weakest drug of interest* from a pilot or a prior atlas.
3. **Pick a tolerance** — the angular precision you need to call two drugs "same direction" (0.1 rad
   ≈ 6° is a sensible default; tighter for fine MoA separation).
4. **`n* = calculate_experimental_cell_quota(sigma2, d, m, tolerance)`** → cells per arm. Multiply by
   (#drugs × #doses × #lines × 2 arms) for the run, and by your per-cell sequencing cost.
5. **Spend cells where they matter:** strong perturbers need far fewer cells (`1/m^2`), so weak
   signatures dominate the budget — either accept a looser tolerance for them or drop them.

```bash
python src/calibrate.py                 # demo on the cached Tahoe array
DATA=/path/to/reference_perturbations.npz python src/calibrate.py
python src/calculator.py                # self-check
```

> **API note (current).** `n_c` is REQUIRED. Ignoring the control-pool size is exactly the error this
> paper documents: against the largest atlas's real shared DMSO pool, ~69% of conditions are
> control-pool-limited (`n* = inf`, unresolvable at any treated depth) — invisible to any equal-arm
> formula. The library calls `src/engine.py`; it never re-implements the quota. (The old
> `calculate_experimental_cell_quota(...)` equal-arm entry point is superseded by the functions below.)

```python
import math, numpy as np
from src.calculator import cell_quota, cell_quota_large_pool, cell_quota_equal_arm
v = np.zeros(50); v[0] = 2.97                    # perturbation vector (m = ||v||); Sigma may be diagonal
Sigma = 0.9567 * np.ones(50)                     # within-condition per-cell variance

cell_quota(v, Sigma, tolerance=0.1, control_pool_size=1_500_000)  # RECOMMENDED two-arm; n_c REQUIRED
cell_quota(v, Sigma, tolerance=0.1, control_pool_size=1514)       # small shared vehicle -> may be math.inf (POOL-LIMITED)
cell_quota_large_pool(v, Sigma, tolerance=0.1)                    # n_c -> inf limit  (= (d-1)sigma^2/(m^2 t^2))
cell_quota_equal_arm(v, Sigma, tolerance=0.1)                     # matched 1:1 vehicle only (= 2x large-pool)
```

## Layout

```
sample_sufficiency_calculator/
├── src/
│   ├── calculator.py   # calculate_experimental_cell_quota(...)        -> n* (isotropic)
│   │                   # calculate_cell_quota_anisotropic(...)         -> n* (anisotropic + tail)
│   │                   # calculate_optimal_resource_allocation(...)    -> {n*, dry-lab reduction, ...}
│   │                   # rms_angular_error(...)
│   └── calibrate.py    # verifiable demo: isotropic, dual-sided, and anisotropic (self-contained)
├── docs/
│   └── THEORY.md       # rigorous anisotropic proof (Delta method, generalized-chi2, tail bound)
├── fixtures/
│   └── tahoe_calibration.json   # per-PC variances + example drug vectors (self-contained demo)
├── tests/
│   ├── verify_theory.py         # standalone SymPy Jacobian proof + Monte-Carlo harness
│   └── test_calculator.py       # pytest suite (units + symbolic + Monte-Carlo + tail coverage)
├── .github/workflows/ci.yml     # GitHub Actions: pytest on py3.10/3.11/3.12 + harness
├── conftest.py                  # puts src/ and tests/ on sys.path for pytest
├── pytest.ini                   # pytest config
├── requirements-dev.txt         # numpy, sympy, pytest
├── README.md
└── .gitignore
```

## Verification & CI

```bash
pip install -r requirements-dev.txt
pytest                       # full suite (12 tests)
python tests/verify_theory.py   # standalone formal+empirical harness
```

The proof is checked three ways (all asserted; run on every push via GitHub Actions over Python
3.10/3.11/3.12 — see `.github/workflows/ci.yml`):

1. **Symbolic (SymPy).** Auto-differentiates `g(x)=x/‖x‖` and asserts the Jacobian equals the analytical
   `(1/m)(I − uuᵀ)` — residual simplifies to the zero matrix (+ numeric cross-check at a concrete `v`).
2. **Empirical (Monte Carlo).** 50,000 noise draws `e ~ N(0, 2Σ/n)` with an anisotropic Σ (decaying
   spectrum in a random basis); the **exact** angle `arccos(⟨v,v̂⟩/(‖v‖‖v̂‖))` matches the closed form
   `2 tr(PΣP)/(n m²)` to **<1%** (last run: 0.13%).
2b. **Tail coverage.** Running the Laurent–Massart confidence quota `n*_δ`, the empirical
   `P(θ>θ*)` stays `≤ δ` (e.g. 0.34% ≤ 10%).
2c. **Second-order term (`verify_theory.py` Layer 4).** The §5 expansion of `E[θ²]` (including the
   arc-vs-tangent correction) matches the exact geometric angle to **<0.1%**, and the older
   lever-arm-only form (`= E[tan²θ]`) is shown to be materially worse — asserted so it can't regress.
3. **Implementation.** `src/calculator.py` returns the analytical `n*` exactly (0.0000%).

Unit tests additionally pin the scaling laws (`1/m²`, `σ²`, `1/θ²`), the anisotropic→isotropic
reduction, the perpendicular-only property (variance along the signal is inert), the dual-allocation
regimes, and input validation.

### Formal verification in Lean 4 / Mathlib (`lean/`)

The **deterministic** backbone is additionally machine-checked in Lean 4 against Mathlib (`v4.31.0`),
with **no `sorry`** and only the three standard axioms (`#print axioms`): Lemma 1 (the normalization
Jacobian `Dg(v)=(1/m)P`), the key trace identity `tr(PΣP)=trΣ−uᵀΣu`, and its isotropic (`→σ²(d−1)`)
and PCA/diagonal (`→Σ_k(1−u_k²)ℓ_k`) reductions. The probabilistic results (E[θ²], the tail) remain
empirical (above). See [`lean/README.md`](lean/README.md) and `docs/THEORY.md` §10.

```bash
cd lean && lake exe cache get && lake build   # checks all Lean proofs (seconds, given the cache)
```

## License / status

Research utility calibrated and validated on two public reference perturbation atlases.
Provided as-is; validate `sigma^2` on your own platform before planning a screen.

## II. Dataset Inventory & Scale Audit

> **Data representation.** Tahoe-100M provides the **full raw gene-expression count matrix**; it is
> not limited to a reduced representation. The **`d = 50` PCA embedding** used throughout this tool is
> a *downstream representation* we compute from that raw data for the geometric power-analysis
> framework (direction estimation lives in PCA space, where `σ²` and the perturbation magnitudes `m`
> are measured). The mathematics and quotas are stated in that fixed embedding, but the underlying
> atlas is the full-resolution raw expression data.

*Scale ledger of two large-scale, independent, multi-line reference perturbation atlases used purely for empirical validation. **GLOBAL** rows are dataset-level (published/metadata); **CAPTURED** rows are computed from the cached arrays of the companion calibration pipeline. Per-well cell-count distributions are reported for CAPTURED data only, since global per-well counts are not in either atlas's metadata. Cell lines are reported by their real identifiers HS-578T, AN3-CA, HEC-1-A, BT-474 and C-33 A (the five established cancer lines common to both atlases).*

### A. Global inventory

| Reference atlas | Total cells (global) | Unique perturbagens | Unique cell lines | Wells (line × condition) |
|---|---:|---:|---:|---:|
| **Tahoe-100M** | ~100,000,000 | 379 | 50 | ~56,850 |
| **EmeraldBay** | ~1,831,648 | 27 molecules (93 conditions) | 52 | 4,992 |

### B. Per-well cell-count distribution (CAPTURED data only)

| Atlas (captured scope) | Wells | Cells | Min | Median | Mean | Max |
|---|---:|---:|---:|---:|---:|---:|
| Tahoe-100M — 5 dose-matched plates | 20,000 | 33,450,029 | 1 | 1,192 | 1,673 | 23,043 |
| EmeraldBay — five cell lines (HS-578T, AN3-CA, HEC-1-A, BT-474, C-33 A) | 430 | 141,720 | 11 | 310 | 330 | 1,978 |

*Tahoe-100M captured median **1,192** cells/well sits far below the median quota n★≈5,794 required at θ=0.1 rad — i.e. most wells are under-sampled at tight tolerance (see `docs/SCALE_AUDIT.md`). Tahoe-100M conditions are essentially unreplicated (R=1 for 96% of drug-line conditions, max R=3), so pooling replicate wells does not change this gating — only ~0.08% of conditions cross UNDER->OVER when pooled (docs/SCALE_AUDIT.md sec 5).*

### C. Cell-line cross-tabulation — vehicle vs active perturbations (CAPTURED)

| Cell line | Atlas | Vehicle wells | Vehicle cells | Active wells | Active cells | Total cells |
|---|---|---:|---:|---:|---:|---:|
| HS-578T | Tahoe-100M | 5 | 12,414 | 395 | 592,831 | 605,245 |
| AN3-CA | Tahoe-100M | 5 | 4,665 | 395 | 205,950 | 210,615 |
| HEC-1-A | Tahoe-100M | 5 | 15,032 | 395 | 808,138 | 823,170 |
| BT-474 | Tahoe-100M | 5 | 10,018 | 395 | 431,939 | 441,957 |
| C-33 A | Tahoe-100M | 5 | 11,593 | 395 | 530,215 | 541,808 |
| HS-578T | EmeraldBay | 2 | 2,584 | 84 | 30,688 | 33,272 |
| AN3-CA | EmeraldBay | 2 | 1,544 | 84 | 15,558 | 17,102 |
| HEC-1-A | EmeraldBay | 2 | 2,296 | 84 | 33,122 | 35,418 |
| BT-474 | EmeraldBay | 2 | 2,439 | 84 | 36,837 | 39,276 |
| C-33 A | EmeraldBay | 2 | 1,010 | 84 | 15,642 | 16,652 |

