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
   with high-variance directions the true `n*` departs from the isotropic estimate (here per-dim
   variance ranges 0.03–31.7, mean 7.66).
2. `tolerance` is an **RMS angular SD in radians**. Sub-degree tolerances on noisy single-cell data
   demand very large `n` — e.g. `tolerance=0.01` (~0.57°) needs ~10^5–10^6 cells/arm. Realistic
   values are typically **0.05–0.2 rad (≈3–11°)**.
3. `sigma^2` is tied to a fixed pipeline (HVG selection, normalization, PCA `d`). Re-derive it if you
   change the embedding.

---

## Calibration on Tahoe-100M (verifiable demonstration)

`src/calibrate.py` reads real perturbation magnitudes from the batch-clean drug-similarity array
(a batch-clean 292-perturbagen reference array) and the per-cell PCA variance
`sigma^2 = 7.66` (mean over 50 dims, calibrated from a reference-atlas subsample of 60k cells
projected through the PCA). It contrasts a **strong** signature — **Resveratrol** (`m = 2.97`), the
moderate-magnitude pathway-modulator profile — against a **weak** signature (`m = 1.33`, 25th
percentile):

| tolerance | Resveratrol (m=2.97) | weak signature (m=1.33) |
|----------:|---------------------:|------------------------:|
| 0.05 rad (2.9°) | ~34,100 cells/arm | ~170,300 |
| 0.10 rad (5.7°) | ~8,500 | ~42,600 |
| 0.20 rad (11.5°) | ~2,100 | ~10,600 |
| 0.30 rad (17.2°) | ~950 | ~4,700 |

The **weak/strong ratio is exactly `(m_strong/m_weak)^2 ≈ 5.0`** at every tolerance — the `m^2` law,
which is the verification check `calibrate.py` prints.

> **Provenance & honesty.** This formula was derived for this tool (Delta method, above); it was not
> taken from a prior design. The numbers in the table are **computed from the cached data**, not
> assumed. An external draft circulated quoting `n*=361` / `n*=2661` is **not reproducible** from the
> real per-cell variance and effect sizes and is deliberately **not** reproduced here.

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

**Honest empirical finding (Tahoe-100M):** on this atlas the anisotropic *mean-quota* correction is
**small (~2–3%; ratio 0.98)** — drug directions carry only ~2–5% of total variance along themselves, so
they are not aligned with the dominant cell-cycle/lineage PCs, and with d=50 no single axis can move the
trace much. The anisotropic machinery's real value here is (i) `d_eff ≈ 28 ≪ 49` and (ii) the rigorous
**tail quota** (95%-confident Resveratrol = 78k vs mean 41k cells/arm, a 1.9× safety factor). On data
where perturbations *do* align with high-variance axes the correction is large — the ratio formula in
§8 of `THEORY.md` says exactly when.

## External validation on an independent atlas (EmeraldBay)

The formula is a theorem (verified symbolically + by Monte Carlo), so a second dataset can't make it
*more* true — but it **can** test whether the CLT/Gaussian-centroid *assumptions* hold on real,
independent single cells, and whether the calibration transfers. We did this on **tahoebio/EmeraldBay**
(a separate 1.8 M-cell, 5-day atlas sharing the 5 representative cell lines), re-estimating σ² and m from
EmeraldBay's *own* cells (within-condition σ² ≈ 2.1, vs Tahoe's marginal 7.66 — re-estimation is
necessary, confirming σ² is the platform-specific input).

**Held-out angular-error curves** (realized RMS angle from n subsampled cells vs the closed form
`θ²(n)=tr(PΣP)/m²·(1/n−1/N)`, finite-population-corrected):

| group (N, m) | mean rel. err | fitted slope vs `tr(PΣP)/m²` | R² |
|---|---|---|---|
| HS-578T (1067, 12.4) | **1.4%** | 0.769 vs 0.763 | **0.998** |
| AN3-CA (752, 6.4) | **0.8%** | 1.379 vs 1.375 | **0.9996** |
| HEC-1-A (1117, 3.0) | 2.4% | 6.66 vs 7.82 | 0.995 |

`realized² ∝ (1/n − 1/N)` is linear-through-origin with **R² > 0.99**, and the fitted slope matches the
anisotropic `tr(PΣP)/m²` to **<1% for strong signatures** (the weaker m≈3 group deviates ~15% at small
n — the expected breakdown of the first-order law at lower SNR, §5 of `THEORY.md`). **Regime gating:** of
101 wells, the 15 predicted OVER-sampled all met the tolerance when downsampled to n* (**100%**).

This is distilled into `fixtures/emeraldbay_calibration.json` and asserted by
`tests/test_emeraldbay_integration.py` (runs in CI, no network). The streaming calibrator that produced
it is the companion calibration pipeline (one-time 58 GB job; figure
`emeraldbay_heldout_validation.png`).

## The dual-sided framework: one threshold, two ledgers

`calculate_optimal_resource_allocation(single_cell_variance, num_dimensions, perturbation_magnitude,
tolerance=0.01, baseline_cells_per_well=1394, complexity="linear")` reads the **same** saturation
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

### Calibration reality check (Tahoe-100M, θ=0.1 rad, N0=1394 cells/well)

| signature | m | n*/well | regime | dry save (lin / quad) | wet multiplex |
|---|---:|---:|:--:|:--:|:--:|
| Resveratrol (moderate-signal modulator) | 2.97 | 8,518 | **UNDER** | 0% / 0% | — |
| weak (25th pct) | 1.33 | 42,571 | **UNDER** | 0% / 0% | — |
| strong cytotoxic (max) | 14.35 | 365 | OVER | **74% / 93%** | **3.8×** |

**Honest finding:** at a tight tolerance this atlas is *under-sampled* for moderate/weak signatures
(those cells are **not** redundant), and only **very strong** perturbers are over-sampled enough to
downsample. The "massive redundant matrix" intuition holds only for high-magnitude drugs or loose
tolerances — the calculator tells you exactly which regime you are in (`regime` field).

## How a wet-lab uses it to cut sequencing cost

1. **Estimate `sigma^2` once** for your platform/pipeline: run a pilot (any condition), embed the
   same way you will analyze (HVG → PCA `d`), take the mean per-dim variance of the per-cell PCA
   coordinates. (For the Tahoe HVG/PCA(50) pipeline this is ~7.66.)
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

```python
from src.calculator import calculate_experimental_cell_quota
n = calculate_experimental_cell_quota(single_cell_variance=7.66, num_dimensions=50,
                                      perturbation_magnitude=2.97, tolerance=0.1)
# -> ~8518 cells per arm
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
3. **Implementation.** `src/calculator.py` returns the analytical `n*` exactly (0.0000%).

Unit tests additionally pin the scaling laws (`1/m²`, `σ²`, `1/θ²`), the anisotropic→isotropic
reduction, the perpendicular-only property (variance along the signal is inert), the dual-allocation
regimes, and input validation.

## License / status

Research utility calibrated and validated on two public reference perturbation atlases.
Provided as-is; validate `sigma^2` on your own platform before planning a screen.

## II. Dataset Inventory & Scale Audit

*Scale ledger of two large-scale, independent, multi-line reference perturbation atlases used purely for empirical validation. **GLOBAL** rows are dataset-level (published/metadata); **CAPTURED** rows are computed from the cached arrays of the companion calibration pipeline. Per-well cell-count distributions are reported for CAPTURED data only, since global per-well counts are not in either atlas's metadata. Cell lines are reported by their real identifiers HS-578T, AN3-CA, HEC-1-A, BT-474 and C-33 A (the five established cancer lines common to both atlases).*

### A. Global inventory

| Reference atlas | Total cells (global) | Unique perturbagens | Unique cell lines | Wells (line × condition) |
|---|---:|---:|---:|---:|
| **Tahoe-100M** | ~100,000,000 | 379 | 50 | ~56,850 |
| **EmeraldBay** | ~1,831,756 | 27 molecules (93 conditions) | 52 | 4,992 |

### B. Per-well cell-count distribution (CAPTURED data only)

| Atlas (captured scope) | Wells | Cells | Min | Median | Mean | Max |
|---|---:|---:|---:|---:|---:|---:|
| Tahoe-100M — 5 dose-matched plates | 20,000 | 33,450,029 | 1 | 1,192 | 1,673 | 23,043 |
| EmeraldBay — five cell lines (HS-578T, AN3-CA, HEC-1-A, BT-474, C-33 A) | 430 | 141,720 | 11 | 310 | 330 | 1,978 |

*Tahoe-100M captured median **1,192** cells/well sits far below the median quota n★≈23,934 required at θ=0.1 rad — i.e. most wells are under-sampled at tight tolerance (see `docs/SCALE_AUDIT.md`). Tahoe-100M conditions are essentially unreplicated (R=1 for 96% of drug-line conditions, max R=3), so pooling replicate wells does not change this gating — only ~0.08% of conditions cross UNDER->OVER when pooled (docs/SCALE_AUDIT.md sec 5).*

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

