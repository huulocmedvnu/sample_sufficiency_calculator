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
(`sig_excl3_corrected.npz`, 292 drugs; RESEARCH_LOG §27/§30) and the per-cell PCA variance
`sigma^2 = 7.66` (mean over 50 dims, calibrated from the plate-6 checkpoint subsample of 60k cells
projected through the PCA). It contrasts a **strong** signature — **Resveratrol** (`m = 2.97`), the
validated functional mTORC1-inhibitor hit from §30 — against a **weak** signature (`m = 1.33`, 25th
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
| Resveratrol (validated mTORi) | 2.97 | 8,518 | **UNDER** | 0% / 0% | — |
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
DATA=/path/to/sig_excl3_corrected.npz python src/calibrate.py
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

Research utility derived from the Tahoe-100M OBGYN drug-similarity work (RESEARCH_LOG §26–§31).
Provided as-is; validate `sigma^2` on your own platform before planning a screen.
