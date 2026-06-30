# Sample-Sufficiency Calculator for Single-Cell Perturbation-Direction Screens

A small, dependency-light power calculator that answers one practical question for single-cell drug
screens that compare **perturbation *directions*** (cosine / correlation of treated-vs-control
centroids in a PCA embedding):

> **How many cells per arm do I need to pin down a drug's transcriptional direction to within a
> given angular tolerance?**

It returns a **cell quota `n*`** so experimentalists can size (and cost) a screen instead of
over-sequencing. Strong perturbations need quadratically fewer cells than weak ones.

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
│   ├── calculator.py   # calculate_experimental_cell_quota(...) + rms_angular_error(...)
│   └── calibrate.py    # empirical calibration / verifiable demo on the cached array
├── README.md
└── .gitignore
```

## License / status

Research utility derived from the Tahoe-100M OBGYN drug-similarity work (RESEARCH_LOG §26–§31).
Provided as-is; validate `sigma^2` on your own platform before planning a screen.
