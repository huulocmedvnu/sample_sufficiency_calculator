# Sample-sufficiency calculator for single-cell perturbation screens

[![CI](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml/badge.svg)](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Lean 4](https://img.shields.io/badge/proofs-Lean%204%20%2B%20Mathlib-green.svg)](lean/)

**How many cells does a perturbation condition need before the *direction* of its transcriptomic effect can be trusted?**
This repository provides a closed-form answer, the code that calibrates it on public single-cell atlases, and the
tests and machine-checked proofs behind it.

It accompanies the manuscript *Closed-form sample size calculation for estimating perturbation directions in single-cell screens* (Tran and Nguyen, 2026, in preparation). The manuscript itself is not
part of this repository.

## The problem

Perturbation screens (Perturb-seq, chemical screens such as Tahoe-100M) summarise each condition by the displacement
of its treated centroid from a control centroid, `v = μ_treated − μ_control`, in a shared low-dimensional embedding.
Mechanism is then read from the **direction** of `v`: two treatments that move the transcriptome the same way are
inferred to act through related programs. No sample-size rule states how many cells that direction needs, and the
usual power calculations do not see the control arm at all.

## The quota

Only the sampling noise **perpendicular** to `v` rotates its direction. With the tangent-space projector
`P = I − uuᵀ` (`u = v/‖v‖`), the per-cell covariance `Σ`, the effect magnitude `m = ‖v‖` and a target
root-mean-square angular error `θ★`, the number of treated cells that holds the error at `θ★` is

```
n_t* = 1 / ( m² θ★² / tr(PΣP)  −  1 / n_c )
```

where `n_c` is the size of the control pool the condition is referenced to. Three consequences follow:

| situation | what the formula gives |
|---|---|
| large shared control, `n_c → ∞` | `n_t* = tr(PΣP) / (m² θ★²)` |
| matched arms, `n_c = n_t` | `n_t* = 2 · tr(PΣP) / (m² θ★²)` |
| effect below the floor `m_min = √( tr(PΣP) / (n_c θ★²) )` | `n_t* = ∞`: the control alone rotates the direction by more than `θ★`, so **no treated depth resolves it** |

The last case is the reason `n_c` is a required argument everywhere in this package. On the largest public chemical
atlas, referenced to its real plate-shared DMSO vehicle, it is the dominant regime.

A Laurent–Massart tail bound extends the mean quota to a confidence quota `n*_δ` with `P(θ > θ★) ≤ δ`. The deterministic
core (the projector and the identity `tr(PΣP) = tr Σ − uᵀΣu`) is proved in Lean 4 against Mathlib, and the probabilistic
statements are checked symbolically and by Monte-Carlo (`tests/verify_theory.py`).

## Quick start

```bash
git clone https://github.com/huulocmedvnu/sample_sufficiency_calculator.git
cd sample_sufficiency_calculator
pip install -r requirements-dev.txt        # numpy, sympy, pytest
export PYTHONPATH=src
python src/calibrate.py                    # worked example on the committed Tahoe-100M calibration
```

```python
import numpy as np
from calculator import cell_quota, cell_quota_report

v = np.zeros(50); v[0] = 2.97            # perturbation vector in the embedding (m = 2.97)
Sigma = 0.9567 * np.ones(50)             # within-condition per-cell variance per PC (diagonal or full 50x50)

cell_quota(v, Sigma, tolerance=0.1, control_pool_size=200_000)   # 533 treated cells (large control pool)
cell_quota(v, Sigma, tolerance=0.1, control_pool_size=3_113)     # 641 (small shared vehicle)
cell_quota(v, Sigma, tolerance=0.1, control_pool_size=300)       # inf: control-pool-limited

cell_quota_report(v, Sigma, 0.1, control_pool_size=3_113, confidence=0.05)
# {'required_cells_treated': 641, 'control_pool_limited': False, 'm_min': 1.23,
#  'effective_dimensions': 49.0, 'required_cells_treated_confident': 1036, ...}
```

`tolerance` is in radians (0.1 rad = 5.7°). `cell_quota_isotropic(σ², d, m, tolerance, control_pool_size)` is the scalar
form when only a per-cell variance is known. Every function raises if `control_pool_size` is missing or not positive.

### What you need to supply

| input | meaning | where it comes from |
|---|---|---|
| `Σ` | within-condition per-cell covariance in the embedding | pooled residual covariance of a pilot, after centring each condition on its own mean (`src/sigma2_canonical.py`) |
| `v` or `m` | perturbation vector, or its bias-corrected magnitude | pilot centroids, `m² = max(0, m_raw² − tr Σ (1/n_t + 1/n_c))` |
| `n_c` | size of the control pool the condition is referenced to | the experimental design |
| `θ★` | angular tolerance | the smallest inter-mechanism angle the analysis must separate, 0.1 rad by convention |

`σ²` depends on the normalisation and embedding. In the fixed recipe used here (normalise to 10⁴, log1p, 2,000 highly
variable genes, PCA to 50) it is close to 1.0 on chemical and genome-wide CRISPRi platforms alike, but it must be
re-estimated for any other pipeline.

## Calibration on public data

The engine was run on four public atlases (six screens, two modalities, three platforms) reprocessed from raw counts
through one embedding. All numbers below come from `fixtures/unified_spectrum.json`, at `θ★ = 0.1` rad.

| screen | conditions | control (`n_c`) | σ² | detectable % | over-sampled % | control-pool-limited % |
|---|--:|---|--:|--:|--:|--:|
| Tahoe-100M | 56,827 | shared DMSO vehicle (≈3,113) | 0.957 | 94.4 | 11.4 | **55.4** |
| EmeraldBay | 4,912 | per-line mean (≈22,800) | 0.917 | 87.6 | 2.7 | 0.5 |
| X-Atlas/Orion HCT116 | 17,585 | pooled NTC (165,562) | 0.913 | 14.7 | 0.0 | 0.0 |
| X-Atlas/Orion HEK293T | 17,856 | pooled NTC (218,838) | 1.033 | 35.3 | 0.0 | 0.0 |
| TRADE Jurkat | 2,187 | pooled NTC (11,514) | 1.498 | 70.6 | 0.8 | 0.7 |
| TRADE HepG2 | 1,883 | pooled NTC (4,380) | 1.867 | 60.9 | 2.7 | 1.2 |

Two findings organise the results. On Tahoe-100M the plate-shared vehicle of about 3,100 cells serves about 94
conditions, its sampling floor `m_min = 1.22` sits above the median effect magnitude (1.10), and 55% of conditions
cannot resolve their direction at any treated depth. Genome-wide CRISPRi screens fail for a different reason: most
knockdowns are too weak to detect, and the detectable ones are sequenced 17 to 35 times too shallow.

The law behind the quota was tested on every atlas by downsampling its own cells: the parameter-free slope
`tr(PΣP)/m²` is recovered in the first-order regime (realised/predicted ratio 0.94 on Tahoe-100M, 0.98 on EmeraldBay,
1.02 and 0.96 on the TRADE screens), and every condition predicted over-sampled met the tolerance when downsampled
to its quota (Tahoe 13,964 of 13,964, EmeraldBay 706 of 706). Details: `docs/FALSIFICATION.md`.

## Repository layout

```
src/
  engine.py             quota_two_arm(): the single definition of the quota; compute(): the six-screen spectrum
  calculator.py         public API (cell_quota, _large_pool, _equal_arm, _isotropic, _report, resource_allocation)
  sigma2_canonical.py   canonical within-condition σ² estimator
  calibrate.py          worked example on the committed calibration
  pipelines/            chemical.py, genetic.py: assemble per-condition sufficient statistics (no math)
scripts/
  run_unified_spectrum.py          engine over all six screens -> fixtures/unified_spectrum.json
  {tahoe,emeraldbay,orion,trade}_recompute/   from-raw streaming, held-out tests, gating (need the raw atlases)
  applications/                    tolerance sweep, budget, reliability, mechanism-label recovery
  make_figures_plotly.py           Figures 2-5 from the committed fixtures (plotly)
fixtures/               committed derived data: spectrum, per-screen calibrations, held-out results, figure inputs
tests/                  26 tests: golden values, calculator suite, tail coverage, symbolic + Monte-Carlo, two guards
lean/                   Lean 4 / Mathlib proofs of the deterministic core
docs/                   THEORY.md (derivation), THEORY_PRIMER.md, ARCHITECTURE.md, SUPPLEMENT.md (constants of record),
                        FALSIFICATION.md, AUDIT_LOG.md, AUDIT_DENOMINATORS.md
```

## Reproducing the results

| what | command | needs |
|---|---|---|
| tests and guards | `PYTHONPATH=src python -m pytest tests/ -q` | Python 3.10 to 3.12 |
| symbolic and Monte-Carlo checks | `python tests/verify_theory.py` | sympy, numpy |
| Lean proofs | `cd lean && lake exe cache get && lake build` | Lean 4 toolchain (pinned in `lean/lean-toolchain`) |
| figures | `python3 scripts/make_figures_plotly.py` | plotly ≥ 5.20 with kaleido 0.2.x, numpy, scipy, pillow |
| six-screen spectrum | `python scripts/run_unified_spectrum.py` | per-condition caches from the recompute pipelines |
| from-raw recompute | `scripts/*_recompute/` | the raw atlases (Tahoe-100M, EmeraldBay, X-Atlas/Orion, TRADE) |

Everything down to the figures runs from the committed fixtures. The recompute pipelines stream the raw atlases
(tens of gigabytes) and are provided for full reproducibility rather than everyday use.

## Design guarantees

- **One definition of the quota.** `src/engine.py::quota_two_arm` is the only place the formula exists. Two CI guards
  fail on any regression: `tests/test_no_duplicate_math.py` scans for re-derived quota expressions and for retired
  constants, and `tests/test_only_engine_imports_math.py` checks the import graph.
- **No dataset-specific branches.** The engine takes one standard per-condition structure (`mu_t, mu_c, n_t, n_c,
  Sigma`) and never branches on dataset or modality.
- **Golden locks.** `tests/test_engine_golden.py` and `tests/test_sigma2_golden.py` pin the published per-screen values.
- **Provenance.** `docs/AUDIT_LOG.md` records the audit that retired an earlier equal-arm version of the analysis and
  how each finding was fixed.

## Limitations

The noise model is the within-condition covariance scaled by `1/n`. Between-well variation, plate effects and
cell-cycle composition do not shrink with cell number and are not included, so the quota understates the requirement
where they dominate. The first-order guarantee holds when the signal-to-noise ratio `ρ² = m² / (uᵀΣu)` is well above
one. Below the detection floor the calculator returns no count. The quota certifies the direction of a pseudobulk
centroid only, not cell-level structure.

## Citation

Please cite the software as described in `CITATION.cff`. The archived releases are on Zenodo (DOI added with the first release).

The accompanying manuscript: Tran T.H.L., Nguyen T.V. (2026). *Closed-form sample size calculation for estimating perturbation directions in single-cell screens.* In preparation.

## License

Code, scripts and derived fixtures are released under the MIT License (see `LICENSE`). The manuscript text and figures
are not part of this repository.

## Contact

Tran Thai Huu Loc, Department of Obstetrics and Gynecology, Hong Hung Hospital, Tay Ninh, Viet Nam.
tranthaihuuloc2@gmail.com. ORCID 0009-0009-1810-010X.
