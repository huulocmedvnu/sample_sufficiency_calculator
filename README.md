# Sample-Sufficiency Calculator for Single-Cell Perturbation-Direction Screens

[![CI](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml/badge.svg)](https://github.com/huulocmedvnu/sample_sufficiency_calculator/actions/workflows/ci.yml)

> **Single source of truth.** Every number in this study is produced by `src/engine.py` (via
> `scripts/run_unified_spectrum.py` → `fixtures/unified_spectrum.json`). The authoritative write-up is
> `docs/MANUSCRIPT_DEEPSEEK.md`; the constants of record are `docs/SUPPLEMENT.md`; the audit narrative
> is `docs/AUDIT_LOG.md`. A prior *equal-arm* thesis (the constants `9,376/m²`, `10.6% over`, `89%
> under`) was retired by the control-pool audit and now survives **only** in `AUDIT_LOG.md` /
> `AUDIT_DENOMINATORS.md` / `ARCHITECTURE.md` / this file — anywhere else, treat those numbers as stale.

A small, dependency-light calculator that answers one question for a direction-based perturbation
screen: **how many cells does it take to resolve a perturbation's transcriptomic *direction* to a target
angular precision — and, just as important, when is that impossible at any depth because the *control*
is too small?**

Mechanism is increasingly read from the *direction* of a treatment's displacement vector
`v = μ_treated − μ_control` in a shared embedding, not its magnitude. Only noise **perpendicular** to `v`
rotates that direction, so the cells needed scale with the tangent-space noise `tr(PΣP)` and inversely
with the squared effect size `m² = ‖v‖²`.

## The quota (two-arm, with the real control-pool size)

The estimator subtracts two cell centroids, so its sampling covariance is `Σ/n_t + Σ/n_c`. Holding the
RMS angular error at a tolerance `θ★` and solving for the treated arm gives

```
              1
n_t* = ───────────────────────────          P = I − uuᵀ,   u = v/‖v‖
       m² θ★² / tr(PΣP)  −  1/n_c
```

- reduces to the **matched equal-arm** value `2·tr(PΣP)/(m²θ★²)` when `n_c = n_t`,
- and to the **large-shared-pool** limit `tr(PΣP)/(m²θ★²)` as `n_c → ∞`;
- returns **`n* = ∞` (CONTROL-POOL-LIMITED)** when `m` is below the floor
  `m_min = √(tr(PΣP)/(n_c θ★²))` — no treated depth resolves the direction; the lever is a bigger
  control pool or a looser tolerance.

**`n_c` is not optional.** Ignoring the control-pool size is the mistake this project documents: against
the largest atlas's real shared vehicle, most conditions are control-pool-limited — invisible to any
equal-arm formula. The quota is defined in exactly one place, `src/engine.py::quota_two_arm`; the public
`src/calculator.py` and the internal `engine.compute` both call it, and two CI guards fail on any file
that recomputes it (`tests/test_no_duplicate_math.py`, `tests/test_only_engine_imports_math.py`).

A tail bound (Laurent–Massart) extends the mean quota to a **confidence quota** `n*_δ` guaranteeing
`P(θ>θ★) ≤ δ`. The deterministic core (the projector Jacobian and `tr(PΣP)=trΣ−uᵀΣu`) is machine-checked
in **Lean 4/Mathlib**; the probabilistic parts are verified symbolically and by Monte-Carlo.

## API

```python
import math, numpy as np
from src.calculator import cell_quota, cell_quota_large_pool, cell_quota_equal_arm, cell_quota_report

v = np.zeros(50); v[0] = 2.97          # perturbation vector (m = ‖v‖)
Sigma = 0.9567 * np.ones(50)           # within-condition per-cell variance (diagonal or dxd)

cell_quota(v, Sigma, tolerance=0.1, control_pool_size=1_500_000)   # RECOMMENDED; n_c REQUIRED
cell_quota(v, Sigma, tolerance=0.1, control_pool_size=3113)        # small shared vehicle → may be math.inf
cell_quota_large_pool(v, Sigma, tolerance=0.1)                     # n_c → ∞ limit
cell_quota_equal_arm(v, Sigma, tolerance=0.1)                      # matched 1:1 vehicle only (= 2× large-pool)
cell_quota_report(v, Sigma, 0.1, control_pool_size=3113, confidence=0.05)  # + m_min, d_eff, tail quota
```

`cell_quota_isotropic(σ², d, m, tolerance, control_pool_size)` is the scalar convenience form. All raise
if `control_pool_size` is missing or ≤ 0. Run the calibrated demo: `python src/calibrate.py`.

## Calibration & validation (all from `fixtures/unified_spectrum.json`, θ★ = 0.1 rad)

Four public atlases (six screens) across two modalities and three platforms, reprocessed from raw through one fixed
embedding (normalize 1e4 → log1p → HVG 2000 → PCA 50). `C` is the large-control bound `(d−1)σ²/θ★²`.

| screen | within σ² | C bound | det % | over % | pool-lim % | control |
|---|--:|--:|--:|--:|--:|---|
| **Tahoe-100M** (56,827 cond.) | 0.9567 | 4,688 | 94.4 | **11.4** | **55.4** | real shared DMSO vehicle (2-3 wells pooled, n_c ≈ 3,113 / ~94 cond.) |
| **EmeraldBay** (4,912) | 0.9174 | 4,495 | 87.6 | 2.7 | 0.5 | per-line mean (large pool) |
| Orion HCT116 (17,585) | 0.9133 | 4,475 | 14.7 | 0.0* | 0.0 | shared NTC pool |
| Orion HEK293T (17,856) | 1.0334 | 5,063 | 35.3 | 0.1* | 0.0 | shared NTC pool |
| TRADE Jurkat (2,187) | 1.4975 | 7,338 | 70.6 | 1.2* | 0.7 | shared NTC pool |
| TRADE HepG2 (1,883) | 1.8667 | 9,147 | 60.9 | 4.4* | 1.2 | shared NTC pool |

*genetic over% is among the *detectable* subset (Table 4).

**Headline (Tahoe, real shared DMSO vehicle).** The plate-shared DMSO pool (the plate's 2-3 vehicle wells
pooled, ~3,113 cells shared across ~94 conditions) sets a floor `m_min = 1.22` above the median bias-corrected
effect (1.096), so **55.4% of conditions are control-pool-limited** — unresolvable at any treated depth — and
only 11.4% are over-sampled (25.2% treated-depth-limited). The binding constraint on the largest single-cell atlas is the size of the shared control
pool, not treated depth. Referencing each condition to a large per-line pool instead (a different
estimand, reported only as a sensitivity) removes the floor and leaves 21.7% over-sampled. Tolerance is a
second lever: loosening 5.7° → 17° frees the control-pool-limited fraction from 55% to 7%.

**σ² transfers across platform and modality** (≈ 1.0 on every chemical and genome-wide-genetic dataset).
**Held-out downsample-and-measure test** (parameter-free slope `tr(PΣP)/m²`, run on every atlas):
in-regime realized/predicted ratio 0.94 (Tahoe) / 0.98 (EmeraldBay), and 1.02 (Jurkat) / 0.96 (HepG2) on
the essential-gene screens (R² = 0.996), with the predicted low-SNR breakdown below the regime.
**Regime-gating** (predicted-OVER conditions meet tolerance when downsampled to n*): Tahoe
**13,964/13,964**, EmeraldBay **706/706** (both 100%). **Budget:** magnitude-adaptive allocation cuts a
600-condition screen from 6.6 M to 2.4 M cells (2.8×).

## Dual-sided use

- **Wet-lab (budget/multiplexing).** `n*` is the treated cells to resolve a direction to `θ★`; surplus
  over `n*` can multiplex more conditions per lane. Under-sampled or control-pool-limited conditions are
  flagged for triage instead of silent under-powering.
- **Dry-lab (safe downsampling).** For **centroid/pseudobulk** analyses (the drug-similarity graph),
  downsampling an over-acquired well to `n*` preserves its direction to `θ★` by construction; savings are
  **polynomial** (`1−(n*/N₀)^p`), never exponential. This does **not** license downsampling for
  cell-resolution structure (per-cell UMAP, rare populations) — those are set by local density.

## Layout

```
src/
  engine.py            # quota_two_arm(...) — THE quota; compute(...) — the 6-screen spectrum
  calculator.py        # standalone API (cell_quota, _equal_arm, _large_pool, _report); calls engine
  pipelines/{chemical,genetic}.py   # assemble the standard per-condition structure (no math)
  sigma2_canonical.py  # canonical within-condition σ²
scripts/
  run_unified_spectrum.py           # runs the engine over all six screens -> fixtures + Table 5
  {tahoe,emeraldbay,orion,trade}_recompute/   # from-raw streaming; held-out falsification & gating
  make_manuscript_figures.py, make_study_flowchart.py, make_heldout_figure.py
  build_manuscript.py (source .md -> submission-format .docx + .pdf), normalize_manuscript.py,
  make_reference_docx.py (Word style sheet -> docs/templates/reference.docx)
docs/   MANUSCRIPT_DEEPSEEK.md (main text source) · SUPPLEMENTARY.md (supplementary source) ·
        templates/reference.docx (journal Word styles) · SUPPLEMENT.md (constants of record)
        THEORY.md · ARCHITECTURE.md · AUDIT_LOG.md · AUDIT_DENOMINATORS.md · HANDOFF.md · FALSIFICATION.md
lean/   Lean 4 / Mathlib proof of the deterministic core
tests/  engine golden · σ² golden · calculator suite · verify_theory (SymPy + Monte-Carlo) · two guards
```

## Verification & CI

```bash
pip install -r requirements-dev.txt
export PYTHONPATH=src
python -m pytest tests/ -q                    # 26 tests (units, golden locks, tail coverage, guards)
python tests/verify_theory.py                 # standalone symbolic + Monte-Carlo harness
cd lean && lake exe cache get && lake build   # machine-checks the Lean proofs
```

CI (GitHub Actions, py3.10–3.12) runs the suite. The two guards keep the quota in one place: the semantic
one (`test_no_duplicate_math.py`) also fails if retired-thesis constants reappear outside the four exempt
docs; the AST one (`test_only_engine_imports_math.py`) fails if any `src/`/`scripts/` file (except
`engine.py`) recomputes a quota formula. Rebuild the manuscript (Briefings in Bioinformatics submission
format: A4, Times 12 pt, double-spaced, line-numbered DOCX plus a PDF rendered from it) with
`python scripts/build_manuscript.py --all` (needs `pip install pypandoc-binary python-docx` and LibreOffice;
keep the `\$` count even, or figures drop silently; verify by rendering pages).

## Status

Research utility, calibrated and validated on four public reference perturbation atlases (six screens). Provided as-is;
**re-estimate σ² and supply the real control-pool size `n_c` for your own platform** before planning a
screen.
