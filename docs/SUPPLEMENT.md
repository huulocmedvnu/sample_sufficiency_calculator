# Constants of record — anisotropic sample-sufficiency calculator

Consolidated, **current** constants for the study after the control-pool audit. Every number here is
produced by the single engine (`src/engine.py`) via `scripts/run_unified_spectrum.py` and lives in
`fixtures/unified_spectrum.json`; do not paraphrase. Authoritative narrative: `docs/MANUSCRIPT_DEEPSEEK.md`.
Audit trail: `docs/AUDIT_LOG.md`, `docs/AUDIT_DENOMINATORS.md`, `docs/ARCHITECTURE.md`.

## 1. The rule (one formula, one owner)

Two-arm cell quota (the only place the math lives is `src/engine.py`):

```
n_t* = 1 / ( m^2 theta*^2 / tr(P Sigma P) - 1/n_c ) ,   P = I - u u^T ,  u = (mu_t - mu_c)/m
```

reducing to the matched equal-arm value `2 tr(P Sigma P)/(m^2 theta*^2)` when `n_c = n_t` and to the
large-control limit `tr(P Sigma P)/(m^2 theta*^2)` as `n_c -> inf`. Bias-corrected magnitude
`m^2 -> max(0, m_raw^2 - tr(Sigma)(1/n_t + 1/n_c))`; detection floor `m_raw/sqrt(tr(S)) > 1.5`;
control-pool floor `m_min = sqrt(tr(P Sigma P)/(n_c theta*^2))` (below it `n* = inf`, control-pool-limited).
Tail-controlled quota `n*_delta = (2/m^2 theta*^2)[tr(PSP) + 2||PSP||_F sqrt(L) + 2||PSP||_op L]`,
`L = log(1/delta)` (Laurent-Massart). Embedding `d = 50`; standard tolerance `theta* = 0.1` rad (5.73 deg).

## 2. Per-screen constants (theta* = 0.1 rad; `fixtures/unified_spectrum.json`)

Regime percentages are over ALL conditions and sum to 100 with not-detectable; `C = (d-1)sigma^2/theta*^2`
is the large-control bound each screen approaches; med. m is the bias-corrected magnitude.

| screen | within sigma^2 | C bound /m^2 | conditions | det % | over % | under % | ghost % | pool-lim % | not-det % | med m | med n* | control |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|---|
| **Tahoe-100M** | **0.9567** | **4,688** | 56,827 | 94.2 | **10.3** | 13.8 | 1.0 | **69.2** | 5.8 | 1.216 | 2,401 | shared DMSO vehicle (n_c≈1,514, m_min=1.75) |
| **EmeraldBay** | **0.9174** | 4,495 | 4,912 | 87.6 | 2.7 | 82.5 | 2.0 | 0.5 | 12.4 | 1.175 | 3,673 | per-line mean (n_c≈22.8k) |
| Orion HCT116 | 0.9133 | 4,475 | 17,585 | 14.7 | 0.0 | 14.7 | 0.0 | 0.0 | 85.3 | 0.141 | 21,376 | NTC pool (165,562) |
| Orion HEK293T | 1.0334 | 5,063 | 17,856 | 35.3 | 0.0 | 35.3 | 0.0 | 0.0 | 64.7 | 0.351 | 7,592 | NTC pool (218,838) |
| TRADE Jurkat | 1.4975 | 7,338 | 2,187 | 70.6 | 0.8 | 68.3 | 0.8 | 0.7 | 29.4 | 1.841 | 1,375 | NTC pool (11,514) |
| TRADE HepG2 | 1.8667 | 9,147 | 1,883 | 60.9 | 2.7 | 56.7 | 0.3 | 1.2 | 39.1 | 2.249 | 620 | NTC pool (4,380) |

The genetic over/under split is also reported among **detectable** knockdowns (Table 4): over(det) =
0.0 (HCT116) / 0.1 (HEK293T) / 1.2 (Jurkat) / 4.4 (HepG2) %.

**Headline (Tahoe, real shared DMSO vehicle).** The plate-shared DMSO_TF pool (median ≈1,514 cells shared
across ≈94 conditions) sets a magnitude floor `m_min = 1.75` above the median bias-corrected effect
(1.216), so **69.2 % of conditions are control-pool-limited** and only 10.3 % over-sampled; median
resolvable `n* = 2,401` cells/arm vs median acquired `N0 = 1,296`. Marginal Tahoe `sigma^2 = 2.406`
(all cells) is retained only as a conservative bound.

**Sensitivity (Tahoe, per-line-mean reference, a large pool).** Removes the floor: 21.7 % over-sampled,
0 % pool-limited, 97.0 % detectable. Reported as a sensitivity, never the headline.

## 3. Tolerance lever, reliability, budget (`fixtures/tahoe_applications.json`)

- **Tolerance (Table 6).** Control-pool-limited fraction 69.2 -> 46.0 -> 29.3 -> 18.7 -> 12.8 % as
  theta* goes 0.10 -> 0.15 -> 0.20 -> 0.25 -> 0.30 rad; over-sampled 10.3 -> 25.6 -> 41.9 -> 55.6 -> 65.9 %.
- **Reliability.** Conditions resolved 10.3 / 41.9 / 65.9 % at theta* = 0.1/0.2/0.3; k-NN (k=10)
  similarity-graph edges with both endpoints resolved **7.0 / 36.5 / 61.7 %**.
- **Budget** (600-condition screen, per-line-mean reference): 90th-percentile quota 10,973 cells ->
  uniform 6.6 M cells vs quota-guided (cap 10k) 2.4 M, a **2.8x** reduction (~$1.98 M vs ~$0.72 M @ $0.30/cell).
- **Multiplex exemplar.** Homoharringtonine 5 uM x NCI-H460 (N0=6,060): m=15.5, n*=19, ~320x multiplexing gain.

## 4. Verification (`fixtures/*`, `tests/`, `outputs/`)

- **Engine gate** (`tests/test_engine_golden.py`): the engine reproduces the genetic ground truth exactly
  (sigma^2 0.9133/1.0334/1.4975/1.8667, detectable 14.7/35.3/70.6/60.9, over(det) 0.0/0.1/1.2/4.4).
- **Single-owner gate** (`tests/test_no_duplicate_math.py`): no quota/bias/detection math outside `engine.py`.
- **Held-out slope** (`docs/FALSIFICATION.md`): in-regime realized/predicted slope 0.94 (Tahoe, 1,790
  conditions) / 0.98 (EmeraldBay, 32 groups), R^2 ≈ 0.999; TRADE strong-knockdown slope 1.02 (Jurkat) /
  0.96 (HepG2), R^2 = 0.996 (`fixtures/trade_{jurkat,hepg2}_falsification.json`).
- **OVER downsample gate**: 5,824/5,851 predicted-over Tahoe conditions (shared vehicle) meet tolerance at
  n* (`fixtures/verification_over.json`); 13,964/13,964 under the per-line-mean gate; EmeraldBay 706/706.
- **Downsampling invariance** (`outputs/emeraldbay_invariance.json`, `pass6_invariance.py`): over-sampled
  DMSO_T0 x HS-578T (n*=92) stays stable when downsampled to its quota (drift ≤0.6 %, cos-SD 0.008); the
  under-sampled Gemcitabine x MIA PaCa-2 (n*=42,666) destabilises at the same depth (+19 %, cos-SD 0.10).
- **Symbolic + Monte-Carlo + Lean 4/Mathlib** proof of the deterministic core (`tests/verify_theory.py`,
  `lean/`): Jacobian residual zero; MC angle matches closed form to 0.13 %; tail coverage 0.34 % <= 10 %.

## 5. Honesty ledger

1. **Variance definition.** Headline uses the within-condition residual sigma^2 (Tahoe 0.9567, EmeraldBay
   0.9174 full-atlas / 0.8673 five-shared-line); marginal Tahoe 2.406 is a conservative bound only.
2. **Tahoe control = the real plate-shared DMSO vehicle**, small (n_c≈1,514) -> control-pool-limited
   majority. The per-line-mean reference is a different estimand, reported only as a sensitivity.
3. **Chemical Sigma is diagonal** (the chemical sufficient statistics store per-PC variances); the genetic
   TRADE path carries the full 50x50 Sigma. tr(PSP) varies only +-1.5-3 % on the chemical atlases.
4. **Held-out geometry, not biology.** Three of the four EmeraldBay illustrative groups are DMSO_T0
   reference populations; the low-SNR slope deficit is the predicted first-order breakdown (THEORY Sec 5).
5. **Scope.** n* governs centroid/direction precision only, not cell-level structure; compute gains are
   polynomial, never exponential; sigma^2 must be re-estimated per platform.
6. **TRADE sigma^2 higher (1.5-1.9), transfer looser**; strong diverse essential-gene perturbations flatten
   the PCA spectrum. TRADE is depth-limited (median 48-85 cells/gene), not magnitude-limited.

## 6. Reproduce

```bash
export PYTHONPATH=src
python3 scripts/run_unified_spectrum.py            # regenerates fixtures/unified_spectrum.json + Table 5
python3 -m pytest tests/ -q                         # engine gate + single-owner guard + golden + verify
python3 scripts/make_manuscript_figures.py && python3 scripts/make_heldout_figure.py \
  && python3 scripts/make_study_flowchart.py        # figures
python3 scripts/normalize_manuscript.py \
  && pandoc docs/manuscript.md -o manuscript.pdf --citeproc --pdf-engine=tectonic
```
