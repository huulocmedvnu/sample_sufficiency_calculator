# Audit log: findings #1-#9, gaps 1-3, and the shared-DMSO thesis

A complete, self-contained record of the audit that overhauled this study. Written so a reader who
was never involved can understand what was wrong, what changed, and why. Authoritative current state:
`docs/MANUSCRIPT.md` + `src/engine.py`. Architecture: `docs/ARCHITECTURE.md`. Original
findings ledger: `docs/AUDIT_DENOMINATORS.md`. Rollback point: git tag `pre-audit-original`.

## 1. Executive summary

**Start:** the manuscript reported that of the 56,827 Tahoe-100M (drug x dose x line) conditions,
**10.6% are over-sampled and 89% under-sampled** at a 5.7 degree tolerance, with the quota
`n* = 9,376/m^2` (a matched-vehicle equal-arm constant). EmeraldBay's spectrum was `4.6/84.0/11.4`
(over/under/ghost). The chemical and genetic spectra were computed by two separate, drifted code paths.

**End:** one engine computes all six screens. Under Tahoe's **real** matched vehicle (the plate-shared
DMSO pool), the headline is **10.3% over-sampled, 13.8% treated-depth-limited, 1.0% ghost, 69.2%
control-pool-limited, 5.8% not-detectable** (94.2% detectable). The thesis changed from *"most
conditions need more treated cells"* to *"on the largest atlas the binding constraint is the size of
the shared control pool, not treated depth"* -- the quota's `m_min` floor, previously a caveat, is
measured to bind on the majority of the largest available atlas. Tolerance is the second lever
(loosening 5.7 to 17 degrees frees the pool-limited fraction 69% to 13%).

The genetic numbers did not change: they were computed correctly all along, and the engine reproduces
them exactly, which is what validates the rewrite.

## 2. Findings and gaps

| # | description | evidence | quantitative impact | fix | commit(s) |
|---|-------------|----------|---------------------|-----|-----------|
| 1 | EmeraldBay spectrum pooled over dose (drug-name x line), mislabelled "condition x line" | `emeraldbay_falsification_full.json::per_group_slope` (27 drug names, not 93 drug x dose) | over 4.6% (dose-pooled, N>=300) -> ~1-2.7% on the true unit | engine uses (sample x line) = (drug x dose x line), 4,912 groups | 2adfe25, 7b9faec |
| 2 | three undeclared cell-count denominators (Tahoe N0>=1, EmeraldBay N>=300, genetic N>=25) | `make_spectrum_unified` / `per_group_slope` / `--min-cells` | inflated EmeraldBay over% ~3x | spectrum over ALL conditions; genetic keeps a >=25-cell measurement QC; declared in Methods | 7b9faec |
| 3 | "4,992 EmeraldBay groups" counted 80 empty combinations | `pseudobulk.npz` (52x96=4,992; 4,912 non-empty) | n 4,992 -> 4,912 | stated | 2adfe25 |
| 4/5 | N>=300 threshold had no statistical basis (held-out fixture-reuse artifact) | `min_cells_slope=300` in the held-out fixture | see #1 | dropped; N>=400 is a held-out-only requirement | 7b9faec |
| 6 | chemical spectrum computed ISOTROPICALLY, not the anisotropic rule the paper is named for | `tahoe_recompute/pass3_quota.py:65` raw isotropic n*=C/m^2, 9,376 reproduces to 0 rel-err | none on aggregate (measured: see gap 1) | engine uses `tr(P.Sigma.P)` with the real Sigma | 9b105b0, 769f146 |
| 7 | genetic bias-corrected the magnitude (m^2 - tr(S)), chemical did not | `orion_recompute/pass3_quota_full.py:58` vs `tahoe .../pass3_quota.py:65` | Tahoe ghost 0.4 -> 2.9% under correction (over ~unchanged) | bias-correction applied to EVERY dataset (decision reversed, section 4) | 7b9faec |
| 8 | genetic detection floor snr>1.5, chemical a baseless factor 1.25 (snr>1.118) | `pass3_quota_full.py:76` (`>1.5`) vs `make_spectrum_unified.py:50` (`1.25`) | EmeraldBay "detectable" 96.8% (1.118) -> 87.6% (1.5) | one snr>1.5 threshold for every dataset; 1.25 removed | 7b9faec |
| 9 (CRITICAL) | Tahoe's control is a small SHARED pool, not matched vehicle | each `DMSO_TF` centroid is shared by a median of 94 conditions at n_c~1,514 cells (Phase-0 stream of `pass4/coords.f32`) | equal-arm C=9,376 and 10.6%/89% invalid; m_min median 1.75 > median m 1.27 -> 69% control-pool-limited | headline references the real shared vehicle; per-line mean kept as a sensitivity | 769f146, 7b9faec |
| gap 1 | chemical Sigma was `sigma2 * I` (isotropic); tr(P.Sigma.P) was a constant 46.88 for all conditions | measured min=max=46.88 | over 21.5 -> 21.7% under real Sigma (negligible) | real diagonal Sigma from per-condition sufficient stats (`chemical_within_cov.json`, `pass7_within_cov.py`); measured finding: Sigma anisotropic 3.6-5x but tr(P.Sigma.P) varies only +-1.5-3% | 769f146 |
| gap 2 | Tahoe estimand: per-line mean vs matched vehicle | Phase-0 (finding #9) | per-line mean 21.7% over / 0 pool-limited; vehicle 10.3% over / 69.2% pool-limited | chose the vehicle (A2 estimand, section 4) | 769f146, 7b9faec |
| gap 3 | sigma2 hardcoded for chemical, recomputed for genetic | `chemical.py` `SIGMA2={...}` (removed) | none (Tahoe recompute reproduces 0.9567) | sigma2 = mean(ell) from the real diagonal Sigma; hardcode dropped | 769f146 |

## 3. Root cause (the most important section)

Every finding #1-#8 is a symptom of ONE structural defect: **the sufficiency calculation was
duplicated across a chemical code path and a genetic code path, which then drifted apart.** The
genetic path (`orion_recompute`, `trade_recompute`) was written carefully and stayed correct; the
chemical path (`make_spectrum_unified.py`, `make_unified_calibration_table.py`,
`tahoe_recompute/pass3_quota.py`) accumulated shortcuts: it hardcoded an isotropic Sigma (#6), skipped
the bias-correction (#7), used a baseless 1.25 detection constant instead of snr>1.5 (#8), hardcoded an
equal-arm quota constant instead of using the real control-pool size (#9/gap 2), and hardcoded sigma2
(gap 3). Two implementations of the same formula could not be kept in agreement by review alone.

**The fix is architectural, not numeric.** `src/engine.py` is now the single place the math lives.
Two thin loaders (`src/pipelines/{chemical,genetic}.py`) only assemble a standard per-condition
structure (`mu_t, mu_c, n_t, n_c, Sigma`) and hand it to the engine; they contain no formula. The
engine takes no `if dataset` / `if modality` branch and no hardcoded arm factor -- every constant
follows from `sigma2`, `theta`, and the REAL `n_c`. `tests/test_no_duplicate_math.py` fails if any
pipeline carries math, if the quota formula appears outside `engine.py`, or if a dead constant
(9,376/9,192) reappears in `src/`. The engine reproduces the genetic ground truth exactly, so a single
gate now protects both modalities. This is why re-duplication cannot silently recur.

## 4. Scientific decisions (who decided, why, what was rejected)

All decisions were made by the study author; the assistant executed and, where the author's stated
expectations were mutually inconsistent, stopped and surfaced the fork with measured numbers.

- **Bias-correction (#7): reversed from (i) to (ii).** Round 1 chose (i) *delete the claim* (the code
  did not bias-correct the chemical magnitudes). On measuring the cost of the alternative, the shift
  was tiny (Tahoe over 10.59 -> 10.42, EmeraldBay 1.00 -> 1.00; only ghost moved), so the author
  reversed to (ii): **bias-correct every dataset**, matching what the Methods always described and the
  genetic path always did. Restored the deleted Statistical-analysis sentence and the GATA1/EIF4A3
  Discussion passage, which are now correct.
- **EmeraldBay unit:** (sample x line) = (drug x dose x line), 4,912 groups. Rejected: the
  `per_group_slope` fixture (dose-pooled drug-name x line), which is for the held-out test only.
- **Spectrum denominator:** no cell-count filter for the chemical spectra; genetic keeps a >=25-cell
  measurement QC (a centroid needs cells to exist). Rejected: the undeclared N>=300 / N>=1 mix.
- **n_c:** use each dataset's REAL control-pool size, never a hardcoded arm factor. The two-arm quota
  `n_t* = 1/(m^2 theta^2 / tr(P.Sigma.P) - 1/n_c)` then yields equal-arm (n_c=n_t) or large-pool
  (n_c large) automatically.
- **Tahoe estimand (gap 2, the biggest decision): shared DMSO vehicle (A), not per-line mean (B).**
  The direction rule is defined on `v = mu_t - mu_vehicle` (estimand A2), the displacement from
  *untreated* cells. Tahoe's real vehicle is the plate-shared DMSO pool (n_c~1,514, shared across ~94
  conditions), giving m_min median 1.75 > median m 1.27 -> **69% control-pool-limited, 10.3% over**.
  Option B (per-line mean) references the average over that line's ~94 treatments -- a different
  quantity ("difference from the average drug"), whose near-noiseless m_min~0.05 claims a control no
  real experiment has, and which gives 21.7% over / 0 pool-limited. B was rejected as the wrong
  estimand and demoted to a sensitivity analysis. This is the decision that changed the paper's thesis.
- **Detection floor (#8):** snr = m_raw/sqrt(tr(S)) > 1.5 for every dataset; the chemical-only 1.25
  constant was removed as baseless (same class of error as the N>=300 threshold).

## 5. Old -> new, all six screens

Regime over ALL conditions (over/under/ghost/pool-limited/not-detectable); genetic over/under are also
reported among detectable (Table 4). Chemical C is the large-control bound `(d-1)sigma^2/theta^2`.

| screen | OLD (pre-audit) | NEW (unified engine) |
|--------|-----------------|----------------------|
| **Tahoe-100M** | 10.6 over / 89.0 under / 0.4 ghost, detectable 98.4, C=9,376 (equal-arm), median n* 5,794 | **10.3 over / 13.8 under / 1.0 ghost / 69.2 pool-limited / 5.8 not-det**, detectable 94.2, C=4,688 bound, median resolvable n* 2,401, m_min 1.75 |
| **EmeraldBay** | 4.6 over / 84.0 under / 11.4 ghost (dose-pooled, N>=300), detectable 96.8, C=9,192 | **2.7 over / 82.5 under / 2.0 ghost / 0.5 pool-limited / 12.4 not-det**, detectable 87.6, C=4,495, sigma^2 0.9174 |
| Orion HCT116 | 14.7 det, 0.0 over(det) | 14.7 det, 0.0 over(det) -- unchanged |
| Orion HEK293T | 35.3 det, 0.1 over(det) | 35.3 det, 0.1 over(det) -- unchanged |
| TRADE Jurkat | 70.6 det, 1.2 over(det) | 70.6 det, 1.2 over(det) -- unchanged |
| TRADE HepG2 | 60.9 det, 4.4 over(det) | 60.9 det, 4.4 over(det) -- unchanged |

Tolerance lever (Tahoe, real vehicle): pool-limited 69.2 -> 46.0 -> 29.3 -> 18.7 -> 12.8% as theta
goes 0.10 -> 0.15 -> 0.20 -> 0.25 -> 0.30 rad (Table 6).

## 6. Verification

- **Engine-vs-genetic gate:** the engine reproduces sigma^2 0.913/1.033/1.50/1.87, detectable
  14.7/35.3/70.6/60.9, and over-among-detectable 0.0/0.1/1.2/4.4 exactly
  (`tests/test_engine_golden.py::test_engine_reproduces_genetic`). A wrong engine fails here.
- **OVER/UNDER downsample-and-measure** on the 5,851 predicted-over-sampled conditions (real vehicle):
  **5,824/5,851 (99.5%)** reach RMS angle <= 1.05 theta. The 27 exceptions are the most aggressively
  compressed strong-effect conditions where the second-order arc-versus-tangent curvature of
  Supplementary Note S1 appears (realized <= 7.0 degrees vs the 5.7 target); accepted as the documented
  regime, written into Results. Fixture: `fixtures/verification_over.json`.
- **Test suite:** 23 tests green, including the engine gate and `test_no_duplicate_math.py`.
- **Build:** `manuscript.pdf` rebuilds to ~1.15 MB with all five figures.

## 7. What is NOT settled / residual risk (do not hide these)

- **EmeraldBay sigma^2 = 0.9174 (recomputed) vs 0.938 (old hardcode).** The ~2% difference is a
  grouping effect: 0.9174 pools the within-condition variance at (sample x line) granularity from the
  committed pseudobulk sumsq; the old 0.938 came from a different pooling. The engine uses **0.9174**
  (the recompute). Tahoe's recompute, by contrast, reproduces the committed 0.9567 exactly. A reviewer
  may ask which EmeraldBay grouping is canonical.
- **Chemical Sigma is DIAGONAL, not full 50x50.** The chemical sufficient statistics store only
  per-PC variances, so `chemical.py` builds `Sigma = diag(ell_within)`; the genetic TRADE path carries
  the full 50x50 Sigma. For the chemical atlases this means `u^T Sigma u` uses only the diagonal, so
  off-diagonal covariance between PCs is ignored in `tr(P.Sigma.P)`. Measured effect is small (tr(PSP)
  varies +-1.5-3%), but it is a genuine approximation, declared in Methods.
- **27/5,851 fail the OVER gate at 1.05 theta.** Explained as the Note-S1 curvature regime (realized
  <= 7.0 degrees), but a reviewer may push on whether the classification should carry a curvature
  safety margin for the strongest-effect conditions.
- **The per-line-mean sensitivity uses a grand-mean baseline** (average over ~94 treatments), which is
  a biased reference; it is reported only as a sensitivity, not the headline.
- **Historical/auxiliary docs** carry pre-audit numbers behind SUPERSEDED banners (not the manuscript).

## 8. Round 2 audit: purge of legacy quota math (branch `rewrite/purge-legacy-math`)

**Trigger.** `scripts/emeraldbay_recompute/pass6_invariance.py` computed the downsampling-invariance
quota with a LOCAL equal-arm formula `2 tr(P.Sigma.P)/(m^2 theta^2)` -- not the engine -- and fed the
result (`n* = 184`, exemplar Galunisertib x AN3-CA) straight into manuscript Table S3 / the
"Downsampling stability" paragraph. `test_no_duplicate_math.py` never caught it: it (a) scanned only
`src/`, never `scripts/`, and (b) matched two exact engine substrings, not the equal-arm spelling.

**New guard.** `tests/test_no_duplicate_math.py` rewritten to scan `src/` AND `scripts/` (every `.py`
except `src/engine.py`) and to flag by MEANING, not spelling: the quota formula (tr(P.Sigma.P) divided
by m^2 theta^2 -- solving for n*, distinct from the held-out measurement law tr/m^2*(1/n-1/N) which is
allowed), the bias-correction `m^2 - tr(S)`, the detection floor `snr > const`, and the dead constants
9376/9192/4688/4495/4475 and the 1.25/2.25/300 thresholds. First run flagged 10 violations in 6 files.

**Migrated to the engine (local quota/bias/detection removed):**
- `scripts/emeraldbay_recompute/pass6_invariance.py` -- now calls `engine.compute` (verified: re-run).
- `scripts/orion_recompute/pass3_quota_full.py` -- re-run on both lines, reproduces the genetic golden
  EXACTLY (sigma^2 0.9133/1.0334, detectable 14.7/35.3, over(det) 0.0/0.1, median n* 21376/7592).
- `scripts/trade_recompute/pass2_falsification.py` -- re-run, reproduces the held-out slopes EXACTLY
  (strong 1.0203 Jurkat / 0.9607 HepG2, R^2 0.996).
- `scripts/trade_recompute/pass1_qc_basis_quota.py` -- quota routed through the engine; streaming/QC/
  basis.npz untouched (needs the raw h5ad to run, so migrated + syntax-checked, not re-run here).
- `scripts/applications/reliability_and_cost.py` -- rewritten engine-based (no atlas cache); reproduces
  reliability 10.3/7.0, 41.9/36.5, 65.9/61.7 and budget 6.6M/2.4M/2.76x.
- `scripts/make_manuscript_figures.py` -- deleted the retired dead functions (`fig4_emeraldbay`,
  `fig6_trade_validation`, `_showcase_curves`) that carried an inline bias-correction.

**Deleted (superseded, outputs unconsumed):** `src/quota_arm.py`,
`scripts/orion_recompute/pass2_quota.py` (pilot), `scripts/recompute_genetic_arm.py` (post-processor;
`genetic_two_arm.json` unconsumed).

**Numbers that were legacy-fed and are now engine-corrected:**
| item | legacy (equal-arm/raw) | engine (two-arm, bias-corrected) |
|---|---|---|
| invariance OVER quota / depth n_d | 184 (130 in one stale spot) | **92** |
| invariance UNDER exemplar | Galunisertib x AN3-CA (n*=22,475) -- now POOL-LIMITED, invalid | **Gemcitabine x MIA PaCa-2**, n*=42,666, UNDER |
| invariance OVER metrics | drift <=0.2%, cos SD 0.005 | drift <=0.6%, cos SD 0.008 |
| invariance UNDER metrics | +11% (0.57->0.63), cos 0.08-0.09, SD 0.047 | +19% (0.78->0.93), cos 0.10-0.12, SD 0.10 |
| Table S1 median m | raw (5.88, 4.24, ...) | bias-corrected (5.79, 4.14, ...); Crizotinib row was stale -> (S)-Crizotinib |
| Table S1 counts | 330 pool-maj / 316 median-below-floor | **331 / 327** |
| Table S2 median m | raw (1.27/1.20/1.35) | bias-corrected (1.213/1.143/1.296) |
| homoharringtonine NCI-H460 5uM n* | 20 | **19** (multiplex ~300 -> ~320-fold) |
| homoharringtonine NCI-H661 5uM | n*=476, UNDER | **POOL-LIMITED** (that line's DMSO pool = 112 cells, m_min=6.5) |
| 90th-pct budget quota | 11,013 | 10,973 (flat/adaptive 6.6M/2.4M unchanged) |

**Fixtures regenerated from the engine:** `orion_{HCT116,HEK293T}_quota.csv` + `_summary.json`,
`trade_{jurkat,hepg2}_falsification.json`, `tahoe_per_dose.csv` (median_m now bias-corrected),
`tahoe_quota_per_condition.csv`, `tahoe_applications.json`, `outputs/emeraldbay_invariance.json`.

**Unchanged (verified NO drift):** `fixtures/unified_spectrum.json` -- every headline conclusion
intact (69.2% pool-limited, 10.3% over, sigma^2 transfer, held-out slopes). Genetic numbers unchanged.

**Verification:** new guard 0 violations; 21 tests green; PDF ~1.15 MB / 35 pp / 5 figures; DOCX built.
Every manuscript number re-derived from a live `engine.compute` stdout (theta-lever Table 6, Table S1/
S2/S3, budget, k-NN edges, resolved fractions, homoharringtonine all reproduced).

**Residual (not blocking):** `fixtures/tahoe_constants.json` still carries pre-audit *thesis* fields
(pct_OVER 10.6, `n*=9,376/m^2`, median 5794) but only its `sigma2=0.9567` is read (by
`test_sigma2_golden`); `scripts/tahoe_recompute/{pass5_within_sigma,pass6_within_recalibrate}.py` are the
old 3-regime tahoe passes (need the 18 GB coords cache to run) and are superseded for the per-condition
CSV by the engine -- they carry no quota formula, so the guard passes them.

## 9. Round 3 audit: the equal-arm recurrence, the AST guard, and the public API

After §8, tightening the (still string-based) guard kept surfacing MORE files that computed the quota
locally -- the same equal-arm `2*tr(PSP)/(m^2 theta^2)`, just spelled `THETA` (uppercase) so the
lowercase-only regex missed it. Rather than patch regexes, a **manual AST-level read of every `.py`**
(src/ + scripts/) found the true set. Six files self-computed the quota/regime outside `engine.py`:
`tahoe_recompute/pass4c_gating`, `emeraldbay_recompute/{pass5_gating_full,pass3_validate}`,
`applications/moa_recovery`, `trade_recompute/pass1` (a diagnostic print), and a dead `_regime_counts`
in `make_manuscript_figures`. All migrated to the engine (or deleted) and re-run where the cache exists.

**Two manuscript numbers were equal-arm-stale** and are now engine-computed (two-arm, real per-line pool):

| gating | old (equal-arm) | new (engine) |
|--------|-----------------|--------------|
| Tahoe full-atlas OVER, downsample-verified (theta=0.1, per-line-mean, N>=400) | 5,503/5,503 | **13,964/13,964 (100%)** |
| EmeraldBay full-atlas OVER, downsample-verified (theta=0.20) | 347/347 | **706/706 (100%)** |

The count roughly doubled because the correct large-pool quota is half the equal-arm one (more
conditions are over-sampled) -- and every one still meets tolerance when downsampled to n*. The paper's
central evidence (predicted-OVER conditions are downsample-stable) holds; only the count changed.
Also updated: k-NN Jaccard null 0.005 -> 0.006 (adaptive vs flat, re-run via engine); "roughly
sixteen-fold" -> "twenty-fold"; EmeraldBay §365 8.7% -> 17.8% OVER.

**Single-owner refactor.** `src/engine.py` now exposes `quota_two_arm(trPSP, m2, n_c, theta)` and
`m_min_floor(...)`; `engine.compute` and `src/calculator.py` both call them (spectrum unchanged, golden
test passes). A second guard, `tests/test_only_engine_imports_math.py`, reasons over the parse tree
(AST) and fails if any src/ or scripts/ file except `engine.py` computes a quota formula locally. The
first guard was made case-insensitive and given the retired-constant check.

**Public API fixed (`src/calculator.py`).** It used to default to the equal-arm quota, so a README user
who omitted the control-pool size silently got a retired number. Now `cell_quota(v, Sigma, tolerance,
control_pool_size)` REQUIRES `n_c`, returns `math.inf` for control-pool-limited, and provides
`cell_quota_equal_arm` / `cell_quota_large_pool` as explicitly-named limits; it delegates to the engine.
Tests, `calibrate.py`, `verify_theory.py`, and the README example were updated to the new signatures.

**Cleanup.** Deleted the retired-thesis hallucination sources -- 9 orphan docs, ~10 stale fixtures,
`agents/` and `dataset_inventory.py` (old numbers in prompt/string literals), `quota_arm.py`, three old
tahoe passes, and figure dead code. Regenerated `docs/SUPPLEMENT.md` (constants of record) and
`fixtures/tahoe_{constants,calibration,within_sigma}.json` + `{orion,trade}_*_quota.csv` from the engine.
`docs/ARCHITECTURE.md` has a "Retired artifacts" section listing every deleted file + fingerprint.

**Final state:** 26 tests green; both guards 0; `unified_spectrum.json` unchanged (69.2% pool-limited /
10.3% over / sigma^2 transfer / held-out slopes 0.94/0.98, 1.02/0.96 all intact); PDF ~1.15 MB / 5
figures (confirmed by rendering pages) + DOCX. Committed on `master` (HEAD `dbf301c`), pushed to origin.

## 9. Round 4 (2026-09-26): the DMSO_TF pool overwrite

`src/pipelines/chemical.py` assigned `dmso_c[(plate, line)] = coords[i]` for each DMSO_TF sample, so on
Tahoe plates that carry 2-3 vehicle wells only the LAST well survived: the headline "n_c ~ 1,514 cells"
was an artifact, not the design. Fixed to pool every DMSO_TF well (cell-weighted). Regenerated:
`run_unified_spectrum.py`, `applications/reliability_and_cost.py` (now also emits the Table 6 sweep),
all figures, `verification_over.json` (new committed generator `pass4d_gating_vehicle.py`; the old
5,824/5,851 number had no script), manuscript, SUPPLEMENT, README, HANDOFF, THEORY_PRIMER.

| Tahoe (shared vehicle, theta 0.1) | before | after |
|---|--:|--:|
| median n_c | 1,514 | 3,113 |
| control-pool-limited % | 69.2 | 55.4 |
| treated-depth-limited % | 13.8 | 25.2 |
| over % / ghost % / not-det % | 10.3 / 1.0 / 5.8 | 11.4 / 2.4 / 5.6 |
| median m_min / m_corr | 1.75 / 1.216 | 1.22 / 1.096 |
| median resolvable n* | 2,401 | 3,622 |
| OVER gate | 5,824/5,851 | 6,442/6,467 |

The other five screens and the per-line-mean sensitivity are unchanged. The 25 gate exceptions are
strong-effect conditions, 23 of them in NCI-H460; the earlier "arc-vs-tangent curvature" explanation
was withdrawn because that term has the opposite sign (a shortfall). The same audit recorded design
issues NOT fixed here: the held-out downsample test is a finite-population identity; the OVER gate
shares the control centroid and so never tests the 1/n_c term; EmeraldBay references the per-line mean
while Tahoe references the vehicle (against its own vehicle EmeraldBay is ~62% pool-limited);
between-well DMSO variance (~2.6x the sampling floor) is absent from Sigma; the MoA-recovery null
(`fixtures/tahoe_moa_recovery.json`) is not reported in the manuscript.
