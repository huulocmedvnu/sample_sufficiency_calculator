# Architecture: two data pipelines, one computation engine

The sufficiency calculation is computed in exactly ONE place. Everything else either feeds it a
standard per-condition structure or consumes its output.

```
raw atlas ──(embedding: normalize 1e4 -> log1p -> HVG 2000 -> PCA 50; cached)──> per-condition
   sufficient statistics
        │
        ├── src/pipelines/chemical.py   (Tahoe, EmeraldBay)
        └── src/pipelines/genetic.py    (Orion, TRADE)
                 │  emit the STANDARD structure, identical shape for every dataset:
                 │    cond_id, mu_t, mu_c, n_t, n_c (REAL control-pool size), Sigma
                 │  they do NOT compute m, tr(S), bias-correction, detectable, or n*
                 ▼
        src/engine.py  ── compute()  ── the ONLY place the math lives:
             sigma2      = tr(Sigma)/d
             m_raw       = ||mu_t - mu_c||
             tr(S)       = tr(Sigma) (1/n_t + 1/n_c)                 real n_c, no hardcoded arm factor
             m2_corr     = max(0, m_raw^2 - tr(S))                   bias-correction, every dataset
             detectable  = m_raw / sqrt(tr(S)) > 1.5                 one snr threshold, every dataset
             n*_t        = 1 / ( m2_corr theta^2 / tr(P Sigma P) - 1/n_c )   full two-arm quota
             m_min       = sqrt( tr(P Sigma P) / (n_c theta^2) )     control-pool floor
             regime      = OVER / UNDER / GHOST / POOL-LIMITED / NOT-DETECTABLE
                 ▼
        scripts/run_unified_spectrum.py  -> fixtures/unified_spectrum.json + Table 5
```

**Contract.** No `if dataset` / `if modality` and no hardcoded arm factor inside the engine. What
differs across the four datasets is only the inputs the pipelines supply (control identity, real
`n_c`, the covariance). `tests/test_no_duplicate_math.py` enforces that pipelines carry no math and
that the quota formula appears only in `engine.py`.

## The nine findings and their single root cause

The math used to be duplicated across a chemical branch and a genetic branch that drifted apart. The
genetic branch was the only one that was right; the engine makes all four go through its logic.

| # | finding | fix |
|---|---------|-----|
| 1 | EmeraldBay spectrum pooled over dose (drug-name x line) | unit is (sample x line) = (drug x dose x line), 4,912 conditions |
| 2 | three undeclared cell-count denominators | spectrum over all conditions, no filter; genetic keeps a >=25-cell measurement QC, declared |
| 3 | "4,992 groups" counted 80 empty combinations | 4,992 = 52x96 cross, 4,912 non-empty |
| 4/5 | N>=300 threshold had no basis | dropped; N>=400 is held-out only |
| 6 | chemical spectrum computed isotropically, not the anisotropic rule | engine uses tr(P Sigma P) with the REAL diagonal Sigma; measured: Sigma anisotropic 3.6-5x but tr(P Sigma P) varies only +-1.5-3% on chemical |
| 7 | genetic bias-corrected, chemical did not | m^2 -> m^2 - tr(S) for every dataset |
| 8 | genetic snr>1.5, chemical a baseless 1.25 | one snr>1.5 threshold for every dataset |
| 9 | Tahoe's control is a small shared DMSO pool, not matched vehicle | headline references the real vehicle (n_c~1,514 shared by ~94 conditions) -> 69% control-pool-limited; per-line-mean kept as a sensitivity |

Root cause of 1-8: duplicated logic in two code paths. #9 is the physics the honest engine then
exposes: with the real control-pool size, the binding constraint on the largest atlas is the shared
control pool, not treated depth.

## What is kept, split, or removed

- **Kept**: the embedding pipeline (`scripts/*_recompute/pass1..pass2`), the cached per-cell PCA-coord
  memmap, the held-out downsample-and-measure falsification (`pass4*`, an independent test of the
  measurement law), the Lean proof, and the verification suite.
- **Removed** (superseded by the engine): `make_spectrum_unified.py`,
  `make_unified_calibration_table.py` (held the dead isotropic constants 9,376 / 9,192 / 1.25).
- **Provenance**: `scripts/tahoe_recompute/pass7_within_cov.py` regenerates the real diagonal Sigma
  (`fixtures/chemical_within_cov.json`) from the cached memmap; `scripts/orion_recompute` and
  `scripts/trade_recompute` remain the streaming + original genetic implementation the engine
  reproduces (the engine-vs-genetic gate in `tests/test_engine_golden.py`).

## Retired artifacts (deleted in the legacy-purge cleanup)

These files were **deleted** because they carried retired-thesis numbers/logic and were not consumed by
the engine, manuscript, tests, or figures. They are listed so that anyone who encounters a copy (in an
old checkout, a cache, or an AI's memory) recognises it as stale. The current story is in
`docs/MANUSCRIPT_DEEPSEEK.md`, `docs/SUPPLEMENT.md`, and `src/engine.py`.

**Retired-thesis fingerprints** (appear now only in AUDIT_LOG / AUDIT_DENOMINATORS / this file / README):
`n★ = 9,376/m²` (equal-arm headline), `10.6% over / 89.0% under / 0.4% ghost` (retired within-condition
spectrum), marginal `σ²=2.406 → 23,577/m²` used as a quota with spectrum `2.5% over / 89.3% under /
8.2% ghost`, median n★ `5,794`, boundary `m=2.69`, budget `18,146` / `10.9M / 3.5M / 3.1×`, multiplex
`~156×`, EmeraldBay `4.6/84.0/11.4`, invariance exemplar `Galunisertib × AN3-CA` with `n★=184` (or
`130`) / `22,475` / `13,800`. The equal-arm form `2·tr(PΣP)/(m²θ²)` is a legitimate *special case*
(matched vehicle) but was retired **as a headline**; the headline is the two-arm engine quota.

**Deleted docs:** `DEVLOG.md`, `CASE_STUDIES.md`, `INVARIANCE.md`, `TECHNICAL_BLUEPRINT.md`,
`EXECUTIVE_SUMMARY.md`, `SCALE_AUDIT.md`, `ORION_GENE_PERTURBATION.md`,
`TRADE_EXTENSION_PLAN.md` (orphan history; superseded by the manuscript + AUDIT_LOG). (`THEORY_PRIMER.md`
was later restored and rewritten to the current two-arm thesis; it builds `theory_primer.pdf`.)

**Deleted fixtures:** `unified_calibration_table.md`, `genetic_two_arm.json`, `emeraldbay_spectrum.json`,
`heldout_summary.{json,md}`, `tahoe_condition_counts.csv`, `tahoe_layout_summary.json`,
`tahoe_per_drug.csv`, `tahoe_per_cell_line.csv` (unconsumed old-thesis artifacts).

**Deleted code:** `src/quota_arm.py`, `scripts/orion_recompute/pass2_quota.py`,
`scripts/recompute_genetic_arm.py`, `scripts/tahoe_recompute/{pass3_quota,pass5_within_sigma,pass6_within_recalibrate}.py`,
`figures/study_flowchart.dot` (duplicated-quota-math or superseded old-tahoe passes).

**Regenerated clean from the engine:** `docs/SUPPLEMENT.md` (constants of record),
`fixtures/tahoe_{constants,calibration,within_sigma}.json`, `fixtures/{orion,trade}_*_quota.csv`,
`fixtures/tahoe_{per_dose,quota_per_condition,applications}.json/.csv`. A guard
(`tests/test_no_duplicate_math.py::test_no_retired_thesis_constants_in_prose_or_code`) fails CI if any of
the fingerprints above reappear in a `.md`/`.py` outside the four exempt docs.
