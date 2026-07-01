# Supplement Index — Anisotropic Sample-Sufficiency Calculator
### Single entry point for manuscript drafting, review, and reproduction

This index cross-links every artifact behind *"Anisotropic Sample-Sufficiency for Single-Cell
Perturbation-Direction Screens"* and consolidates the **constants of record**. Repository:
`github.com/huulocmedvnu/sample_sufficiency_calculator` (frozen at commit recorded below). The external
calibrator lives in the companion data repository.

---

## 1. Document & code map

| Artifact | Role | Use it for |
|---|---|---|
| [`docs/THEORY.md`](THEORY.md) | Full anisotropic proof | Methods: derivation, distribution, tail bound, validity regime, plug-in estimation |
| [`../lean/`](../lean/) | Lean 4 / Mathlib formalization | Machine-checked deterministic core: Lemma 1 Jacobian + `tr(PΣP)=trΣ−uᵀΣu` + iso/diagonal reductions (no `sorry`, standard axioms; see `lean/README.md`, THEORY §10) |
| [`docs/TECHNICAL_BLUEPRINT.md`](TECHNICAL_BLUEPRINT.md) | Audited 5-section blueprint (+ §0 corrections) | Master draft scaffold; the **authoritative** narrative + Auditor's corrections |
| [`docs/CASE_STUDIES.md`](CASE_STUDIES.md) | 5 real-drug sufficiency spectrum | Results: empirical case studies table + analysis paragraph |
| [`docs/INVARIANCE.md`](INVARIANCE.md) | Downsampling-invariance demo (OVER vs UNDER) | Results: "Downstream Functional Invariance" subsection |
| [`docs/MANUSCRIPT_DRAFT.md`](MANUSCRIPT_DRAFT.md) | Full assembled manuscript draft | Submission scaffold (Claude-drafted from the source-of-truth) |
| [`docs/MANUSCRIPT_DEEPSEEK.md`](MANUSCRIPT_DEEPSEEK.md) | Agent-team manuscript (DeepSeek writers + Claude audit) | Independent machine-drafted scaffold; cross-check against `MANUSCRIPT_DRAFT.md` |
| [`docs/manuscript.md`](manuscript.md) | **Compilable** manuscript (unified `$…$` math, pure-ASCII prose) | `pandoc docs/manuscript.md -o manuscript.pdf --pdf-engine=xelatex`; regenerate via `scripts/normalize_manuscript.py` |
| [`../agents/`](../agents/) | Manuscript agent team (`team.py`, `ds_client.py`) | Reproducible fact-gated drafting pipeline; see `agents/README.md` |
| [`docs/DEVLOG.md`](DEVLOG.md) | Manuscript development log | Full build history + pipeline + rebuild commands (`../manuscript.pdf`, `../manuscript.docx`) |
| [`docs/REFERENCES.md`](REFERENCES.md) | Verified bibliography (34 refs) thematic ledger | Crossref-verified citations; `scripts/verify_references.py` re-checks all 32 DOIs |
| [`../references.bib`](../references.bib) | Canonical BibTeX (34 entries) | Cited as pandoc `[@key]`; rendered by citeproc + `csl/vancouver.csl` (Vancouver, BiB style) |
| [`docs/SCALE_AUDIT.md`](SCALE_AUDIT.md) | n★ is per-arm; project-budget ledger; Methods definition | Methods: "Cell Quota per Perturbation" def. + scale-budgeting |
| [`../src/calculator.py`](../src/calculator.py) | Estimator | `calculate_experimental_cell_quota` (isotropic), `calculate_cell_quota_anisotropic` (+ tail), `calculate_optimal_resource_allocation` (dual-sided), `rms_angular_error` |
| [`../src/calibrate.py`](../src/calibrate.py) | Tahoe calibration demo | Reproduce the σ²/N₀/quota numbers from cached arrays |
| [`../tests/verify_theory.py`](../tests/verify_theory.py) | Standalone formal+empirical harness | §3 verification (SymPy Jacobian + Monte-Carlo + implementation) |
| [`../tests/test_calculator.py`](../tests/test_calculator.py) | 12 unit/property tests | Scaling laws, reductions, tail coverage, validation |
| [`../tests/test_emeraldbay_integration.py`](../tests/test_emeraldbay_integration.py) | Gated external-validation test | Asserts the EmeraldBay result (CI-safe) |
| [`../fixtures/tahoe_calibration.json`](../fixtures/tahoe_calibration.json) | Tahoe ℓ_k + drug vectors | Self-contained calibration demo |
| [`../fixtures/emeraldbay_calibration.json`](../fixtures/emeraldbay_calibration.json) | EmeraldBay held-out curves + gating | §4 external-validation numbers |
| [`../.github/workflows/ci.yml`](../.github/workflows/ci.yml) | CI (py3.10/3.11/3.12) | Continuous verification of all asserts |
| [`../README.md`](../README.md) | Overview + usage | Quick start; abstract-level summary |

Companion (companion data repository): `scripts/calibrate_emeraldbay.py` (streaming calibrator),
`outputs/emeraldbay_heldout_validation.png` (Fig.)(provenance retained in the companion repository).

---

## 2. Constants of record (single authoritative source — do not paraphrase)

**Core formula.** `n★ = 2·tr(PΣP)/(m²θ★²) = 2(d−1)σ²/(m²θ★²)` (isotropic special case); `P = I − uuᵀ`.
Tail-controlled: `n★_δ = (2/m²θ★²)[tr(PΣP) + 2‖PΣP‖_F√L + 2‖PΣP‖_op·L]`, `L=log(1/δ)` (Laurent–Massart).

| Symbol / quantity | Value | Source |
|---|---|---|
| Embedding dimension d | 50 | shared PCA(50) |
| **Tahoe** per-cell σ² | **2.406** | fresh recompute; `tahoe_calibration.json`, `scripts/tahoe_recompute/` |
| **Tahoe** baseline N₀ | **1,296 cells / (drug×dose×line) condition** (median, post-filter) | `obs_metadata`; `tahoe_condition_counts.csv` |
| Standard tolerance θ★ | 0.1 rad (5.73°) | case studies / resource demo |
| Tahoe quota law @ standard config | **n★ = 23,577 / m²** | derived |
| Over/under boundary | **m = 4.27** (n★ = N₀ at median depth) | derived |
| **Regime split** (379 drugs × 3 doses × 50 lines = 56,827 conditions, θ=0.1) | **2.5% OVER · 89.3% UNDER · 8.2% Ghost** (97.5% under-or-ghost); median n★ 14,570, median m 1.27 | `tahoe_per_cell_line.csv`, `SCALE_AUDIT.md` |
| **Per-drug spectrum** (of 150 conditions/drug) | 132/379 drugs OVER in 0 conditions; strongest Panobinostat 73/150, Homoharringtonine 72/150 | `tahoe_per_drug.csv` |
| Depth-fixed resolution | θ(N₀) = **0.427 / m** rad (at median N₀) | derived |
| **Study-design layout** | **100,648,790 cells** (95,624,334 pass `full`, 5.0% filtered); **1,344 wells = 14 plates × 96**; cells/well median 71,092 (pre) / 67,212 (post) | `obs_metadata`; `tahoe_layout_summary.json` |
| Plate QC variation | plate3 11.63% filter loss (operationalizes "excl3") vs 3–5% typical | `per_plate` in `tahoe_layout_summary.json` |
| **EmeraldBay** within-condition σ² | **≈ 0.963** (full atlas, 52 lines, 1.83M cells; 5-shared-line slice = 0.896) | `scripts/emeraldbay_recompute/` (own PCA(50); ≠ marginal Tahoe 2.406) |
| EmeraldBay cells streamed (5 shared cell lines) | 141,720 (of 1.83M; 57.7 GB / 116 shards) | `scripts/emeraldbay_recompute/` |
| EmeraldBay groups gated | 132 (**5 shared lines**, ≥100 cells); predicted OVER = 10 → 100% met tol | gating, θ★ = 0.20 rad |
| Symbolic Jacobian residual | zero matrix; max float diff 2.8×10⁻¹⁷ | `verify_theory.py` L1 |
| Monte-Carlo check | K=50,000; seed=0; rel. err 0.13%; ρ²=3.82×10⁵; d_eff=12.2 | `verify_theory.py` L2 |
| Laurent–Massart tail coverage | empirical Pr(θ>θ★) = 0.34% ≤ δ = 10% | `test_calculator.py` |
| Test suite | 13 tests pass (12 unit + EmeraldBay integration) | CI |

**Held-out angular-error curves (EmeraldBay, from-raw recompute `scripts/emeraldbay_recompute/`;
prediction θ²(n)=tr(PΣP)/m²·(1/n−1/N); embedding + σ² fit on the FULL atlas (52 lines); the four
validation groups and the 132 gated groups are restricted to the 5 lines shared with Tahoe-100M):**

| Group (actual identity) | N | m | mean rel. err | slope vs tr(PΣP)/m² | R² |
|---|---:|---:|---:|---:|---:|
| DMSO_T0 × HS-578T | 1067 | 6.46 | 0.6% | 0.789 vs 0.787 | 0.9999 |
| DMSO_T0 × HEC-1-A | 435 | 3.63 | 1.2% | 1.819 vs 1.854 | 0.9990 |
| DMSO_T0 × BT-474 | 461 | 3.56 | 1.2% | 1.561 vs 1.623 | 0.9998 |
| Encorafenib × HEC-1-A (drug) | 1117 | 1.78 | 2.9% | 9.67 vs 11.17 (low-SNR) | 0.9965 |

**Case-study spectrum (Tahoe; n★ = 23,577/m², per-drug median m and OVER out of 150 conditions):**
Panobinostat m=4.24 OVER in 73/150; Homoharringtonine m=5.88 OVER in 72/150; Harringtonine m=3.67 OVER
in 64/150; Palbociclib m=1.41 UNDER in 150/150; Crizotinib m=0.99 Ghost in 28/150 (OVER in 0). The unit
is the **(drug × dose × line) condition**, modulated by BOTH dose and line: homoharringtonine at 5 µM
reaches m=15.5 (n★=98, OVER) in NCI-H460 but only m=3.1 (n★=2,390, UNDER) in NCI-H661; the same drug in
NCI-H460 rises from n★=222 (0.05 µM) to n★=98 (5 µM). Full tables → [`CASE_STUDIES.md`](CASE_STUDIES.md),
`tahoe_per_drug.csv`, `tahoe_per_dose.csv`, `tahoe_quota_per_condition.csv`.

---

## 3. Honesty ledger (consolidated — reviewers will probe these)

1. **Variance-definition mismatch.** EmeraldBay σ²≈0.96 is *within-condition* over the full atlas (52 lines; the 5-shared-line validation slice = 0.90); Tahoe σ²≈2.41 is *marginal* —
   not a clean platform head-to-head. Both establish σ² as a platform/pipeline-specific plug-in.
2. **DMSO_T0 group identity.** The two high-m validation groups are time-zero reference populations, not
   drug effects; only Encorafenib×HEC-1-A is a drug. The held-out test validates geometry, not biology.
3. **Low-SNR breakdown is expected, not hidden.** Encorafenib (m=3.0) shows a ~15% slope deficit at small
   n — the first-order Delta breakdown when ρ²=m²/(uᵀSu) is not ≫1 (THEORY §5).
4. **Gating tolerance = 0.20 rad** in the EmeraldBay experiment (not 0.1).
5. **Resource figures (3.8×, 74%/93%) are Tahoe-derived** (N₀=1394, θ=0.1), not EmeraldBay.
6. **Paclitaxel unavailable** in the batch-clean array (plate3-excluded) → Homoharringtonine is the
   cytotoxic exemplar. Magnitudes are 24 h survivor transcriptional norms (not viability).
7. **Scope limit.** n★ governs centroid/direction (pseudobulk) only — **not** cell-level UMAP local
   structure or rare-population detection; compute gains are polynomial, never exponential.
8. **Resveratrol MoA.** metadata `moa-fine`="unclear"; it is a moderate-magnitude exemplar only — no
   mechanistic/MoA claim is attached.

---

## 4. Reproduce

```bash
pip install -r ../requirements-dev.txt
pytest                          # 13 asserts (unit + symbolic + Monte-Carlo + EmeraldBay)
python ../tests/verify_theory.py    # standalone formal+empirical harness
python ../src/calibrate.py           # Tahoe calibration demo (reads cached arrays if present)
```

External calibration (one-time, 58 GB, companion data repository): `python scripts/calibrate_emeraldbay.py`.

**Frozen at:** `sample_sufficiency_calculator` repository HEAD (see `git log`); companion the companion data repository

