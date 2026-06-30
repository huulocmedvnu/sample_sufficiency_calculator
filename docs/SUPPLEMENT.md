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
| [`docs/TECHNICAL_BLUEPRINT.md`](TECHNICAL_BLUEPRINT.md) | Audited 5-section blueprint (+ §0 corrections) | Master draft scaffold; the **authoritative** narrative + Auditor's corrections |
| [`docs/CASE_STUDIES.md`](CASE_STUDIES.md) | 5 real-drug sufficiency spectrum | Results: empirical case studies table + analysis paragraph |
| [`docs/INVARIANCE.md`](INVARIANCE.md) | Downsampling-invariance demo (OVER vs UNDER) | Results: "Downstream Functional Invariance" subsection |
| [`docs/MANUSCRIPT_DRAFT.md`](MANUSCRIPT_DRAFT.md) | Full assembled manuscript draft | Submission scaffold (Claude-drafted from the source-of-truth) |
| [`docs/MANUSCRIPT_DEEPSEEK.md`](MANUSCRIPT_DEEPSEEK.md) | Agent-team manuscript (DeepSeek writers + Claude audit) | Independent machine-drafted scaffold; cross-check against `MANUSCRIPT_DRAFT.md` |
| [`docs/manuscript.md`](manuscript.md) | **Compilable** manuscript (unified `$…$` math, pure-ASCII prose) | `pandoc docs/manuscript.md -o manuscript.pdf --pdf-engine=xelatex`; regenerate via `scripts/normalize_manuscript.py` |
| [`../agents/`](../agents/) | Manuscript agent team (`team.py`, `ds_client.py`) | Reproducible fact-gated drafting pipeline; see `agents/README.md` |
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
| Embedding dimension d | 50 | shared PCA |
| **Tahoe** marginal σ² | 7.66 | `tahoe_calibration.json`, `calibrate.py` |
| **Tahoe** baseline N₀ | 1,394 cells/well (median) | reference atlas, 15,200 groups |
| Standard tolerance θ★ | 0.1 rad (5.73°) | case studies / resource demo |
| Tahoe quota law @ standard config | n★ = 75,068 / m² | derived |
| Over/under boundary | m = 7.34 (n★ = N₀) | derived |
| **n★ distribution** (292 drugs, θ=0.1) | median 23,934; only 2% over-sampled; 64% need 10k–50k; 17% ghost >50k | `SCALE_AUDIT.md` |
| 100×3×2 screen budget | flat 39.2M vs adaptive(cap10k) 5.7M vs θ=0.2 3.8M cells | `SCALE_AUDIT.md` |
| Depth-fixed resolution | θ(N₀) = 0.734 / m rad | derived |
| **EmeraldBay** within-condition σ² | ≈ 2.12 | `calibrate_emeraldbay.py` (≠ marginal 7.66; see caveats) |
| EmeraldBay cells streamed (5 shared cell lines) | 142,883 (of 58 GB / 116 shards) | companion pipeline |
| EmeraldBay wells gated | 101; predicted OVER = 15 → 100% met tol | gating, θ★ = 0.20 rad |
| Symbolic Jacobian residual | zero matrix; max float diff 2.8×10⁻¹⁷ | `verify_theory.py` L1 |
| Monte-Carlo check | K=50,000; seed=0; rel. err 0.13%; ρ²=3.82×10⁵; d_eff=12.2 | `verify_theory.py` L2 |
| Laurent–Massart tail coverage | empirical Pr(θ>θ★) = 0.34% ≤ δ = 10% | `test_calculator.py` |
| Test suite | 13 tests pass (12 unit + EmeraldBay integration) | CI |

**Held-out angular-error curves (EmeraldBay, prediction θ²(n)=tr(PΣP)/m²·(1/n−1/N)):**

| Group (actual identity) | N | m | mean rel. err | slope vs tr(PΣP)/m² | R² |
|---|---:|---:|---:|---:|---:|
| DMSO_T0 × HS-578T | 1067 | 12.41 | 1.4% | 0.769 vs 0.763 | 0.9979 |
| DMSO_T0 × AN3-CA | 752 | 6.39 | 0.8% | 1.379 vs 1.375 | 0.9996 |
| Encorafenib × HEC-1-A | 1117 | 3.02 | 2.4% | 6.66 vs 7.82 | 0.9954 |

**Case-study spectrum (Tahoe; n★ = 75,068/m²):** Homoharringtonine m=14.4 n★=365 OVER (3.8× / 74%·93%);
Idarubicin m=10.6 n★=664 OVER (2.1× / 52%·77%); Dinaciclib m=6.5 n★=1,760 UNDER; Resveratrol m=3.0
n★=8,518 UNDER; Ribociclib m=0.84 n★=107,547 UNDER (ghost). Full table → [`CASE_STUDIES.md`](CASE_STUDIES.md).

---

## 3. Honesty ledger (consolidated — reviewers will probe these)

1. **Variance-definition mismatch.** EmeraldBay σ²≈2.1 is *within-condition*; Tahoe 7.66 is *marginal* —
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

