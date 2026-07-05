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
| [`../scripts/applications/`](../scripts/applications/) | Applied analyses (reliability audit, cost, adaptive-vs-flat) | Results §2.7; `tahoe_applications.json`, `tahoe_moa_recovery.json` |
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
| **Tahoe** within-condition σ² (headline) | **0.9567** (residual after removing each condition's mean; the noise that enters a centroid) over all 95,624,334 cells | fresh recompute; `tahoe_within_sigma.json`, `pass5_within_sigma.py` |
| **Tahoe** baseline N₀ | **1,296 cells / (drug×dose×line) condition** (median, post-filter) | `obs_metadata`; `tahoe_condition_counts.csv` |
| Standard tolerance θ★ | 0.1 rad (5.73°) | case studies / resource demo |
| Tahoe quota law @ standard config | **n★ = 9,376 / m²** | derived |
| Over/under boundary | **m = 2.69** (n★ = N₀ at median depth) | derived |
| **Regime split** (379 drugs × 3 doses × 50 lines = 56,827 conditions, θ=0.1) | **10.6% OVER · 89.0% UNDER · 0.4% Ghost** (89.4% under-or-ghost); median n★ 5,794, median m 1.27 | `tahoe_per_cell_line.csv`, `SCALE_AUDIT.md` |
| Marginal σ² (conservative bound) | **σ²=2.406** (all cells; additionally absorbs between-condition/line structure) → n★=23,577/m²; spectrum **2.5% OVER · 89.3% UNDER · 8.2% Ghost**, median n★ 14,570, boundary m 4.27 | `tahoe_calibration.json`, `scripts/tahoe_recompute/` |
| **Per-drug spectrum** (of 150 conditions/drug) | 28/379 drugs OVER in 0 conditions; strongest Panobinostat 120/150, Homoharringtonine 100/150 | `tahoe_per_drug.csv` |
| Depth-fixed resolution | θ(N₀) = **0.269 / m** rad (at median N₀) | derived |
| **Study-design layout** | **100,648,790 cells** (95,624,334 pass `full`, 5.0% filtered); **1,344 wells = 14 plates × 96**; cells/well median 71,092 (pre) / 67,212 (post) | `obs_metadata`; `tahoe_layout_summary.json` |
| Plate QC variation | plate3 11.63% filter loss (operationalizes "excl3") vs 3–5% typical | `per_plate` in `tahoe_layout_summary.json` |
| **EmeraldBay** within-condition σ² | **≈ 0.938** (full atlas, 52 lines, 1.83M cells; 5-shared-line slice = 0.848) | `scripts/emeraldbay_recompute/` (own **frozen** HVG(2000)+PCA(50); ≈ Tahoe within-condition 0.9567 to ~2%) |
| EmeraldBay basis reproducibility | own 50-d subspace reproduces across disjoint-shard refits at mean cos² **0.79 → 0.96** once the HVG set is frozen (per-fit HVG churn removed); Tahoe 0.93 → 0.995 | `scripts/check_basis_stability.py`, `confirm_frozen_hvg.py`; `fixtures/emeraldbay_frozen_hvg.json` |
| EmeraldBay cells streamed (5 shared cell lines) | 141,720 (of 1.83M; 57.7 GB / 116 shards) | `scripts/emeraldbay_recompute/` |
| EmeraldBay gating (full atlas) | **3,971 groups, 52 lines** (≥100 cells): 347 OVER (8.7%), 91.3% UNDER; downsample-verified **347/347 → 100% over all 52 lines** (`pass5_gating_full.py`) | gating, θ★ = 0.20 rad |
| Tahoe gating (full atlas, direct) | downsample-verified **5,503/5,503 → 100%** of predicted-OVER conditions over all 44 lines with any OVER condition (per-line-mean baseline; `pass4c_gating.py`) | θ★ = 0.1 rad; ≠ DMSO-referenced 10.6% spectrum |
| Symbolic Jacobian residual | zero matrix; max float diff 2.8×10⁻¹⁷ | `verify_theory.py` L1 |
| Monte-Carlo check | K=50,000; seed=0; rel. err 0.13%; ρ²=3.82×10⁵; d_eff=12.2 | `verify_theory.py` L2 |
| Laurent–Massart tail coverage | empirical Pr(θ>θ★) = 0.34% ≤ δ = 10% | `test_calculator.py` |
| Test suite | 13 tests pass (12 unit + EmeraldBay integration) | CI |
| **Application: reliability audit** | conditions resolved 10.6% / 42.4% / 67.1% at θ★=0.1/0.2/0.3 rad; similarity-graph edges with both endpoints resolved 10.0% / 39.2% / 65.0% | `tahoe_applications.json`, `scripts/applications/` |
| **Application: cost** (600-condition screen) | uniform-safe 10.9M cells (~$3.27M) vs quota-guided 3.5M (~$1.06M) → 3.1×, ~$2.21M saved @ $0.30/cell; 10.6% over-sampled, median multiplex 1.6× | `tahoe_applications.json` |
| **Application: adaptive vs flat** | budget-matched — indistinguishable in this under-sampled regime (k-NN graph Jaccard within 0.005) | `tahoe_moa_recovery.json` |

**Held-out angular-error curves (EmeraldBay, from-raw recompute `scripts/emeraldbay_recompute/`;
prediction θ²(n)=tr(PΣP)/m²·(1/n−1/N); embedding + σ² fit on the FULL atlas (52 lines); the four
validation groups and the 132 gated groups are restricted to the 5 lines shared with Tahoe-100M):**

| Group (actual identity) | N | m | mean rel. err | slope vs tr(PΣP)/m² | R² |
|---|---:|---:|---:|---:|---:|
| DMSO_T0 × HS-578T | 1067 | 5.89 | 0.8% | 0.906 vs 0.918 | 0.9999 |
| DMSO_T0 × HEC-1-A | 435 | 3.58 | 0.6% | 1.831 vs 1.845 | 0.9999 |
| DMSO_T0 × BT-474 | 461 | 3.35 | 0.8% | 1.709 vs 1.780 | 0.9977 |
| Irinotecan × HS-578T (drug) | 563 | 1.96 | 5.3% | 10.08 vs 12.81 (low-SNR) | 0.9786 |

**Case-study spectrum (Tahoe; n★ = 9,376/m², per-drug median m and OVER out of 150 conditions):**
Panobinostat m=4.24 OVER in 120/150; Homoharringtonine m=5.88 OVER in 100/150; Harringtonine m=3.67 OVER
in 78/150; Palbociclib m=1.41 UNDER in 146/150 (OVER 4); Crizotinib m=0.99 Ghost in 3/150 (OVER 1). The unit
is the **(drug × dose × line) condition**, modulated by BOTH dose and line: homoharringtonine at 5 µM
reaches m=15.5 (n★=39, OVER) in NCI-H460 but only m=3.1 (n★=951, UNDER) in NCI-H661; the same drug in
NCI-H460 rises from n★=88 (0.05 µM) to n★=39 (5 µM). Full tables → [`CASE_STUDIES.md`](CASE_STUDIES.md),
`tahoe_per_drug.csv`, `tahoe_per_dose.csv`, `tahoe_quota_per_condition.csv`.

---

## 2b. Constants of record — genetic-perturbation scale-up (X-Atlas/Orion, CRISPRi Perturb-seq)

Cross-modality extension: the same calibrated framework applied to genome-wide gene knockdown
([`ORION_GENE_PERTURBATION.md`](ORION_GENE_PERTURBATION.md); `scripts/orion_recompute/`; fixtures
`orion_{HCT116,HEK293T}_{quota.csv,summary.json}`). Condition = (target gene × line); control arm = pooled
**Non-Targeting (NTC)** centroid; m = ‖μ_knockdown − μ_NTC‖; quota theorem unchanged. Full-atlas streaming
recompute (all 17.4B HCT116 + 29.1B HEK293T expression nonzeros re-projected; PCA projection validated to 1e-4
vs scanpy; `total_counts` == raw row sum exactly). θ★ = 0.1 rad, d = 50.

| quantity | HCT116 | HEK293T | source |
|---|---|---|---|
| cells (post-QC) / knockdowns scored | 3,404,169 / 17,585 | 4,534,299 / 17,856 | `orion_<LINE>_summary.json` |
| NTC control cells | 165,562 | 218,838 | " |
| **within-condition σ²** (headline plug-in) | **0.9133** | **1.0334** | " (≈ Tahoe 0.9567 / EmeraldBay 0.938 across modality) |
| marginal σ² | 0.9167 | 1.0396 | " (marginal≈within ⇒ most knockdowns transcriptionally weak) |
| **quota law** | **n★ = 8,950/m²** | **n★ = 10,127/m²** | derived |
| median magnitude m | 0.141 | 0.351 | " (vs Tahoe 1.27; knockdowns ~4–10× weaker) |
| **regime spectrum (OVER/UNDER/Ghost)** | **0.0% / 28.3% / 71.7%** | **0.0% / 49.4% / 50.6%** | " |
| median n★ / median acquired N₀ | 409,983 / 155 | 52,908 / 204 | " |
| **median deficit n★/N₀** | **~2,260×** | **~248×** | " (Tahoe 4.5×) |
| detectable fraction (signal>1.5×floor); OVER among them | 14.7%; **0.0%** | 35.3%; **0.0%** | " |
| **aniso/iso ratio** | 0.925 | **0.660** | " (HEK293T = where the anisotropic form departs from isotropy; THEORY §8) |
| tolerance robustness (% OVER @ θ=0.2 / 0.3) | 0.0% / 0.4% | 0.2% / 2.1% | derived |
| strongest-m knockdowns (validation) | EIF2S3, RPL31/18/7/30, MED22, NUP93, MAK16, UTP15, NOL8, EXOSC9, HEATR1 | RPL38, NUP93, POLR2C/I, RPS20, SF3B5, CPSF4, SNIP1 | `orion_<LINE>_quota.csv` (essential machinery lands at top) |

**Headline:** no gene knockdown is over-sampled for direction estimation at θ★=0.1 on either line; ~50–72% are
"ghosts." Same theorem, same σ² (~1.0) and n★ (~9,000/m²) law as chemical — but knockdowns' small magnitudes
push `n★ ∝ 1/m²` up so the atlas is near-universally direction-under-sampled. KD efficiency (HCT116 75.4% /
HEK293T 51.5%) does not drive magnitude here — HEK293T (lower KD) has the *larger* median m.

---

## 2c. Constants of record — essential-gene extension (TRADE, GEO GSE264667), **validated** (manuscript Phase F)

The genetic counterpart of Phase C: because the TRADE essential-gene CRISPRi screens are ~2.5e5 cells/line (not
1e8), per-cell PCA coordinates are retained and the parameter-free downsample-and-measure falsification runs
directly. Dual-sgRNA CRISPRi, 2,393 DepMap common-essential genes/line, pooled Non-Targeting control. QC per
TRADE: low-UMI < frac×median (Jurkat 0.14 / HepG2 0.18), high-mito mito_UMI>1750/3000 (mito_UMI =
mitopercent×UMI_count), single-or-dual-same-gene guides. Full 50×50 within-condition Σ; θ★=0.1, d=50.
Pipeline `scripts/trade_recompute/`; docs `TRADE_EXTENSION_PLAN.md`; fixtures `trade_{jurkat,hepg2}_{quota.csv,
summary.json,falsification.json}`.

| quantity | Jurkat-Essential | HepG2-Essential | source |
|---|---|---|---|
| cells (raw → QC-pass) | 262,956 → 252,896 (96.2%) | 145,473 → 129,772 (89.2%) | `trade_<line>_summary.json` |
| knockdowns scored / NTC cells | 2,187 / 11,514 | 1,883 / 4,380 | " |
| **within-condition σ²** | **1.498** | **1.867** | " (higher than ~1.0 band; strong perturbations flatten PCA spectrum) |
| quota law (equal-arm / large-NTC-pool) | 14,676 / **7,338** /m² | 18,294 / **9,147** /m² | " (large-pool ≈ 9,000 band) |
| median magnitude m | 1.84 | 2.25 | " (essential genes STRONG; > Tahoe 1.27, ≫ Orion 0.14) |
| **spectrum OVER/UNDER/Ghost** | **0.2 / 85.5 / 14.3%** | **0.6 / 81.1 / 18.2%** | " |
| median n★ / median depth / deficit | 4,113 / 85 / 48× | 3,346 / 48 / 70× | " (**depth**-limited, not magnitude-limited) |
| detectable (signal>1.5×floor); OVER among | 71%; 0% | 61%; 1% | " |
| aniso/iso ratio | 0.964 | 0.931 | " |
| **Phase-C falsification, strong m>4** | **slope ratio 1.02, R² 0.996** | **0.96, R² 0.996** | `trade_<line>_falsification.json` (parameter-free) |
| in-regime (ρ²≥3) count; breakdown Spearman(1−ratio,1/ρ²) | 547; 0.92 | 120; 0.79 | " (predicted low-SNR breakdown reproduced) |
| strongest-m knockdowns (validation) | EXOSC2–9 (RNA exosome), SLU7, BUD13, MED21/27 | EIF2S1, HSPA9, PSMB5/PSMC5, RPL* | `trade_<line>_quota.csv` |

**Headline:** first end-to-end **validated** gene-perturbation extension — the parameter-free angular-error slope
matches to 0.96–1.02 (R² 0.996) on strong essential-gene knockdowns, with the predicted low-SNR breakdown. New
regime: perturbations strong but shallowly sampled (median 48–85 cells) → depth-limited, yet still ~0% over-sampled.

---

## 3. Honesty ledger (consolidated — reviewers will probe these)

1. **Variance definition.** The headline now uses the *within-condition* residual σ² on both atlases — Tahoe 0.9567 and EmeraldBay ≈0.938 (5-shared-line validation slice = 0.85) — a like-for-like match to ~2%. The coarser *marginal* Tahoe value (2.41, all cells) is retained only as a conservative bound; it is not directly comparable to a within-condition estimate. Both establish σ² as a platform/pipeline-specific plug-in.
2. **DMSO_T0 group identity.** Three of the four validation groups are DMSO_T0 time-zero reference
   populations (HS-578T, HEC-1-A, BT-474), not drug effects; only Irinotecan×HS-578T is a drug. The
   held-out test validates geometry, not biology.
3. **Low-SNR breakdown is expected, not hidden.** Irinotecan (m=1.96) shows a ~21% slope deficit at small
   n — the first-order Delta breakdown when ρ²=m²/(uᵀSu) is not ≫1 (THEORY §5).
4. **EmeraldBay HVG is frozen for reproducibility.** The EmeraldBay embedding uses a frozen 2000-HVG set
   (`fixtures/emeraldbay_frozen_hvg.json`) rather than a per-fit selection: per-fit HVG churn dropped the
   50-d subspace agreement between disjoint-shard refits to cos²≈0.79, and freezing the HVG set restores it
   to ≈0.96 (Tahoe's basis was already stable at 0.93→0.995 and is left per-fit). The validation conclusions
   are unchanged by the switch; σ², the held-out slope match, and the 100% gating guarantee all persist.
4. **Gating tolerance = 0.20 rad** in the EmeraldBay experiment (not 0.1).
5. **Resource figures (~156× multiplex, 99%/~100% compute) are Tahoe-derived** (over-sampled exemplar
   homoharringtonine 5 µM × NCI-H460, N₀=6,060, θ=0.1), not EmeraldBay.
6. **Paclitaxel unavailable** in the batch-clean array (plate3-excluded) → Homoharringtonine is the
   cytotoxic exemplar. Magnitudes are 24 h survivor transcriptional norms (not viability).
7. **Scope limit.** n★ governs centroid/direction (pseudobulk) only — **not** cell-level UMAP local
   structure or rare-population detection; compute gains are polynomial, never exponential.
8. **Resveratrol MoA.** metadata `moa-fine`="unclear"; it is a moderate-magnitude exemplar only — no
   mechanistic/MoA claim is attached.
9. **Orion is an application of the proven framework, not an independent re-proof.** The geometry is a
   theorem (Lean §10); Orion tests whether the σ² *calibration* transfers (it does, ~1.0 on both lines) and
   what the quota says at genome scale. A dedicated Orion downsample-and-measure falsification (the
   parameter-free slope test of Phase C) is the next step and is **not yet run**.
10. **Orion weak-majority magnitudes are order-of-magnitude.** With median ~155–204 acquired cells/knockdown,
    the per-condition floor on m is non-negligible; m is bias-corrected (`−tr(Σ)(1/n+1/n_NTC)`), but for the
    weak majority n★ means "unresolvable at this depth," not an exact count. ~15% (HCT116) / 35% (HEK293T)
    clear the floor; even among those, 0% are over-sampled. Magnitudes are survivor transcriptional norms.
11. **TRADE σ² is higher (1.5–1.9) and the transfer is looser.** Strong, diverse essential-gene perturbations
    flatten the PCA spectrum (top-PC EVR 0.04–0.07 vs Orion 0.13–0.20), raising per-PC variance; the large-NTC-
    pool quota constant (7.3k–9.1k/m²) still lands near ~9,000, but σ² is a platform plug-in and does not claim
    a tight cross-dataset match here. TRADE is depth-limited (median 48–85 cells), not magnitude-limited.
12. **TRADE QC interpretation.** Low-UMI read as `UMI < frac×median` (Jurkat 0.14/HepG2 0.18; in the deposited
    h5ad this filter is near-vacuous — mito is the binding filter); mito threshold stored as `mitopercent`, so
    absolute mito UMIs reconstructed as `mitopercent×UMI_count` (>1750/3000). Baseline is the Non-Targeting
    centroid (both quota and falsification), so m-values are NTC-referenced.

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

