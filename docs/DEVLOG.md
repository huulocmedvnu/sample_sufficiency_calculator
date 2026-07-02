# Manuscript development log

Chronological record of how the manuscript and its tooling were built. Newest phase last. Every
quantitative claim is gated against the Constants of Record ([`SUPPLEMENT.md`](SUPPLEMENT.md)); every
reference is Crossref-verified; no language model was trusted as a source of fact.

## Pipeline (current)

```
docs/MANUSCRIPT_DEEPSEEK.md      (source: prose + \(..\)/\[..\] math + [@cite] keys, BiB section order)
        │  scripts/normalize_manuscript.py   (unify math to $..$, ASCII prose, fix pandoc edge cases)
        ▼
docs/manuscript.md               (compilable; YAML carries bibliography + Vancouver CSL)
        │  pandoc --citeproc      (csl/vancouver.csl + references.bib)
        ├─► manuscript.pdf        (pandoc + tectonic; 14 pp)
        └─► manuscript.docx       (native OMML equations + Word tables; editable)
```
Rebuild from repo root: `pandoc docs/manuscript.md -o manuscript.pdf --citeproc --pdf-engine=xelatex`
and `pandoc docs/manuscript.md -o manuscript.docx --citeproc`. Reference check: `python
scripts/verify_references.py` (32/32 DOIs).

## Phases

1. **Agent team** (`a23603e`). Built `agents/` — a fact-gated drafting pipeline: a deterministic
   SourcePacker injects the Constants of Record, DeepSeek `deepseek-v4-pro` writes each section, an
   adversarial Fact-Checker labels every claim SUPPORTED/UNSUPPORTED/DRIFTED against the pack, a Reviser
   repairs flagged claims, an Editor stitches. The run auto-caught real drifts (a wrong MSE formula, a
   per-well-vs-per-condition median); one section whose JSON came back empty was escalated to manual
   (Claude) audit. Output: `docs/MANUSCRIPT_DEEPSEEK.md`.

2. **Verified bibliography** (`11e4757`, later `d102493`). Every reference checked against the Crossref
   REST API (DOI + first author + year) by `scripts/verify_references.py` — never generated from model
   recall. Tahoe-100M (bioRxiv) and EmeraldBay (Hugging Face) verified against their primary records.

3. **Compilable + typeset** (`bd9ee41`, `d0adb91`, `9b07356`, `a7a7edc`). `normalize_manuscript.py`
   strips the editor's outer code fence, unifies three competing math notations to `$…$`/`$$…$$`, and
   removes breaking Unicode. Real compiles with pandoc + tectonic; rendered pages inspected and math-
   rendering bugs fixed (`$…$` adjacent to a digit, trailing space before a closing `$`, accented author
   names, `{,}` leaking into prose).

4. **Briefings in Bioinformatics format + 34 refs** (`d102493`). Reordered to Abstract → Introduction →
   Materials and methods → Results → Discussion → Conflicts/Funding/Data availability → References.
   Switched to pandoc `[@key]` + citeproc + a Vancouver CSL so citations auto-number in order of
   appearance and the reference list is generated; expanded 17 → **34** verified references (Connectivity
   Map/L1000, Perturb-seq variants, scPerturb, demuxlet/DoubletFinder, sctransform/Seurat/best-practices,
   scVI/CPA/GEARS, Johnstone, Vershynin, Mardia–Jupp).

5. **De-anonymization** (`f16a3a2`, `8afa36f`). Reversed the earlier generalization at the source: removed
   the `ALIAS` map from `agents/team.py` and the `REF_BY_*` relabelling in `scripts/dataset_inventory.py`
   (the actual cause the agent team mirrored), then propagated real names repo-wide — **Tahoe-100M**,
   **EmeraldBay**, and the cell lines **HS-578T, AN3-CA, HEC-1-A, BT-474, C-33 A**.

6. **Prose rewrites** (`f16a3a2`, `aa7c335`, `1d647b2`). Abstract + Introduction, then Results/Methods/
   Discussion, rewritten in neutral scientific register and made more application-oriented (a Methods
   "how to apply" recipe; a Discussion three-stage workflow; Results read as planning rules). The
   Introduction was then trimmed to motivation only — context → gap → high-level aim — with the
   derivation, the closed-form equation, and all result numbers moved out to Results/Methods.

7. **Editable Word output** (`8c375e3`). `manuscript.docx` built via pandoc; all math is **native Office
   Math (OMML)** equation objects (editable in Word, not images), plus native Word tables and the
   numbered reference list. Fixed a normalizer bug that double-escaped `\arccos` (LaTeX tolerated it; the
   docx math parser did not).

8. **Key Points + Limitations** (`646f1c1`). Added a Briefings in Bioinformatics-style **Key Points** box
   (5 high-level bullets) between the Abstract and Introduction. Rewrote the Limitations in the neutral
   register and removed the redundant inline "Limitations." paragraph the Discussion still carried,
   folding its unique points (Gaussian-centroid assumption; additive vs multiplicative/rotational batch)
   into the list — now 10 numbered items.

9. **Proof revision + Lean formalization.** Audited the main derivation (`docs/THEORY.md`). Found and
   fixed a real error in the §5 second-order expansion: the boxed formula was the expansion of
   $\mathbb E[\tan^2\theta]$, not $\mathbb E[\theta^2]$ — it dropped the arc-vs-tangent term
   $-\tfrac23\,\mathbb E\lVert Pe\rVert^4/(m^2\operatorname{tr}(PSP))$, which enters at the *same*
   second order. As written it made the prediction *worse* than leading order (Monte-Carlo: leading
   +0.9%, §5-as-written +3.4% — it matched $\mathbb E[\tan^2\theta]$ to 0.07% — corrected eq. (6)
   +0.06%). The leading-order law, the quota, and the $O(\rho^{-2})$ validity conclusion are all
   unaffected. Locked the fix in `tests/verify_theory.py` (new Layer 4). Then formally verified the
   **deterministic core** in Lean 4 / Mathlib (`lean/`): Lemma 1 (`Dg(v)=(1/m)P`), the key trace
   identity `tr(PΣP)=trΣ−uᵀΣu`, and its isotropic/diagonal reductions — no `sorry`, standard axioms
   only. The probabilistic core (E[θ²], tail) stays empirical (`verify_theory.py`). See
   `docs/THEORY.md` §10 and `lean/README.md`.

10. **Manuscript: Lean verification + rebuild.** Propagated the Lean formalization into the manuscript
    source (`MANUSCRIPT_DEEPSEEK.md`): Abstract and Key Points note the deterministic core is
    machine-checked in Lean 4/Mathlib; Methods "Verification and CI" gained a formal-proof bullet; and
    Results "Three-way verification" → **"Four-way verification"** with a new mode 4 (Lemma 1, the
    trace identity, iso/diagonal reductions; `no sorry`, standard axioms only) plus a scope note that
    the probabilistic E[θ²]/tail results stay empirical. Extended `scripts/normalize_manuscript.py`
    `CODE_RE` so Lean code tokens (`sorry`, `#print axioms`, `propext`, `Classical.choice`,
    `Quot.sound`, `.lean`, `Mathlib`) stay monospace instead of being turned into `$…$` math (which
    also avoided `#` breaking LaTeX). Regenerated `manuscript.md` and rebuilt `manuscript.pdf`
    (pandoc + tectonic, 14 pp) and `manuscript.docx`; verified the Lean text renders with no
    code-as-math leak. No factual correction was needed — the manuscript never carried the §5 formula.

11. **Introduction rewrite (motivation, symbol-free).** Replaced the Introduction's formalism-heavy
    prose with a four-paragraph narrative — (1) single-cell screens and MoA-by-direction, (2) angular
    error as the binding constraint of a direction-based screen, (3) why isotropic power calculations
    fail on strongly anisotropic single-cell noise, (4) the closed-form anisotropic quota as a
    budget-optimized triage instrument. Removed all explicit math (\(\hat\mu\), \(\hat v\),
    \(\mathbb R^d\), \(\arccos\), \(P=I-uu^\top\), \(\theta_\star\)); the formalism now lives only in
    Methods. Regenerated `manuscript.md`, rebuilt PDF + DOCX.

12. **From-raw, dose-resolved Tahoe-100M recompute.** Discarded the old cached calibration and
    reanalysed from the raw 100,648,790-cell dataset (`scripts/tahoe_recompute/`, following the theislab
    vevo_100m recipe: `normalize_total(1e4)`→`log1p`→HVG(2000)→PCA(50)), streamed over all 3,388
    expression shards (337 GB) with a resumable, parallel-prefetch pipeline. **Caught and fixed a
    dose-pooling error** mid-way: the first pass keyed pseudobulk on (drug×line), silently averaging
    over the ~3 doses per drug and summing their cells; re-streamed keyed on (sample×line) and pooled to
    the correct **(drug × dose × line)** unit with plate-matched DMSO controls. New constants:
    **σ²=2.406** (was 7.66), **n★=23,577/m²** (was 75,068), **N₀=1,296** cells/condition (was 1,394);
    panel **56,827 conditions** with **2.5% OVER / 89.3% UNDER / 8.2% Ghost** at θ=0.1. Also derived the
    study-design layout from `obs_metadata` (1,344 wells = 14 plates × 96; 5.0% QC-filtered; plate3
    11.6% loss). Propagated repo-wide: `calculator.py`/`calibrate.py`/tests, `fixtures/tahoe_*`,
    `SUPPLEMENT.md`, the manuscript (§2.4 rewritten to the drug×dose×line spectrum; Tables 2–4; rebuilt
    PDF/DOCX), README; scaffold docs (`CASE_STUDIES`, `SCALE_AUDIT`, `TECHNICAL_BLUEPRINT`,
    `MANUSCRIPT_DRAFT`) carry a superseded-numbers banner and defer to `fixtures/tahoe_*.csv`.

13. **From-raw EmeraldBay recompute (independent validation).** Reprocessed the second atlas
    (tahoebio/EmeraldBay, 1.83M cells, 116 shards, 57.7 GB) from the raw counts with the same recipe
    (`scripts/emeraldbay_recompute/`, parallel-prefetch stream). **Scope (exact):** the HVG(2000)+PCA(50)
    embedding and the within-condition σ² ≈ 0.963 are computed over the **full atlas** (all 52 lines,
    1,831,648 cells, 4,912 condition-line groups); the held-out per-cell validation is
    restricted to the 5 shared lines (141,720 cells; 5-line σ²=0.896). Regime gating classifies all 52 lines from the per-condition moments (3,971 groups, 217 OVER = 5.5%); the downsample-and-measure accuracy check (10/10 met) is on the shared 5. **Held-out angular-error curves** (subsample n cells,
    realized RMS angle vs the closed form θ²(n)=tr(PΣP)/m²·(1/n−1/N)): all four (condition×line) groups
    give R² > 0.996; the strong DMSO_T0×HS-578T (m=6.46) matches the predicted slope to 0.2%, while the
    weak Encorafenib×HEC-1-A drug group (m=1.78) shows the expected ~13% low-SNR slope deficit.
    Regime gating at θ=0.20: 10/132 groups predicted OVER, **100% met tolerance** when downsampled to
    n★. Regenerated `fixtures/emeraldbay_calibration.json`; the integration test passes on it. Updated
    SUPPLEMENT, manuscript §2.3 + Methods, README. The formula is a theorem (Lean/§10), so the
    independent atlas tests the CLT/Gaussian-centroid *assumptions*, which hold to <1.2% on real cells.

14. **Applied-value analyses (Results §2.7).** Tested whether acting on n★ improves real outputs, from committed data (`scripts/applications/`). (i) Reliability audit: at θ★=0.1 only 2.5% of conditions and 2.1% of similarity-graph k-NN edges are resolved (18.5/17.7% at 0.2; 39.0/36.1% at 0.3) — most pairwise MoA calls rest on under-powered directions. (ii) Cost: a 600-condition screen needs 27.4M cells under uniform-safe loading vs 4.9M quota-guided (5.6×, ~$6.7M @ $0.30/cell). (iii) Honest null: budget-matched adaptive ≈ flat for graph recovery in this under-sampled regime (Jaccard within 0.005) — the value is triage and avoiding over-provisioning, not reallocation. Note: fine-grained MoA labels are NOT recoverable from single-condition directions here (full-depth k-NN ≈ chance), so no such claim is made. Fixtures `tahoe_applications.json`, `tahoe_moa_recovery.json`.

15. **Theory review + two corrections (independent re-derivation).** Re-derived the probabilistic core
    from scratch, independently of the Lean proof (which only certifies the *deterministic* algebra), and
    re-ran `tests/verify_theory.py` (all layers pass; Monte-Carlo E[θ²] to 0.13%, eq.(6) to 0.05%). The
    headline results are correct, but two statements were fixed. **(i) Confidence-quota caveat.** Added an
    explicit note (THEORY.md §4; manuscript Methods) that Laurent–Massart rigorously bounds the
    perpendicular-noise energy ‖Pe‖², and only becomes a statement about the angle via the first-order
    step θ≈‖Pe‖/m; the exact geometry tanθ=‖Pe‖/(m+uᵀe) also depends on the along-signal fluctuation in the
    denominator, so n★_δ is an *approximate* (1−δ) guarantee (exact as ρ→∞), conservative in practice
    (empirical coverage 0.34% vs nominal 10%). **(ii) Plug-in robustness claim was wrong.** The old text
    said the trace functional is "1-Lipschitz ⟹ direction error enters at second order" — a non-sequitur.
    Verified numerically that tr(PΣP) is *first-order* sensitive to direction error (δtr(PΣP)=−2δuᵀΣu+…,
    nonzero unless u is a Σ-eigenvector); what is second order is the *bias*, because the pilot direction
    is unbiased to first order and the mean-zero first-order term averages out. Reworded in THEORY.md §7 +
    manuscript Methods. Also relabelled the §5 second-order relative correction O(ρ⁻²)→O(θ★²) (the
    arc-vs-tangent term dominates the lever-arm term by ~d_eff). Rebuilt manuscript PDF/DOCX.

16. **Plain-math primer (`docs/THEORY_PRIMER.md` → `theory_primer.pdf`).** Wrote a self-contained
    companion that derives the same headline quota n★=2 tr(PΣP)/(m²θ★²) using only basic matrix algebra
    (transpose, trace, projection, eigenvalues) and one integral (E[z²]=1 for a standard normal). Path:
    triangle (only sideways noise tilts a direction, θ≈ℓ⊥/m) → projection P=I−uuᵀ picks out the sideways
    part → E[eᵀPe]=tr(PS) by averaging term-by-term → solve for n. Includes the isotropic reduction
    σ²(d−1), the Tahoe number 23,577/m², a worked table by signal strength, and an honest scope section.
    Hand-written (no citations); built with `pandoc --pdf-engine=tectonic` (5 pp).

17. **Readability rewrite of the manuscript (all sections).** Reworded every section of
    `docs/MANUSCRIPT_DEEPSEEK.md` for clarity while keeping the formal Briefings-in-Bioinformatics
    register: crisper finding-first topic sentences in Results, and each Methods derivation step now
    opens with the plain-language idea (e.g. "only noise perpendicular to the signal can rotate the
    direction") before the formal expression, mirroring the new `THEORY_PRIMER.md`. **Content was frozen:**
    an automated audit confirmed the 39 body `[@cite]` keys are identical, every table *value* is
    byte-identical (only U+2011→ASCII hyphens changed, which the normalizer does anyway), and no numeric
    result was altered. **One stale figure fixed:** §2.3 honesty-caveat 1 read "the two groups with m=12.4
    and m=6.4 are DMSO_T0" — contradicting Table 1, the Limitations, and the Constants of Record; corrected
    to the three DMSO_T0 reference groups at m=6.46/3.63/3.56 (drug = Encorafenib, m=1.78). Also unified
    "1,296"→"1{,}296" and de-duplicated one redundant 5-line name listing. Regenerated `docs/manuscript.md`
    (0 flagged chars, math parity even) and rebuilt PDF (16 pp) / DOCX; banner audit-ledger updated.

## Honesty ledger (carried in the manuscript banner + `SUPPLEMENT.md`)

- The two high-magnitude validation groups are DMSO time-zero populations, not drug effects; only one
  group (Encorafenib × HEC-1-A) is a drug — the held-out test validates geometry, not biology.
- EmeraldBay σ² ≈ 2.12 (within-condition) and Tahoe-100M σ² = 2.406 (marginal, fresh recompute) are not
  a clean platform head-to-head; both only establish σ² as a plug-in.
- Gating tolerance in the EmeraldBay experiment is 0.20 rad (not the standard 0.1).
- Resource figures are Tahoe-100M-derived; compute gains are polynomial, never exponential.
- Author lists beyond the third name are abbreviated and should be expanded before submission.
