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
    1,831,648 cells, 4,992 condition-line groups); the held-out per-cell validation is
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

18. **Readability treatment of the scaffold docs (`CASE_STUDIES.md`, `TECHNICAL_BLUEPRINT.md`).** Same
    pass as phase 17, applied to the two supplement docs. **CASE_STUDIES** was already at dose resolution,
    so it got prose-only clarity edits (broke the single-block Analysis into strong/weak-end paragraphs;
    finding-first topic sentences); an automated audit confirmed every table value and number unchanged.
    **TECHNICAL_BLUEPRINT was a pre-recompute scaffold** — its numbers were stale throughout (EmeraldBay
    σ²=2.12, validation table m=12.41/6.39/3.02 with the wrong 3-group set, 142,883 cells, "101 wells",
    "39 M/5.7 M/85-100 drugs", compliance ≤1.15·θ★). At the user's direction these were **reconciled to
    the Constants of Record**, every value sourced from `SUPPLEMENT.md` + fixtures: EmeraldBay σ²≈0.963
    (0.896 5-line); the validation table rebuilt to the current four groups (HS-578T 6.46, HEC-1-A 3.63,
    BT-474 3.56, Encorafenib 1.78) with the manuscript's slopes/R²; 141,720 cells; 3,971 gated groups
    (217 OVER, 10/10 downsample-verified on the 5 shared lines); ≤1.05·θ★; budget 27.4 M→4.9 M (5.6×);
    §0 auditor items, §4, §5, and appendix all updated. Also added the §4 THEORY-consistent tail-bound
    caveat. Verified: no stale token remains, all current values present and matching. (These docs are
    not part of the pandoc build; no PDF rebuild needed.)

19. **Reconciled the stale honesty-ledger lines in `SUPPLEMENT.md` (Constants of Record).** The main
    constants table was already current, but three prose items in the §3 honesty ledger predated the
    dose-resolved recompute: item 2 said "two high-m validation groups" (now three DMSO_T0 groups:
    HS-578T, HEC-1-A, BT-474; drug = Encorafenib); item 3 gave Encorafenib m=3.0 / ~15% deficit (now
    m=1.78 / ~13%); item 5 gave resource figures "3.8×, 74%/93%, N₀=1394" (now ~60× multiplex, 98%/>99%
    compute, over-sampled exemplar homoharringtonine 5 µM × NCI-H460, N₀=6,060). All sourced from the
    current constants; verified no stale token remains and the exemplar N₀=6,060 matches CASE_STUDIES and
    TECHNICAL_BLUEPRINT.

20. **References expanded 34 → 50 (all Crossref-verified) + logic/citation pass.** Added 16 real
    references to strengthen the manuscript's logical scaffolding, not as padding: pseudobulk-vs-per-cell
    aggregation and pseudoreplication (Squair 2021, Crowell 2020, Zimmerman 2021) to justify the centroid
    summary; perturbation-screen breadth (Jaitin 2016, Gasperini 2019, Frangieh 2021, McFarland 2020);
    perturbation-response modeling (Lotfollahi scGen 2019); directional statistics (Fisher 1953); the
    Hanson–Wright parent inequality behind Laurent–Massart (Rudelson–Vershynin 2013); embedding-fidelity
    critiques for the scope caveat (Kobak 2019, Chari 2023); batch integration (Haghverdi 2018, Korsunsky
    2019); replicate/power precedent (Schurch 2016); and the Human Cell Atlas context (Regev 2017).
    **Integrity:** every DOI was fetched from the registrar via content negotiation and Crossref-verified
    (author+year) by `scripts/verify_references.py` — 48/48 DOIs pass, none hand-entered. Rebuilt PDF (17 pp,
    50 refs rendered) / DOCX; audit confirmed zero numeric or table drift (edits are additive citations +
    clauses). Banner reference line updated 34→50.

21. **§2.3 validation methodology made explicit + a sampling-scheme fix.** Added a "Subsampling protocol
    (exact)" paragraph to Results §2.3 spelling out (i) the ground-truth direction
    \(u=(\hat\mu_N-\bar\mu_{\text{line}})/\lVert\cdot\rVert\) built from all \(N\) condition cells against
    the per-line-mean reference (distinct from the DMSO reference of §2.4), (ii) the angular error
    \(\theta=\arccos(\hat u_n^\top u)\) as the arccosine of the cosine similarity (clipped), RMS over
    \(R=300\), and (iii) the finite-population handling — without-replacement subsampling makes
    \(\operatorname{Var}(\hat\mu_n-\hat\mu_N)=\Sigma(1/n-1/N)\), the fit passes through the origin at
    \(n=N\), and the grid caps at \(n\le N/2\) to avoid the degenerate boundary. Prompted by this
    transparency pass, **corrected the Tahoe stage-B curve script** (`pass4b_curves.py`) from
    with-replacement (a bootstrap, variance \(\Sigma/n\), no \(-1/N\)) to without-replacement, matching
    EmeraldBay and the FPC prediction; the fix ships with the falsification bundle but the running
    orchestrator already reads the corrected script. Rebuilt PDF (17 pp) / DOCX; tables and all data
    numbers unchanged.

22. **Editorial restructure for a broad bioinformatics/genomics audience.** Reorganized the manuscript
    per four editorial directions while keeping every constant intact. (1) **Introduction** rewritten as a
    biology-first onboarding: what single-cell perturbation screens are and why they matter → the standard
    workflow (treat → embed → aggregate to pseudobulk centroids → read direction for MoA) → the unanswered
    question "how many cells to trust the direction?" → the anisotropic quota as the answer. (2)
    **Methods total consolidation**: moved all procedure out of Results/Appendix — the exact subsampling
    protocol (all-N reference, without-replacement, FPC (1/n−1/N)), the preprocessing pipeline
    (log1p/HVG2000/PCA50), full dataset layouts, and the 13-test CI + Lean scope now live in Methods; the
    standalone Results §2.2 "four-way verification" was folded into the Methods verification subsection (no
    unique number lost — its figures already appeared there). (3) **Results retold as a linear arc** —
    Phase A (the theory + geometric insight), Phase B (the Tahoe audit / sample-sufficiency spectrum),
    Phase C (EmeraldBay out-of-distribution falsification), Phase D (utility: downsampling-invariance +
    adaptive budgeting). Tables renumbered by new order (1 per-drug, 2 by-dose, 3 EmeraldBay, 4 invariance,
    5 budget) with all in-text refs updated. (4) **Tone** softened around the math (e.g. tangent-space
    projector → "filter out the noise aligned with the effect vector, which only changes its length").
    **Integrity:** automated set-difference audit confirms all 50 distinct [@cite] keys preserved, all
    table-row values identical (multisets equal), and no data number lost or added (only §2.2/2.3/2.6
    section labels dropped, replaced by Phase A–D). Rebuilt PDF (16 pp, 50 refs) / DOCX; 0 flagged chars.

23. **Discussion rewritten as a scholarly dialogue.** Replaced the inward-looking, manual-like opening
    with four framing paragraphs that position the work in the field: (i) *A design rule, not a detector* —
    what the quota is and its three-stage use; (ii) *Detection versus direction* — a constructive contrast
    with scPower [@schmid2021scpower], bulk replicate planning [@schurch2016replicates], and Hotelling
    detection [@hotelling1931], arguing that existing tools collapse or invert the covariance to solve
    *detection* power (Mahalanobis \(v^\top\Sigma^{-1}v\), rewarding low-noise directions) whereas the
    directional problem is governed by \(\operatorname{tr}(P\Sigma P)\) and *hurt* by perpendicular noise —
    hence an anisotropic quota is not just conservative but wrong-signed for some conditions under
    isotropic assumptions; (iii) *A design standard* — the single-cell counterpart of bulk power analysis,
    with Lean-checked core; (iv) *The ghost regime in perspective* — reframing the 97.5% under-sampled/
    ghost audit through the reproducibility lens (raw cell count has masked directional resolution; the
    tool shifts the question from "how many cells?" to "is the mechanism resolvable?"). Technical
    paragraphs (scope, plug-in, dual-use, budget) and the Limitations honesty-ledger kept verbatim; tone
    constructive throughout ("evolving the standard, not criticizing the past"). Audit: all 50 [@cite]
    keys and the full number set preserved, tables unchanged. Rebuilt PDF (17 pp) / DOCX.

24. **Full-atlas falsification of the angular-error law (both atlases, from raw).** Executed the
    parameter-free predicted-vs-realized slope test at full scale. **Tahoe-100M:** re-streamed all 3,388
    shards from HuggingFace (337 GB), projected all **95,624,334** cells through the same PCA(50) basis to
    a 20 GB disk memmap, and ran without-replacement subsample curves on **56,195** of 67,018 (sample×line)
    conditions — the first direct subsample-and-measure test on the primary atlas (previously only the
    derived spectrum had been applied to it). **EmeraldBay:** re-projected the full **1,831,648**-cell /
    52-line atlas and tested 1,064 groups. Results: (1) **σ² confirmed on every cell** — Tahoe marginal
    σ²=2.4158 vs 2.406 (14-shard fit 2.4058), EmeraldBay within-condition 0.9745 vs 0.963, both within
    rounding; the calibration fit on ~0.4% of data holds on the whole atlas. (2) **In-regime (ρ²≥3) the
    parameter-free slope matches** — Tahoe median ratio 0.942 (1,790 conditions, R²=0.999), EmeraldBay
    0.981 (33 groups, R²=0.999). (3) **The predicted ρ⁻² breakdown is reproduced** (Spearman 0.57 Tahoe,
    0.79 EmeraldBay); whole-population ratios (0.68, 0.37) are the expected low-SNR tail, consistent with
    the "most conditions under-sampled/ghost" finding. (4) Two-independent-halves √2 test 1.352 (EmeraldBay,
    642 groups). Writeup `docs/FALSIFICATION.md`; scripts `scripts/*/pass4*.py` (+ `pass2b_project_all.py`);
    fixtures `emeraldbay_falsification{,_full}.json`, `tahoe_direct_curves_full.json`; figures under
    `outputs/`. A with-replacement bug in the Tahoe stage-B script was corrected to without-replacement
    (matching the FPC prediction) before the run. Superseded prototype `pass4_direct_curves.py` removed.

25. **Manuscript Phase C: added the full-atlas population-scale confirmation.** Folded the phase-24
    falsification results into Results Phase C as a "Population-scale confirmation on both atlases"
    paragraph + a 4th honesty caveat: the parameter-free slope test run across the full populations (all
    95,624,334 Tahoe cells / 56,195 conditions; full 52-line EmeraldBay / 1,064 groups), in-regime
    (ρ²≥3) median ratio 0.94 (1,790 Tahoe) and 0.98 (33 EmeraldBay) at R²≈0.999, the ρ⁻² breakdown
    (Spearman 0.57/0.79), and the σ² reconfirmation on every cell (marginal 2.4158 vs 2.406; within-
    condition 0.9745 vs 0.963). All added numbers sourced from `fixtures/tahoe_direct_curves_full.json`
    and `fixtures/emeraldbay_falsification_full.json` and cross-checked; audit confirms no existing number,
    table, or citation changed. Rebuilt PDF (18 pp) / DOCX.

26. **Full-atlas gating verification (all 52 EmeraldBay lines).** The gating *classification* was already
    atlas-wide, but the downsample-and-measure *check* had been limited to the 5 Tahoe-shared lines (10
    groups) because only their per-cell coordinates were retained. With `pass2b`'s full re-projection
    (`all_cells.npz`, 52 lines, 1,831,648 cells) available, `pass5_gating_full.py` reruns the check over
    the whole atlas: it reproduces the published classification exactly (3,971 groups ≥100 cells, **217
    predicted OVER**) and verifies **217/217 (100%)** meet ≤1.05·θ★ when downsampled to n★ — no longer a
    subset. Removed the "restricted to five shared lines" / "10-group subset" language from the manuscript
    (Datasets + Phase C gating paragraph); the five shared lines now serve only the four showcased curves
    (Table 3) and Tahoe cross-comparison. Updated `fixtures/emeraldbay_calibration.json` and the Constants
    of Record. Audit: no manuscript number/table/citation changed (217 and 10 both already present).
    Rebuilt PDF (18 pp) / DOCX.

27. **Full-atlas gating verification for Tahoe (the primary-atlas counterpart of phase 26).** Tahoe had
    the full-atlas *classification* (Phase B spectrum) and the full-atlas direct *curves* (Phase C), but no
    empirical downsample-and-measure *gating* check like EmeraldBay's 217/217. `pass4c_gating.py` supplies
    it from the 20 GB per-cell memmap: classify every (sample×line) condition at θ★=0.1 (per-line-mean
    baseline), and for each predicted OVER, downsample to n★ and check the realized RMS angle ≤ 1.05·θ★.
    Result: **5,503/5,503 (100 %)** of predicted-OVER conditions met the tolerance, across all 44 lines in
    which any over-sampled condition occurs — the EmeraldBay result reproduced at ~25× scale. Added one
    sentence to Phase C's population-scale paragraph, explicitly flagging that the per-line-mean baseline
    means these 5,503 are **not** the DMSO-referenced 2.5% OVER of Phase B (avoids a spurious
    contradiction). Fixture `tahoe_direct_curves_full.json` gains a `gating` block; SUPPLEMENT Constants of
    Record updated. Audit: only new number is 5,503; no existing number/table/citation changed. Rebuilt
    PDF (18 pp) / DOCX. With this, the empirical validation is atlas-wide on both sides — classification,
    direct curves, and gating verification all span the full atlas, not any cell-line subset.

28. **Reconciled `calibrate.py` with the authoritative recompute.** Running the demo surfaced a real
    inconsistency: its stored example vectors in `fixtures/tahoe_calibration.json` gave homoharringtonine
    5 µM × NCI-H460 **m=14.88 / n★=106**, but the authoritative `tahoe_quota_per_condition.csv` (and
    CASE_STUDIES, and the manuscript) say **m=15.53 / n★=98** — the demo vectors were partly stale
    (looked dose-pooled). Regenerated all three example vectors (Homoharringtonine/Panobinostat/Trametinib)
    directly from the dose-resolved `out_dose` pseudobulk with the exact `pass3_quota` computation
    (plate-matched DMSO_TF baseline); magnitudes now match the CSV to 4 dp (15.5333 / 9.9671 / 4.3947).
    Also fixed the honest-read print rounding (`.0f`→`.1f`) so it reports **2.5% OVER / 97.5% under-or-
    ghost** instead of 2% / 98%. `calibrate.py` now agrees end-to-end with the CSV, CASE_STUDIES, the
    manuscript, and the refreshed README; the 13-test suite still passes. (The earlier README edit that
    replaced 14.9→15.5 was in the right direction; this fixes the fixture it should have matched.)

29. **Mathematical-rigor audit of the problem statement (Methods + Results Phase A).** Fixed two
    conceptual ambiguities a computational-biology reviewer would flag, both verified against the source
    pipeline first. (1) **Matrix-vs-vector notation.** (A1)–(A3) now state the raw inputs are matrices
    \(X_t\in\mathbb{R}^{n_t\times d}\), \(X_c\in\mathbb{R}^{n_c\times d}\) (rows = individual \(d\)-vector
    cells); since \(n_t\neq n_c\) they cannot be subtracted, so we collapse the cell dimension to
    centroids and subtract the \(d\)-vectors, with the explicit estimator equation
    \(\hat v=\hat\mu_t-\hat\mu_c\), \(\hat\mu_t=\tfrac1{n_t}\sum x_i^t\), \(\hat\mu_c=\tfrac1{n_c}\sum x_j^c\).
    (2) **PCA→pseudobulk ordering.** Added a "Pipeline topology" paragraph: global embedding first (PCA on
    the pooled single-cell matrix), then cell-level covariance \(\Sigma\in\mathbb{R}^{50\times50}\) from
    individual embedded cells, then *late* pseudobulk as a downstream coordinate average within the
    embedding — explicitly *not* an upstream gene-count aggregation before PCA (which would collapse
    \(\Sigma\to0\) and destroy \(\operatorname{tr}(P\Sigma P)\)). Verified against `pass1_basis.py`
    (PCA on `vstack` of pooled cells), the per-cell `np.cov(C.T)`/within-condition residual σ², and the
    linear projected-mean centroid (`(sums/counts−pca_mean)@compsᵀ`). Notation propagated into Phase A.
    Audit: no number/table/citation changed. Clean compile, PDF (18 pp) / DOCX.

30. **Reviewer-facing revision: logic flags, style, biologist-friendliness.** A three-part pass, all
    prose/structure (no result number changed; audit confirms 50 [@cite] keys and the 5 data tables
    unchanged). **Logic:** (i) added a "Marginal-versus-within-condition variance" caveat with a
    `TODO(author)` — the theory calls for the within-condition residual variance but the headline Tahoe
    calibration uses the marginal σ²=2.406 (a *conservative* over-estimate), flagged for recompute/bias-
    quantification before submission (23,577 unchanged for now); (ii) softened the abstract's n★_δ claim
    from "controls Pr(θ>θ★)≤δ" to "approximate (1−δ) guarantee, exact as ρ→∞"; (iii) abstract now says the
    validation confirms measurement *geometry* (3/4 groups are DMSO_T0 vehicle populations, not drugs);
    (iv) new Discussion "Validity regime" paragraph elevating the ρ²≥3 first-order-breakdown caveat;
    (v) abstract clarifies Lean checks only the deterministic core. **Style:** abstract trimmed to ~265
    words; removed self-praising/literary language ("breakthrough", "principled answer", "four movements",
    "tell its story", "cashes out", "the crux"); defined "condition" and MOSAIC at first use; verified
    equation rendering (μ̂ hats/subscripts clean). **Biologist-friendly:** added a symbol glossary table
    and a "Quickstart (worked number)" box to Methods; marked the Laurent–Massart/Hanson–Wright derivation
    "For the mathematically inclined"; glossed "tangent space" as "the directions u can rotate into";
    inserted three figure placeholders with suggested captions (F1 geometry, F2 n★ histogram, F3 EmeraldBay
    θ²(n) curve). Rebuilt PDF (19 pp) / DOCX.

31. **Acted on the variance TODO: computed Tahoe within-condition σ² (documented as a refinement).**
    `pass5_within_sigma.py` computed the within-condition per-cell variance over all **95,624,334** Tahoe
    cells from the pass4 memmap: **σ² = 0.9567** --- almost exactly EmeraldBay's independent within-condition
    0.963, evidence the quantity transfers across platforms (the marginal 2.406 is inflated by
    between-condition/between-line structure). Rescaling the (unchanged) magnitudes gives n★=9,376/m²,
    median n★ 5,794, boundary m 2.69, and spectrum 10.6% OVER / 89.0% UNDER / 0.4% Ghost (vs 2.5/89.3/8.2).
    Per the author's decision, the conservative marginal figures stay the headline; the within-condition
    recompute is documented as a refinement (one Phase B paragraph + Methods pointer + a SUPPLEMENT row),
    and the TODO(author) is removed. Fixture `tahoe_within_sigma.json`. Audit: no existing headline number,
    citation, or data table changed; only the refinement numbers added. Rebuilt PDF (19 pp) / DOCX.

32. **Handoff.** Wrote `docs/HANDOFF.md` — current state (HEAD, Constants of Record, validation status,
    deliverable docs, build/reproduce commands, and the open items for the human author: three figure
    placeholders to realize, the marginal-vs-within-condition σ² decision, figure legends, author-list
    expansion, banner removal at finalization). Working tree clean; all pushed to origin/master.

33. **Switched fully to within-condition variance and regenerated all tables.** Per the author's decision,
    the headline calibration is now the theory-preferred within-condition residual σ² = **0.9567**
    (n★ = **9,376/m²**), and the marginal σ² = 2.406 (n★ = 23,577/m²) is retained everywhere only as a
    labelled *conservative bound*. `pass6_within_recalibrate.py` rescaled every Tahoe n★ by 0.3976 (m and
    N₀ are σ²-independent) and regenerated `tahoe_quota_per_condition.csv`, `tahoe_per_drug.csv`,
    `tahoe_per_dose.csv`, `tahoe_per_cell_line.csv`, `tahoe_constants.json`; `reliability_and_cost.py`
    (SIGMA2=0.9567) regenerated `tahoe_applications.json`. New spectrum **10.6% OVER / 89.0% UNDER /
    0.4% Ghost** (median n★ 5,794, boundary m 2.69); per-drug now shows majorities over-sampled
    (Panobinostat 120/150, Homoharringtonine 100/150) and only 28/379 drugs over-sampled in no condition;
    budget for the 600-condition screen 10.9M→3.5M cells (~$3.27M→$1.06M, 3.1×); reliability 10.6/42.4/67.1
    resolved, edges 10.0/39.2/65.0; dual-use exemplar homoharringtonine 5µM×NCI-H460 n★=39 (~156× multiplex).
    Propagated through the manuscript, SUPPLEMENT (Constants of Record + honesty ledger), README,
    CASE_STUDIES, HANDOFF, TECHNICAL_BLUEPRINT, THEORY_PRIMER, FALSIFICATION, SCALE_AUDIT, and the
    calibration fixture (annotated). **EmeraldBay results unchanged** (held-out curves, gating 217/217,
    slope-test ratios) — they use EmeraldBay's own Σ; the variance-comparison caveat is reframed from
    "not comparable" to a like-for-like within-condition match (0.9567 vs 0.963, ~1%). Audit: no citation
    keys or validation results changed. Rebuilt PDF / DOCX.

34. **Cross-modality scale-up to genome-wide gene perturbation (X-Atlas/Orion).** Extended the study from
    chemical to *genetic* perturbation on the largest public CRISPRi Perturb-seq atlas (slaf-project/X-Atlas-Orion,
    SLAF/Lance): ~8M cells, 18,903 gene knockdowns × 2 lines (HCT116, HEK293T). Built `scripts/orion_recompute/`
    (pass0 head-block download + basis; **pass2 streams the full expression table** — all 17.4B HCT116 + 29.1B
    HEK293T nonzeros re-projected into PCA(50) on the fly, resumable, per-gene pseudobulk sufficient stats;
    pass3 quota). Projection validated to 1e-4 vs scanpy; `total_counts` == raw row sum exactly. Condition =
    (gene × line), control = pooled **Non-Targeting** centroid. **Results:** within-condition σ² = 0.913 (HCT116)
    / 1.033 (HEK293T) — transfers across modality (≈ Tahoe 0.957 / EmeraldBay 0.963); n* = 8,950 / 10,127 per m²;
    median m 0.14 / 0.35 (≈10× smaller than drugs); spectrum **0.0% OVER** / 28.3% / 71.7% (HCT116) and
    0.0% / 49.4% / 50.6% (HEK293T); median deficit ~2,260× / ~248×. Marginal σ² ≈ within (weak-perturbation
    fingerprint). Self-validations: strongest-m knockdowns are core essential genes (RPL*, POLR2*, SF3B5, NUP93,
    EIF2S3…); HEK293T aniso/iso ratio 0.66 (first dataset where the anisotropic form departs from isotropy;
    THEORY §8). Added `docs/ORION_GENE_PERTURBATION.md`, SUPPLEMENT §2b + honesty items 9–10, README section,
    manuscript **Results Phase E + Table 6** (+ Abstract/Key-Points/Discussion/Limitations/Data-availability),
    citation `huang2025xatlasorion` (Crossref-verified; 49/49 DOIs pass). Fixed `normalize_manuscript.py` to keep
    `.csv/.sh/.npz/.bib` paths monospace (a `.csv` brace path had become math → double-subscript LaTeX error).
    Rebuilt PDF (21 pp) / DOCX. **Honest scope:** this applies the proven framework + transferred calibration;
    a dedicated Orion downsample-falsification (genetic Phase C) is the next step, and weak-majority magnitudes
    (median ~155–204 cells) are order-of-magnitude near the sampling floor.

35. **Validated essential-gene extension (TRADE; manuscript Phase F).** Extended to a third source/regime: the
    TRADE essential-gene CRISPRi screens (Nadig/Replogle/Weissman, Nat Genet 2025; GEO **GSE264667**),
    Jurkat-Essential + HepG2-Essential (2,393 DepMap common-essential genes/line, dual-sgRNA, pooled NTC). Unlike
    the 1e8-cell Orion streams, these are ~2.5e5 cells/line → load in RAM, per-cell coords retained → the
    parameter-free **downsample-and-measure falsification runs directly** (the genetic Phase C Orion couldn't do).
    Built `scripts/trade_recompute/` (pass0 GEO download; pass1 QC+global-embedding+full-50×50-Σ+quota; pass2
    falsification; driver) + `docs/TRADE_EXTENSION_PLAN.md`. Topology-faithful; QC per TRADE (low-UMI frac×median
    0.14/0.18, mito_UMI>1750/3000, single-or-dual-same-gene guides). **Results:** σ²_within 1.50 (Jurkat) / 1.87
    (HepG2) — higher than the ~1.0 band (strong perturbations flatten PCA spectrum); large-NTC-pool constant
    7,338 / 9,147 per m² (≈9,000 band). median m 1.84 / 2.25 (essential genes STRONG); spectrum 0.2/85.5/14.3 and
    0.6/81.1/18.2 % OVER/UNDER/Ghost; median deficit 48× / 70×; **depth-limited, not magnitude-limited** (median
    48–85 cells/gene). **Falsification (headline):** strong-knockdown (m>4) a-priori vs realized slope ratio
    **1.02 / 0.96 at R²=0.996**, parameter-free, with predicted low-SNR breakdown (Spearman 0.92/0.79) — first
    end-to-end *validated* gene-perturbation result. Biology recovered per line (Jurkat: RNA exosome EXOSC2-9,
    splicing, Mediator; HepG2: EIF*/PSM*/RPL*). Added manuscript **Results Phase F + Table 7** (+ Abstract, Key
    Point, Discussion, Limitations 11, Data availability GSE264667), SUPPLEMENT §2c + honesty 11–12, README,
    citation `nadig2025trade` (Crossref-verified; 50/50 DOIs). Fixtures `trade_{jurkat,hepg2}_{quota,summary,
    falsification}`. Rebuilt PDF (23 pp) / DOCX. Honest: σ² transfer looser here; QC low-UMI filter near-vacuous
    in the deposited h5ad (mito is binding).

36. **Manuscript figures generated from fixtures + embedded** (`427b306`). `scripts/make_manuscript_figures.py`
    builds five vector figures (PDF + PNG) from committed fixtures, Okabe-Ito colourblind-safe palette:
    F1 geometry schematic, F2 Tahoe n★ spectrum, F3 EmeraldBay held-out validation, F4 cross-modality summary
    (magnitude distributions + regime bars), F5 TRADE falsification. Replaced the three manuscript placeholders
    with real vector figures and added F4/F5. Caught a units bug (F3 fixture stores θ_RMS, squared to θ²) and a
    regime-% source (F2 now reads the per-condition `regime` column = 10.6/89.0/0.4). `figures/*.pdf` tracked via
    a `.gitignore` exception. PDF -> 25 pp.

37. **Read the four source papers; rebuilt Discussion into the standard five-part structure** (`6181288`,
    `4b08906`). Author added tahoe.pdf / Orion.pdf / nadig.pdf (+ EmeraldBay = same Vevo platform, 5-day).
    Discussion reorganized into Summary of key findings / Interpretation and explanation / Comparison with
    previous studies / Implications / Limitations and recommendations. Added the mechanistic "why" grounded in
    TRADE's effect-size distribution (typical genome-wide KD affects ~45 genes vs essential 500+; response-vs-
    experiment separation; GATA1-vs-EIF4A3 = our m bias-correction). NOTE: the reorganization script had a bug
    that silently dropped the two continuation lines of the "design rule" paragraph -- fixed in phase 41.

38. **Results restructured per the Results-section guideline** (`967a19b`). Author supplied an SJSU Writing-
    Center handout: Results should state facts (not interpret), open with an intro tied to the research question,
    close with a summary. Added the intro + a "Summary of findings" closing paragraph, and moved interpretation
    from the six analyses into the Discussion (kept every number/table/figure/citation).

39. **Completed the Methods Datasets section + fixed EmeraldBay platform framing** (`e3dec22`, `468de60`,
    `849d332`). Cross-checked all dataset numbers against the source papers (Tahoe "379 distinct drugs" verbatim;
    Orion 18,903 genes; TRADE 2,393 essential). Added **X-Atlas/Orion** and **TRADE** dataset paragraphs to
    Methods. **EmeraldBay is the same Vevo Mosaic platform as Tahoe (5-day, not a different platform)** -- reframed
    every "cross-platform" claim for Tahoe<->EmeraldBay as "across timepoints on the same platform", reserving
    cross-platform/cross-modality for the genetic atlases. Fixed stale integration-test numbers (four groups
    0.6/1.2/1.2/2.9%), an inaccurate "eight cell lines", an unexplained "(15,200 groups)", the overloaded "five
    datasets", TRADE depth 45->48, and added in-text Figure 1-5 callouts.

40. **Style pass: no semicolons, no "Phase A/B/C" scheme** (`33ff02b`). Per author direction: (i) removed all
    117 prose semicolons (comma before a coordinating conjunction, else period), keeping the 36 required pandoc
    `[@a; @b]` citation separators and fixing 4 semicolons baked into figure text; (ii) renamed the six Results
    section headers to descriptive titles and rewrote ~30 in-text "Phase X"/"Phases E--F" cross-references as
    descriptive pointers ("the sample-sufficiency spectrum", "the out-of-distribution validation", "the genetic
    atlases below", etc.). **The manuscript no longer uses "Phase A-F" labels** -- DEVLOG entries above keep them
    only as historical shorthand for the six analyses. Caught a bug where the semicolon pass corrupted LaTeX `\;`
    thin-space commands (`\;=\;` -> `\.` = U+0307 dot accent) and broke the tectonic build; restored the estimator
    equation.

41. **EmeraldBay added to Figure 4; Discussion de-bolded + lost content restored** (`d9c519b`, `76efbe9`). Added
    EmeraldBay to Figure 4 (a second chemical dataset was missing) using its 1,064 per-group m/N from
    `emeraldbay_falsification_full.json` (median m 0.71 per-line-mean baseline; regime 4.3/80.0/15.7 at θ=0.1).
    Per author direction, removed the 19 bold paragraph-lead "titles" from the Discussion so it reads as flowing
    prose (kept the five section headers and the numbered limitations list). While verifying, found the phase-37
    restructure bug had dropped the "design rule" paragraph's tail-bound confidence-quota equation (n★_δ) and its
    three use points -- restored from commit 6181288, semicolon-free. Final state: 0 "Phase" labels, 0 visible
    semicolons, 0 unresolved citations, Tables 1-7 + Figures 1-5 all referenced, PDF 28 pp.

    **Known open item (author decision pending):** Figure 4(b) shows EmeraldBay's θ=0.1 regime (4.3% over), a
    number not narrated in the text, while the text gives its gating at θ=0.20 (5.5% over) -- different tolerance,
    not contradictory, but two "% over" figures now coexist. Also a style choice left open: Methods/Results still
    use ~11 bold run-in paragraph headings (e.g. "Pipeline topology.", "Regime-gating accuracy.") whereas the
    Discussion no longer does.

42. **First-order-limit / second-order-deficit passages elevated in three sections.** Following a geometric
    stress-test of the Delta expansion, added formal treatment of the first-order approximation boundary so the
    ~13% low-SNR slope deficit reads as a forecast, not a flaw. (i) **Methods:** new validity-regime paragraph
    closing the Delta subsection -- exact geometry `tanθ = ‖Pe‖/(m+uᵀe)`, the along-signal `uᵀe` as a fluctuating
    lever arm, the `Dg(v)=P/m` annihilation of parallel noise, regime `ρ²≫1`, and the emergent convexity coupling
    `+3‖Pe‖²(uᵀe)²/m⁴` plus arc term `−⅔‖Pe‖⁴/m⁴`. (ii) **Results:** kept the Encorafenib sentence verbatim and
    appended the deterministic-curvature explanation -- `(1/n−1/N)²` concavity, "slope shift not loss of fit"
    (R²>0.996), sign+`1/ρ²`-growth as theory-specified. (iii) **Discussion:** new robustness paragraph -- a
    percentage-level quota error cannot flip the OVER/UNDER decision (4.5-fold median gap, boundary m=2.69), the
    low-SNR error lands on already-under-sampled conditions, gated OVER conditions sit at high SNR (<1% correction),
    confirmed by 5,503/5,503 and 217/217 gating. De-duplicated the `tanθ` equation from the confidence-quota caveat
    (now cross-references the Methods derivation). **Audit vs HEAD:** citation keys 54→54 (0 added/removed), table
    rows 69→69, 0 numeric tokens dropped, 0 new unsourced numbers (every figure reused), 0 unresolved citations in
    the rebuilt PDF, 0 prose semicolons. PDF 28→29 pp; DOCX rebuilt. Deliberately did NOT put the raw Table-3
    slopes (9.67/11.17) in the body -- not in SUPPLEMENT; the Discussion claims decision-robustness, not low-SNR
    quota-accuracy, keeping the ghost/unresolvable framing intact.

43. **Normalizer fix: directory paths kept as code, not math.** A rendering glitch surfaced during the
    phase-42 PDF verification -- inline paths `scripts/tahoe_recompute/` and `scripts/emeraldbay_recompute/`
    typeset with the `_r` as a subscript ("tahoe_r ecompute"). Root cause was in `scripts/normalize_manuscript.py`,
    not the manuscript: its `CODE_RE` keeps a backtick span as `code` only if it matches a known pattern
    (file extensions, Lean names, ...), else converts it to `$...$` math. Directory paths (trailing slash, no
    extension) matched nothing, so the underscore entered math mode. Added `(scripts|fixtures|docs|agents|src|
    tests|lean)/` to `CODE_RE`. Scanned the normalized output for all mathified code identifiers -- exactly the
    two directory paths, both now render as monospace with literal underscores. `manuscript.md` diff is 2 lines
    (`$...$`->`` `...` ``); numeric-token and citation hashes byte-identical to HEAD (0 drift). PDF/DOCX rebuilt,
    29 pp. Source `MANUSCRIPT_DEEPSEEK.md` was already correct (backticks); only the normalizer was patched.

44. **Methods restructured for biologist readers (two-layer): heavy derivations relocated to Supplementary
    Notes S1/S2.** Per author direction (audience = biologists, not mathematicians; keep rigor). The two densest
    Methods subsections were lightened in the main text and their derivations moved -- not deleted -- to a new
    `# Supplementary information` section placed after Data availability, before References. (i) **Delta-method
    subsection:** rewrote the intro around a plain-language lever-arm intuition (sideways jitter / effect length),
    kept the key formulas (`E[θ²]`, isotropic reduction) for quantitative reviewers, and replaced the ~250-word
    second-order paragraph with a one-sentence pointer to S1. (ii) **Distribution/tail-control subsection:**
    replaced the whole chi-square + Laurent-Massart + confidence-quota block (~35 lines) with a 6-line "Beyond the
    average error: a confidence quota" paragraph pointing to S2. (iii) **New S1** (Jacobian, tangent space, exact
    `tanθ` geometry, second-order correction -- source of the 13% deficit) and **S2** (generalized chi-square,
    `d_eff`, Laurent-Massart/Hanson-Wright, `n★_δ`, caveat) carry the relocated text verbatim. **Audit vs HEAD:**
    citation keys 54->54 (identical), table rows 69->69, every scientific number preserved. Two intentional
    numeric-token deltas: `@mardia2000directional` de-duplicated (doubled cite 2->1, key still present) and the
    relocated confidence quota's `\tag{5}`/`eq. (5)` dropped (now untagged in S2). 0 unresolved citations, 0 prose
    semicolons. Rebuilt PDF (29 pp) and DOCX; visually verified the lightened Methods (pp 4-5) and S1/S2 (pp 24-26)
    render correctly with proper math typesetting and correct S1-before-S2 ordering. Rigor is unchanged -- proofs
    relocated, not removed, and still backed by `docs/THEORY.md` + the Lean core.

## Honesty ledger (carried in the manuscript banner + `SUPPLEMENT.md`)

- The two high-magnitude validation groups are DMSO time-zero populations, not drug effects; only one
  group (Encorafenib × HEC-1-A) is a drug — the held-out test validates geometry, not biology.
- The headline uses the within-condition residual σ² on both atlases (Tahoe 0.9567 ≈ EmeraldBay 0.963,
  a like-for-like match to ~1%); the marginal Tahoe σ² = 2.406 is kept only as a labelled conservative
  bound and is not directly comparable to a within-condition estimate. Both establish σ² as a plug-in.
- Gating tolerance in the EmeraldBay experiment is 0.20 rad (not the standard 0.1).
- Resource figures are Tahoe-100M-derived; compute gains are polynomial, never exponential.
- Author lists beyond the third name are abbreviated and should be expanded before submission.
