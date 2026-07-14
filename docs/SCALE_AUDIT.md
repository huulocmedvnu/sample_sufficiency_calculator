> ⚠️ **SUPERSEDED.** Numbers here predate the audit (findings #1-#9, gaps 1-3). The authoritative source is `docs/MANUSCRIPT_DEEPSEEK.md` + `src/engine.py`. See `docs/AUDIT_DENOMINATORS.md` and `docs/ARCHITECTURE.md`.

# Scale-Logic Self-Audit & Alignment Ledger

## 0. AUDITOR CORRECTION (read first)

The drafting premise — *"our calculated n★ values sit in the range of hundreds to a few thousand
cells"* — is **true only for the strong-perturbation tail**, not the typical condition. At the standard
configuration (σ² = 0.9567 within-condition [headline], d = 50, **θ★ = 0.1 rad**; the marginal σ² = 2.406 is retained only as a conservative bound), over the 56,827-condition panel (379 drugs ×
3 doses × 50 lines):

| n★ percentile | 10th | 25th | **50th (median)** | 75th | 90th | 95th |
|---|---:|---:|---:|---:|---:|---:|
| cells/arm | 1,136 | 2,581 | **5,794** | 10,998 | 18,146 | 24,214 |

Only **10.6%** of conditions are over-sampled (n★ < N₀ = 1,296); ~13% are < 1.4 k, ~58% need 1.4 k–10 k,
**28% need 10 k–50 k**, and just 0.4% are ghosts (> 50 k). The "hundreds-of-cells" regime (e.g.
homoharringtonine at 5 µM in NCI-H460, n★ = 39) is the high-magnitude / high-dose / responsive-line
minority (m > 2.69). Because n★ ∝ 1/θ★², the "few thousand" regime for *typical* conditions requires a
looser tolerance (e.g. θ★ = 0.2–0.3 rad). **Manuscript framing must therefore be: at a tight ~5.7° MoA
tolerance the majority (89.4%) of (drug × dose × line) conditions are *under*-sampled — the atlas
is large by aggregating many conditions, not because any single one is cheap.**

---

## 1. Self-audit — n★ is strictly per-arm; the global atlas size does not enter

A condition's centroid is the mean over **that condition's** cells (aggregated over its replicate wells —
§4): μ̂_t = (1/n_t) Σ_{i ∈ condition} x_i. By the CLT its sampling covariance is

> Cov(μ̂_t) = Σ_t / n_t,   n_t = cells in the treated arm of the condition (not the atlas).

Hence the perturbation-vector noise is **S = Σ_t/n_t + Σ_c/n_c**, and the quota solves
tr(P S P)/m² = θ★² for n (per arm). The variance of a sample mean depends only on the per-cell
covariance Σ and the **number of cells averaged**; the count of cells residing elsewhere in the atlas is
irrelevant. The global size N_atlas = Σ_wells n_well is an **output aggregate**, never an input — it
appears in **no** term of S.

**Where a population size *does* enter — the finite-population correction.** If a condition's centroid is
estimated by sub-sampling n cells *without replacement* from the **N₀ cells actually acquired for that
condition** (aggregated over all its replicate wells; see §4), the sample-mean covariance is
Σ/n · (N₀ − n)/(N₀ − 1) ≈ Σ (1/n − 1/N₀). Thus

> E[θ²] = tr(PΣP)/m² · (1/n − 1/N₀).

Here **N₀ is the per-condition acquired pool**, not the atlas. As N₀ → ∞ (super-population / fresh
acquisition) the FPC → 1 and E[θ²] → tr(PΣP)/(n m²). **Conclusion: n★ ⟂ atlas size; the only "N" in the
mathematics is the per-condition acquired pool, and even that vanishes for large N₀.**

---

## 2. Scale-Alignment Ledger — per-condition quota → project budget

**Aggregation identity.** For a screen of D drugs × Z doses × L lines = K treated **conditions** (each
condition may be realized across several replicate wells — §4):

> Cells_to_load = f_recovery · [ Σ_{c=1}^{K} n★_c  +  Cells_control ]
> Reads = Cells_to_load · depth_reads_per_cell ;   Cost = Reads · price_per_read

— a **sum over independent per-condition quotas**, each n★_c fixed by that condition's magnitude m_c
(n★ = 9,376/m² at the standard config) and counting the *total* cells of the condition however many
wells supply them (§4). f_recovery ≈ 1.4 (droplet capture + doublet/QC loss); the control/vehicle pool is
sized once per batch and amortized (it is not multiplied by K).

**Worked example — 100 drugs × 3 doses × 2 lines = 600 treated conditions** (conditions sampled from the
real Tahoe magnitude distribution; per-drug median n★ 5,990, IQR [3,888–8,564]):

| Budgeting policy | treated cells | notes |
|---|---:|---|
| **A. Flat at 90th-pctile n★ = 18,146** | **10.9 M** | uniform over-loading; reviewer-naïve; wasteful |
| **B. Magnitude-adaptive, cap 10,000/condition** | **3.5 M** | 54/379 drugs hit the cap → resolved at coarser-than-0.1-rad θ |
| **B′. Adaptive, cap 20,000/condition** | **4.5 M** | fewer capped |
| **C. Adaptive + relaxed θ★ = 0.2 rad, cap 10 k** | **1.2 M** | n★ ÷ 4; most conditions now within cap |

Add f_recovery (×1.4) and a control pool (~10 k cells/batch): Policy B ≈ **4.9 M cells to load**, which at
~20 k reads/cell is ~1.0×10¹¹ reads. **This is how "small" per-condition numbers compose into a
multi-million-cell project: the magnitude is driven not by any single quota but by (i) the number of
conditions K and (ii) the long tail of weak conditions whose individual n★ is large.**

**The calculator's translational role** is therefore *not* to make screens cheap but to make the budget
**rational and explicit**: (1) magnitude-adaptive allocation replaces flat loading (Policy A → B saves
~3.1×); (2) per-drug **triage** — drugs whose n★ exceeds the per-well cap are flagged *unresolvable at the
chosen θ★* (ghosts) rather than silently under-powered; (3) an explicit **tolerance/depth trade** (Policy
B → C: relaxing θ★ from 5.7° to 11.5° quarters every quota).

---

## 3. Methods-section definition — "Cell Quota per Perturbation"

> **Cell quota per perturbation (n★).** The number of single cells that must be acquired *and pass
> quality control* for one **perturbation condition** — a single (drug × dose × cell-line) treated arm —
> such that the root-mean-square angular error of that condition's estimated perturbation-direction
> vector v̂ = μ̂_treated − μ̂_control, measured against its reference (vehicle) centroid in the shared
> d-dimensional embedding, does not exceed a pre-specified tolerance θ★. Formally
> n★ = 2·tr(PΣP)/(m²θ★²) (isotropic special case 2(d−1)σ²/(m²θ★²)), where σ²/Σ is the per-cell
> *within-condition* covariance, m = ‖v‖ the perturbation magnitude, and P = I − uuᵀ the tangent-space
> projector. n★ is a **per-arm, per-condition** quantity counting the *total* cells aggregated into the
> condition's centroid — **not** the cells of any single physical well; it is independent of the number of
> conditions screened and of the total atlas size, which enter budgeting only through Σ_c n★_c. The
> control/vehicle arm is sized separately and may be an amortized pool; where a condition is estimated by
> sub-sampling from N₀ already-acquired cells the achievable error carries the finite-population factor
> (1/n − 1/N₀) with N₀ the per-condition acquired pool. **Replicate pooling is valid only under the
> batch-exchangeability conditions of §4: pooling reduces only the within-condition variance Σ_cell/n;
> any plate/batch variance Σ_batch is reduced by the number of replicate wells R, not by cells.**

*(Constants of record and provenance: `docs/SUPPLEMENT.md`. Numbers above generated from
the reference perturbation array via the calculator.)*

---

## 4. The experimental unit: replicate pooling and batch variance

### 4.1 What n★ counts (audit answer, Q1)
The code returns *"cells required per arm"* and the proof's assumption **(A1)** is that the arm's cells
x₁,…,x_n are **i.i.d.** with covariance Σ; the estimator is their mean μ̂ = (1/n)Σ xᵢ. The quantity n is
therefore the **number of cells aggregated into the condition's centroid**, with no reference to physical
wells. **n★ is the quota per unique perturbation condition (treated arm) — the total cells whose mean
forms the centroid — not per individual physical well.** The earlier "per well" phrasing is exact only in
the degenerate case where one condition = one well, or where replicate wells are perfectly exchangeable.
A lab may legitimately reach n★ by aggregating cells of the *exact same (drug × dose × line)* condition
across replicate wells — **provided assumption (A1) survives that aggregation**, which §4.2 makes precise.

### 4.2 Does pooling R under-sampled replicates preserve validity? (audit answer, Q2)
Pooling is valid **iff the replicate wells are exchangeable draws from one distribution.** A plate/batch
effect breaks this. Model each cell in replicate well r as

> xᵣᵢ = μ + bᵣ + εᵣᵢ ,  bᵣ ~ (0, Σ_batch) i.i.d. across wells,  εᵣᵢ ~ (0, Σ_cell) i.i.d. within,

where **bᵣ is a per-well/batch mean shift** shared by all cells of well r (the technical/plate confounder)
and Σ_cell is the within-condition biological+technical noise the calculator already uses. For R wells of
nᵣ cells (n = Σ nᵣ) the pooled-centroid covariance is

> Cov(μ̂) = (1/n²) Σᵣ nᵣ² · Σ_batch + (1/n) · Σ_cell   →(balanced nᵣ=n/R)→   **Σ_batch / R + Σ_cell / n .**

The two terms scale **differently**: the within-condition term Σ_cell/n is driven down by **total cells
n** (what n★ controls), but the batch term Σ_batch/R is driven down **only by the number of replicate
wells R**, *independently of how many cells each well contributes*. Hence the angular error of the
perturbation vector (un-cancelled batch, equal arms) is the **two-component**

> E[θ²] ≈ (2/n)·tr(P Σ_cell P)/m²  +  (2/R)·tr(P Σ_batch P)/m² ,

with an **irreducible floor** θ²_floor = (2/R)·tr(P Σ_batch P)/m² that *no amount of cell-pooling within a
fixed set of wells can cross.* **Therefore: simple cell-pooling does NOT in general make an
under-sampled condition sufficient. It suffices only when Σ_batch ≈ 0 or is removed.**

### 4.3 How batch is removed: matched-control (vehicle) referencing
The perturbation vector v̂ = μ̂_treated − μ̂_control **cancels the additive bᵣ** when treated and control
share the same wells/batches (co-plated, common offset): computing v̂ᵣ = (treated − co-plated vehicle)
**within each well/batch** and then aggregating ∑ᵣ wᵣ v̂ᵣ (wᵣ ∝ cells) eliminates the Σ_batch/R term and
restores Cov(v̂) = Σ_cell·(1/n_t + 1/n_c) — i.e. **the i.i.d. formula with n = total pooled cells holds
again.** This is precisely the role of the calculator's control arm.
*Limitation (honest):* differencing cancels only **additive (location)** batch. **Multiplicative or
rotational** batch — a per-batch rescaling or basis-rotation of Σ — is *not* removed by referencing and
requires explicit integration (Harmony / scVI / ComBat); residual non-additive batch is documented to
occur in real perturbation atlases and is out of scope of the first-order quota.

### 4.4 Replication-strategy guidelines for wet-labs
Let **R★ = 2·tr(P Σ_batch P)/(m² θ★²)** — the *replicate-count* analogue of n★ (same form, batch
covariance and R in place of Σ_cell and n). Estimate Σ_batch from a pilot as the between-replicate scatter
of per-well centroids (after referencing); Σ_cell from the within-well residual.

1. **Interpret n★ as a per-condition cell budget**, reachable by pooling replicate wells — *not* a mandate
   to put n★ cells in one physical well.
2. **Always reference within batch.** Compute the perturbation vector per well against a **co-plated
   vehicle**, then pool the per-well vectors weighted by cells. This cancels additive batch and is what
   makes n★ (total cells) the correct, sufficient target.
3. **Prefer splitting n★ across several co-plated/multiplexed replicate wells over concentrating it in
   one well.** Splitting (i) supplies the R replicates that alone control the batch floor, (ii) de-risks
   single-well/plate failure, and (iii) with within-well referencing still reaches n★ total cells while
   suppressing batch. The only cost is a per-well vehicle overhead.
4. **Do not chase n★ by deeper sequencing of a single well when batch effects are plausible.** That
   reduces Σ_cell/n but leaves Σ_batch/1 (R = 1) — the estimate floors at the batch level.
5. **If matched referencing is impossible** (e.g. a single global control from a different batch), the
   condition needs **both** n ≥ n★ cells **and** R ≥ R★ independent replicate batches; report whichever
   binds. R★ ≪ n★ whenever Σ_batch ≪ Σ_cell (the usual regime), so a handful of replicate batches
   typically suffices once referencing is in place.

**Bottom line.** n★ is a per-condition *cell* quota and pooling replicate wells to reach it is valid —
but only after within-batch vehicle referencing removes the additive plate/batch offset; otherwise the
governing unit becomes two-dimensional, (cells n ≥ n★) **and** (replicate batches R ≥ R★), because cell
count and replicate count suppress orthogonal variance components.

---

## 5. Empirical pooled-condition re-gating: does replication rescue these atlases? (audit answer)

§4 establishes that n★ is a **per-condition** quota reachable by pooling replicate wells. A natural
hypothesis is that under-sampled single wells are rescued once replicates are pooled
(n_pooled = R·n_well). **We tested this directly on the data; in these two atlases it is false**, because
the conditions are essentially **unreplicated**.

**Replication structure (measured).** In the full Tahoe atlas (56,827 drug×dose×line conditions),
each drug is screened at ~3 doses, but a given (drug, dose) occupies **a single plate for 86.3%** of
combinations, two plates for 13.4%, and ≥ 3 for 0.4%: same-dose replicate plates are largely absent.
Median cells per (drug × dose × line) condition is **1,296** (post-filter). EmeraldBay mirrors this —
one pooled sample per (drug, dose, line); no independent same-dose replicate wells exist in the data.

**Re-gate (per-condition magnitude and cells, θ★ = 0.1 rad).** Over the 56,827 Tahoe conditions:

| Gating model | OVER-sampled (n ≥ n★) |
|---|---:|
| Per (drug × dose × line) condition | **10.6%** |
| Under-sampled or Ghost | **89.4%** |

The taxonomy is stark: the median condition needs n★ ≈ 5.8 k cells but holds ≈ 1.3 k, an **~4.5×
deficit** that the sparse replicate plates cannot close. **Replication does not rescue under-sampled
conditions here — not because the pooling logic is wrong (§4), but because the data architecture supplies
≈ one plate per drug–dose.** Dose, by contrast, is a real lever (the over-sampled fraction rises from
8.8% at 0.05 µM to 14.2% at 5 µM), but it does not move the bulk of the panel out of under-sampling.

**Case-study profiles resolved across all 150 conditions per drug (3 doses × 50 lines):**

| Drug (profile) | median m | median n★ | conditions OVER-sampled |
|---|---:|---:|---:|
| Panobinostat (HDAC) | 4.24 | 523 | **120 / 150** |
| Homoharringtonine (protein synthesis) | 5.88 | 271 | 100 / 150 |
| Idarubicin (anthracycline; TOP2A) | 2.93 | 1,091 | 70 / 150 |
| Palbociclib (CDK4/6) | 1.41 | 4,684 | 4 / 150 |
| Crizotinib (ALK/MET; ghost tail) | 0.99 | 9,632 | 1 / 150 (3 ghost) |

This **confirms** the single-well profiles (high-magnitude agents over-sampled in most lines; moderate and
weak agents under-sampled in nearly all), and shows the regime is **line-dependent** (a drug can be
saturated in its high-response lines yet under-sampled in the rest). **Bottom line: the pooled-condition
unit is the mathematically correct one (§4), but in Tahoe and EmeraldBay it coincides with the
single-condition unit, so the project-level conclusion stands — ≈ 10.6% of conditions are safely saturated
at a 5.7° tolerance; replication as actually practised in these atlases rescues < 0.1%.** A screen
*designed* with R replicate batches and within-batch referencing would benefit (n★ is per-condition); a
single-sample-per-condition atlas does not.