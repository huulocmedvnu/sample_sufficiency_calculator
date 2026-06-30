# Scale-Logic Self-Audit & Alignment Ledger

## 0. AUDITOR CORRECTION (read first)

The drafting premise — *"our calculated n★ values sit in the range of hundreds to a few thousand
cells"* — is **true only for the strong-perturbation tail**, not the typical drug. At the standard
configuration (σ² = 7.66, d = 50, **θ★ = 0.1 rad**), over the 292-drug Tahoe panel:

| n★ percentile | 10th | 25th | **50th (median)** | 75th | 90th | 95th |
|---|---:|---:|---:|---:|---:|---:|
| cells/arm | 7,029 | 12,087 | **23,934** | 42,571 | 65,392 | 79,196 |

Only **2%** of drugs are over-sampled (n★ < N₀ = 1,394); 17% need 1.4 k–10 k, **64% need 10 k–50 k**,
17% are ghosts (> 50 k). The "hundreds-to-thousands" figure (e.g. Homoharringtonine n★ = 365) is the
m > 7.34 minority. Because n★ ∝ 1/θ★², the "few thousand" regime for *typical* drugs requires a looser
tolerance (e.g. θ★ = 0.2–0.3 rad). **Manuscript framing must therefore be: at a tight ~5.7° MoA
tolerance most atlas wells are *under*-sampled — the atlas is large by aggregating many wells, not
because any single well is cheap.**

---

## 1. Self-audit — n★ is strictly per-arm; the global atlas size does not enter

A condition's centroid is the mean over **that well's** cells:
μ̂_t = (1/n_t) Σ_{i ∈ well} x_i. By the CLT its sampling covariance is

> Cov(μ̂_t) = Σ_t / n_t,   n_t = cells in the treated well (not the atlas).

Hence the perturbation-vector noise is **S = Σ_t/n_t + Σ_c/n_c**, and the quota solves
tr(P S P)/m² = θ★² for n (per arm). The variance of a sample mean depends only on the per-cell
covariance Σ and the **number of cells averaged**; the count of cells residing elsewhere in the atlas is
irrelevant. The global size N_atlas = Σ_wells n_well is an **output aggregate**, never an input — it
appears in **no** term of S.

**Where a population size *does* enter — the finite-population correction.** If a condition's centroid is
estimated by sub-sampling n cells *without replacement* from the **N₀ cells actually acquired for that
well**, the sample-mean covariance is Σ/n · (N₀ − n)/(N₀ − 1) ≈ Σ (1/n − 1/N₀). Thus

> E[θ²] = tr(PΣP)/m² · (1/n − 1/N₀).

Here **N₀ is the per-well acquired count**, not the atlas. As N₀ → ∞ (super-population / fresh
acquisition) the FPC → 1 and E[θ²] → tr(PΣP)/(n m²). **Conclusion: n★ ⟂ atlas size; the only "N" in the
mathematics is the per-well acquired pool, and even that vanishes for large N₀.**

---

## 2. Scale-Alignment Ledger — per-condition quota → project budget

**Aggregation identity.** For a screen of D drugs × R doses × L lines = K treated wells:

> Cells_to_load = f_recovery · [ Σ_{w=1}^{K} n★_w  +  Cells_control ]
> Reads = Cells_to_load · depth_reads_per_cell ;   Cost = Reads · price_per_read

— a **sum over independent per-well quotas**, each n★_w fixed by that well's magnitude m_w (n★ = 75,068/m²
at the standard config). f_recovery ≈ 1.4 (droplet capture + doublet/QC loss); the control/vehicle pool
is sized once per batch/line and amortized across wells (it is not multiplied by K).

**Worked example — 100 drugs × 3 doses × 2 lines = 600 treated wells** (drugs sampled from the real
Tahoe magnitude distribution; per-drug n★ median 24,858, IQR [14,527–46,190]):

| Budgeting policy | treated cells | notes |
|---|---:|---|
| **A. Flat at 90th-pctile n★ = 65,392** | **39.2 M** | uniform over-loading; reviewer-naïve; wasteful |
| **B. Magnitude-adaptive, cap 10,000/well** | **5.7 M** | 85/100 drugs hit the cap → resolved at coarser-than-0.1-rad θ |
| **B′. Adaptive, cap 20,000/well** | **10.1 M** | 59/100 capped |
| **C. Adaptive + relaxed θ★ = 0.2 rad, cap 10 k** | **3.8 M** | n★ ÷ 4; most drugs now within cap |

Add f_recovery (×1.4) and a control pool (~10 k cells/batch): Policy B ≈ **8.0 M cells to load**, which at
~20 k reads/cell is ~1.6×10¹¹ reads. **This is how "small" per-condition numbers compose into a
multi-million-cell project: the magnitude is driven not by any single quota but by (i) the number of
conditions K and (ii) the long tail of weak drugs whose individual n★ is large.**

**The calculator's translational role** is therefore *not* to make screens cheap but to make the budget
**rational and explicit**: (1) magnitude-adaptive allocation replaces flat loading (Policy A → B saves
~7×); (2) per-drug **triage** — drugs whose n★ exceeds the per-well cap are flagged *unresolvable at the
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
> projector. n★ is a **per-arm, per-condition** quantity: it is independent of the number of conditions
> screened and of the total atlas size, which enter project budgeting only through the sum
> Σ_w n★_w over wells. The control/vehicle arm is sized separately and may be a single pool amortized
> across conditions; where a condition is estimated by sub-sampling from N₀ already-acquired cells, the
> achievable error carries the finite-population factor (1/n − 1/N₀) with N₀ the per-well acquired count.

*(Constants of record and provenance: `docs/SUPPLEMENT.md`. Numbers above generated from
the reference perturbation array via the calculator.)*
