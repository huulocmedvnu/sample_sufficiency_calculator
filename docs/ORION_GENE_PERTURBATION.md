> ⚠️ **SUPERSEDED.** Numbers here predate the audit (findings #1-#9, gaps 1-3). The authoritative source is `docs/MANUSCRIPT_DEEPSEEK.md` + `src/engine.py`. See `docs/AUDIT_DENOMINATORS.md` and `docs/ARCHITECTURE.md`.

# Cross-modality scale-up: genome-wide gene perturbation (X-Atlas/Orion)

**Question.** The sample-sufficiency quota `n★ = 2·tr(PΣP)/(m²θ★²)` was derived and calibrated on *chemical*
perturbation (Tahoe-100M, EmeraldBay). Does it transfer to *genetic* perturbation — a different modality, a
different platform, a different lab? The geometry is a proven theorem (Lean core + Monte-Carlo), so data cannot
re-prove it; what an independent genetic atlas *can* test is whether the calibrated plug-in (σ², and the
resulting `n★ ≈ 9,000/m²` law) transfers, and what the quota says when applied at genome scale.

**Answer.** It transfers. The within-condition per-cell PCA variance σ² lands in the same ~1.0 band on both
X-Atlas/Orion cell lines as on both chemical atlases; the strongest-magnitude knockdowns are the canonical
essential genes; and applied genome-wide the quota shows that **no gene knockdown is over-sampled for direction
estimation at a 5.7° tolerance** — a far more extreme version of the chemical finding, driven by the ~4–10×
smaller transcriptional magnitudes of gene knockdowns.

---

## 1. The dataset and how the method maps onto it

**X-Atlas/Orion** [@huang2025xatlasorion] is the largest public genome-wide **CRISPRi Perturb-seq** atlas: two
Fix-Cryopreserve-scRNAseq (FiCS) screens targeting all 18,903 human protein-coding genes, ~8M cells across two
lines (HCT116 colorectal, HEK293T embryonic kidney), median ≥140 cells per perturbation, ~16,000 UMIs/cell,
with a large pooled **Non-Targeting (NTC)** control (165k / 219k cells). Median on-target knockdown is 75.4%
(HCT116) and 51.5% (HEK293T). We read it from the SLAF/Lance re-release (`slaf-project/X-Atlas-Orion`).

The framework maps one-to-one, with the genetic analog of each chemical quantity:

| framework quantity | chemical (Tahoe/EmeraldBay) | genetic (X-Atlas/Orion) |
|---|---|---|
| condition | drug × dose × cell line | **target gene × cell line** |
| control arm | plate-matched DMSO centroid | **pooled Non-Targeting (NTC) centroid** |
| effect magnitude `m` | ‖μ_drug − μ_DMSO‖ | **‖μ_knockdown − μ_NTC‖** |
| direction `u = v/m` | drug mechanism axis | **knockdown transcriptional response axis** |
| quota `n★` | cells/arm for direction to θ★ | identical |

The quota theorem, tail bound, and Lean-checked deterministic core are unchanged — they are indifferent to what
caused the displacement. Grouping knockdowns by the *direction* of their response (pathway/complex membership) is
the central analysis Perturb-seq is built for, so the "how many cells to resolve the direction" question is, if
anything, more native here than for drugs.

## 2. Pipeline (`scripts/orion_recompute/`)

Same fixed recipe as the chemical recomputes: `normalize_total(1e4)` → `log1p` → `highly_variable_genes(2000)`
→ `PCA(50)`. `pass0` downloads a representative head-block per line and `pass1` fits the HVG/PCA basis; **`pass2`
streams the *entire* expression table** (HCT116 17.4B + HEK293T 29.1B nonzeros) in cell-order chunks, projecting
each cell into PCA(50) on the fly and accumulating per-`gene_target` sufficient statistics
(count, Σcoords, Σcoords²) — resumable, no full in-memory hold; `pass3` derives σ² and the quota spectrum. The
streaming projection was validated to `1e-4` against `scanpy`, and `total_counts` is the exact normalize
denominator (0 diff vs the raw row sum).

## 3. Result 1 — σ² transfers across modality (the calibration is portable)

The within-condition per-cell variance — the noise that actually enters a centroid — is nearly identical across
four datasets, two modalities, two labs, three platforms:

| dataset | modality | within-condition σ² | quota law @ θ=0.1 |
|---|---|---:|---:|
| Tahoe-100M | chemical | 0.9567 | n★ = 9,376/m² |
| EmeraldBay | chemical | 0.938 | — |
| **X-Atlas/Orion HCT116** | **genetic (CRISPRi)** | **0.9133** | **n★ = 8,950/m²** |
| **X-Atlas/Orion HEK293T** | **genetic (CRISPRi)** | **1.0334** | **n★ = 10,127/m²** |

σ² is a platform-specific plug-in that *must* be re-estimated (the whole point of running it on a new modality),
and it lands at ~1.0 every time under this fixed recipe. The calibrated quota constant is therefore ~9,000/m²
across both modalities.

**A fingerprint of weak perturbations.** On both Orion lines the *marginal* σ² (all cells pooled) nearly equals
the *within-condition* σ² (HCT116 0.9167 vs 0.9133; HEK293T 1.0396 vs 1.0334), whereas on Tahoe the marginal is
2.5× larger than within (2.406 vs 0.9567). Marginal = within + between-condition variance, so this near-equality
is a direct, quantitative statement that **most gene knockdowns barely move the transcriptome** — there is very
little between-condition spread to add.

## 4. Result 2 — the genome-wide sufficiency spectrum

Applying `n★` to every knockdown against its true acquired cell count:

| | Tahoe (chemical) | **Orion HCT116** | **Orion HEK293T** |
|---|---:|---:|---:|
| perturbations | 56,827 conditions | 17,585 genes | 17,856 genes |
| median magnitude `m` | 1.27 | **0.14** | **0.35** |
| **% OVER-sampled** | 10.6% | **0.0%** | **0.0%** |
| % UNDER | 89.0% | 28.3% | 49.4% |
| % Ghost (n★>50k) | 0.4% | **71.7%** | **50.6%** |
| median n★ | 5,794 | 409,983 | 52,908 |
| median acquired depth | 1,296 | 155 | 204 |
| **median deficit (n★/N₀)** | 4.5× | **~2,260×** | **~248×** |
| aniso/iso ratio | 0.93–0.96 | 0.925 | **0.660** |

**Not a single knockdown is over-sampled** for direction at θ★ = 0.1 rad on either line. The mechanism is
`n★ ∝ 1/m²`: gene knockdowns produce ~4–10× smaller transcriptional displacements than drugs (median m 0.14–0.35
vs 1.27), so the quota inflates 20–80× while acquired depth is ~6–9× lower — the deficit reaches ~2,260× (HCT116)
and ~248× (HEK293T). The result is robust to tolerance: even at a loose θ★ = 0.3 rad (17°) only 0.4% (HCT116) and
2.1% (HEK293T) of knockdowns become over-sampled. HEK293T is the milder of the two because its knockdown
responses are larger (median m 0.35) and its depth higher (204 cells), not because its knockdown efficiency is
higher — it is in fact lower (51% vs 75%), so magnitude here is set by cell biology, not guide potency.

## 5. Result 3 — the pipeline recovers known biology (a validation the theorem can't give)

A theorem cannot tell you the *magnitudes* are right; biology can. The strongest-magnitude knockdowns on both
lines are dominated by **core essential machinery** — ribosomal proteins, ribosome biogenesis, RNA polymerase II,
splicing/mRNA processing, nucleoporins, translation initiation:

- **HCT116 top:** EIF2S3, RPL31, RPL18, MED22, NUP93, MAK16, UTP15, RPL7, RPL30, NOL8, EXOSC9, HEATR1
  (20 of the top 50 are ribosomal / Pol II / splicing / nucleoporin / translation genes).
- **HEK293T top:** RPL38, NUP93, ITGB1, POLR2C, POLR2I, RPS20, SF3B5, USPL1, SNIP1, CPSF4, PGAM2
  (13 of the top 50 in the same families).

Full-depth pseudobulk puts the knockdowns that *should* produce the largest transcriptional shifts exactly at the
top of the magnitude spectrum — evidence the streaming projection, the NTC-referenced magnitudes, and the
directions are correct.

## 6. Result 4 — HEK293T is where the anisotropic form earns its keep

On Tahoe, EmeraldBay, and Orion HCT116 the anisotropic/isotropic quota ratio is ~0.93–0.95: perturbation
directions carry little variance along themselves, so the anisotropic correction is small. **HEK293T is
different — ratio 0.66.** Its noise is far more single-axis-dominated (PC1 alone carries ~20% of the variance vs
~13% in HCT116), and knockdown directions partly align with it, so `tr(PΣP)` drops well below the isotropic
`σ²(d−1)` and the quota is ~34% smaller than an isotropic calculation would give. This is the first dataset in the
study where the anisotropic machinery departs materially from isotropy — exactly the regime `THEORY.md` §8
identifies (signal aligned with high-variance axes → easier than isotropy predicts), observed here on real data.

## 7. Honest scope and caveats

1. **This is an application of a proven framework, not an independent re-proof.** The geometry is a theorem
   (Lean §10); Orion tests whether the *calibration* transfers (it does) and what the quota says at genome scale.
   A dedicated downsample-and-measure falsification on Orion (the parameter-free slope test of Phase C) is the
   natural next step and is not yet run here.
2. **Weak-majority magnitudes are order-of-magnitude, not precise.** With a median of ~155–204 acquired cells per
   knockdown, the per-condition sampling floor on `m` is non-negligible; magnitudes are bias-corrected
   (subtracting `tr(Σ)(1/n+1/n_NTC)`), but for the weak majority — whose true `m` sits near that floor — `n★`
   should be read as "unresolvable at this depth," consistent with the first-order / ρ²≫1 validity regime, not as
   an exact cell count. About 15% (HCT116) / 35% (HEK293T) of knockdowns clear the floor with a clearly
   detectable direction; even among those, 0% are over-sampled.
3. **Knockdown efficiency as a continuous dose.** Orion's sgRNA-UMI abundance is a per-cell proxy for knockdown
   strength, giving a *continuous* within-condition dose axis (unlike the chemical atlas's discrete doses). Using
   it to trace `m` vs knockdown level is a genuine extension the drug data could not support; not exploited here.
4. **Scope unchanged.** `n★` governs centroid-direction precision only, not cell-level structure; magnitudes are
   transcriptional displacements of surviving cells, not viability; compute savings are polynomial.

## Provenance

Scripts `scripts/orion_recompute/{pass0_download_head,pass1_basis_sigma,pass2_stream_pseudobulk,pass3_quota_full}.py`
and driver `run_full_both_lines.sh`. Full-atlas fixtures `fixtures/orion_{HCT116,HEK293T}_{quota.csv,summary.json}`
(all ~17,600 knockdowns per line). Work dir `/mnt/hdd2/loc-tran/orion_work/full` (outside the repo; regenerable).
Constants consolidated in `docs/SUPPLEMENT.md`.
