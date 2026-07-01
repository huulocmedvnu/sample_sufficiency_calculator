# Lean 4 / Mathlib formalization of the angular-error core

Machine-checked proofs of the **deterministic linear-algebra backbone** of the cell-quota theory in
[`../docs/THEORY.md`](../docs/THEORY.md). Everything here is checked by the Lean kernel against
Mathlib — no `sorry`, no axioms beyond Mathlib's.

## What is proved

`SampleSufficiency/AngularError.lean` — the tangent-space projector `P = I − u uᵀ` and the trace
identity the quota is built from:

| Lean name | Statement | THEORY.md |
|---|---|---|
| `proj_transpose` | `Pᵀ = P` | §1 |
| `proj_mul_self` | `P P = P` (idempotent, unit `u`) | §1 |
| `proj_mulVec_self` | `P u = 0` (kills the radial direction) | §2 |
| `trace_proj` | `tr P = d − ‖u‖²` (= `d−1`, rank of the projector) | §1 |
| **`trace_proj_conj`** | **`tr(P Σ P) = tr Σ − uᵀ Σ u`** (only ⟂ noise survives) | §2–3 |
| `trace_proj_iso` | `tr(P (σ² I) P) = σ² (d−1)` (recovers the v1 isotropic formula) | §3 |
| `trace_proj_diag` | `tr(P diag(ℓ) P) = Σ_k (1 − u_k²) ℓ_k` (PCA basis) | §7, §9 |

`SampleSufficiency/Jacobian.lean` — Lemma 1, the normalization Jacobian, on any real inner-product
space:

| Lean name | Statement | THEORY.md |
|---|---|---|
| `hasFDerivAt_norm'` | `D‖·‖(v) = ‖v‖⁻¹ ⟪v, ·⟫` (`v ≠ 0`) | §2 |
| **`hasFDerivAt_normalize`** | **`D(x ↦ x/‖x‖)(v) = ‖v‖⁻¹ (I − ‖v‖⁻² v vᵀ) = (1/‖v‖) P`** | §2, Lemma 1 |
| `hasFDerivAt_normalize_apply_self` | `Dg(v) v = 0` (radial direction ∈ ker) | §2 |

Together these are exactly the algebraic facts the isotropic/anisotropic quota
`n* = 2 tr(PΣP)/(m² θ*²)` rests on: the derivative of the direction map is `(1/m)P` (so only
perpendicular noise rotates `û`), and `tr(PΣP)` collapses to `σ²(d−1)` in the isotropic case and to
`Σ_k(1−u_k²)ℓ_k` in a PCA basis.

## Scope (honest)

Lean here certifies the **deterministic** algebra. The **probabilistic** statements —
`E[θ²] = tr(PSP)/m²` (Isserlis/Delta method), the §5 second-order term, and the Laurent–Massart tail
`Pr(θ>θ*) ≤ δ` — are **not** formalized (they need Gaussian-quadratic-form machinery Mathlib does not
package). Those are checked empirically by [`../tests/verify_theory.py`](../tests/verify_theory.py)
(symbolic Jacobian + Monte-Carlo, incl. the corrected §5 term in Layer 4).

## Build / reproduce

Pinned to **Lean 4.31.0 / Mathlib v4.31.0** (see `lean-toolchain`, `lake-manifest.json`).

```bash
cd lean
lake exe cache get      # download prebuilt Mathlib oleans (~GBs); avoids compiling Mathlib
lake build              # checks AngularError.lean + Jacobian.lean  (seconds, given the cache)
```

`.lake/` (build artifacts, multi-GB) is git-ignored; `lake-manifest.json` is committed so the exact
Mathlib revision is reproducible.

To re-confirm the proofs rest only on the three standard Mathlib axioms:

```bash
lake env lean SampleSufficiency/Axioms.lean   # prints `#print axioms` for the headline theorems
```

Expected: each depends on `[propext, Classical.choice, Quot.sound]` and nothing else.
