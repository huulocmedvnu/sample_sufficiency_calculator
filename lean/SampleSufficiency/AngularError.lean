import Mathlib

/-!
# Formal verification of the deterministic core of the angular-error / cell-quota proof

This file machine-checks, in Lean 4 + Mathlib, the **linear-algebra backbone** of
`docs/THEORY.md` — the part on which the cell quota

```
n* = 2 · tr(P Σ P) / (m² θ²)
```

is built (`P = I - u uᵀ`, the tangent-space projector onto the directions perpendicular to the
unit signal `u`).  Concretely we prove:

* `proj_mul_self` / `proj_transpose` / `proj_mulVec_self` / `trace_proj`
  — `P` is a symmetric idempotent with `P u = 0` and `tr P = d - 1` (orthogonal projector of
  rank `d-1`);
* `trace_proj_conj` — the key identity `tr(P Σ P) = tr Σ - uᵀ Σ u`, i.e. variance *along* the
  signal is subtracted out and only perpendicular noise survives (THEORY.md §2, §3);
* `trace_proj_iso` — the isotropic reduction `tr(P (σ² I) P) = σ² (d-1)`, recovering the v1
  formula (THEORY.md §3 "isotropic reduction");
* `trace_proj_diag` — the PCA / diagonal form `tr(P diag(ℓ) P) = Σ_k (1 - u_k²) ℓ_k`
  (THEORY.md §7, §9).

The probabilistic statements (E[θ²], the Laurent–Massart tail, the §5 second-order term) are out
of scope for a Lean proof and are covered empirically by `tests/verify_theory.py`.
-/

namespace SampleSufficiency

open Matrix Finset
open scoped BigOperators

variable {d : ℕ}

/-- The tangent-space projector `P = I − u uᵀ`. -/
noncomputable def proj (u : Fin d → ℝ) : Matrix (Fin d) (Fin d) ℝ :=
  1 - Matrix.vecMulVec u u

variable (u : Fin d → ℝ)

/-- `tr (a bᵀ) = a ⬝ b`. -/
lemma trace_vecMulVec (a b : Fin d → ℝ) : (Matrix.vecMulVec a b).trace = a ⬝ᵥ b := by
  simp only [Matrix.trace, Matrix.diag_vecMulVec, Pi.mul_apply, dotProduct]

/-- `tr (M · (a bᵀ)) = b ⬝ M a`.  The one computational lemma everything else reuses. -/
lemma trace_mul_vecMulVec (M : Matrix (Fin d) (Fin d) ℝ) (a b : Fin d → ℝ) :
    (M * Matrix.vecMulVec a b).trace = b ⬝ᵥ M.mulVec a := by
  rw [Matrix.mul_vecMulVec, trace_vecMulVec, dotProduct_comm]

/-- For a unit vector, `(u uᵀ)(u uᵀ) = u uᵀ` (the outer product is idempotent). -/
lemma vecMulVec_mul_self (h : u ⬝ᵥ u = 1) :
    Matrix.vecMulVec u u * Matrix.vecMulVec u u = Matrix.vecMulVec u u := by
  ext i j
  simp only [Matrix.mul_apply, Matrix.vecMulVec_apply]
  have hfac : ∑ k, (u i * u k) * (u k * u j) = (u i * u j) * ∑ k, u k * u k := by
    rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun k _ => by ring
  have hsum : ∑ k, u k * u k = u ⬝ᵥ u := rfl
  rw [hfac, hsum, h, mul_one]

/-- `P` is symmetric. -/
lemma proj_transpose : (proj u)ᵀ = proj u := by
  unfold proj
  rw [Matrix.transpose_sub, Matrix.transpose_one, Matrix.transpose_vecMulVec]

/-- `P` is idempotent (an orthogonal projection), for a unit `u`. -/
lemma proj_mul_self (h : u ⬝ᵥ u = 1) : proj u * proj u = proj u := by
  unfold proj
  have hAA := vecMulVec_mul_self u h
  rw [sub_mul, one_mul, mul_sub, mul_one, hAA]
  abel

/-- `P u = 0`: the signal direction is killed by the perpendicular projector. -/
lemma proj_mulVec_self (h : u ⬝ᵥ u = 1) : (proj u).mulVec u = 0 := by
  unfold proj
  rw [Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.vecMulVec_mulVec, h]
  simp

/-- `tr P = d − ‖u‖²` (so `= d − 1` for a unit vector: rank-`(d-1)` projector). -/
lemma trace_proj : (proj u).trace = (d : ℝ) - u ⬝ᵥ u := by
  unfold proj
  rw [Matrix.trace_sub, Matrix.trace_one, trace_vecMulVec]
  simp [Fintype.card_fin]

/-- **The key identity** `tr(P S P) = tr S − uᵀ S u`.

Variance along the signal (`uᵀ S u`) is subtracted out; only the perpendicular trace remains.
This is the quantity the cell quota `n* = 2 tr(PSP)/(m²θ²)` is built from (THEORY.md §2–3). -/
theorem trace_proj_conj (S : Matrix (Fin d) (Fin d) ℝ) (h : u ⬝ᵥ u = 1) :
    (proj u * S * proj u).trace = S.trace - u ⬝ᵥ S.mulVec u := by
  unfold proj
  set A := Matrix.vecMulVec u u with hA
  have hAA : A * A = A := vecMulVec_mul_self u h
  have hX : (S * A).trace = u ⬝ᵥ S.mulVec u := trace_mul_vecMulVec S u u
  have hm : (1 - A) * S * (1 - A) = S - S * A - A * S + A * S * A := by noncomm_ring
  rw [hm, Matrix.trace_add, Matrix.trace_sub, Matrix.trace_sub]
  have h1 : (A * S).trace = (S * A).trace := Matrix.trace_mul_comm A S
  have h2 : (A * S * A).trace = (S * A).trace := by
    rw [Matrix.trace_mul_comm (A * S) A, ← Matrix.mul_assoc, hAA]
    exact Matrix.trace_mul_comm A S
  rw [h1, h2, hX]; ring

/-- **Isotropic reduction**: `tr(P (σ² I) P) = σ² (d − 1)` — recovers the v1 formula. -/
theorem trace_proj_iso (σ2 : ℝ) (h : u ⬝ᵥ u = 1) :
    (proj u * (σ2 • (1 : Matrix (Fin d) (Fin d) ℝ)) * proj u).trace = σ2 * ((d : ℝ) - 1) := by
  rw [trace_proj_conj u _ h, Matrix.trace_smul, Matrix.trace_one,
    Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, h, smul_eq_mul, smul_eq_mul,
    Fintype.card_fin]
  ring

/-- **PCA / diagonal form**: `tr(P diag(ℓ) P) = Σ_k (1 − u_k²) ℓ_k`. -/
theorem trace_proj_diag (ℓ : Fin d → ℝ) (h : u ⬝ᵥ u = 1) :
    (proj u * Matrix.diagonal ℓ * proj u).trace = ∑ k, (1 - (u k) ^ 2) * ℓ k := by
  rw [trace_proj_conj u _ h, Matrix.trace_diagonal]
  have hquad : u ⬝ᵥ (Matrix.diagonal ℓ).mulVec u = ∑ k, (u k) ^ 2 * ℓ k := by
    simp only [Matrix.mulVec_diagonal, dotProduct]
    exact Finset.sum_congr rfl fun k _ => by ring
  rw [hquad, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun k _ => by ring

end SampleSufficiency
