import Mathlib

/-!
# Part B — the normalization Jacobian  (Lemma 1 of `docs/THEORY.md`)

For the direction map `g x = x / ‖x‖` on a real inner-product space, the Fréchet derivative at
`v ≠ 0` is

```
Dg(v) = (1/‖v‖) · P ,   P = I − u uᵀ ,   u = v/‖v‖ ,
```

i.e. `Dg(v) y = ‖v‖⁻¹ · (y − ⟪u,y⟫ u)`.  Only the component of a perturbation **perpendicular** to
`v` changes the direction; the radial part lies in `ker Dg(v)`.  This is the analytic heart behind
the angular-error law (`E[θ²] = tr(PΣP)/m²`).

Mathlib provides the derivative of `‖·‖²`; we lift it to `‖·‖` through `√`, then to `‖·‖⁻¹`, then to
`g = ‖·‖⁻¹ • id` via the product rule.
-/

namespace SampleSufficiency

open scoped RealInnerProductSpace
open ContinuousLinearMap

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

/-- Fréchet derivative of the norm at `v ≠ 0`:  `D‖·‖(v) = ‖v‖⁻¹ ⟪v, ·⟫`. -/
theorem hasFDerivAt_norm' (v : E) (hv : v ≠ 0) :
    HasFDerivAt (fun x => ‖x‖) (‖v‖⁻¹ • innerSL ℝ v) v := by
  have hne : (‖v‖ : ℝ) ≠ 0 := norm_ne_zero_iff.mpr hv
  have h2 : ((‖v‖ : ℝ) ^ 2) ≠ 0 := pow_ne_zero 2 hne
  have hsqrt := (hasStrictFDerivAt_norm_sq v).hasFDerivAt.sqrt (by simpa using h2)
  simp only [Real.sqrt_sq_eq_abs, abs_norm] at hsqrt
  have hD : (‖v‖⁻¹ • innerSL ℝ v) = (1 / (2 * ‖v‖)) • (2 • innerSL ℝ v) := by
    ext y
    simp only [smul_apply, innerSL_apply_apply, smul_eq_mul, nsmul_eq_mul,
      Nat.cast_ofNat]
    field_simp
  rw [hD]; exact hsqrt

/-- **Lemma 1 (normalization Jacobian).**  For `v ≠ 0`, the direction map `g x = x/‖x‖` has

`Dg(v) = ‖v‖⁻¹ · (I − ‖v‖⁻² · v vᵀ) = (1/‖v‖) · P`,

with `P = I − u uᵀ` the tangent-space projector (`u = v/‖v‖`).  The radial direction `u` lies in
`ker Dg(v)` (see `hasFDerivAt_normalize_apply_self`). -/
theorem hasFDerivAt_normalize (v : E) (hv : v ≠ 0) :
    HasFDerivAt (fun x => ‖x‖⁻¹ • x)
      (‖v‖⁻¹ • (ContinuousLinearMap.id ℝ E
        - (‖v‖ ^ 2)⁻¹ • (innerSL ℝ v).smulRight v)) v := by
  have hne : (‖v‖ : ℝ) ≠ 0 := norm_ne_zero_iff.mpr hv
  have hnorm := hasFDerivAt_norm' v hv
  have hinv : HasFDerivAt (fun x => (‖x‖)⁻¹)
      (-(‖v‖ ^ 2)⁻¹ • (‖v‖⁻¹ • innerSL ℝ v)) v :=
    (hasDerivAt_inv hne).comp_hasFDerivAt v hnorm
  have hg := hinv.smul (hasFDerivAt_id v)
  have hD : (‖v‖⁻¹ • (ContinuousLinearMap.id ℝ E
        - (‖v‖ ^ 2)⁻¹ • (innerSL ℝ v).smulRight v))
      = (‖v‖⁻¹ • ContinuousLinearMap.id ℝ E)
        + (-(‖v‖ ^ 2)⁻¹ • (‖v‖⁻¹ • innerSL ℝ v)).smulRight v := by
    ext y
    simp only [add_apply, smul_apply,
      sub_apply, ContinuousLinearMap.id_apply,
      ContinuousLinearMap.smulRight_apply, innerSL_apply_apply]
    module
  rw [hD]
  exact hg

/-- The radial direction is in the kernel: `Dg(v) v = 0` (moving `v` along itself does not rotate
the direction).  Concretely `Dg(v) v = ‖v‖⁻¹ (v − ‖v‖⁻² ⟪v,v⟫ v) = 0`. -/
theorem hasFDerivAt_normalize_apply_self (v : E) (hv : v ≠ 0) :
    (‖v‖⁻¹ • (ContinuousLinearMap.id ℝ E
        - (‖v‖ ^ 2)⁻¹ • (innerSL ℝ v).smulRight v)) v = 0 := by
  have hne : (‖v‖ : ℝ) ≠ 0 := norm_ne_zero_iff.mpr hv
  simp only [smul_apply, sub_apply,
    ContinuousLinearMap.id_apply, ContinuousLinearMap.smulRight_apply, innerSL_apply_apply,
    real_inner_self_eq_norm_sq]
  rw [smul_smul, inv_mul_cancel₀ (pow_ne_zero 2 hne), one_smul, sub_self, smul_zero]

end SampleSufficiency
