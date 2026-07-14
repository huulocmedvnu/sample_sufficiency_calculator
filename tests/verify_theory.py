#!/usr/bin/env python
"""
verify_theory.py — dual-layer formal + empirical verification of the anisotropic angular-error theory
in docs/THEORY.md.

LAYER 1 (symbolic, SymPy): prove the Delta-method Jacobian of g(x)=x/||x|| equals (1/m)·P, with
    P = I - u u^T,  m = ||x||,  u = x/m  — by automatic differentiation, compared to the analytical form.

LAYER 2 (empirical, Monte Carlo): confirm the closed-form mean-squared angular error
    E[theta^2] = 2·tr(P Σ P) / (n·m^2)
against the EXACT geometric angle arccos(<v, v_hat>/(||v|| ||v_hat||)) sampled from N(v, 2Σ/n), under a
high-SNR regime, with an anisotropic Σ (exponentially decaying spectrum in a random basis).

Run:  python tests/verify_theory.py
Requires: numpy, sympy. (Optionally cross-checks the repo implementation if importable.)
"""
import numpy as np


# ====================================================================================
# LAYER 1 — SYMBOLIC VERIFICATION OF THE JACOBIAN  Dg(v) = (1/m) P
# ====================================================================================
def symbolic_jacobian_verification(d: int = 3) -> None:
    import sympy as sp
    print("=" * 78)
    print("LAYER 1 — Symbolic verification:  d/dx [ x/||x|| ]  ?=  (1/m)(I - u u^T)")
    print("=" * 78)

    # symbolic vector x = (x1, ..., xd)
    x = sp.Matrix(sp.symbols(f"x1:{d + 1}", real=True))
    norm = sp.sqrt((x.T * x)[0])                 # ||x|| = sqrt(sum x_i^2)
    g = x / norm                                 # normalization map g(x) = x/||x||

    # (1) SymPy automatic differentiation -> d x d Jacobian
    J_auto = g.jacobian(x)

    # (2) our analytical claim:  (1/m) * P,   m = ||x||,  u = x/m,  P = I - u u^T
    m = norm
    u = x / m
    P = sp.eye(d) - u * u.T
    J_claim = P / m

    # (3) elementwise symbolic equality (simplify each residual to 0)
    residual = sp.simplify(J_auto - J_claim)
    all_zero = all(sp.simplify(residual[i, j]) == 0 for i in range(d) for j in range(d))

    print("SymPy autodiff Jacobian g.jacobian(x):")
    sp.pprint(sp.simplify(J_auto))
    print("\nAnalytical (1/m)(I - u u^T):")
    sp.pprint(sp.simplify(J_claim))
    print(f"\nSymbolic residual (autodiff - analytical) simplifies to zero matrix:  {all_zero}")
    assert all_zero, f"SYMBOLIC MISMATCH — residual not identically zero:\n{residual}"

    # (4) independent numeric cross-check at a concrete non-zero v (exact rational arithmetic)
    v_subs = {x[0]: sp.Integer(3), x[1]: sp.Integer(-1), x[2]: sp.Integer(2)} if d == 3 else \
             {x[k]: sp.Integer(k + 1) for k in range(d)}
    J_auto_num = np.array(J_auto.subs(v_subs).evalf(), dtype=float)
    # analytical built directly in numpy from the same v
    vv = np.array([float(v_subs[x[k]]) for k in range(d)])
    mm = np.linalg.norm(vv); uu = vv / mm
    J_claim_num = (np.eye(d) - np.outer(uu, uu)) / mm
    max_abs = np.max(np.abs(J_auto_num - J_claim_num))
    print(f"Numeric cross-check at v={vv.tolist()}:  max|autodiff - analytical| = {max_abs:.2e}")
    assert max_abs < 1e-12, "NUMERIC MISMATCH at the evaluated vector"
    print(">>> LAYER 1 PASSED: the Jacobian is exactly (1/m)(I - u u^T).\n")


# ====================================================================================
# LAYER 2 — EMPIRICAL MONTE-CARLO VERIFICATION OF  E[theta^2] = 2 tr(P Σ P)/(n m^2)
# ====================================================================================
def monte_carlo_verification(d: int = 50, n: int = 50_000, K: int = 50_000,
                             m: float = 3.0, seed: int = 0, tol: float = 0.01) -> None:
    print("=" * 78)
    print("LAYER 2 — Monte-Carlo verification of  E[theta^2] = 2 tr(P Σ P)/(n m^2)")
    print("=" * 78)
    rng = np.random.default_rng(seed)

    # --- anisotropic covariance: exponentially decaying eigenvalues in a RANDOM basis ---
    # (random rotation => Σ is non-diagonal, testing the coordinate-free trace formula)
    eig = 0.85 ** np.arange(d)                       # dominant early "PCs", long tail
    eig = eig / eig.sum() * d                        # scale so tr(Σ) = d (mean variance 1)
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))  # random orthonormal basis
    Sigma = (Q * eig) @ Q.T
    Sigma = 0.5 * (Sigma + Sigma.T)                  # symmetrize against float drift

    # --- ground-truth perturbation vector v (random direction, fixed magnitude m) ---
    v = rng.standard_normal(d)
    v = v / np.linalg.norm(v) * m
    u = v / m
    P = np.eye(d) - np.outer(u, u)                   # tangent-space projector

    # --- sampling covariance of v_hat for EQUAL arms: S = 2 Σ / n ---
    S = 2.0 * Sigma / n

    # --- analytical prediction (first-order Delta method) ---
    theo = float(np.trace(P @ Sigma @ P)) * 2.0 / (n * m ** 2)   # = tr(P S P)/m^2

    # validity / interpretability diagnostics
    uSu = float(u @ S @ u)                            # noise variance ALONG the signal
    rho2 = m ** 2 / uSu                               # along-signal SNR^2 (need >> 1)
    PSP = P @ S @ P
    d_eff = float(np.trace(PSP) ** 2 / np.trace(PSP @ PSP))      # effective noise dof

    # --- brute-force Monte Carlo: exact geometric angle ---
    L = np.linalg.cholesky(S)                         # S = L L^T
    z = rng.standard_normal((K, d))
    e = z @ L.T                                       # e ~ N(0, S), shape (K, d)
    v_hat = v + e                                     # perturbed estimates
    cos = (v_hat @ v) / (m * np.linalg.norm(v_hat, axis=1))
    theta = np.arccos(np.clip(cos, -1.0, 1.0))        # EXACT angle, clipped for float safety
    emp = float(np.mean(theta ** 2))

    rel_err = abs(emp - theo) / theo
    mc_se = float(np.std(theta ** 2, ddof=1) / np.sqrt(K)) / theo   # MC relative std-error of the mean

    print(f"d={d}, n={n:,} cells/arm, K={K:,} MC iters, m={m}, tr(Σ)={np.trace(Sigma):.1f} "
          f"(anisotropic, random basis)")
    print(f"validity: along-signal SNR rho^2 = {rho2:,.0f} (>>1 required); effective noise dof d_eff = {d_eff:.1f}")
    print(f"RMS angle ~ {np.degrees(np.sqrt(theo)):.3f} deg  (small-angle / high-SNR regime)\n")
    print(f"  Analytical  E[theta^2] = {theo:.6e}")
    print(f"  Monte-Carlo E[theta^2] = {emp:.6e}   (MC rel. std-error {mc_se*100:.3f}%)")
    print(f"  relative error         = {rel_err*100:.4f}%   (threshold {tol*100:.1f}%)")
    assert rel_err < tol, f"EMPIRICAL MISMATCH: relative error {rel_err*100:.3f}% exceeds {tol*100:.1f}%"
    print(">>> LAYER 2 PASSED: closed-form E[theta^2] matches the exact geometric angle.\n")

    # --- (optional) verify the SHIPPED implementation equals the analytical n* ---
    try:
        import os, sys
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
        from calculator import cell_quota_equal_arm
        theta_star = float(np.sqrt(theo))                       # tolerance s.t. n* should return n
        n_calc = cell_quota_equal_arm(v, Sigma, theta_star)     # equal-arm form (matches E[theta^2] above)
        impl_rel = abs(n_calc - n) / n
        print(f"LAYER 3 (implementation cross-check): calculator returns n*={n_calc:.1f} "
              f"for tolerance=sqrt(E[theta^2]); target n={n} -> rel. error {impl_rel*100:.4f}%")
        assert impl_rel < 1e-6, "IMPLEMENTATION MISMATCH vs analytical n*"
        print(">>> LAYER 3 PASSED: src/calculator.py matches the analytical formula exactly.\n")
    except ImportError:
        print("(repo calculator not importable — skipping implementation cross-check)\n")


# ====================================================================================
# LAYER 4 — SECOND-ORDER TERM (THEORY.md §5, eq. 6): arc-vs-tangent must be kept
# ====================================================================================
def second_order_verification(d: int = 6, m: float = 4.0, denom: float = 50.0,
                              K: int = 6_000_000, seed: int = 0) -> None:
    """Verify the §5 second-order expansion of E[theta^2] (eq. 6) and show that the
    lever-arm-only form (= the expansion of E[tan^2 theta]) over-predicts at this order.

        E[theta^2] = tr(PSP)/m^2 * (1 + (1/m^2)[ 3 uSu + 6 uSPSu/tr(PSP)
                                                  - (2/3) tr(PSP) - (4/3) tr((PSP)^2)/tr(PSP) ])
    The last two bracket terms are the arc-vs-tangent correction theta^2 = tan^2 - (2/3) tan^4 + ...
    """
    print("=" * 78)
    print("LAYER 4 — second-order term  E[theta^2]  (THEORY.md §5, eq. 6)")
    print("=" * 78)
    rng = np.random.default_rng(seed)
    # anisotropic SPD Sigma (random basis), unit direction u, magnitude-m signal v
    A = rng.standard_normal((d, d))
    Sigma = A @ A.T
    u = rng.standard_normal(d); u /= np.linalg.norm(u)
    v = m * u
    P = np.eye(d) - np.outer(u, u)
    S = Sigma / denom                                   # sampling covariance of v_hat

    trPSP = float(np.trace(P @ S @ P))
    uSu = float(u @ S @ u)
    uSPSu = float(u @ S @ P @ S @ u)
    trPSP2 = float(np.trace(P @ S @ P @ S))             # tr((PSP)^2)  (P idempotent)
    lead = trPSP / m ** 2
    # lever-arm only  == expansion of E[tan^2 theta]  (the pre-fix THEORY.md formula)
    lever_only = lead * (1.0 + (3.0 / m ** 2) * (uSu + 2.0 * uSPSu / trPSP))
    # full eq. (6): add the arc-vs-tangent term  -(2/3) E||Pe||^4 / m^4
    corrected = (lead
                 + (1.0 / m ** 4) * (3.0 * (trPSP * uSu + 2.0 * uSPSu))
                 - (2.0 / 3.0) * (trPSP ** 2 + 2.0 * trPSP2) / m ** 4)

    L = np.linalg.cholesky(S)
    e = rng.standard_normal((K, d)) @ L.T
    v_hat = v + e
    cos = (v_hat @ v) / (m * np.linalg.norm(v_hat, axis=1))
    theta = np.arccos(np.clip(cos, -1.0, 1.0))
    emp = float(np.mean(theta ** 2))
    se = float(np.std(theta ** 2, ddof=1) / np.sqrt(K)) / emp

    re_lead = abs(emp - lead) / emp
    re_lever = abs(emp - lever_only) / emp
    re_corr = abs(emp - corrected) / emp
    rms_deg = np.degrees(np.sqrt(emp))
    print(f"d={d}, m={m}, S=Sigma/{denom:.0f}, K={K:,}, RMS angle ~ {rms_deg:.2f} deg "
          f"(MC rel. std-error {se*100:.3f}%)")
    print(f"  MC      E[theta^2] = {emp:.6e}")
    print(f"  leading order      -> rel err {re_lead*100:6.3f}%")
    print(f"  lever-arm only (E[tan^2]) -> rel err {re_lever*100:6.3f}%   (pre-fix §5 form)")
    print(f"  eq.(6) corrected   -> rel err {re_corr*100:6.3f}%")
    assert re_corr < 0.5e-2, f"eq.(6) mismatch: {re_corr*100:.3f}% (> 0.5%)"
    assert re_lever > 1.5e-2, "expected lever-arm-only form to be materially worse (sanity check)"
    assert re_corr < re_lead < re_lever, "ordering: corrected < leading < lever-only expected"
    print(">>> LAYER 4 PASSED: arc-vs-tangent term is required; eq.(6) matches MC to <0.5%.\n")


def main():
    symbolic_jacobian_verification(d=3)
    monte_carlo_verification(d=50, n=50_000, K=50_000, m=3.0, seed=0, tol=0.01)
    second_order_verification(d=6, m=4.0, denom=50.0, K=6_000_000, seed=0)
    print("=" * 78)
    print("ALL VERIFICATION LAYERS PASSED.")
    print("=" * 78)


if __name__ == "__main__":
    main()
