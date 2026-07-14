"""pytest suite for the standalone sample-sufficiency calculator (src/calculator.py).

Covers: the two-arm quota and its equal-arm / large-pool limits, the control-pool-limited (inf) case,
scaling laws, anisotropic<->isotropic reduction, the perpendicular-only property, resource-allocation
regimes, that n_c is REQUIRED, and the Laurent-Massart tail-quota coverage. The quota itself is defined
once in src/engine.py; these tests exercise the library that calls it.
"""
import math
import numpy as np
import pytest

from calculator import (
    cell_quota, cell_quota_equal_arm, cell_quota_large_pool, cell_quota_isotropic,
    cell_quota_report, resource_allocation, rms_angular_error,
)


# ---------------------------------------------------------------- limits & their relationship
def test_equal_arm_and_large_pool_closed_form():
    d, s2, m, th = 50, 2.406, 3.0, 0.1
    v = np.zeros(d); v[0] = m
    lp = cell_quota_large_pool(v, s2 * np.ones(d), th)
    eq = cell_quota_equal_arm(v, s2 * np.ones(d), th)
    assert lp == pytest.approx((d - 1) * s2 / (m ** 2 * th ** 2))     # large-pool = tr/(m^2 t^2)
    assert eq == pytest.approx(2 * lp)                                # equal-arm = twice large-pool
    # isotropic convenience with n_c=inf reproduces the large-pool value
    assert cell_quota_isotropic(s2, d, m, th, control_pool_size=math.inf) == pytest.approx(lp)


def test_two_arm_interpolates_and_pool_limits():
    d, s2, m, th = 50, 2.406, 3.0, 0.1
    v = np.zeros(d); v[0] = m
    ell = s2 * np.ones(d)
    lp = cell_quota_large_pool(v, ell, th)
    eq = cell_quota_equal_arm(v, ell, th)
    # n_c = n_t (=n*) gives the equal-arm value; a large n_c approaches the large-pool value
    assert cell_quota(v, ell, th, control_pool_size=eq) == pytest.approx(eq, rel=1e-6)
    assert cell_quota(v, ell, th, control_pool_size=1e9) == pytest.approx(lp, rel=1e-3)
    assert lp < cell_quota(v, ell, th, control_pool_size=5 * lp) < eq


def test_control_pool_limited_returns_inf():
    """Below the control-pool floor m_min the quota is +inf -- unresolvable at ANY treated depth."""
    d, s2, th = 50, 2.406, 0.1
    v = np.zeros(d); v[0] = 0.3                                       # weak effect
    r = cell_quota_report(v, s2 * np.ones(d), th, control_pool_size=200)
    assert r["control_pool_limited"] is True
    assert r["required_cells_treated"] == math.inf
    assert 0.3 < r["m_min"]                                           # effect sits below the floor
    assert cell_quota(v, s2 * np.ones(d), th, control_pool_size=200) == math.inf


def test_scaling_laws():
    d, th = 50, 0.1
    v = np.zeros(d); v[0] = 3.0
    base = cell_quota_large_pool(v, 2.406 * np.ones(d), th)
    v2 = np.zeros(d); v2[0] = 6.0
    assert cell_quota_large_pool(v2, 2.406 * np.ones(d), th) == pytest.approx(base / 4)   # n* ~ 1/m^2
    assert cell_quota_large_pool(v, 4.812 * np.ones(d), th) == pytest.approx(base * 2)    # n* ~ sigma^2
    assert cell_quota_large_pool(v, 2.406 * np.ones(d), 0.2) == pytest.approx(base / 4)   # n* ~ 1/theta^2


def test_rms_angular_error_is_inverse():
    d, s2, m = 50, 2.406, 3.0
    n = cell_quota_equal_arm(np.eye(d)[0] * m, s2 * np.ones(d), 0.1)
    assert rms_angular_error(s2, d, m, n) == pytest.approx(0.1, rel=1e-9)


# ---------------------------------------------------------------- anisotropic
def test_anisotropic_reduces_to_isotropic():
    d, s2, m, th = 50, 2.406, 3.0, 0.1
    v = np.zeros(d); v[0] = m
    assert cell_quota_large_pool(v, s2 * np.ones(d), th) == pytest.approx(
        cell_quota_isotropic(s2, d, m, th, control_pool_size=math.inf))


def test_perpendicular_variance_only():
    d, th = 10, 0.1
    v = np.zeros(d); v[0] = 2.0
    ell = np.ones(d)
    base = cell_quota_large_pool(v, ell, th)
    ell_par = ell.copy(); ell_par[0] = 100.0
    assert cell_quota_large_pool(v, ell_par, th) == pytest.approx(base)          # along-signal inert
    ell_perp = ell.copy(); ell_perp[1] = 100.0
    assert cell_quota_large_pool(v, ell_perp, th) > 1.5 * base                   # perpendicular bites


def test_effective_dimensions_bounds():
    d = 20
    v = np.zeros(d); v[0] = 1.0
    r = cell_quota_report(v, np.ones(d), 0.1, control_pool_size=math.inf)
    assert r["effective_dimensions"] == pytest.approx(d - 1)


# ---------------------------------------------------------------- dual-sided allocation
def test_resource_allocation_under_sampled():
    d = 50
    v = np.zeros(d); v[0] = 1.22
    weak = resource_allocation(v, 2.406 * np.ones(d), 0.1, control_pool_size=math.inf, baseline_cells_per_well=1296)
    assert weak["dry_lab_compute_reduction_ratio"] == 0.0
    assert weak["wet_lab_multiplex_gain"] == pytest.approx(1.0)
    assert "UNDER" in weak["regime"]


def test_resource_allocation_over_sampled():
    d = 50
    v = np.zeros(d); v[0] = 12.07
    strong = resource_allocation(v, 2.406 * np.ones(d), 0.1, control_pool_size=math.inf,
                                 baseline_cells_per_well=1296, complexity="quadratic")
    lin = resource_allocation(v, 2.406 * np.ones(d), 0.1, control_pool_size=math.inf,
                              baseline_cells_per_well=1296, complexity="linear")
    assert 0.0 < lin["dry_lab_compute_reduction_ratio"] < 1.0
    assert strong["dry_lab_compute_reduction_ratio"] >= lin["dry_lab_compute_reduction_ratio"]
    assert strong["wet_lab_multiplex_gain"] > 1.0


def test_n_c_is_required():
    d = 50
    v = np.zeros(d); v[0] = 3.0
    with pytest.raises(TypeError):
        cell_quota(v, np.ones(d), 0.1)                                # missing control_pool_size
    with pytest.raises(ValueError):
        cell_quota(v, np.ones(d), 0.1, control_pool_size=0)           # n_c must be > 0


def test_input_validation():
    d = 5
    with pytest.raises(ValueError):
        cell_quota(np.zeros(d), np.ones(d), 0.1, control_pool_size=1e6)   # zero perturbation
    with pytest.raises(ValueError):
        cell_quota_isotropic(-1.0, 50, 3.0, 0.1, control_pool_size=1e6)   # sigma2 <= 0
    with pytest.raises(ValueError):
        cell_quota_report(np.ones(d), np.ones(d), 0.1, control_pool_size=1e6, confidence=1.5)


# ---------------------------------------------------------------- formal + empirical theory
def test_symbolic_jacobian():
    pytest.importorskip("sympy")
    from verify_theory import symbolic_jacobian_verification
    symbolic_jacobian_verification(d=3)


def test_monte_carlo_expected_angle():
    from verify_theory import monte_carlo_verification
    monte_carlo_verification(d=50, n=50_000, K=50_000, m=3.0, seed=0, tol=0.01)


def test_tail_quota_coverage():
    rng = np.random.default_rng(1)
    d, m, theta_star, delta = 50, 3.0, 0.05, 0.10
    eig = 0.85 ** np.arange(d); eig = eig / eig.sum() * d
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    Sigma = (Q * eig) @ Q.T; Sigma = 0.5 * (Sigma + Sigma.T)
    v = rng.standard_normal(d); v = v / np.linalg.norm(v) * m
    n = cell_quota_report(v, Sigma, theta_star, control_pool_size=math.inf,
                          confidence=delta)["required_cells_treated_confident"]
    S = Sigma / n                      # large-pool design (control arm noiseless): Cov(v_hat) = Sigma/n
    L = np.linalg.cholesky(S)
    v_hat = v + rng.standard_normal((200_000, d)) @ L.T
    cos = (v_hat @ v) / (m * np.linalg.norm(v_hat, axis=1))
    coverage = float(np.mean(np.arccos(np.clip(cos, -1.0, 1.0)) > theta_star))
    assert coverage <= delta, f"tail coverage {coverage:.3f} exceeds delta={delta}"
