"""pytest suite for the sample-sufficiency calculator.

Covers: closed-form value & scaling laws, anisotropic<->isotropic reduction, the perpendicular-only
property, resource-allocation regimes, input validation, and the formal/empirical theory checks
(symbolic Jacobian, Monte-Carlo E[theta^2], and Laurent-Massart tail-quota coverage).
"""
import numpy as np
import pytest

from calculator import (
    calculate_experimental_cell_quota as quota,
    calculate_cell_quota_anisotropic as aquota,
    calculate_optimal_resource_allocation as alloc,
    rms_angular_error,
)


# ---------------------------------------------------------------- isotropic core
def test_isotropic_closed_form():
    d, s2, m, th = 50, 2.406, 3.0, 0.1
    assert quota(s2, d, m, th) == pytest.approx(2 * (d - 1) * s2 / (m ** 2 * th ** 2))


def test_scaling_laws():
    base = quota(2.406, 50, 3.0, 0.1)
    assert quota(2.406, 50, 6.0, 0.1) == pytest.approx(base / 4)     # n* ~ 1/m^2
    assert quota(4.812, 50, 3.0, 0.1) == pytest.approx(base * 2)    # n* ~ sigma^2
    assert quota(2.406, 50, 3.0, 0.2) == pytest.approx(base / 4)     # n* ~ 1/theta^2


def test_rms_angular_error_is_inverse():
    n = quota(2.406, 50, 3.0, 0.1)
    assert rms_angular_error(2.406, 50, 3.0, n) == pytest.approx(0.1, rel=1e-9)


# ---------------------------------------------------------------- anisotropic
def test_anisotropic_reduces_to_isotropic():
    d, s2, m, th = 50, 2.406, 3.0, 0.1
    v = np.zeros(d); v[0] = m
    r = aquota(s2 * np.ones(d), v, th)
    assert r["required_cells_per_arm"] == pytest.approx(quota(s2, d, m, th))
    assert r["isotropic_ratio"] == pytest.approx(1.0)


def test_perpendicular_variance_only():
    """Variance ALONG the signal is irrelevant; perpendicular variance increases the quota."""
    d, m, th = 10, 2.0, 0.1
    v = np.zeros(d); v[0] = m
    ell = np.ones(d)
    base = aquota(ell, v, th)["required_cells_per_arm"]

    ell_par = ell.copy(); ell_par[0] = 100.0            # huge variance along the signal
    assert aquota(ell_par, v, th)["required_cells_per_arm"] == pytest.approx(base)

    ell_perp = ell.copy(); ell_perp[1] = 100.0          # huge variance perpendicular
    assert aquota(ell_perp, v, th)["required_cells_per_arm"] > 1.5 * base


def test_effective_dimensions_bounds():
    d = 20
    v = np.zeros(d); v[0] = 1.0
    r = aquota(np.ones(d), v, 0.1)
    # isotropic perpendicular noise -> d_eff = d-1 exactly
    assert r["effective_dimensions"] == pytest.approx(d - 1)


# ---------------------------------------------------------------- dual-sided allocation
def test_resource_allocation_under_sampled():
    weak = alloc(2.406, 50, 1.22, 0.1, baseline_cells_per_well=1296)
    assert weak["dry_lab_compute_reduction_ratio"] == 0.0
    assert weak["wet_lab_multiplex_gain"] == pytest.approx(1.0)
    assert "UNDER" in weak["regime"]


def test_resource_allocation_over_sampled():
    strong = alloc(2.406, 50, 12.07, 0.1, baseline_cells_per_well=1296, complexity="quadratic")
    lin = alloc(2.406, 50, 12.07, 0.1, baseline_cells_per_well=1296, complexity="linear")
    assert 0.0 < lin["dry_lab_compute_reduction_ratio"] < 1.0
    # quadratic-cost ops save at least as much as linear
    assert strong["dry_lab_compute_reduction_ratio"] >= lin["dry_lab_compute_reduction_ratio"]
    assert strong["wet_lab_multiplex_gain"] > 1.0


def test_input_validation():
    with pytest.raises(ValueError):
        quota(-1.0, 50, 3.0, 0.1)
    with pytest.raises(ValueError):
        quota(2.406, 1, 3.0, 0.1)
    with pytest.raises(ValueError):
        quota(2.406, 50, 0.0, 0.1)
    with pytest.raises(ValueError):
        aquota(np.ones(5), np.zeros(5), 0.1)            # zero perturbation vector
    with pytest.raises(ValueError):
        aquota(np.ones(5), np.ones(5), 0.1, confidence=1.5)


# ---------------------------------------------------------------- formal + empirical theory
def test_symbolic_jacobian():
    pytest.importorskip("sympy")
    from verify_theory import symbolic_jacobian_verification
    symbolic_jacobian_verification(d=3)                 # asserts internally


def test_monte_carlo_expected_angle():
    from verify_theory import monte_carlo_verification
    monte_carlo_verification(d=50, n=50_000, K=50_000, m=3.0, seed=0, tol=0.01)


def test_tail_quota_coverage():
    """Empirical P(theta > theta*) must not exceed delta when running the confidence quota
    (Laurent-Massart is an upper bound, so coverage should be comfortably below delta)."""
    rng = np.random.default_rng(1)
    d, m, theta_star, delta = 50, 3.0, 0.05, 0.10
    eig = 0.85 ** np.arange(d); eig = eig / eig.sum() * d
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    Sigma = (Q * eig) @ Q.T; Sigma = 0.5 * (Sigma + Sigma.T)
    v = rng.standard_normal(d); v = v / np.linalg.norm(v) * m

    n = aquota(Sigma, v, tolerance=theta_star, confidence=delta)["required_cells_per_arm_confident"]
    S = 2.0 * Sigma / n
    L = np.linalg.cholesky(S)
    K = 200_000
    v_hat = v + rng.standard_normal((K, d)) @ L.T
    cos = (v_hat @ v) / (m * np.linalg.norm(v_hat, axis=1))
    theta = np.arccos(np.clip(cos, -1.0, 1.0))
    coverage = float(np.mean(theta > theta_star))
    assert coverage <= delta, f"tail coverage {coverage:.3f} exceeds delta={delta}"
