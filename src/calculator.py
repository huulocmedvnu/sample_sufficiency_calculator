"""calculator.py -- standalone cell-quota API for single-cell perturbation-direction screens.

THE quota formula lives in ONE place, src/engine.py (`quota_two_arm`); this module only assembles its
inputs (the perpendicular-noise trace tr(P Sigma P) and the squared magnitude) and calls it. It never
re-implements the formula -- so a public caller of this library gets exactly the number the paper's
engine produces.

Derivation (Delta method on the angular error). A perturbation vector v = mu_treated - mu_control in a
d-dim embedding is estimated from cell centroids; only noise PERPENDICULAR to v rotates the direction
u = v/||v||. With P = I - u u^T, the treated-arm cell quota to hold the RMS angular error at `tolerance`
is the TWO-ARM form

    n_t* = 1 / ( m^2 tolerance^2 / tr(P Sigma P)  -  1/n_c ) ,      m = ||v|| ,  n_c = control-pool size,

which returns +inf (CONTROL-POOL-LIMITED) when m is below the floor m_min = sqrt(tr(PSP)/(n_c tol^2)):
no treated depth resolves the direction. It reduces to the matched EQUAL-ARM special case
2 tr(PSP)/(m^2 tol^2) when n_c = n_t and to the LARGE-POOL limit tr(PSP)/(m^2 tol^2) as n_c -> inf.

**n_c is a REQUIRED argument of the recommended entry points.** Ignoring the control-pool size is
exactly the error this paper documents: against the largest atlas's real shared DMSO pool, ~55% of
conditions are control-pool-limited (unresolvable at any treated depth) -- invisible to any equal-arm
formula. Use `cell_quota(v, Sigma, tolerance, control_pool_size)`; the `_equal_arm` / `_large_pool`
variants are the two limits, provided only for the cases where they genuinely apply.
"""
from __future__ import annotations
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from engine import quota_two_arm, m_min_floor          # THE quota is defined in engine.py; we call it


def _perp(covariance, perturbation_vector):
    """Assemble (m, tr(P Sigma P), P Sigma P) from the covariance and the perturbation vector.
    This is linear algebra on the inputs (like estimating sigma^2), NOT the quota formula."""
    v = np.asarray(perturbation_vector, dtype=float)
    d = v.size
    m = float(np.linalg.norm(v))
    if m <= 0:
        raise ValueError("perturbation_vector must be nonzero")
    u = v / m
    cov = np.asarray(covariance, dtype=float)
    Sigma = np.diag(cov) if cov.ndim == 1 else cov
    if Sigma.shape != (d, d):
        raise ValueError("covariance must be length-d (diagonal) or d x d")
    P = np.eye(d) - np.outer(u, u)
    M = P @ Sigma @ P
    return m, float(np.trace(M)), M, Sigma


def cell_quota(perturbation_vector, covariance, tolerance, control_pool_size):
    """RECOMMENDED. Treated-arm cell quota n* for the GENERAL two-arm design.

    Parameters
    ----------
    perturbation_vector : array (length d)   v = mu_treated - mu_control; m = ||v||, u = v/m.
    covariance          : array              within-condition per-cell Sigma (length-d diagonal or d x d).
    tolerance           : float              target RMS angular error in radians (paper standard: 0.1).
    control_pool_size   : float  (REQUIRED)  the real number of control cells n_c. There is no default:
        the whole point of the paper is that ignoring n_c is wrong (55% of the largest atlas is
        control-pool-limited). Pass math.inf explicitly for a truly unbounded shared pool.

    Returns
    -------
    float : cells needed in the TREATED arm, or math.inf when the effect is below the control-pool floor
            m_min = sqrt(tr(P Sigma P)/(n_c tolerance^2)) -- i.e. CONTROL-POOL-LIMITED, unresolvable at any
            treated depth. (Enlarge n_c or loosen the tolerance; adding treated cells will not help.)
    """
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    if control_pool_size is None or (control_pool_size <= 0):
        raise ValueError("control_pool_size (n_c) is REQUIRED and must be > 0 (use math.inf for an "
                         "unbounded shared pool). Ignoring the control-pool size is the error this paper documents.")
    m, tr, _, _ = _perp(covariance, perturbation_vector)
    return quota_two_arm(tr, m ** 2, control_pool_size, tolerance)


def cell_quota_large_pool(perturbation_vector, covariance, tolerance):
    """Large shared control pool (n_c -> inf): n* = tr(P Sigma P)/(m^2 tolerance^2). Use ONLY when the
    control really is a large pool amortised across conditions (e.g. a big NTC or per-line reference);
    for a small shared vehicle use cell_quota(...) with the real n_c (it may be control-pool-limited)."""
    return cell_quota(perturbation_vector, covariance, tolerance, math.inf)


def cell_quota_equal_arm(perturbation_vector, covariance, tolerance):
    """Matched EQUAL-ARM special case (n_c = n_t): n* = 2 tr(P Sigma P)/(m^2 tolerance^2), i.e. twice the
    large-pool value. Valid ONLY when the control is a 1:1 co-plated vehicle of the SAME size as the
    treated arm. For a shared control pool this OVER-states resolvability -- use cell_quota(...)."""
    return 2.0 * cell_quota_large_pool(perturbation_vector, covariance, tolerance)


def cell_quota_isotropic(single_cell_variance, num_dimensions, magnitude, tolerance, control_pool_size):
    """Isotropic convenience: tr(P Sigma P) = (d-1) sigma^2 (no direction needed). n_c REQUIRED; returns
    math.inf when control-pool-limited. Use the anisotropic cell_quota(...) when Sigma is available."""
    if single_cell_variance <= 0:
        raise ValueError("single_cell_variance must be > 0")
    if num_dimensions < 2:
        raise ValueError("num_dimensions must be >= 2")
    if magnitude <= 0:
        raise ValueError("magnitude must be > 0")
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    if control_pool_size is None or control_pool_size <= 0:
        raise ValueError("control_pool_size (n_c) is REQUIRED and must be > 0 (use math.inf for a large pool)")
    tr = (num_dimensions - 1) * single_cell_variance
    return quota_two_arm(tr, magnitude ** 2, control_pool_size, tolerance)


def rms_angular_error(single_cell_variance, num_dimensions, perturbation_magnitude, n_cells_per_arm):
    """Inverse (equal-arm measurement law): the RMS angular error (radians) achieved with n cells per arm.
    This predicts theta given n (it does not solve for n), so it is the measurement direction, not a quota."""
    return math.sqrt(2.0 * (num_dimensions - 1) * single_cell_variance
                     / (n_cells_per_arm * perturbation_magnitude ** 2))


def cell_quota_report(perturbation_vector, covariance, tolerance, control_pool_size, confidence=None):
    """Full diagnostics for cell_quota. n_c REQUIRED. Returns a dict with the treated quota (or inf),
    whether it is control-pool-limited, the floor m_min, the perpendicular-noise trace, the effective
    noise dimension, and -- if `confidence`=delta -- a tail-controlled quota (Laurent-Massart) that
    guarantees P(theta>tolerance) <= delta, inflating the mean quota by a factor set by the noise shape."""
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    if control_pool_size is None or control_pool_size <= 0:
        raise ValueError("control_pool_size (n_c) is REQUIRED and must be > 0 (use math.inf for a large pool)")
    m, tr, M, Sigma = _perp(covariance, perturbation_vector)
    n_mean = quota_two_arm(tr, m ** 2, control_pool_size, tolerance)
    floor = float(m_min_floor(tr, control_pool_size, tolerance))
    tr_sq = float(np.trace(M @ M))
    out = {
        "required_cells_treated": n_mean,
        "control_pool_limited": not math.isfinite(n_mean),
        "m_min": floor,
        "perp_noise_trace": tr,
        "effective_dimensions": (tr ** 2 / tr_sq) if tr_sq > 0 else 0.0,
        "magnitude": m,
    }
    if confidence is not None:
        if not (0.0 < confidence < 1.0):
            raise ValueError("confidence (delta) must be in (0,1)")
        # Laurent-Massart (2000) tail for X = sum nu_i z_i^2 inflates the perpendicular-energy target from
        # tr to tr + 2||M||_F sqrt(L) + 2||M||_op L; scale the mean quota by that same ratio.
        fro = math.sqrt(tr_sq); op = float(np.linalg.eigvalsh(M)[-1]); L = math.log(1.0 / confidence)
        factor = (tr + 2.0 * math.sqrt(L) * fro + 2.0 * L * op) / tr if tr > 0 else float("inf")
        out["required_cells_treated_confident"] = n_mean * factor if math.isfinite(n_mean) else float("inf")
        out["confidence"] = confidence
    return out


# complexity exponents for the dry-lab cost model (cost ~ N**p, or N*log N for 'nlogn')
_COMPLEXITY = {"linear": 1.0, "ram": 1.0, "storage": 1.0, "pca": 1.0,
               "quadratic": 2.0, "pairwise": 2.0, "kernel": 2.0, "distance": 2.0}


def resource_allocation(perturbation_vector, covariance, tolerance, control_pool_size,
                        baseline_cells_per_well, complexity="linear"):
    """Dual-sided sample-sufficiency: one threshold n* read two ways (n_c REQUIRED).

    WET-LAB: n* treated cells resolve this direction to `tolerance`; if you run more you can multiplex.
    DRY-LAB: if a well already holds baseline_cells_per_well >= n*, downsampling it to n* preserves the
    direction (hence the drug-drug graph) for centroid/pseudobulk analysis; compute saved scales with the
    op's complexity in N. Scope: n* governs centroid/direction only, NOT cell-resolution structure; gains
    are polynomial, never exponential. If the condition is CONTROL-POOL-LIMITED (n* = inf) it is
    unresolvable at any treated depth -- neither downsampling nor multiplexing applies."""
    n_star = cell_quota(perturbation_vector, covariance, tolerance, control_pool_size)
    N0 = float(baseline_cells_per_well)
    if not math.isfinite(n_star):
        return {"required_cells_per_well": math.inf, "dry_lab_compute_reduction_ratio": 0.0,
                "regime": "CONTROL-POOL-LIMITED (n* = inf): unresolvable at any treated depth; enlarge n_c "
                          "or loosen tolerance", "wet_lab_multiplex_gain": 0.0,
                "baseline_cells_per_well": N0, "complexity": complexity}
    over_sampled = n_star < N0
    if not over_sampled:
        ratio, multiplex = 0.0, 1.0
        regime = "UNDER-sampled (n* >= N0): cells are not redundant -- acquire MORE, do not downsample"
    else:
        key = complexity.lower()
        if key == "nlogn":
            ratio = 1.0 - (n_star * math.log(max(n_star, 2))) / (N0 * math.log(max(N0, 2)))
        else:
            p = _COMPLEXITY.get(key)
            if p is None:
                raise ValueError(f"unknown complexity '{complexity}'; use {sorted(_COMPLEXITY)} or 'nlogn'")
            ratio = 1.0 - (n_star / N0) ** p
        multiplex = N0 / n_star
        regime = "over-sampled: downsampling to n* is safe for centroid/pseudobulk analysis"
    return {"required_cells_per_well": n_star, "dry_lab_compute_reduction_ratio": max(0.0, ratio),
            "regime": regime, "wet_lab_multiplex_gain": multiplex,
            "baseline_cells_per_well": N0, "complexity": complexity}


if __name__ == "__main__":
    # Tahoe-100M within-condition sigma^2 = 0.9567, d = 50 (SUPPLEMENT.md). n_c is required: a small
    # shared vehicle can be control-pool-limited (inf); a large pool is the finite large-pool quota.
    Sig = 0.9567 * np.ones(50)
    for m, lab in [(15.5, "Homoharringtonine 5uM (strong)"), (1.22, "median signature")]:
        v = np.zeros(50); v[0] = m
        big = cell_quota(v, Sig, 0.1, control_pool_size=1e9)
        small = cell_quota(v, Sig, 0.1, control_pool_size=3113)          # Tahoe's real shared DMSO pool
        print(f"{lab:30s} m={m:.2f}  large-pool n*={big:,.0f}  |  shared-pool(n_c=3113) n*={small}")
