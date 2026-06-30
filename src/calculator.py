"""
calculator.py — sample-size (cell-quota) calculator for single-cell perturbation-direction screens.

Derivation (Delta method on the angular error of a perturbation vector)
----------------------------------------------------------------------
A drug's perturbation vector in a d-dimensional embedding (e.g. PCA) is

    v = mu_treated - mu_control

where each centroid mu is a MEAN over single cells. By the CLT, a centroid estimated from n cells
has covariance Sigma_cell / n, so for n cells per arm:

    Cov(v_hat) ~ Sigma_t/n + Sigma_c/n  ~  (2/n) * Sigma_cell          (similar per-cell covariance)

We care about the DIRECTION u = v/||v||. The angular error theta between the estimated and true
direction comes from the component of the noise PERPENDICULAR to v. With P = I - u u^T projecting
onto the (d-1)-dim tangent space, the perpendicular noise has

    E[|| P (v_hat - v) ||^2] = tr(P Cov(v_hat) P) ~ (2 sigma^2 / n) * (d - 1)

(isotropic approximation, sigma^2 = mean per-cell variance per embedding dim). Small-angle,
theta ~ ||perp noise|| / ||v||, so the RMS angular error is

    E[theta^2] ~ 2 (d - 1) sigma^2 / (n * m^2),     m = ||v|| = perturbation magnitude.

Setting the RMS angular error equal to a tolerance and solving for n gives the cell quota:

    n* = 2 (d - 1) sigma^2 / (m^2 * tolerance^2)        [cells PER ARM]

Interpretation: the "angular error ball" is the (d-1)-dim cap of RMS radius `tolerance` (radians)
around the true direction on the unit sphere; n* is the number of cells per arm needed to shrink
that ball to the requested radius. n* scales as sigma^2/m^2 — a STRONG perturbation (large m) needs
quadratically fewer cells than a weak one.

Caveats (honest): (1) isotropic-variance approximation — the exact form uses tr(P Sigma P) for the
specific v; when v aligns with high-variance dims the true n* differs. (2) `tolerance` is an RMS
angular SD in radians; sub-degree tolerances on noisy single-cell data demand very large n.
(3) Assumes the embedding (HVG + PCA) and normalization are fixed; changing them changes sigma^2.
"""
from __future__ import annotations


def calculate_experimental_cell_quota(single_cell_variance: float,
                                      num_dimensions: int,
                                      perturbation_magnitude: float,
                                      tolerance: float = 0.01) -> float:
    """Cells required PER ARM to estimate a perturbation's direction within `tolerance` radians (RMS).

    n* = 2 * (num_dimensions - 1) * single_cell_variance
            / (perturbation_magnitude**2 * tolerance**2)

    Parameters
    ----------
    single_cell_variance : float
        Mean per-cell variance per embedding dimension (sigma^2), in the SAME space as the magnitude
        (e.g. the shared PCA space). Must be > 0.
    num_dimensions : int
        Embedding dimensionality d (e.g. 50 PCA components). Must be >= 2.
    perturbation_magnitude : float
        m = ||v|| of the drug's perturbation vector in that space. Must be > 0.
    tolerance : float, default 0.01
        Target RMS angular error in radians (0.01 rad ~= 0.57 deg). NOTE: tight tolerances on noisy
        single-cell data require very large n; realistic values are often 0.05-0.2 rad.

    Returns
    -------
    float
        n* = cells required per arm (treated and control each). Round UP for planning.
    """
    if single_cell_variance <= 0:
        raise ValueError("single_cell_variance must be > 0")
    if num_dimensions < 2:
        raise ValueError("num_dimensions must be >= 2")
    if perturbation_magnitude <= 0:
        raise ValueError("perturbation_magnitude must be > 0")
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    return (2.0 * (num_dimensions - 1) * single_cell_variance
            / (perturbation_magnitude ** 2 * tolerance ** 2))


def rms_angular_error(single_cell_variance: float, num_dimensions: int,
                      perturbation_magnitude: float, n_cells_per_arm: int) -> float:
    """Inverse: the RMS angular error (radians) achieved with `n_cells_per_arm` cells per arm."""
    import math
    return math.sqrt(2.0 * (num_dimensions - 1) * single_cell_variance
                     / (n_cells_per_arm * perturbation_magnitude ** 2))


def calculate_cell_quota_anisotropic(covariance, perturbation_vector, tolerance=0.01,
                                     large_control_pool=False, confidence=None, hw_constant=2.0):
    """ANISOTROPIC cell quota (peer-review form) — no isotropic simplification. See docs/THEORY.md.

    Implements  n* = 2 * tr(P Sigma P) / (m^2 * tolerance^2),  P = I - u u^T,  u = v/||v||,
    where only the noise PERPENDICULAR to the signal rotates the direction. Variance ALONG the signal
    is subtracted out. Reduces to the isotropic 2*(d-1)*sigma^2/(m^2 t^2) when Sigma = sigma^2 I.

    Parameters
    ----------
    covariance : array
        Per-cell residual covariance Sigma in the SAME basis as `perturbation_vector`. Either a 1-D
        array of per-dimension variances (e.g. PCA explained-variances ell_k; Sigma=diag) or a full
        d x d matrix. Use WITHIN-condition (mean-centred per group) variance, not total.
    perturbation_vector : array
        v = mu_treated - mu_control (length d). m=||v||, u=v/m.
    tolerance : float
        Target RMS angular error (radians).
    large_control_pool : bool
        True -> shared/huge DMSO pool (control arm noiseless, factor 1). False -> equal arms (factor 2).
    confidence : float or None
        If set to delta in (0,1), also return a tail-controlled quota guaranteeing P(theta>tol)<=delta
        via a Hanson-Wright bound (generalized chi-square). `hw_constant` is the absolute constant.

    Returns
    -------
    dict: required_cells_per_arm (mean), [required_cells_per_arm_confident], perp_noise_trace,
          effective_dimensions, magnitude, isotropic_ratio (n*_aniso / n*_iso).
    """
    import math
    import numpy as np
    v = np.asarray(perturbation_vector, dtype=float)
    d = v.size
    m = float(np.linalg.norm(v))
    if m <= 0:
        raise ValueError("perturbation_vector must be nonzero")
    if tolerance <= 0:
        raise ValueError("tolerance must be > 0")
    u = v / m
    cov = np.asarray(covariance, dtype=float)
    Sigma = np.diag(cov) if cov.ndim == 1 else cov
    if Sigma.shape != (d, d):
        raise ValueError("covariance must be length-d (diagonal) or d x d")
    P = np.eye(d) - np.outer(u, u)
    M = P @ Sigma @ P                       # P Sigma P
    tr = float(np.trace(M))                 # = sum_k (1 - u_k^2) ell_k  for diagonal Sigma
    arm = 1.0 if large_control_pool else 2.0
    n_mean = arm * tr / (m ** 2 * tolerance ** 2)
    tr_sq = float(np.trace(M @ M))
    out = {
        "required_cells_per_arm": n_mean,
        "perp_noise_trace": tr,
        "effective_dimensions": (tr ** 2 / tr_sq) if tr_sq > 0 else 0.0,
        "magnitude": m,
        # ratio vs the isotropic estimate that uses the *mean* variance sigma^2 = tr(Sigma)/d
        "isotropic_ratio": tr / ((d - 1) * (np.trace(Sigma) / d)),
    }
    if confidence is not None:
        if not (0.0 < confidence < 1.0):
            raise ValueError("confidence (delta) must be in (0,1)")
        fro = math.sqrt(tr_sq)
        op = float(np.linalg.eigvalsh(M)[-1])
        L = math.log(1.0 / confidence)
        n_conf = arm * (tr + math.sqrt(hw_constant * L) * fro + hw_constant * L * op) \
            / (m ** 2 * tolerance ** 2)
        out["required_cells_per_arm_confident"] = n_conf
        out["confidence"] = confidence
    return out


# complexity exponents for the dry-lab cost model (cost ~ N**p, or N*log N for 'nlogn')
_COMPLEXITY = {"linear": 1.0, "ram": 1.0, "storage": 1.0, "pca": 1.0,
               "quadratic": 2.0, "pairwise": 2.0, "kernel": 2.0, "distance": 2.0}


def calculate_optimal_resource_allocation(single_cell_variance: float,
                                          num_dimensions: int,
                                          perturbation_magnitude: float,
                                          tolerance: float = 0.01,
                                          baseline_cells_per_well: float = 1394.0,
                                          complexity: str = "linear") -> dict:
    """Dual-sided sample-sufficiency: one info-saturation threshold n*, read two ways.

    WET-LAB (budget gating): n* = cells you must sequence per arm to resolve this drug's perturbation
    DIRECTION to `tolerance` radians (RMS). If n* < cells you currently run, you are over-sequencing
    and can multiplex more conditions per lane instead.

    DRY-LAB (safe downsampling): IF a well already has `baseline_cells_per_well` >= n*, you may
    downsample that well to n* cells before CENTROID/PSEUDOBULK analysis and still preserve its
    perturbation direction (hence the drug-drug similarity graph) to `tolerance`. The compute you
    save scales with the algorithm's complexity in N:
        linear (RAM, storage, streaming, PCA-fit):  saved = 1 - (n*/N0)
        quadratic (cell-cell pairwise / kernels):    saved = 1 - (n*/N0)**2
    If n* >= N0 the well is UNDER-sampled (no redundancy) -> ratio 0; do NOT downsample.

    SCOPE / honesty: n* governs centroid/direction-based (pseudobulk) analyses only. It does NOT
    license downsampling for CELL-RESOLUTION tasks (per-cell UMAP local structure, rare-population
    or cell-type detection) — those are governed by local density, not by this angular threshold, and
    can be distorted by downsampling. Runtime gains are POLYNOMIAL (linear-to-quadratic), never
    exponential.

    Parameters
    ----------
    single_cell_variance, num_dimensions, perturbation_magnitude, tolerance
        As in `calculate_experimental_cell_quota`.
    baseline_cells_per_well : float, default 1394
        Cells currently acquired per well (the over/under-sampling reference). Default = the Tahoe-100M
        excl3 atlas median (RESEARCH_LOG §27); override with your platform's number.
    complexity : str, default 'linear'
        Cost model for the dry-lab op: 'linear'/'ram'/'pca' (p=1), 'quadratic'/'pairwise' (p=2),
        or 'nlogn' (N log N).

    Returns
    -------
    dict with at least `required_cells_per_well` and `dry_lab_compute_reduction_ratio`, plus
    `regime`, `wet_lab_multiplex_gain`, `baseline_cells_per_well`, `complexity`.
    """
    import math
    n_star = calculate_experimental_cell_quota(single_cell_variance, num_dimensions,
                                               perturbation_magnitude, tolerance)
    N0 = float(baseline_cells_per_well)
    over_sampled = n_star < N0

    if not over_sampled:
        ratio = 0.0
        regime = "UNDER-sampled (n* >= N0): cells are not redundant — acquire MORE, do not downsample"
        multiplex = 1.0
    else:
        key = complexity.lower()
        if key == "nlogn":
            ratio = 1.0 - (n_star * math.log(max(n_star, 2))) / (N0 * math.log(max(N0, 2)))
        else:
            p = _COMPLEXITY.get(key)
            if p is None:
                raise ValueError(f"unknown complexity '{complexity}'; use {sorted(_COMPLEXITY)} or 'nlogn'")
            ratio = 1.0 - (n_star / N0) ** p
        regime = "over-sampled: downsampling to n* is safe for centroid/pseudobulk analysis"
        multiplex = N0 / n_star          # extra conditions per fixed per-lane read budget

    return {
        "required_cells_per_well": n_star,
        "dry_lab_compute_reduction_ratio": max(0.0, ratio),
        "regime": regime,
        "wet_lab_multiplex_gain": multiplex,
        "baseline_cells_per_well": N0,
        "complexity": complexity,
    }


if __name__ == "__main__":
    # quick self-check with the calibrated Tahoe-100M constant (see calibrate.py / README)
    for m, lab in [(2.97, "Resveratrol (strong)"), (1.33, "weak signature")]:
        for theta in (0.05, 0.1, 0.2):
            n = calculate_experimental_cell_quota(7.66, 50, m, theta)
            print(f"{lab:22s} m={m:.2f} theta={theta} rad -> n*={n:,.0f} cells/arm")
