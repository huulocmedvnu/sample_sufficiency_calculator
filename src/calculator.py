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


if __name__ == "__main__":
    # quick self-check with the calibrated Tahoe-100M constant (see calibrate.py / README)
    for m, lab in [(2.97, "Resveratrol (strong)"), (1.33, "weak signature")]:
        for theta in (0.05, 0.1, 0.2):
            n = calculate_experimental_cell_quota(7.66, 50, m, theta)
            print(f"{lab:22s} m={m:.2f} theta={theta} rad -> n*={n:,.0f} cells/arm")
