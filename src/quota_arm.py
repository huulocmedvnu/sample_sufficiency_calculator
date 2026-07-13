"""
Cell quota with an EXPLICIT two-arm sampling model — one formula, no hard-coded factor 1 or 2.

The direction estimator v_hat = mu_t - mu_c has  Cov = Sigma/n_t + Sigma/n_c.
Requiring  tr(P Cov P)/m^2 = theta^2  and solving for the TREATED arm n_t:

    base := tr(P Sigma P) / (m^2 theta^2)          # large-pool / factor-1 core
    1/n_t = 1/base - 1/n_c
    n_t*  = 1/(1/base - 1/n_c)   if base < n_c      # resolvable at finite treated depth
          = +inf                 if base >= n_c      # CONTROL-POOL-LIMITED (given this n_c)

Limits handled automatically:
  * n_c = n_t (equal arms, matched vehicle)  ->  n_t* = 2*base   (the factor-2 result)
  * n_c -> inf (large pooled NTC)            ->  n_t* -> base     (the factor-1 result)

"Control-pool-limited" is a limit of THIS dataset as acquired (a design parameter n_c the
experimenter controls), NOT a physical limit: enlarging the NTC pool removes it.
"""
import numpy as np


def quota_two_arm(trPSP, m, theta, n_control):
    """Treated-arm quota n_t* (float, may be np.inf). Scalar or array inputs."""
    trPSP = np.asarray(trPSP, float); m = np.asarray(m, float)
    base = trPSP / (m ** 2 * theta ** 2)
    inv = 1.0 / base - 1.0 / n_control
    with np.errstate(divide="ignore"):
        n_t = np.where(inv > 0, 1.0 / np.where(inv > 0, inv, np.nan), np.inf)
    return n_t if n_t.ndim else float(n_t)


def quota_equal_arm(trPSP, m, theta):
    """Matched-vehicle (n_c = n_t) limit = 2*tr(PSP)/(m^2 theta^2)  (chemical convention)."""
    return 2.0 * np.asarray(trPSP, float) / (np.asarray(m, float) ** 2 * theta ** 2)


def m_min_resolvable(trPSP, n_control, theta):
    """Smallest magnitude resolvable at all given this control pool:
       m_min = sqrt(tr(PSP)/(n_c theta^2)).  m <= m_min  =>  n_t* = inf."""
    return np.sqrt(np.asarray(trPSP, float) / (n_control * theta ** 2))


def classify(n_t, N0, ghost=50000.0):
    """Four-class regime for a treated-arm quota n_t vs acquired depth N0.
       OVER: N0>=n_t (compress/multiplex) | UNDER: finite n_t>N0 (add treated cells) |
       GHOST: finite n_t>ghost (add treated, infeasible at routine depth) |
       POOL-LIMITED: n_t=inf (adding treated cells is useless; enlarge the NTC pool)."""
    if not np.isfinite(n_t):
        return "POOL-LIMITED"
    if N0 >= n_t:
        return "OVER"
    return "GHOST" if n_t > ghost else "UNDER"
