"""
Unified sufficiency engine — ONE computation for every dataset (chemical AND genetic).

This is the single source of truth that fixes audit findings #1-#9 (docs/AUDIT_DENOMINATORS.md):
the calculation used to be duplicated across a chemical branch and a genetic branch that drifted
apart (genetic bias-corrected, chemical did not; genetic used snr>1.5, chemical a baseless 1.25;
genetic used the real control pool size, chemical hardcoded an equal-arm factor of 2; genetic used
the full two-arm quota, chemical hardcoded the equal-arm constant). The genetic branch was the only
one that was right. This module makes all four datasets go through its logic.

NO if-dataset, NO if-modality, NO hardcoded arm factor. Everything follows from the REAL n_c that the
pipeline supplies:
    sigma2      = tr(Sigma)/d                                   (within-condition residual variance)
    m_raw       = ||mu_t - mu_c||
    tr(S)       = tr(Sigma) * (1/n_t + 1/n_c)                   real n_c, no hardcoded factor
    m2_corr     = max(0, m_raw^2 - tr(S))                       bias-correction, every dataset
    detectable  = m_raw / sqrt(tr(S)) > 1.5                     one snr threshold, every dataset
    n*_t        = 1 / ( m2_corr*theta^2 / tr(P.Sigma.P) - 1/n_c )   full two-arm quota;
                  n_c ~ n_t  -> ~equal-arm (factor 2);  n_c >> n_t -> ~large-pool (factor 1),
                  automatically, with no manual regime switch.
    m_min       = sqrt( tr(P.Sigma.P) / (n_c*theta^2) )         control-pool floor, every dataset
    regime      = OVER / UNDER / GHOST / POOL-LIMITED / NOT-DETECTABLE
                  over the FULL condition set (no cell-count filter), for every dataset.

P = I - u u^T with u = (mu_t - mu_c)/m_raw, so tr(P.Sigma.P) = tr(Sigma) - u^T Sigma u.
"""
import numpy as np

THETA_DEFAULT = 0.1
GHOST = 50_000.0
SNR_DETECT = 1.5


def quota_two_arm(trPSP, m2, n_c, theta):
    """THE cell quota (the only place it is defined). Treated-arm cells n* to hold the RMS angular error
    at `theta`, given the perpendicular-noise trace tr(P.Sigma.P), the squared effect magnitude `m2`
    (bias-corrected where the caller has finite samples), the REAL control-pool size `n_c`, and `theta`:

        n_t* = 1 / ( m2 * theta^2 / tr(P.Sigma.P)  -  1/n_c )

    Returns +inf (CONTROL-POOL-LIMITED) when m2 <= tr(P.Sigma.P)/(n_c*theta^2), i.e. below the floor
    m_min: no treated depth resolves the direction. Reduces to the equal-arm value 2*tr/(m2*theta^2)
    when n_c = n_t and to the large-pool limit tr/(m2*theta^2) as n_c -> inf. Scalar or array inputs."""
    trPSP = np.asarray(trPSP, float); m2 = np.asarray(m2, float); n_c = np.asarray(n_c, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = m2 * theta ** 2 / trPSP - 1.0 / n_c
        out = np.where(denom > 0, 1.0 / denom, np.inf)
    return out if out.ndim else float(out)


def m_min_floor(trPSP, n_c, theta):
    """Smallest magnitude resolvable at all given this control pool: m_min = sqrt(tr(P.Sigma.P)/(n_c theta^2)).
    m <= m_min  =>  quota_two_arm returns +inf (control-pool-limited)."""
    return np.sqrt(np.asarray(trPSP, float) / (np.asarray(n_c, float) * theta ** 2))


def compute(mu_t, mu_c, n_t, n_c, Sigma, theta=THETA_DEFAULT, ghost=GHOST, snr_detect=SNR_DETECT):
    """Per-condition sufficiency spectrum. All inputs come from a pipeline; no dataset knowledge here.

    mu_t (N,d), mu_c (N,d): treated and control centroids in the d-dim embedding.
    n_t (N,), n_c (N,): REAL treated and control cell counts (n_c is the control POOL size).
    Sigma (d,d): within-condition covariance (diagonal for the chemical sufficient-stat pipeline,
                 full for the genetic pipeline; the engine only uses tr(Sigma) and u^T Sigma u).
    Returns a dict of per-condition arrays plus scalars.
    """
    mu_t = np.asarray(mu_t, float); mu_c = np.asarray(mu_c, float)
    n_t = np.asarray(n_t, float); n_c = np.asarray(n_c, float)
    Sigma = np.asarray(Sigma, float)
    d = Sigma.shape[0]
    trSig = float(np.trace(Sigma))
    sigma2 = trSig / d

    v = mu_t - mu_c
    m_raw = np.linalg.norm(v, axis=1)
    safe = m_raw > 0
    u = np.zeros_like(v); u[safe] = v[safe] / m_raw[safe, None]
    uSu = np.einsum('ij,jk,ik->i', u, Sigma, u)          # u^T Sigma u
    trPSP = trSig - uSu                                   # tr(P Sigma P), P = I - uu^T

    trS = trSig * (1.0 / n_t + 1.0 / n_c)                 # two-arm sampling floor, REAL n_c
    snr = np.divide(m_raw, np.sqrt(trS), out=np.zeros_like(m_raw), where=trS > 0)
    detectable = snr > snr_detect

    m2_corr = np.maximum(m_raw ** 2 - trS, 0.0)           # bias-correction (every dataset)
    m_corr = np.sqrt(m2_corr)

    nstar = quota_two_arm(trPSP, m2_corr, n_c, theta)     # THE quota (defined once, above)
    m_min = m_min_floor(trPSP, n_c, theta)                # control-pool floor

    # regime, precedence: not-detectable -> pool-limited -> ghost -> over -> under
    regime = np.empty(len(m_raw), dtype=object)
    pool_lim = ~np.isfinite(nstar)
    ghost_m = np.isfinite(nstar) & (nstar > ghost)
    over_m = np.isfinite(nstar) & (nstar < n_t)
    regime[:] = "UNDER"
    regime[over_m] = "OVER"
    regime[ghost_m] = "GHOST"
    regime[pool_lim] = "POOL-LIMITED"
    regime[~detectable] = "NOT-DETECTABLE"

    return dict(m_raw=m_raw, m_corr=m_corr, snr=snr, detectable=detectable,
                trPSP=trPSP, trS=trS, n_star=nstar, m_min=m_min, regime=regime,
                n_t=n_t, n_c=n_c, sigma2=sigma2, tr_Sigma=trSig, theta=theta)


def summarize(res):
    """Aggregate spectrum over ALL conditions (no cell-count filter). Percentages sum to 100 across
    the five regimes. Also reports the detectable-subset over/under split (the genetic-style report)."""
    reg = res["regime"]; N = len(reg); det = res["detectable"]
    finite = np.isfinite(res["n_star"])
    pc = lambda mask: round(100.0 * np.count_nonzero(mask) / N, 1) if N else 0.0
    counts = {k: pc(reg == k) for k in ("OVER", "UNDER", "GHOST", "POOL-LIMITED", "NOT-DETECTABLE")}
    ndet = int(det.sum())
    over_det = round(100.0 * np.count_nonzero((reg == "OVER") & det) / ndet, 1) if ndet else 0.0
    under_det = round(100.0 * np.count_nonzero((reg == "UNDER") & det) / ndet, 1) if ndet else 0.0
    med_nstar = float(np.median(res["n_star"][finite])) if finite.any() else float("inf")
    return dict(
        n_conditions=N, sigma2=round(res["sigma2"], 4),
        detectable_pct=round(100.0 * det.mean(), 1), not_detectable_pct=counts["NOT-DETECTABLE"],
        pct_over=counts["OVER"], pct_under=counts["UNDER"], pct_ghost=counts["GHOST"],
        pct_pool_limited=counts["POOL-LIMITED"],
        det_pct_over=over_det, det_pct_under=under_det,
        median_m_raw=round(float(np.median(res["m_raw"])), 3),
        median_m_corr=round(float(np.median(res["m_corr"])), 3),
        median_n_star=round(med_nstar) if np.isfinite(med_nstar) else None,
        median_n_t=int(np.median(res["n_t"])), median_n_c=int(np.median(res["n_c"])),
        median_m_min=round(float(np.median(res["m_min"])), 3),
        # effective quota constant for Table 1: large-pool bound (d-1)sigma^2/theta^2 (approached when
        # n_c >> n_t); equals half the equal-arm 2(d-1)sigma^2/theta^2 reached when n_c ~ n_t.
        C_largepool=round(_c_bound(res)),
    )


def _c_bound(res):
    d = int(round(res["tr_Sigma"] / res["sigma2"]))
    return (d - 1) * res["sigma2"] / res["theta"] ** 2
