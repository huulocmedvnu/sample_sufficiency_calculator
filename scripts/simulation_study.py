#!/usr/bin/env python3
"""Monte Carlo simulation study with known truth (manuscript Section 3.2, Figure 3, Supplementary Note S7).

Four blocks, each varied one factor at a time around a base case (d = 50, the real Tahoe-100M within-condition
spectrum as a diagonal covariance, effect magnitude m = 2 x the control-pool floor, n_c = 3,100 control cells,
tolerance 0.1 rad, Gaussian cells):

  A  correctness   magnitude x control-pool size, tolerance, isotropic vs anisotropic covariance, effect direction
                   along the top / a middle / the bottom eigenvector. At n_t = quota the realized root-mean-square
                   angle is compared with the tolerance (target ratio 1); at the confidence quota (delta = 10%)
                   the exceedance probability is compared with delta.
  B  noise laws    cell-level Gaussian, multivariate t with 5 and 3 degrees of freedom (same covariance), a
                   two-component Gaussian mixture inside each condition (quota with the full mixture covariance and
                   with the within-component covariance only), and negative-binomial counts in gene space pushed
                   through the paper's own normalise -> log1p -> fixed-PCA projection (needs the cached Tahoe basis;
                   skipped, and recorded as skipped, when the cache is absent).
  C  plug-in       covariance and effect estimated from a pilot of 50, 200 or 1,000 cells per arm (sample covariance
                   or Ledoit-Wolf shrinkage), the plug-in quota, the realized angle at that quota, and how often a
                   truly finite condition is declared control-pool-limited (and vice versa) near the floor.
  D  baselines     the equal-arm, large-pool and isotropic rules of src/calculator.py at the same scenarios.

The quota, the floor and the confidence quota are obtained ONLY through src/calculator.py / src/engine.py; this
script simulates data and measures angles. Centroids of Gaussian cells are drawn from their exact sampling law
(the mean of n Gaussian cells is Gaussian with covariance Sigma/n); every other law is simulated cell by cell.
Fixed seeds, numpy only (scipy is imported lazily for the negative-binomial block), bit-for-bit reproducible for a
given seed whatever the number of worker processes.

    python3 scripts/simulation_study.py            # full run -> outputs/simulation_study.json,
                                                   #             fixtures/simulation_summary.json
    python3 scripts/simulation_study.py --quick    # CI mode: few replicates, no count-level block
"""
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):   # one BLAS thread per worker
    os.environ.setdefault(_v, "1")
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from calculator import (cell_quota, cell_quota_equal_arm, cell_quota_isotropic,  # noqa: E402
                        cell_quota_large_pool, cell_quota_report)
import engine  # noqa: E402

D = 50
THETA = 0.1                     # base tolerance (rad)
NC_BASE = 3100.0                # Tahoe median shared-vehicle pool
DELTA = 0.1                     # confidence level of the tail-controlled quota
K_BASE = 2.0                    # base magnitude, in units of the control-pool floor
K_GRID = [1.1, 1.25, 1.5, 2.0, 4.0, 10.0]
NC_GRID = [300.0, 1000.0, 3100.0, 10000.0, math.inf]
TOL_GRID = [0.05, 0.2, 0.3]     # plus the base 0.1
PILOT_GRID = [50, 200, 1000]
K_PILOT = [0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 4.0]
K_NOISE = [1.25, 2.0, 4.0, 10.0]
K_COUNTS = [1.25, 2.0, 4.0]
TAHOE_BASIS = "/mnt/hdd2/loc-tran/tahoe_work/out_dose/basis.npz"
TAHOE_OBS = "/mnt/hdd2/loc-tran/tahoe_work/meta/metadata/obs_metadata.parquet"

with open(os.path.join(ROOT, "fixtures", "chemical_within_cov.json")) as fh:
    ELL = np.array(json.load(fh)["tahoe"]["ell"], float)          # Tahoe per-PC within-condition variances
SIGMA_REAL = np.diag(ELL)
SIGMA_ISO = np.eye(D) * ELL.mean()


# ------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------
def unit(x):
    x = np.asarray(x, float)
    return x / np.linalg.norm(x)


def base_direction():
    """A fixed direction drawn uniformly on the sphere (seed 0): a generic drug displacement that does not
    concentrate on any eigen-axis."""
    return unit(np.random.default_rng(0).standard_normal(D))


def floor_of(u, Sigma, theta, n_c):
    """Control-pool floor m_min for direction u (from the engine, via the public report)."""
    return cell_quota_report(u, Sigma, theta, n_c)["m_min"]


def angles(vhat, u):
    c = (vhat @ u) / np.linalg.norm(vhat, axis=1)
    return np.arccos(np.clip(c, -1.0, 1.0))


def angle_stats(th, theta):
    """Realized RMS angle / tolerance, its Monte Carlo SE, and the exceedance P(theta > tolerance)."""
    t2 = th * th
    R = t2.size
    mean = float(t2.mean())
    se_mean = float(t2.std(ddof=1) / math.sqrt(R))
    rms = math.sqrt(mean)
    exc = float((th > theta).mean())
    return dict(reps=int(R), ratio=rms / theta, ratio_se=se_mean / (2.0 * rms) / theta,
                rms_rad=rms, exceed=exc, exceed_se=math.sqrt(exc * (1.0 - exc) / R))


def gauss_centroids(rng, R, n, L):
    """Exact sampling law of the centroid of n Gaussian cells with covariance L L^T."""
    if not math.isfinite(n):
        return np.zeros((R, D))
    return (rng.standard_normal((R, D)) @ L.T) / math.sqrt(n)


def cell_centroids(rng, R, n, law, L, mix=None):
    """Centroid of n cells drawn cell by cell under a non-Gaussian law (chunked over replicates)."""
    n = int(n)
    out = np.empty((R, D))
    step = max(1, int(4e7 // (n * D)))
    for a in range(0, R, step):
        b = min(R, a + step)
        z = rng.standard_normal((b - a, n, D)) @ L.T
        if law.startswith("t"):
            nu = float(law[1:])
            w = rng.chisquare(nu, (b - a, n, 1))
            z *= np.sqrt((nu - 2.0) / w)
        elif law == "mixture":
            pi, delta = mix
            z += (rng.random((b - a, n, 1)) < pi) * delta - pi * delta
        out[a:b] = z.mean(axis=1)
    return out


def ledoit_wolf(X):
    """Ledoit-Wolf (2004) shrinkage toward a scaled identity, the form used by scikit-learn
    (assume_centered: X holds residuals already centred within arm)."""
    n, d = X.shape
    emp = X.T @ X / n
    trace_vec = np.sum(X * X, axis=0) / n
    mu = trace_vec.sum() / d
    X2 = X * X
    beta_ = float(np.sum(X2.T @ X2))
    delta_ = float(np.sum((X.T @ X) ** 2)) / n ** 2
    beta = (beta_ / n - delta_) / (d * n)
    delta = (delta_ - 2.0 * mu * trace_vec.sum() + d * mu * mu) / d
    beta = min(beta, delta)
    shrink = 0.0 if beta <= 0 else beta / delta
    return (1.0 - shrink) * emp + shrink * mu * np.eye(d), shrink


# ------------------------------------------------------------------------------------------------
# block A / D: Gaussian, exact centroids
# ------------------------------------------------------------------------------------------------
def run_gauss(task):
    """One Gaussian scenario: simulate at n_t = the rule's quota (and at the confidence quota) and measure."""
    ss, rec, R = task["ss"], dict(task["rec"]), task["R"]
    rng = np.random.default_rng(ss)
    u, Sigma, theta, n_c, m = task["u"], task["Sigma"], task["theta"], task["n_c"], task["m"]
    v = m * u
    L = np.linalg.cholesky(Sigma)
    rep = cell_quota_report(v, Sigma, theta, n_c, confidence=DELTA)
    rec.update(m=m, m_floor=rep["m_min"], k=m / rep["m_min"] if rep["m_min"] > 0 else math.inf,
               n_two_arm=rep["required_cells_treated"], n_conf=rep["required_cells_treated_confident"],
               d_eff=rep["effective_dimensions"])
    n_rule = task.get("n_rule", rep["required_cells_treated"])
    rec["n_rule"] = n_rule
    if not math.isfinite(n_rule):
        rec.update(reps=0, ratio=None, ratio_se=None, exceed=None, exceed_se=None)
        return rec
    n_t = max(1, int(round(n_rule)))
    vhat = v + gauss_centroids(rng, R, n_t, L) - gauss_centroids(rng, R, n_c, L)
    rec.update(n_t=n_t, **angle_stats(angles(vhat, u), theta))
    if task.get("confidence", False) and math.isfinite(rep["required_cells_treated_confident"]):
        n_d = max(1, int(round(rep["required_cells_treated_confident"])))
        vhat = v + gauss_centroids(rng, R, n_d, L) - gauss_centroids(rng, R, n_c, L)
        st = angle_stats(angles(vhat, u), theta)
        rec.update(n_t_conf=n_d, ratio_conf=st["ratio"], exceed_conf=st["exceed"], exceed_conf_se=st["exceed_se"])
    return rec


def tasks_block_a(ss_list, R):
    u0 = base_direction()
    tasks = []
    k_iter = iter(ss_list)

    def add(rec, u, Sigma, theta, n_c, k, m_abs=None):
        fl = floor_of(u, Sigma, theta, n_c)
        m = m_abs if m_abs is not None else k * fl
        tasks.append(dict(ss=next(k_iter), rec=dict(block="A", **rec), R=R, u=u, Sigma=Sigma, theta=theta,
                          n_c=n_c, m=m, confidence=True))

    # magnitude x control-pool size (for the infinite pool the floor is zero: reuse the base-pool magnitudes)
    base_fl = floor_of(u0, SIGMA_REAL, THETA, NC_BASE)
    for n_c in NC_GRID:
        for k in K_GRID:
            if math.isfinite(n_c):
                add(dict(factor="magnitude x pool", n_c=n_c, k_nominal=k, tol=THETA, cov="real", dir="random"),
                    u0, SIGMA_REAL, THETA, n_c, k)
            else:
                add(dict(factor="magnitude x pool", n_c=n_c, k_nominal=k, tol=THETA, cov="real", dir="random"),
                    u0, SIGMA_REAL, THETA, n_c, k, m_abs=k * base_fl)
    for tol in TOL_GRID:
        add(dict(factor="tolerance", n_c=NC_BASE, k_nominal=K_BASE, tol=tol, cov="real", dir="random"),
            u0, SIGMA_REAL, tol, NC_BASE, K_BASE)
    add(dict(factor="covariance", n_c=NC_BASE, k_nominal=K_BASE, tol=THETA, cov="isotropic", dir="random"),
        u0, SIGMA_ISO, THETA, NC_BASE, K_BASE)
    order = np.argsort(-ELL)
    for name, idx in (("top", order[0]), ("middle", order[D // 2]), ("bottom", order[-1])):
        e = np.zeros(D); e[idx] = 1.0
        add(dict(factor="direction", n_c=NC_BASE, k_nominal=K_BASE, tol=THETA, cov="real", dir=name),
            e, SIGMA_REAL, THETA, NC_BASE, K_BASE)
    return tasks


def tasks_block_d(ss_list, R):
    u0 = base_direction()
    tasks = []
    k_iter = iter(ss_list)
    base_fl = floor_of(u0, SIGMA_REAL, THETA, NC_BASE)
    sigma2 = float(ELL.mean())
    for n_c in NC_GRID:
        fl = floor_of(u0, SIGMA_REAL, THETA, n_c) if math.isfinite(n_c) else base_fl
        for k in K_GRID:
            m = k * fl
            v = m * u0
            rules = {"equal-arm": cell_quota_equal_arm(v, SIGMA_REAL, THETA),
                     "large-pool": cell_quota_large_pool(v, SIGMA_REAL, THETA),
                     "isotropic": cell_quota_isotropic(sigma2, D, m, THETA, n_c)}
            for rule, n_rule in rules.items():
                tasks.append(dict(ss=next(k_iter), R=R, u=u0, Sigma=SIGMA_REAL, theta=THETA, n_c=n_c, m=m,
                                  n_rule=float(n_rule),
                                  rec=dict(block="D", rule=rule, n_c=n_c, k_nominal=k, tol=THETA)))
    return tasks


# ------------------------------------------------------------------------------------------------
# block B: cell-level noise laws
# ------------------------------------------------------------------------------------------------
def run_noise(task):
    ss, rec, R = task["ss"], dict(task["rec"]), task["R"]
    rng = np.random.default_rng(ss)
    u, theta, n_c, m, law = task["u"], task["theta"], task["n_c"], task["m"], task["law"]
    Sigma_gen, Sigma_quota, mix = task["Sigma_gen"], task["Sigma_quota"], task.get("mix")
    v = m * u
    L = np.linalg.cholesky(Sigma_gen)
    rep = cell_quota_report(v, Sigma_quota, theta, n_c, confidence=DELTA)
    rec.update(m=m, m_floor=rep["m_min"], k=m / rep["m_min"], n_two_arm=rep["required_cells_treated"],
               n_conf=rep["required_cells_treated_confident"])
    if not math.isfinite(rep["required_cells_treated"]):
        rec.update(reps=0, ratio=None)
        return rec
    n_t = max(1, int(round(rep["required_cells_treated"])))
    vhat = v + cell_centroids(rng, R, n_t, law, L, mix) - cell_centroids(rng, R, n_c, law, L, mix)
    rec.update(n_t=n_t, **angle_stats(angles(vhat, u), theta))
    if math.isfinite(rep["required_cells_treated_confident"]):
        n_d = max(1, int(round(rep["required_cells_treated_confident"])))
        vhat = v + cell_centroids(rng, R, n_d, law, L, mix) - cell_centroids(rng, R, n_c, law, L, mix)
        st = angle_stats(angles(vhat, u), theta)
        rec.update(n_t_conf=n_d, ratio_conf=st["ratio"], exceed_conf=st["exceed"], exceed_conf_se=st["exceed_se"])
    return rec


def tasks_block_b(ss_list, R):
    u0 = base_direction()
    k_iter = iter(ss_list)
    tasks = []
    # subpopulation shift for the mixture: a fixed direction, 30% of cells, |delta| = 6 (a strong within-condition
    # subpopulation such as a cell-cycle state); the mixture covariance is Sigma + pi (1 - pi) delta delta^T
    pi = 0.3
    delta = 6.0 * unit(np.random.default_rng(1).standard_normal(D))
    Sigma_mix = SIGMA_REAL + pi * (1.0 - pi) * np.outer(delta, delta)
    laws = [("gaussian", SIGMA_REAL, SIGMA_REAL, None), ("t5", SIGMA_REAL, SIGMA_REAL, None),
            ("t3", SIGMA_REAL, SIGMA_REAL, None),
            ("mixture (full covariance)", SIGMA_REAL, Sigma_mix, (pi, delta)),
            ("mixture (within-component covariance)", SIGMA_REAL, SIGMA_REAL, (pi, delta))]
    for k in K_NOISE:
        for name, Sg, Sq, mix in laws:
            law = "mixture" if name.startswith("mixture") else name
            # truth: the magnitude is k x the floor under the TRUE cell covariance (mixture: the full covariance)
            Strue = Sigma_mix if mix is not None else SIGMA_REAL
            m = k * floor_of(u0, Strue, THETA, NC_BASE)
            tasks.append(dict(ss=next(k_iter), R=R, u=u0, theta=THETA, n_c=NC_BASE, m=m, law=law,
                              Sigma_gen=Sg, Sigma_quota=Sq, mix=mix,
                              rec=dict(block="B", law=name, k_nominal=k, n_c=NC_BASE, tol=THETA)))
    return tasks


# ------------------------------------------------------------------------------------------------
# block C: plug-in estimation from a pilot
# ------------------------------------------------------------------------------------------------
def run_pilot(task):
    ss, rec, R = task["ss"], dict(task["rec"]), task["R"]
    rng = np.random.default_rng(ss)
    u, theta, n_c, m, n_p, est = task["u"], task["theta"], task["n_c"], task["m"], task["n_p"], task["est"]
    Sigma = task["Sigma"]
    v = m * u
    L = np.linalg.cholesky(Sigma)
    rep = cell_quota_report(v, Sigma, theta, n_c)
    n_true = rep["required_cells_treated"]
    rec.update(m=m, m_floor=rep["m_min"], k=m / rep["m_min"], n_two_arm=n_true, n_pilot=n_p, estimator=est)
    n_hat = np.empty(R); reg = np.empty(R, dtype=object); shrink = np.empty(R); rms_hat = np.full(R, np.nan)
    R_ang = task["R_ang"]
    step = max(1, int(2e7 // (2 * n_p * D)))
    for a in range(0, R, step):
        b = min(R, a + step)
        xt = v + rng.standard_normal((b - a, n_p, D)) @ L.T
        xc = rng.standard_normal((b - a, n_p, D)) @ L.T
        mt, mc = xt.mean(axis=1), xc.mean(axis=1)
        rt, rc = xt - mt[:, None, :], xc - mc[:, None, :]
        for i in range(b - a):
            res_mat = np.vstack([rt[i], rc[i]])
            if est == "sample":
                S_hat = res_mat.T @ res_mat / (2 * n_p - 2); sh = 0.0
            else:
                S_hat, sh = ledoit_wolf(res_mat)
            out = engine.compute(mt[i:i + 1], mc[i:i + 1], np.array([n_p]), np.array([n_p]), S_hat, theta=theta)
            reg[a + i] = out["regime"][0]
            shrink[a + i] = sh
            if out["regime"][0] == "NOT-DETECTABLE":
                n_hat[a + i] = math.nan
                continue
            # the plug-in quota for the PLANNED control pool, with the pilot's bias-corrected magnitude and the
            # pilot's perpendicular-noise trace (both from the engine)
            n_hat[a + i] = engine.quota_two_arm(out["trPSP"][0], np.square(out["m_corr"][0]), n_c, theta)
            if math.isfinite(n_hat[a + i]):
                n_t = max(1, int(round(n_hat[a + i])))
                vh = v + gauss_centroids(rng, R_ang, n_t, L) - gauss_centroids(rng, R_ang, n_c, L)
                rms_hat[a + i] = math.sqrt(float(np.mean(np.square(angles(vh, u)))))
    fin = np.isfinite(n_hat)
    nd = np.array([r == "NOT-DETECTABLE" for r in reg])
    pl = np.isinf(n_hat)
    q = lambda x, p: [float(np.percentile(x, pp)) for pp in p] if x.size else None
    ratio_fin = rms_hat[fin] / theta
    rec.update(reps=int(R), frac_not_detectable=float(nd.mean()), frac_pool_limited=float(pl.mean()),
               frac_finite=float(fin.mean()), mean_shrinkage=float(shrink.mean()),
               n_hat_over_true_q=(q(n_hat[fin] / n_true, [5, 25, 50, 75, 95]) if math.isfinite(n_true) else None),
               ratio_q=q(ratio_fin, [5, 25, 50, 75, 95]),
               ratio_pooled=(float(np.sqrt(np.mean(np.square(ratio_fin)))) if fin.any() else None),
               frac_ratio_within_1p05=(float(np.mean(ratio_fin <= 1.05)) if fin.any() else None),
               frac_ratio_within_1p2=(float(np.mean(ratio_fin <= 1.2)) if fin.any() else None),
               n_hat_median=(float(np.median(n_hat[fin])) if fin.any() else None))
    return rec


def tasks_block_c(ss_list, R, R_ang):
    u0 = base_direction()
    k_iter = iter(ss_list)
    tasks = []
    fl = floor_of(u0, SIGMA_REAL, THETA, NC_BASE)
    for n_p in PILOT_GRID:
        for est in ("sample", "ledoit-wolf"):
            for k in K_PILOT:
                tasks.append(dict(ss=next(k_iter), R=R, R_ang=R_ang, u=u0, theta=THETA, n_c=NC_BASE, m=k * fl,
                                  n_p=n_p, est=est, Sigma=SIGMA_REAL,
                                  rec=dict(block="C", k_nominal=k, n_c=NC_BASE, tol=THETA)))
    return tasks


# ------------------------------------------------------------------------------------------------
# block B (counts): negative-binomial counts in gene space through the paper's embedding
# ------------------------------------------------------------------------------------------------
_NB = {}


def _nb_init(comps, pmean, p_gene, lib_mu, lib_sd, phi, target):
    _NB.update(comps=comps, pmean=pmean, p=p_gene, lib_mu=lib_mu, lib_sd=lib_sd, phi=phi, target=target)


def _nb_embed(rng, n, fold):
    """n cells: library size ~ lognormal (fitted to Tahoe UMI totals), HVG counts ~ NB(mean = L p_g fold_g,
    dispersion phi), the remaining counts ~ Poisson, then normalise to the target sum, log1p, centre on the
    frozen PCA mean and project on the frozen components. Returns (n, D)."""
    p, comps, pmean, phi, target = _NB["p"], _NB["comps"], _NB["pmean"], _NB["phi"], _NB["target"]
    out = np.empty((n, D))
    rest = 1.0 - p.sum()
    r = 1.0 / phi
    for a in range(0, n, 2000):
        b = min(n, a + 2000)
        lib = rng.lognormal(_NB["lib_mu"], _NB["lib_sd"], b - a)
        mu = lib[:, None] * (p * fold)[None, :]
        cnt = rng.negative_binomial(r, r / (r + mu)).astype(np.float64)
        tot = cnt.sum(axis=1) + rng.poisson(lib * rest)
        x = np.log1p(cnt * (target / tot)[:, None])
        out[a:b] = (x - pmean) @ comps.T
    return out


def _nb_moments(args):
    """Sum and scatter of n embedded cells (for the truth: means and pooled within-arm covariance)."""
    ss, n, fold = args
    rng = np.random.default_rng(ss)
    x = _nb_embed(rng, n, fold)
    return x.sum(axis=0), x.T @ x, n


def _nb_replicates(args):
    """Centroid differences for a batch of replicates: n_t treated cells and n_c control cells each."""
    ss, reps, n_t, n_c, fold = args
    rng = np.random.default_rng(ss)
    out = np.empty((reps, D))
    for i in range(reps):
        out[i] = _nb_embed(rng, n_t, fold).mean(axis=0) - _nb_embed(rng, n_c, np.ones_like(fold)).mean(axis=0)
    return out


def _calibrate_gene_means(pmean, lib_mu, lib_sd, target):
    """Per-gene proportions p_g such that E[log1p(c target / L)] with c ~ Poisson(L p_g) matches the atlas's
    mean log-normalised expression (bisection on a log grid; Poisson pmf, library sizes on a 51-point grid)."""
    lg = np.exp(lib_mu + lib_sd * np.linspace(-2.5, 2.5, 51))
    lo, hi = np.full(pmean.size, -20.0), np.full(pmean.size, 0.0)
    ks = np.arange(0, 40)
    from scipy.special import gammaln
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lam = lg[:, None] * np.exp(mid)[None, :]                        # (lib, gene)
        logpmf = ks[:, None, None] * np.log(lam)[None] - lam[None] - gammaln(ks + 1)[:, None, None]
        val = np.log1p(ks[:, None, None] * (target / lg)[None, :, None])
        ex = (np.exp(logpmf) * val).sum(axis=0).mean(axis=0)            # E over k, then over library sizes
        hi = np.where(ex > pmean, mid, hi)
        lo = np.where(ex > pmean, lo, mid)
    return np.exp(0.5 * (lo + hi))


def run_counts_block(seed, R, workers, quick):
    if quick or not (os.path.exists(TAHOE_BASIS) and os.path.exists(TAHOE_OBS)):
        return dict(skipped=True, reason="quick mode" if quick else "Tahoe basis / metadata cache absent"), []
    import pyarrow.parquet as pq
    b = np.load(TAHOE_BASIS)
    comps, pmean, target = b["components"].astype(np.float64), b["pca_mean"].astype(np.float64), float(b["target_sum"])
    tc = pq.ParquetFile(TAHOE_OBS).read_row_group(0, columns=["tscp_count"]).to_pandas()["tscp_count"].values
    tc = tc[tc > 0].astype(np.float64)
    lib_mu, lib_sd = float(np.mean(np.log(tc))), float(np.std(np.log(tc)))
    phi = 0.3
    p_gene = _calibrate_gene_means(pmean, lib_mu, lib_sd, target)
    G = pmean.size
    rng0 = np.random.default_rng(seed)
    lfc = rng0.normal(0.0, 0.5, G)             # a global programme over all HVGs (log fold change sd 0.5 x scale)
    meta = dict(skipped=False, phi=phi, lib_mu=lib_mu, lib_sd=lib_sd, lib_median=float(np.exp(lib_mu)),
                hvg_count_fraction=float(p_gene.sum()), programme_genes=int(G), lfc_sd=0.5,
                truth_cells_per_arm=1_000_000, calibration_cells_per_arm=100_000)
    ss = np.random.SeedSequence(seed + 7)
    records = []
    with Pool(workers, initializer=_nb_init, initargs=(comps, pmean, p_gene, lib_mu, lib_sd, phi, target)) as pool:
        def moments(n_total, fold, chunk=10_000):
            parts = [(c, chunk, fold) for c in ss.spawn(n_total // chunk)]
            s = np.zeros(D); sc = np.zeros((D, D)); n = 0
            for a, b_, c in pool.map(_nb_moments, parts):
                s += a; sc += b_; n += c
            mean = s / n
            return mean, sc - n * np.outer(mean, mean), n

        # control truth and covariance (1e6 cells)
        mu_c, sc_c, n_c_truth = moments(1_000_000, np.ones(G))
        mu_c_cal, _, _ = moments(100_000, np.ones(G))

        def m_of(s_):
            """Embedded effect magnitude and floor for programme scale s_ (100,000 treated cells)."""
            mu_t_cal, sc_t_cal, n_cal = moments(100_000, np.exp(s_ * lfc))
            S_pool = (sc_c + sc_t_cal) / (n_c_truth + n_cal - 2)
            v_cal = mu_t_cal - mu_c_cal
            return float(np.linalg.norm(v_cal)), floor_of(unit(v_cal), S_pool, THETA, NC_BASE)

        for k in K_COUNTS:
            # calibrate the programme scale s so that m = k x floor: m(s) is close to linear in s, so two
            # proportional updates from s = 1 land within a few percent of the target
            s = 1.0
            for _ in range(3):
                m_s, fl = m_of(s)
                s *= k * fl / m_s
            fold = np.exp(s * lfc)
            mu_t, sc_t, n_t_truth = moments(1_000_000, fold)
            Sigma = (sc_c + sc_t) / (n_c_truth + n_t_truth - 2)
            v = mu_t - mu_c
            u = unit(v); m = float(np.linalg.norm(v))
            rep = cell_quota_report(v, Sigma, THETA, NC_BASE, confidence=DELTA)
            rec = dict(block="B", law="negative binomial counts", k_nominal=k, n_c=NC_BASE, tol=THETA,
                       programme_scale=s, m=m, m_floor=rep["m_min"], k=m / rep["m_min"],
                       n_two_arm=rep["required_cells_treated"], n_conf=rep["required_cells_treated_confident"],
                       sigma2_embedded=float(np.trace(Sigma) / D), d_eff=rep["effective_dimensions"],
                       sigma2_treated=float(np.trace(sc_t / (n_t_truth - 1)) / D),
                       sigma2_control=float(np.trace(sc_c / (n_c_truth - 1)) / D))
            for key, n_rule in (("", rep["required_cells_treated"]), ("_conf", rep["required_cells_treated_confident"])):
                if not math.isfinite(n_rule):
                    continue
                n_t = max(1, int(round(n_rule)))
                batch = 5
                parts = [(c, batch, n_t, int(NC_BASE), fold) for c in ss.spawn(R // batch)]
                vhat = np.vstack(pool.map(_nb_replicates, parts))          # centroid differences carry v already
                st = angle_stats(angles(vhat, u), THETA)
                if key == "":
                    rec.update(n_t=n_t, **st)
                else:
                    rec.update(n_t_conf=n_t, ratio_conf=st["ratio"], exceed_conf=st["exceed"],
                               exceed_conf_se=st["exceed_se"])
            records.append(rec)
            print(f"  [counts] k={k}: s={s:.3f} m={m:.3f} floor={rep['m_min']:.3f} n*={rep['required_cells_treated']:.0f} "
                  f"ratio={rec['ratio']:.4f} exceed_conf={rec.get('exceed_conf')}  sigma2_emb={rec['sigma2_embedded']:.3f}",
                  flush=True)
    return meta, records


# ------------------------------------------------------------------------------------------------
# summary for the figure and the text
# ------------------------------------------------------------------------------------------------
def summarize(records, meta):
    A = [r for r in records if r["block"] == "A"]
    B = [r for r in records if r["block"] == "B"]
    C = [r for r in records if r["block"] == "C"]
    Dd = [r for r in records if r["block"] == "D"]
    base = next(r for r in A if r["factor"] == "magnitude x pool" and r["n_c"] == NC_BASE and r["k_nominal"] == K_BASE)
    fin = lambda rs: [r for r in rs if r.get("ratio") is not None]
    worst = lambda rs: max(fin(rs), key=lambda r: abs(r["ratio"] - 1.0)) if fin(rs) else None
    pick = lambda r, keys: {k: r.get(k) for k in keys}
    akeys = ["factor", "n_c", "k_nominal", "k", "tol", "cov", "dir", "m", "m_floor", "n_two_arm", "n_t", "n_conf",
             "n_t_conf", "ratio", "ratio_se", "ratio_conf", "exceed", "exceed_conf", "exceed_conf_se", "reps", "d_eff"]
    bkeys = ["law", "k_nominal", "k", "m", "m_floor", "n_two_arm", "n_t", "n_conf", "n_t_conf", "ratio", "ratio_se",
             "ratio_conf", "exceed", "exceed_conf", "exceed_conf_se", "reps", "sigma2_embedded", "d_eff"]
    ckeys = ["n_pilot", "estimator", "k_nominal", "k", "n_two_arm", "frac_not_detectable", "frac_pool_limited",
             "frac_finite", "n_hat_over_true_q", "n_hat_median", "ratio_q", "ratio_pooled", "frac_ratio_within_1p05",
             "frac_ratio_within_1p2", "mean_shrinkage", "reps"]
    dkeys = ["rule", "n_c", "k_nominal", "k", "m", "n_two_arm", "n_rule", "n_t", "ratio", "ratio_se", "reps"]
    panel_a = {str(n_c): [pick(r, akeys) for r in A if r["factor"] == "magnitude x pool" and r["n_c"] == n_c]
               for n_c in NC_GRID}
    return dict(
        meta=meta,
        base_case=pick(base, akeys),
        block_a=dict(panel=panel_a, other=[pick(r, akeys) for r in A if r["factor"] != "magnitude x pool"],
                     worst=pick(worst(A), akeys), ratio_range=[min(r["ratio"] for r in fin(A)), max(r["ratio"] for r in fin(A))],
                     exceed_conf_max=max(r["exceed_conf"] for r in A if r.get("exceed_conf") is not None),
                     exceed_conf_min=min(r["exceed_conf"] for r in A if r.get("exceed_conf") is not None)),
        block_b=dict(rows=[pick(r, bkeys) for r in B], worst=pick(worst(B), bkeys) if fin(B) else None,
                     counts=meta.get("counts")),
        block_c=dict(rows=[pick(r, ckeys) for r in C]),
        block_d=dict(rows=[pick(r, dkeys) for r in Dd],
                     base_pool=[pick(r, dkeys) for r in Dd if r["n_c"] == NC_BASE]),
    )


def clean(o):
    """JSON-safe copy: numpy scalars to floats, infinities and NaN to strings."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, float):
        return ("inf" if o == math.inf else ("-inf" if o == -math.inf else ("nan" if math.isnan(o) else o)))
    if isinstance(o, (np.floating, np.integer)):
        return clean(float(o))
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=20261009)
    ap.add_argument("--quick", action="store_true", help="CI mode: few replicates, no count-level block")
    ap.add_argument("--reps", type=int, default=20000, help="replicates per Gaussian scenario (blocks A, D)")
    ap.add_argument("--reps-cells", type=int, default=4000, help="replicates per cell-level scenario (block B)")
    ap.add_argument("--reps-counts", type=int, default=2000, help="replicates per count-level scenario")
    ap.add_argument("--pilots", type=int, default=2000, help="pilot replicates per scenario (block C)")
    ap.add_argument("--angles-per-pilot", type=int, default=50)
    ap.add_argument("--workers", type=int, default=max(1, min(64, os.cpu_count() or 1)))
    ap.add_argument("--counts-only", action="store_true",
                    help="re-run only the count-level block and splice it into an existing --out file")
    ap.add_argument("--resummarize", action="store_true", help="rewrite --summary from an existing --out file")
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "simulation_study.json"))
    ap.add_argument("--summary", default=os.path.join(ROOT, "fixtures", "simulation_summary.json"))
    a = ap.parse_args()
    if a.quick:
        a.reps, a.reps_cells, a.pilots, a.angles_per_pilot = 400, 200, 100, 20
    t0 = time.time()

    def unclean(o):
        """Inverse of clean(): JSON 'inf'/'nan' strings back to floats so reloaded records compare as computed."""
        if isinstance(o, dict):
            return {k: unclean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [unclean(v) for v in o]
        return {"inf": math.inf, "-inf": -math.inf, "nan": math.nan}.get(o, o) if isinstance(o, str) else o

    if a.resummarize:
        with open(a.out) as fh:
            prev = unclean(json.load(fh))
        prev["meta"]["n_c_grid"] = [("inf" if not math.isfinite(x) else x) for x in NC_GRID]
        with open(a.summary, "w") as fh:
            json.dump(clean(summarize(prev["scenarios"], prev["meta"])), fh, indent=1)
        print(f"[simulation] rewrote {a.summary} from {a.out}")
        return
    if a.counts_only:
        with open(a.out) as fh:
            prev = unclean(json.load(fh))
        records = [r for r in prev["scenarios"] if r.get("law") != "negative binomial counts"]
        t_prev = prev["meta"]["runtime_s"]
    else:
        ss = np.random.SeedSequence(a.seed)
        children = ss.spawn(2000)
        tasks = (tasks_block_a(children[0:200], a.reps) + tasks_block_d(children[200:600], a.reps)
                 + tasks_block_b(children[600:800], a.reps_cells)
                 + tasks_block_c(children[800:1200], a.pilots, a.angles_per_pilot))
        print(f"[simulation] {len(tasks)} scenarios (A/D/B/C), workers={a.workers}, seed={a.seed}", flush=True)
        with Pool(a.workers) as pool:
            records = pool.map(_run_one, tasks, chunksize=1)
        t_prev = time.time() - t0
        print(f"[simulation] cell/centroid blocks done in {t_prev:.0f} s", flush=True)
    t1 = time.time()
    counts_meta, counts_records = run_counts_block(a.seed, a.reps_counts, a.workers, a.quick)
    records += counts_records
    runtime = t_prev + (time.time() - t1)
    counts_meta["runtime_s"] = round(time.time() - t1, 1)
    meta = dict(seed=a.seed, quick=a.quick, reps_gaussian=a.reps, reps_cells=a.reps_cells, reps_counts=a.reps_counts,
                pilots=a.pilots, angles_per_pilot=a.angles_per_pilot, scenarios=len(records), runtime_s=round(runtime, 1),
                d=D, theta_base=THETA, n_c_base=NC_BASE, delta=DELTA, k_base=K_BASE, k_grid=K_GRID,
                n_c_grid=[("inf" if not math.isfinite(x) else x) for x in NC_GRID], tolerance_grid=[THETA] + TOL_GRID,
                pilot_grid=PILOT_GRID, k_pilot=K_PILOT, k_noise=K_NOISE, k_counts=K_COUNTS,
                sigma2=float(ELL.mean()), tr_sigma=float(ELL.sum()), numpy=np.__version__, counts=counts_meta,
                date=time.strftime("%Y-%m-%d"))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(clean(dict(meta=meta, scenarios=records)), fh, indent=1)
    summary = clean(summarize(records, meta))
    with open(a.summary, "w") as fh:
        json.dump(summary, fh, indent=1)
    bc = summary["base_case"]
    print(f"[simulation] base case: n*={bc['n_two_arm']:.0f}, ratio={bc['ratio']:.4f} (SE {bc['ratio_se']:.4f}), "
          f"exceed at n*_delta={bc['exceed_conf']:.4f}; runtime {runtime:.0f} s; wrote {a.out} and {a.summary}")


def _run_one(t):
    return {"A": run_gauss, "D": run_gauss, "B": run_noise, "C": run_pilot}[t["rec"]["block"]](t)


if __name__ == "__main__":
    main()
