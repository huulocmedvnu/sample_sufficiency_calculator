"""
Canonical within-condition per-cell variance sigma^2 — ONE definition for the whole repo.

Root cause of the 0.963/0.848/0.9174 leaks: three subtly different sigma^2 computations
coexisted (per-cell vs pseudobulk-moment, and different group keys). This module fixes ONE
definition and everything must call it:

    sigma^2 = mean_over_PCs(  sum_groups sum_cells (x - mean_group)^2  /  sum_groups (n_group - 1)  )

pooled WITHIN each (cell line x condition) group, over per-cell embedding coordinates.
This is the definition that reproduces the headline Tahoe-100M 0.9567 and EmeraldBay-full 0.9380.

The moment form (from streamed sufficient statistics sum, sumsq, count per group) is exact-equal
to the per-cell form IFF the group keys are the SAME (line x condition). The historical EmeraldBay
0.9174 differed only because it pooled a different grouping (cond_keys), not because moments differ.
"""
import numpy as np


def within_condition_sigma2(coords, group_keys, min_cells=2, return_per_pc=False):
    """Pooled within-(line x condition) per-cell variance, mean over PCs.

    coords     : (n_cells, d) per-cell embedding coordinates.
    group_keys : (n_cells,) label per cell = the (cell line x condition) group it belongs to.
    min_cells  : groups with fewer cells contribute no residual (need >=2 for a variance).
    """
    coords = np.asarray(coords, dtype=float)
    keys = np.asarray(group_keys)
    d = coords.shape[1]
    ss = np.zeros(d)
    dof = 0
    for g in np.unique(keys):
        Cg = coords[keys == g]
        if len(Cg) < min_cells:
            continue
        ss += ((Cg - Cg.mean(0)) ** 2).sum(0)
        dof += len(Cg) - 1
    if dof <= 0:
        raise ValueError("no group had >= min_cells cells")
    per_pc = ss / dof
    return (float(per_pc.mean()), per_pc) if return_per_pc else float(per_pc.mean())


def within_condition_sigma2_moments(sums, sumsq, counts, min_cells=2, return_per_pc=False):
    """Same statistic from streamed sufficient statistics (exact-equal to the per-cell form
    under the SAME group keys). sums, sumsq: (n_groups, d); counts: (n_groups,)."""
    sums = np.asarray(sums, float); sumsq = np.asarray(sumsq, float); counts = np.asarray(counts, float)
    ok = counts >= min_cells
    ss = (sumsq[ok] - sums[ok] ** 2 / counts[ok, None]).sum(0)
    dof = (counts[ok] - 1).sum()
    if dof <= 0:
        raise ValueError("no group had >= min_cells cells")
    per_pc = ss / dof
    return (float(per_pc.mean()), per_pc) if return_per_pc else float(per_pc.mean())
