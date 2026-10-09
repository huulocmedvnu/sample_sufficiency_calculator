"""The Monte Carlo simulation study (scripts/simulation_study.py, Supplementary Note S7) in its CI mode.

Runs the whole scenario grid with few replicates (no count-level block, which needs the Tahoe cache) and checks
the three facts the paper relies on: at the two-arm quota the realized root-mean-square angle sits at the
tolerance, at the confidence quota the exceedance stays below delta, and the plug-in block declares a truly
pool-limited base-pool condition pool-limited from a 1,000-cell pilot. Monte Carlo error with 400 replicates is
about 0.6% on the ratio, so the tolerances below are loose."""
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))


def test_simulation_study_quick(tmp_path):
    out, summ = tmp_path / "sim.json", tmp_path / "summary.json"
    cmd = [sys.executable, os.path.join(ROOT, "scripts", "simulation_study.py"), "--quick",
           "--workers", "2", "--out", str(out), "--summary", str(summ)]
    subprocess.run(cmd, check=True, cwd=ROOT, timeout=1800)
    S = json.load(open(summ))
    assert S["meta"]["quick"] and S["meta"]["counts"]["skipped"]

    # block A: realized RMS angle / tolerance within 5% of one in every finite scenario
    rows = [r for rows in S["block_a"]["panel"].values() for r in rows] + S["block_a"]["other"]
    ratios = np.array([r["ratio"] for r in rows if r["ratio"] is not None])
    assert ratios.size >= 30 and np.all(np.abs(ratios - 1.0) < 0.05), ratios

    # confidence quota (delta = 10%): exceedance below delta wherever the confidence quota is finite
    exc = [r["exceed_conf"] for r in rows if r.get("exceed_conf") is not None]
    assert exc and max(exc) <= S["meta"]["delta"], exc

    # block D: the equal-arm rule overspends and the large-pool rule undershoots at the base case
    base = {r["rule"]: r for r in S["block_d"]["base_pool"] if r["k_nominal"] == S["meta"]["k_base"]}
    assert base["equal-arm"]["ratio"] < 0.95 and base["large-pool"]["ratio"] > 1.05

    # block C: a 1,000-cell pilot declares a condition at 0.8 x the floor pool-limited almost always
    c = next(r for r in S["block_c"]["rows"] if r["n_pilot"] == 1000 and r["estimator"] == "sample"
             and r["k_nominal"] == 0.8)
    assert c["frac_pool_limited"] >= 0.9
