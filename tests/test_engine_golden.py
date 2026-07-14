"""
Golden lock on the unified engine (src/engine.py) and its six-screen spectrum after the round-2
rewrite (docs/AUDIT_DENOMINATORS.md #1-#9).

Three checks:
  1. The engine REPRODUCES the committed genetic ground truth from the genetic pipeline (the whole
     rewrite is validated against the branch that was always right). Skipped if raw data is absent.
  2. The committed fixtures/unified_spectrum.json headline + finding-#9 numbers are locked.
  3. Per-row self-consistency: the five regime percentages sum to 100 for every screen.
"""
import os, json
import numpy as np
import pytest

FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
SRC = os.path.join(os.path.dirname(__file__), "..", "src")
US = os.path.join(FX, "unified_spectrum.json")


def _by_name(js):
    return {s["name"]: s for s in js["headline"]}


# ---- 1. engine reproduces genetic ground truth ----
@pytest.mark.skipif(not os.path.exists("/mnt/hdd2/loc-tran/orion_work/full/HCT116/pseudobulk.npz"),
                    reason="genetic raw intermediates not present")
def test_engine_reproduces_genetic():
    import sys; sys.path.insert(0, SRC)
    from engine import compute, summarize
    from pipelines import genetic
    gold = {"orion_HCT116": (0.9133, 14.7, 0.0), "orion_HEK293T": (1.0334, 35.3, 0.1),
            "trade_jurkat": (1.4975, 70.6, 1.2), "trade_hepg2": (1.8667, 60.9, 4.4)}
    for k, (gs, gd, go) in gold.items():
        d = genetic.load(k); s = summarize(compute(d["mu_t"], d["mu_c"], d["n_t"], d["n_c"], d["Sigma"]))
        assert abs(s["sigma2"] - gs) < 5e-4, f"{k} sigma2 {s['sigma2']} != {gs}"
        assert s["detectable_pct"] == gd, f"{k} detectable {s['detectable_pct']} != {gd}"
        assert s["det_pct_over"] == go, f"{k} over(det) {s['det_pct_over']} != {go}"


# ---- 2. committed spectrum fixtures locked ----
def test_headline_spectrum_golden():
    js = json.load(open(US)); by = _by_name(js)
    gold = {  # detectable, over, under, ghost, pool-limited, C bound  (real diagonal Sigma)
        # HEADLINE: Tahoe on its real shared DMSO vehicle -> control-pool-limited majority
        "Tahoe-100M":   (94.2, 10.3, 13.8, 1.0, 69.2, 4688),
        "EmeraldBay":   (87.6, 2.7, 82.5, 2.0, 0.5, 4495),
        "Orion HCT116": (14.7, 0.0, 14.7, 0.0, 0.0, 4475),
        "Orion HEK293T": (35.3, 0.0, 35.3, 0.0, 0.0, 5063),
        "TRADE Jurkat": (70.6, 0.8, 68.3, 0.8, 0.7, 7338),
        "TRADE HepG2":  (60.9, 2.7, 56.7, 0.3, 1.2, 9147),
    }
    for name, (det, ov, un, gh, pl, C) in gold.items():
        s = by[name]
        assert (s["detectable_pct"], s["pct_over"], s["pct_under"], s["pct_ghost"], s["pct_pool_limited"]) \
            == (det, ov, un, gh, pl), f"{name} spectrum drifted: {s}"
        assert s["C_largepool"] == C, f"{name} C {s['C_largepool']} != {C}"


def test_tahoe_per_line_mean_sensitivity_golden():
    s = json.load(open(US))["sensitivity"]
    # referencing the per-line mean (large pool) instead of the real vehicle removes the pool limit
    assert s["pct_pool_limited"] == 0.0 and s["pct_over"] == 21.7 and s["detectable_pct"] == 97.0


# ---- 3. self-consistency ----
def test_regimes_sum_to_100():
    for s in json.load(open(US))["headline"]:
        tot = s["pct_over"] + s["pct_under"] + s["pct_ghost"] + s["pct_pool_limited"] + s["not_detectable_pct"]
        assert abs(tot - 100) < 0.2, f"{s['name']} regimes sum to {tot}"
        assert abs(s["detectable_pct"] + s["not_detectable_pct"] - 100) < 0.2, s["name"]
