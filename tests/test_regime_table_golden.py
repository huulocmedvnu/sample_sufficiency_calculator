"""
Golden lock on the sufficiency-regime spectrum (Table 5 / Table 1) after the audit-round-1
denominator + unit fixes (docs/AUDIT_DENOMINATORS.md #1-#7).

Pins three things so they cannot silently drift apart again:
  1. The committed aggregate spectrum numbers for every screen (the Table 5 / Table 1 cells).
  2. Internal self-consistency of each row (percent blocks sum to 100; median n* = C / median m^2;
     the EmeraldBay aggregate reproduces from its own per-group m/N under the isotropic rule).
  3. The manuscript text agrees with the fixtures for the two chemical rows -- numbers are the
     source of truth, prose is made to match (not the other way round).

Tahoe MUST stay 10.6 / 89.0 / 0.4 and EmeraldBay MUST be 1.0 / 98.8 / 0.2 on the
(sample x line) = (drug x dose x line) unit with no cell-count filter, isotropic C = 9,192/m^2.
"""
import json, os, re
import numpy as np
import pytest

FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
DOC = os.path.join(os.path.dirname(__file__), "..", "docs", "MANUSCRIPT_DEEPSEEK.md")
d = 50


def _load(n):
    return json.load(open(os.path.join(FX, n)))


# ---- golden aggregate spectrum, one row per screen (over / under / ghost / detectable) ----
def test_tahoe_spectrum_golden():
    c = _load("tahoe_constants.json")
    assert (c["pct_OVER"], c["pct_UNDER"], c["pct_Ghost"]) == (10.6, 89.0, 0.4), \
        "Tahoe spectrum must stay 10.6 / 89.0 / 0.4 (stop condition b)"


def test_emeraldbay_spectrum_golden():
    e = _load("emeraldbay_spectrum.json")
    assert e["n_groups"] == 4912, "EmeraldBay unit is (sample x line), 4,912 non-empty groups"
    assert e["n_lines"] == 52
    assert e["C_isotropic"] == 9192
    assert (e["over"], e["under"], e["ghost"]) == (1.0, 98.8, 0.2)
    assert e["detectable"] == 96.8
    assert e["median_m"] == 1.29 and e["median_N0"] == 240
    assert abs(e["sigma2"] - 0.938) < 5e-4


def test_genetic_spectrum_golden():
    g = _load("genetic_two_arm.json")
    golden = {  # over / under among detectable, detectable %
        "orion_HCT116": (0.0, 100.0, 14.7),
        "orion_HEK293T": (0.1, 99.9, 35.3),
        "trade_jurkat": (1.2, 98.8, 70.6),
        "trade_hepg2": (4.4, 95.6, 60.9),
    }
    for k, (ov, un, det) in golden.items():
        v = g[k]
        assert (v["det_pct_over"], v["det_pct_under"]) == (ov, un), k
        assert v["pct_detectable"] == det, k


# ---- self-consistency of each row ----
def test_percent_blocks_sum_to_100():
    e = _load("emeraldbay_spectrum.json")
    assert abs(e["over"] + e["under"] + e["ghost"] - 100) < 0.15
    assert abs(e["detectable"] + e["not_detectable"] - 100) < 0.15
    c = _load("tahoe_constants.json")
    assert abs(c["pct_OVER"] + c["pct_UNDER"] + c["pct_Ghost"] - 100) < 0.15
    g = _load("genetic_two_arm.json")
    for k in ("orion_HCT116", "orion_HEK293T", "trade_jurkat", "trade_hepg2"):
        v = g[k]
        assert abs(v["det_pct_over"] + v["det_pct_under"] - 100) < 0.15, k
        assert abs(v["pct_detectable"] + v["pct_not_detectable"] - 100) < 0.15, k


def test_emeraldbay_aggregate_reproduces_from_per_group():
    # the isotropic rule on the committed per-group m/N must reproduce the aggregate over/ghost/detectable
    e = _load("emeraldbay_spectrum.json")
    m = np.array(e["per_group"]["m"]); N = np.array(e["per_group"]["N"])
    C = 2 * (d - 1) * e["sigma2"] / 0.1 ** 2
    nstar = C / np.maximum(m ** 2, 1e-12)
    over = 100 * np.mean(nstar < N)
    ghost = 100 * np.mean(nstar > 50000)
    det = 100 * np.mean(m ** 2 > 1.25 * d * e["sigma2"] * (1.0 / N))
    assert round(over, 1) == e["over"], (over, e["over"])
    assert round(ghost, 1) == e["ghost"], (ghost, e["ghost"])
    assert round(det, 1) == e["detectable"], (det, e["detectable"])
    assert round(C) == e["C_isotropic"]


def test_median_nstar_matches_C_over_median_m2():
    # median n* implied by the isotropic constant and the reported median m (order-of-magnitude lock)
    e = _load("emeraldbay_spectrum.json")
    implied = e["C_isotropic"] / e["median_m"] ** 2
    assert 5000 < implied < 6500, implied  # 9192 / 1.29^2 = 5524


# ---- manuscript text must match the fixtures (numbers are the source of truth) ----
def test_manuscript_chemical_rows_match_fixtures():
    txt = open(DOC).read()
    e = _load("emeraldbay_spectrum.json")
    # Table 1 + Table 5 EmeraldBay regime block
    assert "1.0 / 98.8 / 0.2" in txt, "manuscript EmeraldBay regime block drifted from fixture"
    # unit declaration and non-empty group count
    assert "4,912" in txt and "(sample × line)" in txt
    # Tahoe headline must be present and unchanged
    assert "10.6% are over-sampled" in txt or "10.6\\% are over-sampled" in txt
    # no stale dose-pooled EmeraldBay spectrum numbers
    for stale in ("4.6 / 84.0 / 11.4", "84.0", "11.4 / "):
        assert stale not in txt, f"stale EmeraldBay number {stale!r} still in manuscript"
