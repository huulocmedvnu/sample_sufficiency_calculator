"""
Golden lock on the within-condition sigma^2 of every dataset, under the ONE canonical
definition (src/sigma2_canonical.py): per-cell variance pooled over (cell line x condition)
groups, mean over PCs. This is the definition that reproduces the headline Tahoe-100M 0.9567
and EmeraldBay-full 0.9380.

Purpose: stop stale/mis-scoped sigma^2 values (the 0.963 / 0.848 / 0.9174 incident) from ever
leaking again. Raw per-cell atlases are not in CI, so this pins the committed fixture values
against their golden constants; the recompute that produced them is in scripts/*_recompute/ and
verified once in scripts (EmeraldBay reproduced 0.9380 / 0.8673 from eb_work per-cell data).
"""
import json, os
import numpy as np
import pytest
from sigma2_canonical import within_condition_sigma2, within_condition_sigma2_moments

FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")

# canonical within-condition sigma^2 (scope) -> golden value
GOLDEN = {
    "Tahoe-100M (full atlas)":       0.9567,
    "EmeraldBay (full 52-line)":     0.9380,
    "EmeraldBay (5-shared-line)":    0.8673,
    "Orion HCT116":                  0.9133,
    "Orion HEK293T":                 1.0334,
    "TRADE Jurkat":                  1.4975,
    "TRADE HepG2":                   1.8667,
}


def _load(name):
    return json.load(open(os.path.join(FX, name)))


def test_sigma2_fixtures_match_golden():
    got = {
        "Tahoe-100M (full atlas)":    _load("tahoe_constants.json")["sigma2"],
        "EmeraldBay (full 52-line)":  _load("emeraldbay_falsification_full.json")["sigma2_confirmation"]["within_condition"],
        "EmeraldBay (5-shared-line)": _load("emeraldbay_calibration.json")["sigma2_mean"],
        "Orion HCT116":               _load("orion_HCT116_summary.json")["sigma2_within"],
        "Orion HEK293T":              _load("orion_HEK293T_summary.json")["sigma2_within"],
        "TRADE Jurkat":               _load("trade_jurkat_summary.json")["sigma2_within"],
        "TRADE HepG2":                _load("trade_hepg2_summary.json")["sigma2_within"],
    }
    for k, gold in GOLDEN.items():
        assert abs(got[k] - gold) < 5e-4, f"{k}: fixture {got[k]} != golden {gold}"


def test_5line_matches_falsification_fixture():
    # the 5-line calibration value must equal the 5-line falsification value (same canonical def)
    a = _load("emeraldbay_calibration.json")["sigma2_mean"]
    b = _load("emeraldbay_falsification.json")["sigma2_confirmation"]["within_condition"]
    assert abs(a - b) < 5e-4, f"5-line calibration {a} != falsification {b}"


def test_canonical_function_pooled_identity():
    # per-cell form == moment form under identical grouping (guards against the 0.9174 regrouping bug)
    rng = np.random.default_rng(0)
    coords = rng.normal(size=(500, 8))
    keys = rng.integers(0, 12, size=500).astype(str)
    a, _ = within_condition_sigma2(coords, keys, return_per_pc=True)
    uk = np.unique(keys)
    sums = np.array([coords[keys == g].sum(0) for g in uk])
    sumsq = np.array([(coords[keys == g] ** 2).sum(0) for g in uk])
    cnt = np.array([(keys == g).sum() for g in uk])
    b, _ = within_condition_sigma2_moments(sums, sumsq, cnt, return_per_pc=True)
    assert abs(within_condition_sigma2(coords, keys) - within_condition_sigma2_moments(sums, sumsq, cnt)) < 1e-9
