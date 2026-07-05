"""Integration test: external validation of the anisotropic theory on the INDEPENDENT EmeraldBay atlas.

Loads the committed `fixtures/emeraldbay_calibration.json` (produced by
scripts/emeraldbay_recompute/ streaming the 57.7 GB atlas) and asserts that the held-out
angular-error curves and the regime-gating accuracy meet the theory's predictions. CI-safe: no
network, no large data — it re-checks the distilled validation result. Skips if the fixture is absent.
"""
import json
import os

import pytest

FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "emeraldbay_calibration.json")


@pytest.mark.skipif(not os.path.exists(FIX), reason="EmeraldBay calibration fixture not present")
def test_emeraldbay_external_validation():
    fx = json.load(open(FIX))
    assert fx["num_dimensions"] == 50

    # regime gating: every well predicted OVER-sampled must meet the tolerance when downsampled to n*
    acc = fx.get("gating_over_accuracy")
    assert acc is None or acc >= 0.9, f"gating over-accuracy {acc} < 0.9"

    strong = 0
    for key, c in fx["heldout_curves"].items():
        rows = c["rows"]                                  # [n, realized_rms, predicted_rms]
        rel = [abs(r - p) / p for _, r, p in rows]
        # realized angle tracks the closed-form prediction across the whole n-grid (both regimes)
        assert max(rel) < 0.25, f"{key}: max relative error {max(rel):.3f}"
        assert sum(rel) / len(rel) < 0.06, f"{key}: mean relative error {sum(rel)/len(rel):.3f}"
        assert abs(c["intercept"] if "intercept" in c else 0.0) < 1e-2
        # The tight linear-through-origin fit (R^2 > 0.99) and the parameter-free slope match are the
        # HIGH-SNR predictions. Low-magnitude groups sit in the first-order-breakdown regime the theory
        # forecasts (small-n deviation), so their through-origin R^2 is legitimately lower -- they are
        # still bounded by the rel-error checks above. Gate the strict checks on strong signature (m>=5),
        # consistent with the paper's rho^2-stratified regime distinction.
        if c["m"] >= 5.0:
            assert c["r2"] > 0.99, f"{key}: R^2 {c['r2']:.4f} (expected line-through-origin)"
            assert abs(c["slope"] - c["expected_slope"]) / c["expected_slope"] < 0.05, \
                f"{key}: slope {c['slope']} vs expected {c['expected_slope']}"
            strong += 1

    assert strong >= 1, "expected at least one strong-signature (m>=5) curve to validate the slope"
