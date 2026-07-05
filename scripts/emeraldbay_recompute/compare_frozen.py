#!/usr/bin/env python
"""Side-by-side comparison of the committed (per-fit HVG) vs frozen-HVG EmeraldBay recompute.

Prints every manuscript-facing EmeraldBay number old vs new so the manuscript/fixtures can be updated
from a single audited table. Reads committed fixtures from a snapshot dir (default: the live committed
fixtures) and the frozen ones written by run_frozen.sh.
"""
import os, json, argparse

HERE = os.path.dirname(__file__)
FIXDIR = os.path.join(HERE, "..", "..", "fixtures")


def load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def row(label, old, new, pct=False):
    def fmt(x):
        if x is None: return "  --  "
        if isinstance(x, float): return f"{x:.4f}"
        return str(x)
    d = ""
    if isinstance(old, (int, float)) and isinstance(new, (int, float)) and old not in (0, None):
        rel = 100 * (new - old) / abs(old)
        d = f"  (Δ {rel:+.1f}%)"
    print(f"  {label:<34} {fmt(old):>12} -> {fmt(new):>12}{d}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--committed-dir", default=FIXDIR,
                    help="dir with committed emeraldbay_*.json (default live fixtures)")
    a = ap.parse_args()
    cd = a.committed_dir

    cal_o = load(os.path.join(cd, "emeraldbay_calibration.json"))
    cal_n = load(os.path.join(FIXDIR, "emeraldbay_calibration_frozen.json"))
    fal_o = load(os.path.join(cd, "emeraldbay_falsification_full.json"))
    fal_n = load(os.path.join(FIXDIR, "emeraldbay_falsification_frozen.json"))

    print("=" * 66)
    print("EMERALDBAY  committed (per-fit HVG)  ->  frozen HVG")
    print("=" * 66)

    if cal_o and cal_n:
        print("\n[calibration] within-condition sigma^2 + held-out curves + gating")
        row("sigma2_mean (within, pass3)", cal_o.get("sigma2_mean"), cal_n.get("sigma2_mean"))
        for g in cal_o.get("heldout_curves", {}):
            co = cal_o["heldout_curves"][g]; cn = cal_n.get("heldout_curves", {}).get(g, {})
            print(f"   -- {g}")
            row("   slope", co.get("slope"), cn.get("slope"))
            row("   expected_slope", co.get("expected_slope"), cn.get("expected_slope"))
            row("   r2", co.get("r2"), cn.get("r2"))
            row("   mean_rel_err", co.get("mean_rel_err"), cn.get("mean_rel_err"))
            row("   m", co.get("m"), cn.get("m"))
        row("gating_n_groups", cal_o.get("gating_n_groups"), cal_n.get("gating_n_groups"))
        row("gating_n_over", cal_o.get("gating_n_over"), cal_n.get("gating_n_over"))
        row("gating_pct_over", cal_o.get("gating_pct_over"), cal_n.get("gating_pct_over"))
        row("gating_over_accuracy", cal_o.get("gating_over_accuracy"), cal_n.get("gating_over_accuracy"))
    else:
        print("\n[calibration] frozen fixture not ready" if not cal_n else "[calibration] committed missing")

    if fal_o and fal_n:
        print("\n[falsification full-atlas] population-scale slope test")
        so, sn = fal_o.get("sigma2_confirmation", {}), fal_n.get("sigma2_confirmation", {})
        row("sigma2 within (52 lines)", so.get("within_condition"), sn.get("within_condition"))
        row("sigma2 marginal", so.get("marginal"), sn.get("marginal"))
        row("n_groups_slope", fal_o.get("n_groups_slope"), fal_n.get("n_groups_slope"))
        eo, en = fal_o.get("exp1_predicted_vs_realized", {}), fal_n.get("exp1_predicted_vs_realized", {})
        row("in-regime n (rho^2>=3)", eo.get("n_in_regime"), en.get("n_in_regime"))
        row("median_ratio_in_regime", eo.get("median_ratio_in_regime"), en.get("median_ratio_in_regime"))
        row("r2_median_in_regime", eo.get("r2_median_in_regime"), en.get("r2_median_in_regime"))
        ho, hn = fal_o.get("exp3_two_halves", {}), fal_n.get("exp3_two_halves", {})
        row("exp3 median_ratio (sqrt2=1.41)", ho.get("median_ratio"), hn.get("median_ratio"))
    else:
        print("\n[falsification] frozen fixture not ready" if not fal_n else "[falsification] committed missing")
    print("=" * 66)


if __name__ == "__main__":
    main()
