"""
Regenerate the genetic (Orion + TRADE) quota fixtures under the exact two-arm model
(src/quota_arm.py) AND the detection floor, from the committed per-gene CSVs (no re-stream).

Detection floor (from scripts/orion_recompute/pass3_quota_full.py): a knockdown is DETECTABLE when
snr_floor = m_raw/sqrt(trS) > 1.5, i.e. its magnitude clears the per-perturbation sampling noise.
Genes below the floor have m ~ noise and are NOT classified further (a detailed regime for them is
meaningless). The HEADLINE spectrum is over/under among DETECTABLE knockdowns only; ghost and
control-pool-limited are reported as secondary diagnostics, because among detectable knockdowns they
sit at 0-2% (the detection floor sits above the pool-limited floor m_min, so the two barely overlap).
"""
import os, json, argparse, numpy as np, pandas as pd, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from quota_arm import quota_two_arm, m_min_resolvable

TH = 0.1; GHOST = 50000.0
DATASETS = [  # (key, csv, summary_json, n_ntc, sigma2_within)
    ("orion_HCT116",  "orion_HCT116_quota.csv",  "orion_HCT116_summary.json",  165562, 0.9133),
    ("orion_HEK293T", "orion_HEK293T_quota.csv", "orion_HEK293T_summary.json", 218838, 1.0334),
    ("trade_jurkat",  "trade_jurkat_quota.csv",  "trade_jurkat_summary.json",  11514,  1.4975),
    ("trade_hepg2",   "trade_hepg2_quota.csv",   "trade_hepg2_summary.json",   4380,   1.8667),
]
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")


def recompute(csv, n_ntc, s2):
    df = pd.read_csv(os.path.join(FX, csv))
    m = df["m"].values.astype(float); N0 = df["n_cells"].values.astype(float)
    base = df["n_star_aniso"].values.astype(float) / 2.0        # = tr(PSP)/(m^2 th^2)
    trPSP = base * m ** 2 * TH ** 2
    n_t = quota_two_arm(trPSP, m, TH, n_ntc)
    det = df["snr_floor"].values > 1.5                          # DETECTION FLOOR (code definition)
    d = det                                                     # shorthand

    # headline: over/under among DETECTABLE only
    over_d = float(np.mean(np.isfinite(n_t[d]) & (N0[d] >= n_t[d])) * 100)
    under_d = round(100 - over_d, 1)
    res_d = det & np.isfinite(n_t)
    out = dict(
        n_perturbations=len(df), n_ntc_cells=int(n_ntc), sigma2_within=s2, theta=TH,
        C_largepool_limit=round(49 * s2 / TH ** 2, 0),
        # (1) detectable fraction — a dataset property, independent of the quota
        pct_detectable=round(100 * float(d.mean()), 1),
        pct_not_detectable=round(100 * float((~d).mean()), 1),
        # (2) headline spectrum among detectable (over/under only; ghost retired for genetic)
        det_pct_over=round(over_d, 1),
        det_pct_under=under_d,
        det_median_nstar=int(np.median(n_t[res_d])) if res_d.any() else None,
        det_median_deficit=round(float(np.median(n_t[res_d] / N0[res_d])), 1) if res_d.any() else None,
        # secondary: control-pool-limited (demoted) — negligible among detectable, large among all
        pct_pool_limited_among_detectable=round(100 * float(np.mean(base[d] >= n_ntc)), 1),
        pct_pool_limited_among_all=round(100 * float(np.mean(base >= n_ntc)), 1),
        m_min=round(float(m_min_resolvable(np.median(trPSP), n_ntc, TH)), 2),
        min_detectable_m=round(float(m[d].min()), 2) if d.any() else None,
    )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--write", action="store_true"); a = ap.parse_args()
    allout = {}
    for key, csv, js, n_ntc, s2 in DATASETS:
        r = recompute(csv, n_ntc, s2); allout[key] = r
        print(f"\n=== {key} (n_ntc={n_ntc:,}) ===")
        print(f"  detectable {r['pct_detectable']}%  (not-detectable {r['pct_not_detectable']}%)")
        print(f"  among detectable: OVER {r['det_pct_over']}% / UNDER {r['det_pct_under']}%  "
              f"(median n* {r['det_median_nstar']}, deficit {r['det_median_deficit']}x)")
        print(f"  [demoted] pool-limited among detectable {r['pct_pool_limited_among_detectable']}%  "
              f"(among all {r['pct_pool_limited_among_all']}%);  m_min {r['m_min']}, min detectable m {r['min_detectable_m']}")
        if a.write:
            p = os.path.join(FX, js); s = json.load(open(p)); s["two_arm"] = r
            json.dump(s, open(p, "w"), indent=1)
    json.dump(allout, open(os.path.join(FX, "genetic_two_arm.json"), "w"), indent=1)
    print("\nwrote fixtures/genetic_two_arm.json" + (" + summaries" if a.write else " (dry run)"))
