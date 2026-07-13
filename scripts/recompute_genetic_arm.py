"""
Regenerate the genetic (Orion + TRADE) quota fixtures under the EXPLICIT two-arm model
(src/quota_arm.py), reading the COMMITTED per-gene CSVs — no atlas re-stream.

Per gene we already have (from the committed CSV):  m, n_cells (=acquired treated depth N0),
n_star_aniso (= 2 tr(PSP)/(m^2 theta^2) = the old factor-2 quota).  Hence
    base = n_star_aniso / 2 = tr(PSP)/(m^2 theta^2),   trPSP = base * m^2 * theta^2.
Apply the exact two-arm quota with n_c = the acquired NTC pool, classify into four regimes,
and emit an updated summary fixture (sigma^2 unchanged — canonical, from streaming).
"""
import os, json, argparse, numpy as np, pandas as pd, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from quota_arm import quota_two_arm, m_min_resolvable, classify

TH = 0.1; GHOST = 50000.0
DATASETS = [  # (key, csv, summary_json, n_ntc, sigma2_within)
    ("orion_HCT116",  "orion_HCT116_quota.csv",  "orion_HCT116_summary.json",  165562, 0.9133),
    ("orion_HEK293T", "orion_HEK293T_quota.csv", "orion_HEK293T_summary.json", 218838, 1.0334),
    ("trade_jurkat",  "trade_jurkat_quota.csv",  "trade_jurkat_summary.json",  11514,  1.4975),
    ("trade_hepg2",   "trade_hepg2_quota.csv",   "trade_hepg2_summary.json",   4380,   1.8667),
]
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")


def recompute(csv, n_ntc):
    df = pd.read_csv(os.path.join(FX, csv))
    m = df["m"].values.astype(float); N0 = df["n_cells"].values.astype(float)
    base = df["n_star_aniso"].values.astype(float) / 2.0        # = tr(PSP)/(m^2 th^2)
    trPSP = base * m ** 2 * TH ** 2
    n_t = quota_two_arm(trPSP, m, TH, n_ntc)                    # exact two-arm, may be inf
    limited = ~np.isfinite(n_t)
    over = (np.isfinite(n_t)) & (N0 >= n_t)
    ghost = (np.isfinite(n_t)) & (n_t > GHOST)
    under = (np.isfinite(n_t)) & (~over) & (~ghost)
    res = np.isfinite(n_t)
    n = len(df)
    med_base = float(np.median(base))                          # NTC needed to un-limit the median gene
    out = dict(
        n_perturbations=n, n_ntc_cells=int(n_ntc), theta=TH,
        C_largepool_limit=round(base_const(df, n_ntc), 1),     # 49*sigma2/theta^2 style, but from data
        pct_over=round(100 * over.mean(), 2),
        pct_under=round(100 * under.mean(), 2),
        pct_ghost=round(100 * ghost.mean(), 2),
        pct_pool_limited=round(100 * limited.mean(), 2),
        n_star_median_resolvable=int(np.median(n_t[res])) if res.any() else None,
        median_deficit_resolvable=round(float(np.median((n_t[res] / N0[res]))), 1) if res.any() else None,
        median_n_star_all=("inf" if limited.mean() >= 0.5 else int(np.median(n_t[np.isfinite(n_t)]))),
        m_min_median_trPSP=round(float(m_min_resolvable(np.median(trPSP), n_ntc, TH)), 3),
        ntc_to_unlimit_median_gene=int(round(med_base)),        # need n_c >= median(base)
        ntc_shortfall_for_median=int(round(max(0, med_base - n_ntc))),
    )
    return out


def base_const(df, n_ntc):
    # dataset "quota law constant" printed = large-pool (factor-1) limit = tr(PSP)/(m^2 th^2) scaled;
    # report the isotropic constant (d-1)*sigma2/theta^2 equivalent = median over genes of base*m^2
    return float(np.median(df["n_star_aniso"].values) / 2.0 * np.median(df["m"].values) ** 2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--write", action="store_true"); a = ap.parse_args()
    allout = {}
    for key, csv, js, n_ntc, s2 in DATASETS:
        r = recompute(csv, n_ntc); r["sigma2_within"] = s2; allout[key] = r
        print(f"\n=== {key}  (n_ntc={n_ntc:,}, sigma^2={s2}) ===")
        print(f"  over {r['pct_over']}%  under {r['pct_under']}%  ghost>50k {r['pct_ghost']}%  "
              f"POOL-LIMITED {r['pct_pool_limited']}%")
        print(f"  median n* (resolvable) {r['n_star_median_resolvable']}  deficit {r['median_deficit_resolvable']}x  "
              f"median-all {r['median_n_star_all']}")
        print(f"  m_min {r['m_min_median_trPSP']}   NTC to un-limit median gene {r['ntc_to_unlimit_median_gene']:,} "
              f"(shortfall {r['ntc_shortfall_for_median']:,})")
        if a.write:
            js_path = os.path.join(FX, js); s = json.load(open(js_path))
            s["two_arm"] = r
            json.dump(s, open(js_path, "w"), indent=1)
    json.dump(allout, open(os.path.join(FX, "genetic_two_arm.json"), "w"), indent=1)
    print("\nwrote fixtures/genetic_two_arm.json" + (" + updated per-dataset summaries" if a.write else " (dry run)"))
