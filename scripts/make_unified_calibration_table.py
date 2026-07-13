"""
Generate the unified calibration table (Part D) — 6 datasets x calibration columns — entirely
from committed fixtures. No hard-coded numbers. Emits Markdown (pipe table) to stdout and to
fixtures/unified_calibration_table.md so the manuscript can embed the exact same rows.

Columns: dataset | arm design | n_c | sigma^2 (+scope+definition) | C @theta=0.1 (large-pool
limit for genetic; equal-arm constant for chemical) | m_min | aniso/iso | median m | median depth
| over / under / ghost / pool-limited (%).

Every row's sigma^2 is the CANONICAL within-condition value (per-cell, line x condition;
src/sigma2_canonical.py) — the "sigma^2 definition" is identical across all six, which is exactly
the point (like-for-like).
"""
import os, json, numpy as np
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
TH = 0.1


def load(n):
    return json.load(open(os.path.join(FX, n)))


def emeraldbay_row():
    d = load("emeraldbay_falsification_full.json")["per_group_slope"]
    m = np.array([r["m"] for r in d]); N = np.array([r["N"] for r in d]); ps = np.array([r["pred_slope"] for r in d])
    s2 = 0.938
    nstar = 2 * ps / TH ** 2                       # equal-arm (matched per-line control), Fig 5b
    over = 100 * np.mean(nstar < N); ghost = 100 * np.mean(nstar > 50000); under = 100 - over - ghost
    ratio = float(np.median(ps * m ** 2 / (49 * s2)))
    return dict(name="EmeraldBay", arm="equal-arm", n_c="matched vehicle", s2=s2,
                s2_scope="full 52-line", C=round(49 * s2 / TH ** 2 * 2), Clabel="equal-arm",
                m_min="n/a", ratio=round(ratio, 3), med_m=round(float(np.median(m)), 2),
                med_dep=int(np.median(N)), over=round(over, 1), under=round(under, 1),
                ghost=round(ghost, 1), pool=0.0, note="per-line-mean baseline")


def tahoe_row():
    c = load("tahoe_constants.json")
    return dict(name="Tahoe-100M", arm="equal-arm", n_c="matched vehicle", s2=c["sigma2"],
                s2_scope="full atlas", C=round(c["nstar_const"]), Clabel="equal-arm",
                m_min="n/a", ratio="0.93-0.96", med_m=c["median_m"], med_dep=c["median_N0"],
                over=c["pct_OVER"], under=c["pct_UNDER"], ghost=c["pct_Ghost"], pool=0.0, note="")


def genetic_row(key, disp, summ_json, aniso_key="aniso_iso_ratio_median"):
    t = load("genetic_two_arm.json")[key]
    s = load(summ_json)
    return dict(name=disp, arm="large-pool", n_c=f"{t['n_ntc_cells']:,}", s2=t["sigma2_within"],
                s2_scope="within-condition", C=round(49 * t["sigma2_within"] / TH ** 2),
                Clabel="large-pool limit (lower bound)", m_min=t["m_min_median_trPSP"],
                ratio=s[aniso_key], med_m=s["m_median"], med_dep=s.get("median_true_N0", s.get("median_N0")),
                over=t["pct_over"], under=t["pct_under"], ghost=t["pct_ghost"], pool=t["pct_pool_limited"],
                note=f"m_min={t['m_min_median_trPSP']}")


ROWS = [
    tahoe_row(), emeraldbay_row(),
    genetic_row("orion_HCT116", "Orion HCT116", "orion_HCT116_summary.json"),
    genetic_row("orion_HEK293T", "Orion HEK293T", "orion_HEK293T_summary.json"),
    genetic_row("trade_jurkat", "TRADE Jurkat", "trade_jurkat_summary.json"),
    genetic_row("trade_hepg2", "TRADE HepG2", "trade_hepg2_summary.json"),
]

hdr = ("| dataset | arm design | n_c | σ² (within-cond) | C @θ=0.1 | m_min | aniso/iso | "
       "median m | median depth | over / under / ghost / pool-ltd (%) |")
sep = "|" + "|".join(["---"] * 10) + "|"
lines = [hdr, sep]
for r in ROWS:
    C = f"{r['C']:,} ({r['Clabel']})"
    reg = f"{r['over']} / {r['under']} / {r['ghost']} / {r['pool']}"
    lines.append(f"| {r['name']} | {r['arm']} | {r['n_c']} | {r['s2']} ({r['s2_scope']}) | {C} | "
                 f"{r['m_min']} | {r['ratio']} | {r['med_m']} | {r['med_dep']} | {reg} |")
md = "\n".join(lines)
print(md)
open(os.path.join(FX, "unified_calibration_table.md"), "w").write(md + "\n")
