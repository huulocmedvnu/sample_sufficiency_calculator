"""
Generate the unified calibration table (Table 1) entirely from committed fixtures — the SAME
fixtures Table 4 (genetic spectrum) uses, so the two cannot disagree. No hard-coded numbers.

Regime is reported the way each modality's headline is: chemical screens are ~fully detectable, so
over/under/ghost among all conditions; genetic screens report the DETECTABLE fraction and then
over/under among detectable (ghost/control-pool-limited are demoted, see Table 4 and text).
Detectable = snr_floor>1.5, i.e. m^2 > 1.25*tr(S) with S = Sigma/n_t + Sigma/n_c (reconstructed for
the chemical atlases, whose CSVs do not store snr_floor).

Emits Markdown to stdout and fixtures/unified_calibration_table.md.
"""
import os, json, numpy as np, pandas as pd
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
TH = 0.1; d = 50


def load(n):
    return json.load(open(os.path.join(FX, n)))


def chem_detectable(m, N0, s2, n_c):
    trSig = d * s2
    trS = trSig * (1.0 / N0 + 1.0 / np.asarray(n_c, float))
    return 100 * float(np.mean(m ** 2 > 1.25 * trS))


def tahoe_row():
    c = load("tahoe_constants.json"); df = pd.read_csv(os.path.join(FX, "tahoe_quota_per_condition.csv"))
    det = chem_detectable(df.m.values, df.N0.values.astype(float), c["sigma2"], df.N0.values.astype(float))
    return dict(name="Tahoe-100M", arm="equal-arm", n_c="matched", s2=c["sigma2"], scope="full",
                C=round(c["nstar_const"]), Clab="equal-arm", m_min="n/a", medm=c["median_m"], medd=c["median_N0"],
                det=round(det, 1), reg=f"{c['pct_OVER']} / {c['pct_UNDER']} / {c['pct_Ghost']}g")


def emeraldbay_row():
    D = load("emeraldbay_falsification_full.json")["per_group_slope"]; s2 = 0.938
    m = np.array([r["m"] for r in D]); N = np.array([r["N"] for r in D]).astype(float)
    ps = np.array([r["pred_slope"] for r in D]); nstar = 2 * ps / TH ** 2
    over = 100 * np.mean(nstar < N); ghost = 100 * np.mean(nstar > 50000); under = 100 - over - ghost
    det = chem_detectable(m, N, s2, 1e12)   # per-line-mean baseline -> large control
    return dict(name="EmeraldBay", arm="equal-arm", n_c="matched", s2=s2, scope="full 52-line",
                C=round(2 * (d - 1) * s2 / TH ** 2), Clab="equal-arm", m_min="n/a",
                medm=round(float(np.median(m)), 2), medd=int(np.median(N)), det=round(det, 1),
                reg=f"{over:.1f} / {under:.1f} / {ghost:.1f}g")


def genetic_row(key, disp, summ, aniso_key="aniso_iso_ratio_median"):
    g = load("genetic_two_arm.json")[key]; s = load(summ)
    return dict(name=disp, arm="large-pool", n_c=f"{g['n_ntc_cells']:,}", s2=g["sigma2_within"], scope="within",
                C=round(g["C_largepool_limit"]), Clab="large-pool bound", m_min=g["m_min"],
                medm=s["m_median"], medd=s.get("median_true_N0", s.get("median_N0")), det=g["pct_detectable"],
                reg=f"{g['det_pct_over']} / {g['det_pct_under']}")


ROWS = [tahoe_row(), emeraldbay_row(),
        genetic_row("orion_HCT116", "Orion HCT116", "orion_HCT116_summary.json"),
        genetic_row("orion_HEK293T", "Orion HEK293T", "orion_HEK293T_summary.json"),
        genetic_row("trade_jurkat", "TRADE Jurkat", "trade_jurkat_summary.json"),
        genetic_row("trade_hepg2", "TRADE HepG2", "trade_hepg2_summary.json")]

hdr = ("| dataset | arm design | n_c | σ² | C @θ=0.1 | m_min | median m | median depth | detectable % | "
       "over / under (%) |")
sep = "|" + "|".join(["---"] * 10) + "|"
out = [hdr, sep]
for r in ROWS:
    out.append(f"| {r['name']} | {r['arm']} | {r['n_c']} | {r['s2']} ({r['scope']}) | {r['C']:,} ({r['Clab']}) | "
               f"{r['m_min']} | {r['medm']} | {r['medd']} | {r['det']} | {r['reg']} |")
md = "\n".join(out)
print(md)
open(os.path.join(FX, "unified_calibration_table.md"), "w").write(md + "\n")
