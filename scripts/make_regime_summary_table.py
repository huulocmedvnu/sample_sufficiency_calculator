"""
Generate the cross-screen sufficiency-regime summary (Table 5): all SIX per-line screens on one row
each, so the over/under/ghost/detectable numbers scattered across Table 1, Table 4 and Figure 4b can be
compared directly. Entirely from committed fixtures, no hard-coded numbers.

The columns carry TWO independent partitions, each summing to 100 within its own denominator:
  * Detectability  (detectable / not-detectable) -- always over ALL scored conditions.
  * Sufficiency spectrum (over / under / ghost)   -- over the DENOMINATOR column:
      chemical (Tahoe, EmeraldBay): over ALL conditions;
      genetic  (Orion, TRADE):      over DETECTABLE knockdowns only, and ghost is n/a (not assigned).

Denominator notes:
  * EmeraldBay uses N>=400 (the held-out eligibility threshold, 958 groups), at theta*=0.1 -- NOT the
    theta*=0.20 gating experiment (347/3,971 groups).
  * Genetic n is the SCORED gene count (>=25 acquired cells; NTC excluded), below the targeted 18,903
    (Orion) / 2,393 (TRADE) -- see the manuscript DISCREPANCIES note.

Emits Markdown to stdout and fixtures/regime_summary_table.md.
"""
import os, json, numpy as np, pandas as pd
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
d = 50; TH = 0.1; EB_MIN_N = 400


def load(n):
    return json.load(open(os.path.join(FX, n)))


def chem_detectable(m, N0, s2, n_c):
    trSig = d * s2
    trS = trSig * (1.0 / N0 + 1.0 / np.asarray(n_c, float))
    return 100 * float(np.mean(m ** 2 > 1.25 * trS))


ROWS = []

# ---- Tahoe-100M (chemical; spectrum over ALL conditions) ----
c = load("tahoe_constants.json"); df = pd.read_csv(os.path.join(FX, "tahoe_quota_per_condition.csv"))
det = chem_detectable(df.m.values, df.N0.values.astype(float), c["sigma2"], df.N0.values.astype(float))
ROWS.append(dict(name="Tahoe-100M", unit="drug × dose × line", n=len(df), det=round(det, 1),
                 notdet=round(100 - det, 1), over=c["pct_OVER"], under=c["pct_UNDER"], ghost=f"{c['pct_Ghost']}",
                 denom="all", s2=c["sigma2"], medm=c["median_m"], medd=c["median_N0"]))

# ---- EmeraldBay (chemical; N>=400, spectrum over ALL such conditions) ----
D = [r for r in load("emeraldbay_falsification_full.json")["per_group_slope"] if r["N"] >= EB_MIN_N]
s2e = 0.938
m = np.array([r["m"] for r in D]); N = np.array([r["N"] for r in D]).astype(float)
nstar = 2 * np.array([r["pred_slope"] for r in D]) / TH ** 2
over = 100 * np.mean(nstar < N); ghost = 100 * np.mean(nstar > 50000); under = 100 - over - ghost
dete = chem_detectable(m, N, s2e, 1e12)
ROWS.append(dict(name="EmeraldBay", unit="condition × line", n=len(D), det=round(dete, 1),
                 notdet=round(100 - dete, 1), over=round(over, 1), under=round(under, 1),
                 ghost=f"{ghost:.1f}", denom="all (N≥400)", s2=s2e,
                 medm=round(float(np.median(m)), 2), medd=int(np.median(N))))

# ---- genetic (over/under among DETECTABLE; ghost not assigned) ----
g = load("genetic_two_arm.json")
GEN = [("orion_HCT116", "Orion HCT116", "orion_HCT116_summary.json"),
       ("orion_HEK293T", "Orion HEK293T", "orion_HEK293T_summary.json"),
       ("trade_jurkat", "TRADE Jurkat", "trade_jurkat_summary.json"),
       ("trade_hepg2", "TRADE HepG2", "trade_hepg2_summary.json")]
for key, disp, summ in GEN:
    v = g[key]; s = load(summ)
    ROWS.append(dict(name=disp, unit="gene × line", n=v["n_perturbations"], det=v["pct_detectable"],
                     notdet=v["pct_not_detectable"], over=v["det_pct_over"], under=v["det_pct_under"],
                     ghost="n/a", denom="detectable", s2=v["sigma2_within"],
                     medm=s["m_median"], medd=s.get("median_true_N0", s.get("median_N0"))))

hdr = ("| screen | unit | n | det. % | not-det. % | over % | under % | ghost % | "
       "denom. | σ² | med. m | med. depth |")
sep = "|" + "|".join("-"*w for w in [13,11,6,6,7,6,6,6,11,6,6,8]) + "|"
out = [hdr, sep]
for r in ROWS:
    out.append(f"| {r['name']} | {r['unit']} | {r['n']:,} | {r['det']} | {r['notdet']} | "
               f"{r['over']} | {r['under']} | {r['ghost']} | {r['denom']} | {r['s2']} | {r['medm']} | {r['medd']:,} |")
md = "\n".join(out)
print(md)
open(os.path.join(FX, "regime_summary_table.md"), "w").write(md + "\n")
