"""
Unified sufficiency-regime spectrum for all six per-line screens, computed under ONE consistent unit and
denominator rule (fixes the EmeraldBay condition-unit error documented in docs/AUDIT_DENOMINATORS.md #1).

UNIT: condition = (perturbation x dose x cell line).
  - Tahoe:      (drug x dose x line), 56,827      -- unchanged, read from the committed CSV.
  - EmeraldBay: (sample x line), 4,912            -- sample = (drug/combo, dose) per sample2cond.json;
                computed here from eb_work/out/pseudobulk.npz (per-condition sufficient statistics).
                NOT per_group_slope, which pools over dose (drug-name x line) and is for the held-out test.
  - Genetic:    (gene x line)                     -- unchanged, from genetic_two_arm.json.

DENOMINATOR: chemical spectra are over ALL conditions (no cell-count filter, N0 >= 1, matching Tahoe).
  Genetic keep the N>=25 scoring QC and report over/under among DETECTABLE only (ghost not assigned).

METHOD (matches Tahoe's committed spectrum): ISOTROPIC quota n* = C/m^2 with C = 2(d-1)sigma^2/theta^2
  (EmeraldBay C = 9,192, its own Table 1 constant; Tahoe C = 9,376). Tahoe's n* reproduces 9,376/m^2 to
  0 relative error, i.e. the committed chemical spectrum is isotropic (see AUDIT #6). EmeraldBay magnitude
  is the raw ||centroid - per-line-mean|| (the prior EmeraldBay convention; no new bias-correction step).

Generates:
  fixtures/emeraldbay_spectrum.json   -- EmeraldBay per-(sample x line) spectrum (aggregate + per-group)
  fixtures/spectrum_summary.md        -- Table 5, all six screens, two-block layout
"""
import os, json, ast, numpy as np, pandas as pd
FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
EB = "/mnt/hdd2/loc-tran/eb_work/out"
d = 50; TH = 0.1; GHOST = 50000.0


def load(n):
    return json.load(open(os.path.join(FX, n)))


# ---------- EmeraldBay: (sample x line) spectrum from the per-condition sufficient statistics ----------
def build_emeraldbay_spectrum():
    z = np.load(os.path.join(EB, "pseudobulk.npz"), allow_pickle=True)
    ck = np.array([str(k) for k in z["cond_keys"]]); sums = z["sums"]; cnt = z["counts"]
    smp = np.array([k.split("|")[0] for k in ck]); line = np.array([k.split("|")[1] for k in ck])
    s2 = 0.938; C = 2 * (d - 1) * s2 / TH ** 2                     # 9,192 (isotropic, EmeraldBay's own C)
    ls, lc = {}, {}
    for i, l in enumerate(line):                                   # per-line mean (large near-noiseless control)
        ls[l] = ls.get(l, 0) + sums[i]; lc[l] = lc.get(l, 0) + cnt[i]
    base = np.array([ls[l] / lc[l] for l in line])
    mu = sums / cnt[:, None]; v = mu - base
    m = np.linalg.norm(v, axis=1)                                  # raw magnitude vs per-line mean
    nstar = C / np.maximum(m ** 2, 1e-12)
    over = 100 * float(np.mean(nstar < cnt))
    ghost = 100 * float(np.mean(nstar > GHOST))
    under = 100 - over - ghost
    det = 100 * float(np.mean(m ** 2 > 1.25 * d * s2 * (1.0 / cnt)))   # SNR floor (one-arm, large control)
    s2c = json.load(open(os.path.join(EB, "sample2cond.json")))
    treatments = set()
    for smp_id in set(smp):
        vraw = s2c.get(smp_id)
        if vraw is None: continue
        try: treatments.add("+".join(sorted(x[0] for x in ast.literal_eval(vraw))))
        except: treatments.add(str(vraw))
    out = dict(
        description="EmeraldBay sufficiency spectrum, (sample x line) = (treatment x line) unit, all "
                    "non-empty groups (N0>=1), isotropic n*=9,192/m^2. From eb_work pseudobulk.npz.",
        unit="(drug x dose x line)  [sample = (drug/combo, dose)]",
        n_groups=int(len(ck)), n_lines=int(len(set(line))), n_samples=int(len(set(smp))),
        n_distinct_treatments=int(len(treatments)), sigma2=s2, C_isotropic=round(C),
        detectable=round(det, 1), not_detectable=round(100 - det, 1),
        over=round(over, 1), under=round(under, 1), ghost=round(ghost, 1),
        median_m=round(float(np.median(m)), 2), median_N0=int(np.median(cnt)),
        per_group=dict(m=[round(float(x), 4) for x in m], N=[int(x) for x in cnt],
                       nstar=[round(float(x), 1) for x in nstar]))
    json.dump(out, open(os.path.join(FX, "emeraldbay_spectrum.json"), "w"), indent=1)
    return out


# ---------- assemble the six-screen table ----------
def rows():
    R = []
    c = load("tahoe_constants.json"); df = pd.read_csv(os.path.join(FX, "tahoe_quota_per_condition.csv"))
    det_t = 100 * float(np.mean(df.m.values ** 2 > 1.25 * d * c["sigma2"] * (2.0 / df.N0.values)))
    R.append(dict(name="Tahoe-100M", unit="drug × dose × line", n=len(df), det=round(det_t, 1),
                  notdet=round(100 - det_t, 1), over=c["pct_OVER"], under=c["pct_UNDER"], ghost=str(c["pct_Ghost"]),
                  denom="all", s2=c["sigma2"], medm=c["median_m"], medd=c["median_N0"]))
    e = load("emeraldbay_spectrum.json")
    R.append(dict(name="EmeraldBay", unit="drug × dose × line", n=e["n_groups"], det=e["detectable"],
                  notdet=e["not_detectable"], over=e["over"], under=e["under"], ghost=str(e["ghost"]),
                  denom="all", s2=e["sigma2"], medm=e["median_m"], medd=e["median_N0"]))
    g = load("genetic_two_arm.json")
    for key, disp, summ in [("orion_HCT116", "Orion HCT116", "orion_HCT116_summary.json"),
                            ("orion_HEK293T", "Orion HEK293T", "orion_HEK293T_summary.json"),
                            ("trade_jurkat", "TRADE Jurkat", "trade_jurkat_summary.json"),
                            ("trade_hepg2", "TRADE HepG2", "trade_hepg2_summary.json")]:
        v = g[key]; s = load(summ)
        R.append(dict(name=disp, unit="gene × line", n=v["n_perturbations"], det=v["pct_detectable"],
                      notdet=v["pct_not_detectable"], over=v["det_pct_over"], under=v["det_pct_under"],
                      ghost="n/a", denom="detectable", s2=v["sigma2_within"],
                      medm=s["m_median"], medd=s.get("median_true_N0", s.get("median_N0"))))
    return R


def write_table(R):
    hdr = ("| screen | unit | n | det. % | not-det. % | over % | under % | ghost % | denom. | σ² | "
           "med. m | med. depth |")
    sep = "|" + "|".join("-" * w for w in [13, 11, 6, 6, 7, 6, 6, 6, 11, 6, 6, 8]) + "|"
    out = [hdr, sep]
    for r in R:
        out.append(f"| {r['name']} | {r['unit']} | {r['n']:,} | {r['det']} | {r['notdet']} | {r['over']} | "
                   f"{r['under']} | {r['ghost']} | {r['denom']} | {r['s2']} | {r['medm']} | {r['medd']:,} |")
    md = "\n".join(out)
    open(os.path.join(FX, "spectrum_summary.md"), "w").write(md + "\n")
    print(md)


if __name__ == "__main__":
    build_emeraldbay_spectrum()
    write_table(rows())
