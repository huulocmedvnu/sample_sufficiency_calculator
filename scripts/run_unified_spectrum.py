"""
Run the unified engine over all six screens and write the canonical fixtures + Table 5.

ONE engine (src/engine.py), TWO data pipelines (src/pipelines/{chemical,genetic}.py). Fixes audit
findings #1-#9: single computation, real n_c, bias-correction + snr>1.5 detection for every dataset,
full two-arm quota (no hardcoded arm factor), spectrum over all conditions (no cell-count filter).

Chemical HEADLINE control = per-line mean (large pool), the same reference the held-out test uses.
Tahoe's real vehicle-matched shared DMSO control is reported separately (finding #9).

Writes:
  fixtures/unified_spectrum.json   aggregate summary for all six screens (+ Tahoe shared-DMSO finding)
  fixtures/unified_spectrum_per_condition.npz   per-condition m_raw/m_corr/n_t/n_c/n_star/regime
  fixtures/spectrum_summary.md     Table 5 (two blocks)
"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from engine import compute, summarize
from pipelines import chemical, genetic

FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
UNIT = {"Tahoe-100M": "drug × dose × line", "EmeraldBay": "sample × line (= drug × dose × line)",
        "Orion HCT116": "gene × line", "Orion HEK293T": "gene × line",
        "TRADE Jurkat": "gene × line", "TRADE HepG2": "gene × line"}


def run(d):
    r = compute(d["mu_t"], d["mu_c"], d["n_t"], d["n_c"], d["Sigma"])
    s = summarize(r); s["name"] = d["dataset"]; s["modality"] = d["modality"]
    s["control_type"] = d["control_type"]; s["unit"] = UNIT.get(d["dataset"], "")
    return s, r


def main():
    # HEADLINE: Tahoe references its REAL vehicle-matched shared DMSO pool (gap 2 decision).
    specs = [chemical.load("Tahoe-100M", control="shared-dmso"), chemical.load("EmeraldBay"),
             genetic.load("orion_HCT116"), genetic.load("orion_HEK293T"),
             genetic.load("trade_jurkat"), genetic.load("trade_hepg2")]
    summaries = []; percond = {}
    for d in specs:
        s, r = run(d); summaries.append(s)
        key = s["name"].replace(" ", "_")
        percond[f"{key}_m_raw"] = r["m_raw"]; percond[f"{key}_m_corr"] = r["m_corr"]
        percond[f"{key}_n_t"] = r["n_t"]; percond[f"{key}_n_c"] = r["n_c"]
        percond[f"{key}_n_star"] = np.where(np.isfinite(r["n_star"]), r["n_star"], -1.0)
        percond[f"{key}_regime"] = np.array([str(x) for x in r["regime"]])

    # sensitivity: Tahoe referenced to the per-line mean (large pool) instead of the true vehicle
    ssens, _ = run(chemical.load("Tahoe-100M", control="per-line-mean"))
    ssens["name"] = "Tahoe-100M (per-line-mean sensitivity)"

    out = dict(headline=summaries, sensitivity=ssens,
               note="Unified engine, REAL diagonal within-Sigma. HEADLINE Tahoe control = real shared "
                    "DMSO vehicle (n_c~1,500, control-pool-limited); sensitivity = per-line mean (large pool).")
    json.dump(out, open(os.path.join(FX, "unified_spectrum.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(FX, "unified_spectrum_per_condition.npz"), **percond)

    # ---- Table 5: two blocks (detectability sums to 100; sufficiency over the denominator) ----
    hdr = ("| screen | unit | n | det. % | not-det. % | over % | under % | ghost % | pool-lim % | denom. | "
           "\\(\\sigma^2\\) | med. \\(m\\) | med. depth | C bound |")
    sep = "|" + "|".join("-" * w for w in [13, 11, 6, 6, 7, 6, 6, 6, 7, 9, 6, 6, 8, 6]) + "|"
    lines = [hdr, sep]
    for s in summaries:
        denom = "all conditions"
        lines.append(f"| {s['name']} | {s['unit']} | {s['n_conditions']:,} | {s['detectable_pct']} | "
                     f"{s['not_detectable_pct']} | {s['pct_over']} | {s['pct_under']} | {s['pct_ghost']} | "
                     f"{s['pct_pool_limited']} | {denom} | {s['sigma2']} | {s['median_m_corr']} | "
                     f"{s['median_n_t']:,} | {s['C_largepool']:,} |")
    md = "\n".join(lines)
    open(os.path.join(FX, "spectrum_summary.md"), "w").write(md + "\n")
    print(md)
    print("\n=== Tahoe per-line-mean sensitivity ===")
    print(f"  det {ssens['detectable_pct']}%  over {ssens['pct_over']}  under {ssens['pct_under']}  "
          f"ghost {ssens['pct_ghost']}  pool {ssens['pct_pool_limited']}")


if __name__ == "__main__":
    main()
