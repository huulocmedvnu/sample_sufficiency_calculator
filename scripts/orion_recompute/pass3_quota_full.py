"""
Pass 3 (FULL SCALE): genome-wide sample-sufficiency spectrum for one X-Atlas/Orion line, computed
THROUGH THE UNIFIED ENGINE (src/engine.py).

Reads the committed pseudobulk sufficient statistics (gt_names, cnt, csum, csq), assembles the standard
per-condition structure -- each knockdown centroid mu_g vs the pooled Non-Targeting centroid mu_NTC,
the real acquired cell counts (n_t = per-gene count, n_c = NTC pool), and the within-condition diagonal
Sigma = diag(ell_k) -- and hands it to engine.compute. The bias-correction m^2 - tr(S), the snr>1.5
detection floor, the two-arm quota n* (large-pool limit for the big NTC pool) and the regime all come
from the engine; there is NO local quota/bias/detection formula here (this file used to duplicate it,
which is exactly the drift the audit removed).

Outputs -> fixtures/orion_<LINE>_quota.csv (per gene) and fixtures/orion_<LINE>_summary.json.
"""
import os, sys, argparse, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute, summarize

NTC = "Non-Targeting"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", required=True)
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/orion_work/full")
    ap.add_argument("--theta", type=float, default=0.1)
    ap.add_argument("--ghost", type=float, default=50000)
    ap.add_argument("--min-cells", type=int, default=25, help="min acquired cells to score a perturbation")
    ap.add_argument("--out-fixtures", default="fixtures")
    args = ap.parse_args()
    d = os.path.join(args.dir, args.line)
    z = np.load(os.path.join(d, "pseudobulk.npz"), allow_pickle=True)
    names = z["gt_names"]; cnt = z["cnt"].astype(np.float64)
    csum = z["csum"].astype(np.float64); csq = z["csq"].astype(np.float64)
    D = int(z["n_comps"]); th = args.theta

    # --- sufficient statistics only (no quota math): within-condition per-PC variance ell_k ---
    ok = cnt >= 2
    ell = (csq[ok] - csum[ok] ** 2 / cnt[ok, None]).sum(0) / (cnt[ok] - 1).sum()
    sigma2_within = float(ell.mean())
    Ncells = cnt.sum(); gsum = csum.sum(0); gsq = csq.sum(0)
    sigma2_marg = float((gsq / Ncells - (gsum / Ncells) ** 2).mean())

    ntc_i = int(np.where(names == NTC)[0][0])
    n_ntc = float(cnt[ntc_i]); mu_ntc = csum[ntc_i] / n_ntc
    sel = [i for i, g in enumerate(names) if g != NTC and cnt[i] >= args.min_cells]
    mu_t = np.array([csum[i] / cnt[i] for i in sel])
    n_t = cnt[list(sel)]
    Sigma = np.diag(ell)                                          # within-condition diagonal Sigma

    # --- ALL sufficiency math via the single engine ---
    r = compute(mu_t, np.tile(mu_ntc, (len(sel), 1)), n_t, np.full(len(sel), n_ntc), Sigma,
                theta=th, ghost=args.ghost)
    s = summarize(r)

    df = pd.DataFrame(dict(
        gene_target=[str(names[i]) for i in sel], n_cells=n_t.astype(int),
        m_raw=np.round(r["m_raw"], 4), m_corr=np.round(r["m_corr"], 4), snr=np.round(r["snr"], 4),
        detectable=r["detectable"], n_star=r["n_star"], regime=r["regime"],
    )).sort_values("m_corr", ascending=False).reset_index(drop=True)
    os.makedirs(args.out_fixtures, exist_ok=True)
    df.to_csv(os.path.join(args.out_fixtures, f"orion_{args.line}_quota.csv"), index=False)

    summ = dict(
        line=args.line, theta=th, format_version="engine-unified",
        n_cells_total=int(Ncells), n_perturbations_scored=len(df), min_cells=args.min_cells,
        n_ntc_cells=int(n_ntc),
        sigma2_within=round(sigma2_within, 4), sigma2_marginal=round(sigma2_marg, 4),
        C_largepool=s["C_largepool"],
        m_median=round(float(np.median(r["m_corr"])), 3),
        detectable_pct=s["detectable_pct"], not_detectable_pct=s["not_detectable_pct"],
        det_pct_over=s["det_pct_over"], det_pct_under=s["det_pct_under"],
        pct_over=s["pct_over"], pct_under=s["pct_under"], pct_ghost=s["pct_ghost"],
        pct_pool_limited=s["pct_pool_limited"],
        median_n_star=s["median_n_star"], median_true_N0=int(np.median(n_t)),
        top_strong=df.head(12)[["gene_target", "m_corr", "n_star", "n_cells", "regime"]].to_dict("records"),
    )
    json.dump(summ, open(os.path.join(args.out_fixtures, f"orion_{args.line}_summary.json"), "w"), indent=1)

    print(f"\n=== ORION FULL-ATLAS SPECTRUM (engine): {args.line} (theta={th}) ===")
    print(f"cells {int(Ncells):,}  scored {len(df):,}  NTC {int(n_ntc):,}  "
          f"sigma^2 within={sigma2_within:.4f} marginal={sigma2_marg:.4f}  C={s['C_largepool']:,}")
    print(f"detectable {s['detectable_pct']}%  over(det) {s['det_pct_over']}%  under(det) {s['det_pct_under']}%  "
          f"median n* {s['median_n_star']}")
    print(f"wrote fixtures/orion_{args.line}_quota.csv + summary.json")


if __name__ == "__main__":
    main()
