"""
Pass 3 (FULL SCALE): from the streamed pseudobulk sufficient statistics, compute the genome-wide
sample-sufficiency quota spectrum for one cell line.

  within-condition sigma^2  = pooled Sum_c (csq_c - csum_c^2/n_c) / Sum_c (n_c - 1)   [per-PC -> ell_k]
  marginal sigma^2          = global per-PC variance over all cells
  m_g = || mu_g - mu_NTC ||,  bias-corrected: m^2 <- max(0, m_raw^2 - tr(Sigma)(1/n_g + 1/n_NTC))
  n*_iso   = 2 (d-1) sigma2_within / (m^2 theta^2)
  n*_aniso = 2 tr(P Sigma P) / (m^2 theta^2),  tr(P Sigma P) = sum_k (1 - u_k^2) ell_k
Regime vs the perturbation's TRUE acquired cell count n_g (= full pseudobulk count). Ghost if n* > 50k.

Outputs -> fixtures/orion_<LINE>_quota.csv (per gene) and fixtures/orion_<LINE>_summary.json.
"""
import os, sys, argparse, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from quota_arm import quota_two_arm, classify   # explicit two-arm quota + 4-class regime

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

    # within-condition pooled per-PC variance (ell_k) over perturbations + NTC with >=2 cells
    ok = cnt >= 2
    ss = (csq[ok] - csum[ok]**2 / cnt[ok, None]).sum(0)
    dof = (cnt[ok] - 1).sum()
    ell = ss / dof
    sigma2_within = float(ell.mean())
    trSig = float(ell.sum())

    # marginal per-PC variance over all cells (global)
    N = cnt.sum(); gsum = csum.sum(0); gsq = csq.sum(0)
    per_pc_marg = gsq / N - (gsum / N)**2
    sigma2_marg = float(per_pc_marg.mean())

    ntc_i = int(np.where(names == NTC)[0][0])
    n_ntc = cnt[ntc_i]; mu_ntc = csum[ntc_i] / n_ntc

    rows = []
    for i, g in enumerate(names):
        if g == NTC or cnt[i] < args.min_cells:
            continue
        n_g = cnt[i]; mu = csum[i] / n_g
        v = mu - mu_ntc; m_raw = float(np.linalg.norm(v))
        trS = trSig * (1.0 / n_g + 1.0 / n_ntc)
        m2 = max(m_raw**2 - trS, 1e-6); m = float(np.sqrt(m2))
        u = v / m_raw
        trPSP = float(ell.sum() - (u**2 * ell).sum())
        n_iso = 2 * (D - 1) * sigma2_within / (m2 * th**2)
        n_ani = 2 * trPSP / (m2 * th**2)                       # equal-arm (matched-vehicle) reference
        # EXACT two-arm quota: uses the large pooled NTC (n_ntc) consistently with the trS above,
        # instead of discarding it. n_c -> inf gives factor 1; n_c = n_t gives the equal-arm value.
        n_two = quota_two_arm(trPSP, m, th, float(n_ntc))
        rows.append(dict(gene_target=g, n_cells=int(n_g), m=m, m_raw=m_raw,
                         n_star_iso=n_iso, n_star_aniso=n_ani, n_star_two_arm=float(n_two),
                         snr_floor=m_raw / np.sqrt(trS),
                         regime=classify(n_two, n_g, args.ghost)))
    df = pd.DataFrame(rows).sort_values("m", ascending=False).reset_index(drop=True)
    os.makedirs(args.out_fixtures, exist_ok=True)
    df.to_csv(os.path.join(args.out_fixtures, f"orion_{args.line}_quota.csv"), index=False)

    ns = df["n_star_aniso"].values; ncell = df["n_cells"].values
    over = float((ns < ncell).mean()*100); ghost = float((ns > args.ghost).mean()*100); under = 100-over-ghost
    det = df["snr_floor"].values > 1.5      # signal clears the per-perturbation sampling floor
    summ = dict(
        line=args.line, theta=th, format_version="full-atlas",
        n_cells_total=int(N), n_gene_targets=int((cnt>=args.min_cells).sum()-0),
        n_perturbations_scored=len(df), min_cells=args.min_cells,
        n_ntc_cells=int(n_ntc),
        sigma2_within=round(sigma2_within,4), sigma2_marginal=round(sigma2_marg,4),
        quota_const_iso=round(2*(D-1)*sigma2_within/th**2,1),
        m_median=round(float(np.median(df.m)),3),
        m_p90=round(float(np.percentile(df.m,90)),3), m_max=round(float(df.m.max()),3),
        n_star_median=int(np.median(ns)),
        pct_over=round(over,1), pct_under=round(under,1), pct_ghost=round(ghost,1),
        pct_detectable=round(100*float(det.mean()),1),
        detectable_pct_over=round(100*float((ns[det]<ncell[det]).mean()),1) if det.any() else None,
        aniso_iso_ratio_median=round(float(np.median(df.n_star_aniso/df.n_star_iso)),4),
        median_true_N0=int(np.median(ncell)),
        median_deficit_x=round(float(np.median(ns/ncell)),2),
        top_strong=df.head(12)[["gene_target","m","n_star_aniso","n_cells","regime"]].to_dict("records"),
    )
    json.dump(summ, open(os.path.join(args.out_fixtures, f"orion_{args.line}_summary.json"), "w"), indent=1)

    print(f"\n=== ORION FULL-ATLAS QUOTA SPECTRUM: {args.line} (theta={th}) ===")
    print(f"cells {int(N):,}  perturbations scored {len(df):,}  NTC cells {int(n_ntc):,}")
    print(f"sigma^2 within={sigma2_within:.4f}  marginal={sigma2_marg:.4f}  -> n* = {2*(D-1)*sigma2_within/th**2:,.0f}/m^2")
    print(f"m: median {np.median(df.m):.3f}  p90 {np.percentile(df.m,90):.2f}  max {df.m.max():.2f}")
    print(f"n*_aniso median {np.median(ns):,.0f}   aniso/iso ratio {np.median(df.n_star_aniso/df.n_star_iso):.3f}")
    print(f"SPECTRUM vs true acquired cells:  {over:.1f}% OVER / {under:.1f}% UNDER / {ghost:.1f}% Ghost")
    print(f"  (median acquired {np.median(ncell):.0f} cells, median deficit {np.median(ns/ncell):.1f}x; "
          f"detectable {100*det.mean():.0f}%, of which OVER {100*(ns[det]<ncell[det]).mean():.0f}%)")
    print("top strong knockdowns:")
    print(df.head(12)[["gene_target","n_cells","m","n_star_aniso","regime"]].to_string(index=False))
    print(f"\nwrote fixtures/orion_{args.line}_quota.csv + summary.json")

if __name__ == "__main__":
    main()
