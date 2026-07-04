"""
Pass 2 (TRADE): parameter-free downsample-and-measure falsification -- the genetic-modality counterpart
of the manuscript's Phase C. Because the TRADE lines are ~10^5 cells, per-cell PCA coordinates fit in
memory (saved by pass1), so we can run the held-out subsampling test the Orion streaming run could not.

For each sufficiently sampled knockdown (N >= --min-N cells):
  * ground-truth direction u from ALL N cells against the NTC centroid: u = (mu_N - mu_NTC)/||.||
  * for a geometric grid of depths n <= N/2, draw R subsamples WITHOUT replacement, form the centroid,
    its direction to NTC, and the angle to u; take the RMS -> realized theta^2(n).
  * the closed form predicts theta^2(n) = tr(P Sigma P)/m^2 * (1/n - 1/N)  (finite-population), a line
    through the origin in (1/n - 1/N) whose slope is the A-PRIORI tr(P Sigma P)/m^2 (NO fitted parameter).
Regress realized theta^2 on (1/n - 1/N); compare the fitted slope to the a-priori slope. Stratify by the
theory's validity criterion rho^2 = m^2/(u^T Sigma u) >= 3 (first-order regime), exactly as Phase C does.

Outputs -> fixtures/trade_<line>_falsification.json  (+ per-knockdown rows).
"""
import os, json, argparse, numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", required=True)
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/trade_work/out")
    ap.add_argument("--out-fixtures", default="fixtures")
    ap.add_argument("--min-N", type=int, default=120)
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--rho2-min", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    b = np.load(os.path.join(args.dir, args.line, "basis.npz"), allow_pickle=True)
    coords = b["coords"].astype(np.float64); gene = b["gene"].astype(str)
    Sig = b["Sigma"]; trSig = float(np.trace(Sig)); ntc = str(b["ntc_label"])
    rng = np.random.default_rng(args.seed)

    idx_by = {}
    for i, g in enumerate(gene):
        idx_by.setdefault(g, []).append(i)
    ntc_idx = np.array(idx_by[ntc]); mu_ntc = coords[ntc_idx].mean(0); n_ntc = ntc_idx.size

    rows = []
    for g, ii in idx_by.items():
        if g == ntc:
            continue
        ii = np.array(ii); N = ii.size
        if N < args.min_N:
            continue
        C = coords[ii]; muN = C.mean(0)
        v = muN - mu_ntc; m_raw = float(np.linalg.norm(v))
        if m_raw == 0:
            continue
        u = v / m_raw
        trPSP = float(trSig - u @ Sig @ u)
        uSu = float(u @ Sig @ u)
        # bias-correct m (full-N centroid vs large NTC pool)
        trS_full = trSig * (1.0 / N + 1.0 / n_ntc)
        m2 = max(m_raw ** 2 - trS_full, 1e-6)
        rho2 = m2 / (uSu * (1.0 / N))          # along-signal SNR at depth N (uses S=Sigma/N ref scale)
        pred_slope = trPSP / m2

        # geometric depth grid capped at N/2
        grid = np.unique(np.geomspace(20, max(21, N // 2), 8).astype(int))
        grid = grid[grid < N]
        xs, ys = [], []
        for n in grid:
            ang2 = np.empty(args.reps)
            for r in range(args.reps):
                sub = ii[rng.choice(N, size=n, replace=False)]
                mu = coords[sub].mean(0)
                w = mu - mu_ntc; nw = np.linalg.norm(w)
                if nw == 0:
                    ang2[r] = 0.0; continue
                cos = np.clip((w / nw) @ u, -1.0, 1.0)
                ang2[r] = np.arccos(cos) ** 2
            xs.append(1.0 / n - 1.0 / N); ys.append(ang2.mean())
        xs = np.array(xs); ys = np.array(ys)
        # slope through origin: minimize ||y - s x|| -> s = <x,y>/<x,x>
        fit_slope = float((xs @ ys) / (xs @ xs))
        ss_res = float(((ys - fit_slope * xs) ** 2).sum())
        ss_tot = float(((ys - ys.mean()) ** 2).sum())
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rows.append(dict(gene_target=g, N=int(N), m=float(np.sqrt(m2)), rho2=float(rho2),
                         pred_slope=pred_slope, fit_slope=fit_slope, ratio=fit_slope / pred_slope, r2=r2))

    R = rows
    inreg = [r for r in R if r["rho2"] >= args.rho2_min]
    strong = [r for r in inreg if r["m"] > 4.0]      # clean first-order regime (strong perturbations)
    def med(xs): return float(np.median(xs)) if xs else None
    # predicted low-SNR breakdown: (1 - ratio) should grow with 1/rho^2
    try:
        from scipy.stats import spearmanr
        spear = float(spearmanr([1 - r["ratio"] for r in inreg], [1.0 / r["rho2"] for r in inreg]).correlation) if len(inreg) > 5 else None
    except Exception:
        spear = None
    summ = dict(
        line=args.line, n_knockdowns_tested=len(R), min_N=args.min_N, reps=args.reps,
        n_in_regime=len(inreg), rho2_min=args.rho2_min,
        median_ratio_in_regime=med([r["ratio"] for r in inreg]),
        median_r2_in_regime=med([r["r2"] for r in inreg]),
        n_strong=len(strong),
        median_ratio_strong=med([r["ratio"] for r in strong]),
        median_r2_strong=med([r["r2"] for r in strong]),
        median_ratio_all=med([r["ratio"] for r in R]),
        spearman_1mratio_vs_inv_rho2=spear,   # positive => predicted low-SNR breakdown reproduced
        examples=sorted(inreg, key=lambda r: -r["m"])[:10],
    )
    os.makedirs(args.out_fixtures, exist_ok=True)
    json.dump(dict(summary=summ, rows=R),
              open(os.path.join(args.out_fixtures, f"trade_{args.line}_falsification.json"), "w"), indent=1)
    print(f"\n=== TRADE FALSIFICATION (Phase C, genetic): {args.line} ===")
    print(f"tested {len(R)} knockdowns (N>={args.min_N}); in-regime (rho^2>={args.rho2_min}): {len(inreg)}")
    print(f"in-regime median fitted/predicted slope ratio = {summ['median_ratio_in_regime']}  "
          f"(R^2 median {summ['median_r2_in_regime']})   [parameter-free]")
    print(f"STRONG (m>4, n={summ['n_strong']}): median slope ratio = {summ['median_ratio_strong']} "
          f"(R^2 {summ['median_r2_strong']}) -- clean first-order regime")
    print(f"whole-population median ratio = {summ['median_ratio_all']}; "
          f"Spearman(1-ratio, 1/rho^2) = {summ['spearman_1mratio_vs_inv_rho2']} (predicted breakdown)")
    if inreg:
        print("strong in-regime examples (gene, N, m, ratio, R^2):")
        for r in sorted(inreg, key=lambda r: -r["m"])[:8]:
            print(f"  {r['gene_target']:10s} N={r['N']:5d} m={r['m']:.2f} ratio={r['ratio']:.3f} R2={r['r2']:.4f}")
    print(f"wrote fixtures/trade_{args.line}_falsification.json")


if __name__ == "__main__":
    main()
