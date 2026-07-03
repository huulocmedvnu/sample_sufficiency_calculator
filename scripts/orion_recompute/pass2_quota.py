"""
Pass 2 (pilot): compute the per-perturbation effect magnitude m = || mu_knockdown - mu_NTC || and the
sample-sufficiency quota n* on the X-Atlas/Orion head-block, then tally the OVER/UNDER/Ghost spectrum.

Control arm = the pooled Non-Targeting (NTC) centroid (the genetic analog of Tahoe's plate-matched DMSO).
Direction u = v/m is the knockdown's transcriptional heading; the quota asks how many cells are needed
to resolve it to tolerance theta*.

  n*_iso   = 2 (d-1) sigma2_within / (m^2 theta*^2)
  n*_aniso = 2 tr(P Sigma P) / (m^2 theta*^2),   tr(P Sigma P) = sum_k (1 - u_k^2) ell_k   (within-cond diag)

Regime is classified against each perturbation's TRUE total cell count N0 (from the full cells table),
not the block count -- the block only supplies enough cells to estimate m for the well-sampled subset.

Outputs -> $OUT/<LINE>/quota_pilot.csv (per perturbation) and quota_pilot_summary.json.
"""
import os, argparse, json, numpy as np, pandas as pd
import lance

HF = "hf://datasets/slaf-project/X-Atlas-Orion"
NTC = "Non-Targeting"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HCT116")
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/orion_work/pilot")
    ap.add_argument("--theta", type=float, default=0.1)
    ap.add_argument("--min-block-cells", type=int, default=30, help="min block cells to trust an m estimate")
    ap.add_argument("--ghost", type=float, default=50000)
    args = ap.parse_args()
    d = os.path.join(args.dir, args.line)
    th = args.theta

    b = np.load(os.path.join(d, "basis.npz"))
    coords = b["coords"]; ell = b["ell_within"]; s2w = float(b["sigma2_within"])
    D = int(b["n_comps"])
    meta = pd.read_parquet(os.path.join(d, "cells_meta_aligned.parquet"))
    gt = meta["gene_target"].values

    # true per-perturbation total cell counts (full atlas), for regime classification
    cds = lance.dataset(f"{HF}/data/{args.line}/cells.lance")
    full = cds.to_table(columns=["gene_target"]).to_pandas()
    N0_true = full.value_counts("gene_target")

    # NTC centroid (control arm)
    ntc_mask = gt == NTC
    mu_c = coords[ntc_mask].mean(0)
    n_ntc_block = int(ntc_mask.sum())
    trSigma = float(ell.sum())     # tr(Sigma_within) = total within-condition per-cell variance
    print(f"[pass2] NTC block cells: {n_ntc_block}  (true {int(N0_true.get(NTC,0)):,})  tr(Sigma)={trSigma:.2f}", flush=True)

    # null check: split NTC in half -> pseudo-magnitude is the sampling floor, should be ~0 vs real effects
    rng = np.random.default_rng(0)
    ni = np.where(ntc_mask)[0]; rng.shuffle(ni); h = ni.size // 2
    m_null = float(np.linalg.norm(coords[ni[:h]].mean(0) - coords[ni[h:]].mean(0)))

    rows = []
    for g, idx in meta.groupby("gene_target").indices.items():
        if g == NTC or idx.size < args.min_block_cells:
            continue
        k = idx.size
        mu_t = coords[idx].mean(0)
        v = mu_t - mu_c; m_raw = float(np.linalg.norm(v))
        if m_raw == 0:
            continue
        # bias-correct: E||v_hat||^2 = m^2 + tr(S), tr(S)=tr(Sigma)(1/k + 1/n_ntc); subtract the floor
        trS = trSigma * (1.0 / k + 1.0 / n_ntc_block)
        m = float(np.sqrt(max(m_raw**2 - trS, 1e-6)))
        u = v / m_raw                              # direction from the (unbiased) mean vector
        trPSP = float(ell.sum() - (u**2 * ell).sum())
        n_iso = 2 * (D - 1) * s2w / (m**2 * th**2)
        n_ani = 2 * trPSP / (m**2 * th**2)
        N0 = float(N0_true.get(g, k))
        rows.append(dict(gene_target=g, block_cells=int(k), N0=N0, m=m, m_raw=m_raw,
                         n_star_iso=n_iso, n_star_aniso=n_ani, trPSP=trPSP,
                         regime=("OVER" if n_ani < N0 else ("Ghost" if n_ani > args.ghost else "UNDER"))))

    df = pd.DataFrame(rows).sort_values("m", ascending=False).reset_index(drop=True)
    df.to_csv(os.path.join(d, "quota_pilot.csv"), index=False)

    ns = df["n_star_aniso"].values; N0 = df["N0"].values
    # pilot detection floor on m: with only ~k block cells/perturbation, |v_hat| under pure noise ~ sqrt(tr(S)).
    # Magnitudes below this floor are unreliable, so the weak-majority spectrum is NOT quantitative from the
    # block -- only knockdowns whose signal clears the floor have a trustworthy m. The robust statement is the
    # detectable subset. The full 18k-gene spectrum needs the full expression stream (per-perturbation ~150 cells).
    floor = np.sqrt(trSigma*(1.0/df["block_cells"].values + 1.0/n_ntc_block))
    det = df["m_raw"].values > 1.5*floor
    dsub = df[det]
    over = float((ns < N0).mean()*100); ghost = float((ns > args.ghost).mean()*100)
    summ = dict(
        line=args.line, theta=th, n_perturbations_scored=len(df),
        min_block_cells=args.min_block_cells,
        sigma2_within=round(s2w, 4), sigma2_marginal=round(float(b["sigma2_marginal"]), 4),
        sigma2_note="within~=marginal (0.897 vs 0.902) => most knockdowns transcriptionally weak; "
                    "matches Tahoe 0.957 / EmeraldBay 0.963 within-condition across modality",
        quota_const_iso=round(2*(D-1)*s2w/th**2, 1),
        pilot_m_detection_floor_median=round(float(np.median(floor)), 3),
        m_null_ntc_split=round(m_null, 4),
        n_detectable=int(det.sum()), pct_detectable=round(100*float(det.mean()),1),
        detectable_m_median=round(float(dsub.m.median()),3),
        detectable_nstar_median=int(dsub.n_star_aniso.median()),
        detectable_N0_median=int(dsub.N0.median()),
        detectable_deficit_x=round(float(dsub.n_star_aniso.median()/dsub.N0.median()),1),
        detectable_pct_over=round(100*float((dsub.n_star_aniso<dsub.N0).mean()),1),
        aniso_iso_ratio_median=round(float(np.median(df.n_star_aniso/df.n_star_iso)),4),
        caveat="block supplies ~7% of each perturbation's cells; weak-majority m unreliable at ~35 cells. "
               "Quantitative full-atlas OVER/UNDER/Ghost spectrum requires the full stream.",
        top_strong=dsub.sort_values("m",ascending=False).head(10)[["gene_target","m","n_star_aniso","N0","regime"]].to_dict("records"),
    )
    json.dump(summ, open(os.path.join(d, "quota_pilot_summary.json"), "w"), indent=1)
    print("\n=== ORION PILOT QUOTA SPECTRUM ({} , theta={}) ===".format(args.line, th))
    print(f"sigma^2 within={s2w:.4f}  marginal={float(b['sigma2_marginal']):.4f}  -> n* = {2*(D-1)*s2w/th**2:,.0f}/m^2")
    print(f"perturbations scored (>= {args.min_block_cells} block cells): {len(df)}; "
          f"pilot m-detection floor median {np.median(floor):.2f} (NTC-split null {m_null:.2f})")
    print(f"clearly-detectable knockdowns (signal>1.5x floor): {int(det.sum())} ({100*det.mean():.0f}%)")
    print(f"  of those: m median {dsub.m.median():.2f}, n*_aniso median {dsub.n_star_aniso.median():,.0f}, "
          f"true N0 median {dsub.N0.median():.0f} -> {dsub.n_star_aniso.median()/dsub.N0.median():.1f}x deficit; "
          f"OVER {100*(dsub.n_star_aniso<dsub.N0).mean():.0f}%")
    print(f"aniso/iso ratio median {np.median(df.n_star_aniso/df.n_star_iso):.3f}")
    print("\ntop detectable knockdowns:")
    print(dsub.sort_values('m',ascending=False).head(10)[["gene_target","block_cells","N0","m","n_star_aniso","regime"]].to_string(index=False))
    print(f"\nwrote {d}/quota_pilot.csv + quota_pilot_summary.json")

if __name__ == "__main__":
    main()
