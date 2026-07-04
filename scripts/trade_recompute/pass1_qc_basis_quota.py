"""
Pass 1 (TRADE essential-gene lines): QC -> GLOBAL embedding -> cell-level within-condition covariance
Sigma -> LATE pseudobulk quota. Faithful to the manuscript's load-bearing pipeline topology
(docs/THEORY.md; "Pipeline topology (order of operations is load-bearing)"):

  (i)  Global embedding FIRST: normalize_total(1e4) -> log1p -> HVG(2000) -> PCA(50) fit on the POOLED
       single-cell matrix (after QC), defining one shared d=50 coordinate space.
  (ii) Cell-level covariance: Sigma in R^{50x50} estimated from INDIVIDUAL embedded cells by centering
       each condition on its own mean and pooling the residuals (the within-condition noise that enters
       a centroid). We keep the FULL 50x50 Sigma (not just its diagonal) for tr(P Sigma P).
  (iii) Late pseudobulk: only afterwards collapse each condition to a centroid (mu_t knockdown, mu_c NTC).

QC (exact, per TRADE / user spec; --line sets the constants):
  * low-UMI:   drop cells with UMI_count <  UMI_FRAC * median(UMI_count)   (Jurkat 0.14, HepG2 0.18)
  * high-mito: drop cells with mito UMIs  >  MITO_UMI_MAX                  (Jurkat 1750, HepG2 3000)
               mito UMIs = mitopercent(as fraction) * UMI_count  (mitopercent auto-detected %/fraction)
  * guide:     retain cells with a single sgRNA OR two sgRNAs targeting the SAME gene (parsed from
               sgID_AB); drop unassigned / multi-gene cells.

Effect magnitude m = || mu_knockdown - mu_NTC ||, bias-corrected for the finite-sample floor.
n*_iso   = 2 (d-1) sigma2_within / (m^2 theta^2);  n*_aniso = 2 tr(P Sigma P) / (m^2 theta^2)
(equal-arms factor 2, matching the Tahoe/Orion headline convention; a large shared NTC pool would
halve it -- reported in the summary as the large-control-pool variant.)

Outputs -> $OUT/<line>/{basis.npz, quota.csv, summary.json}
  basis.npz keeps per-cell PCA coords + gene labels so pass2 (falsification) needs no re-embed.
"""
import os, re, json, argparse, numpy as np, scipy.sparse as sp
import scanpy as sc, anndata as ad, pandas as pd

# per-line QC constants (TRADE)
LINE_QC = {
    "jurkat": dict(umi_frac=0.14, mito_umi_max=1750.0),
    "hepg2":  dict(umi_frac=0.18, mito_umi_max=3000.0),
}
TARGET_SUM = 1e4
N_COMPS = 50
NTC_RE = re.compile(r"non.?target|^ntc$|control|safe.?harbor|scramble|negative", re.I)


def resolve_cols(obs):
    """Map obs to the fields we need, tolerant to naming."""
    cols = {c.lower(): c for c in obs.columns}
    def pick(*names):
        for n in names:
            if n in obs.columns: return n
            if n.lower() in cols: return cols[n.lower()]
        return None
    return dict(
        umi=pick("UMI_count", "total_counts", "n_counts", "num_umis"),
        mito=pick("mitopercent", "percent_mito", "pct_counts_mt", "mito_frac"),
        gene=pick("gene", "gene_target", "target_gene", "target"),
        guide=pick("sgID_AB", "guide_identity", "sgID", "sgRNA", "protospacer"),
        gem=pick("gem_group", "gemgroup", "batch", "lane"),
    )


def guide_ok(sgid_series, gene_series):
    """Retain single-sgRNA cells or dual-sgRNA cells targeting the SAME gene.
    Dual guides in these libraries are encoded like 'GENE_+_1|GENE_+_2' (two tokens split on |, ; or _AB_).
    We keep a cell if its guide tokens map to <=1 distinct target gene. Falls back to keeping any cell that
    already carries a single assigned `gene` when sgID_AB is uninformative."""
    keep = np.ones(len(sgid_series), dtype=bool)
    def target_of(tok):
        # strip trailing direction/index, e.g. 'RPL7_+_2' -> 'RPL7'
        return re.split(r"[_|]", tok.strip())[0] if tok else ""
    for i, s in enumerate(sgid_series.astype(str).values):
        toks = re.split(r"[|;]", s) if any(d in s for d in "|;") else [s]
        tgs = {target_of(t) for t in toks if t and t.lower() != "nan"}
        tgs.discard("")
        if len(tgs) > 1:
            keep[i] = False
    return keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", required=True, choices=["jurkat", "hepg2"])
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", default="/mnt/hdd2/loc-tran/trade_work/out")
    ap.add_argument("--theta", type=float, default=0.1)
    ap.add_argument("--min-cells", type=int, default=25)
    ap.add_argument("--ghost", type=float, default=50000)
    args = ap.parse_args()
    qc = LINE_QC[args.line]
    d = os.path.join(args.out, args.line); os.makedirs(d, exist_ok=True)
    th = args.theta

    print(f"[pass1:{args.line}] loading {args.h5ad}", flush=True)
    A = sc.read_h5ad(args.h5ad)
    A.var_names_make_unique()
    col = resolve_cols(A.obs)
    print(f"[pass1:{args.line}] {A.n_obs:,} cells x {A.n_vars:,} genes; resolved cols: {col}", flush=True)
    for k in ("umi", "mito", "gene"):
        if col[k] is None:
            raise SystemExit(f"could not resolve obs column for '{k}'; obs has {list(A.obs.columns)}")

    umi = np.asarray(A.obs[col["umi"]], dtype=np.float64)
    mito = np.asarray(A.obs[col["mito"]], dtype=np.float64)
    mito_frac = mito / 100.0 if np.nanmax(mito) > 1.5 else mito          # detect % vs fraction
    mito_umi = mito_frac * umi
    gene = A.obs[col["gene"]].astype(str).values

    # --- QC ---
    umi_min = qc["umi_frac"] * np.median(umi)
    keep_umi = umi >= umi_min
    keep_mito = mito_umi <= qc["mito_umi_max"]
    if col["guide"] is not None:
        keep_guide = guide_ok(A.obs[col["guide"]], A.obs[col["gene"]])
    else:
        keep_guide = np.ones(A.n_obs, bool)
        print(f"[pass1:{args.line}] WARNING: no guide column; skipping guide-assignment filter", flush=True)
    keep = keep_umi & keep_mito & keep_guide
    print(f"[pass1:{args.line}] QC: umi>= {umi_min:.0f} keeps {keep_umi.mean():.1%}; "
          f"mito_umi<= {qc['mito_umi_max']:.0f} keeps {keep_mito.mean():.1%}; "
          f"guide keeps {keep_guide.mean():.1%}; ALL -> {int(keep.sum()):,}/{A.n_obs:,} ({keep.mean():.1%})",
          flush=True)
    A = A[keep].copy()
    gene = gene[keep]

    # --- (i) GLOBAL embedding first (pooled, post-QC) ---
    if not sp.issparse(A.X):
        A.X = sp.csr_matrix(A.X)
    sc.pp.normalize_total(A, target_sum=TARGET_SUM)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="seurat")
    A = A[:, A.var["highly_variable"].values].copy()
    sc.pp.pca(A, n_comps=N_COMPS, zero_center=True, svd_solver="arpack")
    coords = A.obsm["X_pca"].astype(np.float64)                          # (n,50) per-cell
    evr = A.uns["pca"]["variance_ratio"].astype(np.float64)
    print(f"[pass1:{args.line}] embedding: {A.n_obs:,} cells, HVG 2000, PCA{N_COMPS}; top5 evr "
          f"{evr[:5].round(3).tolist()}", flush=True)

    # NTC label
    cats = pd.unique(gene)
    ntc_labels = [c for c in cats if NTC_RE.search(str(c))]
    if not ntc_labels:
        raise SystemExit(f"no non-targeting label found among gene categories (examples: {cats[:10]})")
    ntc_label = max(ntc_labels, key=lambda c: (gene == c).sum())
    is_ntc = gene == ntc_label
    print(f"[pass1:{args.line}] NTC label='{ntc_label}' ({int(is_ntc.sum()):,} cells); "
          f"{len(cats)-len(ntc_labels)} knockdown genes", flush=True)

    # --- (ii) cell-level within-condition covariance Sigma (full 50x50) + marginal ---
    order = np.argsort(gene, kind="stable")
    gs, starts = np.unique(gene[order], return_index=True)
    ends = np.r_[starts[1:], len(gene)]
    Sig = np.zeros((N_COMPS, N_COMPS)); dof = 0
    cent_sums = {}; cent_cnt = {}
    for g, s, e in zip(gs, starts, ends):
        idx = order[s:e]; n = idx.size
        C = coords[idx]; mu = C.mean(0)
        cent_sums[g] = C.sum(0); cent_cnt[g] = n
        if n >= 2:
            R = C - mu
            Sig += R.T @ R
            dof += n - 1
    Sig /= dof                                                           # pooled within-condition Sigma
    ell = np.diag(Sig).copy()
    sigma2_within = float(ell.mean())
    trSig = float(np.trace(Sig))
    gsum = coords.sum(0); gsq = (coords * coords).sum(0); N = coords.shape[0]
    sigma2_marg = float((gsq / N - (gsum / N) ** 2).mean())
    print(f"[pass1:{args.line}] sigma^2 within={sigma2_within:.4f} marginal={sigma2_marg:.4f}  "
          f"tr(Sigma)={trSig:.2f}  n* = {2*(N_COMPS-1)*sigma2_within/th**2:,.0f}/m^2", flush=True)

    # --- (iii) late pseudobulk quota ---
    mu_ntc = cent_sums[ntc_label] / cent_cnt[ntc_label]; n_ntc = cent_cnt[ntc_label]
    rows = []
    for g in gs:
        if g == ntc_label or cent_cnt[g] < args.min_cells:
            continue
        n_g = cent_cnt[g]; mu = cent_sums[g] / n_g
        v = mu - mu_ntc; m_raw = float(np.linalg.norm(v))
        if m_raw == 0:
            continue
        trS = trSig * (1.0 / n_g + 1.0 / n_ntc)              # equal-arm sampling-noise floor
        m2 = max(m_raw ** 2 - trS, 1e-6); m = float(np.sqrt(m2))
        u = v / m_raw
        trPSP = float(trSig - u @ Sig @ u)                   # full-Sigma anisotropic functional
        n_iso = 2 * (N_COMPS - 1) * sigma2_within / (m2 * th ** 2)
        n_ani = 2 * trPSP / (m2 * th ** 2)
        rows.append(dict(gene_target=g, n_cells=int(n_g), m=m, m_raw=m_raw,
                         n_star_iso=n_iso, n_star_aniso=n_ani, trPSP=trPSP,
                         snr_floor=m_raw / np.sqrt(trS),
                         regime=("OVER" if n_ani < n_g else ("Ghost" if n_ani > args.ghost else "UNDER"))))
    df = pd.DataFrame(rows).sort_values("m", ascending=False).reset_index(drop=True)
    df.to_csv(os.path.join(d, "quota.csv"), index=False)

    ns = df["n_star_aniso"].values; ncell = df["n_cells"].values
    det = df["snr_floor"].values > 1.5
    over = float((ns < ncell).mean() * 100); ghost = float((ns > args.ghost).mean() * 100)
    summ = dict(
        line=args.line, theta=th, qc=qc, ntc_label=str(ntc_label),
        n_cells_qc=int(N), n_cells_raw=int(keep.size), qc_pass_frac=round(float(keep.mean()), 4),
        n_knockdowns_scored=len(df), n_ntc_cells=int(n_ntc),
        sigma2_within=round(sigma2_within, 4), sigma2_marginal=round(sigma2_marg, 4),
        quota_const_iso=round(2 * (N_COMPS - 1) * sigma2_within / th ** 2, 1),
        quota_const_iso_largepool=round((N_COMPS - 1) * sigma2_within / th ** 2, 1),
        m_median=round(float(np.median(df.m)), 3), m_p90=round(float(np.percentile(df.m, 90)), 3),
        m_max=round(float(df.m.max()), 3),
        n_star_median=int(np.median(ns)),
        pct_over=round(over, 1), pct_under=round(100 - over - ghost, 1), pct_ghost=round(ghost, 1),
        pct_detectable=round(100 * float(det.mean()), 1),
        detectable_pct_over=round(100 * float((ns[det] < ncell[det]).mean()), 1) if det.any() else None,
        median_deficit_x=round(float(np.median(ns / ncell)), 1),
        aniso_iso_ratio_median=round(float(np.median(df.n_star_aniso / df.n_star_iso)), 4),
        median_N0=int(np.median(ncell)),
        top_strong=df.head(12)[["gene_target", "m", "n_star_aniso", "n_cells", "regime"]].to_dict("records"),
    )
    json.dump(summ, open(os.path.join(d, "summary.json"), "w"), indent=1)
    np.savez(os.path.join(d, "basis.npz"),
             coords=coords.astype(np.float32), gene=gene.astype(str),
             Sigma=Sig, ell=ell, sigma2_within=sigma2_within, sigma2_marginal=sigma2_marg,
             ntc_label=str(ntc_label), n_comps=N_COMPS, evr=evr)

    print(f"\n=== TRADE QUOTA SPECTRUM: {args.line} (theta={th}) ===")
    print(f"QC-passed {N:,} cells; {len(df):,} knockdowns; NTC {n_ntc:,}")
    print(f"sigma^2 within={sigma2_within:.4f} -> n* = {2*(N_COMPS-1)*sigma2_within/th**2:,.0f}/m^2  "
          f"(large-NTC-pool variant {(N_COMPS-1)*sigma2_within/th**2:,.0f}/m^2)")
    print(f"m median {np.median(df.m):.3f} (p90 {np.percentile(df.m,90):.2f}, max {df.m.max():.2f})")
    print(f"SPECTRUM: {over:.1f}% OVER / {100-over-ghost:.1f}% UNDER / {ghost:.1f}% Ghost  "
          f"(median n* {np.median(ns):,.0f} vs median depth {np.median(ncell):.0f}; "
          f"deficit {np.median(ns/ncell):.1f}x; detectable {100*det.mean():.0f}%, OVER {100*(ns[det]<ncell[det]).mean():.0f}%)")
    print(f"aniso/iso ratio median {np.median(df.n_star_aniso/df.n_star_iso):.3f}")
    print("top strong knockdowns:")
    print(df.head(12)[["gene_target", "n_cells", "m", "n_star_aniso", "regime"]].to_string(index=False))
    print(f"wrote {d}/quota.csv + summary.json + basis.npz")


if __name__ == "__main__":
    main()
