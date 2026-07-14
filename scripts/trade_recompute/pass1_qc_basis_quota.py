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
import os, sys, re, json, argparse, numpy as np, scipy.sparse as sp
import scanpy as sc, anndata as ad, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from engine import compute, summarize   # the SINGLE owner of the quota / bias / detection / regime

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
          f"tr(Sigma)={trSig:.2f}  (quota via src/engine.py)", flush=True)

    # --- (iii) late pseudobulk: assemble the standard structure, hand it to the SINGLE engine ---
    mu_ntc = cent_sums[ntc_label] / cent_cnt[ntc_label]; n_ntc = float(cent_cnt[ntc_label])
    sel = [g for g in gs if g != ntc_label and cent_cnt[g] >= args.min_cells]
    mu_t = np.array([cent_sums[g] / cent_cnt[g] for g in sel])
    n_t = np.array([cent_cnt[g] for g in sel], float)
    r = compute(mu_t, np.tile(mu_ntc, (len(sel), 1)), n_t, np.full(len(sel), n_ntc), Sig,
                theta=th, ghost=args.ghost)                  # bias / detection / quota / regime all engine
    s = summarize(r)
    df = pd.DataFrame(dict(
        gene_target=sel, n_cells=n_t.astype(int),
        m_raw=np.round(r["m_raw"], 4), m_corr=np.round(r["m_corr"], 4), snr=np.round(r["snr"], 4),
        detectable=r["detectable"], n_star=r["n_star"], regime=r["regime"],
    )).sort_values("m_corr", ascending=False).reset_index(drop=True)
    df.to_csv(os.path.join(d, "quota.csv"), index=False)

    summ = dict(
        line=args.line, theta=th, qc=qc, ntc_label=str(ntc_label),
        n_cells_qc=int(N), n_cells_raw=int(keep.size), qc_pass_frac=round(float(keep.mean()), 4),
        n_knockdowns_scored=len(df), n_ntc_cells=int(n_ntc),
        sigma2_within=round(sigma2_within, 4), sigma2_marginal=round(sigma2_marg, 4),
        C_largepool=s["C_largepool"],
        m_median=round(float(np.median(r["m_corr"])), 3),
        detectable_pct=s["detectable_pct"], not_detectable_pct=s["not_detectable_pct"],
        det_pct_over=s["det_pct_over"], det_pct_under=s["det_pct_under"],
        pct_over=s["pct_over"], pct_under=s["pct_under"], pct_ghost=s["pct_ghost"],
        pct_pool_limited=s["pct_pool_limited"],
        median_n_star=s["median_n_star"], median_N0=int(np.median(n_t)),
        top_strong=df.head(12)[["gene_target", "m_corr", "n_star", "n_cells", "regime"]].to_dict("records"),
    )
    json.dump(summ, open(os.path.join(d, "summary.json"), "w"), indent=1)
    np.savez(os.path.join(d, "basis.npz"),
             coords=coords.astype(np.float32), gene=gene.astype(str),
             Sigma=Sig, ell=ell, sigma2_within=sigma2_within, sigma2_marginal=sigma2_marg,
             ntc_label=str(ntc_label), n_comps=N_COMPS, evr=evr)

    print(f"\n=== TRADE SPECTRUM (engine): {args.line} (theta={th}) ===")
    print(f"QC-passed {N:,} cells; {len(df):,} knockdowns; NTC {int(n_ntc):,}  "
          f"sigma^2 within={sigma2_within:.4f}  C={s['C_largepool']:,}")
    print(f"detectable {s['detectable_pct']}%  over(det) {s['det_pct_over']}%  under(det) {s['det_pct_under']}%  "
          f"median n* {s['median_n_star']}")
    print("top strong knockdowns:")
    print(df.head(12)[["gene_target", "n_cells", "m_corr", "n_star", "regime"]].to_string(index=False))
    print(f"wrote {d}/quota.csv + summary.json + basis.npz")


if __name__ == "__main__":
    main()
