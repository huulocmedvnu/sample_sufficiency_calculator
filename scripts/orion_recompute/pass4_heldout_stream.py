"""
Pass 4 (held-out falsification, streamed): run the direct downsample-and-measure angular-error test
DIRECTLY on the X-Atlas/Orion genome-wide CRISPRi atlas, independently per line, on EVERY gene
knockdown with N >= 400 cells. This is the genetic-modality analogue of the EmeraldBay held-out
test (scripts/emeraldbay_recompute/pass4_falsification.py) and uses the IDENTICAL routine
(group_geometry + rms_angle + through-origin slope fit), so the two modalities are directly
comparable. No "shared cell line" restriction: the test runs on the atlas's own eligible groups.

Why a re-stream: pass2 kept only per-gene sufficient statistics (count, sum, sumsq), which suffice
for the quota but NOT for held-out subsampling (that needs per-cell coordinates). Here we re-stream
the full expression table once per line and RETAIN the per-cell PCA(50) coordinates ONLY for cells
whose gene_target has >= MIN_N cells (a bounded subset), then run the held-out fit per gene against
the pooled Non-Targeting centroid (the genetic reference, from the committed pseudobulk.npz).

Projection is byte-identical to pass2: coord = log1p(count/total_counts*1e4)[HVG] @ components.T
- (pca_mean @ components.T). Resumable: checkpoints (offset + retained coords) to $OUT/<LINE>/heldout_ckpt.npz.

Run:  python scripts/orion_recompute/pass4_heldout_stream.py --line HCT116 \
        --basis /mnt/hdd2/loc-tran/orion_work/full/HCT116/basis.npz \
        --pseudobulk /mnt/hdd2/loc-tran/orion_work/full/HCT116/pseudobulk.npz
Output -> fixtures/orion_<LINE>_heldout.json
"""
import os, sys, time, json, argparse
import numpy as np, scipy.sparse as sp, pandas as pd, lance

HF = "hf://datasets/slaf-project/X-Atlas-Orion"
NTC = "Non-Targeting"
TH = 0.1
MIN_N = 400          # eligibility: knockdowns with >= 400 acquired cells
REPS = 200           # resamples per depth (matches EmeraldBay pass4)
NGRID = 12           # geometric depth-grid points
rng = np.random.default_rng(0)


def group_geometry(C, base, d):
    """m, u, trPSP, rho2 for a cell block C vs baseline (identical to EmeraldBay pass4)."""
    v = C.mean(0) - base
    m = float(np.linalg.norm(v))
    if m <= 0:
        return None
    u = v / m
    P = np.eye(d) - np.outer(u, u)
    Sig = np.cov(C.T)
    trPSP = float(np.trace(P @ Sig @ P))
    uSu = float(u @ Sig @ u)
    rho2 = m ** 2 / uSu if uSu > 0 else np.inf
    return m, u, trPSP, rho2


def rms_angle(C, base, u, n, reps):
    N = len(C); acc = np.empty(reps)
    for j in range(reps):
        idx = rng.choice(N, n, replace=False)
        vn = C[idx].mean(0) - base
        un = vn / np.linalg.norm(vn)
        acc[j] = np.arccos(np.clip(un @ u, -1, 1)) ** 2
    return acc.mean()


def held_out(C, base, d):
    N = len(C)
    g = group_geometry(C, base, d)
    if g is None:
        return None
    m, u, trPSP, rho2 = g
    if m <= 0 or trPSP <= 0:
        return None
    pred_slope = trPSP / m ** 2
    ns = np.unique(np.geomspace(20, max(N // 2, 40), NGRID).round().astype(int))
    ns = ns[ns < N]
    if ns.size < 3:
        return None
    xs = np.array([1.0 / n - 1.0 / N for n in ns])
    ys = np.array([rms_angle(C, base, u, int(n), REPS) for n in ns])
    fit_slope = float((xs @ ys) / (xs @ xs))
    ss_tot = float(((ys - ys.mean()) ** 2).sum())
    r2 = 1.0 - float(((ys - fit_slope * xs) ** 2).sum()) / ss_tot if ss_tot > 0 else 1.0
    return dict(N=int(N), m=round(m, 4), rho2=round(float(rho2), 2),
                pred_slope=round(pred_slope, 6), fit_slope=round(fit_slope, 6),
                ratio=round(fit_slope / pred_slope, 4), r2=round(r2, 5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HCT116", choices=["HCT116", "HEK293T"])
    ap.add_argument("--basis", required=True)
    ap.add_argument("--pseudobulk", required=True, help="pass2 pseudobulk.npz (for NTC centroid)")
    ap.add_argument("--out", default=os.environ.get("OUT", "/mnt/hdd2/loc-tran/orion_work/full"))
    ap.add_argument("--fixtures", default=os.path.join(os.path.dirname(__file__), "..", "..", "fixtures"))
    ap.add_argument("--chunk", type=int, default=8_000_000)
    ap.add_argument("--ckpt-every", type=int, default=40)
    a = ap.parse_args()
    d = 50
    outdir = os.path.join(a.out, a.line); os.makedirs(outdir, exist_ok=True)
    ckpt = os.path.join(outdir, "heldout_ckpt.npz")

    b = np.load(a.basis)
    comps = b["components"].astype(np.float64); hvg = b["hvg_ids"]
    mean_proj = (b["pca_mean"].astype(np.float64) @ comps.T)
    n_genes_full = 38606
    hvg_col = np.full(n_genes_full, -1, np.int64); hvg_col[hvg] = np.arange(hvg.size)

    cds = lance.dataset(f"{HF}/data/{a.line}/cells.lance")
    cm = cds.to_table(columns=["cell_integer_id", "total_counts", "gene_target"]).to_pandas()
    cm = cm.sort_values("cell_integer_id").reset_index(drop=True)
    assert (cm["cell_integer_id"].values == np.arange(len(cm))).all()
    tot = cm["total_counts"].values.astype(np.float64)
    gt_names, gt_code = np.unique(cm["gene_target"].values, return_inverse=True)
    gt_code = gt_code.astype(np.int32); G = gt_names.size

    # eligible target genes: >= MIN_N cells (exclude NTC), from the atlas cell table itself
    per_gene_n = np.bincount(gt_code, minlength=G)
    ntc_idx = int(np.where(gt_names == NTC)[0][0])
    is_target = (per_gene_n >= MIN_N)
    is_target[ntc_idx] = False
    tgt_codes = np.where(is_target)[0]
    print(f"[p4:{a.line}] {len(cm):,} cells, {G:,} genes; {tgt_codes.size} target genes "
          f"(N>={MIN_N}); retaining ~{int(per_gene_n[is_target].sum()):,} cells", flush=True)

    # NTC centroid (base) from committed pseudobulk sufficient stats
    z = np.load(a.pseudobulk, allow_pickle=True)
    pnames = z["gt_names"]; pcnt = z["cnt"].astype(np.float64); pcsum = z["csum"].astype(np.float64)
    pn_ntc = int(np.where(pnames == NTC)[0][0])
    base_ntc = pcsum[pn_ntc] / pcnt[pn_ntc]

    ex = lance.dataset(f"{HF}/data/{a.line}/expression.lance")
    total_rows = ex.count_rows()

    if os.path.exists(ckpt):
        zz = np.load(ckpt, allow_pickle=True)
        hcoords = zz["hcoords"].astype(np.float64); hcode = zz["hcode"].astype(np.int32)
        offset = int(zz["offset_done"])
        buf_c = [hcoords]; buf_g = [hcode]
        print(f"[p4:{a.line}] RESUME offset {offset:,}/{total_rows:,} "
              f"({100*offset/total_rows:.1f}%), retained {hcode.size:,} cells", flush=True)
    else:
        buf_c, buf_g, offset = [], [], 0

    def save(off):
        allc = np.concatenate(buf_c) if buf_c else np.zeros((0, d))
        allg = np.concatenate(buf_g) if buf_g else np.zeros(0, np.int32)
        np.savez(ckpt, hcoords=allc, hcode=allg, offset_done=off)
        # collapse buffers to keep the list short
        buf_c[:] = [allc]; buf_g[:] = [allg]

    def read_chunk(off):
        last = None
        for att in range(6):
            try:
                return ex.scanner(columns=["cell_integer_id", "gene_integer_id", "value"],
                                  limit=a.chunk, offset=off).to_table()
            except Exception as e:
                last = e; print(f"[p4:{a.line}] retry {att+1} @ {off:,}", flush=True)
                time.sleep(min(60, 5 * 2 ** att))
        raise last

    t0 = time.time(); ch = 0
    while offset < total_rows:
        tab = read_chunk(offset)
        n = tab.num_rows
        if n == 0:
            break
        cid = tab.column("cell_integer_id").to_numpy()
        last = cid[-1]
        end_of_table = (n < a.chunk) or (offset + n >= total_rows)
        if end_of_table:
            complete = np.ones(n, dtype=bool); consumed = n
        else:
            complete = cid < last; consumed = int(complete.sum())
        gid = tab.column("gene_integer_id").to_numpy()
        val = tab.column("value").to_numpy().astype(np.float64)
        cid = cid[complete]; gid = gid[complete]; val = val[complete]
        keep = hvg_col[gid] >= 0
        cid = cid[keep]; col = hvg_col[gid[keep]]; val = val[keep]
        uc, inv = np.unique(cid, return_inverse=True)
        ln = np.log1p(val / tot[cid] * 1e4)
        S = sp.csr_matrix((ln, (inv, col)), shape=(uc.size, hvg.size))
        coords = S @ comps.T - mean_proj                 # (uc,50)
        gc = gt_code[uc]
        m_t = is_target[gc]
        if m_t.any():
            buf_c.append(coords[m_t]); buf_g.append(gc[m_t].astype(np.int32))
        offset += consumed; ch += 1
        if ch % a.ckpt_every == 0 or end_of_table:
            save(offset)
            el = time.time() - t0; rate = offset / el if el > 0 else 0
            eta = (total_rows - offset) / rate / 3600 if rate > 0 else 0
            nret = sum(len(x) for x in buf_g)
            print(f"[p4:{a.line}] {100*offset/total_rows:5.1f}%  retained {nret:,} cells  "
                  f"ETA {eta:.1f}h", flush=True)

    # ---- held-out per target gene ----
    allc = np.concatenate(buf_c); allg = np.concatenate(buf_g)
    print(f"[p4:{a.line}] stream done: {allg.size:,} retained cells over {tgt_codes.size} genes; "
          f"running held-out ({REPS} reps)...", flush=True)
    recs = []
    for k, code in enumerate(tgt_codes):
        C = allc[allg == code]
        if len(C) < MIN_N:
            continue
        r = held_out(C, base_ntc, d)
        if r is None:
            continue
        r["gene_target"] = str(gt_names[code]); recs.append(r)
        if (k + 1) % 200 == 0:
            print(f"[p4:{a.line}] held-out {k+1}/{tgt_codes.size}", flush=True)

    ratio = np.array([r["ratio"] for r in recs]); rho2 = np.array([r["rho2"] for r in recs])
    r2s = np.array([r["r2"] for r in recs]); inreg = rho2 >= 3
    deficit = 1.0 - ratio
    from scipy.stats import spearmanr
    sp_def = float(spearmanr(1.0 / np.maximum(rho2, 1e-9), deficit).correlation) if len(recs) > 2 else None
    summary = dict(
        line=a.line, min_cells=MIN_N, reps=REPS, n_groups=len(recs),
        median_ratio_all=round(float(np.median(ratio)), 4),
        median_r2_all=round(float(np.median(r2s)), 5),
        n_in_regime=int(inreg.sum()),
        median_ratio_in_regime=round(float(np.median(ratio[inreg])), 4) if inreg.any() else None,
        median_r2_in_regime=round(float(np.median(r2s[inreg])), 5) if inreg.any() else None,
        spearman_deficit_vs_inv_rho2=round(sp_def, 4) if sp_def is not None else None,
        m_range=[round(float(min(r["m"] for r in recs)), 3), round(float(max(r["m"] for r in recs)), 3)],
    )
    out = dict(description=f"Orion {a.line} held-out angular-error falsification, all N>={MIN_N} "
               f"knockdowns vs pooled-NTC centroid (streamed per-cell). Same routine as EmeraldBay pass4.",
               **summary, per_group=recs)
    op = os.path.join(a.fixtures, f"orion_{a.line}_heldout.json")
    json.dump(out, open(op, "w"), indent=1)
    print(f"[p4:{a.line}] WROTE {op}\n  {json.dumps(summary)}", flush=True)


if __name__ == "__main__":
    main()
