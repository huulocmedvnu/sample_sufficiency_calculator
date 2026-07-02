#!/usr/bin/env python
"""Stage B of the full Tahoe direct-curve falsification: subsample-and-measure angular-error curves
per condition, read from the memmap written by pass4a_stream.py.

For each (sample x line) condition with >= MIN_CELLS cells: estimate the direction from all cells
(vs the per-line-mean baseline), then for a grid of n subsample (REPS reps) and measure realized RMS
angle; fit the through-origin slope of theta^2 vs (1/n - 1/N) and compare to the parameter-free
prediction tr(PSP)/m^2. Resumable (checkpoints results every 2000 conditions).

Writes fixtures/tahoe_direct_curves_full.json and outputs/tahoe_direct_curves_full.png.
Env: MIN_CELLS(400), REPS(120), NGRID(7), NMAX(5000).
"""
import os, json, numpy as np
BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"
_dose = "/mnt/hdd2/loc-tran/tahoe_work/out_dose"
BASIS = _dose if os.path.exists(os.path.join(_dose, "basis.npz")) else "/mnt/hdd2/loc-tran/tahoe_work/out"
HERE = os.path.dirname(__file__)
FIX = os.path.join(HERE, "..", "..", "fixtures", "tahoe_direct_curves_full.json")
FIGDIR = os.path.join(HERE, "..", "..", "outputs"); os.makedirs(FIGDIR, exist_ok=True)
PROG = os.path.join(BASE, "curves_progress.json")
MIN_CELLS = int(os.environ.get("MIN_CELLS", "400")); REPS = int(os.environ.get("REPS", "120"))
NGRID = int(os.environ.get("NGRID", "7")); NMAX = int(os.environ.get("NMAX", "5000"))
rng = np.random.default_rng(0)

NC = int(np.load(os.path.join(BASIS, "basis.npz"))["n_comps"])
z = np.load(os.path.join(BASE, "ckpt.npz"), allow_pickle=True)
cursor = int(z["cursor"]); cond_keys = list(z["cond_keys"]); cond_drug = dict(z["cond_drug"].item())
line_sum = dict(z["line_sum"].item()); line_cnt = dict(z["line_cnt"].item())
baseline = {l: line_sum[l] / line_cnt[l] for l in line_sum}
print(f"[curves] {cursor:,} cells, {len(cond_keys)} conditions", flush=True)

coords = np.memmap(os.path.join(BASE, "coords.f32"), dtype=np.float32, mode="r", shape=(int(os.path.getsize(os.path.join(BASE,'coords.f32'))//(4*NC)), NC))[:cursor]
condid = np.memmap(os.path.join(BASE, "condid.i32"), dtype=np.int32, mode="r", shape=(cursor,))

# --- sigma^2 full-atlas confirmation: Tahoe headline 2.406 is the MARGINAL per-cell variance ---
# (mean over PCs of the per-PC variance across ALL projected cells). Chunked to stay memory-safe.
_bf = np.load(os.path.join(BASIS, "basis.npz"))
sigma2_fit = float(_bf["per_pc_var"].mean()) if "per_pc_var" in _bf else float(_bf["sigma2"])
_psum = np.zeros(NC); _pssq = np.zeros(NC); _CH = 2_000_000
for _s in range(0, cursor, _CH):
    _blk = np.asarray(coords[_s:_s + _CH], dtype=np.float64)
    _psum += _blk.sum(0); _pssq += (_blk * _blk).sum(0)
_mean = _psum / cursor; _var = _pssq / cursor - _mean ** 2
sigma2_marginal_full = float(_var.mean())
sigma2_confirmation = dict(
    scope=f"{cursor} cells (full atlas)", marginal_full_atlas=round(sigma2_marginal_full, 4),
    fit_subsample_estimate=round(sigma2_fit, 4), manuscript_marginal=2.406,
    matches_manuscript=bool(abs(sigma2_marginal_full - 2.406) < 0.05))
print(f"[sigma2] Tahoe full-atlas marginal={sigma2_marginal_full:.4f} (manuscript 2.406; "
      f"14-shard fit estimate {sigma2_fit:.4f}) over {cursor:,} cells", flush=True)

# group rows by condition id
order = np.argsort(condid, kind="stable")
cid_sorted = condid[order]
bnd = np.searchsorted(cid_sorted, np.arange(len(cond_keys) + 1))

records = []
if os.path.exists(PROG):
    records = json.load(open(PROG))["records"]; done_ids = {r["cid"] for r in records}
    print(f"[resume] {len(records)} conditions already done", flush=True)
else:
    done_ids = set()


def curve(C, base):
    N = len(C); v = C.mean(0) - base; m = float(np.linalg.norm(v))
    if m <= 0: return None
    u = v / m; P = np.eye(NC) - np.outer(u, u); Sig = np.cov(C.T.astype(np.float64))
    trPSP = float(np.trace(P @ Sig @ P)); uSu = float(u @ Sig @ u)
    if trPSP <= 0 or uSu <= 0: return None
    rho2 = m ** 2 / uSu
    ns = np.unique(np.geomspace(20, min(N // 2, NMAX), NGRID).round().astype(int)); ns = ns[ns < N]
    xs, ys = [], []
    for n in ns:
        # WITHOUT replacement from the acquired pool of N (matches EmeraldBay pass3/pass4 and the
        # finite-population prediction theta^2 = tr(PSP)/m^2 * (1/n - 1/N)); vectorised over REPS.
        idx = np.argsort(rng.random((REPS, N)), axis=1)[:, :int(n)]
        cent = C[idx].mean(1) - base                            # (REPS, 50)
        un = cent / np.linalg.norm(cent, axis=1, keepdims=True)
        ang2 = np.arccos(np.clip(un @ u, -1, 1)) ** 2
        ys.append(float(ang2.mean())); xs.append(1.0 / n - 1.0 / N)
    xs, ys = np.array(xs), np.array(ys)
    fit = float((xs @ ys) / (xs @ xs))
    sst = float(((ys - ys.mean()) ** 2).sum())
    r2 = 1.0 - float(((ys - fit * xs) ** 2).sum()) / sst if sst > 0 else 1.0
    return dict(N=int(N), m=round(m, 4), rho2=round(float(rho2), 3), pred_slope=round(trPSP / m ** 2, 6),
                fit_slope=round(fit, 6), ratio=round(fit / (trPSP / m ** 2), 4), r2=round(r2, 5))


t_last = 0
for j in range(len(cond_keys)):
    if j in done_ids: continue
    a, b_ = bnd[j], bnd[j + 1]
    if b_ - a < MIN_CELLS: continue
    key = cond_keys[j]; line = key.rsplit("|", 1)[1]
    if line not in baseline: continue
    C = np.asarray(coords[order[a:b_]], dtype=np.float64)
    r = curve(C, baseline[line])
    if r is None: continue
    r.update(cid=j, sample=key.rsplit("|", 1)[0], line=line, drug=cond_drug.get(key, "?"))
    records.append(r)
    if len(records) % 2000 == 0:
        json.dump({"records": records}, open(PROG, "w"))
        print(f"[curves] {len(records)} conditions done (scan {j}/{len(cond_keys)})", flush=True)

json.dump({"records": records}, open(PROG, "w"))

pred = np.array([r["pred_slope"] for r in records]); fit = np.array([r["fit_slope"] for r in records])
rho2 = np.array([r["rho2"] for r in records]); ratio = fit / pred; ms = np.array([r["m"] for r in records])
inreg = rho2 >= 3.0
summary = dict(
    description="Tahoe-100M FULL direct subsample-and-measure angular-error test (all shards re-streamed; "
                "same PCA(50) basis; per-line-mean baseline; condition = sample x line). Parameter-free "
                "slope tr(PSP)/m^2 vs realized, across every condition with >=%d cells." % MIN_CELLS,
    n_conditions=len(records), min_cells=MIN_CELLS, reps=REPS,
    sigma2_confirmation=sigma2_confirmation,
    pearson_loglog=round(float(np.corrcoef(np.log10(pred), np.log10(fit))[0, 1]), 4),
    median_ratio_all=round(float(np.median(ratio)), 4),
    n_in_regime=int(inreg.sum()),
    median_ratio_in_regime=round(float(np.median(ratio[inreg])), 4) if inreg.any() else None,
    median_r2_in_regime=round(float(np.median([records[i]["r2"] for i in range(len(records)) if inreg[i]])), 5) if inreg.any() else None,
    m_range=[round(float(ms.min()), 3), round(float(ms.max()), 3)],
    spearman_deficit_vs_inv_rho2=round(float(np.corrcoef(
        np.argsort(np.argsort(1.0 / rho2)), np.argsort(np.argsort(1.0 - ratio)))[0, 1]), 4),
    per_group=records)
json.dump(summary, open(FIX, "w"), indent=1)
print(f"[curves] DONE {len(records)} conds; Pearson(loglog)={summary['pearson_loglog']}; "
      f"in-regime n={summary['n_in_regime']} ratio={summary['median_ratio_in_regime']} "
      f"R2={summary['median_r2_in_regime']}; wrote {FIX}", flush=True)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    lo, hi = min(pred.min(), fit.min()), max(pred.max(), fit.max())
    ax[0].loglog([lo, hi], [lo, hi], "k--", lw=1, label="y = x")
    ax[0].scatter(pred, fit, c=np.log10(rho2), cmap="viridis", s=6, alpha=0.4)
    ax[0].set_xlabel("predicted slope tr(PΣP)/m²"); ax[0].set_ylabel("realized slope")
    ax[0].set_title(f"Tahoe FULL — {len(records):,} conditions"); ax[0].legend(fontsize=8)
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].scatter(rho2, ratio, c=np.log10(ms), cmap="plasma", s=6, alpha=0.4)
    ax[1].set_xscale("log"); ax[1].set_xlabel("ρ² = m²/(uᵀΣu)"); ax[1].set_ylabel("realized/predicted slope")
    ax[1].set_title("slope ratio vs SNR")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "tahoe_direct_curves_full.png"), dpi=130)
    print("[fig] wrote tahoe_direct_curves_full.png", flush=True)
except Exception as e:
    print(f"[fig] skipped ({e})", flush=True)
