"""
Chemical data pipeline (Tahoe-100M, EmeraldBay). Emits the STANDARD per-condition structure the
unified engine consumes; it does NOT compute m, tr(S), bias-correction, detectable or n* (engine.py).

condition = (drug x dose x line). Reads committed per-condition sufficient statistics (no atlas
re-stream):
  Tahoe:       <tahoe_work>/out_dose/{pseudobulk,basis}.npz + sample metadata  (raw-feature sums;
               projected through the committed PCA basis; pooled over plate-replicates)
  EmeraldBay:  <eb_work>/out/pseudobulk.npz                                     (already in PCA space)

Control design (finding #9, gap 2): the Tahoe HEADLINE control is its real vehicle-matched shared
DMSO_TF pool (all 2-3 vehicle wells on the plate pooled; n_c ~ 3,100 cells, shared across ~94
conditions), which is control-pool-limited for most conditions. The per-line mean (a large pool, n_c = all cells in the line) is reported as a sensitivity
analysis. EmeraldBay references the per-line mean.

No cell-count filter on the spectrum (finding #2/#4): every non-empty condition is included.
Within-condition Sigma is the REAL diagonal covariance diag(ell_within) computed from per-condition
sufficient statistics (fixtures/chemical_within_cov.json; gap 1/3), not an isotropic scalar. The
per-PC variance spans 3.6-5x but tr(P.Sigma.P) varies only +-1.5-3% across conditions.
"""
import os, ast, json, numpy as np, pyarrow.parquet as pq

TAHOE = "/mnt/hdd2/loc-tran/tahoe_work/out_dose"
TAHOE_META = "/mnt/hdd2/loc-tran/tahoe_work/meta/metadata"
EB = "/mnt/hdd2/loc-tran/eb_work/out"
CTRL = "DMSO_TF"

# REAL within-condition covariance (diagonal), computed from per-condition sufficient statistics
# (fixtures/chemical_within_cov.json). No hardcoded scalar sigma2; the engine derives sigma2 = tr(Sigma)/d.
_COV = json.load(open(os.path.join(os.path.dirname(__file__), "..", "..", "fixtures", "chemical_within_cov.json")))
ELL = {"Tahoe-100M": np.array(_COV["tahoe"]["ell"]), "EmeraldBay": np.array(_COV["emeraldbay"]["ell"])}


def _std(name, cond_id, mu_t, mu_c, n_t, n_c, ell, control_type):
    return dict(dataset=name, modality="chemical", control_type=control_type,
                cond_id=np.asarray(cond_id, dtype=object), mu_t=np.asarray(mu_t, float),
                mu_c=np.asarray(mu_c, float), n_t=np.asarray(n_t, float), n_c=np.asarray(n_c, float),
                Sigma=np.diag(np.asarray(ell, float)))          # REAL diagonal within-Sigma (anisotropic)


def load_tahoe(control="per-line-mean"):
    b = np.load(f"{TAHOE}/basis.npz"); comps = b["components"].astype(np.float64); pmean = b["pca_mean"].astype(np.float64)
    pb = np.load(f"{TAHOE}/pseudobulk.npz", allow_pickle=True); keys = [str(k) for k in pb["cond_keys"]]
    sums = pb["sums"].astype(np.float64); counts = pb["counts"].astype(np.float64)
    coords = (sums / counts[:, None] - pmean) @ comps.T
    sm = pq.read_table(f"{TAHOE_META}/sample_metadata.parquet").to_pandas()
    sm["dose"] = sm.drugname_drugconc.map(lambda s: (lambda t: float(t[0][1]))(ast.literal_eval(s)) if s else np.nan)
    S2 = {r.sample: (r.drug, r.dose, r.plate) for r in sm.itertuples()}
    lsum = {}; lcnt = {}; dmso_c = {}; dmso_n = {}; cond = {}
    for i, k in enumerate(keys):
        s, line = k.rsplit("|", 1); mt = S2.get(s)
        if mt is None: continue
        drug, ds, plate = mt
        lsum[line] = lsum.get(line, 0) + sums[i]; lcnt[line] = lcnt.get(line, 0) + counts[i]
        if drug == CTRL:
            # POOL every DMSO_TF well on this (plate, line): Tahoe plates carry 2-3 vehicle wells each.
            # (A plain assignment here silently kept only the last well -- halving n_c; audit 2026-09.)
            k2 = (plate, line); n_prev = dmso_n.get(k2, 0.0)
            dmso_c[k2] = (dmso_c.get(k2, 0.0) * n_prev + coords[i] * counts[i]) / (n_prev + counts[i])
            dmso_n[k2] = n_prev + counts[i]
        else:
            cond.setdefault((drug, ds, line), []).append((coords[i], counts[i], plate))
    linebase = {l: (lsum[l] / lcnt[l] - pmean) @ comps.T for l in lsum}
    cid, mu_t, mu_c, n_t, n_c = [], [], [], [], []
    for (drug, ds, line), items in cond.items():
        tot = sum(c for _, c, _ in items)
        if tot <= 0: continue
        cen = sum(c * co for co, c, _ in items) / tot
        if control == "per-line-mean":
            mc = linebase[line]; nc = lcnt[line]
        elif control == "shared-dmso":
            if not all((p, line) in dmso_c for _, _, p in items): continue
            mc = sum(c * dmso_c[(p, line)] for _, c, p in items) / tot                 # treated-weighted matched DMSO
            nc = tot * tot / sum(c * c / dmso_n[(p, line)] for _, c, p in items)        # effective shared-pool size
        else:
            raise ValueError(control)
        cid.append(f"{drug}|{ds}|{line}"); mu_t.append(cen); mu_c.append(mc); n_t.append(tot); n_c.append(nc)
    ct = "per-line mean (large pool)" if control == "per-line-mean" else "shared vehicle DMSO_TF pool"
    return _std("Tahoe-100M", cid, np.array(mu_t), np.array(mu_c), np.array(n_t), np.array(n_c),
                ELL["Tahoe-100M"], ct)


def load_emeraldbay():
    z = np.load(f"{EB}/pseudobulk.npz", allow_pickle=True)
    ck = np.array([str(k) for k in z["cond_keys"]]); sums = z["sums"].astype(np.float64); cnt = z["counts"].astype(np.float64)
    mu = sums / cnt[:, None]; line = np.array([k.split("|")[1] for k in ck])          # sums already in PCA space
    lsum = {}; lcnt = {}
    for i, l in enumerate(line): lsum[l] = lsum.get(l, 0) + sums[i]; lcnt[l] = lcnt.get(l, 0) + cnt[i]
    base = np.array([lsum[l] / lcnt[l] for l in line]); nc = np.array([lcnt[l] for l in line])
    return _std("EmeraldBay", list(ck), mu, base, cnt, nc, ELL["EmeraldBay"], "per-line mean (large pool)")


DATASETS = {"Tahoe-100M": load_tahoe, "EmeraldBay": load_emeraldbay}


def load(key, **kw):
    return DATASETS[key](**kw)
