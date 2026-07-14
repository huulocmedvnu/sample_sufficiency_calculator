"""
Genetic data pipeline (X-Atlas/Orion, TRADE). Emits the STANDARD per-condition structure the unified
engine consumes; it does NOT compute m, tr(S), bias-correction, detectable or n* (that is engine.py).

condition = (gene knockdown), control = the shared non-targeting (NTC) pool, n_c = NTC cell count.
Reads the committed per-condition sufficient statistics / covariance (no atlas re-stream):
  Orion:  <orion_work>/full/<line>/pseudobulk.npz   (gt_names, cnt, csum, csq; within Sigma diagonal)
  TRADE:  <trade_work>/out/<line>/basis.npz          (per-cell coords + gene labels + full 50x50 Sigma)

Scored-condition QC: a knockdown must have >= MIN_CELLS acquired cells to have an estimable centroid
(this is a measurement requirement, not the arbitrary cell-count denominator filter of finding #2).
"""
import os, re, numpy as np

ORION = "/mnt/hdd2/loc-tran/orion_work/full"
TRADE = "/mnt/hdd2/loc-tran/trade_work/out"
MIN_CELLS = 25
NTC_RE = re.compile(r"non.?target|^ntc$|control|safe.?harbor|scramble|negative", re.I)


def _std(cond_id, mu_t, mu_c, n_t, n_c, Sigma, name, control_type):
    return dict(dataset=name, modality="genetic", control_type=control_type,
                cond_id=np.asarray(cond_id), mu_t=np.asarray(mu_t, float), mu_c=np.asarray(mu_c, float),
                n_t=np.asarray(n_t, float), n_c=np.asarray(n_c, float), Sigma=np.asarray(Sigma, float))


def load_orion(line):
    z = np.load(os.path.join(ORION, line, "pseudobulk.npz"), allow_pickle=True)
    names = z["gt_names"].astype(str); cnt = z["cnt"].astype(np.float64)
    csum = z["csum"].astype(np.float64); csq = z["csq"].astype(np.float64)
    D = int(z["n_comps"])
    ok = cnt >= 2                                              # within-Sigma over conditions with >=2 cells
    ss = (csq[ok] - csum[ok] ** 2 / cnt[ok, None]).sum(0); dof = (cnt[ok] - 1).sum()
    ell = ss / dof; Sigma = np.diag(ell)                      # Orion within-Sigma is diagonal (per-PC)
    ntc = [i for i, g in enumerate(names) if NTC_RE.search(g)]
    ni = max(ntc, key=lambda i: cnt[i]); mu_ntc = csum[ni] / cnt[ni]; n_ntc = cnt[ni]
    sel = [i for i, g in enumerate(names) if i != ni and cnt[i] >= MIN_CELLS]
    mu_t = np.array([csum[i] / cnt[i] for i in sel]); n_t = np.array([cnt[i] for i in sel])
    return _std([names[i] for i in sel], mu_t, np.tile(mu_ntc, (len(sel), 1)), n_t,
                np.full(len(sel), n_ntc), Sigma, f"Orion {line}", "shared NTC pool")


def load_trade(line):
    z = np.load(os.path.join(TRADE, line, "basis.npz"), allow_pickle=True)
    coords = z["coords"].astype(np.float64); gene = z["gene"].astype(str)
    Sigma = z["Sigma"].astype(np.float64); ntc_label = str(z["ntc_label"])
    genes, idx = np.unique(gene, return_inverse=True)
    csum = np.zeros((len(genes), coords.shape[1])); cnt = np.zeros(len(genes))
    np.add.at(csum, idx, coords); np.add.at(cnt, idx, 1.0)
    ni = int(np.where(genes == ntc_label)[0][0]); mu_ntc = csum[ni] / cnt[ni]; n_ntc = cnt[ni]
    sel = [i for i in range(len(genes)) if i != ni and cnt[i] >= MIN_CELLS]
    mu_t = np.array([csum[i] / cnt[i] for i in sel]); n_t = np.array([cnt[i] for i in sel])
    disp = "TRADE " + ("Jurkat" if "jurkat" in line else "HepG2")
    return _std([genes[i] for i in sel], mu_t, np.tile(mu_ntc, (len(sel), 1)), n_t,
                np.full(len(sel), n_ntc), Sigma, disp, "shared NTC pool")


DATASETS = {"orion_HCT116": lambda: load_orion("HCT116"), "orion_HEK293T": lambda: load_orion("HEK293T"),
            "trade_jurkat": lambda: load_trade("jurkat"), "trade_hepg2": lambda: load_trade("hepg2")}


def load(key):
    return DATASETS[key]()
