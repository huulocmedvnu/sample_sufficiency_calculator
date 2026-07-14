"""Pass 7 -- real within-condition diagonal covariance for the chemical atlases (gap 1/3).

Streams the cached per-cell PCA-coord memmap from pass4a_stream.py (no raw re-stream), accumulates
per-condition csum/csq/cnt, and writes the pooled within-condition per-PC variance vector ell_within
(diagonal Sigma) + sigma2 = mean(ell) to fixtures/chemical_within_cov.json, plus the per-condition
sufficient statistics to fixtures/tahoe_within_sumsq.npz. This replaces the hardcoded scalar sigma2:
the engine now sees the REAL anisotropic diagonal Sigma. Tahoe sigma2 reproduces the committed 0.9567.
EmeraldBay's ell is computed directly from its pseudobulk sumsq (already stored).
"""
import os, json, numpy as np
BASE = "/mnt/hdd2/loc-tran/tahoe_work/pass4"; NC = 50
FX = os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")


def tahoe_ell():
    z = np.load(f"{BASE}/ckpt.npz", allow_pickle=True); cursor = int(z["cursor"]); K = len(z["cond_keys"])
    coords = np.memmap(f"{BASE}/coords.f32", dtype=np.float32, mode="r", shape=(cursor, NC))
    condid = np.memmap(f"{BASE}/condid.i32", dtype=np.int32, mode="r", shape=(cursor,))
    csum = np.zeros((K, NC)); csq = np.zeros((K, NC)); cnt = np.zeros(K)
    for s in range(0, cursor, 8_000_000):
        e = min(s + 8_000_000, cursor); cc = np.asarray(condid[s:e]); X = np.asarray(coords[s:e], np.float64)
        cnt += np.bincount(cc, minlength=K)
        for k in range(NC):
            csum[:, k] += np.bincount(cc, weights=X[:, k], minlength=K)
            csq[:, k] += np.bincount(cc, weights=X[:, k] * X[:, k], minlength=K)
    ok = cnt >= 2; ell = (csq[ok] - csum[ok] ** 2 / cnt[ok, None]).sum(0) / (cnt[ok] - 1).sum()
    np.savez_compressed(f"{FX}/tahoe_within_sumsq.npz", cond_keys=z["cond_keys"],
                        csum=csum.astype(np.float32), csq=csq.astype(np.float32), cnt=cnt.astype(np.int32))
    return ell


def eb_ell():
    z = np.load("/mnt/hdd2/loc-tran/eb_work/out/pseudobulk.npz", allow_pickle=True)
    s, sq, c = z["sums"], z["sumsq"], z["counts"]; ok = c >= 2
    return (sq[ok] - s[ok] ** 2 / c[ok, None]).sum(0) / (c[ok] - 1).sum()


if __name__ == "__main__":
    te = tahoe_ell(); ee = eb_ell()
    out = dict(note="Within-condition per-PC variance (diagonal Sigma) for the chemical atlases, from "
                    "per-condition sufficient statistics. sigma2 = mean(ell). Replaces the hardcoded scalar.",
               tahoe=dict(sigma2=round(float(te.mean()), 4), ell=[round(float(x), 6) for x in te],
                          per_pc_min=round(float(te.min()), 4), per_pc_max=round(float(te.max()), 4)),
               emeraldbay=dict(sigma2=round(float(ee.mean()), 4), ell=[round(float(x), 6) for x in ee],
                               per_pc_min=round(float(ee.min()), 4), per_pc_max=round(float(ee.max()), 4)))
    json.dump(out, open(f"{FX}/chemical_within_cov.json", "w"), indent=1)
    print(f"Tahoe sigma2={out['tahoe']['sigma2']} (committed 0.9567)  EmeraldBay sigma2={out['emeraldbay']['sigma2']}")
