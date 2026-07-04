"""
Pass 0: fetch the TRADE (Nadig/Replogle/Weissman, Nat Genet 2025) essential-gene Perturb-seq h5ad files
from GEO GSE264667. Resumable (curl -C -). These are raw single-cell UMI counts with obs metadata
(UMI_count, mitopercent, z_gemgroup_UMI, gem_group, gene [target], sgID_AB [dual-guide identity]).

  GSE264667_jurkat_raw_singlecell_01.h5ad  (~8.7 GB, Jurkat-Essential)
  GSE264667_hepg2_raw_singlecell_01.h5ad   (~5.2 GB, HepG2-Essential)

Usage:  python pass0_download.py [--dir /mnt/hdd2/loc-tran/trade_work]
(uses the system `curl`; the two files are small enough to keep locally.)
"""
import os, subprocess, argparse

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl"
FILES = {
    "jurkat": "GSE264667_jurkat_raw_singlecell_01.h5ad",
    "hepg2":  "GSE264667_hepg2_raw_singlecell_01.h5ad",
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="/mnt/hdd2/loc-tran/trade_work")
    args = ap.parse_args()
    os.makedirs(args.dir, exist_ok=True)
    for line, fn in FILES.items():
        dst = os.path.join(args.dir, fn)
        print(f"[pass0] {line}: {fn}", flush=True)
        subprocess.run(["curl", "-sS", "-C", "-", "-o", dst, f"{BASE}/{fn}"], check=True)
        print(f"[pass0]   -> {os.path.getsize(dst)/1e9:.2f} GB", flush=True)
    print("[pass0] done")

if __name__ == "__main__":
    main()
