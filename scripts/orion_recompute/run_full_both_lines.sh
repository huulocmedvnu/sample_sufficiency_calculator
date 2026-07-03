#!/usr/bin/env bash
# Full-scale genome-wide quota spectrum for BOTH X-Atlas/Orion lines, fully automated.
#   1. fit a head-block HVG/PCA(50) basis per line (reuse HCT116 pilot basis; fit HEK293T fresh)
#   2. stream the ENTIRE expression table per line -> per-gene pseudobulk in PCA space (resumable)
#   3. compute the genome-wide n* spectrum per line -> fixtures/orion_<LINE>_{quota.csv,summary.json}
# Idempotent: re-running skips completed bases and resumes pass2 from its checkpoint.
set -u
cd "$(dirname "$0")/../.." || exit 1          # repo root
PY=python
FULL=/mnt/hdd2/loc-tran/orion_work/full
PILOT=/mnt/hdd2/loc-tran/orion_work/pilot
mkdir -p "$FULL/HCT116" "$FULL/HEK293T"
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

# --- HCT116 basis: reuse the validated pilot basis ---
if [ ! -f "$FULL/HCT116/basis.npz" ]; then
  cp "$PILOT/HCT116/basis.npz" "$FULL/HCT116/basis.npz" && log "HCT116 basis: reused pilot basis"
fi

# --- HEK293T basis: download a representative head-block and fit ---
if [ ! -f "$FULL/HEK293T/basis.npz" ]; then
  if [ ! -f "$FULL/HEK293T/expr_coo.npz" ]; then
    log "HEK293T: downloading 250k head-block for basis..."
    OUT="$FULL" $PY scripts/orion_recompute/pass0_download_head.py --line HEK293T --max-cells 250000 || exit 2
  fi
  log "HEK293T: fitting HVG/PCA(50) basis + sigma^2..."
  $PY scripts/orion_recompute/pass1_basis_sigma.py --line HEK293T --dir "$FULL" || exit 3
fi

# --- Full streaming pseudobulk (resumable) + quota, per line ---
for LINE in HCT116 HEK293T; do
  log "$LINE: streaming full expression table (resumable)..."
  $PY scripts/orion_recompute/pass2_stream_pseudobulk.py \
        --line "$LINE" --dir "$FULL" --basis "$FULL/$LINE/basis.npz" || { log "$LINE pass2 FAILED"; exit 4; }
  log "$LINE: computing genome-wide quota spectrum..."
  $PY scripts/orion_recompute/pass3_quota_full.py --line "$LINE" --dir "$FULL" --theta 0.1 || exit 5
done

log "ALL DONE. Fixtures: fixtures/orion_HCT116_*.{csv,json}, fixtures/orion_HEK293T_*.{csv,json}"
