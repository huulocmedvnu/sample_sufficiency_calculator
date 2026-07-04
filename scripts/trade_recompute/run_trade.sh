#!/usr/bin/env bash
# TRADE essential-gene sample-sufficiency extension, both lines, end to end.
#   pass0 download (GSE264667) -> pass1 QC+global-embedding+within-Sigma+quota -> pass2 falsification.
# Idempotent; pass1/2 are fast (in-memory, ~10^5 cells) so no checkpointing needed.
set -u
cd "$(dirname "$0")/../.." || exit 1
PY=python
WORK=/mnt/hdd2/loc-tran/trade_work
declare -A H5=( [hepg2]="$WORK/GSE264667_hepg2_raw_singlecell_01.h5ad"
                [jurkat]="$WORK/GSE264667_jurkat_raw_singlecell_01.h5ad" )
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

# pass0: download if missing
for line in hepg2 jurkat; do
  [ -f "${H5[$line]}" ] || { log "downloading $line..."; $PY scripts/trade_recompute/pass0_download.py --dir "$WORK" || exit 2; }
done

# process smaller line (hepg2) first so results land sooner
for line in hepg2 jurkat; do
  log "$line: QC + global embedding + within-Sigma + quota..."
  $PY scripts/trade_recompute/pass1_qc_basis_quota.py --line "$line" --h5ad "${H5[$line]}" --out "$WORK/out" --theta 0.1 || exit 3
  log "$line: downsample-and-measure falsification (Phase C)..."
  $PY scripts/trade_recompute/pass2_falsification.py --line "$line" --dir "$WORK/out" || exit 4
done
log "ALL DONE. fixtures/trade_{hepg2,jurkat}_falsification.json ; $WORK/out/<line>/{quota.csv,summary.json}"
