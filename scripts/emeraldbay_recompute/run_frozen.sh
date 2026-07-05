#!/bin/bash
# Full EmeraldBay recompute with a FROZEN HVG set, written to an ISOLATED output dir + separate fixtures
# so the committed per-fit-HVG results are preserved until the frozen numbers are reviewed and promoted.
#
# Prereq: fixtures/emeraldbay_frozen_hvg.json must exist (make_frozen_hvg_fixture.py).
# Env isolation:
#   OUT_EB        -> eb_work/out_frozen           (basis, pseudobulk, shared_cells, all_cells)
#   FROZEN_HVG    -> the frozen HVG token-id list (pass1 skips per-fit seurat selection)
#   EB_CALIB_FIX  -> fixtures/emeraldbay_calibration_frozen.json   (pass3 + pass5)
#   EB_FIX        -> fixtures/emeraldbay_falsification_frozen.json (pass4)
set -euo pipefail
cd /mnt/hdd2/loc-tran/sample_sufficiency_calculator
ROOT=/mnt/hdd2/loc-tran/sample_sufficiency_calculator
export OUT_EB=/mnt/hdd2/loc-tran/eb_work/out_frozen
export FROZEN_HVG=$ROOT/fixtures/emeraldbay_frozen_hvg.json
export EB_CALIB_FIX=$ROOT/fixtures/emeraldbay_calibration_frozen.json
mkdir -p "$OUT_EB"

if [ ! -f "$FROZEN_HVG" ]; then echo "FATAL: $FROZEN_HVG missing (run make_frozen_hvg_fixture.py)"; exit 1; fi

echo "=== [1/6] pass1_basis (frozen HVG) -> $OUT_EB/basis.npz ==="
python scripts/emeraldbay_recompute/pass1_basis.py

echo "=== [2/6] pass2_project (pseudobulk + shared_cells) ==="
python scripts/emeraldbay_recompute/pass2_project.py

echo "=== [3/6] pass2b_project_all (52 lines -> all_cells.npz) ==="
python scripts/emeraldbay_recompute/pass2b_project_all.py

echo "=== [4/6] pass3_validate -> $EB_CALIB_FIX ==="
python scripts/emeraldbay_recompute/pass3_validate.py

echo "=== [5/6] pass4_falsification (full atlas) -> falsification_frozen.json ==="
EB_COORDS=$OUT_EB/all_cells.npz \
EB_FIX=$ROOT/fixtures/emeraldbay_falsification_frozen.json \
python scripts/emeraldbay_recompute/pass4_falsification.py

echo "=== [6/6] pass5_gating_full -> updates $EB_CALIB_FIX ==="
python scripts/emeraldbay_recompute/pass5_gating_full.py

echo "=== FROZEN EMERALDBAY PIPELINE DONE ==="
