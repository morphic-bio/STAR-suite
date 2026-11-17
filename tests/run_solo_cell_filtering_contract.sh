#!/usr/bin/env bash
# Contract test of cell filtering on a raw matrix (--runMode soloCellFiltering)
# and of --soloCellFilter TopCells: every barcode of the matrix is offered to
# the caller, TopCells N reports N barcodes plus ties, and its number has to
# be a whole number. The matrices are generated; no genome and no external
# data are read.
#
#   STAR_BIN=/path/to/STAR OUT_ROOT=/fresh/dir tests/run_solo_cell_filtering_contract.sh
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_solo_cell_filtering_contract_$(date -u +%Y%m%dT%H%M%SZ)_$$}"

command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }
[[ -x "$STAR_BIN" ]] || { echo "ERROR: STAR binary not found at $STAR_BIN" >&2; exit 1; }
[[ ! -e "$OUT_ROOT" ]] || { echo "ERROR: use a fresh OUT_ROOT: $OUT_ROOT" >&2; exit 1; }
mkdir -p "$OUT_ROOT"

cd "$ROOT_DIR"
python3 tests/test_solo_cell_filtering_last_barcode.py --star "$STAR_BIN" --out "${OUT_ROOT}/cell_filtering"
echo "PASS: cell filtering contract test (${OUT_ROOT})"
