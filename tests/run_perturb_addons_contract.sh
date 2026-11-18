#!/usr/bin/env bash
# Contract tests of the experimental STAR-perturb add-ons: the split-pool RNA
# mode (--soloSplitseqMode Fixed), the PIP RNA mode (--soloPipseqMode Illumina)
# and their list and reference tools. All inputs are generated; no external
# data is read.
#
#   STAR_BIN=/path/to/STAR OUT_ROOT=/fresh/dir tests/run_perturb_addons_contract.sh
#
# STAR_BIN defaults to this repository's build, OUT_ROOT to a fresh folder
# under /tmp. STAR_SUITE_SHARE names the share/star-suite folder that holds
# barcode_lists; it defaults to this repository's and then to
# /usr/local/share/star-suite (the container image).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_perturb_addons_contract_$(date -u +%Y%m%dT%H%M%SZ)_$$}"
SHARE="${STAR_SUITE_SHARE:-${ROOT_DIR}/share/star-suite}"
[[ -d "${SHARE}/barcode_lists" ]] || SHARE=/usr/local/share/star-suite

command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }
[[ -x "$STAR_BIN" ]] || { echo "ERROR: STAR binary not found at $STAR_BIN" >&2; exit 1; }
[[ -d "${SHARE}/barcode_lists/splitseq_v3" ]] || { echo "ERROR: no barcode_lists under $SHARE" >&2; exit 1; }
[[ ! -e "$OUT_ROOT" ]] || { echo "ERROR: use a fresh OUT_ROOT: $OUT_ROOT" >&2; exit 1; }
mkdir -p "$OUT_ROOT"

step() {
    local name="$1"; shift
    echo "=== ${name} ==="
    "$@" > "${OUT_ROOT}/${name}.log" 2>&1 || { tail -n 40 "${OUT_ROOT}/${name}.log" >&2; echo "FAIL: ${name}" >&2; exit 1; }
    tail -n 1 "${OUT_ROOT}/${name}.log"
}

cd "$ROOT_DIR"
step splitseq_smoke python3 tests/test_splitseq_smoke.py --star "$STAR_BIN" --out "${OUT_ROOT}/splitseq_smoke"
step splitseq_rounds python3 tests/test_splitseq_rounds.py --star "$STAR_BIN" --out "${OUT_ROOT}/splitseq_rounds"
step pipseq_smoke python3 tests/test_pipseq_smoke.py --star "$STAR_BIN" --out "${OUT_ROOT}/pipseq_smoke"
step pipseq_edges python3 tests/test_pipseq_edges.py --star "$STAR_BIN" --fixture "${OUT_ROOT}/pipseq_smoke" --out "${OUT_ROOT}/pipseq_edges"
step pipseq_reference_gtf python3 tests/test_make_pipseq_reference_gtf.py --out "${OUT_ROOT}/pipseq_reference_gtf"
step pipseq_tier_lists python3 tests/test_make_pipseq_tier_lists.py --out "${OUT_ROOT}/pipseq_tier_lists"
step splitseq_lists python3 tests/test_make_splitseq_lists.py --bundle "${SHARE}/barcode_lists/splitseq_v3" --out "${OUT_ROOT}/splitseq_lists"
echo "PASS: perturb add-on contract tests (${OUT_ROOT})"
