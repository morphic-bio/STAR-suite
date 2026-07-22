#!/bin/bash
# Phase 2 gate script (runbook):
#   Gate A: SE fixture — STAR-Slam vs GRAND-SLAM (run_slam_fixture_parity.sh)
#   Gate B: PE symmetry on cloned fixture (regenerate_slam_fixture.sh --pe-clone)
#
# Gate B requires:
#   - STAR_BIN (default: build/STAR)
#   - STAR_INDEX (genome directory) — required for a real gate (no silent skip)
#   - Cloned PE FASTQs from --pe-clone (see PE_R1 / PE_R2 below)
#
# Optional:
#   RUN_SE_GATE (default 1), RUN_PE_GATE (default 1)
#   SLAM_PE_GATE_ALLOW_SKIP=1  — if STAR_INDEX missing, skip Gate B but do not fail the script
#
# Usage:
#   bash tests/run_slam_fixture_pe_parity.sh [--max-trim-delta N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
export STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
export SLAM_FIXTURE_ROOT="${SLAM_FIXTURE_ROOT:-$ROOT_DIR/test/fixtures/slam}"

RUN_SE_GATE="${RUN_SE_GATE:-1}"
RUN_PE_GATE="${RUN_PE_GATE:-1}"
FAILED=0
MAX_TRIM_DELTA=2
SLAM_PE_GATE_ALLOW_SKIP="${SLAM_PE_GATE_ALLOW_SKIP:-0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --max-trim-delta)
      MAX_TRIM_DELTA="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

PE_R1="${PE_R1:-$SLAM_FIXTURE_ROOT/raw/slam_pe_cloned_SRR32576116_R1_001.fastq.gz}"
PE_R2="${PE_R2:-$SLAM_FIXTURE_ROOT/raw/slam_pe_cloned_SRR32576116_R2_001.fastq.gz}"

echo "======================================================================"
echo "run_slam_fixture_pe_parity.sh"
echo "ROOT_DIR=$ROOT_DIR  MAX_TRIM_DELTA=$MAX_TRIM_DELTA"
echo "======================================================================"

if [[ "$RUN_SE_GATE" == "1" ]]; then
  echo ""
  echo "--- Gate A: SE no-regression (fixture parity) ---"
  if bash "$SCRIPT_DIR/run_slam_fixture_parity.sh"; then
    echo "Gate A: PASS"
  else
    echo "Gate A: FAIL" >&2
    FAILED=1
  fi
else
  echo "Gate A: skipped (RUN_SE_GATE=0)"
fi

if [[ "$RUN_PE_GATE" != "1" ]]; then
  echo "Gate B: skipped (RUN_PE_GATE=0)"
  exit "$FAILED"
fi

echo ""
echo "--- Gate B: PE trim symmetry (cloned fixture) ---"

if [[ ! -x "$STAR_BIN" ]]; then
  echo "Gate B: FAIL — STAR binary missing or not executable: $STAR_BIN" >&2
  exit 1
fi

if [[ ! -f "$PE_R1" || ! -f "$PE_R2" ]]; then
  echo "Gate B: FAIL — PE clone fixture missing." >&2
  echo "  Expected: $PE_R1" >&2
  echo "  Expected: $PE_R2" >&2
  echo "  Run: bash tests/regenerate_slam_fixture.sh --pe-clone" >&2
  exit 1
fi

STAR_INDEX="${STAR_INDEX:-}"
if [[ -z "$STAR_INDEX" || ! -d "$STAR_INDEX" ]]; then
  if [[ "$SLAM_PE_GATE_ALLOW_SKIP" == "1" ]]; then
    echo "SKIP Gate B: STAR_INDEX not set or not a directory (SLAM_PE_GATE_ALLOW_SKIP=1)"
    exit "$FAILED"
  fi
  echo "Gate B: FAIL — set STAR_INDEX to a genome directory (STAR --genomeDir)." >&2
  exit 1
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/slam_pe_gate.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

OUT="$WORK/out_"
QC_JSON="$WORK/slam_qc_gate.json"

echo "  STAR_BIN=$STAR_BIN"
echo "  STAR_INDEX=$STAR_INDEX"
echo "  readFilesIn: $PE_R1 $PE_R2"

set +e
"$STAR_BIN" \
  --runThreadN 2 \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$PE_R1" "$PE_R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT" \
  --outSAMtype None \
  --slamQuantMode 1 \
  --autoTrim variance \
  --autoTrimMaxReads 8000 \
  --autoTrimMinReads 50 \
  --autoTrimDetectionReads 12000 \
  --slamAutoTrimPerMate Yes \
  --slamQcJson "$QC_JSON" \
  >"$WORK/star.log" 2>&1
st=$?
set -e

if [[ "$st" -ne 0 ]]; then
  echo "Gate B: FAIL — STAR exit $st (log: $WORK/star.log)" >&2
  tail -40 "$WORK/star.log" >&2 || true
  exit 1
fi

if [[ ! -f "$QC_JSON" ]]; then
  ALT="${OUT}slam_qc.json"
  if [[ -f "$ALT" ]]; then
    QC_JSON=$ALT
  else
    echo "Gate B: FAIL — QC JSON not produced under $WORK" >&2
    exit 1
  fi
fi

python3 - "$QC_JSON" "$MAX_TRIM_DELTA" <<'PY'
import json
import sys

path, dmax_s = sys.argv[1], sys.argv[2]
dmax = int(dmax_s)
with open(path) as f:
    d = json.load(f)

if d.get("variance_histogram_mode") != "per_mate_separate":
    print("FAIL: variance_histogram_mode expected per_mate_separate, got", d.get("variance_histogram_mode"), file=sys.stderr)
    sys.exit(1)

def limpair(k1, k2):
    a = int(d.get(k1, 0))
    b = int(d.get(k2, 0))
    return abs(a - b)

pairs = [
    ("trim5p_mate1", "trim5p_mate2"),
    ("trim3p_mate1", "trim3p_mate2"),
]
worst = 0
for k1, k2 in pairs:
    dd = limpair(k1, k2)
    if dd > worst:
        worst = dd
    print(f"  |{k1} - {k2}| = {dd}")

if worst > dmax:
    print(f"FAIL: trim delta {worst} > --max-trim-delta {dmax}", file=sys.stderr)
    sys.exit(1)

mates = d.get("mates") or []
if len(mates) != 2:
    print("FAIL: mates array length", len(mates), file=sys.stderr)
    sys.exit(1)
print("Gate B: PASS")
PY

echo "======================================================================"
exit "$FAILED"
