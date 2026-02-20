#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="${ROOT_DIR}/tests/run_dynamic_threads_tiny_fixture.sh"

READ_MAP_NUMBER="${READ_MAP_NUMBER:-30000}"
OUT_BASE="${OUT_BASE:-/tmp/dynamic_threads_variable_sequences_$(date +%Y%m%d_%H%M%S)_$$}"

if [[ ! -x "${HARNESS}" ]]; then
    echo "ERROR: harness not found or not executable: ${HARNESS}" >&2
    exit 1
fi

run_case() {
    local case_name="$1"
    local run_threads="$2"
    local map_permits="$3"
    local sequence="$4"
    local expected_seq_len="$5"
    local case_out="${OUT_BASE}/${case_name}"

    echo "=== Variable sequence case: ${case_name} ==="
    RUN_THREADS="${run_threads}" \
    MAP_PERMITS="${map_permits}" \
    READ_MAP_NUMBER="${READ_MAP_NUMBER}" \
    VARIABLE_THREADS=1 \
    VARIABLE_THREADS_RETUNE_EVERY_ACQUIRES=1 \
    VARIABLE_THREADS_PERMIT_SEQUENCE="${sequence}" \
    MIN_RETUNES=1 \
    OUT_BASE="${case_out}" \
        "${HARNESS}"

    local report_json="${case_out}/dynamic_thread_report.json"
    python3 - "$report_json" "$run_threads" "$map_permits" "$expected_seq_len" "$sequence" <<'PY'
import json
import re
import sys
from pathlib import Path

report_path = Path(sys.argv[1])
expected_run_threads = int(sys.argv[2])
expected_initial_permits = int(sys.argv[3])
expected_seq_len = int(sys.argv[4])
expected_sequence = [int(x) for x in re.split(r"[\s,]+", sys.argv[5].strip()) if x]

if not report_path.is_file():
    raise SystemExit(f"missing report JSON: {report_path}")

report = json.loads(report_path.read_text(encoding="utf-8"))
interface = report.get("interface") or {}
telemetry = report.get("telemetry") or {}

assert interface.get("run_thread_n") == expected_run_threads, interface
assert interface.get("map_permits") == expected_initial_permits, interface
assert interface.get("variable_threads_enabled") is True, interface
assert interface.get("retune_every_acquires") == 1, interface
assert interface.get("retune_sequence_length") == expected_seq_len, interface

assert telemetry.get("retunes", 0) >= 1, telemetry
assert telemetry.get("retune_every_acquires") == 1, telemetry
assert telemetry.get("retune_sequence_length") == expected_seq_len, telemetry
assert telemetry.get("target_permits") is not None, telemetry
assert telemetry.get("configured_permits") is not None, telemetry
trace = telemetry.get("retune_trace") or []
assert len(trace) >= 1, telemetry
assert telemetry.get("retune_trace_dropped") == 0, telemetry
for idx, observed in enumerate(trace):
    expected = expected_sequence[idx % len(expected_sequence)]
    assert observed == expected, (idx, observed, expected, trace, expected_sequence)
PY
}

echo "=== Dynamic Threads Variable Sequence Smoke ==="
echo "HARNESS=${HARNESS}"
echo "OUT_BASE=${OUT_BASE}"
echo "READ_MAP_NUMBER=${READ_MAP_NUMBER}"

# Requested scenario: start at 3, then 2, then 4 (repeat)
run_case "sequence_3_2_4" 4 3 "2 4" 2

# Requested edge scenario: start at 1, then 2, then 1 (repeat)
run_case "sequence_1_2_1" 2 1 "2 1" 2

echo "PASS: variable sequence smoke completed"
echo "Outputs: ${OUT_BASE}"
