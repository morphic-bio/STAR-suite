#!/bin/bash
# Build and run the SlamCompat unit tests.
#
# Usage:
#   ./tests/run_slam_compat_test.sh
#
# Or from project root:
#   bash tests/run_slam_compat_test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CXX="${CXX:-g++}"
OUT_BIN="${OUT_BIN:-/tmp/test_slam_compat}"

SRC_TEST="$ROOT_DIR/tests/slam/test_slam_compat.cpp"
SRC_COMPAT="$ROOT_DIR/source/SlamCompat.cpp"
SRC_TRANSCRIPT="$ROOT_DIR/source/Transcript.cpp"

if ! command -v "$CXX" >/dev/null 2>&1; then
    echo "FAIL: compiler not found: $CXX"
    exit 1
fi
if [[ ! -f "$SRC_TEST" ]]; then
    echo "FAIL: test source not found: $SRC_TEST"
    exit 1
fi
if [[ ! -f "$SRC_COMPAT" ]]; then
    echo "FAIL: SlamCompat source not found: $SRC_COMPAT"
    exit 1
fi
if [[ ! -f "$SRC_TRANSCRIPT" ]]; then
    echo "FAIL: Transcript source not found: $SRC_TRANSCRIPT"
    exit 1
fi

echo "=== Building SlamCompat test ==="
# Compile with minimal dependencies (test uses mock data via test constructor)
"$CXX" -std=c++17 -I"$ROOT_DIR/source" -o "$OUT_BIN" \
    "$SRC_TEST" "$SRC_COMPAT" "$SRC_TRANSCRIPT" 2>&1

echo "=== Running SlamCompat test ==="
"$OUT_BIN"
