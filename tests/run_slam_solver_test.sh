#!/bin/bash
# Build and run the SlamSolver unit test.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SLAM_DIR="$ROOT_DIR/slam"
CORE_DIR="$ROOT_DIR/core/legacy"

CXX="${CXX:-g++}"
OUT_BIN="${OUT_BIN:-/tmp/test_slam_solver}"

SRC_TEST="$ROOT_DIR/tests/slam/test_slam_solver.cpp"
SRC_SOLVER="$SLAM_DIR/source/SlamSolver.cpp"

if ! command -v "$CXX" >/dev/null 2>&1; then
    echo "FAIL: compiler not found: $CXX"
    exit 1
fi
if [[ ! -f "$SRC_TEST" ]]; then
    echo "FAIL: test source not found: $SRC_TEST"
    exit 1
fi
if [[ ! -f "$SRC_SOLVER" ]]; then
    echo "FAIL: solver source not found: $SRC_SOLVER"
    exit 1
fi

echo "=== Building SlamSolver test ==="
"$CXX" -std=c++11 -I"$SLAM_DIR/source" -I"$CORE_DIR/source" \
    -o "$OUT_BIN" "$SRC_TEST" "$SRC_SOLVER"

echo "=== Running SlamSolver test ==="
"$OUT_BIN"
