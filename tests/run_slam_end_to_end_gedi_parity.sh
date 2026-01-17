#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Wrapper to run the documented GEDI parity workflow without changing flags.
# See: tests/SLAM_E2E_IMPLEMENTATION_SUMMARY.md

WORK_DIR="${WORK_DIR:-/storage/slam_e2e_$(date +%Y%m%d)}"
export WORK_DIR

exec "$ROOT_DIR/tests/run_slam_end_to_end.sh" "$@"
