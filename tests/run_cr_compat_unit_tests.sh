#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SUPPRESSIONS="${ROOT_DIR}/tests/compat/lsan_suppressions.txt"
TEST_BIN="${ROOT_DIR}/core/legacy/source/compat_unit_tests"

make -C "${ROOT_DIR}/core/legacy/source" -j ASAN=1 compat-unit-tests

cd "${ROOT_DIR}"
ASAN_OPTIONS="detect_leaks=1" \
LSAN_OPTIONS="detect_leaks=1:suppressions=${SUPPRESSIONS}" \
  "${TEST_BIN}"
