#!/usr/bin/env bash
# OCM streaming materializer memory smoke (tiny fixture + optional profile log parse).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
FIXTURE_ROOT="${OCM_TEST_FIXTURE_ROOT:-${REPO_ROOT}/tests/fixtures/ocm_multi_tiny}"
UNIT_BIN="${REPO_ROOT}/core/legacy/source/ocm_multi_unit_tests"
LOG="${OCM_MATERIALIZER_MEM_LOG:-/tmp/ocm_materializer_mem.log}"

export STAR_SOLO_MEMORY_PROFILE=1
export OCM_TEST_FIXTURE_ROOT="${FIXTURE_ROOT}"
export OCM_TEST_LOG="${LOG}"
unset OCM_TEST_RUN_DIR OCM_TEST_CONFIG

rm -rf "${FIXTURE_ROOT}/run/outs" "${FIXTURE_ROOT}/samples"
[[ -x "${UNIT_BIN}" ]] || make -C "${REPO_ROOT}/core/legacy/source" -j8 ocm-multi-unit-tests

"${UNIT_BIN}" materialize
"${REPO_ROOT}/tests/test_ocm_mex_materializer_tiny.sh"

echo "=== OCM materializer memory checkpoints (${LOG}) ==="
rg '^Solo memory: ocm_' "${LOG}" || true
echo
"${REPO_ROOT}/tests/run_solo_memory_profile_harness.sh" --parse-log "${LOG}" 2>/dev/null | rg 'ocm_' || true

echo "PASS: OCM streaming materializer memory smoke"
