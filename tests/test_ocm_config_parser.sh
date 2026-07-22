#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
UNIT_BIN="${REPO_ROOT}/core/legacy/source/ocm_multi_unit_tests"

if [[ ! -x "${UNIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" -j8 ocm-multi-unit-tests
fi

unset OCM_TEST_RUN_DIR OCM_TEST_CONFIG OCM_TEST_LOG
export OCM_TEST_FIXTURE_ROOT="${REPO_ROOT}/tests/fixtures/ocm_multi_tiny"
"${UNIT_BIN}" config sample_id
