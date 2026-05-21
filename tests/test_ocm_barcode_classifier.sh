#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
UNIT_BIN="${REPO_ROOT}/core/legacy/source/ocm_multi_unit_tests"

if [[ ! -x "${UNIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" -j8 ocm-multi-unit-tests
fi

"${UNIT_BIN}" classifier
