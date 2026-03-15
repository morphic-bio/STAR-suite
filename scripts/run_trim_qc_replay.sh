#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
BIN="${TRIM_QC_REPLAY_BIN:-${REPO_ROOT}/core/legacy/source/trim_qc_replay}"

if [[ ! -x "${BIN}" ]]; then
  echo "ERROR: missing executable ${BIN}" >&2
  echo "Build it with: make -C ${REPO_ROOT}/core/legacy/source trim_qc_replay" >&2
  exit 1
fi

exec "${BIN}" "$@"
