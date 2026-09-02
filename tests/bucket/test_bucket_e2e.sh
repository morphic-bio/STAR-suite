#!/usr/bin/env bash
set -euo pipefail

CASE="${1:-}"
OUT_ROOT="${2:-/tmp/star_suite_cb_bucket_tests}"
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"

case "${CASE}" in B3|B4|B5|B6) ;; *) echo "usage: $0 B3|B4|B5|B6 [out-root]" >&2; exit 2 ;; esac

# Phase 2 replaces this guard with the STAR fixture/JAX-downsample comparisons.
if [[ ! -x "${ROOT_DIR}/core/legacy/source/cb_bucket_harness" ]]; then
    echo "FAIL: ${CASE} requires the bucket harness and its enabling implementation" >&2
    exit 1
fi

"${ROOT_DIR}/core/legacy/source/cb_bucket_harness" \
    --mode "${CASE}" --scratch "${OUT_ROOT}/${CASE}"
