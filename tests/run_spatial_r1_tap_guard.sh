#!/usr/bin/env bash
set -euo pipefail

STAR_BIN="${STAR_BIN:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../core/legacy/source" && pwd)/STAR}"
work="$(mktemp -d /tmp/star-spatial-r1-tap-guard.XXXXXX)"
trap 'rm -rf "${work}"' EXIT

set +e
"${STAR_BIN}" --soloSpatialR1FastqTap "${work}/not-opened.fifo" \
  >"${work}/stdout.log" 2>"${work}/stderr.log"
status=$?
set -e

[[ ${status} -ne 0 ]]
grep -q -- '--soloSpatialR1FastqTap requires --soloSpatialFeatureSidecar' \
  "${work}/stdout.log" "${work}/stderr.log"
if grep -q 'could not open Fastx input module' "${work}/stdout.log" "${work}/stderr.log"; then
  echo 'ERROR: fused-only option was rejected after opening ordinary read input' >&2
  exit 1
fi

echo 'spatial R1 tap fail-closed guard passed'
