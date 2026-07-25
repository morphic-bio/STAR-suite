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

set +e
"${STAR_BIN}" --soloSpatialGexIntegrated maybe \
  >"${work}/integrated-value.stdout.log" \
  2>"${work}/integrated-value.stderr.log"
status=$?
set -e

[[ ${status} -ne 0 ]]
grep -q -- '--soloSpatialGexIntegrated accepts only yes or no' \
  "${work}/integrated-value.stdout.log" "${work}/integrated-value.stderr.log"
if grep -q 'could not open Fastx input module' \
    "${work}/integrated-value.stdout.log" "${work}/integrated-value.stderr.log"; then
  echo 'ERROR: invalid integrated mode was rejected after opening ordinary read input' >&2
  exit 1
fi

set +e
"${STAR_BIN}" --soloSpatialGexIntegrated yes \
  --soloSpatialFeatureSidecar "${work}/diagnostic-sidecar" \
  --soloSpatialR1FastqTap "${work}/not-opened-integrated.fifo" \
  >"${work}/integrated-tap.stdout.log" \
  2>"${work}/integrated-tap.stderr.log"
status=$?
set -e

[[ ${status} -ne 0 ]]
grep -q -- '--soloSpatialR1FastqTap is not part of integrated spatial GEX mode' \
  "${work}/integrated-tap.stdout.log" "${work}/integrated-tap.stderr.log"
if grep -q 'could not open Fastx input module' \
    "${work}/integrated-tap.stdout.log" "${work}/integrated-tap.stderr.log"; then
  echo 'ERROR: integrated FIFO tap was rejected after opening ordinary read input' >&2
  exit 1
fi

echo 'spatial R1 tap and integrated-mode fail-closed guards passed'
