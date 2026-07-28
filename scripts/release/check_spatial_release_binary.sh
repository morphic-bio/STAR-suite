#!/usr/bin/env bash

set -euo pipefail

BINARY=""
EXPECTED_VERSION=""
EXPECTED_COMMIT=""

usage() {
  cat <<USAGE
Usage: $0 --binary <STAR> --expected-version <version> --expected-commit <sha>

Checks immutable identity and the integrated spatial/Flex feature inventory in
an installed release binary.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2 ;;
    --expected-version) EXPECTED_VERSION="$2"; shift 2 ;;
    --expected-commit) EXPECTED_COMMIT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -x "$BINARY" ]] || { echo "ERROR: STAR binary is missing: $BINARY" >&2; exit 1; }
[[ -n "$EXPECTED_VERSION" ]] || { echo "ERROR: --expected-version is required" >&2; exit 1; }
[[ "$EXPECTED_COMMIT" =~ ^[0-9a-fA-F]{40}$ ]] || {
  echo "ERROR: --expected-commit must be a full 40-character SHA" >&2
  exit 1
}

[[ "$($BINARY --version)" == "$EXPECTED_VERSION" ]] || {
  echo "ERROR: STAR Suite version mismatch" >&2
  exit 1
}
source_revision="$($BINARY --source-revision)"
[[ "$source_revision" == "${EXPECTED_COMMIT,,}" ]] || {
  echo "ERROR: STAR source revision $source_revision does not match ${EXPECTED_COMMIT,,}" >&2
  exit 1
}

help_text="$($BINARY --help 2>&1)"
for parameter in \
  soloSpatialFlexIntegrated \
  soloSpatialGexIntegrated \
  soloFlexGdna \
  soloCrMultimapRescueEvidence \
  soloSpatialOverflowPolicy; do
  grep -q "$parameter" <<< "$help_text" || {
    echo "ERROR: release binary is missing parameter $parameter" >&2
    exit 1
  }
done

command -v nm >/dev/null 2>&1 || {
  echo "ERROR: nm is required for release feature-inventory validation" >&2
  exit 1
}
symbols="$(nm -C "$BINARY")"
for symbol in \
  'spatial_gex::Pipeline::completeCurrentThread' \
  'spatial_gex::Pipeline::finalize' \
  'flexResolveGeneIdx15_inlineResolver' \
  'FlexGdnaProbeMetadata'; do
  grep -q "$symbol" <<< "$symbols" || {
    echo "ERROR: release binary is missing symbol $symbol" >&2
    exit 1
  }
done

echo "PASS: STAR spatial release identity and feature inventory"
