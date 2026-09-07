#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_SOURCE_DIR="${STAR_SOURCE_DIR:-${REPO_ROOT}/core/legacy/source}"
SYNTH_SOURCE="${STAR_SOURCE_DIR}/ReadAlign_hashCacheSynth.cpp"
FLEX_RESOLVER_SOURCE="${STAR_SOURCE_DIR}/flex/SoloReadFeature_record_flex.cpp"

if [[ ! -f "${SYNTH_SOURCE}" ]]; then
    echo "ERROR: missing synthetic H1 verifier source: ${SYNTH_SOURCE}" >&2
    exit 1
fi

if [[ ! -f "${FLEX_RESOLVER_SOURCE}" ]]; then
    echo "ERROR: missing Flex resolver source: ${FLEX_RESOLVER_SOURCE}" >&2
    exit 1
fi

if grep -Eq 'P\.outFilter(ScoreMinOverLread|MatchNminOverLread)[[:space:]]*=' \
    "${SYNTH_SOURCE}"; then
    echo "ERROR: synthetic H1 verification writes shared alignment thresholds" >&2
    exit 1
fi

if ! grep -Fq 'mappedFilter(false);' "${SYNTH_SOURCE}"; then
    echo "ERROR: synthetic H1 verification does not select absolute-only filtering" >&2
    exit 1
fi

if grep -Fq 'gProbeIndexLoaded' "${FLEX_RESOLVER_SOURCE}"; then
    echo "ERROR: H1 resolver publishes the global probe index through an unsynchronized flag" >&2
    exit 1
fi

if ! grep -Fq 'std::call_once(gProbeIndexInitOnce' "${FLEX_RESOLVER_SOURCE}"; then
    echo "ERROR: H1 resolver does not initialize the global probe index exactly once" >&2
    exit 1
fi

echo "PASS: parallel H1 verification uses immutable thresholds and thread-safe probe-index initialization"
