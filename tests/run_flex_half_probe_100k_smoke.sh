#!/usr/bin/env bash
# Current release gate: H0/H1X2, matching model/filter axes, no genomic alignment.
# The March-2026 hash-versus-alignment test is a separate historical diagnostic.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export FLEX_SINGLE_RUN=1
export FLEX_REFERENCE_OUTPUT="${FLEX_REFERENCE_OUTPUT:-/mnt/pikachu/star_suite_v1100_gates_20260928/flex_modern_reference_v195a}"
export TEST_WORKDIR="${TEST_WORKDIR:-${OUT_ROOT:-/tmp/star_flex_half_probe_100k_$(date -u +%Y%m%d_%H%M%S)_$$}}"
exec bash "${ROOT}/tests/test_flex_v194_default_route.sh"
