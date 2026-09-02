#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOL_DIR="${ROOT_DIR}/flex/tools/hash_screen_replay"
CACHE="${CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin}"
DUMP="${DUMP:-/storage/downsampled_100K/SC2300771/results/flex_hash_screen_dump_20260320_005618/hash_screen_dump_v2.bin}"
CACHE_SHA256="${CACHE_SHA256:-4e363bf902c46f6d3a5045343081cfab3b6ebc74c2efe4be2cfd86e7ce0d4c48}"
DUMP_SHA256="${DUMP_SHA256:-c1cea03a9307e78b8ef2058a69015734f52385b1bf0ca9948f5347b03bd3276a}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_flex_hash_screen_replay}"

verify_sha256() {
    local expected="$1"
    local path="$2"
    local observed
    observed="$(sha256sum "${path}" | awk '{print $1}')"
    [[ "${observed}" == "${expected}" ]] || {
        echo "ERROR: SHA-256 mismatch for ${path}: expected ${expected}, got ${observed}" >&2
        return 1
    }
}

[[ -f "${CACHE}" ]] || { echo "ERROR: missing Flex hash cache: ${CACHE}" >&2; exit 1; }
[[ -f "${DUMP}" ]] || { echo "ERROR: missing Flex decision dump: ${DUMP}" >&2; exit 1; }
verify_sha256 "${CACHE_SHA256}" "${CACHE}"
verify_sha256 "${DUMP_SHA256}" "${DUMP}"

mkdir -p "${OUT_ROOT}"
make -C "${TOOL_DIR}" all
"${TOOL_DIR}/test_hash_screen" | tee "${OUT_ROOT}/unit.log"
"${TOOL_DIR}/hash_screen_replay" \
    "${CACHE}" "${DUMP}" --mode both --stats --diff "${OUT_ROOT}/diff.tsv" \
    | tee "${OUT_ROOT}/replay.log"

[[ "$(wc -l < "${OUT_ROOT}/diff.tsv")" -eq 1 ]] || {
    echo "ERROR: Flex replay emitted decision mismatches" >&2
    exit 1
}
grep -Fq 'Flat vs truth:        0 mismatches' "${OUT_ROOT}/replay.log"
grep -Fq 'Tiered vs truth:      0 mismatches' "${OUT_ROOT}/replay.log"
grep -Fq 'Flat vs tiered:       0 mismatches' "${OUT_ROOT}/replay.log"
grep -Fq 'Result: PASS' "${OUT_ROOT}/replay.log"

printf 'cache_sha256\t%s\n' "${CACHE_SHA256}" > "${OUT_ROOT}/provenance.tsv"
printf 'dump_sha256\t%s\n' "${DUMP_SHA256}" >> "${OUT_ROOT}/provenance.tsv"
echo "Flex hash-screen replay regression passed"
