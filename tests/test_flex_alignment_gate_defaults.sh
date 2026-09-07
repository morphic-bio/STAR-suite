#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"

if [[ ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: STAR binary is not executable: ${STAR_BIN}" >&2
    exit 1
fi

TEST_ROOT="$(mktemp -d /tmp/star_flex_alignment_gate.XXXXXX)"
trap 'rm -rf "${TEST_ROOT}"' EXIT

printf '@read1\nACGT\n+\nIIII\n' >"${TEST_ROOT}/r2.fastq"
printf '@read1\nACGT\n+\nIIII\n' >"${TEST_ROOT}/r1.fastq"

run_parameter_probe() {
    local name="$1"
    shift
    local out_dir="${TEST_ROOT}/${name}"
    mkdir -p "${out_dir}"

    # Minimal reads plus a missing genome allow parameter initialization to
    # complete without starting an alignment. The assertions below concern
    # only the effective Flex defaults recorded during initialization.
    if "${STAR_BIN}" \
        --genomeDir "${TEST_ROOT}/missing-genome" \
        --readFilesIn "${TEST_ROOT}/r2.fastq" "${TEST_ROOT}/r1.fastq" \
        --outFileNamePrefix "${out_dir}/" \
        --flex yes \
        "$@" \
        >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"; then
        echo "ERROR: parameter probe unexpectedly completed an alignment" >&2
        return 1
    fi

    if [[ ! -f "${out_dir}/Log.out" ]]; then
        echo "ERROR: STAR did not write ${out_dir}/Log.out" >&2
        return 1
    fi
}

run_parameter_probe defaults
grep -Fq "outFilterScoreMinOverLread=0 (Flex default)" "${TEST_ROOT}/defaults/Log.out"
grep -Fq "outFilterMatchNminOverLread=0 (Flex default)" "${TEST_ROOT}/defaults/Log.out"

run_parameter_probe explicit \
    --outFilterScoreMinOverLread 0.75 \
    --outFilterMatchNminOverLread 0.8
grep -Fq "outFilterScoreMinOverLread=0.75 (explicit)" "${TEST_ROOT}/explicit/Log.out"
grep -Fq "outFilterMatchNminOverLread=0.8 (explicit)" "${TEST_ROOT}/explicit/Log.out"

echo "PASS: Flex alignment length gates default safely and preserve explicit overrides"
