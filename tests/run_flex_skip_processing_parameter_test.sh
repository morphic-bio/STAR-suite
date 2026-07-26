#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
star_bin="${STAR_BIN:-${repo_dir}/core/legacy/source/STAR}"
fastq="${script_dir}/fixtures/trim_qc_fastq_tiny.fastq"
probe_list="${script_dir}/fixtures/flex_probe_gene_list_tiny.txt"
test_dir="$(mktemp -d /tmp/star-flex-skip-processing-parameters.XXXXXX)"
trap 'rm -rf "${test_dir}"' EXIT

common_args=(
    --genomeDir /nonexistent
    --readFilesIn "${fastq}" "${fastq}"
    --soloType CB_UMI_Simple
    --soloCBwhitelist None
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloFeatures Gene
    --soloProbeList "${probe_list}"
    --flex yes
)

set +e
"${star_bin}" "${common_args[@]}" \
    --soloSkipProcessing yes \
    --outFileNamePrefix "${test_dir}/skip/" \
    >"${test_dir}/skip.stdout" 2>"${test_dir}/skip.stderr"
skip_rc=$?

"${star_bin}" "${common_args[@]}" \
    --soloSkipProcessing no \
    --outFileNamePrefix "${test_dir}/process/" \
    >"${test_dir}/process.stdout" 2>"${test_dir}/process.stderr"
process_rc=$?

"${star_bin}" "${common_args[@]}" \
    --soloSkipProcessing yes \
    --soloRunFlexFilter no \
    --soloMapqMode off \
    --outFileNamePrefix "${test_dir}/explicit/" \
    >"${test_dir}/explicit.stdout" 2>"${test_dir}/explicit.stderr"
explicit_rc=$?
set -e

if [[ ${skip_rc} -eq 0 ]]; then
    echo "FAIL: fixture unexpectedly completed with a nonexistent genome" >&2
    exit 1
fi
if grep -Fq "FlexFilter requires expected cells count" "${test_dir}/skip.stderr"; then
    echo "FAIL: soloSkipProcessing=yes still requires expected cells" >&2
    exit 1
fi
if ! grep -Fq "could not open genome file" "${test_dir}/skip.stderr"; then
    echo "FAIL: skip-processing invocation did not advance beyond Solo validation" >&2
    sed -n '1,80p' "${test_dir}/skip.stderr" >&2
    exit 1
fi

if [[ ${process_rc} -eq 0 ]]; then
    echo "FAIL: processing fixture unexpectedly completed without expected cells" >&2
    exit 1
fi
if ! grep -Fq "FlexFilter requires expected cells count" "${test_dir}/process.stderr"; then
    echo "FAIL: active Flex processing did not require expected cells" >&2
    sed -n '1,80p' "${test_dir}/process.stderr" >&2
    exit 1
fi

if [[ ${explicit_rc} -eq 0 ]]; then
    echo "FAIL: explicit-override fixture unexpectedly completed with a nonexistent genome" >&2
    exit 1
fi
if ! grep -Fq "could not open genome file" "${test_dir}/explicit.stderr"; then
    echo "FAIL: explicit-override invocation did not advance beyond Solo validation" >&2
    sed -n '1,80p' "${test_dir}/explicit.stderr" >&2
    exit 1
fi
if ! grep -Fq "soloRunFlexFilter=no" "${test_dir}/explicit/Log.out"; then
    echo "FAIL: explicit --soloRunFlexFilter no was overridden" >&2
    exit 1
fi
if ! grep -Fq "soloMapqMode=off (enum=0)" "${test_dir}/explicit/Log.out"; then
    echo "FAIL: explicit --soloMapqMode off was overridden" >&2
    exit 1
fi

echo "PASS: skip-processing bypasses the expected-cell requirement; active processing retains it; explicit Flex policy overrides are preserved"
