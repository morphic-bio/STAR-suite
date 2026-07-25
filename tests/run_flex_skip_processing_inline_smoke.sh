#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
star_bin="${STAR_BIN:-${repo_dir}/core/legacy/source/STAR}"
genome_dir="${FLEX_SKIP_PROCESSING_GENOME_DIR:-}"
fastq="${script_dir}/fixtures/trim_qc_fastq_tiny.fastq"
probe_list="${script_dir}/fixtures/flex_probe_gene_list_tiny.txt"
test_dir="$(mktemp -d /tmp/star-flex-skip-processing-inline.XXXXXX)"
trap 'rm -rf "${test_dir}"' EXIT

if [[ -z "${genome_dir}" ]]; then
    echo "FAIL: set FLEX_SKIP_PROCESSING_GENOME_DIR to a STAR Flex index" >&2
    exit 2
fi
if [[ ! -x "${star_bin}" ]]; then
    echo "FAIL: STAR binary is not executable: ${star_bin}" >&2
    exit 2
fi
if [[ ! -f "${genome_dir}/Genome" ]]; then
    echo "FAIL: Genome is absent from Flex index: ${genome_dir}" >&2
    exit 2
fi

"${star_bin}" \
    --runThreadN 1 \
    --genomeDir "${genome_dir}" \
    --readFilesIn "${fastq}" "${fastq}" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist None \
    --soloCBstart 1 \
    --soloCBlen 3 \
    --soloUMIstart 4 \
    --soloUMIlen 2 \
    --soloBarcodeReadLength 0 \
    --soloFeatures Gene \
    --soloProbeList "${probe_list}" \
    --flex yes \
    --flexPipeline no \
    --soloSkipProcessing yes \
    --outSAMtype None \
    --outFileNamePrefix "${test_dir}/"

if grep -Fq "SoloReadInfoLoader::load called with null streamReads" "${test_dir}/Log.out"; then
    echo "FAIL: inline skip-processing attempted legacy stream replay" >&2
    exit 1
fi
grep -Fq "skipping legacy readInfo stream replay" "${test_dir}/Log.out"
grep -Fq "skipping Solo counting (soloSkipProcessing=yes)" "${test_dir}/Log.out"
test -f "${test_dir}/Log.final.out"

echo "PASS: inline Flex skip-processing completes without legacy stream replay"
