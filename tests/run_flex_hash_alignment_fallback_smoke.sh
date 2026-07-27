#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
star_bin="${STAR_BIN:-${repo_dir}/core/legacy/source/STAR}"
genome_dir="${FLEX_HASH_FALLBACK_GENOME_DIR:-}"
cache="${FLEX_HASH_FALLBACK_CACHE:-${genome_dir}/flex_h01_sequence_cache.bin}"
probe_list="${FLEX_HASH_FALLBACK_PROBE_LIST:-${genome_dir}/probe_gene_list.txt}"
test_dir="$(mktemp -d /tmp/star-flex-hash-alignment-fallback.XXXXXX)"
trap 'rm -rf "${test_dir}"' EXIT

if [[ -z "${genome_dir}" ]]; then
    echo "FAIL: set FLEX_HASH_FALLBACK_GENOME_DIR to a STAR Flex index" >&2
    exit 2
fi
for required in "${star_bin}" "${genome_dir}/Genome" "${cache}" "${probe_list}"; do
    if [[ ! -e "${required}" ]]; then
        echo "FAIL: required input is absent: ${required}" >&2
        exit 2
    fi
done
command -v samtools >/dev/null 2>&1 || {
    echo "FAIL: samtools is required" >&2
    exit 2
}

STAR_INLINE_REJECT_LOG="${test_dir}/inline_trace.tsv" "${star_bin}" \
    --runThreadN 1 \
    --genomeDir "${genome_dir}" \
    --readFilesIn "${script_dir}/fixtures/flex_hash_fallback_r2.fastq" "${script_dir}/fixtures/flex_hash_fallback_r1.fastq" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist None \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloFeatures Gene \
    --soloProbeList "${probe_list}" \
    --soloInlineHashMode yes \
    --soloHashScreenFile "${cache}" \
    --flex yes \
    --flexPipeline no \
    --flexNoAlign 0 \
    --soloSkipProcessing yes \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped Within \
    --outSAMattributes NH HI AS nM \
    --outFileNamePrefix "${test_dir}/"

metric() {
    local label="$1"
    awk -F '|' -v label="${label}" '$1 ~ label "[[:space:]]*$" {gsub(/[[:space:]]/, "", $2); print $2}' "${test_dir}/Log.final.out"
}

[[ "$(metric 'Number of input reads')" == "3" ]]
[[ "$(metric 'Hash screen: KEEP')" == "2" ]]
[[ "$(metric 'Hash screen: KEEP, no barcode')" == "1" ]]
[[ "$(metric 'Hash screen: DENY')" == "0" ]]
[[ "$(metric 'Hash screen: PASS')" == "1" ]]
samtools quickcheck "${test_dir}/Aligned.out.bam"
[[ "$(samtools view -c "${test_dir}/Aligned.out.bam")" == "1" ]]
if [[ "$(awk -F '\t' 'NR > 1 && ($9 == "KEEP_HASH" || $9 == "RESOLVER_DROP") { n++ } END { print n+0 }' \
    "${test_dir}/inline_trace.tsv")" != "1" ]]; then
    echo "FAIL: the cache miss did not reach exactly one terminal feature-resolver decision" >&2
    exit 1
fi
if awk -F '\t' 'NR > 1 && $9 == "UNMATCHED_TAG" { found=1 } END { exit !found }' \
    "${test_dir}/inline_trace.tsv"; then
    echo "FAIL: disabled sample whitelist sentinel was treated as sample tagging" >&2
    exit 1
fi

echo "PASS: no-sample H0 and no-barcode H0 stay on hash path; only the cache miss reaches alignment and resolves feature evidence before barcode eligibility"
