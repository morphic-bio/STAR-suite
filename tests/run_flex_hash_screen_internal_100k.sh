#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
STAR_BIN="${STAR_BIN:-${REPO_DIR}/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_$(date -u +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-16}"
TMP_ROOT="${TMP_ROOT:-/storage/100K/tmp}"
FLEX_PIPELINE_MODE="${FLEX_PIPELINE_MODE:-no}"
WRITE_BAM="${WRITE_BAM:-no}"

GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
HASH_CACHE="${HASH_CACHE:-}"
HASH_CACHE_SOURCE="${HASH_CACHE_SOURCE:-/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin}"
HASH_CACHE_SOURCE_SHA256="${HASH_CACHE_SOURCE_SHA256:-4e363bf902c46f6d3a5045343081cfab3b6ebc74c2efe4be2cfd86e7ce0d4c48}"
HASH_CACHE_SOURCE_SAMPLE="${HASH_CACHE_SOURCE_SAMPLE:-BC004}"
HASH_CACHE_INCLUDE_CLASSES="${HASH_CACHE_INCLUDE_CLASSES:-}"
HASH_CACHE_STAMP_CLASSES="${HASH_CACHE_STAMP_CLASSES:-0}"
RESTAMP_TOOL="${RESTAMP_TOOL:-${REPO_DIR}/scripts/restamp_flex_hash_cache.py}"
MEX_COMPARATOR="${MEX_COMPARATOR:-${SCRIPT_DIR}/compare_flex_hash_screen_mex.py}"
PARITY_MODE="${PARITY_MODE:-bounded}"
MAX_MISMATCH_FRACTION="${MAX_MISMATCH_FRACTION:-0.002}"
MAX_COUNT_DELTA_FRACTION="${MAX_COUNT_DELTA_FRACTION:-0.002}"
MAX_BARCODE_DIFFERENCE_FRACTION="${MAX_BARCODE_DIFFERENCE_FRACTION:-0.001}"
MAX_COORDINATE_DELTA="${MAX_COORDINATE_DELTA:-1}"
READS_R2="${READS_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz}"
READS_R1="${READS_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz}"

run_star() {
    local label="$1"
    shift
    local out_dir="${OUT_ROOT}/${label}"
    local tmp_dir="${TMP_ROOT}/flex_hash_screen_internal_${label}"
    local sam_args=(--outSAMtype None --outSAMattributes None)

    if [[ "${WRITE_BAM}" == "yes" ]]; then
        sam_args=(--outSAMtype BAM Unsorted --outSAMattributes None)
    elif [[ "${WRITE_BAM}" != "no" ]]; then
        echo "ERROR: WRITE_BAM must be yes or no" >&2
        return 2
    fi

    rm -rf "${out_dir}" "${tmp_dir}"
    mkdir -p "${out_dir}"

    echo "=== ${label} ==="
    echo "Output: ${out_dir}"
    echo "Tmp:    ${tmp_dir}"

    STAR_FLEX_HASH_SCREEN_CACHE="${HASH_CACHE}" \
    /usr/bin/time -f 'elapsed=%E user=%U sys=%S maxrss_kb=%M' \
      "${STAR_BIN}" \
      --runThreadN "${THREADS}" \
      --outTmpDir "${tmp_dir}" \
      --genomeDir "${GENOME_DIR}" \
      --soloType CB_UMI_Simple \
      --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
      --soloCBwhitelist "${CB_WHITELIST}" \
      --flex yes \
      --flexPipeline "${FLEX_PIPELINE_MODE}" \
      --soloFlexExpectedCellsPerTag 3000 \
      --soloSampleWhitelist "${SAMPLE_WHITELIST}" \
      --soloProbeList "${PROBE_LIST}" \
      --soloSampleProbes "${SAMPLE_PROBES}" \
      --soloSampleProbeOffset 68 \
      --soloFlexAllowedTags "${SAMPLE_WHITELIST}" \
      --soloFlexOutputPrefix "${out_dir}/per_sample" \
      --limitIObufferSize 50000000 50000000 \
      --outSJtype None \
      --soloMultiMappers Rescue \
      --alignIntronMax 500000 \
      --outFilterMismatchNmax 6 \
      --outFilterMismatchNoverReadLmax 1.0 \
      --outFilterMatchNmin 25 \
      --outSAMunmapped None \
      --outFilterMatchNminOverLread 0 \
      --outFilterMultimapNmax 10000 \
      --outFilterMultimapScoreRange 4 \
      --outSAMmultNmax 10000 \
      --winAnchorMultimapNmax 200 \
      --outSAMprimaryFlag AllBestScore \
      --outFilterScoreMin 0 \
      --outFilterScoreMinOverLread 0 \
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
      --soloUMIfiltering MultiGeneUMI_CR \
      --soloUMIdedup 1MM_CR \
      --soloCellFilter None \
      --clipAdapterType CellRanger4 \
      --soloFeatures Gene \
      --alignEndsType Local \
      --soloStrand Unstranded \
      --chimSegmentMin 1000000 \
      --soloKeysCompat cr \
      "${sam_args[@]}" \
      --soloSampleSearchNearby no \
      --readFilesIn "${READS_R2}" "${READS_R1}" \
      --outFileNamePrefix "${out_dir}/" \
      "$@" \
      >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"
}

prepare_hash_cache() {
    local cache_mode="external"
    local sample_index=""

    if [[ -z "${HASH_CACHE}" ]]; then
        local -a restamp_args
        cache_mode="restamped_v1"
        [[ -f "${HASH_CACHE_SOURCE}" ]] || {
            echo "ERROR: legacy Flex hash cache source not found: ${HASH_CACHE_SOURCE}" >&2
            return 1
        }
        [[ -f "${RESTAMP_TOOL}" ]] || {
            echo "ERROR: Flex hash-cache restamp tool not found: ${RESTAMP_TOOL}" >&2
            return 1
        }
        sample_index="$(awk -F '\t' -v sample="${HASH_CACHE_SOURCE_SAMPLE}" '$1 == sample { print NR }' "${SAMPLE_WHITELIST}")"
        [[ "${sample_index}" =~ ^[1-9][0-9]*$ ]] || {
            echo "ERROR: source sample ${HASH_CACHE_SOURCE_SAMPLE} is not uniquely present in ${SAMPLE_WHITELIST}" >&2
            return 1
        }
        HASH_CACHE="${OUT_ROOT}/cache/sequence_cache_${HASH_CACHE_SOURCE_SAMPLE}_full_whitelist.bin"
        mkdir -p "$(dirname "${HASH_CACHE}")"
        restamp_args=(
          --input "${HASH_CACHE_SOURCE}"
          --output "${HASH_CACHE}"
          --sample-index "${sample_index}"
          --stamp-classes "${HASH_CACHE_STAMP_CLASSES}"
          --expected-input-sha256 "${HASH_CACHE_SOURCE_SHA256}"
        )
        if [[ -n "${HASH_CACHE_INCLUDE_CLASSES}" ]]; then
            restamp_args+=(--include-classes "${HASH_CACHE_INCLUDE_CLASSES}")
        fi
        python3 "${RESTAMP_TOOL}" "${restamp_args[@]}" \
          >"${OUT_ROOT}/cache/restamp.log"
    fi

    [[ -s "${HASH_CACHE}" ]] || {
        echo "ERROR: Flex hash cache is missing or empty: ${HASH_CACHE}" >&2
        return 1
    }
    [[ "$(LC_ALL=C head -c 8 "${HASH_CACHE}")" == "FH01SEQ1" ]] || {
        echo "ERROR: Flex hash cache has invalid magic: ${HASH_CACHE}" >&2
        return 1
    }

    {
        printf 'cache_mode\t%s\n' "${cache_mode}"
        printf 'cache_path\t%s\n' "${HASH_CACHE}"
        printf 'cache_sha256\t%s\n' "$(sha256sum "${HASH_CACHE}" | awk '{print $1}')"
        printf 'source_path\t%s\n' "${HASH_CACHE_SOURCE}"
        printf 'source_expected_sha256\t%s\n' "${HASH_CACHE_SOURCE_SHA256}"
        printf 'source_sample\t%s\n' "${HASH_CACHE_SOURCE_SAMPLE}"
        printf 'source_sample_index\t%s\n' "${sample_index:-external}"
        printf 'cache_include_classes\t%s\n' "${HASH_CACHE_INCLUDE_CLASSES:-all}"
        printf 'cache_stamp_classes\t%s\n' "${HASH_CACHE_STAMP_CLASSES}"
        printf 'sample_whitelist\t%s\n' "${SAMPLE_WHITELIST}"
        printf 'sample_whitelist_sha256\t%s\n' "$(sha256sum "${SAMPLE_WHITELIST}" | awk '{print $1}')"
        printf 'star_binary\t%s\n' "${STAR_BIN}"
        printf 'star_sha256\t%s\n' "$(sha256sum "${STAR_BIN}" | awk '{print $1}')"
        printf 'flex_pipeline_mode\t%s\n' "${FLEX_PIPELINE_MODE}"
        printf 'write_bam\t%s\n' "${WRITE_BAM}"
    } >"${OUT_ROOT}/cache_manifest.tsv"
}

mkdir -p "${OUT_ROOT}" "${TMP_ROOT}"
prepare_hash_cache

run_star hash_on
run_star legacy --no-hash-screen yes

if [[ "${PARITY_MODE}" == "exact" ]]; then
    MAX_MISMATCH_FRACTION=0
    MAX_COUNT_DELTA_FRACTION=0
    MAX_BARCODE_DIFFERENCE_FRACTION=0
    MAX_COORDINATE_DELTA=0
elif [[ "${PARITY_MODE}" != "bounded" ]]; then
    echo "ERROR: PARITY_MODE must be bounded or exact" >&2
    exit 2
fi

python3 "${MEX_COMPARATOR}" \
    "${OUT_ROOT}/legacy" "${OUT_ROOT}/hash_on" \
    --max-mismatch-fraction "${MAX_MISMATCH_FRACTION}" \
    --max-count-delta-fraction "${MAX_COUNT_DELTA_FRACTION}" \
    --max-barcode-difference-fraction "${MAX_BARCODE_DIFFERENCE_FRACTION}" \
    --max-coordinate-delta "${MAX_COORDINATE_DELTA}" \
    | tee "${OUT_ROOT}/compare.log"

{
    echo "label\telapsed\tuser\tsys\tmaxrss_kb"
    for label in legacy hash_on; do
        if [ -f "${OUT_ROOT}/${label}/stderr.log" ]; then
            tail -n 1 "${OUT_ROOT}/${label}/stderr.log" | sed "s/^/${label}\t/"
        fi
    done
} > "${OUT_ROOT}/timing.tsv"

echo ""
echo "Outputs: ${OUT_ROOT}"
