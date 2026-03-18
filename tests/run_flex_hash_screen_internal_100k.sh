#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${SCRIPT_DIR}/.."
STAR_BIN="${STAR_BIN:-${REPO_DIR}/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_$(date -u +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-16}"
TMP_ROOT="${TMP_ROOT:-/storage/100K/tmp}"

GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
HASH_CACHE="${HASH_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin}"
READS_R2="${READS_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz}"
READS_R1="${READS_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz}"

run_star() {
    local label="$1"
    shift
    local out_dir="${OUT_ROOT}/${label}"
    local tmp_dir="${TMP_ROOT}/flex_hash_screen_internal_${label}"

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
      --outSAMattributes None \
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
      --outSAMtype None \
      --soloSampleSearchNearby no \
      --readFilesCommand zcat \
      --readFilesIn "${READS_R2}" "${READS_R1}" \
      --outFileNamePrefix "${out_dir}/" \
      "$@" \
      >"${out_dir}/stdout.log" 2>"${out_dir}/stderr.log"
}

mkdir -p "${OUT_ROOT}"

run_star hash_on
run_star legacy --no-hash-screen yes

"${SCRIPT_DIR}/compare_mex_outputs.sh" "${OUT_ROOT}/legacy" "${OUT_ROOT}/hash_on" | tee "${OUT_ROOT}/compare.log"

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
