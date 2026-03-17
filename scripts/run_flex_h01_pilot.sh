#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

THREADS="${THREADS:-8}"
if (( THREADS > 16 )); then
  echo "THREADS=${THREADS} exceeds requested cap; set THREADS<=16" >&2
  exit 1
fi

STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_DIR:-/storage/downsampled_100K/SC2300771/results/flex_h01_pilot_${STAMP}}"
STAR_SRC_DIR="${STAR_SRC_DIR:-${REPO_ROOT}/core/legacy/source}"
ISOLATED_REPO_ROOT="${OUT_DIR}/isolated_repo"
STAR_BUILD_DIR="${ISOLATED_REPO_ROOT}/core/legacy/source"
STAR_BIN="${STAR_BUILD_DIR}/STAR"

GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
PROBES_FASTA="${PROBES_FASTA:-/storage/flex_filtered_reference/filtered_reference/probes_only.fa}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
FEATURES_TSV="${FEATURES_TSV:-/storage/downsampled_100K/SC2300771/alignment/replay/flex_replay_features.tsv}"
MATRIX_MTX="${MATRIX_MTX:-/storage/downsampled_100K/SC2300771/alignment/replay/flex_replay_matrix.mtx}"
READS_R2="${READS_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz}"
READS_R1="${READS_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz}"

mkdir -p "${OUT_DIR}" "${OUT_DIR}/synthetic" "${OUT_DIR}/baseline"
rm -rf "${ISOLATED_REPO_ROOT}"
mkdir -p "${ISOLATED_REPO_ROOT}/core/legacy"
mkdir -p "${STAR_BUILD_DIR}"
cp -a "${STAR_SRC_DIR}/." "${STAR_BUILD_DIR}/"
ln -s "${REPO_ROOT}/flex" "${ISOLATED_REPO_ROOT}/flex"
ln -s "${REPO_ROOT}/slam" "${ISOLATED_REPO_ROOT}/slam"
mkdir -p "${ISOLATED_REPO_ROOT}/core/features"
ln -s "${REPO_ROOT}/core/features/vbem" "${ISOLATED_REPO_ROOT}/core/features/vbem"
ln -s "${REPO_ROOT}/core/features/yremove_fastq" "${ISOLATED_REPO_ROOT}/core/features/yremove_fastq"
ln -s "${REPO_ROOT}/core/features/bamsort" "${ISOLATED_REPO_ROOT}/core/features/bamsort"
ln -s "${REPO_ROOT}/core/features/process_features" "${ISOLATED_REPO_ROOT}/core/features/process_features"
ln -s "${REPO_ROOT}/core/features/libscrna" "${ISOLATED_REPO_ROOT}/core/features/libscrna"
make -C "${STAR_BUILD_DIR}" clean
make -C "${STAR_BUILD_DIR}" -j"${THREADS}" STAR

echo "Pilot output: ${OUT_DIR}"
echo "Isolated STAR build: ${STAR_BIN}"
echo "Threads: ${THREADS}"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" select-probes \
  --probe-list "${PROBE_LIST}" \
  --probes-fasta "${PROBES_FASTA}" \
  --features-tsv "${FEATURES_TSV}" \
  --matrix-mtx "${MATRIX_MTX}" \
  --count 100 \
  --output "${OUT_DIR}/pilot_probes.tsv"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" make-synth-fastq \
  --selected-probes "${OUT_DIR}/pilot_probes.tsv" \
  --cb-whitelist "${CB_WHITELIST}" \
  --sample-whitelist "${SAMPLE_WHITELIST}" \
  --probe-offset 0 \
  --sample-offset 68 \
  --read-length 90 \
  --fill-base A \
  --r1-fastq "${OUT_DIR}/synthetic/pilot_R1.fastq.gz" \
  --r2-fastq "${OUT_DIR}/synthetic/pilot_R2.fastq.gz" \
  --manifest-tsv "${OUT_DIR}/synthetic/pilot_manifest.tsv"

mkdir -p "${OUT_DIR}/synthetic/alignment"
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${CB_WHITELIST}" \
  --flex yes \
  --soloRunFlexFilter no \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist "${SAMPLE_WHITELIST}" \
  --soloSampleProbes "${SAMPLE_PROBES}" \
  --soloSampleProbeOffset 68 \
  --soloSampleSearchNearby no \
  --soloProbeList "${PROBE_LIST}" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
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
  --outSAMattributes NH HI AS nM NM GX GN ZG \
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
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "${OUT_DIR}/synthetic/alignment/" \
  --readFilesCommand zcat \
  --readFilesIn "${OUT_DIR}/synthetic/pilot_R2.fastq.gz" "${OUT_DIR}/synthetic/pilot_R1.fastq.gz"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" build-cache-from-mex \
  --manifest-tsv "${OUT_DIR}/synthetic/pilot_manifest.tsv" \
  --probe-list "${PROBE_LIST}" \
  --barcodes-tsv "${OUT_DIR}/synthetic/alignment/Solo.out/Gene/raw/barcodes.tsv" \
  --features-tsv "${OUT_DIR}/synthetic/alignment/Solo.out/Gene/raw/features.tsv" \
  --matrix-mtx "${OUT_DIR}/synthetic/alignment/Solo.out/Gene/raw/matrix.mtx" \
  --qname-cache-tsv "${OUT_DIR}/synthetic/qname_cache.tsv" \
  --sequence-cache-tsv "${OUT_DIR}/sequence_cache.tsv"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" scan-fastq \
  --sequence-cache-tsv "${OUT_DIR}/sequence_cache.tsv" \
  --offsets "0,1" \
  --output-tsv "${OUT_DIR}/scan_decisions.tsv" \
  --summary-tsv "${OUT_DIR}/scan_summary.tsv" \
  "${READS_R2}"

mkdir -p "${OUT_DIR}/baseline/alignment"
STAR_INLINE_REJECT_LOG="${OUT_DIR}/baseline/reject.tsv" \
STAR_INLINE_TRACE_QNAME=1 \
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
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
  --soloFlexOutputPrefix "${OUT_DIR}/baseline/per_sample" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
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
  --outSAMattributes NH HI AS nM NM GX GN ZG \
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
  --outSAMtype BAM Unsorted \
  --soloSampleSearchNearby no \
  --readFilesCommand zcat \
  --readFilesIn "${READS_R2}" "${READS_R1}" \
  --outFileNamePrefix "${OUT_DIR}/baseline/alignment/"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" compare-baseline \
  --probe-list "${PROBE_LIST}" \
  --scan-decisions-tsv "${OUT_DIR}/scan_decisions.tsv" \
  --baseline-reject-log "${OUT_DIR}/baseline/reject.tsv" \
  --summary-tsv "${OUT_DIR}/comparison_summary.tsv" \
  --mismatches-tsv "${OUT_DIR}/comparison_mismatches.tsv"

echo "Pilot completed: ${OUT_DIR}"
