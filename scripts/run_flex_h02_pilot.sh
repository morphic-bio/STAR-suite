#!/usr/bin/env bash
# One-shot H2 pilot: synthetic FASTQ from H0 seeds in an FH01SEQ1 cache → STAR-Flex →
# MEX-derived H2 sequence TSV → optional FH01SEQ1 binary (class 3 = H2 KEEP).
#
# Sharding: if unique_H0 * 11025 > max_reads_per_shard, increase --num-shards and run
# once per shard-id, then merge sequence_cache TSVs (same columns) before h2-write-binary-cache.
#
# Env:
#   PARENT_LIMIT   — cap unique H0 sequences (default 40 → ≤441k variants < 600k whitelist)
#   STAR_BIN       — STAR binary (default: repo core/legacy/source/STAR)
#   THREADS        — default 8
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
THREADS="${THREADS:-8}"
PARENT_LIMIT="${PARENT_LIMIT:-40}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"

STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_DIR:-/storage/downsampled_100K/SC2300771/results/flex_h02_pilot_${STAMP}}"
GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
BINARY_CACHE="${BINARY_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"

SHARD_ID="${SHARD_ID:-0}"
NUM_SHARDS="${NUM_SHARDS:-1}"
MAX_READS_PER_SHARD="${MAX_READS_PER_SHARD:-600000}"

mkdir -p "${OUT_DIR}/synthetic" "${OUT_DIR}/alignment"

echo "H2 pilot output: ${OUT_DIR}"
echo "STAR: ${STAR_BIN}  threads=${THREADS}  parent_limit=${PARENT_LIMIT}  shard=${SHARD_ID}/${NUM_SHARDS}"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-make-synth-fastq \
  --binary-cache "${BINARY_CACHE}" \
  --probe-list "${PROBE_LIST}" \
  --cb-whitelist "${CB_WHITELIST}" \
  --sample-whitelist "${SAMPLE_WHITELIST}" \
  --probe-offset 0 \
  --sample-offset 68 \
  --read-length 90 \
  --fill-base A \
  --num-shards "${NUM_SHARDS}" \
  --shard-id "${SHARD_ID}" \
  --max-reads-per-shard "${MAX_READS_PER_SHARD}" \
  --parent-limit "${PARENT_LIMIT}" \
  --r1-fastq "${OUT_DIR}/synthetic/h2_R1.fastq.gz" \
  --r2-fastq "${OUT_DIR}/synthetic/h2_R2.fastq.gz" \
  --manifest-tsv "${OUT_DIR}/synthetic/h2_manifest.tsv"

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
  --outFileNamePrefix "${OUT_DIR}/alignment/" \
  --readFilesCommand zcat \
  --readFilesIn "${OUT_DIR}/synthetic/h2_R2.fastq.gz" "${OUT_DIR}/synthetic/h2_R1.fastq.gz"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-build-cache-from-mex \
  --manifest-tsv "${OUT_DIR}/synthetic/h2_manifest.tsv" \
  --probe-list "${PROBE_LIST}" \
  --barcodes-tsv "${OUT_DIR}/alignment/Solo.out/Gene/raw/barcodes.tsv" \
  --features-tsv "${OUT_DIR}/alignment/Solo.out/Gene/raw/features.tsv" \
  --matrix-mtx "${OUT_DIR}/alignment/Solo.out/Gene/raw/matrix.mtx" \
  --qname-cache-tsv "${OUT_DIR}/h2_qname_cache.tsv" \
  --sequence-cache-tsv "${OUT_DIR}/h2_sequence_cache.tsv"

python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-write-binary-cache \
  --sequence-cache-tsv "${OUT_DIR}/h2_sequence_cache.tsv" \
  --output-bin "${OUT_DIR}/h2_keep_only.bin" \
  --sample-idx 0 \
  --keep-only

echo "Done. Sequence TSV: ${OUT_DIR}/h2_sequence_cache.tsv"
echo "H2 KEEP-only binary: ${OUT_DIR}/h2_keep_only.bin"
