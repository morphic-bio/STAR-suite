#!/usr/bin/env bash
# Smoke: STAR + libchromap + MACS3 concurrent on 100K fixture, full lo-mem stack.
# Validates that --chromapAtacLowMem 1 + --chromapAtacMacs3FragLowMem 1 produce
# byte-identical narrowPeak vs the sidecar-native peak/MEX reference.
#
# Defaults:
#   STAR    := /mnt/pikachu/STAR-suite/core/legacy/source/STAR (Chromap-enabled build)
#   PEAK_MEX:= core/features/libchromap_contract/star_multiome_atac_peak_mex
set -euo pipefail

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
PEAK_MEX="${PEAK_MEX:-/mnt/pikachu/STAR-suite/core/features/libchromap_contract/star_multiome_atac_peak_mex}"

GENOME_DIR="${GENOME_DIR:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star}"
REF_FA="${REF_FA:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
CHROMAP_INDEX="${CHROMAP_INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
ATAC_WL="${ATAC_WL:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
GEX_WL="${GEX_WL:-/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k/open_source_full_20260424_015259/refs/737K-arc-v1_gex.txt}"
ATAC_DIR="${ATAC_DIR:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
GEX_DIR="${GEX_DIR:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/gex}"

OUT="${OUT:-$(mktemp -d /tmp/star_chromap_macs3_lowmem_smoke_100k.XXXXXX)}"
mkdir -p "${OUT}/out" "${OUT}/tmp"

if [[ ! -x "${PEAK_MEX}" ]]; then
  make -C "$(dirname "${PEAK_MEX}")" star_multiome_atac_peak_mex
fi

ATAC_R1="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
ATAC_R2="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
ATAC_BC="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"
GEX_R1="${GEX_DIR}/pbmc_unsorted_3k_S01_L003_R1_001.fastq.gz,${GEX_DIR}/pbmc_unsorted_3k_S01_L004_R1_001.fastq.gz"
GEX_R2="${GEX_DIR}/pbmc_unsorted_3k_S01_L003_R2_001.fastq.gz,${GEX_DIR}/pbmc_unsorted_3k_S01_L004_R2_001.fastq.gz"

echo "[smoke] outputs in ${OUT}"
echo "[smoke] STAR=${STAR_BIN}"

MACS3_QVALUE="${MACS3_QVALUE:-}"
MACS3_PVALUE="${MACS3_PVALUE:-}"
if [[ -n "${MACS3_QVALUE}" && -n "${MACS3_PVALUE}" ]]; then
  echo "ERROR: set only one of MACS3_QVALUE or MACS3_PVALUE" >&2
  exit 2
fi

star_macs3_threshold_args=()
peak_macs3_threshold_args=()
if [[ -n "${MACS3_QVALUE}" ]]; then
  star_macs3_threshold_args+=(--chromapAtacMacs3FragQvalue "${MACS3_QVALUE}")
  peak_macs3_threshold_args+=(--macs3-frag-qvalue "${MACS3_QVALUE}")
  echo "[smoke] MACS3 threshold: q=${MACS3_QVALUE}"
elif [[ -n "${MACS3_PVALUE}" ]]; then
  star_macs3_threshold_args+=(--chromapAtacMacs3FragPvalue "${MACS3_PVALUE}")
  peak_macs3_threshold_args+=(--macs3-frag-pvalue "${MACS3_PVALUE}")
  echo "[smoke] MACS3 threshold: p=${MACS3_PVALUE}"
else
  echo "[smoke] MACS3 threshold: default p=1e-5"
fi

/usr/bin/time -v "${STAR_BIN}" \
  --runThreadN 8 --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${GEX_R2}" "${GEX_R1}" --readFilesCommand zcat \
  --outFileNamePrefix "${OUT}/out/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CR UR CB UB GX GN \
  --limitBAMsortRAM 8000000000 --outBAMsortingThreadN 2 \
  --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 --soloCBwhitelist "${GEX_WL}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique --soloCbUbRequireTogether no \
  --soloCellFilter None --soloStrand Unstranded --soloFeatures Gene GeneFull \
  --chromapAtacEnable 1 --chromapAtacStartMode concurrent \
  --chromapAtacReferenceFasta "${REF_FA}" --chromapAtacIndex "${CHROMAP_INDEX}" \
  --chromapAtacRead1 "${ATAC_R1}" --chromapAtacRead2 "${ATAC_R2}" --chromapAtacBarcode "${ATAC_BC}" \
  --chromapAtacBarcodeWhitelist "${ATAC_WL}" \
  --chromapAtacOutputFragments "${OUT}/out/atac_possorted_bam.bam" \
  --chromapAtacSecondaryFragments "${OUT}/out/atac_fragments.bin" \
  --chromapAtacOutputFormat BAM --chromapAtacSummary "${OUT}/out/chromap_summary.csv" \
  --chromapAtacThreads 8 --chromapAtacHtsThreads 2 \
  --chromapAtacSortBam 1 --chromapAtacWriteIndex 1 \
  --chromapAtacSortBamRam 2147483648 --chromapAtacTempDir "${OUT}/tmp" \
  --chromapAtacTn5ShiftMode classical \
  --chromapAtacLowMem 1 \
  --chromapAtacCallMacs3FragPeaks 1 \
  --chromapAtacMacs3FragPeaksOutput "${OUT}/out/atac_peaks.narrowPeak" \
  --chromapAtacMacs3FragSummitsOutput "${OUT}/out/atac_summits.bed" \
  "${star_macs3_threshold_args[@]}" \
  --chromapAtacMacs3FragLowMem 1 \
  > "${OUT}/star.stdout.log" 2> "${OUT}/star.stderr.log"

echo "[smoke] STAR finished, validating peaks against sidecar reference..."

# Independent reference: production sidecar-native peak/MEX path on the same
# binary sidecar that STAR wrote during the Chromap run.
"${PEAK_MEX}" \
  --sidecar "${OUT}/out/atac_fragments.bin" \
  --call-peaks-from-sidecar \
  --peaks "${OUT}/ref.narrowPeak" \
  --summits-out "${OUT}/ref.summits.bed" \
  --out-dir "${OUT}/ref_peak_mex" \
  --metrics-tsv "${OUT}/ref_atac_metrics.tsv" \
  --threads 8 \
  "${peak_macs3_threshold_args[@]}" \
  --temp-dir "${OUT}/tmp" \
  --force \
  > "${OUT}/sidecar_peak_mex.log" 2>&1

if cmp -s "${OUT}/out/atac_peaks.narrowPeak" "${OUT}/ref.narrowPeak" \
   && cmp -s "${OUT}/out/atac_summits.bed" "${OUT}/ref.summits.bed" ; then
  np_md5=$(md5sum "${OUT}/out/atac_peaks.narrowPeak" | cut -d' ' -f1)
  n=$(wc -l < "${OUT}/out/atac_peaks.narrowPeak")
  echo "PASS: lo-mem stack narrowPeak byte-identical to sidecar reference (md5=${np_md5}, peaks=${n})"
else
  echo "FAIL: lo-mem stack output diverges from sidecar reference"
  diff <(head -5 "${OUT}/out/atac_peaks.narrowPeak") <(head -5 "${OUT}/ref.narrowPeak") | head -20
  exit 1
fi
