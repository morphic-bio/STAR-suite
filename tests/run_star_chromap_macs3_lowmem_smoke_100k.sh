#!/usr/bin/env bash
# Smoke: STAR + libchromap + MACS3 concurrent on 100K fixture, full lo-mem stack.
# Validates that --chromapAtacLowMem 1 + --chromapAtacMacs3FragLowMem 1 produce
# byte-identical narrowPeak vs the standalone chromap_callpeaks reference.
#
# Defaults:
#   STAR     := /mnt/pikachu/STAR-suite/core/legacy/source/STAR (built WITH_CHROMAP=1)
#   CALLPEAKS:= /tmp/chromap-suite-relink/chromap_callpeaks (or libMACS3/bin/macs3frag)
set -euo pipefail

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
CALLPEAKS="${CALLPEAKS:-/tmp/chromap-suite-relink/chromap_callpeaks}"

GENOME_DIR="${GENOME_DIR:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star}"
REF_FA="${REF_FA:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
CHROMAP_INDEX="${CHROMAP_INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
ATAC_WL="${ATAC_WL:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
GEX_WL="${GEX_WL:-/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k/open_source_full_20260424_015259/refs/737K-arc-v1_gex.txt}"
ATAC_DIR="${ATAC_DIR:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
GEX_DIR="${GEX_DIR:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/gex}"

OUT="${OUT:-$(mktemp -d /tmp/star_chromap_macs3_lowmem_smoke_100k.XXXXXX)}"
mkdir -p "${OUT}/out" "${OUT}/tmp"

ATAC_R1="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
ATAC_R2="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
ATAC_BC="${ATAC_DIR}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${ATAC_DIR}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"
GEX_R1="${GEX_DIR}/pbmc_unsorted_3k_S01_L003_R1_001.fastq.gz,${GEX_DIR}/pbmc_unsorted_3k_S01_L004_R1_001.fastq.gz"
GEX_R2="${GEX_DIR}/pbmc_unsorted_3k_S01_L003_R2_001.fastq.gz,${GEX_DIR}/pbmc_unsorted_3k_S01_L004_R2_001.fastq.gz"

echo "[smoke] outputs in ${OUT}"
echo "[smoke] STAR=${STAR_BIN}"

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
  --chromapAtacSecondaryFragments "${OUT}/out/atac_fragments.tsv.gz" \
  --chromapAtacOutputFormat BAM --chromapAtacSummary "${OUT}/out/chromap_summary.csv" \
  --chromapAtacThreads 8 --chromapAtacHtsThreads 2 \
  --chromapAtacSortBam 1 --chromapAtacWriteIndex 1 \
  --chromapAtacSortBamRam 2147483648 --chromapAtacTempDir "${OUT}/tmp" \
  --chromapAtacTn5ShiftMode classical \
  --chromapAtacLowMem 1 \
  --chromapAtacCallMacs3FragPeaks 1 \
  --chromapAtacMacs3FragPeaksOutput "${OUT}/out/atac_peaks.narrowPeak" \
  --chromapAtacMacs3FragSummitsOutput "${OUT}/out/atac_summits.bed" \
  --chromapAtacMacs3FragPeaksSource memory \
  --chromapAtacMacs3FragLowMem 1 \
  > "${OUT}/star.stdout.log" 2> "${OUT}/star.stderr.log"

echo "[smoke] STAR finished, validating peaks against standalone reference..."

# Standalone reference: chromap_callpeaks on the same fragments.
"${CALLPEAKS}" -i "${OUT}/out/atac_fragments.tsv.gz" \
  --frag-pileup-macs3-uint8-counts \
  --macs3-frag-narrowpeak "${OUT}/ref.narrowPeak" \
  --macs3-frag-summits "${OUT}/ref.summits.bed" \
  --bdgpeakcall-cutoff 5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
  --frag-score-pseudocount 0 \
  > "${OUT}/callpeaks.log" 2>&1

if cmp -s "${OUT}/out/atac_peaks.narrowPeak" "${OUT}/ref.narrowPeak" \
   && cmp -s "${OUT}/out/atac_summits.bed" "${OUT}/ref.summits.bed" ; then
  np_md5=$(md5sum "${OUT}/out/atac_peaks.narrowPeak" | cut -d' ' -f1)
  n=$(wc -l < "${OUT}/out/atac_peaks.narrowPeak")
  echo "PASS: lo-mem stack narrowPeak byte-identical to standalone (md5=${np_md5}, peaks=${n})"
else
  echo "FAIL: lo-mem stack output diverges from standalone reference"
  diff <(head -5 "${OUT}/out/atac_peaks.narrowPeak") <(head -5 "${OUT}/ref.narrowPeak") | head -20
  exit 1
fi
