#!/usr/bin/env bash
set -euo pipefail

ROOT="${CR_MULTI_ROOT:-/storage/A375}"
FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${CR_MULTI_GEX_DIR:-${FASTQ_ROOT}/gex}"
NREADS="10000"
DOWN_SCRIPT="${CR_MULTI_DOWNSAMPLE_SCRIPT:-/mnt/pikachu/process_features/scripts/downsample_fastq_directory.sh}"

WHITELIST_GZ="${CR_MULTI_WHITELIST_GZ:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz}"
WHITELIST="${CR_MULTI_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
OUTPREFIX="${CR_MULTI_GEX_OUTPREFIX:-${ROOT}/star_gex_smoke/}"

if [[ ! -x "${DOWN_SCRIPT}" ]]; then
  echo "Missing downsample script: ${DOWN_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${WHITELIST}" ]]; then
  if [[ ! -f "${WHITELIST_GZ}" ]]; then
    echo "Missing whitelist: ${WHITELIST_GZ}" >&2
    exit 1
  fi
  zcat "${WHITELIST_GZ}" > "${WHITELIST}"
fi

if [[ ! -d "${GEX_DIR}/downsampled" ]] || ! ls "${GEX_DIR}/downsampled"/*.fastq* >/dev/null 2>&1; then
  "${DOWN_SCRIPT}" "${NREADS}" "${GEX_DIR}"
fi

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
  exit 1
fi

mkdir -p "${OUTPREFIX}"
chmod 775 "${OUTPREFIX}"

R1_FILES=$(ls "${GEX_DIR}/downsampled/"*R1*.fastq.gz | paste -sd, -)
R2_FILES=$(ls "${GEX_DIR}/downsampled/"*R2*.fastq.gz | paste -sd, -)

"${STAR_BIN}" \
  --runThreadN 4 \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_FILES}" "${R1_FILES}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUTPREFIX}" \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Rescue \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --soloStrand Unstranded \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --outSAMtype BAM Unsorted

OUT_GEX="${OUTPREFIX}/Solo.out/Gene/raw"
echo "GEX outputs in: ${OUT_GEX}"

OUT_GEX="${OUT_GEX}" python3 - <<'PY'
import os

out = os.environ["OUT_GEX"]
matrix = os.path.join(out, "matrix.mtx")
features = os.path.join(out, "features.tsv")
barcodes = os.path.join(out, "barcodes.tsv")
for path in (matrix, features, barcodes):
    if not os.path.exists(path):
        raise SystemExit(f"Missing output file: {path}")

sum_counts = 0
with open(matrix, "r", encoding="utf-8") as f:
    header_seen = False
    for line in f:
        if line.startswith("%"):
            continue
        if not header_seen:
            header_seen = True
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        sum_counts += int(float(parts[2]))

if sum_counts == 0:
    raise SystemExit("Gene counts are zero in GEX raw matrix")

print(f"Gene counts OK: sum={sum_counts}")
PY
