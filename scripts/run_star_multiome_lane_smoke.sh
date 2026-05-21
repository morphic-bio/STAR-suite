#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"
NORMALIZE_ATAC_BC="${REPO_ROOT}/scripts/normalize_multiome_atac_barcode_fastq.py"
PACKAGE_GENEFULL="${REPO_ROOT}/scripts/package_star_genefull_mex.py"
PREPARE_VELOCYTO="${REPO_ROOT}/scripts/prepare_velocyto_mex.py"
ALLOW_LEGACY_PREPARE_VELOCYTO="${ALLOW_LEGACY_PREPARE_VELOCYTO:-0}"
REMOTE_POST_MEX="${REPO_ROOT}/scripts/run_remote_multiome_post_mex_rsync.sh"
LOCAL_DOWNSTREAM="${REPO_ROOT}/scripts/run_scrna_downstream_gene_full_velocyto.sh"
BUILD_ATAC_MEX_NATIVE="${REPO_ROOT}/core/features/libchromap_contract/star_multiome_atac_peak_mex"
BUILD_ATAC_MEX_PY="${REPO_ROOT}/scripts/build_atac_peak_matrix_from_fragments.py"
BUILD_MUDATA="${REPO_ROOT}/scripts/build_multiome_mudata.py"

GEX_R1=""
GEX_R2=""
ATAC_R1=""
ATAC_BARCODE=""
ATAC_R2=""
OUT_DIR=""
THREADS="${STAR_MULTIOME_THREADS:-16}"
CHROMAP_THREADS="${STAR_MULTIOME_CHROMAP_THREADS:-8}"
GENOME_DIR="${STAR_MULTIOME_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
GEX_WHITELIST="${STAR_MULTIOME_GEX_WHITELIST:-/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k/open_source_full_20260424_015259/refs/737K-arc-v1_gex.txt}"
CHROMAP_REF="${STAR_MULTIOME_CHROMAP_REF:-/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa}"
CHROMAP_INDEX="${STAR_MULTIOME_CHROMAP_INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
ATAC_WHITELIST="${STAR_MULTIOME_ATAC_WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
ATAC_TO_GEX="${STAR_MULTIOME_ATAC_TO_GEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv}"
CHROMAP_LOW_MEM="${STAR_MULTIOME_CHROMAP_LOW_MEM:-0}"
CHROMAP_LOW_MEM_RAM="${STAR_MULTIOME_CHROMAP_LOW_MEM_RAM:-0}"
CHROMAP_MACS3_FRAG_LOW_MEM="${STAR_MULTIOME_CHROMAP_MACS3_FRAG_LOW_MEM:-0}"
CHROMAP_START_MODE="${STAR_MULTIOME_CHROMAP_START_MODE:-concurrent}"
SOLO_STRAND="${STAR_MULTIOME_SOLO_STRAND:-Forward}"
ATAC_BARCODE_START="${STAR_MULTIOME_ATAC_BARCODE_START:-9}"
ATAC_BARCODE_LENGTH="${STAR_MULTIOME_ATAC_BARCODE_LENGTH:-16}"
ATAC_BARCODE_READ_FORMAT="${STAR_MULTIOME_ATAC_READ_FORMAT:-}"
REMOTE_HOST="${MULTIOME_REMOTE_HOST:-}"
REMOTE_ROOT="${MULTIOME_REMOTE_ROOT:-}"
REMOTE_OUTPUT_NAME="${MULTIOME_REMOTE_OUTPUT_NAME:-downstream_genefull_velocyto_cellbender}"
CELLBENDER_CPU_CORES="${MULTIOME_CELLBENDER_CPU_CORES:-}"
CELLBENDER_GPU="0"
NO_SYNC_IMAGES="0"
KEEP_REMOTE="0"
ALLOW_LOCAL_DOWNSTREAM="0"
FORCE="0"
FORCE_ATAC_BARCODE="0"
SKIP_BUILD="0"
STOP_AFTER_LOCAL_MEX="0"
SOLO_INLINE_HASH="0"
USE_NATIVE_ATAC_BARCODE="${STAR_MULTIOME_USE_NATIVE_ATAC_BARCODE:-1}"
USE_NATIVE_ATAC_MEX="${STAR_MULTIOME_USE_NATIVE_ATAC_MEX:-1}"

usage() {
  cat <<'EOF'
Usage:
  run_star_multiome_lane_smoke.sh --gex-r1 PATH --gex-r2 PATH \
    --atac-r1 PATH --atac-barcode PATH --atac-r2 PATH --out-dir PATH [options]

Runs a 10x Multiome STAR-suite sample workflow; FASTQ arguments may be a single
path or a comma-separated multi-lane path list:
  1. run Chromap against the raw ATAC i5/barcode read using native read-format support
  2. run STAR GEX with GeneFull+Velocyto, Y/noY outputs, and in-process Chromap ATAC
  3. package GeneFull/Velocyto MEX
  4. materialize ATAC peak MEX locally from the Chromap binary sidecar
  5. run post-MEX downstream h5ad, CellBender, and MuData remotely when configured

Required:
  --gex-r1 PATH             GEX R1 FASTQ containing CB/UMI
  --gex-r2 PATH             GEX R2 FASTQ containing cDNA
  --atac-r1 PATH            ATAC genomic read 1 FASTQ
  --atac-barcode PATH       ATAC 24 bp i5/barcode read FASTQ
  --atac-r2 PATH            ATAC genomic read 2 / R3 FASTQ
  --out-dir PATH            Fresh output directory

Reference/defaults:
  --genome-dir PATH         STAR GEX genomeDir (default: /storage/autoindex_110_44/bulk_index)
  --gex-whitelist PATH      737K ARC GEX whitelist
  --chromap-ref PATH        Chromap FASTA, default from STAR reference set
  --chromap-index PATH      Chromap index
  --atac-whitelist PATH     737K ARC ATAC whitelist
  --atac-to-gex PATH        two-column ATAC->GEX barcode translation table
  --solo-strand STR         STARsolo strand (default: Forward)

Remote downstream:
  --remote-host HOST
  --remote-root PATH
  --remote-output-name NAME
  --cellbender-cpu-cores N
  --cellbender-gpu
  --no-sync-images
  --keep-remote
  --allow-local-downstream  Permit local downstream without CellBender
  --stop-after-local-mex    Stop after local STAR MEX packaging and ATAC peak MEX

Other:
  --threads N
  --chromap-threads N
  --chromap-low-mem        Enable STAR/Chromap low-memory overflow-spill mode
  --chromap-low-mem-ram N  RAM threshold in bytes for low-memory spill; 0 uses Chromap defaults
  --chromap-macs3-frag-low-mem
                           Enable low-memory MACS3 fragment peak workspace
  --chromap-start-mode MODE
                           STAR/Chromap scheduling: postMapping or concurrent
                           (default: concurrent)
  --atac-barcode-start N    1-based barcode window start in ATAC barcode read (default: 9)
  --atac-barcode-length N   barcode length (default: 16)
  --chromap-atac-read-format FORMAT
                            Native Chromap read format (default derived as bc:8:23:-)
  --normalize-atac-barcode  Use the legacy Python barcode FASTQ normalizer fallback
  --python-atac-mex         Use the legacy Python/bedtools ATAC peak-MEX fallback
  --skip-build              Do not rebuild STAR WITH_CHROMAP=1
  --force                   Regenerate outputs
  --force-atac-barcode      Regenerate normalized ATAC barcode FASTQ
  --solo-inline-hash        Enable STARsolo inline hash mode (off by default)
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gex-r1) GEX_R1="$2"; shift 2 ;;
    --gex-r2) GEX_R2="$2"; shift 2 ;;
    --atac-r1) ATAC_R1="$2"; shift 2 ;;
    --atac-barcode) ATAC_BARCODE="$2"; shift 2 ;;
    --atac-r2) ATAC_R2="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --chromap-threads) CHROMAP_THREADS="$2"; shift 2 ;;
    --chromap-low-mem) CHROMAP_LOW_MEM="1"; shift ;;
    --chromap-low-mem-ram) CHROMAP_LOW_MEM_RAM="$2"; shift 2 ;;
    --chromap-macs3-frag-low-mem) CHROMAP_MACS3_FRAG_LOW_MEM="1"; shift ;;
    --chromap-start-mode) CHROMAP_START_MODE="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --gex-whitelist) GEX_WHITELIST="$2"; shift 2 ;;
    --chromap-ref) CHROMAP_REF="$2"; shift 2 ;;
    --chromap-index) CHROMAP_INDEX="$2"; shift 2 ;;
    --atac-whitelist) ATAC_WHITELIST="$2"; shift 2 ;;
    --atac-to-gex) ATAC_TO_GEX="$2"; shift 2 ;;
    --solo-strand) SOLO_STRAND="$2"; shift 2 ;;
    --atac-barcode-start) ATAC_BARCODE_START="$2"; shift 2 ;;
    --atac-barcode-length) ATAC_BARCODE_LENGTH="$2"; shift 2 ;;
    --chromap-atac-read-format) ATAC_BARCODE_READ_FORMAT="$2"; shift 2 ;;
    --remote-host) REMOTE_HOST="$2"; shift 2 ;;
    --remote-root) REMOTE_ROOT="$2"; shift 2 ;;
    --remote-output-name) REMOTE_OUTPUT_NAME="$2"; shift 2 ;;
    --cellbender-cpu-cores) CELLBENDER_CPU_CORES="$2"; shift 2 ;;
    --cellbender-gpu) CELLBENDER_GPU="1"; shift ;;
    --no-sync-images) NO_SYNC_IMAGES="1"; shift ;;
    --keep-remote) KEEP_REMOTE="1"; shift ;;
    --allow-local-downstream) ALLOW_LOCAL_DOWNSTREAM="1"; shift ;;
    --stop-after-local-mex) STOP_AFTER_LOCAL_MEX="1"; shift ;;
    --skip-build) SKIP_BUILD="1"; shift ;;
    --force) FORCE="1"; shift ;;
    --force-atac-barcode) FORCE_ATAC_BARCODE="1"; shift ;;
    --solo-inline-hash) SOLO_INLINE_HASH="1"; shift ;;
    --normalize-atac-barcode) USE_NATIVE_ATAC_BARCODE="0"; shift ;;
    --python-atac-mex) USE_NATIVE_ATAC_MEX="0"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

for required_name in GEX_R1 GEX_R2 ATAC_R1 ATAC_BARCODE ATAC_R2 OUT_DIR; do
  [[ -n "${!required_name}" ]] || { echo "ERROR: --${required_name,,} is required" >&2; exit 1; }
done

case "${CHROMAP_START_MODE}" in
  postMapping|concurrent) ;;
  *) echo "ERROR: --chromap-start-mode must be postMapping or concurrent" >&2; exit 1 ;;
esac

realpath_csv() {
  local csv="$1"
  local out=()
  local field
  IFS=',' read -r -a fields <<< "${csv}"
  for field in "${fields[@]}"; do
    [[ -n "${field}" ]] || { echo "ERROR: empty path in CSV: ${csv}" >&2; exit 1; }
    out+=("$(realpath "${field}")")
  done
  local joined
  joined="$(IFS=','; echo "${out[*]}")"
  echo "${joined}"
}

check_csv_paths_exist() {
  local csv="$1"
  local field
  IFS=',' read -r -a fields <<< "${csv}"
  for field in "${fields[@]}"; do
    [[ -e "${field}" ]] || { echo "ERROR: missing ${field}" >&2; exit 1; }
  done
}

first_csv_path() {
  local csv="$1"
  local first
  IFS=',' read -r first _ <<< "${csv}"
  echo "${first}"
}

GEX_R1="$(realpath_csv "${GEX_R1}")"
GEX_R2="$(realpath_csv "${GEX_R2}")"
ATAC_R1="$(realpath_csv "${ATAC_R1}")"
ATAC_BARCODE="$(realpath_csv "${ATAC_BARCODE}")"
ATAC_R2="$(realpath_csv "${ATAC_R2}")"
GENOME_DIR="$(realpath "${GENOME_DIR}")"
GEX_WHITELIST="$(realpath "${GEX_WHITELIST}")"
CHROMAP_REF="$(realpath "${CHROMAP_REF}")"
CHROMAP_INDEX="$(realpath "${CHROMAP_INDEX}")"
ATAC_WHITELIST="$(realpath "${ATAC_WHITELIST}")"
ATAC_TO_GEX="$(realpath "${ATAC_TO_GEX}")"
OUT_DIR="$(realpath -m "${OUT_DIR}")"

for csv_paths in "${GEX_R1}" "${GEX_R2}" "${ATAC_R1}" "${ATAC_BARCODE}" "${ATAC_R2}"
do
  check_csv_paths_exist "${csv_paths}"
done
for path in "${GENOME_DIR}" "${GEX_WHITELIST}" "${CHROMAP_REF}" "${CHROMAP_INDEX}" \
  "${ATAC_WHITELIST}" "${ATAC_TO_GEX}" "${PACKAGE_GENEFULL}" \
  "${BUILD_MUDATA}"
do
  [[ -e "${path}" ]] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done
if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
  [[ -e "${PREPARE_VELOCYTO}" ]] || { echo "ERROR: missing ${PREPARE_VELOCYTO}" >&2; exit 1; }
fi
if [[ "${USE_NATIVE_ATAC_BARCODE}" != "1" ]]; then
  [[ -e "${NORMALIZE_ATAC_BC}" ]] || { echo "ERROR: missing ${NORMALIZE_ATAC_BC}" >&2; exit 1; }
fi
if [[ "${USE_NATIVE_ATAC_MEX}" != "1" ]]; then
  [[ -e "${BUILD_ATAC_MEX_PY}" ]] || { echo "ERROR: missing ${BUILD_ATAC_MEX_PY}" >&2; exit 1; }
fi

mkdir -p "${OUT_DIR}/logs" "${OUT_DIR}/atac" "${OUT_DIR}/mudata"
SAMPLE_DIR="${OUT_DIR}/star_sample"
RUN_DIR="${SAMPLE_DIR}/run"
ATAC_BC_NORM="${OUT_DIR}/atac/$(basename "$(first_csv_path "${ATAC_BARCODE}")" .fastq.gz).arc_atac_bc.fastq.gz"
ATAC_BC_FOR_CHROMAP="${ATAC_BARCODE}"
ATAC_BAM="${RUN_DIR}/atac_possorted.bam"
ATAC_SIDECAR="${RUN_DIR}/atac_fragments.bin"
ATAC_PEAKS="${RUN_DIR}/atac_peaks.narrowPeak"
ATAC_SUMMITS="${RUN_DIR}/atac_summits.bed"
ATAC_MEX="${OUT_DIR}/atac/peak_mex"
ATAC_METRICS="${OUT_DIR}/atac/atac_metrics.tsv"
DOWNSTREAM_DIR="${SAMPLE_DIR}/${REMOTE_OUTPUT_NAME}"
UNFILTERED_H5MU="${OUT_DIR}/mudata/star_chromap_unfiltered_multiome.h5mu"
FILTERED_H5MU="${OUT_DIR}/mudata/star_chromap_filtered_multiome.h5mu"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

ensure_mudata_python() {
  if python3 - <<'PY' >/dev/null 2>&1
import mudata
PY
  then
    echo python3
    return 0
  fi
  local venv="${OUT_DIR}/mudata_venv"
  if [[ ! -x "${venv}/bin/python" ]]; then
    python3 -m venv --system-site-packages "${venv}"
  fi
  "${venv}/bin/python" - <<'PY' >/dev/null 2>&1 || "${venv}/bin/python" -m pip install --quiet mudata
import mudata
PY
  echo "${venv}/bin/python"
}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  log "Clean rebuilding STAR with Chromap support"
  make -C "${REPO_ROOT}/core/legacy/source" clean > "${OUT_DIR}/logs/build_clean.log" 2>&1
  make -C "${REPO_ROOT}/core/legacy/source" -j8 STAR WITH_CHROMAP=1 > "${OUT_DIR}/logs/build_star_with_chromap.log" 2>&1
  if [[ "${USE_NATIVE_ATAC_MEX}" == "1" ]]; then
    log "Building native ATAC peak-MEX tool"
    make -C "${REPO_ROOT}/core/features/libchromap_contract" star_multiome_atac_peak_mex > "${OUT_DIR}/logs/build_star_multiome_atac_peak_mex.log" 2>&1
  fi
fi

if [[ "${USE_NATIVE_ATAC_MEX}" == "1" && ! -x "${BUILD_ATAC_MEX_NATIVE}" ]]; then
  log "Building native ATAC peak-MEX tool"
  make -C "${REPO_ROOT}/core/features/libchromap_contract" star_multiome_atac_peak_mex > "${OUT_DIR}/logs/build_star_multiome_atac_peak_mex.log" 2>&1
fi

if [[ "${USE_NATIVE_ATAC_BARCODE}" == "1" ]]; then
  if [[ -z "${ATAC_BARCODE_READ_FORMAT}" ]]; then
    if ! [[ "${ATAC_BARCODE_START}" =~ ^[0-9]+$ && "${ATAC_BARCODE_LENGTH}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: --atac-barcode-start/--atac-barcode-length must be positive integers" >&2
      exit 1
    fi
    if (( ATAC_BARCODE_START < 1 || ATAC_BARCODE_LENGTH < 1 )); then
      echo "ERROR: --atac-barcode-start/--atac-barcode-length must be positive integers" >&2
      exit 1
    fi
    atac_zero_start=$((ATAC_BARCODE_START - 1))
    atac_zero_end=$((atac_zero_start + ATAC_BARCODE_LENGTH - 1))
    ATAC_BARCODE_READ_FORMAT="bc:${atac_zero_start}:${atac_zero_end}:-"
  fi
  log "Using native Chromap ATAC barcode read format: ${ATAC_BARCODE_READ_FORMAT}"
else
  ATAC_BC_FOR_CHROMAP="${ATAC_BC_NORM}"
fi

for numeric_name in THREADS CHROMAP_THREADS CHROMAP_LOW_MEM CHROMAP_LOW_MEM_RAM CHROMAP_MACS3_FRAG_LOW_MEM; do
  if ! [[ "${!numeric_name}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: ${numeric_name} must be a non-negative integer" >&2
    exit 1
  fi
done

if [[ "${USE_NATIVE_ATAC_BARCODE}" != "1" && ( "${FORCE_ATAC_BARCODE}" == "1" || ! -f "${ATAC_BC_NORM}" ) ]]; then
  log "Normalizing ATAC barcode FASTQ"
  normalize_args=(
    python3 "${NORMALIZE_ATAC_BC}"
    --input "${ATAC_BARCODE}"
    --output "${ATAC_BC_NORM}"
    --start "${ATAC_BARCODE_START}"
    --length "${ATAC_BARCODE_LENGTH}"
    --reverse-complement
  )
  [[ "${FORCE_ATAC_BARCODE}" == "1" ]] && normalize_args+=(--force)
  "${normalize_args[@]}" > "${OUT_DIR}/logs/normalize_atac_barcode.log"
elif [[ "${USE_NATIVE_ATAC_BARCODE}" != "1" ]]; then
  log "Reusing normalized ATAC barcode FASTQ"
fi

if [[ "${FORCE}" == "1" ]]; then
  rm -rf "${RUN_DIR}" "${SAMPLE_DIR}/tmp"
fi
mkdir -p "${RUN_DIR}"

STAR_COMMAND="${OUT_DIR}/RUN_STAR_MULTIOME.sh"
solo_inline_hash_mode="no"
[[ "${SOLO_INLINE_HASH}" == "1" ]] && solo_inline_hash_mode="yes"
chromap_read_format_block=""
if [[ "${USE_NATIVE_ATAC_BARCODE}" == "1" ]]; then
  printf -v chromap_read_format_block '  --chromapAtacReadFormat "%s" \\\n' "${ATAC_BARCODE_READ_FORMAT}"
fi
cat > "${STAR_COMMAND}" <<EOF
#!/usr/bin/env bash
set -euo pipefail

"${STAR_BIN}" \\
  --runThreadN "${THREADS}" \\
  --genomeDir "${GENOME_DIR}" \\
  --readFilesIn "${GEX_R2}" "${GEX_R1}" \\
  --readFilesCommand zcat \\
  --outFileNamePrefix "${RUN_DIR}/" \\
  --outTmpDir "${SAMPLE_DIR}/tmp" \\
  --outSAMtype BAM Unsorted \\
  --outSAMattributes NH HI AS nM NM GX GN \\
  --emitNoYBAM yes \\
  --emitYNoYFastq yes \\
  --emitYNoYFastqCompression gz \\
  --clipAdapterType CellRanger4 \\
  --clip3pPolyG yes \\
  --alignEndsType Local \\
  --chimSegmentMin 1000000 \\
  --soloType CB_UMI_Simple \\
  --soloCBstart 1 \\
  --soloCBlen 16 \\
  --soloUMIstart 17 \\
  --soloUMIlen 12 \\
  --soloBarcodeReadLength 0 \\
  --soloCBwhitelist "${GEX_WHITELIST}" \\
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
  --soloUMIfiltering MultiGeneUMI_CR \\
  --soloUMIdedup 1MM_CR \\
  --soloMultiMappers Unique \\
  --soloCellFilter EmptyDrops_CR \\
  --soloCbUbRequireTogether no \\
  --soloStrand "${SOLO_STRAND}" \\
  --soloFeatures GeneFull Velocyto \\
  --soloCrGexFeature genefull \\
  --soloCrMultimapRescue yes \\
  --soloInlineHashMode "${solo_inline_hash_mode}" \\
  --chromapAtacEnable 1 \\
  --chromapAtacStartMode "${CHROMAP_START_MODE}" \\
  --chromapAtacReferenceFasta "${CHROMAP_REF}" \\
  --chromapAtacIndex "${CHROMAP_INDEX}" \\
  --chromapAtacRead1 "${ATAC_R1}" \\
  --chromapAtacRead2 "${ATAC_R2}" \\
  --chromapAtacBarcode "${ATAC_BC_FOR_CHROMAP}" \\
${chromap_read_format_block}  --chromapAtacBarcodeWhitelist "${ATAC_WHITELIST}" \\
  --chromapAtacBarcodeTranslate "${ATAC_TO_GEX}" \\
  --chromapAtacBarcodeTranslateFromFirst 1 \\
  --chromapAtacOutputFormat BAM \\
  --chromapAtacOutputFragments "${ATAC_BAM}" \\
  --chromapAtacSecondaryFragments "${ATAC_SIDECAR}" \\
  --chromapAtacSortBam 1 \\
  --chromapAtacSummary "${RUN_DIR}/chromap_summary.csv" \\
  --chromapAtacThreads "${CHROMAP_THREADS}" \\
  --chromapAtacLowMem "${CHROMAP_LOW_MEM}" \\
  --chromapAtacLowMemRam "${CHROMAP_LOW_MEM_RAM}" \\
  --chromapAtacMacs3FragLowMem "${CHROMAP_MACS3_FRAG_LOW_MEM}" \\
  --chromapAtacTempDir "${SAMPLE_DIR}/chromap_tmp" \\
  --chromapAtacTn5ShiftMode classical
EOF
chmod +x "${STAR_COMMAND}"

if [[ "${FORCE}" == "1" || ! -f "${RUN_DIR}/Log.final.out" || ! -f "${ATAC_BAM}" || ! -f "${ATAC_SIDECAR}" || ! -f "${ATAC_SIDECAR}.chroms.tsv" ]]; then
  log "Running STAR-suite GEX + Chromap ATAC"
  rm -rf "${SAMPLE_DIR}/tmp" "${SAMPLE_DIR}/chromap_tmp"
  mkdir -p "${SAMPLE_DIR}/chromap_tmp"
  bash "${STAR_COMMAND}" > "${OUT_DIR}/logs/star_multiome.stdout.log" 2> "${OUT_DIR}/logs/star_multiome.stderr.log"
else
  log "Reusing STAR/Chromap run"
fi

if [[ "${FORCE}" == "1" \
  || ! -f "${RUN_DIR}/outs/gene_full_feature_bc_matrix_manifest.json" \
  || ! -f "${RUN_DIR}/outs/raw_feature_bc_matrix/matrix.mtx.gz" \
  || ! -f "${RUN_DIR}/outs/filtered_feature_bc_matrix/matrix.mtx.gz" ]]; then
  log "Packaging GeneFull MEX"
  python3 "${PACKAGE_GENEFULL}" --run-dir "${RUN_DIR}" > "${OUT_DIR}/logs/package_star_genefull_mex.log"
else
  log "Reusing packaged GeneFull MEX"
fi

if [[ ! -f "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json" \
  || ! -f "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix/unspliced.mtx.gz" \
  || ! -f "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/unspliced.mtx.gz" ]]; then
  if [[ "${ALLOW_LEGACY_PREPARE_VELOCYTO}" == "1" ]]; then
    log "WARNING: using legacy prepare_velocyto_mex.py fallback"
    python3 "${PREPARE_VELOCYTO}" --run-dir "${RUN_DIR}" > "${OUT_DIR}/logs/prepare_velocyto_mex.log"
  else
    echo "ERROR: Native STAR Velocyto MEX outputs are missing under ${RUN_DIR}/outs" >&2
    exit 1
  fi
else
  log "Native Velocyto MEX outputs present"
fi

if [[ "${FORCE}" == "1" || ! -f "${ATAC_MEX}/matrix.mtx.gz" || ! -f "${ATAC_METRICS}" ]]; then
  log "Building local ATAC peak MEX from Chromap binary sidecar"
  rm -rf "${ATAC_MEX}"
  mkdir -p "${SAMPLE_DIR}/chromap_tmp"
  if [[ "${USE_NATIVE_ATAC_MEX}" == "1" ]]; then
    atac_mex_args=(
      "${BUILD_ATAC_MEX_NATIVE}"
      --sidecar "${ATAC_SIDECAR}"
      --barcode-translate "${ATAC_TO_GEX}"
      --barcode-translate-from-first
      --call-peaks-from-sidecar
      --peaks "${ATAC_PEAKS}" \
      --summits-out "${ATAC_SUMMITS}" \
      --out-dir "${ATAC_MEX}" \
      --metrics-tsv "${ATAC_METRICS}" \
      --threads "${CHROMAP_THREADS}" \
      --temp-dir "${SAMPLE_DIR}/chromap_tmp"
    )
    [[ "${FORCE}" == "1" ]] && atac_mex_args+=(--force)
    "${atac_mex_args[@]}" | tee "${OUT_DIR}/logs/build_atac_peak_matrix.log"
  else
    echo "ERROR: --python-atac-mex is not compatible with the binary-sidecar production boundary" >&2
    exit 1
  fi
else
  log "Reusing ATAC peak MEX"
fi

if [[ "${STOP_AFTER_LOCAL_MEX}" == "1" ]]; then
  {
    printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'run_dir=%s\n' "${RUN_DIR}"
    printf 'rna_raw_mex=%s\n' "${RUN_DIR}/outs/raw_feature_bc_matrix"
    printf 'rna_filtered_mex=%s\n' "${RUN_DIR}/outs/filtered_feature_bc_matrix"
    printf 'velocyto_raw_mex=%s\n' "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix"
    printf 'atac_bam=%s\n' "${ATAC_BAM}"
    printf 'atac_sidecar=%s\n' "${ATAC_SIDECAR}"
    printf 'atac_peaks=%s\n' "${ATAC_PEAKS}"
    printf 'atac_summits=%s\n' "${ATAC_SUMMITS}"
    printf 'atac_mex=%s\n' "${ATAC_MEX}"
    printf 'atac_metrics=%s\n' "${ATAC_METRICS}"
    printf 'post_mex_boundary=local_star_chromap_matrices_ready\n'
  } > "${OUT_DIR}/LOCAL_MEX_READY.txt"
  log "PASS: local STAR/Chromap MEX boundary complete"
  echo "Output dir: ${OUT_DIR}"
  exit 0
fi

POST_MEX_BUILT_MUDATA="0"
if [[ "${FORCE}" != "1" \
  && -f "${UNFILTERED_H5MU}" \
  && -f "${FILTERED_H5MU}" ]]; then
  log "Reusing post-MEX MuData outputs"
  POST_MEX_BUILT_MUDATA="1"
elif [[ -n "${REMOTE_HOST}" && -n "${REMOTE_ROOT}" ]]; then
  remote_args=(
    "${REMOTE_POST_MEX}"
    --sample-dir "${SAMPLE_DIR}"
    --remote-host "${REMOTE_HOST}"
    --remote-root "${REMOTE_ROOT}"
    --output-name "${REMOTE_OUTPUT_NAME}"
    --run-cellbender
    --adaptive-filter
    --local-log "${OUT_DIR}/logs/remote_post_mex.log"
  )
  [[ -n "${CELLBENDER_CPU_CORES}" ]] && remote_args+=(--cellbender-cpu-cores "${CELLBENDER_CPU_CORES}")
  [[ "${CELLBENDER_GPU}" == "1" ]] && remote_args+=(--cellbender-gpu)
  [[ "${NO_SYNC_IMAGES}" == "1" ]] && remote_args+=(--no-sync-images)
  [[ "${KEEP_REMOTE}" == "1" ]] && remote_args+=(--keep-remote)
  log "Running remote post-MEX downstream, CellBender, and MuData"
  "${remote_args[@]}" | tee "${OUT_DIR}/logs/remote_post_mex.wrapper.log"
  POST_MEX_BUILT_MUDATA="1"
elif [[ "${ALLOW_LOCAL_DOWNSTREAM}" == "1" ]]; then
  DOWNSTREAM_DIR="${SAMPLE_DIR}/downstream_genefull_velocyto"
  log "Running local RNA downstream without CellBender"
  "${LOCAL_DOWNSTREAM}" \
    --run-dir "${RUN_DIR}" \
    --output-dir "${DOWNSTREAM_DIR}" \
    --adaptive-filter | tee "${OUT_DIR}/logs/local_downstream.log"
else
  echo "ERROR: provide --remote-host/--remote-root for post-MEX CellBender/MuData, or pass --allow-local-downstream" >&2
  exit 1
fi

RNA_UNFILTERED="${DOWNSTREAM_DIR}/final_counts.h5ad"
[[ -f "${RNA_UNFILTERED}" ]] || RNA_UNFILTERED="${DOWNSTREAM_DIR}/unfiltered_counts.h5ad"
RNA_FILTERED="${DOWNSTREAM_DIR}/filtered_counts.h5ad"
[[ -f "${RNA_UNFILTERED}" ]] || { echo "ERROR: missing RNA h5ad in ${DOWNSTREAM_DIR}" >&2; exit 1; }
[[ -f "${RNA_FILTERED}" ]] || { echo "ERROR: missing filtered RNA h5ad ${RNA_FILTERED}" >&2; exit 1; }

MUDATA_PYTHON="$(ensure_mudata_python)"

if [[ "${POST_MEX_BUILT_MUDATA}" != "1" ]]; then
  log "Building MuData outputs locally"
  "${MUDATA_PYTHON}" "${BUILD_MUDATA}" \
    --rna-h5ad "${RNA_UNFILTERED}" \
    --atac-mex-dir "${ATAC_MEX}" \
    --per-barcode-metrics "${ATAC_METRICS}" \
    --metrics-barcode-column barcode \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_h5ad_chromap_atac \
    --rna-source "${RNA_UNFILTERED}" \
    --atac-source "${ATAC_MEX}" \
    --fragments-source "${ATAC_SIDECAR}" \
    --peaks-source "${ATAC_PEAKS}" \
    --evidence-source "${ATAC_METRICS}" \
    --y-removal-enabled true \
    --output-h5mu "${UNFILTERED_H5MU}" | tee "${OUT_DIR}/logs/build_unfiltered_h5mu.log"

  "${MUDATA_PYTHON}" "${BUILD_MUDATA}" \
    --rna-h5ad "${RNA_FILTERED}" \
    --atac-mex-dir "${ATAC_MEX}" \
    --per-barcode-metrics "${ATAC_METRICS}" \
    --metrics-barcode-column barcode \
    --all-barcodes-are-cells \
    --require-rna-velocyto-layers \
    --cell-call-source star_downstream_filtered_h5ad_chromap_atac \
    --rna-source "${RNA_FILTERED}" \
    --atac-source "${ATAC_MEX}" \
    --fragments-source "${ATAC_SIDECAR}" \
    --peaks-source "${ATAC_PEAKS}" \
    --evidence-source "${ATAC_METRICS}" \
    --y-removal-enabled true \
    --output-h5mu "${FILTERED_H5MU}" | tee "${OUT_DIR}/logs/build_filtered_h5mu.log"
else
  log "Using remote-built MuData outputs"
fi

"${MUDATA_PYTHON}" - <<PY | tee "${OUT_DIR}/logs/validate_h5mu.log"
import mudata as md
for path in ["${UNFILTERED_H5MU}", "${FILTERED_H5MU}"]:
    m = md.read_h5mu(path)
    rna = m.mod["rna"]
    atac = m.mod["atac"]
    required_layers = {"counts", "spliced", "unspliced", "ambiguous"}
    missing = sorted(required_layers - set(rna.layers))
    if missing:
        raise SystemExit(f"{path}: missing RNA layers {missing}")
    if "counts" not in atac.layers:
        raise SystemExit(f"{path}: missing ATAC counts layer")
    print(path)
    print(f"  obs={m.n_obs} rna_vars={rna.n_vars} atac_vars={atac.n_vars}")
    print(f"  rna_layers={sorted(rna.layers.keys())}")
    print(f"  atac_layers={sorted(atac.layers.keys())}")
    print(f"  y_removal={m.uns['multiome'].get('y_removal_enabled')}")
PY

{
  printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'gex_r1=%s\n' "${GEX_R1}"
  printf 'gex_r2=%s\n' "${GEX_R2}"
  printf 'atac_r1=%s\n' "${ATAC_R1}"
  printf 'atac_barcode_raw=%s\n' "${ATAC_BARCODE}"
  printf 'atac_barcode_for_chromap=%s\n' "${ATAC_BC_FOR_CHROMAP}"
  printf 'atac_barcode_normalized=%s\n' "$([[ "${USE_NATIVE_ATAC_BARCODE}" == "1" ]] && echo "-" || echo "${ATAC_BC_NORM}")"
  printf 'atac_r2=%s\n' "${ATAC_R2}"
  printf 'genome_dir=%s\n' "${GENOME_DIR}"
  printf 'gex_whitelist=%s\n' "${GEX_WHITELIST}"
  printf 'chromap_ref=%s\n' "${CHROMAP_REF}"
  printf 'chromap_index=%s\n' "${CHROMAP_INDEX}"
  printf 'chromap_threads=%s\n' "${CHROMAP_THREADS}"
  printf 'chromap_low_mem=%s\n' "${CHROMAP_LOW_MEM}"
  printf 'chromap_low_mem_ram=%s\n' "${CHROMAP_LOW_MEM_RAM}"
  printf 'chromap_macs3_frag_low_mem=%s\n' "${CHROMAP_MACS3_FRAG_LOW_MEM}"
  printf 'atac_whitelist=%s\n' "${ATAC_WHITELIST}"
  printf 'atac_to_gex=%s\n' "${ATAC_TO_GEX}"
  printf 'solo_strand=%s\n' "${SOLO_STRAND}"
  printf 'solo_inline_hash=%s\n' "${SOLO_INLINE_HASH}"
  printf 'native_atac_barcode=%s\n' "${USE_NATIVE_ATAC_BARCODE}"
  printf 'native_atac_mex=%s\n' "${USE_NATIVE_ATAC_MEX}"
  printf 'atac_barcode_window=%s:%s_rc\n' "${ATAC_BARCODE_START}" "${ATAC_BARCODE_LENGTH}"
  printf 'chromap_atac_read_format=%s\n' "$([[ "${USE_NATIVE_ATAC_BARCODE}" == "1" ]] && echo "${ATAC_BARCODE_READ_FORMAT}" || echo "-")"
  printf 'run_dir=%s\n' "${RUN_DIR}"
  printf 'downstream_dir=%s\n' "${DOWNSTREAM_DIR}"
  printf 'atac_bam=%s\n' "${ATAC_BAM}"
  printf 'atac_sidecar=%s\n' "${ATAC_SIDECAR}"
  printf 'atac_peaks=%s\n' "${ATAC_PEAKS}"
  printf 'atac_summits=%s\n' "${ATAC_SUMMITS}"
  printf 'atac_mex=%s\n' "${ATAC_MEX}"
  printf 'atac_metrics=%s\n' "${ATAC_METRICS}"
  printf 'unfiltered_h5mu=%s\n' "${UNFILTERED_H5MU}"
  printf 'filtered_h5mu=%s\n' "${FILTERED_H5MU}"
} > "${OUT_DIR}/RUN_MANIFEST.txt"

log "PASS: STAR multiome lane smoke complete"
echo "Output dir: ${OUT_DIR}"
