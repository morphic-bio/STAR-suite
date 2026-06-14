#!/usr/bin/env bash
# Downsampled CAT-ATAC trimodal STAR smoke: GEX Solo + Chromap ATAC + guide pf-multi.
# This is a validation/debug harness, not a production benchmark.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
OUT_DIR="${CATATAC_TRIMODAL_SMOKE_OUT:-${REPO_ROOT}/tests/catatac_trimodal_downsample_smoke_output}"
MAX_READS="${CATATAC_TRIMODAL_MAX_READS:-100000}"
THREADS="${CATATAC_TRIMODAL_THREADS:-16}"
INLINE_ATAC_PEAK_MEX="${CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX:-no}"
REQUIRE_NO_STAR="${CATATAC_TRIMODAL_REQUIRE_NO_STAR:-1}"

GENOME="/mnt/pikachu/autoindex_98_32/pe_index"
GEX_R1="/mnt/pikachu/catatac_gse288996/fastq/GEX/SRR32265752_1.fastq.gz"
GEX_R2="/mnt/pikachu/catatac_gse288996/fastq/GEX/SRR32265752_2.fastq.gz"
ATAC_R1_SRC="/mnt/pikachu/catatac_gse288996/fastq/ATAC/SRR32265760_1.fastq.gz"
ATAC_R2_SRC="/mnt/pikachu/catatac_gse288996/fastq/ATAC/SRR32265760_2.fastq.gz"
ATAC_R3_SRC="/mnt/pikachu/catatac_gse288996/fastq/ATAC/SRR32265760_3.fastq.gz"
GUIDE_DIR="${CATATAC_TRIMODAL_GUIDE_DIR:-/mnt/pikachu/catatac_gse288996/guide_redump/fixture}"
STAGED_ROOT="${CATATAC_TRIMODAL_STAGED_ROOT:-/mnt/pikachu/catatac_gse288996/full_bench/catatac_trimodal_staged}"
STAGED_ATAC="${CATATAC_TRIMODAL_STAGED_ATAC:-${STAGED_ROOT}/atac_${MAX_READS}}"
FEATURE_REF="${REPO_ROOT}/core/features/process_features/feature_lists/catatac_crispri_guide_capture.csv"
GEX_WL="/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt"
ATAC_WL="/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt"
ATAC2GEX="/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv"
CHROMAP_FASTA="/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
CHROMAP_IDX="/mnt/pikachu/catatac_gse288996/refs/GRCh38-arc.chromap.idx"
LIBRARY_ID="${CATATAC_GUIDE_LIBRARY_ID:-catatac_guide_fixture}"
CS1="CAAGTTGATAACGGACTAGCC"
CS2="CAAGTTGTAAACGGACTAGCC"

case "${INLINE_ATAC_PEAK_MEX}" in
  yes|no) ;;
  1|true|TRUE|True) INLINE_ATAC_PEAK_MEX="yes" ;;
  0|false|FALSE|False) INLINE_ATAC_PEAK_MEX="no" ;;
  *) echo "ERROR: CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX must be yes/no" >&2; exit 1 ;;
esac

if [[ ! "${MAX_READS}" =~ ^[0-9]+$ ]] || [[ "${MAX_READS}" -le 0 ]]; then
  echo "ERROR: CATATAC_TRIMODAL_MAX_READS must be a positive integer for smoke runs" >&2
  exit 1
fi

resolve_guide_fastq() {
  local read_idx="$1"
  local read_token="R${read_idx}"
  local path
  for path in \
    "${GUIDE_DIR}/guide_${read_token}.fastq.gz" \
    "${GUIDE_DIR}/SRR32265756_${read_idx}.fastq.gz" \
    "${GUIDE_DIR}"/*_"${read_idx}".fastq.gz \
    "${GUIDE_DIR}"/*_"${read_token}"_*.fastq.gz; do
    if [[ -f "${path}" ]]; then
      printf '%s\n' "${path}"
      return 0
    fi
  done
  return 1
}

GUIDE_R1="$(resolve_guide_fastq 1)" || {
  echo "ERROR: could not resolve guide R1 in ${GUIDE_DIR}" >&2
  exit 1
}
GUIDE_R2="$(resolve_guide_fastq 2)" || {
  echo "ERROR: could not resolve guide R2 in ${GUIDE_DIR}" >&2
  exit 1
}
GUIDE_R3="$(resolve_guide_fastq 3)" || {
  echo "ERROR: could not resolve guide R3 in ${GUIDE_DIR}" >&2
  exit 1
}

mkdir -p "${OUT_DIR}"
LOCK_FILE="${OUT_DIR}/.trimodal_smoke.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  echo "FATAL: another trimodal downsample smoke holds ${LOCK_FILE}" >&2
  exit 1
fi

for path in "${STAR_BIN}" "${GENOME}/Genome" "${GEX_R1}" "${GEX_R2}" \
  "${ATAC_R1_SRC}" "${ATAC_R2_SRC}" "${ATAC_R3_SRC}" "${GUIDE_R1}" "${GUIDE_R2}" "${GUIDE_R3}" \
  "${FEATURE_REF}" "${GEX_WL}" "${ATAC_WL}" "${ATAC2GEX}" "${CHROMAP_FASTA}" "${CHROMAP_IDX}"; do
  [[ -r "${path}" ]] || { echo "ERROR: missing ${path}" >&2; exit 1; }
done

if [[ "${REQUIRE_NO_STAR}" != "0" ]] && pgrep -x STAR >/dev/null 2>&1; then
  echo "ERROR: competing STAR job already running (serialize benchmark runs)" >&2
  pgrep -x STAR -a | head -3 >&2
  exit 1
fi

stage_atac_downsample() {
  local out_dir="$1" limit="$2"
  mkdir -p "${out_dir}"
  local manifest="${out_dir}/MANIFEST.txt"
  local -a srcs=("${ATAC_R1_SRC}" "${ATAC_R2_SRC}" "${ATAC_R3_SRC}")
  local -a dsts=(
    "${out_dir}/ATAC_R1.fastq.gz"
    "${out_dir}/ATAC_R2.fastq.gz"
    "${out_dir}/ATAC_R3.fastq.gz"
  )
  local ok=1
  for dst in "${dsts[@]}"; do
    [[ -s "${dst}" ]] || ok=0
  done
  if [[ ! -s "${manifest}" ]] || ! grep -qx "read_count=${limit}" "${manifest}"; then
    ok=0
  fi
  if [[ "${ok}" -eq 1 ]]; then
    echo "Using cached staged ATAC: ${out_dir}"
    return 0
  fi
  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"
  echo "Staging ${limit} ATAC read triplets -> ${out_dir}"
  python3 - <<PY "${limit}" "${srcs[0]}" "${srcs[1]}" "${srcs[2]}" "${dsts[0]}" "${dsts[1]}" "${dsts[2]}"
import gzip, sys
limit = int(sys.argv[1])
triplets = list(zip(sys.argv[2:5], sys.argv[5:8]))
for src, dst in triplets:
    n = 0
    with gzip.open(src, "rt") as inp, gzip.open(dst, "wt") as out:
        while n < limit:
            header = inp.readline()
            if not header:
                break
            out.write(header)
            for _ in range(3):
                line = inp.readline()
                if not line:
                    break
                out.write(line)
            n += 1
    print(f"staged {n} reads -> {dst}")
PY
  {
    printf 'read_count=%s\n' "${limit}"
    printf 'source_r1=%s\n' "${ATAC_R1_SRC}"
    printf 'source_r2=%s\n' "${ATAC_R2_SRC}"
    printf 'source_r3=%s\n' "${ATAC_R3_SRC}"
  } > "${manifest}"
}

# Preserve flock inode; clear prior run artifacts only.
find "${OUT_DIR}" -mindepth 1 -maxdepth 1 ! -name '.trimodal_smoke.lock' -exec rm -rf {} +
mkdir -p "${OUT_DIR}/logs"
exec > >(tee -a "${OUT_DIR}/logs/smoke.log") 2>&1
STAR_RUN="${OUT_DIR}/star_run"
STAR_TMP="${OUT_DIR}/star_tmp"
CHROMAP_TMP="${OUT_DIR}/chromap_tmp"
PF_CONFIG="${OUT_DIR}/pf_multi_config.csv"

stage_atac_downsample "${STAGED_ATAC}" "${MAX_READS}"
ATAC_R1="${STAGED_ATAC}/ATAC_R1.fastq.gz"
ATAC_R2="${STAGED_ATAC}/ATAC_R2.fastq.gz"
ATAC_R3="${STAGED_ATAC}/ATAC_R3.fastq.gz"

cat > "${PF_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_layout,star_barcode_read,star_barcode_format,star_umi_read,star_umi_start,star_umi_length,star_feature_read,star_capture_read,star_capture_sequences,star_capture_max_hamming,star_barcode_output_map,star_feature_search_mode,star_max_hamming
/mnt/pikachu/catatac_gse288996/fastq/GEX,DMSO1,Gene Expression,Gene Expression,TRU,gex_dmso1,,,,
${GUIDE_DIR},DMSO1,CRISPR Guide Capture,CRISPR Guide Capture,,${LIBRARY_ID},${FEATURE_REF},${ATAC_WL},catatac_guide,R2,bc:8:23:-,R1,0,12,R3,R1,${CS1}|${CS2},0,${ATAC2GEX},free,1
EOF

cat > "${OUT_DIR}/PREFLIGHT.txt" <<EOF
trimodal_downsample_smoke
max_reads=${MAX_READS}
threads=${THREADS}
inline_atac_peak_mex=${INLINE_ATAC_PEAK_MEX}
guide_dir=${GUIDE_DIR}
guide_r1=${GUIDE_R1}
guide_r2=${GUIDE_R2}
guide_r3=${GUIDE_R3}
staged_atac=${STAGED_ATAC}
star_binary=${STAR_BIN}
EOF

echo "=== CAT-ATAC trimodal downsample smoke ==="
echo "  OUT_DIR=${OUT_DIR}"
echo "  MAX_READS=${MAX_READS}"
echo "  THREADS=${THREADS}"
echo "  INLINE_ATAC_PEAK_MEX=${INLINE_ATAC_PEAK_MEX}"
echo "  GUIDE_DIR=${GUIDE_DIR}"

rm -rf "${STAR_TMP}" "${CHROMAP_TMP}"
mkdir -p "${STAR_RUN}" "${CHROMAP_TMP}"

STAR_CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME}" \
  --readFilesIn "${GEX_R2}" "${GEX_R1}" \
  --readFilesCommand zcat \
  --readMapNumber "${MAX_READS}" \
  --outFileNamePrefix "${STAR_RUN}/" \
  --outTmpDir "${STAR_TMP}" \
  --outSAMtype None \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${GEX_WL}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Forward \
  --soloFeatures GeneFull \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --soloInlineHashMode no \
  --defaultCrCompat yes \
  --pfMultiConfig "${PF_CONFIG}" \
  --crWhitelist "${GEX_WL}" \
  --crFeatureRef "${FEATURE_REF}" \
  --dynamicThreadInterface 1 \
  --dynamicThreadTelemetry 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  --chromapAtacEnable 1 \
  --chromapAtacStartMode concurrent \
  --chromapAtacReferenceFasta "${CHROMAP_FASTA}" \
  --chromapAtacIndex "${CHROMAP_IDX}" \
  --chromapAtacRead1 "${ATAC_R1}" \
  --chromapAtacRead2 "${ATAC_R3}" \
  --chromapAtacBarcode "${ATAC_R2}" \
  --chromapAtacReadFormat "bc:8:23:-" \
  --chromapAtacBarcodeWhitelist "${ATAC_WL}" \
  --chromapAtacBarcodeTranslate "${ATAC2GEX}" \
  --chromapAtacBarcodeTranslateFromFirst 1 \
  --chromapAtacOutputFormat BAM \
  --chromapAtacOutputFragments "${STAR_RUN}/atac_possorted.bam" \
  --chromapAtacSecondaryFragments "${STAR_RUN}/atac_fragments.bin" \
  --chromapAtacSortBam 1 \
  --chromapAtacSummary "${STAR_RUN}/chromap_summary.csv" \
  --chromapAtacThreads "${THREADS}" \
  --chromapAtacLowMem 1 \
  --chromapAtacLowMemRam 0 \
  --chromapAtacMacs3FragLowMem 1 \
  --chromapAtacTempDir "${CHROMAP_TMP}" \
  --chromapAtacTn5ShiftMode classical
)

if [[ "${INLINE_ATAC_PEAK_MEX}" == "yes" ]]; then
  STAR_CMD+=(
    --multiomeAtacPeakMexInline yes
    --multiomeAtacPeakMexOutDir "${STAR_RUN}/atac/peak_mex"
    --multiomeAtacPeakNarrowPeak "${STAR_RUN}/atac/atac_peaks.narrowPeak"
    --multiomeAtacPeakSummits "${STAR_RUN}/atac/atac_summits.bed"
    --multiomeAtacPeakMetricsTsv "${STAR_RUN}/atac/atac_metrics.tsv"
    --multiomeAtacPeakThreads "${THREADS}"
  )
fi

{
  printf '#!/usr/bin/env bash\n'
  printf 'set -euo pipefail\n'
  printf '/usr/bin/time -v -o %q \\\n' "${OUT_DIR}/logs/time_trimodal.log"
  printf '  %q' "${STAR_CMD[0]}"
  for ((i = 1; i < ${#STAR_CMD[@]}; ++i)); do
    printf ' \\\n  %q' "${STAR_CMD[$i]}"
  done
  printf '\n'
} > "${OUT_DIR}/RUN_STAR_TRIMODAL_SMOKE.sh"
chmod +x "${OUT_DIR}/RUN_STAR_TRIMODAL_SMOKE.sh"

set +e
/usr/bin/time -v -o "${OUT_DIR}/logs/time_trimodal.log" \
  "${STAR_CMD[@]}" \
  2>&1 | tee "${OUT_DIR}/logs/star_stdout.log"

STAR_EXIT=${PIPESTATUS[0]}
set -e
echo "${STAR_EXIT}" > "${OUT_DIR}/logs/star_exit_code.txt"

python3 "${SCRIPT_DIR}/catatac_trimodal_downsample_verify.py" \
  --star-run "${STAR_RUN}" \
  --gex-whitelist "${GEX_WL}" \
  --guide-library-id "${LIBRARY_ID}" \
  --expect-atac-peak-mex "${INLINE_ATAC_PEAK_MEX}" \
  --report "${OUT_DIR}/trimodal_verify.json" \
  || VERIFY_EXIT=$?

VERIFY_EXIT="${VERIFY_EXIT:-0}"
if [[ "${STAR_EXIT}" -ne 0 ]]; then
  echo "FATAL: STAR exit ${STAR_EXIT}" >&2
  exit "${STAR_EXIT}"
fi
if [[ "${VERIFY_EXIT}" -ne 0 ]]; then
  exit "${VERIFY_EXIT}"
fi

echo "CAT-ATAC trimodal downsample smoke PASSED: ${OUT_DIR}"
