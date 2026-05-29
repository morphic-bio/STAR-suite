#!/usr/bin/env bash
set -euo pipefail

# SLAM parity smoke:
#   1. Run STAR-SLAM on the host-local 100K fixture.
#   2. Gate correlation parity against the noSNP GRAND-SLAM reference.
#   3. Replay the STAR dump through slam_requant and gate STAR-vs-requant parity.
#
# Missing external fixtures are treated as SKIP so this remains safe on hosts
# that do not have the private SLAM fixture cache.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -f "${SCRIPT_DIR}/external_fixtures_env.sh" ]]; then
  # shellcheck source=/dev/null
  source "${SCRIPT_DIR}/external_fixtures_env.sh"
fi

DEFAULT_REQ_BIN="${ROOT_DIR}/slam/tools/slam_requant/slam_requant"
if [[ ! -x "${DEFAULT_REQ_BIN}" && -x "${ROOT_DIR}/tools/slam_requant/slam_requant" ]]; then
  DEFAULT_REQ_BIN="${ROOT_DIR}/tools/slam_requant/slam_requant"
fi

OUT_BASE="${OUT_BASE:-/tmp/slam_parity_smoke_${STAMP}}"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
REQ_BIN="${REQ_BIN:-${DEFAULT_REQ_BIN}}"
FASTQ="${FASTQ:-${SLAM_FIXTURE_FASTQ:-${ROOT_DIR}/test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz}}"
STAR_INDEX="${STAR_INDEX:-${SLAM_FIXTURE_STAR_INDEX:-${ROOT_DIR}/test/fixtures/slam/ref/star_index}}"
SNPS_BED="${SNPS_BED:-${SLAM_FIXTURE_SNPS_BED:-${ROOT_DIR}/test/fixtures/slam/ref/snps.bed}}"
REF_TSV="${REF_TSV:-${SLAM_FIXTURE_NOSNP_REF_TSV:-${SLAM_FIXTURE_REF_TSV:-${ROOT_DIR}/test/fixtures/slam/expected/from_nosnp.tsv.gz}}}"

THREADS="${THREADS:-4}"
REF_CORR_MIN="${REF_CORR_MIN:-0.99}"
EXACT_MIN_PEARSON="${EXACT_MIN_PEARSON:-0.999999}"
EXACT_MAX_ABS_DELTA="${EXACT_MAX_ABS_DELTA:-1e-6}"
ALIGN_MIN_PEARSON="${ALIGN_MIN_PEARSON:-0.98}"
AUTO_BUILD_REQUANT="${AUTO_BUILD_REQUANT:-1}"
RUN_REFERENCE="${RUN_REFERENCE:-1}"
RUN_REQUANT="${RUN_REQUANT:-1}"
COMPARE_ONLY=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Runs STAR-SLAM parity gates on the small SLAM fixture.

Options:
  --outdir DIR             Output root (default: ${OUT_BASE})
  --star-bin PATH          STAR binary (default: ${STAR_BIN})
  --requant-bin PATH       slam_requant binary (default: ${REQ_BIN})
  --fastq PATH             Fixture FASTQ (default: ${FASTQ})
  --genome-dir DIR         STAR genomeDir (default: ${STAR_INDEX})
  --snp-bed PATH           SNP mask BED (default: ${SNPS_BED})
  --reference PATH         noSNP GRAND-SLAM reference TSV.GZ (default: ${REF_TSV})
  --threads N              STAR threads (default: ${THREADS})
  --ref-corr-min X         Minimum STAR-vs-GRAND-SLAM NTR Pearson (default: ${REF_CORR_MIN})
  --exact-min-pearson X    Minimum exact replay Pearson (default: ${EXACT_MIN_PEARSON})
  --exact-max-abs-delta X  Maximum exact replay per-gene delta (default: ${EXACT_MAX_ABS_DELTA})
  --align-min-pearson X    Minimum recomputed-alignment replay Pearson (default: ${ALIGN_MIN_PEARSON})
  --no-reference           Skip GRAND-SLAM reference correlation gate
  --no-requant             Skip slam_requant replay gates
  --compare-only           Reuse existing STAR/requant outputs in --outdir
EOF
}

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

skip() {
  log "SKIP: $*"
  exit 0
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

quote_cmd() {
  printf '%q ' "$@"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir|--out-base) OUT_BASE="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --requant-bin) REQ_BIN="$2"; shift 2 ;;
    --fastq) FASTQ="$2"; shift 2 ;;
    --genome-dir) STAR_INDEX="$2"; shift 2 ;;
    --snp-bed) SNPS_BED="$2"; shift 2 ;;
    --reference) REF_TSV="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --ref-corr-min) REF_CORR_MIN="$2"; shift 2 ;;
    --exact-min-pearson) EXACT_MIN_PEARSON="$2"; shift 2 ;;
    --exact-max-abs-delta) EXACT_MAX_ABS_DELTA="$2"; shift 2 ;;
    --align-min-pearson) ALIGN_MIN_PEARSON="$2"; shift 2 ;;
    --no-reference) RUN_REFERENCE=0; shift ;;
    --no-requant) RUN_REQUANT=0; shift ;;
    --compare-only) COMPARE_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "--threads must be a positive integer"
[[ -f "${SCRIPT_DIR}/slam/compare_fixture.py" ]] || die "Missing compare_fixture.py"
[[ -f "${SCRIPT_DIR}/slam/compare_star_outputs.py" ]] || die "Missing compare_star_outputs.py"

if [[ "${COMPARE_ONLY}" == "0" ]]; then
  [[ -x "${STAR_BIN}" ]] || skip "STAR binary not found: ${STAR_BIN}"
  [[ -f "${FASTQ}" ]] || skip "SLAM fixture FASTQ not found: ${FASTQ}"
  [[ -d "${STAR_INDEX}" ]] || skip "SLAM fixture STAR index not found: ${STAR_INDEX}"
  [[ -f "${SNPS_BED}" ]] || skip "SLAM SNP BED not found: ${SNPS_BED}"
  if [[ "${RUN_REFERENCE}" == "1" ]]; then
    [[ -f "${REF_TSV}" ]] || skip "SLAM noSNP reference not found: ${REF_TSV}"
  fi
fi

if [[ "${RUN_REQUANT}" == "1" && ! -x "${REQ_BIN}" && "${AUTO_BUILD_REQUANT}" == "1" ]]; then
  if [[ -f "${ROOT_DIR}/slam/tools/slam_requant/Makefile" ]]; then
    log "Building slam_requant"
    make -C "${ROOT_DIR}/slam/tools/slam_requant" -j"${THREADS}"
  fi
fi
if [[ "${RUN_REQUANT}" == "1" && ! -x "${REQ_BIN}" ]]; then
  skip "slam_requant binary not found: ${REQ_BIN}"
fi

mkdir -p "${OUT_BASE}/logs" "${OUT_BASE}/compare"

SLAM_PREFIX="${OUT_BASE}/star_slam_"
SLAM_OUT="${SLAM_PREFIX}SlamQuant.out"
DUMP_PATH="${OUT_BASE}/star_slam.dump.bin"
WEIGHT_PATH="${OUT_BASE}/star_slam.weights.bin"

write_command() {
  local path="$1"
  shift
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\n\n'
    quote_cmd "$@"
    printf '\n'
  } > "${path}"
  chmod +x "${path}"
}

if [[ "${COMPARE_ONLY}" == "0" ]]; then
  cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --genomeDir "${STAR_INDEX}"
    --readFilesIn "${FASTQ}"
    --readFilesCommand zcat
    --outFileNamePrefix "${SLAM_PREFIX}"
    --outSAMtype None
    --clip3pAdapterSeq AGATCGGAAGAG
    --clip3pAdapterMMp 0.1
    --slamQuantMode 1
    --slamSnpMaskIn "${SNPS_BED}"
    --slamGrandSlamOut 1
    --slamDumpBinary "${DUMP_PATH}"
    --slamDumpWeights "${WEIGHT_PATH}"
  )
  write_command "${OUT_BASE}/RUN_STAR_SLAM.sh" "${cmd[@]}"

  log "Running STAR-SLAM fixture parity source"
  /usr/bin/time -v -o "${OUT_BASE}/logs/star_slam.time.log" \
    "${cmd[@]}" \
    > "${OUT_BASE}/logs/star_slam.stdout.log" \
    2> "${OUT_BASE}/logs/star_slam.stderr.log"
fi

[[ -s "${SLAM_OUT}" ]] || die "Missing STAR-SLAM output: ${SLAM_OUT}"

if [[ "${RUN_REFERENCE}" == "1" ]]; then
  log "Comparing STAR-SLAM against noSNP GRAND-SLAM reference"
  python3 "${SCRIPT_DIR}/slam/compare_fixture.py" \
    --reference "${REF_TSV}" \
    --test "${SLAM_OUT}" \
    --correlation-only \
    --corr-min "${REF_CORR_MIN}" \
    > "${OUT_BASE}/compare/grandslam_nosnp_reference.log" \
    2>&1
fi

if [[ "${RUN_REQUANT}" == "1" ]]; then
  [[ -s "${DUMP_PATH}" ]] || die "Missing STAR-SLAM dump: ${DUMP_PATH}"
  [[ -s "${WEIGHT_PATH}" ]] || die "Missing STAR-SLAM weight sidecar: ${WEIGHT_PATH}"

  log "Replaying STAR-SLAM dump through slam_requant"
  "${REQ_BIN}" \
    --dump "${DUMP_PATH}" \
    --out "${OUT_BASE}/requant_dump_" \
    --slamSnpMaskIn "${SNPS_BED}" \
    > "${OUT_BASE}/logs/requant_dump.log" 2>&1

  "${REQ_BIN}" \
    --dump "${DUMP_PATH}" \
    --out "${OUT_BASE}/requant_align_" \
    --slamSnpMaskIn "${SNPS_BED}" \
    --slamWeightMode alignments \
    > "${OUT_BASE}/logs/requant_align.log" 2>&1

  "${REQ_BIN}" \
    --dump "${DUMP_PATH}" \
    --out "${OUT_BASE}/requant_weight_" \
    --slamSnpMaskIn "${SNPS_BED}" \
    --slamWeightFile "${WEIGHT_PATH}" \
    > "${OUT_BASE}/logs/requant_weight.log" 2>&1

  python3 "${SCRIPT_DIR}/slam/compare_star_outputs.py" \
    --reference "${SLAM_OUT}" \
    --test "${OUT_BASE}/requant_dump_SlamQuant.out" \
    --min-pearson "${EXACT_MIN_PEARSON}" \
    --max-abs-delta "${EXACT_MAX_ABS_DELTA}" \
    --require-same-genes \
    > "${OUT_BASE}/compare/requant_dump.log" 2>&1

  python3 "${SCRIPT_DIR}/slam/compare_star_outputs.py" \
    --reference "${SLAM_OUT}" \
    --test "${OUT_BASE}/requant_weight_SlamQuant.out" \
    --min-pearson "${EXACT_MIN_PEARSON}" \
    --max-abs-delta "${EXACT_MAX_ABS_DELTA}" \
    --require-same-genes \
    > "${OUT_BASE}/compare/requant_weight.log" 2>&1

  python3 "${SCRIPT_DIR}/slam/compare_star_outputs.py" \
    --reference "${SLAM_OUT}" \
    --test "${OUT_BASE}/requant_align_SlamQuant.out" \
    --min-pearson "${ALIGN_MIN_PEARSON}" \
    > "${OUT_BASE}/compare/requant_align.log" 2>&1
fi

{
  printf 'SLAM parity smoke: PASS\n'
  printf 'out_base=%s\n' "${OUT_BASE}"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'fastq=%s\n' "${FASTQ}"
  printf 'genome_dir=%s\n' "${STAR_INDEX}"
  printf 'snp_bed=%s\n' "${SNPS_BED}"
  printf 'reference=%s\n' "${REF_TSV}"
  printf 'ref_corr_min=%s\n' "${REF_CORR_MIN}"
  printf 'exact_min_pearson=%s\n' "${EXACT_MIN_PEARSON}"
  printf 'exact_max_abs_delta=%s\n' "${EXACT_MAX_ABS_DELTA}"
  printf 'align_min_pearson=%s\n' "${ALIGN_MIN_PEARSON}"
} > "${OUT_BASE}/PARITY_SUMMARY.txt"

cat "${OUT_BASE}/PARITY_SUMMARY.txt"
log "Done"
