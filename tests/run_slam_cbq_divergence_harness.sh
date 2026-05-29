#!/usr/bin/env bash
set -euo pipefail

# SLAM FASTQ-vs-CBQ divergence harness.
#
# This is intentionally stricter than the GRAND-SLAM correlation smoke. It first
# proves the FASTQ path is deterministic, then compares ordered CBQ input against
# the same FASTQ at progressively later boundaries:
#   input payload -> SAM body/SJ/log metrics -> SLAM diagnostics -> pre-NTR stats.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -f "${SCRIPT_DIR}/external_fixtures_env.sh" ]]; then
  # shellcheck source=/dev/null
  source "${SCRIPT_DIR}/external_fixtures_env.sh"
fi

OUT_BASE="${OUT_BASE:-/tmp/slam_cbq_divergence_harness_${STAMP}}"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
REQ_BIN="${REQ_BIN:-${ROOT_DIR}/slam/tools/slam_requant/slam_requant}"
FASTX_BIN="${FASTX_INPUT_HARNESS_BIN:-${ROOT_DIR}/core/legacy/source/fastx_input_harness}"
CBQ_BIN="${CBQ_READER_HARNESS_BIN:-${ROOT_DIR}/core/legacy/source/cbq_reader_harness}"
ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-${ROOT_DIR}/core/legacy/source/cbq_ordered_encoder}"

FASTQ="${FASTQ:-${SLAM_FIXTURE_FASTQ:-${ROOT_DIR}/test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz}}"
STAR_INDEX="${STAR_INDEX:-${SLAM_FIXTURE_STAR_INDEX:-${ROOT_DIR}/test/fixtures/slam/ref/star_index}}"
SNPS_BED="${SNPS_BED:-${SLAM_FIXTURE_SNPS_BED:-${ROOT_DIR}/test/fixtures/slam/ref/snps.bed}}"

THREADS="${THREADS:-4}"
CONTROL_THREADS="${CONTROL_THREADS:-1}"
RUN_THREAD_CONTROL="${RUN_THREAD_CONTROL:-1}"
EXACT_MAX_ABS_DELTA="${EXACT_MAX_ABS_DELTA:-1e-6}"
EXACT_MIN_PEARSON="${EXACT_MIN_PEARSON:-0.999999}"
ENCODER_COMPRESSION_LEVEL="${ENCODER_COMPRESSION_LEVEL:-0}"
ENCODER_BLOCK_SIZE="${ENCODER_BLOCK_SIZE:-1048576}"
AUTO_BUILD="${AUTO_BUILD:-1}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Runs a staged SLAM FASTQ-vs-CBQ divergence harness on a single-end fixture.

Options:
  --outdir DIR              Output root (default: ${OUT_BASE})
  --star-bin PATH           STAR binary (default: ${STAR_BIN})
  --requant-bin PATH        slam_requant binary (default: ${REQ_BIN})
  --fastx-bin PATH          fastx_input_harness binary (default: ${FASTX_BIN})
  --cbq-bin PATH            cbq_reader_harness binary (default: ${CBQ_BIN})
  --encoder-bin PATH        cbq_ordered_encoder binary (default: ${ENCODER_BIN})
  --fastq PATH              Single-end FASTQ[.gz] (default: ${FASTQ})
  --genome-dir DIR          STAR genomeDir (default: ${STAR_INDEX})
  --snp-bed PATH            SNP mask BED (default: ${SNPS_BED})
  --threads N               Main STAR threads (default: ${THREADS})
  --control-threads N       FASTQ thread-control threads (default: ${CONTROL_THREADS})
  --no-thread-control       Skip FASTQ 1-thread vs N-thread control
  --compression-level N     Ordered CBQ compression level (default: ${ENCODER_COMPRESSION_LEVEL})
  --block-size N            Ordered CBQ block size (default: ${ENCODER_BLOCK_SIZE})
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
    --fastx-bin) FASTX_BIN="$2"; shift 2 ;;
    --cbq-bin) CBQ_BIN="$2"; shift 2 ;;
    --encoder-bin) ENCODER_BIN="$2"; shift 2 ;;
    --fastq) FASTQ="$2"; shift 2 ;;
    --genome-dir) STAR_INDEX="$2"; shift 2 ;;
    --snp-bed) SNPS_BED="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --control-threads) CONTROL_THREADS="$2"; shift 2 ;;
    --no-thread-control) RUN_THREAD_CONTROL=0; shift ;;
    --compression-level) ENCODER_COMPRESSION_LEVEL="$2"; shift 2 ;;
    --block-size) ENCODER_BLOCK_SIZE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "--threads must be a positive integer"
[[ "${CONTROL_THREADS}" =~ ^[0-9]+$ && "${CONTROL_THREADS}" -gt 0 ]] || die "--control-threads must be a positive integer"
[[ -f "${FASTQ}" ]] || skip "SLAM FASTQ not found: ${FASTQ}"
[[ -d "${STAR_INDEX}" ]] || skip "SLAM STAR index not found: ${STAR_INDEX}"
[[ -f "${SNPS_BED}" ]] || skip "SLAM SNP BED not found: ${SNPS_BED}"

if [[ "${AUTO_BUILD}" == "1" ]]; then
  if [[ ! -x "${STAR_BIN}" ]]; then
    make -C "${ROOT_DIR}/core/legacy/source" -j"${THREADS}" STAR
  fi
  if [[ ! -x "${FASTX_BIN}" || ! -x "${CBQ_BIN}" || ! -x "${ENCODER_BIN}" ]]; then
    make -C "${ROOT_DIR}/core/legacy/source" -j"${THREADS}" \
      fastx-input-harness cbq-reader-harness cbq-ordered-encoder
  fi
  if [[ ! -x "${REQ_BIN}" && -f "${ROOT_DIR}/slam/tools/slam_requant/Makefile" ]]; then
    make -C "${ROOT_DIR}/slam/tools/slam_requant" -j"${THREADS}"
  fi
fi

[[ -x "${STAR_BIN}" ]] || skip "STAR binary not found: ${STAR_BIN}"
[[ -x "${FASTX_BIN}" ]] || skip "fastx_input_harness not found: ${FASTX_BIN}"
[[ -x "${CBQ_BIN}" ]] || skip "cbq_reader_harness not found: ${CBQ_BIN}"
[[ -x "${ENCODER_BIN}" ]] || skip "cbq_ordered_encoder not found: ${ENCODER_BIN}"
[[ -x "${REQ_BIN}" ]] || skip "slam_requant binary not found: ${REQ_BIN}"
[[ -f "${SCRIPT_DIR}/slam/compare_star_outputs.py" ]] || die "Missing compare_star_outputs.py"

mkdir -p "${OUT_BASE}/inputs" "${OUT_BASE}/logs" "${OUT_BASE}/compare" "${OUT_BASE}/stages"
REPORT="${OUT_BASE}/DIVERGENCE_REPORT.txt"
STAGE_TSV="${OUT_BASE}/stage_status.tsv"
FIRST_STAGE=""
FASTQ_CONTROL_FAILED=0

{
  printf 'SLAM CBQ divergence harness\n'
  printf 'timestamp=%s\n' "${STAMP}"
  printf 'fastq=%s\n' "${FASTQ}"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'requant_bin=%s\n' "${REQ_BIN}"
  printf 'genome_dir=%s\n' "${STAR_INDEX}"
  printf 'snp_bed=%s\n' "${SNPS_BED}"
  printf 'threads=%s\n' "${THREADS}"
  printf 'control_threads=%s\n' "${CONTROL_THREADS}"
  printf 'encoder_compression_level=%s\n' "${ENCODER_COMPRESSION_LEVEL}"
  printf 'encoder_block_size=%s\n' "${ENCODER_BLOCK_SIZE}"
  printf '\n'
} > "${REPORT}"
printf 'stage\tstatus\trole\tdetail\n' > "${STAGE_TSV}"

record_stage() {
  local stage="$1"
  local role="$2"
  local rc="$3"
  local detail="$4"
  local status="PASS"
  if [[ "${rc}" != "0" ]]; then
    status="FAIL"
    if [[ -z "${FIRST_STAGE}" ]]; then
      FIRST_STAGE="${stage}"
    fi
    if [[ "${role}" == fastq-control ]]; then
      FASTQ_CONTROL_FAILED=1
    fi
  fi
  printf '%s\t%s\t%s\t%s\n' "${stage}" "${status}" "${role}" "${detail}" | tee -a "${STAGE_TSV}" >> "${REPORT}"
}

cmp_stage() {
  local stage="$1"
  local role="$2"
  local ref="$3"
  local test="$4"
  local diff_out="${OUT_BASE}/compare/${stage}.diff"
  if cmp -s "${ref}" "${test}"; then
    record_stage "${stage}" "${role}" 0 "exact match"
  else
    diff -u "${ref}" "${test}" > "${diff_out}" 2>&1 || true
    record_stage "${stage}" "${role}" 1 "mismatch; see ${diff_out}"
  fi
}

compare_star_stage() {
  local stage="$1"
  local role="$2"
  local ref="$3"
  local test="$4"
  local out="${OUT_BASE}/compare/${stage}.log"
  if python3 "${SCRIPT_DIR}/slam/compare_star_outputs.py" \
      --reference "${ref}" \
      --test "${test}" \
      --min-pearson "${EXACT_MIN_PEARSON}" \
      --max-abs-delta "${EXACT_MAX_ABS_DELTA}" \
      --require-same-genes \
      > "${out}" 2>&1; then
    record_stage "${stage}" "${role}" 0 "exact pre-NTR/NTR match; see ${out}"
  else
    record_stage "${stage}" "${role}" 1 "pre-NTR/NTR mismatch; see ${out}"
  fi
}

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

normalize_outputs() {
  local label="$1"
  local run_dir="${OUT_BASE}/${label}"
  local prefix="${run_dir}/star_"
  local sam_body="${OUT_BASE}/stages/${label}.sam.body"
  local sam_sorted="${OUT_BASE}/stages/${label}.sam.body.sorted"
  local log_metrics="${OUT_BASE}/stages/${label}.Log.final.metrics"

  if [[ -s "${prefix}Aligned.out.sam" ]]; then
    grep -v '^@' "${prefix}Aligned.out.sam" > "${sam_body}" || true
    LC_ALL=C sort "${sam_body}" > "${sam_sorted}"
  fi

  if [[ -s "${prefix}Log.final.out" ]]; then
    grep '|' "${prefix}Log.final.out" \
      | grep -Ev 'Started job|Started mapping|Finished on|Mapping speed' \
      | grep -Ev 'FLEX STAGE TIMING|PIPELINE DIAGNOSTICS| calls | total ms | mean us/call ' \
      | grep -Ev 'Chunks processed|thread-seconds|throughput|Avg mutex wait|Avg chunk read|Avg map chunk' \
      | sed 's/[[:space:]]*$//' \
      > "${log_metrics}" || true
  fi
}

run_star_fastq() {
  local label="$1"
  local threads="$2"
  local run_dir="${OUT_BASE}/${label}"
  local prefix="${run_dir}/star_"
  mkdir -p "${run_dir}/logs"
  rm -rf "${run_dir}/tmp"

  local read_cmd=()
  if [[ "${FASTQ}" == *.gz ]]; then
    read_cmd=(--readFilesCommand zcat)
  fi

  local cmd=(
    "${STAR_BIN}"
    --runThreadN "${threads}"
    --genomeDir "${STAR_INDEX}"
    --readFilesType Fastx
    --readFilesIn "${FASTQ}"
    "${read_cmd[@]}"
    --outFileNamePrefix "${prefix}"
    --outTmpDir "${run_dir}/tmp"
    --outSAMtype SAM
    --clip3pAdapterSeq AGATCGGAAGAG
    --clip3pAdapterMMp 0.1
    --slamQuantMode 1
    --slamSnpMaskIn "${SNPS_BED}"
    --slamGrandSlamOut 1
    --slamDumpBinary "${run_dir}/star_slam.dump.bin"
    --slamDumpWeights "${run_dir}/star_slam.weights.bin"
  )
  write_command "${run_dir}/RUN_STAR.sh" "${cmd[@]}"
  log "Running ${label} (${threads} threads, FASTQ)"
  /usr/bin/time -v -o "${run_dir}/logs/star.time.log" \
    "${cmd[@]}" > "${run_dir}/logs/star.stdout.log" 2> "${run_dir}/logs/star.stderr.log"
  normalize_outputs "${label}"
}

run_star_cbq() {
  local label="$1"
  local threads="$2"
  local cbq="$3"
  local run_dir="${OUT_BASE}/${label}"
  local prefix="${run_dir}/star_"
  mkdir -p "${run_dir}/logs"
  rm -rf "${run_dir}/tmp"

  local cmd=(
    "${STAR_BIN}"
    --runThreadN "${threads}"
    --genomeDir "${STAR_INDEX}"
    --readFilesType Binseq SE
    --readFilesIn "${cbq}"
    --outFileNamePrefix "${prefix}"
    --outTmpDir "${run_dir}/tmp"
    --outSAMtype SAM
    --clip3pAdapterSeq AGATCGGAAGAG
    --clip3pAdapterMMp 0.1
    --slamQuantMode 1
    --slamSnpMaskIn "${SNPS_BED}"
    --slamGrandSlamOut 1
    --slamDumpBinary "${run_dir}/star_slam.dump.bin"
    --slamDumpWeights "${run_dir}/star_slam.weights.bin"
  )
  write_command "${run_dir}/RUN_STAR.sh" "${cmd[@]}"
  log "Running ${label} (${threads} threads, CBQ)"
  /usr/bin/time -v -o "${run_dir}/logs/star.time.log" \
    "${cmd[@]}" > "${run_dir}/logs/star.stdout.log" 2> "${run_dir}/logs/star.stderr.log"
  normalize_outputs "${label}"
}

run_requant() {
  local label="$1"
  local run_dir="${OUT_BASE}/${label}"
  [[ -s "${run_dir}/star_slam.dump.bin" ]] || die "Missing dump for ${label}"
  [[ -s "${run_dir}/star_slam.weights.bin" ]] || die "Missing weights for ${label}"

  log "Replaying ${label} dump"
  "${REQ_BIN}" \
    --dump "${run_dir}/star_slam.dump.bin" \
    --out "${run_dir}/requant_dump_" \
    --slamSnpMaskIn "${SNPS_BED}" \
    > "${run_dir}/logs/requant_dump.log" 2>&1

  "${REQ_BIN}" \
    --dump "${run_dir}/star_slam.dump.bin" \
    --out "${run_dir}/requant_weight_" \
    --slamSnpMaskIn "${SNPS_BED}" \
    --slamWeightFile "${run_dir}/star_slam.weights.bin" \
    > "${run_dir}/logs/requant_weight.log" 2>&1
}

CBQ_FILE="${OUT_BASE}/inputs/slam_fixture.cbq"
log "Encoding ordered CBQ fixture"
"${ENCODER_BIN}" \
  --readFilesIn "${FASTQ}" \
  --outFile "${CBQ_FILE}" \
  --compressionLevel "${ENCODER_COMPRESSION_LEVEL}" \
  --blockSize "${ENCODER_BLOCK_SIZE}" \
  > "${OUT_BASE}/logs/encode_cbq.stdout.log" \
  2> "${OUT_BASE}/logs/encode_cbq.stderr.log"
[[ -s "${CBQ_FILE}" ]] || die "Ordered CBQ encoder did not produce ${CBQ_FILE}"

log "Comparing FASTQ and CBQ input payloads"
"${FASTX_BIN}" --readNameSeparator space --readFilesIn "${FASTQ}" --dump-fastq \
  > "${OUT_BASE}/stages/input.fastq.dump"
"${CBQ_BIN}" --readNameSeparator space --mateCount 1 --readFilesIn "${CBQ_FILE}" --dump-fastq \
  > "${OUT_BASE}/stages/input.cbq.dump"
cmp_stage input_payload_fastq_vs_cbq cbq-specific \
  "${OUT_BASE}/stages/input.fastq.dump" \
  "${OUT_BASE}/stages/input.cbq.dump"

run_star_fastq fastq_a "${THREADS}"
run_star_fastq fastq_b "${THREADS}"
if [[ "${RUN_THREAD_CONTROL}" == "1" ]]; then
  run_star_fastq fastq_control "${CONTROL_THREADS}"
fi
run_star_cbq cbq_a "${THREADS}" "${CBQ_FILE}"

for label in fastq_a fastq_b cbq_a; do
  run_requant "${label}"
done
if [[ "${RUN_THREAD_CONTROL}" == "1" ]]; then
  run_requant fastq_control
fi

cmp_stage fastq_repeat_sam_body_exact fastq-control \
  "${OUT_BASE}/stages/fastq_a.sam.body" \
  "${OUT_BASE}/stages/fastq_b.sam.body"
cmp_stage fastq_repeat_sam_body_sorted fastq-control \
  "${OUT_BASE}/stages/fastq_a.sam.body.sorted" \
  "${OUT_BASE}/stages/fastq_b.sam.body.sorted"
cmp_stage fastq_repeat_sj_exact fastq-control \
  "${OUT_BASE}/fastq_a/star_SJ.out.tab" \
  "${OUT_BASE}/fastq_b/star_SJ.out.tab"
cmp_stage fastq_repeat_log_metrics fastq-control \
  "${OUT_BASE}/stages/fastq_a.Log.final.metrics" \
  "${OUT_BASE}/stages/fastq_b.Log.final.metrics"
cmp_stage fastq_repeat_slam_diagnostics fastq-control \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out.diagnostics" \
  "${OUT_BASE}/fastq_b/star_SlamQuant.out.diagnostics"
compare_star_stage fastq_repeat_pre_ntr fastq-control \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out" \
  "${OUT_BASE}/fastq_b/star_SlamQuant.out"

if [[ "${RUN_THREAD_CONTROL}" == "1" ]]; then
  cmp_stage fastq_thread_sam_body_sorted fastq-control \
    "${OUT_BASE}/stages/fastq_control.sam.body.sorted" \
    "${OUT_BASE}/stages/fastq_a.sam.body.sorted"
  cmp_stage fastq_thread_log_metrics fastq-control \
    "${OUT_BASE}/stages/fastq_control.Log.final.metrics" \
    "${OUT_BASE}/stages/fastq_a.Log.final.metrics"
  cmp_stage fastq_thread_slam_diagnostics fastq-control \
    "${OUT_BASE}/fastq_control/star_SlamQuant.out.diagnostics" \
    "${OUT_BASE}/fastq_a/star_SlamQuant.out.diagnostics"
  compare_star_stage fastq_thread_pre_ntr fastq-control \
    "${OUT_BASE}/fastq_control/star_SlamQuant.out" \
    "${OUT_BASE}/fastq_a/star_SlamQuant.out"
fi

compare_star_stage fastq_a_dump_replay_pre_ntr internal-replay \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out" \
  "${OUT_BASE}/fastq_a/requant_dump_SlamQuant.out"
compare_star_stage fastq_a_weight_replay_pre_ntr internal-replay \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out" \
  "${OUT_BASE}/fastq_a/requant_weight_SlamQuant.out"
compare_star_stage cbq_a_dump_replay_pre_ntr internal-replay \
  "${OUT_BASE}/cbq_a/star_SlamQuant.out" \
  "${OUT_BASE}/cbq_a/requant_dump_SlamQuant.out"
compare_star_stage cbq_a_weight_replay_pre_ntr internal-replay \
  "${OUT_BASE}/cbq_a/star_SlamQuant.out" \
  "${OUT_BASE}/cbq_a/requant_weight_SlamQuant.out"

cmp_stage fastq_vs_cbq_sam_body_exact cbq-specific \
  "${OUT_BASE}/stages/fastq_a.sam.body" \
  "${OUT_BASE}/stages/cbq_a.sam.body"
cmp_stage fastq_vs_cbq_sam_body_sorted cbq-specific \
  "${OUT_BASE}/stages/fastq_a.sam.body.sorted" \
  "${OUT_BASE}/stages/cbq_a.sam.body.sorted"
cmp_stage fastq_vs_cbq_sj_exact cbq-specific \
  "${OUT_BASE}/fastq_a/star_SJ.out.tab" \
  "${OUT_BASE}/cbq_a/star_SJ.out.tab"
cmp_stage fastq_vs_cbq_log_metrics cbq-specific \
  "${OUT_BASE}/stages/fastq_a.Log.final.metrics" \
  "${OUT_BASE}/stages/cbq_a.Log.final.metrics"
cmp_stage fastq_vs_cbq_slam_diagnostics cbq-specific \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out.diagnostics" \
  "${OUT_BASE}/cbq_a/star_SlamQuant.out.diagnostics"
cmp_stage fastq_vs_cbq_transitions cbq-specific \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out.transitions.tsv" \
  "${OUT_BASE}/cbq_a/star_SlamQuant.out.transitions.tsv"
compare_star_stage fastq_vs_cbq_pre_ntr cbq-specific \
  "${OUT_BASE}/fastq_a/star_SlamQuant.out" \
  "${OUT_BASE}/cbq_a/star_SlamQuant.out"

{
  printf '\nConclusion\n'
  if [[ -z "${FIRST_STAGE}" ]]; then
    printf 'first_divergence=none\n'
    printf 'classification=pass\n'
    printf 'interpretation=FASTQ controls, CBQ input payload, alignment surfaces, SLAM diagnostics, and pre-NTR stats are exact within configured tolerance.\n'
  elif [[ "${FASTQ_CONTROL_FAILED}" == "1" ]]; then
    printf 'first_divergence=%s\n' "${FIRST_STAGE}"
    printf 'classification=not_cbq_specific\n'
    printf 'interpretation=A FASTQ-only control failed before or alongside the CBQ comparison; investigate FASTQ determinism/threading first.\n'
  else
    printf 'first_divergence=%s\n' "${FIRST_STAGE}"
    printf 'classification=cbq_specific\n'
    printf 'interpretation=FASTQ controls passed; first failing CBQ boundary above is the earliest observed divergence point.\n'
  fi
  printf 'stage_status=%s\n' "${STAGE_TSV}"
} >> "${REPORT}"

cat "${REPORT}"

if [[ -n "${FIRST_STAGE}" ]]; then
  exit 1
fi
