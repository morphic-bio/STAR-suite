#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
VALIDATOR="${SCRIPT_DIR}/validate_partition_manifest.py"

MANIFEST=""
OUTPUT_DIR=""
STAR_BIN="${STAR_BIN:-STAR}"
GENOME_DIR=""
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
THREADS_PER_PARTITION="${THREADS_PER_PARTITION:-8}"
MAX_PARALLEL="${MAX_PARALLEL:-1}"
GATHER_THREADS="${GATHER_THREADS:-8}"
READ_FILES_COMMAND="${READ_FILES_COMMAND:-zcat}"
EXECUTE="0"
SKIP_GATHER="0"
declare -a EXTRA_STAR_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  run_partitioned_bulk_pe_alignment.sh --partition-manifest FILE \
    --output-dir DIR --star-bin PATH --genome-dir DIR [options]

Consumes a provider-neutral biodepot.partition_manifest/v1 document. The
recipe does not create partitions. It validates contiguous logical coverage,
runs one paired-end STAR alignment per declared partition, records completion,
and gathers coordinate-sorted BAMs only after every partition succeeds.

Options:
  --partition-manifest FILE   JSON partition manifest (required)
  --output-dir DIR            fresh or resumable output root (required)
  --star-bin PATH             STAR executable
  --genome-dir DIR            STAR genomeDir (required)
  --samtools-bin PATH         samtools executable
  --threads-per-partition N   STAR threads per partition (default: 8)
  --max-parallel N            concurrent STAR processes (default: 1)
  --gather-threads N          samtools merge/index threads (default: 8)
  --read-files-command CMD    STAR --readFilesCommand value; empty omits it
  --extra-star-arg ARG        additional single STAR argument; repeatable
  --skip-gather               retain partition BAMs without merging
  --execute                   execute; otherwise print a read-only plan
  --help                      show this help
EOF
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

print_cmd() {
  printf '  '
  printf '%q ' "$@"
  printf '\n'
}

while (($#)); do
  case "$1" in
    --partition-manifest) MANIFEST="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --samtools-bin) SAMTOOLS_BIN="$2"; shift 2 ;;
    --threads-per-partition) THREADS_PER_PARTITION="$2"; shift 2 ;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
    --gather-threads) GATHER_THREADS="$2"; shift 2 ;;
    --read-files-command) READ_FILES_COMMAND="$2"; shift 2 ;;
    --extra-star-arg) EXTRA_STAR_ARGS+=("$2"); shift 2 ;;
    --skip-gather) SKIP_GATHER="1"; shift ;;
    --execute) EXECUTE="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "${MANIFEST}" ]] || die "--partition-manifest is required"
[[ -n "${OUTPUT_DIR}" ]] || die "--output-dir is required"
[[ -n "${GENOME_DIR}" ]] || die "--genome-dir is required"
for value in THREADS_PER_PARTITION MAX_PARALLEL GATHER_THREADS; do
  [[ "${!value}" =~ ^[1-9][0-9]*$ ]] || die "${value} must be a positive integer"
done
[[ -x "${VALIDATOR}" ]] || die "missing manifest validator: ${VALIDATOR}"

VALIDATE_ARGS=("${VALIDATOR}" "${MANIFEST}")
[[ "${EXECUTE}" == "1" ]] && VALIDATE_ARGS+=(--check-files)
mapfile -t PARTITIONS < <("${VALIDATE_ARGS[@]}" --emit-tsv)
((${#PARTITIONS[@]} > 0)) || die "manifest contains no partitions"

if [[ "${EXECUTE}" == "1" ]]; then
  [[ -d "${GENOME_DIR}" ]] || die "missing genome directory: ${GENOME_DIR}"
  if [[ "${STAR_BIN}" != */* ]]; then STAR_BIN="$(command -v "${STAR_BIN}" || true)"; fi
  [[ -x "${STAR_BIN}" ]] || die "STAR executable not found: ${STAR_BIN:-unset}"
  if [[ "${SKIP_GATHER}" != "1" ]]; then
    if [[ "${SAMTOOLS_BIN}" != */* ]]; then SAMTOOLS_BIN="$(command -v "${SAMTOOLS_BIN}" || true)"; fi
    [[ -x "${SAMTOOLS_BIN}" ]] || die "samtools executable not found: ${SAMTOOLS_BIN:-unset}"
  fi
  mkdir -p "${OUTPUT_DIR}/partitions"
else
  printf 'DRY RUN: validated %s partitions; no files will be written.\n' "${#PARTITIONS[@]}"
fi

run_partition() {
  local ordinal="$1" part_id="$2" mate1="$3" mate2="$4" start="$5" end="$6"
  local part_dir="${OUTPUT_DIR}/partitions/$(printf '%06d' "${ordinal}")_${part_id}"
  local complete="${part_dir}/COMPLETE"
  local -a command=(
    "${STAR_BIN}" --runThreadN "${THREADS_PER_PARTITION}"
    --genomeDir "${GENOME_DIR}" --readFilesIn "${mate1}" "${mate2}"
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${part_dir}/"
  )
  [[ -n "${READ_FILES_COMMAND}" ]] && command+=(--readFilesCommand "${READ_FILES_COMMAND}")
  command+=("${EXTRA_STAR_ARGS[@]}")

  if [[ "${EXECUTE}" != "1" ]]; then
    printf 'partition %s [%s,%s):\n' "${part_id}" "${start}" "${end}"
    print_cmd "${command[@]}"
    return 0
  fi
  if [[ -f "${complete}" && -s "${part_dir}/Aligned.sortedByCoord.out.bam" ]]; then
    printf 'resume: %s already complete\n' "${part_id}"
    return 0
  fi
  mkdir -p "${part_dir}"
  print_cmd "${command[@]}" > "${part_dir}/command.txt"
  "${command[@]}" > "${part_dir}/stdout.log" 2> "${part_dir}/stderr.log"
  [[ -s "${part_dir}/Aligned.sortedByCoord.out.bam" ]] || return 1
  [[ -s "${part_dir}/Log.final.out" ]] || return 1
  {
    printf 'schema=biodepot.partition_completion/v1\n'
    printf 'ordinal=%s\n' "${ordinal}"
    printf 'partition_id=%s\n' "${part_id}"
    printf 'logical_start=%s\n' "${start}"
    printf 'logical_end=%s\n' "${end}"
    printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "${complete}"
}

if [[ "${EXECUTE}" != "1" ]]; then
  for row in "${PARTITIONS[@]}"; do
    IFS=$'\t' read -r ordinal part_id mate1 mate2 start end <<< "${row}"
    run_partition "${ordinal}" "${part_id}" "${mate1}" "${mate2}" "${start}" "${end}"
  done
  if [[ "${SKIP_GATHER}" != "1" ]]; then
    print_cmd "${SAMTOOLS_BIN}" merge -@ "${GATHER_THREADS}" -o "${OUTPUT_DIR}/Aligned.sortedByCoord.out.bam" '<partition BAMs>'
  fi
  exit 0
fi

declare -a BATCH_PIDS=()
wait_batch() {
  local pid failed=0
  for pid in "${BATCH_PIDS[@]}"; do
    wait "${pid}" || failed=1
  done
  BATCH_PIDS=()
  ((failed == 0)) || die "one or more partitions failed; gather was not started"
}

for row in "${PARTITIONS[@]}"; do
  IFS=$'\t' read -r ordinal part_id mate1 mate2 start end <<< "${row}"
  run_partition "${ordinal}" "${part_id}" "${mate1}" "${mate2}" "${start}" "${end}" &
  BATCH_PIDS+=("$!")
  ((${#BATCH_PIDS[@]} < MAX_PARALLEL)) || wait_batch
done
((${#BATCH_PIDS[@]} == 0)) || wait_batch

if [[ "${SKIP_GATHER}" == "1" ]]; then
  printf 'PASS: all partitions complete; gather skipped\n'
  exit 0
fi

declare -a PARTITION_BAMS=()
for row in "${PARTITIONS[@]}"; do
  IFS=$'\t' read -r ordinal part_id _ <<< "${row}"
  part_dir="${OUTPUT_DIR}/partitions/$(printf '%06d' "${ordinal}")_${part_id}"
  [[ -f "${part_dir}/COMPLETE" ]] || die "missing completion evidence for ${part_id}"
  PARTITION_BAMS+=("${part_dir}/Aligned.sortedByCoord.out.bam")
done
FINAL_BAM="${OUTPUT_DIR}/Aligned.sortedByCoord.out.bam"
"${SAMTOOLS_BIN}" merge -f -@ "${GATHER_THREADS}" -o "${FINAL_BAM}" "${PARTITION_BAMS[@]}"
"${SAMTOOLS_BIN}" index -@ "${GATHER_THREADS}" "${FINAL_BAM}"
printf 'PASS: gathered %s partitions into %s\n' "${#PARTITION_BAMS[@]}" "${FINAL_BAM}"
