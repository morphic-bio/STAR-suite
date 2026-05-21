#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

RAW_DIR="${RAW_DIR:-/mnt/pikachu/JAX_scRNAseq02/raw}"
ORACLE_DIR="${ORACLE_DIR:-/mnt/pikachu/JAX_scRNAseq02/cellranger-logs}"
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_$(date -u +%Y%m%dT%H%M%SZ)}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt}"
SAMPLE_STEM="${SAMPLE_STEM:-25E32-L3_GT25-03394_ACCTCGAGCT-ATCGAACACA_S44}"
SAMPLE_ID="${SAMPLE_ID:-25E32-L3}"
DOWNSAMPLE_READ_PAIRS="${DOWNSAMPLE_READ_PAIRS:-2000000}"
THREADS="${THREADS:-16}"
SOLO_INLINE_HASH_MODE="${SOLO_INLINE_HASH_MODE:-no}"
RUN_STAR="0"
RUN_VALIDATE="0"
FORCE="0"
KEEP_TMP="0"
FULL_FASTQS="0"
EXTRA_STAR_ARGS=()

usage() {
  cat <<'EOF'
Usage:
  run_jax_scrnaseq02_ocm_oracle_smoke.sh [options] [-- <extra STAR args>]

Stages a 2M read-pair downsample, or full-depth FASTQ symlinks with
--full-fastqs, of the JAX scRNAseq02 OCM oracle library 25E32-L3 and renders a
STAR smoke command. It does not compile STAR.

By default the script stages FASTQs and writes RUN_STAR.sh only. Use --run-star
after the current production mapping job is done. Use --validate to run the
Cell Ranger multi layout validator after STAR completes or against an existing
run directory.

Options:
  --raw-dir PATH              Raw FASTQ directory
  --oracle-dir PATH           Cell Ranger oracle/log directory
  --out-root PATH             Smoke output root
  --star-bin PATH             Existing STAR binary to use
  --genome-dir PATH           STAR genomeDir
  --solo-cb-whitelist PATH    GEM-X TRU whitelist
  --downsample-read-pairs N   Total read pairs across both lanes (default: 2000000)
  --full-fastqs               Symlink the complete 25E32-L3 FASTQs instead of
                              staging a downsample
  --threads N                 STAR threads (default: 16)
  --solo-inline-hash-mode M   STAR --soloInlineHashMode value (default: no).
                              Use no for this BAM/Y-removal + Velocyto MEX
                              smoke; yes is only safe on the direct hash-bridge
                              path that does not write BAM.
  --run-star                  Execute RUN_STAR.sh after staging
  --validate                  Run oracle validator after staging or STAR run
  --force                     Recreate staged FASTQs and scripts
  --keep-tmp                  Keep STAR tmp directory after a successful run
  -h, --help                  Show help

Extra STAR arguments after -- are appended verbatim for development-only
overrides. The native OCM compatibility flags are passed by default.
EOF
}

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-dir) RAW_DIR="$2"; shift 2 ;;
    --oracle-dir) ORACLE_DIR="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --solo-cb-whitelist) SOLO_CB_WHITELIST="$2"; shift 2 ;;
    --downsample-read-pairs) DOWNSAMPLE_READ_PAIRS="$2"; shift 2 ;;
    --full-fastqs) FULL_FASTQS="1"; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --solo-inline-hash-mode) SOLO_INLINE_HASH_MODE="$2"; shift 2 ;;
    --run-star) RUN_STAR="1"; shift ;;
    --validate) RUN_VALIDATE="1"; shift ;;
    --force) FORCE="1"; shift ;;
    --keep-tmp) KEEP_TMP="1"; shift ;;
    --) shift; EXTRA_STAR_ARGS+=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ "${DOWNSAMPLE_READ_PAIRS}" =~ ^[0-9]+$ ]] || die "--downsample-read-pairs must be an integer"
[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "--threads must be a positive integer"
[[ "${SOLO_INLINE_HASH_MODE}" =~ ^(yes|no|auto)$ ]] || die "--solo-inline-hash-mode must be yes, no, or auto"
[[ -d "${RAW_DIR}" ]] || die "Missing raw FASTQ directory: ${RAW_DIR}"
[[ -d "${ORACLE_DIR}" ]] || die "Missing oracle directory: ${ORACLE_DIR}"
[[ -x "${STAR_BIN}" ]] || die "Missing executable STAR binary: ${STAR_BIN}"
[[ -d "${GENOME_DIR}" ]] || die "Missing genomeDir: ${GENOME_DIR}"
[[ -f "${SOLO_CB_WHITELIST}" ]] || die "Missing whitelist: ${SOLO_CB_WHITELIST}"
[[ -f "${ORACLE_DIR}/config.csv" ]] || die "Missing oracle config.csv"
[[ -f "${ORACLE_DIR}/cells_per_tag.json" ]] || die "Missing oracle cells_per_tag.json"
[[ -f "${ORACLE_DIR}/25E32-L3_Day4-pool-1.mri.tgz" ]] || die "Missing oracle MRI archive"

OUT_ROOT="$(realpath -m "${OUT_ROOT}")"
STAGE_DIR="${OUT_ROOT}/stage_fastqs"
SAMPLE_DIR="${OUT_ROOT}/samples/${SAMPLE_ID}"
RUN_DIR="${SAMPLE_DIR}/run"
LOG_DIR="${OUT_ROOT}/logs"
TMP_DIR="${SAMPLE_DIR}/tmp"
MANIFEST="${OUT_ROOT}/downsample_manifest.tsv"
RUN_SCRIPT="${OUT_ROOT}/RUN_STAR.sh"
VALIDATION_JSON="${OUT_ROOT}/ocm_oracle_validation.json"

mkdir -p "${STAGE_DIR}" "${RUN_DIR}" "${LOG_DIR}"

LANES=(L007 L008)

source_fastq() {
  local lane="$1"
  local read="$2"
  printf '%s/%s_%s_%s_001.fastq.gz' "${RAW_DIR}" "${SAMPLE_STEM}" "${lane}" "${read}"
}

stage_fastq() {
  local src_r1="$1"
  local src_r2="$2"
  local dst_r1="$3"
  local dst_r2="$4"
  local read_pairs="$5"

  python3 - "$src_r1" "$src_r2" "$dst_r1" "$dst_r2" "$read_pairs" <<'PY'
import gzip
import sys
from pathlib import Path

src_r1, src_r2, dst_r1, dst_r2, n_s = sys.argv[1:]
n = int(n_s)
Path(dst_r1).parent.mkdir(parents=True, exist_ok=True)

def copy_pair():
    with gzip.open(src_r1, "rt") as r1_in, gzip.open(src_r2, "rt") as r2_in, \
         gzip.open(dst_r1, "wt", compresslevel=1) as r1_out, \
         gzip.open(dst_r2, "wt", compresslevel=1) as r2_out:
        written = 0
        for _ in range(n):
            rec1 = [r1_in.readline() for _ in range(4)]
            rec2 = [r2_in.readline() for _ in range(4)]
            if not all(rec1) or not all(rec2):
                raise RuntimeError(f"Source FASTQs ended before {n} read pairs; wrote {written}")
            r1_out.writelines(rec1)
            r2_out.writelines(rec2)
            written += 1
        return written

written = copy_pair()
print(written)
PY
}

stage_inputs() {
  if [[ -f "${MANIFEST}" && "${FORCE}" != "1" ]]; then
    log "Input manifest exists; reusing staged FASTQs: ${MANIFEST}"
    return
  fi

  if [[ "${FORCE}" == "1" ]]; then
    rm -rf "${STAGE_DIR}" "${MANIFEST}"
    mkdir -p "${STAGE_DIR}"
  fi

  if [[ "${FULL_FASTQS}" == "1" ]]; then
    printf 'sample_id\tlane\tread_pairs\tr1\tr2\n' > "${MANIFEST}"
    local full_lane full_src_r1 full_src_r2 full_dst_r1 full_dst_r2
    for full_lane in "${LANES[@]}"; do
      full_src_r1="$(source_fastq "${full_lane}" R1)"
      full_src_r2="$(source_fastq "${full_lane}" R2)"
      [[ -f "${full_src_r1}" ]] || die "Missing source FASTQ: ${full_src_r1}"
      [[ -f "${full_src_r2}" ]] || die "Missing source FASTQ: ${full_src_r2}"
      full_dst_r1="${STAGE_DIR}/$(basename "${full_src_r1}")"
      full_dst_r2="${STAGE_DIR}/$(basename "${full_src_r2}")"
      ln -sfn "${full_src_r1}" "${full_dst_r1}"
      ln -sfn "${full_src_r2}" "${full_dst_r2}"
      printf '%s\t%s\t%s\t%s\t%s\n' "${SAMPLE_ID}" "${full_lane}" "full" "${full_dst_r1}" "${full_dst_r2}" >> "${MANIFEST}"
      log "Linked full ${full_lane} FASTQs"
    done
    return
  fi

  local lane_count="${#LANES[@]}"
  local base_count=$((DOWNSAMPLE_READ_PAIRS / lane_count))
  local remainder=$((DOWNSAMPLE_READ_PAIRS % lane_count))

  printf 'sample_id\tlane\tread_pairs\tr1\tr2\n' > "${MANIFEST}"
  local i lane lane_pairs src_r1 src_r2 dst_r1 dst_r2 written
  for i in "${!LANES[@]}"; do
    lane="${LANES[$i]}"
    lane_pairs="${base_count}"
    if (( i < remainder )); then
      lane_pairs=$((lane_pairs + 1))
    fi

    src_r1="$(source_fastq "${lane}" R1)"
    src_r2="$(source_fastq "${lane}" R2)"
    [[ -f "${src_r1}" ]] || die "Missing source FASTQ: ${src_r1}"
    [[ -f "${src_r2}" ]] || die "Missing source FASTQ: ${src_r2}"

    dst_r1="${STAGE_DIR}/$(basename "${src_r1}")"
    dst_r2="${STAGE_DIR}/$(basename "${src_r2}")"
    log "Staging ${lane}: ${lane_pairs} read pairs"
    written="$(stage_fastq "${src_r1}" "${src_r2}" "${dst_r1}" "${dst_r2}" "${lane_pairs}")"
    [[ "${written}" == "${lane_pairs}" ]] || die "Wrote ${written} read pairs for ${lane}, expected ${lane_pairs}"
    printf '%s\t%s\t%s\t%s\t%s\n' "${SAMPLE_ID}" "${lane}" "${lane_pairs}" "${dst_r1}" "${dst_r2}" >> "${MANIFEST}"
  done
}

join_by_comma() {
  local IFS=,
  printf '%s' "$*"
}

render_star_script() {
  local r1_files=()
  local r2_files=()
  local lane
  for lane in "${LANES[@]}"; do
    r1_files+=("${STAGE_DIR}/$(basename "$(source_fastq "${lane}" R1)")")
    r2_files+=("${STAGE_DIR}/$(basename "$(source_fastq "${lane}" R2)")")
  done

  local r1_csv r2_csv
  r1_csv="$(join_by_comma "${r1_files[@]}")"
  r2_csv="$(join_by_comma "${r2_files[@]}")"

  local cmd=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --dynamicThreadInterface 1
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${r2_csv}" "${r1_csv}"
    --readFilesCommand zcat
    --outFileNamePrefix "${RUN_DIR}/"
    --outTmpDir "${TMP_DIR}"
    --outSAMtype BAM Unsorted
    --emitNoYBAM yes
    --emitYNoYFastq yes
    --emitYNoYFastqCompression gz
    --clipAdapterType CellRanger4
    --clip3pPolyG yes
    --alignEndsType Local
    --chimSegmentMin 1000000
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloBarcodeReadLength 0
    --soloCBwhitelist "${SOLO_CB_WHITELIST}"
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloInlineCBCorrection yes
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloMultiMappers Unique
    --soloCellFilter EmptyDrops_CR
    --soloCbUbRequireTogether no
    --soloStrand Forward
    --soloFeatures GeneFull Velocyto
    --soloCrGexFeature genefull
    --soloCrMultimapRescue yes
    --soloInlineHashMode "${SOLO_INLINE_HASH_MODE}"
    --ocmMultiEnable yes
    --ocmMultiConfig "${ORACLE_DIR}/config.csv"
    --ocmMultiBarcodeMode flex
    --ocmMultiOutputCompat cellranger
  )
  if (( ${#EXTRA_STAR_ARGS[@]} > 0 )); then
    cmd+=("${EXTRA_STAR_ARGS[@]}")
  fi

  {
    printf '#!/usr/bin/env bash\n'
    printf 'set -euo pipefail\n\n'
    printf 'cd %q\n' "${REPO_ROOT}"
    printf ': "${STAR_VELOCYTO_LOW_MEM:=1}"\n'
    printf ': "${STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS:=8192}"\n'
    printf ': "${STAR_VELOCYTO_UMI_RESERVE_CAP:=32}"\n'
    printf ': "${STAR_SOLO_BINARY_SPOOL:=1}"\n'
    printf ': "${MALLOC_ARENA_MAX:=2}"\n'
    printf ': "${MALLOC_TRIM_THRESHOLD_:=131072}"\n'
    printf 'export STAR_VELOCYTO_LOW_MEM STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS STAR_VELOCYTO_UMI_RESERVE_CAP STAR_SOLO_BINARY_SPOOL MALLOC_ARENA_MAX MALLOC_TRIM_THRESHOLD_\n'
    printf 'mkdir -p %q %q\n' "${RUN_DIR}" "${LOG_DIR}"
    printf 'printf "started_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/STAR_STARTED.txt"
    printf 'printf "STAR_VELOCYTO_LOW_MEM=%%s\\nSTAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS=%%s\\nSTAR_VELOCYTO_UMI_RESERVE_CAP=%%s\\nSTAR_SOLO_BINARY_SPOOL=%%s\\nMALLOC_ARENA_MAX=%%s\\nMALLOC_TRIM_THRESHOLD_=%%s\\n" "$STAR_VELOCYTO_LOW_MEM" "$STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS" "$STAR_VELOCYTO_UMI_RESERVE_CAP" "$STAR_SOLO_BINARY_SPOOL" "$MALLOC_ARENA_MAX" "$MALLOC_TRIM_THRESHOLD_" > %q\n' "${LOG_DIR}/star.env.txt"
    printf 'cmd=('
    local arg
    for arg in "${cmd[@]}"; do
      printf ' %q' "${arg}"
    done
    printf ' )\n'
    printf 'printf "STAR command:" | tee %q\n' "${LOG_DIR}/star.command.txt"
    printf 'printf " %%q" "${cmd[@]}" | tee -a %q\n' "${LOG_DIR}/star.command.txt"
    printf 'printf "\\n" | tee -a %q\n' "${LOG_DIR}/star.command.txt"
    printf '"${cmd[@]}" 2>&1 | tee %q\n' "${LOG_DIR}/star.log"
    if [[ "${KEEP_TMP}" != "1" ]]; then
      printf 'rm -rf %q\n' "${TMP_DIR}"
    fi
    printf 'printf "completed_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/STAR_COMPLETED.txt"
  } > "${RUN_SCRIPT}"
  chmod +x "${RUN_SCRIPT}"
}

run_validator() {
  python3 "${REPO_ROOT}/scripts/validate_jax_scrnaseq02_ocm_oracle.py" \
    --star-run-dir "${RUN_DIR}" \
    --oracle-dir "${ORACLE_DIR}" \
    --report-json "${VALIDATION_JSON}"
}

stage_inputs
render_star_script

cat > "${OUT_ROOT}/README.txt" <<EOF
JAX scRNAseq02 OCM oracle smoke

sample_id=${SAMPLE_ID}
sample_stem=${SAMPLE_STEM}
downsample_read_pairs=${DOWNSAMPLE_READ_PAIRS}
full_fastqs=${FULL_FASTQS}
stage_dir=${STAGE_DIR}
run_dir=${RUN_DIR}
oracle_dir=${ORACLE_DIR}
star_bin=${STAR_BIN}
genome_dir=${GENOME_DIR}
solo_cb_whitelist=${SOLO_CB_WHITELIST}
solo_inline_hash_mode=${SOLO_INLINE_HASH_MODE}

Run STAR later:
  ${RUN_SCRIPT}

Validate existing or completed run:
  ${REPO_ROOT}/scripts/validate_jax_scrnaseq02_ocm_oracle.py --star-run-dir ${RUN_DIR} --oracle-dir ${ORACLE_DIR}
EOF

log "Prepared OCM oracle smoke harness: ${OUT_ROOT}"
log "STAR script: ${RUN_SCRIPT}"

if [[ "${RUN_STAR}" == "1" ]]; then
  log "Running STAR smoke"
  "${RUN_SCRIPT}"
fi

if [[ "${RUN_VALIDATE}" == "1" ]]; then
  log "Running Cell Ranger multi layout validation"
  run_validator
fi
