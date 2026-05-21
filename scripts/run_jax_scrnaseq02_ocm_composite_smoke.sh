#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

RAW_DIR="${RAW_DIR:-/mnt/pikachu/JAX_scRNAseq02/raw}"
CONFIG="${CONFIG:-/mnt/pikachu/JAX_scRNAseq02/cellranger-logs/config.csv}"
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/JAX_scRNAseq02_processed/ocm_composite_smoke_$(date -u +%Y%m%dT%H%M%SZ)}"
SAMPLE_ID="${SAMPLE_ID:-25E32-L3}"
SAMPLE_STEM="${SAMPLE_STEM:-25E32-L3_GT25-03394_ACCTCGAGCT-ATCGAACACA_S44}"
LANES="${LANES:-L007,L008}"
READ_PAIRS="${READ_PAIRS:-100000}"
THREADS="${THREADS:-16}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt}"
STAR_OCM_BARCODE_MODE="${STAR_OCM_BARCODE_MODE:-flex}"
STAR_NATIVE_OCM_BAM_SPLIT="${STAR_NATIVE_OCM_BAM_SPLIT:-yes}"
STAR_NATIVE_OCM_MEX="${STAR_NATIVE_OCM_MEX:-auto}"
STAR_YREMOVE="${STAR_YREMOVE:-yes}"
STAR_BAM_CBUB_TAGS="${STAR_BAM_CBUB_TAGS:-no}"
STAR_BAM_GXGN_TAGS="${STAR_BAM_GXGN_TAGS:-no}"
STAR_SOLO_FEATURES="${STAR_SOLO_FEATURES:-GeneFull Velocyto}"
STAR_COMPARE_FEATURE="${STAR_COMPARE_FEATURE:-GeneFull}"
CR_REFERENCE="${CR_REFERENCE:-/mnt/pikachu/CR-references/refdata-gex-GRCh38-2024-A}"
CELLRANGER_BIN="${CELLRANGER_BIN:-cellranger}"
CR_LOCALMEM="${CR_LOCALMEM:-80}"
CR_CREATE_BAM="${CR_CREATE_BAM:-true}"
CR_REUSE_RUN_DIR="${CR_REUSE_RUN_DIR:-auto}"
CR_REUSE_SEARCH_ROOT="${CR_REUSE_SEARCH_ROOT:-/mnt/pikachu/JAX_scRNAseq02_processed}"
OCM_PREP_REUSE_ROOT="${OCM_PREP_REUSE_ROOT:-}"
SIMPLEED_SIM_N="${SIMPLEED_SIM_N:-100000}"
SIMPLEED_MODE="${SIMPLEED_MODE:-full}"
OCM_MATERIALIZE_THREADS="${OCM_MATERIALIZE_THREADS:-4}"
STAR_OUTSAMTYPE="${STAR_OUTSAMTYPE:-BAM Unsorted}"

RUN_PREPARE=0
RUN_STAR=0
RUN_SPLIT_BAM=0
RUN_MATERIALIZE=0
RUN_CR=0
RUN_MERGE_CR_BAM=0
RUN_COMPARE=0
RUN_ALL_REQUESTED=0
RUN_CR_EXPLICIT=0
FORCE=0

usage() {
  cat <<'EOF'
Usage:
  run_jax_scrnaseq02_ocm_composite_smoke.sh [options]

Stages a JAX scRNAseq02 OCM downsample, runs the native STAR OCM-Flex
composite-barcode smoke, and compares STAR per-sample MEX outputs against a
completed Cell Ranger 9 multi reference for the same deterministic downsample.

Options:
  --read-pairs N             Total read pairs across lanes (default: 100000)
  --out-root PATH            Output root
  --threads N                STAR / CR local cores (default: 16)
  --run-all                  Run prepare, STAR, materialize STAR comparison feature, and compare; reuse CR by default
  --prepare                  Stage downsample and render run scripts
  --run-star                 Execute rendered STAR command
  --split-bam                Split STAR pooled BAM into 4 per-OCM BAMs
  --materialize              Materialize per-sample STAR MEX outputs
  --run-cr                   Force execution of rendered Cell Ranger multi command
  --merge-cr-bam             Merge CR per-sample/unassigned BAMs into one BAM
  --compare                  Compare STAR materialization to CR9 multi
  --force                    Recreate staged/scripted outputs for selected steps
  --star-out-samtype VALUE   None or "BAM Unsorted" (default: "BAM Unsorted")
  --star-yremove yes|no      Emit STAR Y/noY BAM and FASTQ sidecars (default: yes)
  --star-bam-cbub yes|no     Emit barcode/UMI tags in BAMs; yes includes final CB/UB and triggers tagged replay (default: no)
  --star-bam-gxgn yes|no     Emit GX/GN gene tags in BAMs (default: no)
  --star-solo-features VALUE STAR solo features (default: "GeneFull Velocyto")
  --star-compare-feature VAL STAR feature to materialize for CR comparison (default: GeneFull)
  -h, --help                 Show help

Environment overrides:
  RAW_DIR CONFIG SAMPLE_ID SAMPLE_STEM LANES STAR_BIN GENOME_DIR STAR_OCM_BARCODE_MODE
  STAR_NATIVE_OCM_BAM_SPLIT STAR_NATIVE_OCM_MEX STAR_YREMOVE STAR_BAM_CBUB_TAGS STAR_BAM_GXGN_TAGS
  STAR_SOLO_FEATURES STAR_COMPARE_FEATURE SOLO_CB_WHITELIST CR_REFERENCE CELLRANGER_BIN CR_LOCALMEM CR_CREATE_BAM
  CR_REUSE_RUN_DIR CR_REUSE_SEARCH_ROOT OCM_PREP_REUSE_ROOT SIMPLEED_SIM_N OCM_MATERIALIZE_THREADS

Reuse behavior:
  OCM_PREP_REUSE_ROOT=/path Reuse staged FASTQs from a previous prepared harness
  STAR_NATIVE_OCM_MEX=auto   Skip Python GeneFull/Velocyto production materialization when native STAR OCM outputs exist
  CR_REUSE_RUN_DIR=auto      Find a completed CR run for this sample/read count (default)
  CR_REUSE_RUN_DIR=/path     Use that completed CR run for comparison
  --run-cr                   Force a new CR run for intentional reference generation
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
    --read-pairs) READ_PAIRS="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --run-all)
      RUN_ALL_REQUESTED=1
      RUN_PREPARE=1
      RUN_STAR=1
      if [[ "${STAR_NATIVE_OCM_BAM_SPLIT}" == "yes" || "${STAR_NATIVE_OCM_BAM_SPLIT}" == "auto" ]]; then
        RUN_SPLIT_BAM=0
      else
        RUN_SPLIT_BAM=1
      fi
      RUN_MATERIALIZE=1
      RUN_CR=0
      RUN_MERGE_CR_BAM=0
      RUN_COMPARE=1
      shift
      ;;
    --prepare) RUN_PREPARE=1; shift ;;
    --run-star) RUN_STAR=1; shift ;;
    --split-bam) RUN_SPLIT_BAM=1; shift ;;
    --materialize) RUN_MATERIALIZE=1; shift ;;
    --run-cr) RUN_CR=1; RUN_CR_EXPLICIT=1; CR_REUSE_RUN_DIR=""; shift ;;
    --merge-cr-bam) RUN_MERGE_CR_BAM=1; shift ;;
    --compare) RUN_COMPARE=1; shift ;;
    --force) FORCE=1; shift ;;
    --star-out-samtype) STAR_OUTSAMTYPE="$2"; shift 2 ;;
    --star-yremove) STAR_YREMOVE="$2"; shift 2 ;;
    --star-bam-cbub) STAR_BAM_CBUB_TAGS="$2"; shift 2 ;;
    --star-bam-gxgn) STAR_BAM_GXGN_TAGS="$2"; shift 2 ;;
    --star-solo-features) STAR_SOLO_FEATURES="$2"; shift 2 ;;
    --star-compare-feature) STAR_COMPARE_FEATURE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

if [[ "${RUN_PREPARE}${RUN_STAR}${RUN_SPLIT_BAM}${RUN_MATERIALIZE}${RUN_CR}${RUN_MERGE_CR_BAM}${RUN_COMPARE}" == "0000000" ]]; then
  RUN_PREPARE=1
fi

[[ "${READ_PAIRS}" =~ ^[0-9]+$ ]] || die "--read-pairs must be an integer"
[[ "${THREADS}" =~ ^[0-9]+$ && "${THREADS}" -gt 0 ]] || die "--threads must be positive"
[[ "${OCM_MATERIALIZE_THREADS}" =~ ^[0-9]+$ && "${OCM_MATERIALIZE_THREADS}" -gt 0 ]] || die "OCM_MATERIALIZE_THREADS must be positive"
[[ -d "${RAW_DIR}" ]] || die "Missing RAW_DIR: ${RAW_DIR}"
[[ -f "${CONFIG}" ]] || die "Missing Cell Ranger multi config: ${CONFIG}"
[[ -x "${STAR_BIN}" ]] || die "Missing STAR binary: ${STAR_BIN}"
[[ -d "${GENOME_DIR}" ]] || die "Missing STAR genomeDir: ${GENOME_DIR}"
[[ -f "${SOLO_CB_WHITELIST}" ]] || die "Missing STAR whitelist: ${SOLO_CB_WHITELIST}"
case "${STAR_OCM_BARCODE_MODE}" in
  flex|legacy17|posthoc) ;;
  *) die "STAR_OCM_BARCODE_MODE must be flex, legacy17, or posthoc" ;;
esac
case "${STAR_NATIVE_OCM_BAM_SPLIT}" in
  yes|no|auto) ;;
  *) die "STAR_NATIVE_OCM_BAM_SPLIT must be yes, no, or auto" ;;
esac
case "${STAR_NATIVE_OCM_MEX}" in
  yes|no|auto) ;;
  *) die "STAR_NATIVE_OCM_MEX must be yes, no, or auto" ;;
esac
case "${STAR_YREMOVE}" in
  yes|no) ;;
  *) die "STAR_YREMOVE must be yes or no" ;;
esac
case "${STAR_BAM_CBUB_TAGS}" in
  yes|no) ;;
  *) die "STAR_BAM_CBUB_TAGS must be yes or no" ;;
esac
case "${STAR_BAM_GXGN_TAGS}" in
  yes|no) ;;
  *) die "STAR_BAM_GXGN_TAGS must be yes or no" ;;
esac
[[ -n "${STAR_SOLO_FEATURES//[[:space:]]/}" ]] || die "STAR_SOLO_FEATURES must contain at least one feature"
[[ -n "${STAR_COMPARE_FEATURE//[[:space:]]/}" ]] || die "STAR_COMPARE_FEATURE must not be empty"

OUT_ROOT="$(realpath -m "${OUT_ROOT}")"
STAGE_ROOT="${OUT_ROOT}/stage"
CR_FASTQS="${STAGE_ROOT}/cr_fastqs"
STAR_FASTQS="${STAGE_ROOT}/star_composite_fastqs"
COMPOSITE_WHITELIST="${STAGE_ROOT}/3M-3pgex-may-2023_TRU_OCM17.txt"
MANIFEST="${OUT_ROOT}/downsample_manifest.tsv"
PREP_STATS="${OUT_ROOT}/prepare_stats.json"
LOG_DIR="${OUT_ROOT}/logs"
STAR_RUN_DIR="${OUT_ROOT}/star_composite/run"
STAR_TMP_DIR="${OUT_ROOT}/star_composite/tmp"
STAR_SCRIPT="${OUT_ROOT}/RUN_STAR_COMPOSITE.sh"
CR_PARENT="${OUT_ROOT}/cellranger"
CR_ID="${SAMPLE_ID}_ocm_composite_${READ_PAIRS}"
CR_COMPARE_RUN_DIR="${CR_PARENT}/${CR_ID}"
CR_CONFIG="${OUT_ROOT}/cellranger_multi_config.csv"
CR_SCRIPT="${OUT_ROOT}/RUN_CELLRANGER_MULTI.sh"
MAT_COMPARE="${OUT_ROOT}/star_materialized/${STAR_COMPARE_FEATURE}_compare"
MAT_GENEFULL="${OUT_ROOT}/star_materialized/GeneFull"
STAR_NATIVE_OUTS="${OUT_ROOT}/star_composite/outs"
STAR_NATIVE_SAMPLES="${OUT_ROOT}/star_composite/samples"
STAR_BAM_SPLIT="${OUT_ROOT}/star_materialized/bam"
CR_MERGED_BAM_DIR="${OUT_ROOT}/cellranger_merged_bam"
CR_MERGED_BAM="${CR_MERGED_BAM_DIR}/merged_alignments.bam"
COMPARE_JSON="${OUT_ROOT}/parity_gene_vs_cr9.json"
COMPARE_TSV="${OUT_ROOT}/parity_gene_vs_cr9.tsv"

mkdir -p "${LOG_DIR}" "${STAGE_ROOT}"
FORCE_ARGS=()
if [[ "${FORCE}" == "1" ]]; then
  FORCE_ARGS=(--force)
fi

link_or_copy_file() {
  local src="$1"
  local dst="$2"
  mkdir -p "$(dirname "${dst}")"
  rm -f "${dst}"
  ln "${src}" "${dst}" 2>/dev/null || cp -p "${src}" "${dst}"
}

link_or_copy_dir_files() {
  local src_dir="$1"
  local dst_dir="$2"
  [[ -d "${src_dir}" ]] || die "Missing reuse source directory: ${src_dir}"
  rm -rf "${dst_dir}"
  mkdir -p "${dst_dir}"
  local src
  while IFS= read -r -d '' src; do
    link_or_copy_file "${src}" "${dst_dir}/$(basename "${src}")"
  done < <(find "${src_dir}" -maxdepth 1 -type f -print0)
}

is_prepared_harness() {
  local root="$1"
  [[ -f "${root}/prepare_stats.json" ]] || return 1
  [[ -d "${root}/stage/cr_fastqs" ]] || return 1
  [[ -d "${root}/stage/star_composite_fastqs" ]] || return 1
  [[ -f "${root}/RUN_STAR_COMPOSITE.sh" ]] || return 1
  [[ -f "${root}/cellranger_multi_config.csv" ]] || return 1
  python3 - "$root/prepare_stats.json" "$SAMPLE_ID" "$SAMPLE_STEM" "$READ_PAIRS" "$STAR_OCM_BARCODE_MODE" <<'PY'
import json
import sys

path, sample_id, sample_stem, read_pairs, barcode_mode = sys.argv[1:]
with open(path, "r", encoding="utf-8") as handle:
    stats = json.load(handle)
expected_star_mode = "native" if barcode_mode == "flex" else "legacy17"
ok = (
    stats.get("sample_id") == sample_id
    and stats.get("sample_stem") == sample_stem
    and int(stats.get("read_pairs", -1)) == int(read_pairs)
    and stats.get("star_mode") == expected_star_mode
)
raise SystemExit(0 if ok else 1)
PY
}

reuse_prepared_harness() {
  local reuse_root="$1"
  reuse_root="$(realpath -m "${reuse_root}")"
  is_prepared_harness "${reuse_root}" || die "OCM_PREP_REUSE_ROOT is not a matching prepared harness: ${reuse_root}"
  log "Reusing prepared downsample artifacts from ${reuse_root}"
  link_or_copy_dir_files "${reuse_root}/stage/cr_fastqs" "${CR_FASTQS}"
  link_or_copy_dir_files "${reuse_root}/stage/star_composite_fastqs" "${STAR_FASTQS}"
  [[ -f "${reuse_root}/downsample_manifest.tsv" ]] && link_or_copy_file "${reuse_root}/downsample_manifest.tsv" "${MANIFEST}"
  [[ -f "${reuse_root}/prepare_stats.json" ]] && link_or_copy_file "${reuse_root}/prepare_stats.json" "${PREP_STATS}"
  if [[ -f "${reuse_root}/stage/$(basename "${COMPOSITE_WHITELIST}")" ]]; then
    link_or_copy_file "${reuse_root}/stage/$(basename "${COMPOSITE_WHITELIST}")" "${COMPOSITE_WHITELIST}"
  fi
  printf 'reused_from=%s\nreused_utc=%s\n' \
    "${reuse_root}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "${OUT_ROOT}/PREP_REUSED.txt"
}

is_completed_cr_run() {
  local run_dir="$1"
  [[ -d "${run_dir}/outs" ]] || return 1
  [[ -d "${run_dir}/outs/per_sample_outs" ]] || return 1
  find "${run_dir}/outs/per_sample_outs" \
    -path '*/count/sample_filtered_feature_bc_matrix/matrix.mtx.gz' \
    -type f -print -quit | grep -q .
}

find_reusable_cr_run() {
  local current_cr_dir
  current_cr_dir="$(realpath -m "${CR_PARENT}/${CR_ID}")"
  [[ -d "${CR_REUSE_SEARCH_ROOT}" ]] || return 1
  local candidate
  while IFS= read -r candidate; do
    candidate="$(realpath -m "${candidate}")"
    [[ "${candidate}" == "${current_cr_dir}" ]] && continue
    if is_completed_cr_run "${candidate}"; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done < <(find "${CR_REUSE_SEARCH_ROOT}" -maxdepth 4 -path "*/cellranger/${CR_ID}" -type d | sort -r)
  return 1
}

if [[ "${RUN_CR_EXPLICIT}" != "1" && "${CR_REUSE_RUN_DIR}" == "auto" && \
      ( "${RUN_COMPARE}" == "1" || "${RUN_CR}" == "1" || "${RUN_ALL_REQUESTED}" == "1" ) ]]; then
  if is_completed_cr_run "${CR_PARENT}/${CR_ID}"; then
    CR_REUSE_RUN_DIR=""
    RUN_CR=0
    log "Reusing completed Cell Ranger reference in current output: ${CR_PARENT}/${CR_ID}"
  elif CR_AUTO_DIR="$(find_reusable_cr_run)"; then
    CR_REUSE_RUN_DIR="${CR_AUTO_DIR}"
    log "Reusing completed Cell Ranger reference: ${CR_REUSE_RUN_DIR}"
  else
    CR_REUSE_RUN_DIR=""
    if [[ "${RUN_COMPARE}" == "1" ]]; then
      die "No completed CR reference found for ${CR_ID}; set CR_REUSE_RUN_DIR or run --run-cr intentionally"
    fi
  fi
fi

if [[ -n "${CR_REUSE_RUN_DIR}" && "${CR_REUSE_RUN_DIR}" != "auto" ]]; then
  CR_REUSE_RUN_DIR="$(realpath -m "${CR_REUSE_RUN_DIR}")"
  is_completed_cr_run "${CR_REUSE_RUN_DIR}" || die "CR_REUSE_RUN_DIR must point to a completed Cell Ranger run with per-sample filtered MEX: ${CR_REUSE_RUN_DIR}"
  RUN_CR=0
  RUN_MERGE_CR_BAM=0
fi
if [[ "${RUN_CR}" == "1" ]]; then
  [[ -d "${CR_REFERENCE}" ]] || die "Missing Cell Ranger reference: ${CR_REFERENCE}"
  command -v "${CELLRANGER_BIN}" >/dev/null 2>&1 || die "Missing cellranger binary: ${CELLRANGER_BIN}"
fi
CR_COMPARE_RUN_DIR="${CR_PARENT}/${CR_ID}"
if [[ -n "${CR_REUSE_RUN_DIR}" && "${CR_REUSE_RUN_DIR}" != "auto" ]]; then
  CR_COMPARE_RUN_DIR="${CR_REUSE_RUN_DIR}"
fi

join_by_comma() {
  local IFS=,
  printf '%s' "$*"
}

source_fastq_name() {
  local lane="$1"
  local read="$2"
  printf '%s_%s_%s_001.fastq.gz' "${SAMPLE_STEM}" "${lane}" "${read}"
}

render_star_script() {
  local -a lane_array r1_files r2_files
  IFS=',' read -r -a lane_array <<< "${LANES}"
  local star_cb_len=17
  local star_umi_start=18
  local star_whitelist="${COMPOSITE_WHITELIST}"
  local star_ocm_param="posthoc"
  local star_inline_cb_correction="no"
  if [[ "${STAR_OCM_BARCODE_MODE}" == "flex" ]]; then
    star_cb_len=16
    star_umi_start=17
    star_whitelist="${SOLO_CB_WHITELIST}"
    star_ocm_param="flex"
    star_inline_cb_correction="yes"
  fi
  local lane
  for lane in "${lane_array[@]}"; do
    r1_files+=("${STAR_FASTQS}/$(source_fastq_name "${lane}" R1)")
    r2_files+=("${STAR_FASTQS}/$(source_fastq_name "${lane}" R2)")
  done
  local r1_csv r2_csv
  r1_csv="$(join_by_comma "${r1_files[@]}")"
  r2_csv="$(join_by_comma "${r2_files[@]}")"

  local -a out_sam_args
  read -r -a out_sam_args <<< "${STAR_OUTSAMTYPE}"

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
    printf 'mkdir -p %q %q\n' "${STAR_RUN_DIR}" "${LOG_DIR}"
    printf 'cmd=(\n'
    printf '  %q\n' "${STAR_BIN}"
    printf '  --runThreadN %q\n' "${THREADS}"
    printf '  --dynamicThreadInterface 1\n'
    printf '  --genomeDir %q\n' "${GENOME_DIR}"
    printf '  --readFilesIn %q %q\n' "${r2_csv}" "${r1_csv}"
    printf '  --readFilesCommand zcat\n'
    printf '  --outFileNamePrefix %q\n' "${STAR_RUN_DIR}/"
    printf '  --outTmpDir %q\n' "${STAR_TMP_DIR}"
    printf '  --outSAMtype'
    local arg
    for arg in "${out_sam_args[@]}"; do
      printf ' %q' "${arg}"
    done
    printf '\n'
    if [[ "${STAR_YREMOVE}" == "yes" ]]; then
      printf '  --emitNoYBAM yes\n'
      printf '  --emitYNoYFastq yes\n'
      printf '  --emitYNoYFastqCompression gz\n'
    fi
    local star_sam_attrs=(NH HI AS nM NM)
    if [[ "${STAR_BAM_CBUB_TAGS}" == "yes" ]]; then
      star_sam_attrs+=(CB UB CR UR)
    fi
    if [[ "${STAR_BAM_GXGN_TAGS}" == "yes" ]]; then
      star_sam_attrs+=(GX GN)
    fi
    printf '  --outSAMattributes'
    for arg in "${star_sam_attrs[@]}"; do
      printf ' %q' "${arg}"
    done
    printf '\n'
    printf '  --clipAdapterType CellRanger4\n'
    printf '  --clip3pPolyG yes\n'
    printf '  --alignEndsType Local\n'
    printf '  --chimSegmentMin 1000000\n'
    printf '  --soloType CB_UMI_Simple\n'
    printf '  --soloCBstart 1\n'
    printf '  --soloCBlen %q\n' "${star_cb_len}"
    printf '  --soloUMIstart %q\n' "${star_umi_start}"
    printf '  --soloUMIlen 12\n'
    printf '  --soloBarcodeReadLength 0\n'
    printf '  --soloCBwhitelist %q\n' "${star_whitelist}"
    printf '  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\n'
    printf '  --soloInlineCBCorrection %q\n' "${star_inline_cb_correction}"
    printf '  --soloUMIfiltering MultiGeneUMI_CR\n'
    printf '  --soloUMIdedup 1MM_CR\n'
    printf '  --soloMultiMappers Unique\n'
    printf '  --soloCellFilter None\n'
    printf '  --soloCbUbRequireTogether no\n'
    printf '  --soloStrand Forward\n'
    local solo_features=()
    read -r -a solo_features <<< "${STAR_SOLO_FEATURES}"
    printf '  --soloFeatures'
    for arg in "${solo_features[@]}"; do
      printf ' %q' "${arg}"
    done
    printf '\n'
    printf '  --soloCrGexFeature genefull\n'
    printf '  --soloCrMultimapRescue yes\n'
    printf '  --soloInlineHashMode no\n'
    printf '  --ocmMultiEnable auto\n'
    printf '  --ocmMultiConfig %q\n' "${CONFIG}"
    printf '  --ocmMultiBarcodeMode %q\n' "${star_ocm_param}"
    printf '  --ocmMultiBamSplit %q\n' "${STAR_NATIVE_OCM_BAM_SPLIT}"
    printf '  --ocmMultiOutputCompat cellranger\n'
    printf ')\n'
    printf 'printf "started_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/STAR_STARTED.txt"
    printf 'printf "STAR command:" | tee %q\n' "${LOG_DIR}/star.command.txt"
    printf 'printf " %%q" "${cmd[@]}" | tee -a %q\n' "${LOG_DIR}/star.command.txt"
    printf 'printf "\\n" | tee -a %q\n' "${LOG_DIR}/star.command.txt"
    printf '"${cmd[@]}" 2>&1 | tee %q\n' "${LOG_DIR}/star.log"
    printf 'rm -rf %q\n' "${STAR_TMP_DIR}"
    printf 'printf "completed_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/STAR_COMPLETED.txt"
  } > "${STAR_SCRIPT}"
  chmod +x "${STAR_SCRIPT}"
}

render_cr_config() {
  cat > "${CR_CONFIG}" <<EOF
[gene-expression]
reference,${CR_REFERENCE}
create-bam,${CR_CREATE_BAM}
include-introns,false

[libraries]
fastq_id,fastqs,feature_types
25E32-L3_GT25-03394_ACCTCGAGCT-ATCGAACACA,${CR_FASTQS},Gene Expression

[samples]
sample_id,ocm_barcode_ids,description
GCM1-Day-4,OB1,iPSCs
GRHL1-Day-4,OB2,iPSCs
OVOL1-Day-4,OB3,iPSCs
WT-PrS-20pct-Day-4,OB4,iPSCs
EOF
}

render_cr_script() {
  mkdir -p "${CR_PARENT}"
  {
    printf '#!/usr/bin/env bash\n'
    printf 'set -euo pipefail\n\n'
    printf 'cd %q\n' "${CR_PARENT}"
    printf 'rm -rf %q\n' "${CR_PARENT}/${CR_ID}"
    printf 'printf "started_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/CELLRANGER_STARTED.txt"
    printf 'cmd=(%q multi --id=%q --csv=%q --localcores=%q --localmem=%q)\n' "${CELLRANGER_BIN}" "${CR_ID}" "${CR_CONFIG}" "${THREADS}" "${CR_LOCALMEM}"
    printf 'printf "Cell Ranger command:" | tee %q\n' "${LOG_DIR}/cellranger.command.txt"
    printf 'printf " %%q" "${cmd[@]}" | tee -a %q\n' "${LOG_DIR}/cellranger.command.txt"
    printf 'printf "\\n" | tee -a %q\n' "${LOG_DIR}/cellranger.command.txt"
    printf '"${cmd[@]}" 2>&1 | tee %q\n' "${LOG_DIR}/cellranger.log"
    printf 'printf "completed_utc=%%s\\n" "$(date -u +%%Y-%%m-%%dT%%H:%%M:%%SZ)" > %q\n' "${OUT_ROOT}/CELLRANGER_COMPLETED.txt"
  } > "${CR_SCRIPT}"
  chmod +x "${CR_SCRIPT}"
}

write_readme() {
  cat > "${OUT_ROOT}/README.txt" <<EOF
JAX scRNAseq02 OCM composite-barcode smoke

sample_id=${SAMPLE_ID}
sample_stem=${SAMPLE_STEM}
read_pairs=${READ_PAIRS}
lanes=${LANES}
out_root=${OUT_ROOT}
star_ocm_barcode_mode=${STAR_OCM_BARCODE_MODE}
star_native_ocm_bam_split=${STAR_NATIVE_OCM_BAM_SPLIT}
star_native_ocm_mex=${STAR_NATIVE_OCM_MEX}
star_yremove=${STAR_YREMOVE}
star_bam_cbub_tags=${STAR_BAM_CBUB_TAGS}
star_bam_gxgn_tags=${STAR_BAM_GXGN_TAGS}
star_solo_features=${STAR_SOLO_FEATURES}
star_compare_feature=${STAR_COMPARE_FEATURE}
cr_reuse_run_dir=${CR_REUSE_RUN_DIR}

STAR composite run:
  ${STAR_SCRIPT}

Cell Ranger 9 multi control:
  ${CR_SCRIPT}
  compare source: ${CR_COMPARE_RUN_DIR}

Materialized STAR outputs:
  ${MAT_COMPARE}
  ${MAT_GENEFULL}
  native STAR OCM outputs: ${STAR_NATIVE_OUTS}, ${STAR_NATIVE_SAMPLES}

STAR BAM outputs:
  native OCM BAMs when enabled: ${OUT_ROOT}/star_composite/outs/per_sample_outs/*/count/sample_alignments.bam
  Y/noY BAMs and FASTQs when STAR_YREMOVE=yes: ${STAR_RUN_DIR}/Aligned.out_{Y,noY}.bam and ${STAR_RUN_DIR}/y_separated/
  post-split BAMs when --split-bam is used: ${STAR_BAM_SPLIT}

Cell Ranger BAM:
  per-sample BAMs under ${CR_PARENT}/${CR_ID}/outs/per_sample_outs/*/count/
  merged BAM: ${CR_MERGED_BAM}

Parity report:
  ${COMPARE_TSV}
EOF
}

if [[ "${RUN_PREPARE}" == "1" ]]; then
  PREP_STAR_MODE="legacy17"
  if [[ "${STAR_OCM_BARCODE_MODE}" == "flex" ]]; then
    PREP_STAR_MODE="native"
  fi
  if [[ "${FORCE}" != "1" && -z "${OCM_PREP_REUSE_ROOT}" ]] && is_prepared_harness "${OUT_ROOT}"; then
    log "Prepared harness already exists; reusing ${OUT_ROOT}"
  else
    if [[ -n "${OCM_PREP_REUSE_ROOT}" ]]; then
      reuse_prepared_harness "${OCM_PREP_REUSE_ROOT}"
    else
      log "Preparing ${READ_PAIRS} read-pair OCM composite downsample"
      python3 "${SCRIPT_DIR}/ocm_composite_adapter.py" prepare \
        --raw-dir "${RAW_DIR}" \
        --sample-id "${SAMPLE_ID}" \
        --sample-stem "${SAMPLE_STEM}" \
        --lanes "${LANES}" \
        --read-pairs "${READ_PAIRS}" \
        --cr-fastq-dir "${CR_FASTQS}" \
        --star-fastq-dir "${STAR_FASTQS}" \
        --star-mode "${PREP_STAR_MODE}" \
        --whitelist "${SOLO_CB_WHITELIST}" \
        --composite-whitelist "${COMPOSITE_WHITELIST}" \
        --manifest "${MANIFEST}" \
        --stats-json "${PREP_STATS}" \
        "${FORCE_ARGS[@]}"
    fi
  fi
  render_star_script
  render_cr_config
  render_cr_script
  write_readme
  log "Prepared harness: ${OUT_ROOT}"
fi

if [[ "${RUN_STAR}" == "1" ]]; then
  if [[ "${FORCE}" != "1" && -f "${OUT_ROOT}/STAR_COMPLETED.txt" && -d "${STAR_RUN_DIR}/Solo.out" ]]; then
    log "STAR output already complete; reusing ${STAR_RUN_DIR}"
  else
    log "Running STAR composite smoke"
    "${STAR_SCRIPT}"
  fi
fi

if [[ "${RUN_SPLIT_BAM}" == "1" ]]; then
  if [[ "${FORCE}" != "1" && -f "${STAR_BAM_SPLIT}/split_bam_summary.json" ]]; then
    log "STAR post-split BAMs already exist; reusing ${STAR_BAM_SPLIT}"
  else
    log "Splitting STAR composite BAM by OCM"
    python3 "${SCRIPT_DIR}/ocm_composite_adapter.py" split-bam \
      --bam "${STAR_RUN_DIR}/Aligned.out.bam" \
      --config "${CONFIG}" \
      --out-dir "${STAR_BAM_SPLIT}" \
      "${FORCE_ARGS[@]}"
  fi
fi

if [[ "${RUN_MATERIALIZE}" == "1" ]]; then
  if [[ "${RUN_COMPARE}" != "1" ]]; then
    log "Skipping STAR ${STAR_COMPARE_FEATURE} comparator materialization because --compare was not requested"
  elif [[ "${FORCE}" != "1" && -f "${MAT_COMPARE}/materialization_summary.json" ]]; then
    log "STAR ${STAR_COMPARE_FEATURE} materialization already exists; reusing ${MAT_COMPARE}"
  else
    log "Materializing STAR ${STAR_COMPARE_FEATURE} outputs for CR parity"
    python3 "${SCRIPT_DIR}/ocm_composite_adapter.py" materialize \
      --repo-root "${REPO_ROOT}" \
      --star-run-dir "${STAR_RUN_DIR}" \
      --config "${CONFIG}" \
      --feature "${STAR_COMPARE_FEATURE}" \
      --out-dir "${MAT_COMPARE}" \
      --simpleed-mode "${SIMPLEED_MODE}" \
      --sim-n "${SIMPLEED_SIM_N}" \
      --threads "${OCM_MATERIALIZE_THREADS}" \
      "${FORCE_ARGS[@]}"
  fi

  native_mex_ready=0
  if [[ -d "${STAR_NATIVE_OUTS}/per_sample_outs" && -d "${STAR_NATIVE_SAMPLES}" ]]; then
    native_mex_ready=1
  fi
  if [[ "${STAR_NATIVE_OCM_MEX}" == "yes" || ( "${STAR_NATIVE_OCM_MEX}" == "auto" && "${native_mex_ready}" == "1" ) ]]; then
    log "Skipping Python GeneFull/Velocyto production materialization; using native STAR OCM outputs at ${STAR_NATIVE_OUTS}"
  elif [[ "${FORCE}" != "1" && -f "${MAT_GENEFULL}/materialization_summary.json" ]]; then
    log "STAR GeneFull/Velocyto materialization already exists; reusing ${MAT_GENEFULL}"
  else
    log "Materializing STAR GeneFull outputs with Velocyto using optimized production router"
    python3 "${SCRIPT_DIR}/ocm_composite_adapter.py" materialize-production \
      --repo-root "${REPO_ROOT}" \
      --star-run-dir "${STAR_RUN_DIR}" \
      --config "${CONFIG}" \
      --out-dir "${MAT_GENEFULL}" \
      --simpleed-mode "${SIMPLEED_MODE}" \
      --sim-n "${SIMPLEED_SIM_N}" \
      --threads "${OCM_MATERIALIZE_THREADS}" \
      "${FORCE_ARGS[@]}"
  fi
fi

if [[ "${RUN_CR}" == "1" ]]; then
  if [[ "${FORCE}" != "1" && -f "${OUT_ROOT}/CELLRANGER_COMPLETED.txt" ]] && is_completed_cr_run "${CR_PARENT}/${CR_ID}"; then
    log "Cell Ranger output already complete; reusing ${CR_PARENT}/${CR_ID}"
  else
    log "Running Cell Ranger 9 multi control"
    "${CR_SCRIPT}"
  fi
fi

if [[ "${RUN_MERGE_CR_BAM}" == "1" ]]; then
  mkdir -p "${CR_MERGED_BAM_DIR}"
  if [[ "${FORCE}" != "1" && -f "${CR_MERGED_BAM}" ]] && samtools quickcheck "${CR_MERGED_BAM}" 2>/dev/null; then
    log "Merged Cell Ranger BAM already exists; reusing ${CR_MERGED_BAM}"
    if [[ ! -f "${CR_MERGED_BAM}.bai" ]]; then
      log "Indexing existing merged Cell Ranger BAM"
      samtools index -@ "${THREADS}" "${CR_MERGED_BAM}"
    fi
  else
    rm -f "${CR_MERGED_BAM}" "${CR_MERGED_BAM}.bai"
    log "Merging Cell Ranger per-sample and unassigned BAMs"
    mapfile -t CR_BAMS < <(
      find "${CR_PARENT}/${CR_ID}/outs/per_sample_outs" -path '*/count/sample_alignments.bam' -type f | sort
      find "${CR_PARENT}/${CR_ID}/outs/multi/count" -name 'unassigned_alignments.bam' -type f | sort
    )
    [[ "${#CR_BAMS[@]}" -gt 0 ]] || die "No Cell Ranger BAMs found to merge"
    samtools merge -f -@ "${THREADS}" "${CR_MERGED_BAM}" "${CR_BAMS[@]}"
    samtools index -@ "${THREADS}" "${CR_MERGED_BAM}"
  fi
fi

if [[ "${RUN_COMPARE}" == "1" ]]; then
  log "Comparing STAR ${STAR_COMPARE_FEATURE} materialization to CR9 multi"
  python3 "${SCRIPT_DIR}/ocm_composite_adapter.py" compare \
    --star-materialized-dir "${MAT_COMPARE}" \
    --cr-run-dir "${CR_COMPARE_RUN_DIR}" \
    --config "${CONFIG}" \
    --out-json "${COMPARE_JSON}" \
    --out-tsv "${COMPARE_TSV}"
fi

log "Done: ${OUT_ROOT}"
