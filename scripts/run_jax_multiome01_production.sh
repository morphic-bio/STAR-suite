#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SMOKE_RUNNER="${REPO_ROOT}/scripts/run_star_multiome_lane_smoke.sh"
REMOTE_POST_MEX="${REPO_ROOT}/scripts/run_remote_multiome_post_mex_rsync.sh"
GLOBUS_UPLOADER="${REPO_ROOT}/scripts/upload_jax_multiome01_large_files_globus.sh"

RAW_DIR="/mnt/pikachu/JAX_Multiome01/raw"
METADATA_XLSX="/mnt/pikachu/DPC_metadata_template_Multiome1-complete.xlsx"
OUTPUT_ROOT="/mnt/pikachu/JAX_Multiome01_processed/star_multiome_$(date -u +%Y%m%dT%H%M%SZ)"
REMOTE_HOST="${MULTIOME_REMOTE_HOST:-10.159.4.53}"
REMOTE_ROOT="${MULTIOME_REMOTE_ROOT:-/home/lhhung/jax_multiome_remote_downstream_production}"
THREADS="${STAR_MULTIOME_THREADS:-16}"
CHROMAP_THREADS="${STAR_MULTIOME_CHROMAP_THREADS:-16}"
CHROMAP_LOW_MEM="${STAR_MULTIOME_CHROMAP_LOW_MEM:-1}"
CHROMAP_LOW_MEM_RAM="${STAR_MULTIOME_CHROMAP_LOW_MEM_RAM:-0}"
CHROMAP_MACS3_FRAG_LOW_MEM="${STAR_MULTIOME_CHROMAP_MACS3_FRAG_LOW_MEM:-1}"
CHROMAP_START_MODE="${STAR_MULTIOME_CHROMAP_START_MODE:-concurrent}"
CELLBENDER_CPU_CORES="${MULTIOME_CELLBENDER_CPU_CORES:-24}"
CELLBENDER_GPU="1"
NO_SYNC_IMAGES="0"
KEEP_REMOTE="0"
FORCE="0"
SKIP_BUILD="0"
MANIFEST_ONLY="0"
GLOBUS_UPLOAD_LARGE_FILES="${JAX_MULTIOME_GLOBUS_UPLOAD_LARGE_FILES:-0}"
GLOBUS_SOURCE_ENDPOINT="${JAX_MULTIOME_GLOBUS_SOURCE_ENDPOINT:-07446cad-33b8-11f0-8c0c-0afffb017b7d}"
GLOBUS_DEST_ENDPOINT="${JAX_MULTIOME_GLOBUS_DEST_ENDPOINT:-61fb8b9a-9b52-456e-928c-30c0fb0140bf}"
GLOBUS_DEST_ROOT="${JAX_MULTIOME_GLOBUS_DEST_ROOT:-/JAX_Multiome01_processed/large_files}"
START_AT=""
SAMPLES=""

usage() {
  cat <<'EOF'
Usage:
  run_jax_multiome01_production.sh [options]

Runs the 9-sample JAX_Multiome01 STAR/Chromap production workflow from the
metadata workbook. Each sample is processed independently with multi-lane FASTQ
CSVs, native ATAC barcode parsing, GeneFull+Velocyto RNA, Y/noY outputs, local
ATAC sidecar/peak-MEX materialization, and remote post-MEX CellBender/MuData
validation. The local STAR/Chromap phase stops at the MEX/sidecar boundary so
the next sample can start while remote post-MEX work runs.

Options:
  --raw-dir PATH
  --metadata-xlsx PATH
  --output-root PATH
  --remote-host HOST
  --remote-root PATH
  --threads N
  --chromap-threads N
  --chromap-low-mem
  --chromap-low-mem-ram N
  --chromap-macs3-frag-low-mem
  --chromap-start-mode MODE
                           STAR/Chromap scheduling: postMapping or concurrent
                           (default: concurrent)
  --cellbender-cpu-cores N
  --no-cellbender-gpu
  --no-sync-images
  --keep-remote
  --skip-build
  --manifest-only
  --dry-run                 Alias for --manifest-only
  --globus-upload-large-files
  --globus-source-endpoint UUID
  --globus-dest-endpoint UUID
  --globus-dest-root PATH
  --force
  --start-at SAMPLE_LABEL_OR_SLUG
  --samples CSV_LABELS_OR_SLUGS
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw-dir) RAW_DIR="$2"; shift 2 ;;
    --metadata-xlsx) METADATA_XLSX="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --remote-host) REMOTE_HOST="$2"; shift 2 ;;
    --remote-root) REMOTE_ROOT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --chromap-threads) CHROMAP_THREADS="$2"; shift 2 ;;
    --chromap-low-mem) CHROMAP_LOW_MEM="1"; shift ;;
    --chromap-low-mem-ram) CHROMAP_LOW_MEM_RAM="$2"; shift 2 ;;
    --chromap-macs3-frag-low-mem) CHROMAP_MACS3_FRAG_LOW_MEM="1"; shift ;;
    --chromap-start-mode) CHROMAP_START_MODE="$2"; shift 2 ;;
    --cellbender-cpu-cores) CELLBENDER_CPU_CORES="$2"; shift 2 ;;
    --no-cellbender-gpu) CELLBENDER_GPU="0"; shift ;;
    --no-sync-images) NO_SYNC_IMAGES="1"; shift ;;
    --keep-remote) KEEP_REMOTE="1"; shift ;;
    --skip-build) SKIP_BUILD="1"; shift ;;
    --manifest-only|--dry-run) MANIFEST_ONLY="1"; shift ;;
    --globus-upload-large-files) GLOBUS_UPLOAD_LARGE_FILES="1"; shift ;;
    --globus-source-endpoint) GLOBUS_SOURCE_ENDPOINT="$2"; shift 2 ;;
    --globus-dest-endpoint) GLOBUS_DEST_ENDPOINT="$2"; shift 2 ;;
    --globus-dest-root) GLOBUS_DEST_ROOT="$2"; shift 2 ;;
    --force) FORCE="1"; shift ;;
    --start-at) START_AT="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

RAW_DIR="$(realpath "${RAW_DIR}")"
METADATA_XLSX="$(realpath "${METADATA_XLSX}")"
OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"

case "${CHROMAP_START_MODE}" in
  postMapping|concurrent) ;;
  *) echo "ERROR: --chromap-start-mode must be postMapping or concurrent" >&2; exit 1 ;;
esac

mkdir -p "${OUTPUT_ROOT}/metadata" "${OUTPUT_ROOT}/logs" "${OUTPUT_ROOT}/samples"
MANIFEST="${OUTPUT_ROOT}/metadata/sample_manifest.tsv"

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${OUTPUT_ROOT}/logs/production.log"
}

log "Writing production sample manifest"
python3 - "${METADATA_XLSX}" "${RAW_DIR}" > "${MANIFEST}" <<'PY'
import collections
import re
import sys
from pathlib import Path

import openpyxl

xlsx = Path(sys.argv[1])
raw = Path(sys.argv[2])
wb = openpyxl.load_workbook(xlsx, data_only=True)

def slug(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")

ws = wb["Library preparation"]
libs = []
for row in ws.iter_rows(min_row=6, values_only=True):
    if not any(v is not None and str(v).strip() for v in row):
        continue
    lib = row[0]
    sample = row[3]
    lib_type = row[4] or "GEX"
    libs.append((lib, sample, lib_type))

by_sample = collections.defaultdict(dict)
sample_order = []
for lib, sample, lib_type in libs:
    by_sample[sample][lib_type] = lib
    if lib_type == "ATAC":
        sample_order.append(sample)

ws = wb["Sequence file"]
headers = [cell.value for cell in ws[4]]
seq = collections.defaultdict(lambda: collections.defaultdict(list))
run_ids = collections.defaultdict(set)
for row in ws.iter_rows(min_row=6, values_only=True):
    if not any(v is not None and str(v).strip() for v in row):
        continue
    record = dict(zip(headers, row))
    lib = record["library_preparation.label"]
    read = record["sequence_file.read_index"]
    label = record["sequence_file.label"]
    seq[lib][read].append(label)
    if record.get("sequence_file.run_id"):
        run_ids[lib].add(str(record["sequence_file.run_id"]))

def paths(lib, read):
    labels = sorted(seq[lib][read])
    if not labels:
        raise SystemExit(f"Missing metadata rows for {lib} {read}")
    missing = [label for label in labels if not (raw / label).exists()]
    if missing:
        raise SystemExit("Missing required FASTQs:\n" + "\n".join(missing))
    return ",".join(str(raw / label) for label in labels)

print("\t".join([
    "sample",
    "sample_slug",
    "atac_library",
    "gex_library",
    "gex_r1",
    "gex_r2",
    "atac_r1",
    "atac_barcode",
    "atac_r2",
    "gex_run_ids",
    "atac_run_ids",
]))
for sample in sample_order:
    atac_lib = by_sample[sample]["ATAC"]
    gex_lib = by_sample[sample]["GEX"]
    print("\t".join([
        sample,
        slug(sample),
        atac_lib,
        gex_lib,
        paths(gex_lib, "read1"),
        paths(gex_lib, "read2"),
        paths(atac_lib, "read1"),
        paths(atac_lib, "read2"),
        paths(atac_lib, "read3"),
        ",".join(sorted(run_ids[gex_lib])),
        ",".join(sorted(run_ids[atac_lib])),
    ]))
PY

sample_count="$(tail -n +2 "${MANIFEST}" | wc -l)"
log "Manifest contains ${sample_count} samples with all required workflow FASTQs present"
log "Production settings: STAR threads=${THREADS}, Chromap threads=${CHROMAP_THREADS}, Chromap low-mem=${CHROMAP_LOW_MEM}, Chromap low-mem RAM=${CHROMAP_LOW_MEM_RAM}, MACS3 fragment low-mem=${CHROMAP_MACS3_FRAG_LOW_MEM}"
if [[ "${MANIFEST_ONLY}" == "1" ]]; then
  log "Manifest-only mode complete"
  echo "Manifest: ${MANIFEST}"
  exit 0
fi

if [[ "${SKIP_BUILD}" != "1" ]]; then
  log "Clean rebuilding STAR with Chromap support once for production"
  make -C "${REPO_ROOT}/core/legacy/source" clean > "${OUTPUT_ROOT}/logs/build_clean.log" 2>&1
  make -C "${REPO_ROOT}/core/legacy/source" -j8 STAR WITH_CHROMAP=1 > "${OUTPUT_ROOT}/logs/build_star_with_chromap.log" 2>&1
  make -C "${REPO_ROOT}/core/features/libchromap_contract" star_multiome_atac_peak_mex > "${OUTPUT_ROOT}/logs/build_star_multiome_atac_peak_mex.log" 2>&1
fi

sample_selected() {
  local sample="$1"
  local slug="$2"
  [[ -z "${SAMPLES}" ]] && return 0
  local wanted
  IFS=',' read -r -a wanted <<< "${SAMPLES}"
  for item in "${wanted[@]}"; do
    [[ "${item}" == "${sample}" || "${item}" == "${slug}" ]] && return 0
  done
  return 1
}

start_reached="0"
[[ -z "${START_AT}" ]] && start_reached="1"
REMOTE_POST_MEX_LOCK="${OUTPUT_ROOT}/logs/remote_post_mex.lock"
post_mex_pids=()
post_mex_samples=()

while IFS=$'\t' read -r sample sample_slug atac_lib gex_lib gex_r1 gex_r2 atac_r1 atac_barcode atac_r2 gex_run_ids atac_run_ids; do
  if [[ "${start_reached}" != "1" ]]; then
    if [[ "${START_AT}" == "${sample}" || "${START_AT}" == "${sample_slug}" ]]; then
      start_reached="1"
    else
      log "Skipping ${sample} until --start-at ${START_AT}"
      continue
    fi
  fi
  if ! sample_selected "${sample}" "${sample_slug}"; then
    log "Skipping ${sample}; not in --samples"
    continue
  fi

  sample_out="${OUTPUT_ROOT}/samples/${sample_slug}"
  sample_log="${OUTPUT_ROOT}/logs/${sample_slug}.wrapper.log"
  if [[ "${FORCE}" != "1" && -f "${sample_out}/mudata/star_chromap_filtered_multiome.h5mu" ]]; then
    log "Reusing completed sample ${sample}"
    continue
  fi

  log "Running ${sample} (ATAC ${atac_lib}; GEX ${gex_lib})"
  args=(
    "${SMOKE_RUNNER}"
    --gex-r1 "${gex_r1}"
    --gex-r2 "${gex_r2}"
    --atac-r1 "${atac_r1}"
    --atac-barcode "${atac_barcode}"
    --atac-r2 "${atac_r2}"
    --out-dir "${sample_out}"
    --threads "${THREADS}"
    --chromap-threads "${CHROMAP_THREADS}"
    --chromap-low-mem-ram "${CHROMAP_LOW_MEM_RAM}"
    --chromap-start-mode "${CHROMAP_START_MODE}"
    --skip-build
    --stop-after-local-mex
  )
  [[ "${CHROMAP_LOW_MEM}" == "1" ]] && args+=(--chromap-low-mem)
  [[ "${CHROMAP_MACS3_FRAG_LOW_MEM}" == "1" ]] && args+=(--chromap-macs3-frag-low-mem)
  [[ "${FORCE}" == "1" ]] && args+=(--force)

  "${args[@]}" 2>&1 | tee "${sample_log}"
  log "Local MEX boundary complete for ${sample}; queueing remote post-MEX work"

  post_log="${OUTPUT_ROOT}/logs/${sample_slug}.post_mex_remote.log"
  (
    set -euo pipefail
    remote_args=(
      "${REMOTE_POST_MEX}"
      --sample-dir "${sample_out}/star_sample"
      --remote-host "${REMOTE_HOST}"
      --remote-root "${REMOTE_ROOT}"
      --output-name downstream_genefull_velocyto_cellbender
      --run-cellbender
      --adaptive-filter
      --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}"
      --local-log "${sample_out}/logs/remote_post_mex.log"
    )
    [[ "${CELLBENDER_GPU}" == "1" ]] && remote_args+=(--cellbender-gpu)
    [[ "${NO_SYNC_IMAGES}" == "1" ]] && remote_args+=(--no-sync-images)
    [[ "${KEEP_REMOTE}" == "1" ]] && remote_args+=(--keep-remote)
    flock "${REMOTE_POST_MEX_LOCK}" "${remote_args[@]}"
    if [[ "${GLOBUS_UPLOAD_LARGE_FILES}" == "1" ]]; then
      "${GLOBUS_UPLOADER}" \
        --run-root "${OUTPUT_ROOT}" \
        --samples "${sample_slug}" \
        --source-endpoint "${GLOBUS_SOURCE_ENDPOINT}" \
        --dest-endpoint "${GLOBUS_DEST_ENDPOINT}" \
        --dest-root "${GLOBUS_DEST_ROOT}"
    fi
  ) > "${post_log}" 2>&1 < /dev/null &
  post_mex_pids+=("$!")
  post_mex_samples+=("${sample_slug}")
  log "Remote post-MEX PID for ${sample}: ${post_mex_pids[-1]} (log ${post_log})"
done < <(tail -n +2 "${MANIFEST}")

post_failures=0
for i in "${!post_mex_pids[@]}"; do
  pid="${post_mex_pids[$i]}"
  slug="${post_mex_samples[$i]}"
  if wait "${pid}"; then
    log "Remote post-MEX complete for ${slug}"
  else
    log "ERROR: remote post-MEX failed for ${slug}; see ${OUTPUT_ROOT}/logs/${slug}.post_mex_remote.log"
    post_failures=$((post_failures + 1))
  fi
done

if [[ "${post_failures}" -gt 0 ]]; then
  echo "ERROR: ${post_failures} remote post-MEX jobs failed" >&2
  exit 1
fi

log "PASS: JAX_Multiome01 production workflow complete"
echo "Output root: ${OUTPUT_ROOT}"
