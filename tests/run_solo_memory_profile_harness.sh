#!/usr/bin/env bash
# Run a Solo post-map path with STAR_SOLO_MEMORY_PROFILE checkpoints in Log.out.
# Use this to see which counting structure dominates RSS on a failing run.
#
# Examples:
#   # JAX OCM 2M smoke (soloInlineHashMode=no → CountingSink + PackedReadInfo + Velocyto):
#   tests/run_solo_memory_profile_harness.sh \
#     --jax-ocm-smoke --downsample-read-pairs 2000000 --run-star
#
#   # Ladder (100k → 500k → 2M) to bisect OOM:
#   tests/run_solo_memory_profile_harness.sh --jax-ocm-smoke --ladder 100000,500000,2000000 --run-star
#
#   # UCSF GeneFull+Velocyto 100k (CountingSink path, no inline hash):
#   tests/run_solo_memory_profile_harness.sh --ucsf-velocyto-100k --threads 8 --run-star
#
#   # Parse an existing log:
#   tests/run_solo_memory_profile_harness.sh --parse-log /path/to/logs/star.log
#
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
MODE="jax_ocm_smoke"
DOWNSAMPLE_ARGS=()
RUN_STAR="0"
PARSE_LOG=""
OUT_ROOT=""
LADDER=""
UCSF_THREADS="${UCSF_THREADS:-8}"
UCSF_OUT_PREFIX=""

usage() {
  cat <<'EOF'
Usage:
  run_solo_memory_profile_harness.sh [options]

Options:
  --jax-ocm-smoke           Use scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh (default)
  --ucsf-velocyto-100k      GeneFull+Velocyto via run_star_velocyto_canonical.sh --profile 100k
  --downsample-read-pairs N Passed to JAX smoke
  --ladder N1,N2,...        Run JAX smoke once per downsample size (implies --run-star)
  --threads N               Thread count for UCSF 100k mode (default: 8)
  --out-prefix PATH         UCSF STAR out prefix (fresh directory)
  --run-star                Execute STAR (otherwise only stages RUN_STAR.sh)
  --out-root PATH           JAX smoke output root
  --parse-log PATH          Print Solo memory checkpoints + RSS deltas from an existing log
  -h, --help

Environment (set automatically unless already exported):
  STAR_SOLO_MEMORY_PROFILE=1
  STAR_SOLO_PHASE_DEBUG=1   (unless already set)

Look for lines:
  Solo memory: <phase> | <counters> | VmPeak: ...; VmRSS: ...; VmHWM: ...
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --jax-ocm-smoke) MODE="jax_ocm_smoke"; shift ;;
    --ucsf-velocyto-100k) MODE="ucsf_velocyto_100k"; shift ;;
    --downsample-read-pairs) DOWNSAMPLE_ARGS=(--downsample-read-pairs "$2"); shift 2 ;;
    --ladder) LADDER="$2"; shift 2 ;;
    --threads) UCSF_THREADS="$2"; shift 2 ;;
    --out-prefix) UCSF_OUT_PREFIX="$2"; shift 2 ;;
    --run-star) RUN_STAR="1"; shift ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --parse-log) PARSE_LOG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

parse_log() {
  local log="$1"
  [[ -f "${log}" ]] || { echo "Missing log: ${log}" >&2; exit 1; }
  echo "=== Solo memory checkpoints from ${log} ==="
  rg '^Solo memory:' "${log}" || true
  echo
  echo "=== VmRSS progression (label → kB → delta kB) ==="
  rg '^Solo memory:' "${log}" | awk '
    {
      label = $3
      rss = 0
      if (match($0, /VmRSS:[[:space:]]+([0-9]+) kB/, m)) rss = m[1] + 0
      delta = (NR == 1) ? 0 : rss - prev
      printf "%3d  %+8d kB  %s\n", NR, delta, label
      prev = rss
    }
  ' || true
  echo
  echo "=== Largest single-step VmRSS jumps ==="
  rg '^Solo memory:' "${log}" | awk '
    {
      label = $3
      if (match($0, /VmRSS:[[:space:]]+([0-9]+) kB/, m)) rss = m[1] + 0
      if (NR > 1) print (rss - prev), label
      prev = rss
    }
  ' | sort -nr | head -8 || true
}

if [[ -n "${PARSE_LOG}" ]]; then
  parse_log "${PARSE_LOG}"
  exit 0
fi

export STAR_SOLO_MEMORY_PROFILE=1
export STAR_SOLO_PHASE_DEBUG="${STAR_SOLO_PHASE_DEBUG:-1}"
export STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"

log_msg() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"; }

run_jax_smoke() {
  local ds="$1"
  local root="${2:-}"
  local cmd=(
    "${REPO_ROOT}/scripts/run_jax_scrnaseq02_ocm_oracle_smoke.sh"
    --downsample-read-pairs "${ds}"
  )
  if [[ -n "${root}" ]]; then
    cmd+=(--out-root "${root}")
  fi
  if [[ "${RUN_STAR}" == "1" ]]; then
    cmd+=(--run-star)
  fi
  log_msg "STAR_SOLO_MEMORY_PROFILE=1 downsample=${ds}"
  log_msg "Running: ${cmd[*]}"
  "${cmd[@]}"
  local star_log log_main
  if [[ -z "${root}" ]]; then
    root="$(ls -td /mnt/pikachu/JAX_scRNAseq02_processed/ocm_oracle_smoke_* 2>/dev/null | head -1)"
  fi
  star_log="${root}/logs/star.log"
  log_main="$(find "${root}/samples" -path '*/run/Log.out' 2>/dev/null | head -1)"
  echo
  if [[ -n "${log_main}" && -f "${log_main}" ]]; then
    parse_log "${log_main}"
  elif [[ -f "${star_log}" ]]; then
    parse_log "${star_log}"
  else
    echo "WARN: could not find ${root}/samples/*/run/Log.out or ${star_log}" >&2
  fi
}

case "${MODE}" in
  jax_ocm_smoke)
    if [[ -n "${LADDER}" ]]; then
      RUN_STAR="1"
      IFS=',' read -r -a ladder_sizes <<< "${LADDER}"
      for ds in "${ladder_sizes[@]}"; do
        ds="${ds//[[:space:]]/}"
        [[ -n "${ds}" ]] || continue
        ladder_root="/mnt/pikachu/JAX_scRNAseq02_processed/ocm_memprof_${ds}_$(date -u +%Y%m%dT%H%M%SZ)"
        run_jax_smoke "${ds}" "${ladder_root}"
        echo
      done
    else
      ds=2000000
      if [[ ${#DOWNSAMPLE_ARGS[@]} -ge 2 ]]; then
        ds="${DOWNSAMPLE_ARGS[1]}"
      fi
      run_jax_smoke "${ds}" "${OUT_ROOT}"
    fi
    ;;
  ucsf_velocyto_100k)
    if [[ "${RUN_STAR}" != "1" ]]; then
      echo "ERROR: --ucsf-velocyto-100k requires --run-star" >&2
      exit 1
    fi
    if [[ -z "${UCSF_OUT_PREFIX}" ]]; then
      UCSF_OUT_PREFIX="/mnt/pikachu/JAX_scRNAseq02_processed/ucsf_velocyto_memprof_100k_$(date -u +%Y%m%dT%H%M%SZ)"
    fi
  # CountingSink path (not inline hash) for apples-to-apples with JAX OCM smoke.
    export STAR_EXTRA_ARGS="${STAR_EXTRA_ARGS:-} --soloInlineHashMode no --soloFeatures GeneFull Velocyto"
    log_msg "UCSF 100k GeneFull+Velocyto out=${UCSF_OUT_PREFIX}"
    "${REPO_ROOT}/scripts/run_star_velocyto_canonical.sh" \
      --profile 100k \
      --threads "${UCSF_THREADS}" \
      --out-prefix "${UCSF_OUT_PREFIX}"
    star_log="${UCSF_OUT_PREFIX}/Log.out"
    if [[ -f "${star_log}" ]]; then
      parse_log "${star_log}"
    else
      echo "WARN: missing ${star_log}" >&2
    fi
    ;;
  *)
    echo "Unknown mode: ${MODE}" >&2
    exit 1
    ;;
esac
