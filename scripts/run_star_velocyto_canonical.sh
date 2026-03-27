#!/usr/bin/env bash
# Canonical UCSF perturb STAR invocation for Gene + GeneFull + Velocyto parity work.
# Pins solo / CR-compat flags so comparisons are not confounded by parameter drift.
#
# Usage:
#   run_star_velocyto_canonical.sh --profile 100k|2m --threads N --out-prefix /path/to/run/ [--prepare-mex]
#
# Optional env:
#   STAR_BIN, STAR_VELOCYTO_DETERMINISTIC_REPLAY, STAR_VELOCYTO_INTEGRATED_HASH (passed through to STAR)
#   INTEGRATED_HASH=1 with DETERMINISTIC_REPLAY=1 selects Stage 2 (default: disk spill shards; see below).
#   STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS (default 128, max 4096) — shard records to bound RAM during merge.
#   STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY=1 — hold all per-CB records in RAM (debug / A–B vs spill only).
#   For profile=2m (full-sample harness; canonical UCSF surface is corrected EBs2_2):
#     UCSF_2M_PFCONFIG — required; pfMultiConfig CSV (no default).
#     UCSF_2M_GEX_DIR, UCSF_2M_GUIDE_DIR — default under UCSF_2M_ROOT (see below).
#     UCSF_2M_FEATURE_REF — CRISPR/feature reference CSV; defaults to UCSF_100K_FEATURE_REF, then repo path.
#     UCSF_2M_CB_WHITELIST — Solo CB whitelist; defaults to UCSF_100K_CB_WHITELIST, then repo path.
#     UCSF_2M_GENOME_DIR — STAR genomeDir; defaults to UCSF_100K_GENOME_DIR, then bulk_index.
#   Export UCSF_2M_FEATURE_REF / UCSF_2M_CB_WHITELIST explicitly when site paths differ from defaults.
#
# Memory: STAR_VELOCYTO_DETERMINISTIC_REPLAY=1 materializes every Velocyto stream record before merge.
# On multi-hundred-million-read surfaces (full 2M), budget peak RSS and read Log.out "RAM after Velocyto sorted-replay" lines.
#
# Output hygiene: --out-prefix must be a fresh directory (no Solo.out / Log.out / prior packaged outs).
# Override only when intentional: UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${EXTERNAL_ENV}"
fi

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
PREPARE_VELOCYTO_MEX="${REPO_ROOT}/scripts/prepare_velocyto_mex.py"

PROFILE=""
THREADS=""
OUT_PREFIX=""
PREPARE_MEX=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      PROFILE="${2:?}"
      shift 2
      ;;
    --threads)
      THREADS="${2:?}"
      shift 2
      ;;
    --out-prefix)
      OUT_PREFIX="${2:?}"
      shift 2
      ;;
    --prepare-mex)
      PREPARE_MEX=1
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -n "${PROFILE}" ]] || die "--profile 100k|2m is required"
[[ -n "${THREADS}" ]] || die "--threads N is required"
[[ -n "${OUT_PREFIX}" ]] || die "--out-prefix is required"

join_by_comma() {
  local IFS=,
  echo "$*"
}

[[ -x "${STAR_BIN}" ]] || die "STAR not executable: ${STAR_BIN}"

if [[ "${PROFILE}" == "100k" ]]; then
  FIXTURE_ROOT="${UCSF_100K_PFCONFIG_ROOT:-/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig}"
  GEX_DIR="${UCSF_100K_GEX_DIR:-${FIXTURE_ROOT}/GEX/iPSC2_1_AALG2}"
  GUIDE_DIR="${UCSF_100K_GUIDE_DIR:-${FIXTURE_ROOT}/guides/iPSC2_1_AALG2}"
  PF_MULTI_CONFIG="${UCSF_100K_PFCONFIG:-${FIXTURE_ROOT}/multi_config_100k.csv}"
  FEATURE_REF="${UCSF_100K_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
  WHITELIST="${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
  GENOME_DIR="${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
elif [[ "${PROFILE}" == "2m" ]]; then
  # Default tree: corrected full-sample EBs2_2 (repo benchmark surface). Override UCSF_2M_ROOT / *_GEX_DIR / *_GUIDE_DIR for other layouts.
  UCSF_2M_ROOT="${UCSF_2M_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2}"
  GEX_DIR="${UCSF_2M_GEX_DIR:-${UCSF_2M_ROOT}/GEX}"
  GUIDE_DIR="${UCSF_2M_GUIDE_DIR:-${UCSF_2M_ROOT}/guides}"
  PF_MULTI_CONFIG="${UCSF_2M_PFCONFIG:?Set UCSF_2M_PFCONFIG to the full-sample pfMultiConfig CSV}"
  FEATURE_REF="${UCSF_2M_FEATURE_REF:-${UCSF_100K_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}}"
  WHITELIST="${UCSF_2M_CB_WHITELIST:-${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}}"
  GENOME_DIR="${UCSF_2M_GENOME_DIR:-${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}}"
else
  die "Unknown --profile ${PROFILE} (use 100k or 2m)"
fi

[[ -d "${GEX_DIR}" ]] || die "GEX dir missing: ${GEX_DIR}"
[[ -d "${GUIDE_DIR}" ]] || die "Guide dir missing: ${GUIDE_DIR}"
[[ -f "${PF_MULTI_CONFIG}" ]] || die "pfMultiConfig missing: ${PF_MULTI_CONFIG}"
[[ -f "${FEATURE_REF}" ]] || die "Feature ref missing: ${FEATURE_REF}"
[[ -f "${WHITELIST}" ]] || die "Whitelist missing: ${WHITELIST}"
[[ -d "${GENOME_DIR}" ]] || die "Genome dir missing: ${GENOME_DIR}"

# -L: follow symlinks (corrected EBs2_2 layout uses symlinked FASTQs in GEX/guides).
mapfile -t GEX_R2_FILES < <(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GEX_R1_FILES < <(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
mapfile -t GUIDE_R2_FILES < <(find -L "${GUIDE_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GUIDE_R1_FILES < <(find -L "${GUIDE_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)

[[ "${#GEX_R2_FILES[@]}" -gt 0 ]] || die "No GEX R2 FASTQs under ${GEX_DIR}"
[[ "${#GUIDE_R2_FILES[@]}" -gt 0 ]] || die "No guide R2 FASTQs under ${GUIDE_DIR}"
[[ "${#GEX_R2_FILES[@]}" -eq "${#GEX_R1_FILES[@]}" ]] || die "GEX R1/R2 mismatch"
[[ "${#GUIDE_R2_FILES[@]}" -eq "${#GUIDE_R1_FILES[@]}" ]] || die "Guide R1/R2 mismatch"

ALL_R2=("${GEX_R2_FILES[@]}" "${GUIDE_R2_FILES[@]}")
ALL_R1=("${GEX_R1_FILES[@]}" "${GUIDE_R1_FILES[@]}")
R2_LIST="$(join_by_comma "${ALL_R2[@]}")"
R1_LIST="$(join_by_comma "${ALL_R1[@]}")"

OUT_PREFIX="${OUT_PREFIX%/}"

if [[ "${UCSF_VELOCYTO_REUSE_STAR_OUTDIR:-0}" != "1" ]]; then
  if [[ -e "${OUT_PREFIX}/Solo.out" || -e "${OUT_PREFIX}/Log.out" || -e "${OUT_PREFIX}/Aligned.out.bam" || -e "${OUT_PREFIX}/Aligned.out.sam" ]]; then
    die "Refusing non-fresh --out-prefix ${OUT_PREFIX} (Solo.out, Log.out, Aligned.out.bam, or Aligned.out.sam present). Use a new directory or UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
  if [[ -d "${OUT_PREFIX}/outs/raw_velocyto_feature_bc_matrix" || -d "${OUT_PREFIX}/outs/filtered_velocyto_feature_bc_matrix" ]]; then
    die "Refusing non-fresh --out-prefix ${OUT_PREFIX} (outs/*velocyto* present). Remove or pick a new path, or UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
fi

mkdir -p "${OUT_PREFIX}"

echo "=== run_star_velocyto_canonical ==="
echo "profile=${PROFILE} threads=${THREADS} out=${OUT_PREFIX}"
if [[ -n "${STAR_VELOCYTO_DETERMINISTIC_REPLAY:-}" ]]; then
  echo "STAR_VELOCYTO_DETERMINISTIC_REPLAY=${STAR_VELOCYTO_DETERMINISTIC_REPLAY}"
fi
if [[ -n "${STAR_VELOCYTO_INTEGRATED_HASH:-}" ]]; then
  echo "STAR_VELOCYTO_INTEGRATED_HASH=${STAR_VELOCYTO_INTEGRATED_HASH}"
fi
if [[ -n "${STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS:-}" ]]; then
  echo "STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS=${STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS}"
fi
if [[ -n "${STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY:-}" ]]; then
  echo "STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY=${STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY}"
fi

# shellcheck disable=SC2086
"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_PREFIX}/" \
  --outSAMtype BAM Unsorted \
  --clipAdapterType CellRanger4 \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Unstranded \
  --soloFeatures Gene GeneFull Velocyto \
  --soloCrGexFeature genefull \
  --pfMultiConfig "${PF_MULTI_CONFIG}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crWhitelist "${WHITELIST}" \
  --crMinUmi 10

if [[ "${PREPARE_MEX}" -eq 1 ]]; then
  [[ -f "${PREPARE_VELOCYTO_MEX}" ]] || die "Missing ${PREPARE_VELOCYTO_MEX}"
  python3 "${PREPARE_VELOCYTO_MEX}" --run-dir "${OUT_PREFIX}"
fi

echo "OK: ${OUT_PREFIX}"
