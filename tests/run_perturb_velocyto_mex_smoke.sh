#!/usr/bin/env bash
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

OUT_DIR="${PERTURB_VELOCYTO_MEX_SMOKE_OUTDIR:-${SCRIPT_DIR}/perturb_velocyto_mex_smoke_output_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${OUT_DIR}/run"
SUMMARY_FILE="${OUT_DIR}/SUMMARY.txt"

FIXTURE_ROOT="${UCSF_100K_PFCONFIG_ROOT:-/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig}"
GEX_DIR="${UCSF_100K_GEX_DIR:-${FIXTURE_ROOT}/GEX/iPSC2_1_AALG2}"
GUIDE_DIR="${UCSF_100K_GUIDE_DIR:-${FIXTURE_ROOT}/guides/iPSC2_1_AALG2}"
PF_MULTI_CONFIG="${UCSF_100K_PFCONFIG:-${FIXTURE_ROOT}/multi_config_100k.csv}"
FEATURE_REF="${UCSF_100K_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
WHITELIST="${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
GENOME_DIR="${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
THREADS="${PERTURB_VELOCYTO_MEX_THREADS:-1}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"
[[ -f "${PREPARE_VELOCYTO_MEX}" ]] || die "Missing helper: ${PREPARE_VELOCYTO_MEX}"
[[ -d "${GEX_DIR}" ]] || die "Missing GEX fixture dir: ${GEX_DIR}"
[[ -d "${GUIDE_DIR}" ]] || die "Missing guide fixture dir: ${GUIDE_DIR}"
[[ -f "${PF_MULTI_CONFIG}" ]] || die "Missing pfMultiConfig: ${PF_MULTI_CONFIG}"
[[ -f "${FEATURE_REF}" ]] || die "Missing feature ref: ${FEATURE_REF}"
[[ -f "${WHITELIST}" ]] || die "Missing whitelist: ${WHITELIST}"
[[ -d "${GENOME_DIR}" ]] || die "Missing genome dir: ${GENOME_DIR}"

mapfile -t GEX_R2_FILES < <(find "${GEX_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GEX_R1_FILES < <(find "${GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
mapfile -t GUIDE_R2_FILES < <(find "${GUIDE_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GUIDE_R1_FILES < <(find "${GUIDE_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)

[[ "${#GEX_R2_FILES[@]}" -gt 0 ]] || die "No GEX R2 FASTQs found under ${GEX_DIR}"
[[ "${#GUIDE_R2_FILES[@]}" -gt 0 ]] || die "No guide R2 FASTQs found under ${GUIDE_DIR}"
[[ "${#GEX_R2_FILES[@]}" -eq "${#GEX_R1_FILES[@]}" ]] || die "GEX R1/R2 lane count mismatch"
[[ "${#GUIDE_R2_FILES[@]}" -eq "${#GUIDE_R1_FILES[@]}" ]] || die "Guide R1/R2 lane count mismatch"

ALL_R2=("${GEX_R2_FILES[@]}" "${GUIDE_R2_FILES[@]}")
ALL_R1=("${GEX_R1_FILES[@]}" "${GUIDE_R1_FILES[@]}")

rm -rf "${OUT_DIR}"
mkdir -p "${RUN_DIR}"

R2_LIST="$(join_by_comma "${ALL_R2[@]}")"
R1_LIST="$(join_by_comma "${ALL_R1[@]}")"

echo "=== Perturb Velocyto MEX Smoke ==="
echo "Output: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo "GEX cDNA lanes: ${#GEX_R2_FILES[@]}"
echo "Guide cDNA lanes: ${#GUIDE_R2_FILES[@]}"

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${RUN_DIR}/" \
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

python3 "${PREPARE_VELOCYTO_MEX}" --run-dir "${RUN_DIR}"

RAW_VELO_DIR="${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix"
FILTERED_VELO_DIR="${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix"
MANIFEST="${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json"

for path in \
  "${RUN_DIR}/Solo.out/Velocyto/raw/spliced.mtx" \
  "${RUN_DIR}/Solo.out/Velocyto/raw/unspliced.mtx" \
  "${RUN_DIR}/Solo.out/Velocyto/raw/ambiguous.mtx" \
  "${RUN_DIR}/Solo.out/Gene/raw/barcodes.tsv" \
  "${RUN_DIR}/Solo.out/Gene/raw/features.tsv" \
  "${RUN_DIR}/Solo.out/Gene/filtered/barcodes.tsv" \
  "${RAW_VELO_DIR}/barcodes.tsv.gz" \
  "${RAW_VELO_DIR}/features.tsv.gz" \
  "${RAW_VELO_DIR}/matrix.mtx.gz" \
  "${RAW_VELO_DIR}/spliced.mtx.gz" \
  "${RAW_VELO_DIR}/unspliced.mtx.gz" \
  "${RAW_VELO_DIR}/ambiguous.mtx.gz" \
  "${FILTERED_VELO_DIR}/barcodes.tsv.gz" \
  "${FILTERED_VELO_DIR}/features.tsv.gz" \
  "${FILTERED_VELO_DIR}/matrix.mtx.gz" \
  "${FILTERED_VELO_DIR}/spliced.mtx.gz" \
  "${FILTERED_VELO_DIR}/unspliced.mtx.gz" \
  "${FILTERED_VELO_DIR}/ambiguous.mtx.gz" \
  "${MANIFEST}"; do
  [[ -f "${path}" ]] || die "Missing expected output: ${path}"
done

readarray -t VALIDATION < <(python3 - <<'PY' "${RAW_VELO_DIR}" "${FILTERED_VELO_DIR}" "${MANIFEST}"
import gzip
import json
import sys
from pathlib import Path

from scipy.io import mmread


def read_lines(path: Path):
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


def read_matrix(path: Path):
    with gzip.open(path, "rb") as handle:
        return mmread(handle).tocsr()


raw_dir = Path(sys.argv[1])
filtered_dir = Path(sys.argv[2])
manifest_path = Path(sys.argv[3])
manifest = json.loads(manifest_path.read_text())

raw_barcodes = read_lines(raw_dir / "barcodes.tsv.gz")
raw_features = read_lines(raw_dir / "features.tsv.gz")
raw_total = read_matrix(raw_dir / "matrix.mtx.gz")
raw_spliced = read_matrix(raw_dir / "spliced.mtx.gz")
raw_unspliced = read_matrix(raw_dir / "unspliced.mtx.gz")
raw_ambiguous = read_matrix(raw_dir / "ambiguous.mtx.gz")
raw_sum = raw_spliced + raw_unspliced + raw_ambiguous

if raw_total.shape != raw_sum.shape:
    raise SystemExit(f"raw total shape mismatch: {raw_total.shape} vs {raw_sum.shape}")
if raw_total.shape != (len(raw_features), len(raw_barcodes)):
    raise SystemExit(
        f"raw axis mismatch: matrix={raw_total.shape}, features={len(raw_features)}, barcodes={len(raw_barcodes)}"
    )
if (raw_total != raw_sum).nnz != 0:
    raise SystemExit("raw total matrix does not equal spliced+unspliced+ambiguous")

filtered_barcodes = read_lines(filtered_dir / "barcodes.tsv.gz")
filtered_features = read_lines(filtered_dir / "features.tsv.gz")
filtered_total = read_matrix(filtered_dir / "matrix.mtx.gz")
filtered_spliced = read_matrix(filtered_dir / "spliced.mtx.gz")
filtered_unspliced = read_matrix(filtered_dir / "unspliced.mtx.gz")
filtered_ambiguous = read_matrix(filtered_dir / "ambiguous.mtx.gz")
filtered_sum = filtered_spliced + filtered_unspliced + filtered_ambiguous

if filtered_total.shape != filtered_sum.shape:
    raise SystemExit(f"filtered total shape mismatch: {filtered_total.shape} vs {filtered_sum.shape}")
if filtered_total.shape != (len(filtered_features), len(filtered_barcodes)):
    raise SystemExit(
        "filtered axis mismatch: "
        f"matrix={filtered_total.shape}, features={len(filtered_features)}, barcodes={len(filtered_barcodes)}"
    )
if (filtered_total != filtered_sum).nnz != 0:
    raise SystemExit("filtered total matrix does not equal spliced+unspliced+ambiguous")

if raw_features != filtered_features:
    raise SystemExit("raw and filtered velocyto feature tables differ")
if not set(filtered_barcodes).issubset(set(raw_barcodes)):
    raise SystemExit("filtered barcodes are not a subset of raw barcodes")

if manifest["raw"]["features"] != len(raw_features) or manifest["raw"]["barcodes"] != len(raw_barcodes):
    raise SystemExit("manifest raw dimensions do not match emitted raw MEX")
if manifest["filtered"]["features"] != len(filtered_features) or manifest["filtered"]["barcodes"] != len(filtered_barcodes):
    raise SystemExit("manifest filtered dimensions do not match emitted filtered MEX")

print(f"RAW_FEATURES={len(raw_features)}")
print(f"RAW_BARCODES={len(raw_barcodes)}")
print(f"FILTERED_BARCODES={len(filtered_barcodes)}")
print(f"RAW_NNZ_TOTAL={raw_total.nnz}")
print(f"FILTERED_NNZ_TOTAL={filtered_total.nnz}")
PY
)

declare -A STATS=()
for line in "${VALIDATION[@]}"; do
  key="${line%%=*}"
  value="${line#*=}"
  STATS["${key}"]="${value}"
done

{
  echo "Perturb velocyto MEX smoke: PASS"
  echo "Output: ${OUT_DIR}"
  echo "Threads: ${THREADS}"
  echo "GEX cDNA lanes: ${#GEX_R2_FILES[@]}"
  echo "Guide cDNA lanes: ${#GUIDE_R2_FILES[@]}"
  echo "Raw velocyto features: ${STATS[RAW_FEATURES]}"
  echo "Raw velocyto barcodes: ${STATS[RAW_BARCODES]}"
  echo "Filtered velocyto barcodes: ${STATS[FILTERED_BARCODES]}"
  echo "Raw velocyto nnz total: ${STATS[RAW_NNZ_TOTAL]}"
  echo "Filtered velocyto nnz total: ${STATS[FILTERED_NNZ_TOTAL]}"
  echo "Raw velocyto MEX: ${RAW_VELO_DIR}"
  echo "Filtered velocyto MEX: ${FILTERED_VELO_DIR}"
  echo "Manifest: ${MANIFEST}"
} > "${SUMMARY_FILE}"

echo "PASS: Perturb velocyto MEX smoke"
echo "Summary: ${SUMMARY_FILE}"
