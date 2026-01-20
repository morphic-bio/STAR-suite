#!/usr/bin/env bash
set -euo pipefail

ROOT="${CR_MULTI_ROOT:-/storage/A375}"
FASTQ_ROOT="${CR_MULTI_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${CR_MULTI_GEX_DIR:-${FASTQ_ROOT}/gex}"
CRISPR_DIR="${CR_MULTI_CRISPR_DIR:-${FASTQ_ROOT}/crispr}"
DOWN_SCRIPT="${CR_MULTI_DOWNSAMPLE_SCRIPT:-/mnt/pikachu/process_features/scripts/downsample_fastq_directory.sh}"
TIER_SCRIPT="${CR_MULTI_TIER_SCRIPT:-/mnt/pikachu/STAR-suite/scripts/a375_make_downsample_tiers.sh}"
TIERS=(10000 50000 100000 500000 1000000)

FEATURE_REF="${CR_MULTI_FEATURE_REF:-${ROOT}/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
WHITELIST_GZ="${CR_MULTI_WHITELIST_GZ:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-5pgex-jan-2023.txt.gz}"
WHITELIST="${CR_MULTI_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

STAR_BIN="${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
OUTPREFIX="${CR_MULTI_OUTPREFIX:-${ROOT}/star_multi_smoke_cpp/}"

if [[ ! -x "${DOWN_SCRIPT}" ]]; then
  echo "Missing downsample script: ${DOWN_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${WHITELIST}" ]]; then
  if [[ ! -f "${WHITELIST_GZ}" ]]; then
    echo "Missing whitelist: ${WHITELIST_GZ}" >&2
    exit 1
  fi
  zcat "${WHITELIST_GZ}" > "${WHITELIST}"
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "Missing genome index: ${GENOME_DIR}" >&2
  exit 1
fi

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
  exit 1
fi

if [[ -x "${TIER_SCRIPT}" ]]; then
  "${TIER_SCRIPT}"
fi

mkdir -p "${OUTPREFIX}"

validate_mex() {
  local mex_dir="$1"
  local mex_name="$2"
  RAW_MEX="${mex_dir}" MEX_NAME="${mex_name}" python3 - <<'PY'
import os
import sys
import gzip
from collections import defaultdict

def open_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")

mex_dir = os.environ["RAW_MEX"]
mex_name = os.environ["MEX_NAME"]

features_path = os.path.join(mex_dir, "features.tsv.gz")
matrix_path = os.path.join(mex_dir, "matrix.mtx.gz")
barcodes_path = os.path.join(mex_dir, "barcodes.tsv.gz")

if not os.path.exists(features_path):
    features_path = os.path.join(mex_dir, "features.tsv")
if not os.path.exists(matrix_path):
    matrix_path = os.path.join(mex_dir, "matrix.mtx")
if not os.path.exists(barcodes_path):
    barcodes_path = os.path.join(mex_dir, "barcodes.tsv")

for path in (features_path, matrix_path, barcodes_path):
    if not os.path.exists(path):
        raise SystemExit(f"Missing {mex_name} output file: {path}")

types = []
with open_file(features_path) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        ftype = parts[2] if len(parts) > 2 else "Unknown"
        types.append(ftype)

barcode_list = []
with open_file(barcodes_path) as f:
    for line in f:
        bc = line.rstrip("\n")
        barcode_list.append(bc)
        if not bc.endswith("-1"):
            raise SystemExit(f"{mex_name} barcode does not end with -1: {bc}")

if barcode_list != sorted(barcode_list):
    raise SystemExit(f"{mex_name} barcodes are not lexicographically sorted")

sums = defaultdict(int)
matrix_rows = 0
matrix_cols = 0
matrix_nnz = 0
with open_file(matrix_path) as f:
    header_seen = False
    for line in f:
        if line.startswith("%"):
            continue
        parts = line.split()
        if not header_seen:
            header_seen = True
            if len(parts) >= 3:
                matrix_rows = int(parts[0])
                matrix_cols = int(parts[1])
                matrix_nnz = int(parts[2])
            continue
        row = int(parts[0]) - 1
        val = int(float(parts[2]))
        ftype = types[row] if row < len(types) else "Unknown"
        sums[ftype] += val

# Debug: Print feature type distribution
type_counts = {}
for t in types:
    type_counts[t] = type_counts.get(t, 0) + 1

print(f"{mex_name} DEBUG: Matrix dimensions: rows={matrix_rows} cols={matrix_cols} nnz={matrix_nnz}", file=sys.stderr)
print(f"{mex_name} DEBUG: Feature type counts in features.tsv: {type_counts}", file=sys.stderr)
print(f"{mex_name} DEBUG: Sums keys found in matrix: {sorted(sums.keys())}", file=sys.stderr)
print(f"{mex_name} DEBUG: Sums values: {dict(sums)}", file=sys.stderr)

gene_sum = sums.get("Gene Expression", 0)
crispr_sum = sums.get("CRISPR Guide Capture", 0)

# Check if Gene Expression features exist in features.tsv but have no counts
gene_features_in_tsv = type_counts.get("Gene Expression", 0)
crispr_features_in_tsv = type_counts.get("CRISPR Guide Capture", 0)

# Check if we should allow missing feature types (for testing/debugging)
allow_missing_gex = os.environ.get("ALLOW_MISSING_GEX", "0") == "1"

if gene_sum == 0 or crispr_sum == 0:
    # Provide detailed diagnostics
    print(f"{mex_name} VALIDATION FAILED:", file=sys.stderr)
    print(f"  Gene Expression: {gene_sum} counts (from {gene_features_in_tsv} features in features.tsv)", file=sys.stderr)
    print(f"  CRISPR Guide Capture: {crispr_sum} counts (from {crispr_features_in_tsv} features in features.tsv)", file=sys.stderr)
    print(f"  Matrix has {matrix_nnz} non-zero entries total", file=sys.stderr)
    
    if gene_features_in_tsv > 0 and gene_sum == 0:
        print(f"  WARNING: {gene_features_in_tsv} Gene Expression features exist in features.tsv but have zero counts in matrix", file=sys.stderr)
        print(f"  This indicates GEX MEX data was empty or not merged correctly into the combined MEX", file=sys.stderr)
        print(f"  Check Solo.out/Gene/raw/matrix.mtx to verify GEX counts exist", file=sys.stderr)
        if allow_missing_gex:
            print(f"  ALLOW_MISSING_GEX=1: Allowing this case for testing/debugging", file=sys.stderr)
        else:
            print(f"  Set ALLOW_MISSING_GEX=1 to allow this case during testing", file=sys.stderr)
    
    if crispr_features_in_tsv > 0 and crispr_sum == 0:
        print(f"  WARNING: {crispr_features_in_tsv} CRISPR Guide Capture features exist in features.tsv but have zero counts in matrix", file=sys.stderr)
        print(f"  This indicates CRISPR MEX data was empty or not merged correctly", file=sys.stderr)
    
    # Allow missing GEX if explicitly enabled (for testing when GEX data is intentionally missing)
    if allow_missing_gex and gene_sum == 0 and crispr_sum > 0:
        print(f"{mex_name} WARNING (allowed): Gene Expression=0, CRISPR Guide Capture={crispr_sum} (ALLOW_MISSING_GEX=1)", file=sys.stderr)
        print(f"{mex_name} PARTIAL OK: barcodes={len(barcode_list)} Gene={gene_sum} CRISPR={crispr_sum} (GEX missing allowed)")
        sys.exit(0)
    else:
        raise SystemExit(
            f"{mex_name} comparable counts missing: Gene Expression={gene_sum}, "
            f"CRISPR Guide Capture={crispr_sum}"
        )

# Success case: both feature types have counts
print(f"{mex_name} OK: barcodes={len(barcode_list)} Gene={gene_sum} CRISPR={crispr_sum}")
PY
}

run_one_tier() {
  local tier="$1"
  local gex_tier_dir="${GEX_DIR}/downsampled_${tier}_v2"
  local crispr_tier_dir="${CRISPR_DIR}/downsampled_${tier}_v2"
  local tier_out="${OUTPREFIX}/tier_${tier}/"

  if [[ ! -d "${gex_tier_dir}" ]] || ! ls "${gex_tier_dir}"/*.fastq* >/dev/null 2>&1; then
    echo "Missing GEX tier: ${gex_tier_dir}" >&2
    return 2
  fi
  if [[ ! -d "${crispr_tier_dir}" ]] || ! ls "${crispr_tier_dir}"/*.fastq* >/dev/null 2>&1; then
    echo "Missing CRISPR tier: ${crispr_tier_dir}" >&2
    return 2
  fi

  local multi_config="${tier_out}/multi_config.csv"
  mkdir -p "${tier_out}"
  cat > "${multi_config}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${gex_tier_dir},A375,Gene Expression,Gene Expression
${crispr_tier_dir},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

  local r1_files
  local r2_files
  r1_files=$(ls "${gex_tier_dir}/"*R1*.fastq.gz | paste -sd, -)
  r2_files=$(ls "${gex_tier_dir}/"*R2*.fastq.gz | paste -sd, -)

  set +e
  "${STAR_BIN}" \
    --runThreadN 4 \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${r2_files}" "${r1_files}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${tier_out}" \
    --outSAMtype BAM Unsorted \
    --clipAdapterType CellRanger4 \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIdedup 1MM_CR \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloMultiMappers Rescue \
    --soloStrand Unstranded \
    --soloAddTagsToUnsorted no \
    --soloFeatures GeneFull \
    --soloCellFilter None \
    --crMultiConfig "${multi_config}" \
    --crFeatureRef "${FEATURE_REF}" \
    --crWhitelist "${WHITELIST}"
  local exit_code=$?
  set -e

  if [[ ${exit_code} -ne 0 ]]; then
    echo "STAR failed for tier ${tier} (exit ${exit_code})"
    return 1
  fi

  local raw_mex="${tier_out}outs/raw_feature_bc_matrix"
  local filtered_mex="${tier_out}outs/filtered_feature_bc_matrix"

  echo "Validating tier ${tier} outputs:"
  echo "  Raw: ${raw_mex}"
  echo "  Filtered: ${filtered_mex}"

  validate_mex "${raw_mex}" "raw"
  validate_mex "${filtered_mex}" "filtered"

  return 0
}

for tier in "${TIERS[@]}"; do
  echo "=== Running tier ${tier} ==="
  if run_one_tier "${tier}"; then
    echo "Tier ${tier} completed without crash"
  else
    echo "Tier ${tier} reproduced the crash"
    echo "First failing tier: ${tier}"
    exit 1
  fi
done

echo "No crash reproduced in tiers: ${TIERS[*]}"
