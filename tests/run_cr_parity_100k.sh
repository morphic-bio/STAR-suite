#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1091
  source "${EXTERNAL_ENV}"
fi

# Parity run configuration (override via env vars as needed)
CR_PARITY_TIER="${CR_PARITY_TIER:-100000}"
CR_PARITY_THREADS="${CR_PARITY_THREADS:-4}"
CR_PARITY_GENOME_DIR="${CR_PARITY_GENOME_DIR:-/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A/star_2.7.11b}"
CR_PARITY_OUT_ROOT="${CR_PARITY_OUT_ROOT:-/tmp/cr_parity_100k}"
CR_PARITY_GEX_DIR="${CR_PARITY_GEX_DIR:-/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/gex/downsampled_${CR_PARITY_TIER}_v2}"

CR_COMPAT_OUT="${CR_PARITY_OUT_ROOT}/cr_compat"
STAR_GEX_OUT="${CR_PARITY_OUT_ROOT}/star_genefull"

CR_COMPAT_MEX="${CR_COMPAT_OUT}/tier_${CR_PARITY_TIER}/outs/raw_feature_bc_matrix"
STAR_MEX="${STAR_GEX_OUT}/Solo.out/GeneFull/raw"
export CR_COMPAT_MEX

echo "== CR parity run (tier ${CR_PARITY_TIER}) =="
echo "  CR_GENOME_DIR=${CR_PARITY_GENOME_DIR}"
echo "  CR_PARITY_GEX_DIR=${CR_PARITY_GEX_DIR}"
echo "  CR_COMPAT_OUT=${CR_COMPAT_OUT}"
echo "  STAR_GEX_OUT=${STAR_GEX_OUT}"

CR_MULTI_TIERS="${CR_PARITY_TIER}" \
CR_MULTI_OUTPREFIX="${CR_COMPAT_OUT}/" \
CR_GENOME_DIR="${CR_PARITY_GENOME_DIR}" \
CR_MULTI_SOLO_MULTIMAPPERS=Unique \
CR_MULTI_SOLO_CELL_FILTER=EmptyDrops_CR \
CR_MULTI_GEX_FEATURE=genefull \
bash "${REPO_ROOT}/tests/cr_multi_smoke/run_cr_multi_smoke.sh"

A375_SKIP_DOWNSAMPLE=1 \
A375_CRLIKE_OUTPREFIX="${STAR_GEX_OUT}/" \
A375_THREADS="${CR_PARITY_THREADS}" \
A375_SOLO_MULTIMAPPERS=Unique \
A375_REQUIRE_CBUB_TOGETHER=no \
A375_WRITE_BAM=1 \
CR_MULTI_GEX_DIR="${CR_PARITY_GEX_DIR}" \
CR_GENOME_DIR="${CR_PARITY_GENOME_DIR}" \
bash "${REPO_ROOT}/tests/run_a375_cr_like_gex.sh"

python3 "${REPO_ROOT}/tests/compare_a375_star_mex.py" "${CR_COMPAT_MEX}" "${STAR_MEX}"

echo "== CR-format checks =="
if [[ ! -f "${CR_COMPAT_MEX}/matrix.mtx.gz" ]] || [[ ! -f "${CR_COMPAT_MEX}/features.tsv.gz" ]] || [[ ! -f "${CR_COMPAT_MEX}/barcodes.tsv.gz" ]]; then
  echo "Missing CR-compat MEX files in ${CR_COMPAT_MEX}" >&2
  exit 1
fi

python3 - <<'PY'
import gzip
import os
import sys

mex_dir = os.environ["CR_COMPAT_MEX"]
barcodes = os.path.join(mex_dir, "barcodes.tsv.gz")
features = os.path.join(mex_dir, "features.tsv.gz")

prev = None
with gzip.open(barcodes, "rt", encoding="utf-8") as fh:
    for line in fh:
        bc = line.strip()
        if not bc.endswith("-1"):
            raise SystemExit(f"Barcode missing -1 suffix: {bc}")
        if prev is not None and bc < prev:
            raise SystemExit("Barcodes are not lexicographically sorted")
        prev = bc

with gzip.open(features, "rt", encoding="utf-8") as fh:
    for i, line in enumerate(fh, start=1):
        if len(line.rstrip("\n").split("\t")) < 3:
            raise SystemExit(f"features.tsv.gz line {i} has < 3 columns")

print("CR-format checks OK")
PY

echo "CR parity run completed."
