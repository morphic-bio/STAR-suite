#!/usr/bin/env bash
# Raw counts retain cells AND an ambient plateau; this is not a 100K-read test.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec python3 "${ROOT}/tests/host_api/run_cellbender_cuda_smoke.py" \
  --raw-mex "${CELLBENDER_SMOKE_RAW_MEX:-/storage/A375/paper_bench_20260326_134444/outs/raw_feature_bc_matrix}" \
  --recipes-root "${MORPHIC_RECIPES_ROOT:-/mnt/pikachu/morphic-recipes}" \
  --outdir "${CELLBENDER_SMOKE_OUT:-/tmp/star_cellbender_cuda_$(date -u +%Y%m%d_%H%M%S)_$$}" \
  --expected-cells "${CELLBENDER_SMOKE_EXPECTED_CELLS:-1188}" \
  --image "${CELLBENDER_IMAGE:-biodepot/cellbender:0.3.2}" \
  --cellbender-gpu
