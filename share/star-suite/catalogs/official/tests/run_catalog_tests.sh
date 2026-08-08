#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

python3 tests/test_partition_manifest.py
python3 tests/validate_catalog.py
python3 tests/scan_publication.py
bash -n scripts/*.sh scripts/docker/*.sh

for script in \
  scripts/run_partitioned_bulk_pe_alignment.sh \
  scripts/run_star_multiome_lane_smoke.sh \
  scripts/run_multiome_minimal.sh \
  scripts/run_remote_multiome_post_mex_rsync.sh \
  scripts/run_scrna_downstream_gene_full_velocyto.sh \
  scripts/run_catatac_trimodal_downsample_smoke.sh \
  scripts/run_hiv_dogma_four_arm_downsample_smoke.sh \
  scripts/run_slam_100k_se_pe_smoke.sh \
  scripts/run_slam_prod_set.sh \
  scripts/docker/build_slam_deseq2_image.sh
do
  "${script}" --help >/dev/null
done

python3 scripts/validate_partition_manifest.py \
  tests/fixtures/partition_manifest.json >/dev/null
echo "PASS: STAR Suite recipe catalog tests"

