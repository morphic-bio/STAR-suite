#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
processing_repo="${VISIUM_HD_PROCESSING_REPO:-/mnt/pikachu/visium-hd-processing}"
resolver="${MOLECULE_FIRST_RESOLVER:-${repo_root}/flex/tools/molecule_first_resolver/molecule_first_resolver}"

if [[ ! -f "${processing_repo}/src/star_spatial/hd_probabilistic_umi.py" ]]; then
  echo "ERROR: frozen processing reference not found: ${processing_repo}" >&2
  exit 2
fi

make -C "${repo_root}/flex/tools/molecule_first_resolver" molecule_first_resolver
python3 "${repo_root}/tests/molecule_first/compare_reference.py" \
  --resolver "${resolver}" \
  --processing-repo "${processing_repo}" \
  --fixture "${repo_root}/tests/molecule_first/candidate_reads.tsv"
