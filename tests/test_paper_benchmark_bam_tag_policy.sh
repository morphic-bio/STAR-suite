#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
paper_scripts=(
  "scripts/paper/run_a375_benchmark.sh"
  "scripts/paper/run_ucsf_ebs2_2_benchmark.sh"
  "scripts/paper/run_msk_30polyko_benchmark.sh"
)

for relative_path in "${paper_scripts[@]}"; do
  script="${repo_root}/${relative_path}"
  grep -q 'PARITY_BAM_TAG_POLICY' "${script}"
  if grep -n -- '--outSAMattributes' "${script}"; then
    echo "ERROR: canonical paper recipe requests BAM attributes: ${relative_path}" >&2
    exit 1
  fi
done

grep -q 'BAM Tag Policy for Comparisons' "${repo_root}/AGENTS.md"
grep -q 'BAM tags are not parity outputs' "${repo_root}/docs/PAPER_BENCHMARK_METHODOLOGY.md"
grep -q 'No `GX`/`UR` for count parity' "${repo_root}/docs/feature_barcodes.md"
grep -q 'Scoped exception' "${repo_root}/docs/RUNBOOK_POST_COLLAPSE_PRODUCTION_INTEGRATION_20260717.md"
grep -q 'GX/UR are raw adapter evidence here' \
  "${repo_root}/flex/tools/molecule_first_resolver/molecule_first_bam_ledger.cpp"

echo "paper benchmark BAM-tag policy: PASS"
