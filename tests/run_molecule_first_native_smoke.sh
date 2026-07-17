#!/usr/bin/env bash
set -euo pipefail

resolver="${MOLECULE_FIRST_RESOLVER:-/usr/local/bin/molecule_first_resolver}"
fixture="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/molecule_first/candidate_reads.tsv"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

if [[ ! -x "${resolver}" ]]; then
  echo "ERROR: molecule-first resolver is not executable: ${resolver}" >&2
  exit 2
fi

"${resolver}" --input "${fixture}" --out-dir "${tmp}/out"

expected='gated_hard_molecules.tsv
hard_call_audit.tsv
hard_molecules.tsv
read_cliques.tsv
resolved_config.tsv
soft_expected_molecules.tsv
strict_molecules.tsv
summary.tsv'
actual="$(find "${tmp}/out" -maxdepth 1 -type f -printf '%f\n' | LC_ALL=C sort)"
[[ "${actual}" == "${expected}" ]]

grep -q $'^prior_application\tonce_per_read_clique$' "${tmp}/out/resolved_config.tsv"
grep -q $'^spatial_lambda\t0$' "${tmp}/out/resolved_config.tsv"
awk -F '\t' '
  $1 == "1mm_cr.soft_occupancy_mass" {
    found = 1
    delta = $2 - 3.3077118714766973
    if (delta < 0) delta = -delta
    if (delta > 1e-12) exit 1
  }
  END { if (!found) exit 1 }
' "${tmp}/out/summary.tsv"

mkdir "${tmp}/nonempty"
touch "${tmp}/nonempty/existing"
if "${resolver}" --input "${fixture}" --out-dir "${tmp}/nonempty" 2>/dev/null; then
  echo "ERROR: resolver accepted a nonempty output directory" >&2
  exit 1
fi

echo "molecule-first native smoke passed"
