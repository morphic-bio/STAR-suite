#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tool_dir="${repo_root}/flex/tools/molecule_first_resolver"
resolver="${MOLECULE_FIRST_RESOLVER:-${tool_dir}/molecule_first_resolver}"
materializer="${MOLECULE_FIRST_MATERIALIZER:-${tool_dir}/molecule_first_materialize}"
fixture="${repo_root}/tests/molecule_first/hd_materializer_candidates.feature_sorted.tsv"
work="$(mktemp -d /tmp/star_suite_molecule_first_materializer_bounded.XXXXXX)"
trap 'rm -rf "${work}"' EXIT

if [[ ! -x "${resolver}" || ! -x "${materializer}" ]]; then
  make -C "${tool_dir}" all >/dev/null
fi

"${resolver}" --input "${fixture}" --out-dir "${work}/resolved" --input-feature-sorted
for run in one two; do
  "${materializer}" \
    --resolved-dir "${work}/resolved" \
    --out-dir "${work}/mex_${run}" \
    --assay visium-hd \
    --umi-mode 1mm_cr \
    --sort-memory-mb 1
done

diff -qr "${work}/mex_one" "${work}/mex_two"

cat > "${work}/expected_barcodes.tsv" <<'EOF'
s_002um_0_1-1
s_002um_1000_1-1
s_002um_10_1-1
s_002um_10_10-1
s_002um_1_1-1
EOF
cmp "${work}/expected_barcodes.tsv" "${work}/mex_one/strict/square_002um/barcodes.tsv"

canonical_features="$(stat -c '%d:%i' "${work}/mex_one/strict/square_002um/features.tsv")"
for policy in strict soft_expected hard gated_hard; do
  for scale in square_002um square_008um square_016um; do
    directory="${work}/mex_one/${policy}/${scale}"
    test -s "${directory}/matrix.mtx"
    [[ "$(stat -c '%d:%i' "${directory}/features.tsv")" == "${canonical_features}" ]]
    [[ "$(stat -c '%d:%i' "${directory}/barcodes.tsv")" == \
       "$(stat -c '%d:%i' "${work}/mex_one/strict/${scale}/barcodes.tsv")" ]]
  done
done

awk -F '\t' '
  NR > 8 {
    rows++
    if (!($1 in mass)) mass[$1]=$6
    delta=$6-mass[$1]
    if (delta < 0) delta=-delta
    if (delta > 1e-12) exit 1
    scales[$1]++
  }
  END {
    if (rows != 12) exit 1
    for (policy in scales) if (scales[policy] != 3) exit 1
  }
' "${work}/mex_one/summary.tsv"

if find "${work}/mex_one" -maxdepth 1 -name '.molecule_first_materialize.*' | grep -q .; then
  echo "ERROR: private materializer spill was not removed" >&2
  exit 1
fi

echo "bounded molecule-first materializer test: PASS"
