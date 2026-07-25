#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tool_dir="${repo_root}/flex/tools/molecule_first_resolver"
ledger="${tool_dir}/molecule_first_bam_ledger"
resolver="${tool_dir}/molecule_first_resolver"
materialize="${tool_dir}/molecule_first_materialize"
fixture="${repo_root}/tests/molecule_first/raw_tag_fixture.sam"
whitelist="${repo_root}/tests/molecule_first/raw_tag_whitelist.txt"
work="$(mktemp -d)"
trap 'rm -rf "${work}"' EXIT

make -C "${tool_dir}" all >/dev/null

"${ledger}" \
  --input "${fixture}" \
  --whitelist "${whitelist}" \
  --output "${work}/candidate_reads.tsv" \
  --summary "${work}/ledger_summary.tsv" \
  --assay scrna

[[ "$(tail -n +2 "${work}/candidate_reads.tsv" | wc -l)" -eq 7 ]]
awk -F '\t' 'NR > 1 && $4 == "AAAA" && $6 != 2 {exit 1}
             NR > 1 && $4 == "AAAC" && $6 != 1 {exit 1}
             END {if (NR != 8) exit 1}' "${work}/candidate_reads.tsv"
[[ "$(awk -F '\t' 'NR > 1 && $1 == "r_ambiguous" {n++} END {print n+0}' "${work}/candidate_reads.tsv")" -eq 2 ]]
[[ "$(awk -F '\t' 'NR > 1 && $1 == "r_nbase" {n++} END {print n+0}' "${work}/candidate_reads.tsv")" -eq 2 ]]
grep -q $'^corrected_tags_used\tfalse$' "${work}/ledger_summary.tsv"
grep -q $'^called_cell_fields_used\tfalse$' "${work}/ledger_summary.tsv"

sed -e 's/CB:Z:[^[:space:]]*/CB:Z:CHANGED-1/g' \
    -e 's/UB:Z:[^[:space:]]*/UB:Z:CHANGED/g' \
    "${fixture}" > "${work}/corrected_tags_changed.sam"
"${ledger}" \
  --input "${work}/corrected_tags_changed.sam" \
  --whitelist "${whitelist}" \
  --output "${work}/candidate_reads_changed.tsv" \
  --summary "${work}/ledger_summary_changed.tsv" \
  --assay scrna
cmp "${work}/candidate_reads.tsv" "${work}/candidate_reads_changed.tsv"

if "${ledger}" --input "${fixture}" --whitelist "${whitelist}" \
    --output "${work}/flex_without_sample.tsv" --summary "${work}/flex_without_sample.summary.tsv" \
    --assay flex >/dev/null 2>&1; then
  echo "ERROR: Flex adapter accepted a ledger without a sample boundary" >&2
  exit 1
fi
"${ledger}" --input "${fixture}" --whitelist "${whitelist}" \
  --output "${work}/flex.tsv" --summary "${work}/flex.summary.tsv" \
  --assay flex --sample-id sample_a

grep -v '^@PG' "${fixture}" > "${work}/not_star.sam"
if "${ledger}" --input "${work}/not_star.sam" --whitelist "${whitelist}" \
    --output "${work}/not_star.tsv" --summary "${work}/not_star.summary.tsv" \
    --assay scrna >/dev/null 2>&1; then
  echo "ERROR: adapter accepted input without STAR provenance" >&2
  exit 1
fi
sed 's/ID:STAR/ID:STARLIKE/; s/PN:STAR/PN:STARLIKE/' "${fixture}" > "${work}/not_exact_star.sam"
if "${ledger}" --input "${work}/not_exact_star.sam" --whitelist "${whitelist}" \
    --output "${work}/not_exact_star.tsv" --summary "${work}/not_exact_star.summary.tsv" \
    --assay scrna >/dev/null 2>&1; then
  echo "ERROR: adapter accepted a non-STAR program-name prefix" >&2
  exit 1
fi

"${resolver}" --input "${work}/candidate_reads.tsv" --out-dir "${work}/resolved"
"${materialize}" --resolved-dir "${work}/resolved" --out-dir "${work}/scrna_mex" \
  --assay scrna --umi-mode 1mm_cr
for policy in strict soft_expected hard gated_hard; do
  test -s "${work}/scrna_mex/${policy}/raw/matrix.mtx"
  cmp "${work}/scrna_mex/strict/raw/features.tsv" "${work}/scrna_mex/${policy}/raw/features.tsv"
  cmp "${work}/scrna_mex/strict/raw/barcodes.tsv" "${work}/scrna_mex/${policy}/raw/barcodes.tsv"
done
grep -q '^%%MatrixMarket matrix coordinate real general$' "${work}/scrna_mex/soft_expected/raw/matrix.mtx"
grep -q '^%%MatrixMarket matrix coordinate integer general$' "${work}/scrna_mex/hard/raw/matrix.mtx"
grep -q $'^soft_cell_calling_allowed\tfalse$' "${work}/scrna_mex/summary.tsv"
if soft_error="$("${repo_root}/scripts/run_molecule_first_cell_calling.sh" \
    --star "${repo_root}/core/legacy/source/STAR" \
    --policy-root "${work}/scrna_mex" --out-root "${work}/forbidden_soft" \
    --policy soft_expected 2>&1)"; then
  echo "ERROR: soft expected counts were accepted for integer EmptyDrops" >&2
  exit 1
fi
grep -q 'soft_expected is real-valued' <<<"${soft_error}"

cat > "${work}/hd_candidates.tsv" <<'HD'
read_id	feature_id	raw_umi	candidate	log_sequence_likelihood	exact_read_count
h1	geneA	U001	s_002um_4_4	0	4
h2	geneA	U002	s_002um_4_4	-0.1	4
h2	geneA	U002	s_002um_5_5	-1.1	2
h3	geneB	U003	s_002um_5_5	0	2
HD
"${resolver}" --input "${work}/hd_candidates.tsv" --out-dir "${work}/hd_resolved"
"${materialize}" --resolved-dir "${work}/hd_resolved" --out-dir "${work}/hd_mex" \
  --assay visium-hd --umi-mode 1mm_cr
for policy in strict soft_expected hard gated_hard; do
  for scale in square_002um square_008um square_016um; do
    test -s "${work}/hd_mex/${policy}/${scale}/matrix.mtx"
  done
  awk -F '\t' -v policy="${policy}" '
    $1 == policy {if (seen++ == 0) mass=$6; else if (($6-mass > 1e-9) || (mass-$6 > 1e-9)) exit 1}
    END {if (seen != 3) exit 1}' "${work}/hd_mex/summary.tsv"
done

echo "molecule-first production adapters: PASS"
