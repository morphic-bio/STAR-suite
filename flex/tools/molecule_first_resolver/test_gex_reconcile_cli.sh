#!/usr/bin/env bash
set -euo pipefail

reconciler="${1:?GEX reconciler path required}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

mkdir "$tmp/source"
printf 'read_clique_id\n' >"$tmp/source/read_cliques.tsv"
printf 'read_id\n' >"$tmp/source/hard_call_audit.tsv"

cat >"$tmp/support.tsv" <<'EOF'
umi_mode	product	candidate	corrected_umi	feature_id	corrected_count	original_at_corrected_count	expected_count	original_expected_count	molecule_id	member_read_count	member_read_ids	read_clique_ids
1mm_cr	hard	A	AAAA	G1	5	3			m1	5	r1;r2;r3;r4;r5	c1
1mm_cr	hard	A	AAAA	G2	3	3			m2	3	r6;r7;r8	c2
1mm_cr	hard	B	AAAA	G1	2	1			m3	2	r9;r10	c3
1mm_cr	hard	B	AAAA	G2	2	1			m4	2	r11;r12	c4
1mm_cr	hard	C	AAAA	G1	4	1			m5	4	r13;r14;r15;r16	c5
1mm_cr	hard	C	AAAA	G2	3	2			m6	3	r17;r18;r19	c6
1mm_cr	soft_expected	A	AAAC	G1			0.8	0.7				c7
1mm_cr	soft_expected	A	AAAC	G2			0.2	0.1				c8
1mm_cr	soft_expected	B	AAAC	G1			0.5	0.4				c9
1mm_cr	soft_expected	B	AAAC	G2			0.5	0.3				c10
1mm_cr	soft_expected	C	AAAC	G1			0.7	0.1				c11
1mm_cr	soft_expected	C	AAAC	G2			0.2	0.2				c12
EOF

"$reconciler" --input "$tmp/support.tsv" --resolved-source "$tmp/source" \
  --out-dir "$tmp/out"

cat >"$tmp/expected-hard.tsv" <<'EOF'
umi_mode	product	molecule_id	feature_id	corrected_umi	candidate	member_read_count	member_read_ids	read_clique_ids
1mm_cr	hard	m1	G1	AAAA	A	5	r1;r2;r3;r4;r5	c1
EOF
diff -u "$tmp/expected-hard.tsv" "$tmp/out/hard_molecules.tsv"

cat >"$tmp/expected-soft.tsv" <<'EOF'
umi_mode	feature_id	corrected_umi	candidate	expected_count	read_clique_ids
1mm_cr	G1	AAAC	A	0.80000000000000004	c7
EOF
diff -u "$tmp/expected-soft.tsv" "$tmp/out/soft_expected_molecules.tsv"

for suffix in \
  multigene_groups:3 \
  multigene_accepted:1 \
  multigene_corrected_count_ties:1 \
  multigene_original_umi_dominance_rejected:1; do
  metric="${suffix%%:*}"
  value="${suffix##*:}"
  grep -q "^1mm_cr.hard.${metric}$(printf '\t')${value}$" "$tmp/out/summary.tsv"
  grep -q "^1mm_cr.soft_expected.${metric}$(printf '\t')${value}$" \
    "$tmp/out/summary.tsv"
done

[[ "$tmp/source/read_cliques.tsv" -ef "$tmp/out/read_cliques.tsv" ]]
[[ "$tmp/source/hard_call_audit.tsv" -ef "$tmp/out/hard_call_audit.tsv" ]]

if "$reconciler" --input "$tmp/support.tsv" --resolved-source "$tmp/source" \
  --out-dir "$tmp/out" 2>/dev/null; then
  echo "nonempty output directory unexpectedly accepted" >&2
  exit 1
fi

echo "bounded GEX reconciliation CLI tests passed"
