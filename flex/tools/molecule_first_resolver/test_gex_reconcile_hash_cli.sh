#!/usr/bin/env bash
set -euo pipefail

stream_reconciler="${1:?streaming GEX reconciler path required}"
hash_reconciler="${2:?hash-partitioned GEX reconciler path required}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

mkdir "$tmp/source" "$tmp/partitions"
printf 'read_clique_id\n' >"$tmp/source/read_cliques.tsv"
printf 'read_id\n' >"$tmp/source/hard_call_audit.tsv"

# This is valid feature order and deliberately invalid reconciliation-key order.
cat >"$tmp/feature-support.tsv" <<'EOF'
umi_mode	product	candidate	corrected_umi	feature_id	corrected_count	original_at_corrected_count	expected_count	original_expected_count	molecule_id	member_read_count	member_read_ids	read_clique_ids
1mm_cr	hard	A	AAAA	G1	5	3			m1	5	r1;r2;r3;r4;r5	c1
1mm_cr	hard	B	AAAA	G1	2	1			m3	2	r9;r10	c3
1mm_cr	hard	C	AAAA	G1	4	1			m5	4	r13;r14;r15;r16	c5
1mm_cr	soft_expected	A	AAAC	G1			0.8	0.7				c7
1mm_cr	soft_expected	B	AAAC	G1			0.5	0.4				c9
1mm_cr	soft_expected	C	AAAC	G1			0.7	0.1				c11
1mm_cr	hard	A	AAAA	G2	3	3			m2	3	r6;r7;r8	c2
1mm_cr	hard	B	AAAA	G2	2	1			m4	2	r11;r12	c4
1mm_cr	hard	C	AAAA	G2	3	2			m6	3	r17;r18;r19	c6
1mm_cr	soft_expected	A	AAAC	G2			0.2	0.1				c8
1mm_cr	soft_expected	B	AAAC	G2			0.5	0.3				c10
1mm_cr	soft_expected	C	AAAC	G2			0.2	0.2				c12
EOF

{
  head -n 1 "$tmp/feature-support.tsv"
  tail -n +2 "$tmp/feature-support.tsv" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4 -k5,5
} >"$tmp/key-support.tsv"

"$stream_reconciler" --input "$tmp/key-support.tsv" --resolved-source "$tmp/source" \
  --out-dir "$tmp/stream"

for repeat in one two; do
  "$hash_reconciler" --input "$tmp/feature-support.tsv" --input-feature-sorted \
    --resolved-source "$tmp/source" --out-dir "$tmp/hash-$repeat" \
    --tmp-dir "$tmp/partitions" --hash-partitions 3 --partition-buffer-kb 1
done

for table in strict_molecules.tsv hard_molecules.tsv gated_hard_molecules.tsv \
  soft_expected_molecules.tsv; do
  {
    head -n 1 "$tmp/stream/$table"
    tail -n +2 "$tmp/stream/$table" | LC_ALL=C sort
  } >"$tmp/stream.$table"
  {
    head -n 1 "$tmp/hash-one/$table"
    tail -n +2 "$tmp/hash-one/$table" | LC_ALL=C sort
  } >"$tmp/hash.$table"
  diff -u "$tmp/stream.$table" "$tmp/hash.$table"
  cmp "$tmp/hash-one/$table" "$tmp/hash-two/$table"
done

# Compare all scientific summary metrics, excluding hash execution metadata.
tail -n +2 "$tmp/stream/summary.tsv" | LC_ALL=C sort >"$tmp/stream.summary"
tail -n +2 "$tmp/hash-one/summary.tsv" | grep -v '^hash\.' | LC_ALL=C sort \
  >"$tmp/hash.summary"
diff -u "$tmp/stream.summary" "$tmp/hash.summary"
cmp "$tmp/hash-one/summary.tsv" "$tmp/hash-two/summary.tsv"
cmp "$tmp/hash-one/resolved_config.tsv" "$tmp/hash-two/resolved_config.tsv"

[[ "$tmp/source/read_cliques.tsv" -ef "$tmp/hash-one/read_cliques.tsv" ]]
[[ "$tmp/source/hard_call_audit.tsv" -ef "$tmp/hash-one/hard_call_audit.tsv" ]]
if find "$tmp/partitions" -mindepth 1 -print -quit | grep -q .; then
  echo "processed hash partitions were not removed" >&2
  exit 1
fi

if "$hash_reconciler" --input "$tmp/feature-support.tsv" --input-feature-sorted \
  --resolved-source "$tmp/source" --out-dir "$tmp/hash-one" \
  --tmp-dir "$tmp/partitions" --hash-partitions 3 2>/dev/null; then
  echo "nonempty output directory unexpectedly accepted" >&2
  exit 1
fi

{
  head -n 1 "$tmp/feature-support.tsv"
  sed -n '8p' "$tmp/feature-support.tsv"
  sed -n '2p' "$tmp/feature-support.tsv"
} >"$tmp/not-feature-sorted.tsv"
if "$hash_reconciler" --input "$tmp/not-feature-sorted.tsv" --input-feature-sorted \
  --resolved-source "$tmp/source" --out-dir "$tmp/out-of-order" \
  --tmp-dir "$tmp/partitions" --hash-partitions 3 2>/dev/null; then
  echo "out-of-order feature input unexpectedly accepted" >&2
  exit 1
fi

echo "hash-partitioned GEX reconciliation CLI tests passed"
