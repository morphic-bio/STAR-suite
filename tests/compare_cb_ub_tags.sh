#!/bin/bash
# Compare CB/UB tags in sorted BAM against authoritative TSV table
# Usage: compare_cb_ub_tags.sh <bam> <tsv>
# Returns 0 on success, 1 on mismatch

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <bam> <tsv>" >&2
    exit 1
fi

BAM="$1"
TSV="$2"

# Check inputs exist
if [ ! -f "$BAM" ]; then
    echo "ERROR: BAM file not found: $BAM" >&2
    exit 1
fi

if [ ! -f "$TSV" ]; then
    echo "ERROR: TSV table not found: $TSV" >&2
    exit 1
fi

# Check TSV is not empty
if [ ! -s "$TSV" ]; then
    echo "ERROR: TSV table is empty: $TSV" >&2
    exit 1
fi

# Check samtools is available
if ! command -v samtools &> /dev/null; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# 1. Extract readId, CB, UB from BAM (treat empty string as "-")
# Note: Duplicate readIds are expected due to multi-mappers (same readId can appear in multiple alignments)
# Use Python for reliable numeric sorting and zero-padding for join compatibility
samtools view "$BAM" | python3 -c "
import sys
data = []
for line in sys.stdin:
    fields = line.strip().split('\t')
    zi = ''
    cb = '-'
    ub = '-'
    for f in fields[11:]:  # Tags start at field 11 (0-indexed)
        if f.startswith('ZI:i:'):
            zi = f[5:]
        elif f.startswith('CB:Z:'):
            cb = f[5:]
        elif f.startswith('UB:Z:'):
            ub = f[5:]
    if zi == '':
        print('ERROR: BAM record missing ZI tag:', fields[0], file=sys.stderr)
        sys.exit(1)
    data.append((int(zi), cb, ub))
data.sort(key=lambda x: x[0])
# Zero-pad to 10 digits for lexicographic sort compatibility with join
for zi, cb, ub in data:
    print(f'{zi:010d}\t{cb}\t{ub}')
" > "$TMPDIR/bam_tags.tsv"

# 2. Sort TSV table by readId (skip header, zero-pad for join compatibility)
python3 -c "
import sys
data = []
with open('$TSV') as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 4:
            data.append((int(parts[0]), parts[1], parts[2], int(parts[3])))
data.sort(key=lambda x: x[0])
# Zero-pad to 10 digits for lexicographic sort compatibility with join
for rid, cb, ub, status in data:
    print(f'{rid:010d}\t{cb}\t{ub}\t{status}')
" > "$TMPDIR/expected_tags.tsv"

# 3. Join and compare
# Expected format: readId cb ub status
# BAM format: readId cb ub
# Note: join will create multiple output lines if readId appears multiple times in BAM (multi-mappers),
#       but each alignment is validated against the same expected table entry
join -t$'\t' -1 1 -2 1 "$TMPDIR/bam_tags.tsv" "$TMPDIR/expected_tags.tsv" | awk -F'\t' '
BEGIN { 
    mismatches=0
    total=0
}
{
    readId=$1+0  # Convert zero-padded back to integer
    bam_cb=$2; bam_ub=$3
    exp_cb=$4; exp_ub=$5; status=$6
    
    total++
    
    # Treat empty BAM tags as "-"
    if (bam_cb == "") bam_cb = "-"
    if (bam_ub == "") bam_ub = "-"
    
    # Validation per status
    if (status == 0) {
        # Both should be absent
        if (bam_cb != "-" || bam_ub != "-") {
            print "MISMATCH readId="readId" status=0 expected CB/UB absent, got CB="bam_cb" UB="bam_ub
            mismatches++
        }
    } else if (status == 1) {
        # Both should match exactly
        if (bam_cb != exp_cb || bam_ub != exp_ub) {
            print "MISMATCH readId="readId" status=1 expected CB="exp_cb" UB="exp_ub" got CB="bam_cb" UB="bam_ub
            mismatches++
        }
    } else if (status == 2) {
        # CB should match, UB should be absent
        if (bam_cb != exp_cb || bam_ub != "-") {
            print "MISMATCH readId="readId" status=2 expected CB="exp_cb" UB=- got CB="bam_cb" UB="bam_ub
            mismatches++
        }
    }
}
END {
    if (mismatches > 0) {
        print "FAILED: "mismatches" mismatches out of "total" records" > "/dev/stderr"
        exit 1
    }
    print "PASSED: All CB/UB tags match expected values ("total" records)"
}
'

# 4. Verify all BAM readIds appear in TSV (check for missing readIds)
missing_readids=$(join -t$'\t' -v1 -1 1 -2 1 "$TMPDIR/bam_tags.tsv" "$TMPDIR/expected_tags.tsv" | cut -f1 | sort -u | wc -l)
if [ "$missing_readids" -gt 0 ]; then
    echo "ERROR: Found $missing_readids unique readIds in BAM that are missing from TSV table" >&2
    join -t$'\t' -v1 -1 1 -2 1 "$TMPDIR/bam_tags.tsv" "$TMPDIR/expected_tags.tsv" | cut -f1 | sort -u | head -10 | while read rid; do
        echo "  Missing readId: $rid" >&2
    done
    if [ "$missing_readids" -gt 10 ]; then
        echo "  ... and $((missing_readids - 10)) more" >&2
    fi
    exit 1
fi

# 5. Report statistics
tsv_count=$(wc -l < "$TMPDIR/expected_tags.tsv")
bam_count=$(wc -l < "$TMPDIR/bam_tags.tsv")
unique_readids=$(cut -f1 "$TMPDIR/bam_tags.tsv" | sort -u | wc -l)
echo "Statistics: TSV has $tsv_count reads, BAM has $bam_count alignments ($unique_readids unique readIds)" >&2
if [ "$bam_count" -gt "$tsv_count" ]; then
    echo "Note: Multi-mappers detected ($((bam_count - tsv_count)) additional alignments)" >&2
fi

