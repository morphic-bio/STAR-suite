#!/bin/bash
# Minimal bulk fixture test for --defaultBulk
# This script passes ONLY dataset/hardware args and relies entirely on the default bundle.
# It verifies that defaults are applied and output is produced.
#
# Usage: ./tests/run_default_bundle_bulk_fixture.sh
#        STAR_BIN=/path/to/STAR ./tests/run_default_bundle_bulk_fixture.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
# Use /tmp for test outputs (not tracked in repo)
TEST_DIR="/tmp/default_bundle_bulk_$$"

# Apply bulk defaults via the default bundle
export STAR_EXTRA_ARGS="--defaultBulk Yes"

echo "=== Minimal Bulk Default Bundle Fixture Test ==="
echo "STAR binary: $STAR_BIN"
echo "Using: $STAR_EXTRA_ARGS"
echo ""

# Ensure dependencies
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
[ -f "$STAR_BIN" ] || { echo "ERROR: STAR binary not found at $STAR_BIN"; exit 1; }

# Setup directories
REF_DIR="${TEST_DIR}/ref"
IDX_DIR="${REF_DIR}/star_index"
FASTQ_DIR="${TEST_DIR}/fastq"
OUT_DIR="${TEST_DIR}/output"
TMP_DIR="${TEST_DIR}/tmp"

rm -rf "$TEST_DIR"
mkdir -p "$REF_DIR" "$IDX_DIR" "$FASTQ_DIR" "$OUT_DIR"

# Create minimal reference
cat > "$REF_DIR/chr1.fa" << 'EOF'
>chr1
AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
EOF

cat > "$REF_DIR/genes.gtf" << 'EOF'
chr1	.	gene	1	60	.	+	.	gene_id "GENE1";
chr1	.	exon	1	60	.	+	.	gene_id "GENE1"; transcript_id "T1";
EOF

# Create minimal FASTQ (bulk RNA-seq style - just cDNA)
READ="AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
QUAL=$(printf 'I%.0s' $(seq 1 ${#READ}))

cat > "$FASTQ_DIR/reads.fastq" << EOF
@read1
$READ
+
$QUAL
@read2
$READ
+
$QUAL
@read3
$READ
+
$QUAL
EOF

# Build STAR index
echo "Building index..."
"$STAR_BIN" \
  --runMode genomeGenerate \
  --genomeDir "$IDX_DIR" \
  --genomeFastaFiles "$REF_DIR/chr1.fa" \
  --sjdbGTFfile "$REF_DIR/genes.gtf" \
  --sjdbOverhang 20 \
  --genomeSAindexNbases 4 >/dev/null

# Run STAR with ONLY dataset/hardware specific params
# All pipeline params should come from --defaultBulk
echo ""
echo "Running STAR with minimal args + --defaultBulk..."
echo "  Dataset args: --genomeDir, --readFilesIn, --outFileNamePrefix, --outTmpDir"
echo "  All other params from default bundle"
echo ""

# shellcheck disable=SC2086
"$STAR_BIN" ${STAR_EXTRA_ARGS:-} \
  --runThreadN 2 \
  --genomeDir "$IDX_DIR" \
  --readFilesIn "$FASTQ_DIR/reads.fastq" \
  --outFileNamePrefix "$OUT_DIR/" \
  --outTmpDir "$TMP_DIR"

# Verify defaults were applied
echo ""
echo "=== Verifying default bundle was applied ==="
if grep -q "Applying --defaultBulk" "$OUT_DIR/Log.out"; then
    echo "PASS: Found 'Applying --defaultBulk' in Log.out"
else
    echo "FAIL: Default bundle header not found in Log.out"
    exit 1
fi

if grep -q "\[default-group\]" "$OUT_DIR/Log.out"; then
    echo "PASS: Found [default-group] parameter lines in Log.out"
    grep "\[default-group\]" "$OUT_DIR/Log.out" | head -5
else
    echo "FAIL: No [default-group] parameters found in Log.out"
    exit 1
fi

# Verify outSAMtype was set by the default bundle
if grep -q "\[default-group\] outSAMtype = BAM SortedByCoordinate" "$OUT_DIR/Log.out"; then
    echo "PASS: outSAMtype BAM SortedByCoordinate was set by default bundle"
else
    echo "FAIL: outSAMtype was not set by default bundle"
    exit 1
fi

# Verify output - defaultBulk should produce sorted BAM
echo ""
echo "=== Verifying output ==="
if [ -f "$OUT_DIR/Aligned.sortedByCoord.out.bam" ]; then
    echo "PASS: Sorted BAM exists (Aligned.sortedByCoord.out.bam)"
    samtools quickcheck "$OUT_DIR/Aligned.sortedByCoord.out.bam"
    echo "PASS: BAM passes quickcheck"
else
    echo "FAIL: Sorted BAM missing"
    ls -la "$OUT_DIR/"
    exit 1
fi

echo ""
echo "=== Minimal Bulk Default Bundle Test PASSED ==="
echo "Defaults from --defaultBulk were successfully applied"

# Cleanup (optional, comment out to inspect outputs)
rm -rf "$TEST_DIR"
