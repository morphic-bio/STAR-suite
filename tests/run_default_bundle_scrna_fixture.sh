#!/bin/bash
# Minimal scRNA fixture test for --defaultCoreScrna
# This script passes ONLY dataset/hardware args and relies entirely on the default bundle.
# It verifies that defaults are applied and output is produced.
#
# Usage: ./tests/run_default_bundle_scrna_fixture.sh
#        STAR_BIN=/path/to/STAR ./tests/run_default_bundle_scrna_fixture.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-${SCRIPT_DIR}/../core/legacy/source/STAR}"
# Use /tmp for test outputs (not tracked in repo)
TEST_DIR="/tmp/default_bundle_scrna_$$"

# Apply scRNA defaults via the default bundle
export STAR_EXTRA_ARGS="--defaultCoreScrna Yes"

echo "=== Minimal scRNA Default Bundle Fixture Test ==="
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
WL="${TEST_DIR}/whitelist.txt"

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

# Create whitelist
cat > "$WL" << 'EOF'
ACGTACGTACGTACGT
TGCATGCATGCATGCA
EOF

# Create minimal FASTQ
CB="ACGTACGTACGTACGT"
UMI="AAAAAAAAAAAA"
READ1="${CB}${UMI}"
READ2="AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"

q1=$(printf 'I%.0s' $(seq 1 ${#READ1}))
q2=$(printf 'I%.0s' $(seq 1 ${#READ2}))

cat > "$FASTQ_DIR/R1.fastq" << EOF
@read1/1
$READ1
+
$q1
EOF

cat > "$FASTQ_DIR/R2.fastq" << EOF
@read1/2
$READ2
+
$q2
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
# All pipeline params should come from --defaultCoreScrna (including soloType Droplet)
echo ""
echo "Running STAR with minimal args + --defaultCoreScrna..."
echo "  Dataset args: --genomeDir, --readFilesIn, --outFileNamePrefix, --outTmpDir"
echo "  Protocol args: --soloCB*, --soloUMI*, --soloCBwhitelist (CB/UMI geometry is dataset-specific)"
echo "  Default bundle provides: soloType, soloFeatures, soloCellFilter, etc."
echo ""

# shellcheck disable=SC2086
"$STAR_BIN" ${STAR_EXTRA_ARGS:-} \
  --runThreadN 2 \
  --genomeDir "$IDX_DIR" \
  --readFilesIn "$FASTQ_DIR/R2.fastq" "$FASTQ_DIR/R1.fastq" \
  --outFileNamePrefix "$OUT_DIR/" \
  --outTmpDir "$TMP_DIR" \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist "$WL"

# Verify defaults were applied
echo ""
echo "=== Verifying default bundle was applied ==="
if grep -q "Applying --defaultCoreScrna" "$OUT_DIR/Log.out"; then
    echo "PASS: Found 'Applying --defaultCoreScrna' in Log.out"
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

# Verify soloType was set by the default bundle (not overridden)
if grep -q "\[default-group\] soloType = Droplet" "$OUT_DIR/Log.out"; then
    echo "PASS: soloType Droplet was set by default bundle"
else
    echo "FAIL: soloType was not set by default bundle"
    exit 1
fi

# Verify output
echo ""
echo "=== Verifying output ==="
if [ -f "$OUT_DIR/Solo.out/Gene/raw/matrix.mtx" ]; then
    echo "PASS: Solo matrix.mtx exists"
else
    echo "FAIL: Solo matrix.mtx missing"
    exit 1
fi

echo ""
echo "=== Minimal scRNA Default Bundle Test PASSED ==="
echo "Defaults from --defaultCoreScrna were successfully applied"

# Cleanup (optional, comment out to inspect outputs)
rm -rf "$TEST_DIR"
