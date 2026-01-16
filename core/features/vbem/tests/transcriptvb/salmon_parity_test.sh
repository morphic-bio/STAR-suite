#!/bin/bash
# salmon_parity_test.sh - Run full Salmon parity test
#
# This test:
# 1. Runs STAR with TranscriptomeSAM to generate BAM
# 2. Runs Salmon on the BAM (alignment mode)
# 3. Runs STAR with TranscriptVB
# 4. Compares Salmon vs TranscriptVB results
#
# Usage: ./salmon_parity_test.sh [GENOME_DIR] [TRANSCRIPTOME_FASTA] [READS_1] [READS_2]

set -e

# Paths
STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-Flex/source/STAR}
SALMON_BIN=${SALMON_BIN:-salmon}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default test data
GENOME_DIR=${1:-/tmp/star_vb_test/star_new_index}
TRANSCRIPTOME=${2:-/mnt/pikachu/test-datasets-rnaseq/reference/transcriptome.fasta}
READS_1=${3:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz}
READS_2=${4:-/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_2.fastq.gz}

OUTDIR=/tmp/salmon_parity_test_$$
mkdir -p $OUTDIR
cd $OUTDIR

echo "========================================"
echo "Salmon Parity Test"
echo "========================================"
echo "STAR: $STAR_BIN"
echo "Salmon: $SALMON_BIN"
echo "Genome: $GENOME_DIR"
echo "Transcriptome: $TRANSCRIPTOME"
echo "Output: $OUTDIR"
echo "========================================"
echo ""

# Step 1: Run STAR with TranscriptomeSAM
echo "=== Step 1: Generate Transcriptome BAM ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --outSAMtype BAM Unsorted \
    --runThreadN 4 \
    --outFileNamePrefix trbam_

if [ ! -f trbam_Aligned.toTranscriptome.out.bam ]; then
    echo "ERROR: Transcriptome BAM not generated"
    exit 1
fi
echo "✓ Transcriptome BAM generated"
echo ""

# Step 2: Run Salmon on BAM
echo "=== Step 2: Run Salmon on BAM ==="
mkdir -p salmon_out
$SALMON_BIN quant \
    -t $TRANSCRIPTOME \
    -l A \
    -a trbam_Aligned.toTranscriptome.out.bam \
    -o salmon_out \
    2>&1 | tail -10

if [ ! -f salmon_out/quant.sf ]; then
    echo "ERROR: Salmon quant.sf not generated"
    exit 1
fi
echo "✓ Salmon quantification complete"
echo ""

# Step 3: Run STAR TranscriptVB
echo "=== Step 3: Run STAR TranscriptVB ==="
$STAR_BIN --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_1 $READS_2 \
    --readFilesCommand zcat \
    --quantMode TranscriptVB \
    --runThreadN 4 \
    --outFileNamePrefix star_vb_

if [ ! -f star_vb_quant.sf ]; then
    echo "ERROR: TranscriptVB quant.sf not generated"
    exit 1
fi
echo "✓ TranscriptVB quantification complete"
echo ""

# Check for gene-level output
if [ -f star_vb_quant.genes.sf ]; then
    echo "✓ quant.genes.sf generated"
fi
echo ""

# Step 4: Compare transcript-level results
echo "=== Step 4: Compare Transcript-Level Results ==="
if [ -f "$SCRIPT_DIR/compare_salmon_star.py" ]; then
    python3 "$SCRIPT_DIR/compare_salmon_star.py" salmon_out/quant.sf star_vb_quant.sf --verbose
    RESULT=$?
else
    # Fallback: basic comparison
    echo "compare_salmon_star.py not found, using basic comparison"
    python3 << 'EOF'
import pandas as pd
from scipy.stats import spearmanr, pearsonr

salmon = pd.read_csv('salmon_out/quant.sf', sep='\t')
star = pd.read_csv('star_vb_quant.sf', sep='\t')
merged = salmon.merge(star, on='Name', suffixes=('_salmon', '_star'))

expressed = merged[(merged['NumReads_salmon'] > 0) & (merged['NumReads_star'] > 0)]
r_spear, _ = spearmanr(expressed['NumReads_salmon'], expressed['NumReads_star'])
r_pear, _ = pearsonr(expressed['NumReads_salmon'], expressed['NumReads_star'])

print(f"Jointly expressed: {len(expressed)}")
print(f"Spearman: {r_spear:.4f}")
print(f"Pearson: {r_pear:.4f}")

if r_spear > 0.99:
    print("✓ PASS")
    exit(0)
else:
    print("✗ FAIL")
    exit(1)
EOF
    RESULT=$?
fi

# Step 5: Compare gene-level results (if available)
if [ -f star_vb_quant.genes.sf ] && [ -f salmon_out/quant.genes.sf ]; then
    echo ""
    echo "=== Step 5: Compare Gene-Level Results ==="
    if [ -f "$SCRIPT_DIR/compare_salmon_star.py" ]; then
        python3 "$SCRIPT_DIR/compare_salmon_star.py" salmon_out/quant.genes.sf star_vb_quant.genes.sf --verbose --gene
        GENE_RESULT=$?
        if [ $GENE_RESULT -ne 0 ]; then
            RESULT=$GENE_RESULT  # Fail if gene comparison fails
        fi
    else
        echo "compare_salmon_star.py not found, skipping gene-level comparison"
    fi
elif [ -f star_vb_quant.genes.sf ]; then
    echo ""
    echo "=== Step 5: Gene-Level Output ==="
    echo "✓ STAR quant.genes.sf generated"
    echo "  (Salmon quant.genes.sf not found - run Salmon with --geneMap to generate)"
fi

echo ""
echo "========================================"
if [ $RESULT -eq 0 ]; then
    echo "✓ Salmon Parity Test PASSED"
else
    echo "✗ Salmon Parity Test FAILED"
fi
echo "========================================"
echo ""
echo "Output files in: $OUTDIR"

exit $RESULT

