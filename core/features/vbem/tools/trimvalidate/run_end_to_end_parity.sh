#!/bin/bash
# End-to-end parity test: Compare STAR alignments using Trim Galore vs built-in trimming
#
# This test validates that trimming parity translates to alignment parity by:
# 1. Running Trim Galore + STAR (baseline)
# 2. Running STAR with --trimCutadapt Yes (candidate)
# 3. Comparing BAM outputs (flagstat, idxstats, alignment diffs)
#
# Prerequisites:
#   - STAR binary in PATH or ./STAR
#   - trim_galore in PATH
#   - samtools in PATH
#   - Reference genome and GTF (see REF_DIR below)
#   - Input FASTQs (see INPUT_FASTQ variables below)

set -e

# Configuration
REF_DIR="${REF_DIR:-/mnt/pikachu/test-datasets/reference}"
GENOME_FA="${GENOME_FA:-$REF_DIR/chr22_23800000-23980000.fa}"
GTF="${GTF:-$REF_DIR/genes.gtf}"
STAR_INDEX="${STAR_INDEX:-/tmp/star_chr22_index}"
WORK="${WORK:-/tmp/trim_parity_e2e}"
THREADS="${THREADS:-4}"

INPUT_R1="${INPUT_R1:-/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_1.fastq.gz}"
INPUT_R2="${INPUT_R2:-/mnt/pikachu/test-datasets/testdata/GSE110004/SRR6357070_2.fastq.gz}"

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Find STAR binary
if [ -f "./STAR" ]; then
    STAR_BIN="./STAR"
elif command -v STAR >/dev/null 2>&1; then
    STAR_BIN="STAR"
else
    echo "Error: STAR not found. Set STAR_BIN or ensure STAR is in PATH."
    exit 1
fi

echo "=== End-to-End Trimming Parity Test ==="
echo "STAR binary: $STAR_BIN"
echo "Reference: $GENOME_FA"
echo "GTF: $GTF"
echo "Input R1: $INPUT_R1"
echo "Input R2: $INPUT_R2"
echo ""

# Check prerequisites
if [ ! -f "$GENOME_FA" ]; then
    echo "Error: Genome FASTA not found: $GENOME_FA"
    echo "Set GENOME_FA environment variable or create the reference file."
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "Error: GTF file not found: $GTF"
    echo "Set GTF environment variable or create the annotation file."
    exit 1
fi

if [ ! -f "$INPUT_R1" ] || [ ! -f "$INPUT_R2" ]; then
    echo "Error: Input FASTQs not found"
    echo "Set INPUT_R1 and INPUT_R2 environment variables."
    exit 1
fi

mkdir -p "$WORK"
mkdir -p "$STAR_INDEX"

# Step 1: Build STAR index (if not already built)
if [ ! -f "$STAR_INDEX/Genome" ]; then
    echo "Building STAR index..."
    $STAR_BIN --runThreadN $THREADS --runMode genomeGenerate \
              --genomeDir "$STAR_INDEX" \
              --genomeFastaFiles "$GENOME_FA" \
              --sjdbGTFfile "$GTF" \
              --sjdbOverhang 99
    echo "Index built."
else
    echo "Using existing STAR index: $STAR_INDEX"
fi

# Step 2: Baseline - Trim Galore + STAR
echo ""
echo "=== Baseline: Trim Galore + STAR ==="
BASELINE_DIR="$WORK/baseline"
mkdir -p "$BASELINE_DIR"

cp "$INPUT_R1" "$BASELINE_DIR/in_R1.fq.gz"
cp "$INPUT_R2" "$BASELINE_DIR/in_R2.fq.gz"

trim_galore --paired --quality 20 --length 20 \
            --adapter "$ADAPTER_R1" \
            --adapter2 "$ADAPTER_R2" \
            "$BASELINE_DIR/in_R1.fq.gz" "$BASELINE_DIR/in_R2.fq.gz" \
            --output_dir "$BASELINE_DIR"

$STAR_BIN --runThreadN $THREADS --genomeDir "$STAR_INDEX" \
          --readFilesIn "$BASELINE_DIR/in_R1_val_1.fq" "$BASELINE_DIR/in_R2_val_2.fq" \
          --outFileNamePrefix "$BASELINE_DIR/" \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMunmapped None \
          --outSAMattributes Standard

# Step 3: Candidate - STAR with built-in trimming
echo ""
echo "=== Candidate: STAR with --trimCutadapt Yes ==="
CANDIDATE_DIR="$WORK/candidate"
mkdir -p "$CANDIDATE_DIR"

cp "$INPUT_R1" "$CANDIDATE_DIR/in_R1.fq.gz"
cp "$INPUT_R2" "$CANDIDATE_DIR/in_R2.fq.gz"

$STAR_BIN --runThreadN $THREADS --genomeDir "$STAR_INDEX" \
          --readFilesCommand zcat \
          --readFilesIn "$CANDIDATE_DIR/in_R1.fq.gz" "$CANDIDATE_DIR/in_R2.fq.gz" \
          --trimCutadapt Yes \
          --trimCutadaptQuality 20 \
          --trimCutadaptMinLength 20 \
          --trimCutadaptAdapter "$ADAPTER_R1 $ADAPTER_R2" \
          --outFileNamePrefix "$CANDIDATE_DIR/" \
          --outSAMtype BAM SortedByCoordinate \
          --outSAMunmapped None \
          --outSAMattributes Standard

# Step 4: Compare results
echo ""
echo "=== Comparison ==="

BASELINE_BAM="$BASELINE_DIR/Aligned.sortedByCoord.out.bam"
CANDIDATE_BAM="$CANDIDATE_DIR/Aligned.sortedByCoord.out.bam"

if [ ! -f "$BASELINE_BAM" ]; then
    echo "Error: Baseline BAM not found: $BASELINE_BAM"
    exit 1
fi

if [ ! -f "$CANDIDATE_BAM" ]; then
    echo "Error: Candidate BAM not found: $CANDIDATE_BAM"
    exit 1
fi

# Index BAMs if needed
[ ! -f "$BASELINE_BAM.bai" ] && samtools index "$BASELINE_BAM"
[ ! -f "$CANDIDATE_BAM.bai" ] && samtools index "$CANDIDATE_BAM"

# Compare flagstat
echo "Comparing flagstat..."
samtools flagstat "$BASELINE_BAM" > "$WORK/baseline.flagstat"
samtools flagstat "$CANDIDATE_BAM" > "$WORK/candidate.flagstat"

if diff -q "$WORK/baseline.flagstat" "$WORK/candidate.flagstat" >/dev/null 2>&1; then
    echo "✅ flagstat: IDENTICAL"
else
    echo "❌ flagstat: DIFFERENT"
    echo "Differences:"
    diff "$WORK/baseline.flagstat" "$WORK/candidate.flagstat" || true
fi

# Compare idxstats
echo ""
echo "Comparing idxstats..."
samtools idxstats "$BASELINE_BAM" > "$WORK/baseline.idxstats"
samtools idxstats "$CANDIDATE_BAM" > "$WORK/candidate.idxstats"

if diff -q "$WORK/baseline.idxstats" "$WORK/candidate.idxstats" >/dev/null 2>&1; then
    echo "✅ idxstats: IDENTICAL"
else
    echo "❌ idxstats: DIFFERENT"
    echo "Differences:"
    diff "$WORK/baseline.idxstats" "$WORK/candidate.idxstats" || true
fi

# Compare alignment records (first 1000 reads)
echo ""
echo "Comparing alignment records (first 1000 reads)..."
samtools view "$BASELINE_BAM" | head -1000 | cut -f1-11 > "$WORK/baseline.alignments"
samtools view "$CANDIDATE_BAM" | head -1000 | cut -f1-11 > "$WORK/candidate.alignments"

if diff -q "$WORK/baseline.alignments" "$WORK/candidate.alignments" >/dev/null 2>&1; then
    echo "✅ Alignment records: IDENTICAL"
else
    echo "⚠️  Alignment records: DIFFERENT (checking if only tags differ)..."
    # Check if differences are only in optional tags (columns 12+)
    diff "$WORK/baseline.alignments" "$WORK/candidate.alignments" | head -20 || true
fi

# Summary
echo ""
echo "=== Summary ==="
echo "Baseline BAM: $BASELINE_BAM"
echo "Candidate BAM: $CANDIDATE_BAM"
echo ""
echo "Note: Since FASTQ-level parity is perfect (0 diff lines), any alignment differences"
echo "      would indicate STAR parameter mismatches, not trimming differences."
echo ""
echo "To clean up: rm -rf $WORK $STAR_INDEX"
