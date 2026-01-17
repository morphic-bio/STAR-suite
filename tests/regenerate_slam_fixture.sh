#!/bin/bash
#
# Regenerate SLAM-seq fixture data for parity testing
#
# This script rebuilds the entire fixture pipeline:
# 1. Download 100k reads from SRR32576116 (SE150)
# 2. Align with STAR
# 3. Run GRAND-SLAM SNP detection
# 4. Create SNP BED from snpdata
# 5. Prefilter BAM by removing reads overlapping SNPs
# 6. Run GRAND-SLAM quant on prefiltered BAM
# 7. Verify output matches existing fixture
#
# Usage:
#   bash regenerate_slam_fixture.sh [--verify-only] [--skip-download]
#
# Requirements:
#   - sra-tools (fastq-dump)
#   - samtools
#   - GEDI/GRAND-SLAM
#   - STAR
#   - Python 3

set -euo pipefail

# Parse arguments
VERIFY_ONLY=0
SKIP_DOWNLOAD=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify-only) VERIFY_ONLY=1; shift ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
FIXTURE_ROOT="$PROJECT_ROOT/test/fixtures/slam"
WORK_DIR="$PROJECT_ROOT/test/tmp_slam_fixture"
STAR_BIN="$PROJECT_ROOT/core/legacy/source/STAR"
GEDI_BIN="$PROJECT_ROOT/gedi"

# Reference paths (from source.txt)
STAR_INDEX="/storage/autoindex_110_44/bulk_index"
GENOME_FA="/storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa"
GENES_GTF="/storage/autoindex_110_44/bulk_index/cellranger_ref/genes.gtf"
GEDI_GENOME="homo_sapiens_110_44"

# Fixture files
RAW_FASTQ="$FIXTURE_ROOT/raw/slam_100000_reads_SRR32576116.fastq.gz"
SNPS_BED="$FIXTURE_ROOT/ref/snps.bed"
REF_TSV="$FIXTURE_ROOT/expected/fixture_ref_human.tsv.gz"
CHECKSUM_FILE="$FIXTURE_ROOT/meta/checksums.sha256"

# SRA accession
SRA_ACCESSION="SRR32576116"
READ_COUNT=100000

# Output prefixes
STAR_PREFIX="$WORK_DIR/fixture_human_"
GEDI_PREFIX="$WORK_DIR/fixture_ref_human"

echo "========================================================================"
echo "SLAM Fixture Regeneration"
echo "========================================================================"
echo "Date: $(date)"
echo "Work directory: $WORK_DIR"
echo "Fixture root: $FIXTURE_ROOT"
echo "STAR index: $STAR_INDEX"
echo "GEDI genome: $GEDI_GENOME"
echo "========================================================================"
echo

# Verify binaries
if [[ ! -x "$STAR_BIN" ]]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi
if [[ ! -x "$GEDI_BIN" ]]; then
    echo "ERROR: GEDI binary not found: $GEDI_BIN"
    exit 1
fi
if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools not found in PATH"
    exit 1
fi

# Verify references
if [[ ! -d "$STAR_INDEX" ]]; then
    echo "ERROR: STAR index not found: $STAR_INDEX"
    exit 1
fi
if [[ ! -f "$GENOME_FA" ]]; then
    echo "ERROR: Genome FASTA not found: $GENOME_FA"
    exit 1
fi
if [[ ! -f "$GENES_GTF" ]]; then
    echo "ERROR: GTF not found: $GENES_GTF"
    exit 1
fi

mkdir -p "$WORK_DIR"

# Step 1: Download/extract FASTQ
if [[ $SKIP_DOWNLOAD -eq 0 ]] && [[ $VERIFY_ONLY -eq 0 ]]; then
    echo "[1/7] Downloading reads from $SRA_ACCESSION..."
    
    if command -v fastq-dump &>/dev/null; then
        # Download with SRA toolkit
        fastq-dump -X "$READ_COUNT" --gzip -O "$WORK_DIR" "$SRA_ACCESSION" || {
            echo "WARNING: fastq-dump failed, using existing FASTQ if available"
        }
        
        # Rename if needed
        if [[ -f "$WORK_DIR/${SRA_ACCESSION}.fastq.gz" ]]; then
            # Take first 100k reads (400k lines)
            zcat "$WORK_DIR/${SRA_ACCESSION}.fastq.gz" | head -n 400000 | gzip > \
                "$WORK_DIR/slam_100000_reads_${SRA_ACCESSION}.fastq.gz"
            echo "✓ Downloaded and extracted $READ_COUNT reads"
        fi
    else
        echo "WARNING: fastq-dump not available, using existing fixture FASTQ"
    fi
fi

# Use fixture FASTQ if work copy doesn't exist
FASTQ_FILE="$WORK_DIR/slam_100000_reads_${SRA_ACCESSION}.fastq.gz"
if [[ ! -f "$FASTQ_FILE" ]]; then
    if [[ -f "$RAW_FASTQ" ]]; then
        FASTQ_FILE="$RAW_FASTQ"
        echo "Using existing fixture FASTQ: $FASTQ_FILE"
    else
        echo "ERROR: No FASTQ available at $FASTQ_FILE or $RAW_FASTQ"
        exit 1
    fi
fi

if [[ $VERIFY_ONLY -eq 1 ]]; then
    echo "=== Verify-only mode: checking existing fixtures ==="
    goto_verify=1
else
    goto_verify=0
fi

if [[ $goto_verify -eq 0 ]]; then
    # Step 2: Align with STAR
    # Use adapter clipping (same as original fixture generation)
    echo "[2/7] Aligning reads with STAR..."
    "$STAR_BIN" \
        --runThreadN 4 \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ_FILE" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$STAR_PREFIX" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --clip3pAdapterMMp 0.1 \
        > "$WORK_DIR/star_align.log" 2>&1
    
    BAM_FILE="${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "ERROR: STAR alignment failed - BAM not found"
        cat "$WORK_DIR/star_align.log"
        exit 1
    fi
    
    # Index BAM
    samtools index "$BAM_FILE"
    echo "✓ STAR alignment complete: $BAM_FILE"

    # Step 3: Run GRAND-SLAM SNP detection
    echo "[3/7] Running GRAND-SLAM SNP detection..."
    SNPDATA_PREFIX="$WORK_DIR/snpdetect"
    
    # Clean up previous GEDI outputs (GEDI skips if output exists)
    rm -f "$SNPDATA_PREFIX".* 2>/dev/null || true
    
    # Run GRAND-SLAM with SNP detection enabled
    # Use -keep to preserve temp files (including snpdata)
    "$GEDI_BIN" -e Slam \
        -reads "$BAM_FILE" \
        -genomic "$GEDI_GENOME" \
        -prefix "$SNPDATA_PREFIX" \
        -strandness Sense \
        -snpConv 0.3 \
        -snppval 0.001 \
        -nthreads 4 \
        -progress \
        -keep \
        > "$WORK_DIR/gedi_snp.log" 2>&1 || {
        echo "WARNING: GEDI SNP detection had errors, check log"
        tail -20 "$WORK_DIR/gedi_snp.log"
    }
    
    echo "✓ GRAND-SLAM SNP detection complete"

    # Step 4: Convert snpdata to BED
    echo "[4/7] Converting SNP data to BED format..."
    SNPDATA_FILE="$SNPDATA_PREFIX.snpdata"
    NEW_SNPS_BED="$WORK_DIR/snps.bed"
    
    if [[ -f "$SNPDATA_FILE" ]]; then
        # snpdata format: Location (chrom:pos 1-based), Coverage, Mismatches, P value
        # Example: 1:11993237   9.0   9.0   1.37e-05
        # Convert to BED3 format (0-based, half-open) with chr prefix
        # Note: Keep duplicates (no uniq) to match original fixture
        tail -n +2 "$SNPDATA_FILE" | awk -F'\t' '{
            split($1, loc, ":");
            chrom = loc[1];
            pos = loc[2];
            # Handle chromosome naming: MT -> chrM, others get chr prefix
            if (chrom == "MT") chrom = "chrM";
            else if (chrom !~ /^chr/) chrom = "chr" chrom;
            # Convert to 0-based BED (pos-1 to pos)
            print chrom "\t" (pos-1) "\t" pos
        }' > "$NEW_SNPS_BED"
        SNP_COUNT=$(wc -l < "$NEW_SNPS_BED")
        echo "✓ Generated SNP BED: $SNP_COUNT entries"
    else
        echo "WARNING: snpdata file not found at $SNPDATA_FILE"
        # Try to use existing snpdata from previous run
        EXISTING_SNPDATA="$WORK_DIR/snpdetect_human.snpdata"
        if [[ -f "$EXISTING_SNPDATA" ]]; then
            echo "Using existing snpdata: $EXISTING_SNPDATA"
            tail -n +2 "$EXISTING_SNPDATA" | awk -F'\t' '{
                split($1, loc, ":");
                chrom = loc[1];
                pos = loc[2];
                if (chrom == "MT") chrom = "chrM";
                else if (chrom !~ /^chr/) chrom = "chr" chrom;
                print chrom "\t" (pos-1) "\t" pos
            }' > "$NEW_SNPS_BED"
            SNP_COUNT=$(wc -l < "$NEW_SNPS_BED")
            echo "✓ Generated SNP BED from existing snpdata: $SNP_COUNT entries"
        else
            echo "ERROR: Cannot generate SNP BED - no snpdata found"
            exit 1
        fi
    fi

    # Step 5: Prefilter BAM by removing reads overlapping SNPs
    echo "[5/7] Prefiltering BAM to remove SNP-overlapping reads..."
    PREFILTERED_BAM="$WORK_DIR/fixture_human_prefiltered.bam"
    
    # Use samtools view -L -U to split BAM into overlapping/non-overlapping
    # -L keeps reads overlapping BED regions, -U outputs reads NOT overlapping
    samtools view -b -h -U "$PREFILTERED_BAM" -L "$NEW_SNPS_BED" "$BAM_FILE" > /dev/null
    
    if [[ ! -f "$PREFILTERED_BAM" ]]; then
        echo "ERROR: Prefiltering failed"
        exit 1
    fi
    
    # Index prefiltered BAM
    samtools index "$PREFILTERED_BAM"
    
    ORIG_READS=$(samtools view -c "$BAM_FILE")
    FILT_READS=$(samtools view -c "$PREFILTERED_BAM")
    REMOVED=$((ORIG_READS - FILT_READS))
    echo "✓ Prefiltered BAM created: $FILT_READS reads ($REMOVED removed for SNP overlap)"

    # Step 6: Run GRAND-SLAM quantification on original BAM (not prefiltered)
    # GEDI handles SNP masking internally via snpdata
    # Use -full to include Conversions/Coverage columns
    # CRITICAL: Use -err 0.001 to match original fixture (without this, correlation drops to ~0.86)
    echo "[6/7] Running GRAND-SLAM quantification..."
    
    # Clean up previous GEDI outputs (GEDI skips if output exists)
    rm -f "$GEDI_PREFIX".* 2>/dev/null || true
    
    "$GEDI_BIN" -e Slam \
        -reads "$BAM_FILE" \
        -genomic "$GEDI_GENOME" \
        -prefix "$GEDI_PREFIX" \
        -strandness Sense \
        -nthreads 4 \
        -full \
        -err 0.001 \
        -progress \
        > "$WORK_DIR/gedi_quant.log" 2>&1 || {
        echo "WARNING: GEDI quantification had errors, check log"
        tail -20 "$WORK_DIR/gedi_quant.log"
    }
    
    NEW_REF_TSV="$GEDI_PREFIX.tsv.gz"
    if [[ ! -f "$NEW_REF_TSV" ]]; then
        echo "ERROR: GRAND-SLAM quantification failed - output not found"
        exit 1
    fi
    echo "✓ GRAND-SLAM quantification complete: $NEW_REF_TSV"
fi

# Step 7: Verify output matches existing fixture
echo "[7/7] Verifying regenerated outputs..."

# Compare SNP BED
if [[ -f "$NEW_SNPS_BED" ]] && [[ -f "$SNPS_BED" ]]; then
    ORIG_SNP_COUNT=$(wc -l < "$SNPS_BED")
    NEW_SNP_COUNT=$(wc -l < "$NEW_SNPS_BED")
    echo ""
    echo "SNP BED comparison:"
    echo "  Original: $ORIG_SNP_COUNT entries"
    echo "  Regenerated: $NEW_SNP_COUNT entries"
    
    # Compare sorted versions (order may differ due to snpdata ordering)
    if diff -q <(sort "$SNPS_BED") <(sort "$NEW_SNPS_BED") > /dev/null 2>&1; then
        echo "  Status: ✅ MATCH (same content)"
    elif [[ "$ORIG_SNP_COUNT" -eq "$NEW_SNP_COUNT" ]]; then
        COMMON=$(comm -12 <(sort "$SNPS_BED") <(sort "$NEW_SNPS_BED") | wc -l)
        echo "  Status: ⚠️ DIFFERENT (same count, $COMMON in common)"
    else
        echo "  Status: ⚠️ DIFFERENT (count mismatch)"
    fi
fi

# Compare reference TSV
if [[ -f "$NEW_REF_TSV" ]] && [[ -f "$REF_TSV" ]]; then
    echo ""
    echo "Reference TSV comparison:"
    
    # Get checksums
    ORIG_SHA=$(sha256sum "$REF_TSV" | cut -d' ' -f1)
    NEW_SHA=$(sha256sum "$NEW_REF_TSV" | cut -d' ' -f1)
    
    echo "  Original SHA256: $ORIG_SHA"
    echo "  Regenerated SHA256: $NEW_SHA"
    
    if [[ "$ORIG_SHA" == "$NEW_SHA" ]]; then
        echo "  Status: ✅ EXACT MATCH (identical checksums)"
    else
        echo "  Status: ⚠️ DIFFERENT (checksums differ)"
        
        # Compare gene counts
        ORIG_GENES=$(zcat "$REF_TSV" | tail -n +2 | wc -l)
        NEW_GENES=$(zcat "$NEW_REF_TSV" | tail -n +2 | wc -l)
        echo "  Original genes: $ORIG_GENES"
        echo "  Regenerated genes: $NEW_GENES"
        
        # Calculate correlation of MAP values if possible
        if command -v python3 &>/dev/null; then
            echo ""
            echo "  Computing NTR correlation..."
            python3 - "$REF_TSV" "$NEW_REF_TSV" << 'PY'
import gzip
import numpy as np
import sys

def load_map(path):
    data = {}
    with gzip.open(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        # Find MAP column
        map_col = None
        for i, h in enumerate(header):
            if 'MAP' in h:
                map_col = i
                break
        if map_col is None:
            return data
        for line in f:
            parts = line.strip().split('\t')
            gene = parts[0]
            try:
                data[gene] = float(parts[map_col])
            except:
                pass
    return data

orig = load_map(sys.argv[1])
regen = load_map(sys.argv[2])

# Common genes
common = set(orig.keys()) & set(regen.keys())
if len(common) < 10:
    print(f"  Too few common genes ({len(common)}) to compute correlation")
    sys.exit(0)

orig_vals = np.array([orig[g] for g in common])
regen_vals = np.array([regen[g] for g in common])

# Filter out zeros for correlation
mask = (orig_vals > 0) | (regen_vals > 0)
if mask.sum() < 10:
    print(f"  Too few non-zero values for correlation")
    sys.exit(0)

corr = np.corrcoef(orig_vals[mask], regen_vals[mask])[0, 1]
print(f"  Pearson correlation (MAP): {corr:.4f}")
print(f"  Common genes: {len(common)}")
print(f"  Non-zero genes compared: {mask.sum()}")
print(f"  Mean absolute difference: {np.mean(np.abs(orig_vals - regen_vals)):.4f}")
PY
        fi
    fi
fi

# Summary
echo ""
echo "========================================================================"
echo "Regeneration Summary"
echo "========================================================================"
echo "Work directory: $WORK_DIR"
echo ""
echo "Generated files:"
[[ -f "$WORK_DIR/slam_100000_reads_${SRA_ACCESSION}.fastq.gz" ]] && \
    echo "  ✓ FASTQ: $WORK_DIR/slam_100000_reads_${SRA_ACCESSION}.fastq.gz"
[[ -f "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" ]] && \
    echo "  ✓ BAM: ${STAR_PREFIX}Aligned.sortedByCoord.out.bam"
[[ -f "$NEW_SNPS_BED" ]] && \
    echo "  ✓ SNP BED: $NEW_SNPS_BED ($(wc -l < "$NEW_SNPS_BED") entries)"
[[ -f "$PREFILTERED_BAM" ]] && \
    echo "  ✓ Prefiltered BAM: $PREFILTERED_BAM"
[[ -f "$NEW_REF_TSV" ]] && \
    echo "  ✓ Reference TSV: $NEW_REF_TSV"
echo ""
echo "To update fixtures with regenerated files:"
echo "  cp $NEW_SNPS_BED $SNPS_BED"
echo "  cp $NEW_REF_TSV $REF_TSV"
echo ""
echo "To verify parity with STAR-Slam:"
echo "  bash $SCRIPT_DIR/run_slam_fixture_parity.sh"
echo ""
echo "========================================================================"
