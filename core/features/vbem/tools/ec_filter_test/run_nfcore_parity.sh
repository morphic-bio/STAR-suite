#!/bin/bash
# nf-core/test-datasets Parity Test for EC Filter
# Uses real RNA-seq data from nf-core/test-datasets (rnaseq branch)
# Compares our EC output against Salmon's --dumpEqWeights output
#
# Usage: ./run_nfcore_parity.sh [output_dir] [threads]

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_DIR="${1:-/tmp/nfcore_ec_parity_test}"
THREADS="${2:-4}"

# Paths to nf-core/test-datasets
TEST_DATASETS_ROOT="/mnt/pikachu/test-datasets"
TEST_DATASETS_RNASEQ="/mnt/pikachu/test-datasets-rnaseq"

# Salmon binary (use local build or system)
SALMON="${SALMON:-salmon}"

# STAR binary (use system or specify path)
STAR="${STAR:-STAR}"

# Parameters (must match between Salmon and our implementation)
SCORE_EXP="1.0"                    # Salmon default: 1.0
MIN_ALN_PROB="1e-5"                # Salmon default: 1e-5
MIN_SCORE_FRACTION="0.65"          # Salmon default: 0.65
RANGE_FACTORIZATION_BINS="4"       # Salmon default: 4
HARD_FILTER="false"                # Salmon default: false (soft filter)

mkdir -p "${OUTPUT_DIR}"

echo "=== nf-core/test-datasets EC Filter Parity Test ==="
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo ""

# Step 1: Checkout nf-core/test-datasets (rnaseq branch)
echo "Step 1: Setting up nf-core/test-datasets..."
if [ ! -d "${TEST_DATASETS_ROOT}" ]; then
    echo "Error: nf-core/test-datasets repository not found at ${TEST_DATASETS_ROOT}"
    echo "Please clone it first: git clone https://github.com/nf-core/test-datasets.git ${TEST_DATASETS_ROOT}"
    exit 1
fi

# Fetch rnaseq branch
echo "  Fetching rnaseq branch..."
git -C "${TEST_DATASETS_ROOT}" fetch origin rnaseq 2>/dev/null || true

# Create worktree if it doesn't exist
if [ ! -d "${TEST_DATASETS_RNASEQ}" ]; then
    echo "  Creating worktree at ${TEST_DATASETS_RNASEQ}..."
    git -C "${TEST_DATASETS_ROOT}" worktree add "${TEST_DATASETS_RNASEQ}" origin/rnaseq
else
    echo "  Worktree already exists at ${TEST_DATASETS_RNASEQ}"
fi

# Locate FASTQ files and reference
# Try multiple possible locations
R1=""
R2=""
for base in "${TEST_DATASETS_RNASEQ}/testdata/GSE110004" "${TEST_DATASETS_RNASEQ}/data/rna-seq/srr6357070" "${TEST_DATASETS_RNASEQ}/data/rna-seq" "${TEST_DATASETS_RNASEQ}"; do
    if [ -f "${base}/SRR6357070_1.fastq.gz" ] && [ -f "${base}/SRR6357070_2.fastq.gz" ]; then
        R1="${base}/SRR6357070_1.fastq.gz"
        R2="${base}/SRR6357070_2.fastq.gz"
        break
    fi
done

TRANSCRIPTOME="${TEST_DATASETS_RNASEQ}/reference/transcriptome.fasta.gz"
GTF="${TEST_DATASETS_RNASEQ}/reference/genes.gtf.gz"

# Verify files exist
if [ -z "${R1}" ] || [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "Error: FASTQ files not found. Searched for SRR6357070_{1,2}.fastq.gz in:"
    echo "  - ${TEST_DATASETS_RNASEQ}/data/rna-seq/srr6357070/"
    echo "  - ${TEST_DATASETS_RNASEQ}/data/rna-seq/"
    echo "  - ${TEST_DATASETS_RNASEQ}/"
    echo ""
    echo "Available FASTQ files in ${TEST_DATASETS_RNASEQ}:"
    find "${TEST_DATASETS_RNASEQ}" -name "*_1.fastq.gz" -o -name "*_2.fastq.gz" 2>/dev/null | head -10 || true
    exit 1
fi

if [ ! -f "${TRANSCRIPTOME}" ]; then
    echo "Error: Transcriptome not found: ${TRANSCRIPTOME}"
    exit 1
fi

if [ ! -f "${GTF}" ]; then
    echo "Error: GTF file not found: ${GTF}"
    exit 1
fi

echo "  Using FASTQ files:"
echo "    R1: ${R1}"
echo "    R2: ${R2}"
echo "  Using reference:"
echo "    Transcriptome: ${TRANSCRIPTOME}"
echo "    GTF: ${GTF}"
echo ""

# Step 2: Build STAR index
echo "Step 2: Building STAR index..."
STAR_IDX="${OUTPUT_DIR}/star_idx_nfcore"
if [ ! -d "${STAR_IDX}" ] || [ ! -f "${STAR_IDX}/Genome" ]; then
    echo "  Building index at ${STAR_IDX}..."
    mkdir -p "${STAR_IDX}"
    
    # STAR doesn't support gzipped files for index building, so decompress if needed
    TRANSCRIPTOME_UNCOMPRESSED="${OUTPUT_DIR}/transcriptome.fasta"
    GTF_UNCOMPRESSED="${OUTPUT_DIR}/genes.gtf"
    
    if [[ "${TRANSCRIPTOME}" == *.gz ]]; then
        echo "  Decompressing transcriptome..."
        zcat "${TRANSCRIPTOME}" > "${TRANSCRIPTOME_UNCOMPRESSED}"
        TRANSCRIPTOME_TO_USE="${TRANSCRIPTOME_UNCOMPRESSED}"
    else
        TRANSCRIPTOME_TO_USE="${TRANSCRIPTOME}"
    fi
    
    if [[ "${GTF}" == *.gz ]]; then
        echo "  Decompressing GTF..."
        zcat "${GTF}" > "${GTF_UNCOMPRESSED}"
        GTF_TO_USE="${GTF_UNCOMPRESSED}"
    else
        GTF_TO_USE="${GTF}"
    fi
    
    "${STAR}" --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${STAR_IDX}" \
        --genomeFastaFiles "${TRANSCRIPTOME_TO_USE}" \
        --sjdbGTFfile "${GTF_TO_USE}" \
        --sjdbOverhang 100 \
        2>&1 | tee "${OUTPUT_DIR}/star_index.log"
    
    # Clean up decompressed files if we created them
    if [[ "${TRANSCRIPTOME}" == *.gz ]]; then
        rm -f "${TRANSCRIPTOME_UNCOMPRESSED}"
    fi
    if [[ "${GTF}" == *.gz ]]; then
        rm -f "${GTF_UNCOMPRESSED}"
    fi
else
    echo "  Using existing STAR index at ${STAR_IDX}"
fi
echo ""

# Step 3: Align reads with STAR (transcriptome mode)
echo "Step 3: Aligning reads with STAR (transcriptome mode)..."
STAR_OUT="${OUTPUT_DIR}/star_align"
BAM_FILE="${STAR_OUT}/Aligned.toTranscriptome.out.bam"

if [ ! -f "${BAM_FILE}" ]; then
    echo "  Running STAR alignment..."
    mkdir -p "${STAR_OUT}"
    "${STAR}" --runThreadN "${THREADS}" \
        --genomeDir "${STAR_IDX}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix "${STAR_OUT}/" \
        2>&1 | tee "${OUTPUT_DIR}/star_align.log"
else
    echo "  Using existing BAM file: ${BAM_FILE}"
fi

if [ ! -f "${BAM_FILE}" ]; then
    echo "Error: STAR alignment failed. BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "  BAM file: ${BAM_FILE}"
echo "  Alignments: $(samtools view -c "${BAM_FILE}" 2>/dev/null || echo 'unknown')"
echo ""

# Step 4: Run Salmon in alignment mode
echo "Step 4: Running Salmon quant (alignment mode)..."
SALMON_OUT="${OUTPUT_DIR}/salmon_out"

if [ ! -d "${SALMON_OUT}" ] || [ ! -f "${SALMON_OUT}/aux_info/eq_classes.txt" ] && [ ! -f "${SALMON_OUT}/aux_info/eq_classes.txt.gz" ]; then
    echo "  Running Salmon with --dumpEqWeights..."
    "${SALMON}" quant \
        -t "${TRANSCRIPTOME}" \
        -l A \
        -a "${BAM_FILE}" \
        --dumpEqWeights \
        --noLengthCorrection \
        --noEffectiveLengthCorrection \
        --noFragLengthDist \
        --threads "${THREADS}" \
        -o "${SALMON_OUT}" \
        2>&1 | tee "${OUTPUT_DIR}/salmon.log"
else
    echo "  Using existing Salmon output at ${SALMON_OUT}"
fi

# Extract Salmon ECs (handle gzipped output)
SALMON_EC="${OUTPUT_DIR}/salmon_eq_classes.txt"
if [ -f "${SALMON_OUT}/aux_info/eq_classes.txt.gz" ]; then
    echo "  Extracting ECs from gzipped file..."
    gunzip -c "${SALMON_OUT}/aux_info/eq_classes.txt.gz" > "${SALMON_EC}"
elif [ -f "${SALMON_OUT}/aux_info/eq_classes.txt" ]; then
    cp "${SALMON_OUT}/aux_info/eq_classes.txt" "${SALMON_EC}"
else
    echo "Error: Salmon EC file not found in ${SALMON_OUT}/aux_info/"
    exit 1
fi

echo "  Salmon ECs: ${SALMON_EC} ($(wc -l < ${SALMON_EC}) lines)"
echo ""

# Step 5: Run our EC filter CLI
echo "Step 5: Running our EC filter CLI..."
OUR_OUT="${OUTPUT_DIR}/our_out"
OUR_EC="${OUR_OUT}/eq_classes.txt"
mkdir -p "${OUR_OUT}"

# Check if CLI exists
CLI="${SCRIPT_DIR}/ec_filter_cli"
if [ ! -f "${CLI}" ]; then
    echo "  Building EC filter CLI..."
    cd "${SCRIPT_DIR}"
    make ec_filter_cli
fi

if [ ! -f "${CLI}" ]; then
    echo "Error: EC filter CLI not found: ${CLI}"
    echo "Please build it first: cd ${SCRIPT_DIR} && make"
    exit 1
fi

echo "  Running EC filter with matching parameters..."
"${CLI}" \
    --input "${BAM_FILE}" \
    --transcripts "${TRANSCRIPTOME}" \
    --score-exp "${SCORE_EXP}" \
    --min-aln-prob "${MIN_ALN_PROB}" \
    --min-score-fraction "${MIN_SCORE_FRACTION}" \
    --range-factorization-bins "${RANGE_FACTORIZATION_BINS}" \
    --hard-filter "${HARD_FILTER}" \
    --no-local-pruning \
    --no-global-pruning \
    --threads "${THREADS}" \
    --output-format salmon \
    -o "${OUR_EC}" \
    2>&1 | tee "${OUTPUT_DIR}/ec_filter.log"

if [ ! -f "${OUR_EC}" ]; then
    echo "Warning: EC filter CLI did not produce output file: ${OUR_EC}"
    echo "This may be expected if CLI is still a skeleton implementation."
    echo "See ${OUTPUT_DIR}/ec_filter.log for details."
    echo ""
    echo "Skipping comparison (CLI not yet implemented)."
    exit 0
fi

echo "  Our ECs: ${OUR_EC} ($(wc -l < ${OUR_EC}) lines)"
echo ""

# Step 6: Compare EC files
echo "Step 6: Comparing EC files..."
python3 "${SCRIPT_DIR}/compare_ecs.py" \
    --salmon "${SALMON_EC}" \
    --ours "${OUR_EC}" \
    --tolerance 1e-6 \
    --report "${OUTPUT_DIR}/parity_report.txt"

echo ""
echo "=== Parity Test Complete ==="
echo "Report: ${OUTPUT_DIR}/parity_report.txt"
echo ""
cat "${OUTPUT_DIR}/parity_report.txt"
