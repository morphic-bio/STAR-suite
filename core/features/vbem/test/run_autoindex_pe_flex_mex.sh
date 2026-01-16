#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR_BIN="${STAR_BIN:-$SCRIPT_DIR/../source/STAR}"

CELLRANGER_REF_RELEASE="${CELLRANGER_REF_RELEASE:-2024-A}"

normalize_release() {
    local release="$1"
    local normalized alias
    normalized="$(echo "$release" | tr '[:upper:]' '[:lower:]')"
    alias="${normalized//-/}"
    alias="${alias// /}"
    case "$normalized" in
        2024-a) echo "2024-A" ;;
        2020-a) echo "2020-A" ;;
        *)
            case "$alias" in
                2024a|2024) echo "2024-A" ;;
                2020a|2020) echo "2020-A" ;;
                *) echo "$release" ;;
            esac
            ;;
    esac
}

CELLRANGER_REF_RELEASE="$(normalize_release "$CELLRANGER_REF_RELEASE")"

case "$CELLRANGER_REF_RELEASE" in
    2024-A)
        DEFAULT_INDEX_ROOT="/storage/autoindex_110_44"
        DEFAULT_FLEX_PROBE_SET="/storage/genome-110-44/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"
        ;;
    2020-A)
        DEFAULT_INDEX_ROOT="/storage/autoindex_98_32"
        DEFAULT_FLEX_PROBE_SET="/storage/genome-98-32/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
        ;;
    *)
        echo "ERROR: Unsupported CELLRANGER_REF_RELEASE: $CELLRANGER_REF_RELEASE"
        echo "Supported: 2024-A (aliases: 2024A, 2024), 2020-A (aliases: 2020A, 2020)"
        exit 1
        ;;
esac

INDEX_ROOT="${INDEX_ROOT:-$DEFAULT_INDEX_ROOT}"
INDEX_PE="$INDEX_ROOT/pe_index"
INDEX_FLEX="$INDEX_ROOT/flex_index"
SHARED_CACHE_DIR="${SHARED_CACHE_DIR:-$INDEX_ROOT/cellranger_ref_cache}"

OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/autoindex_pe_flex_mex_output}"
TMP_DIR="${TMP_DIR:-/storage/100K/tmp/autoindex_pe_flex_mex}"

BASELINE_PE_DIR="${BASELINE_PE_DIR:-$SCRIPT_DIR/../tests/baseline_test_output}"
BASELINE_FLEX_DIR="${BASELINE_FLEX_DIR:-$SCRIPT_DIR/../tests/gold_standard}"

FLEX_PROBE_SET="${FLEX_PROBE_SET:-$DEFAULT_FLEX_PROBE_SET}"

CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"

READS_R2="${READS_R2:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R2_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R2_001.fastq.gz}"
READS_R1="${READS_R1:-/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L002_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L003_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L004_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L005_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L006_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L007_R1_001.fastq.gz,/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L008_R1_001.fastq.gz}"

SJDB_OVERHANG="${SJDB_OVERHANG:-}"
RUN_THREADS="${RUN_THREADS:-8}"
CLEAN_INDEX="${CLEAN_INDEX:-0}"
AUTO_CKSUM_UPDATE="${AUTO_CKSUM_UPDATE:-yes}"
REMOVE_DEPRECATED="${REMOVE_DEPRECATED:-Yes}"
LIMIT_GENOME_GENERATE_RAM="${LIMIT_GENOME_GENERATE_RAM:-100000000000}"
FORCE_BUILD="${FORCE_BUILD:-1}"
FORCE_ALL_INDEX="${FORCE_ALL_INDEX:-0}"

if [ ! -x "$STAR_BIN" ]; then
    echo "ERROR: STAR binary not found: $STAR_BIN"
    exit 1
fi

for path in "$FLEX_PROBE_SET" "$CB_WHITELIST" "$SAMPLE_WHITELIST" "$SAMPLE_PROBES"; do
    if [ ! -f "$path" ]; then
        echo "ERROR: Missing required file: $path"
        exit 1
    fi
done

mkdir -p "$OUT_DIR" "$TMP_DIR"

if [ "$CLEAN_INDEX" = "1" ]; then
    rm -rf "$INDEX_ROOT"
fi

run_star() {
    local label="$1"
    local log_file="$2"
    shift 2
    echo ">>> $label"
    "$@" > "$log_file" 2>&1 || {
        echo "FAIL: $label"
        echo "Log: $log_file"
        tail -n 200 "$log_file"
        exit 1
    }
}

compare_mex_dir() {
    local baseline="$1"
    local candidate="$2"
    local label="$3"
    for f in matrix.mtx barcodes.tsv features.tsv; do
        if [ ! -f "$baseline/$f" ]; then
            echo "FAIL: Missing baseline $baseline/$f"
            exit 1
        fi
        if [ ! -f "$candidate/$f" ]; then
            echo "FAIL: Missing candidate $candidate/$f"
            exit 1
        fi
        if ! cmp -s "$baseline/$f" "$candidate/$f"; then
            echo "FAIL: $label differs in $f"
            diff -u "$baseline/$f" "$candidate/$f" || true
            exit 1
        fi
    done
    echo "PASS: $label MEX matches"
}

build_index() {
    local label="$1"
    local genome_dir="$2"
    local extra_args="$3"
    local sjdb_args=()
    local ram_args=()
    local force_args=()

    if [ "$FORCE_BUILD" = "0" ] && [ -f "$genome_dir/Genome" ] && [ -f "$genome_dir/SA" ] && [ -f "$genome_dir/SAindex" ]; then
        echo "Skipping $label index build (already exists): $genome_dir"
        return
    fi

    mkdir -p "$genome_dir"
    if [ -n "$SJDB_OVERHANG" ]; then
        sjdb_args+=(--sjdbOverhang "$SJDB_OVERHANG")
    fi
    if [ -n "$LIMIT_GENOME_GENERATE_RAM" ]; then
        ram_args+=(--limitGenomeGenerateRAM "$LIMIT_GENOME_GENERATE_RAM")
    fi
    if [ "$FORCE_BUILD" = "1" ]; then
        force_args+=(--forceIndex Yes)
        if [ "$FORCE_ALL_INDEX" = "1" ]; then
            force_args+=(--forceAllIndex Yes)
        fi
    fi

    run_star "$label index build" "$genome_dir/build.log" \
        "$STAR_BIN" --runMode genomeGenerate \
        --autoIndex Yes \
        --cellrangerStyleIndex Yes \
        --cellrangerRefRelease "$CELLRANGER_REF_RELEASE" \
        --autoCksumUpdate "$AUTO_CKSUM_UPDATE" \
        --removeDeprecated "$REMOVE_DEPRECATED" \
        --cellrangerStyleCacheDir "$SHARED_CACHE_DIR" \
        --genomeDir "$genome_dir" \
        --runThreadN "$RUN_THREADS" \
        "${sjdb_args[@]}" \
        "${ram_args[@]}" \
        "${force_args[@]}" \
        $extra_args
}

echo "=============================================="
echo "AutoIndex default (${CELLRANGER_REF_RELEASE})"
echo "=============================================="
echo "STAR binary: $STAR_BIN"
echo "CellRanger release: $CELLRANGER_REF_RELEASE"
echo "Index root: $INDEX_ROOT"
echo "Output dir: $OUT_DIR"
echo "Baseline PE: $BASELINE_PE_DIR"
echo "Baseline Flex: $BASELINE_FLEX_DIR"
echo "Flex probe set: $FLEX_PROBE_SET"
echo "Remove deprecated: $REMOVE_DEPRECATED"
echo "Reads R2: $READS_R2"
echo "Reads R1: $READS_R1"
echo ""

build_index "PE" "$INDEX_PE" ""
build_index "Flex" "$INDEX_FLEX" "--flexGeneProbeSet $FLEX_PROBE_SET"

rm -rf "$OUT_DIR/pe_autoindex" "$OUT_DIR/flex_autoindex"
mkdir -p "$OUT_DIR/pe_autoindex" "$OUT_DIR/flex_autoindex"

# PE run (flex disabled)
run_star "PE autoIndex" "$OUT_DIR/pe_autoindex/run.log" \
    "$STAR_BIN" \
    --runThreadN "$RUN_THREADS" \
    --outTmpDir "$TMP_DIR/pe" \
    --genomeDir "$INDEX_PE" \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$CB_WHITELIST" \
    --limitIObufferSize 50000000 50000000 \
    --outSJtype None \
    --outBAMcompression 6 \
    --alignIntronMax 500000 \
    --outFilterMismatchNmax 6 \
    --outFilterMismatchNoverReadLmax 1.0 \
    --outFilterMatchNmin 25 \
    --outSAMunmapped None \
    --outFilterMatchNminOverLread 0 \
    --outFilterMultimapNmax 10000 \
    --outFilterMultimapScoreRange 4 \
    --outSAMmultNmax 10000 \
    --winAnchorMultimapNmax 200 \
    --outSAMprimaryFlag AllBestScore \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0 \
    --outSAMattributes NH HI AS nM NM GX GN \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter None \
    --clipAdapterType CellRanger4 \
    --soloFeatures Gene \
    --alignEndsType Local \
    --soloStrand Unstranded \
    --chimSegmentMin 1000000 \
    --outSAMtype BAM Unsorted \
    --readFilesCommand zcat \
    --readFilesIn "$READS_R2" "$READS_R1" \
    --outFileNamePrefix "$OUT_DIR/pe_autoindex/"

# Flex run (inline + per-sample)
run_star "Flex autoIndex" "$OUT_DIR/flex_autoindex/run.log" \
    "$STAR_BIN" \
    --runThreadN "$RUN_THREADS" \
    --outTmpDir "$TMP_DIR/flex" \
    --genomeDir "$INDEX_FLEX" \
    --soloType CB_UMI_Simple \
    --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "$CB_WHITELIST" \
    --flex yes \
    --soloRemoveDeprecated "$REMOVE_DEPRECATED" \
    --soloFlexExpectedCellsPerTag 3000 \
    --soloSampleWhitelist "$SAMPLE_WHITELIST" \
    --soloSampleProbes "$SAMPLE_PROBES" \
    --soloSampleProbeOffset 68 \
    --soloFlexAllowedTags "$SAMPLE_WHITELIST" \
    --soloFlexOutputPrefix "$OUT_DIR/flex_autoindex/per_sample" \
    --limitIObufferSize 50000000 50000000 \
    --outSJtype None \
    --outBAMcompression 6 \
    --soloMultiMappers Rescue \
    --alignIntronMax 500000 \
    --outFilterMismatchNmax 6 \
    --outFilterMismatchNoverReadLmax 1.0 \
    --outFilterMatchNmin 25 \
    --outSAMunmapped None \
    --outFilterMatchNminOverLread 0 \
    --outFilterMultimapNmax 10000 \
    --outFilterMultimapScoreRange 4 \
    --outSAMmultNmax 10000 \
    --winAnchorMultimapNmax 200 \
    --outSAMprimaryFlag AllBestScore \
    --outFilterScoreMin 0 \
    --outFilterScoreMinOverLread 0 \
    --outSAMattributes NH HI AS nM NM GX GN ZG \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter None \
    --clipAdapterType CellRanger4 \
    --soloFeatures Gene \
    --alignEndsType Local \
    --soloAddTagsToUnsorted no \
    --soloStrand Unstranded \
    --chimSegmentMin 1000000 \
    --soloKeysCompat cr \
    --outSAMtype BAM Unsorted \
    --soloSampleSearchNearby no \
    --readFilesCommand zcat \
    --readFilesIn "$READS_R2" "$READS_R1" \
    --outFileNamePrefix "$OUT_DIR/flex_autoindex/"

# Compare PE raw MEX
PE_BASE_RAW="$BASELINE_PE_DIR/Solo.out/Gene/raw"
PE_CUR_RAW="$OUT_DIR/pe_autoindex/Solo.out/Gene/raw"
compare_mex_dir "$PE_BASE_RAW" "$PE_CUR_RAW" "PE"

# Compare Flex raw MEX against gold_standard/raw
FLEX_BASE_RAW="$BASELINE_FLEX_DIR/raw"
FLEX_CUR_RAW="$OUT_DIR/flex_autoindex/Solo.out/Gene/raw"
compare_mex_dir "$FLEX_BASE_RAW" "$FLEX_CUR_RAW" "Flex raw"

# Compare per-sample MEX
for SAMPLE in BC004 BC006 BC007 BC008; do
    BASE_SAMPLE="$BASELINE_FLEX_DIR/per_sample/$SAMPLE/Gene/filtered"
    CUR_SAMPLE="$OUT_DIR/flex_autoindex/per_sample/$SAMPLE/Gene/filtered"
    compare_mex_dir "$BASE_SAMPLE" "$CUR_SAMPLE" "Flex per-sample $SAMPLE"
done

echo ""
echo "=============================================="
echo "PASS: AutoIndex PE + Flex MEX parity"
echo "=============================================="
