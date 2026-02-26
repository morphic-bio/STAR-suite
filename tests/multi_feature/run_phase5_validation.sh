#!/usr/bin/env bash
# Phase 5: Multi-feature validation suite
# Runs STAR with a downsampled multi-feature fixture (10K reads, 500 LARRY features)
# and validates: E2E smoke, per-library provenance, merged MEX, chemistry normalization,
# scheduler logging, and edge cases.
#
# Usage:
#   bash run_phase5_validation.sh [fixture_dir] [output_dir]
#
# Prerequisites:
#   - STAR binary built (make -C core/legacy/source -j8 STAR)
#   - Downsampled fixture created (create_fixture_downsampled.sh)
#   - Genome index at /storage/autoindex_110_44/.../star
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR="$REPO_ROOT/core/legacy/source/STAR"
GENOME=/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star
WHITELIST=/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt
FIXTURE=${1:-/tmp/msk_multi_fixture_ds}
OUT_DIR=${2:-/tmp/msk_phase5_validation_$(date +%Y%m%d_%H%M%S)}

if [ ! -x "$STAR" ]; then echo "FATAL: STAR binary not found at $STAR"; exit 1; fi
if [ ! -d "$FIXTURE/mRNA" ]; then echo "FATAL: fixture not found — run create_fixture_downsampled.sh first"; exit 1; fi
if [ ! -d "$GENOME" ]; then echo "FATAL: genome index not found at $GENOME"; exit 1; fi

echo "=== Phase 5: Multi-Feature Validation Suite ==="
echo "  STAR:     $STAR"
echo "  Genome:   $GENOME"
echo "  Fixture:  $FIXTURE"
echo "  Output:   $OUT_DIR"
echo ""

mkdir -p "$OUT_DIR"

GRNA_REF="$FIXTURE/ref_feature_geneBC.csv"
LARRY_REF="$FIXTURE/ref_feature_larryBC.csv"

PASS=0
FAIL=0

pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

# =========================================================================
# Part A: Baseline multi-library run (mRNA + PolyIII + LARRY, 8 threads)
# =========================================================================
echo "--- Part A: Baseline multi-library run ---"

RUN_A="$OUT_DIR/run_a_baseline"
mkdir -p "$RUN_A"

cat > "$RUN_A/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$GRNA_REF,grna_de
$FIXTURE/LARRY,DE_30KO,Custom,Custom,TRU,$LARRY_REF,larry_de
EOF

GEX_R1=$(ls "$FIXTURE"/mRNA/*_R1_*.fastq.gz | head -n1)
GEX_R2=$(ls "$FIXTURE"/mRNA/*_R2_*.fastq.gz | head -n1)

PHASE5_PERF_TOL_PCT=${PHASE5_PERF_TOL_PCT:-10}

echo "  Running STAR (8 threads, no dynamic permits)..."
START_A=$(date +%s%N)

"$STAR" \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$RUN_A/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$RUN_A/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    > "$RUN_A/star_stdout.log" 2>&1

END_A=$(date +%s%N)
ELAPSED_A=$(( (END_A - START_A) / 1000000 ))
STAR_EXIT=$?

echo "  STAR exit: $STAR_EXIT (${ELAPSED_A}ms)"

if [ "$STAR_EXIT" -ne 0 ]; then
    echo "FATAL: STAR baseline run failed"
    tail -20 "$RUN_A/Log.out" 2>/dev/null
    exit 1
fi

# --- A1: Output directories exist ---
echo ""
echo "--- A1: Output directory structure ---"
if [ -d "$RUN_A/cr_assign/CRISPR_Guide_Capture/grna_de" ]; then
    pass "grna_de assign output directory exists"
else
    fail "grna_de assign output directory missing"
fi

if [ -d "$RUN_A/cr_assign/Custom/larry_de" ]; then
    pass "larry_de assign output directory exists"
else
    fail "larry_de assign output directory missing"
fi

# --- A2: Per-library provenance manifests ---
echo ""
echo "--- A2: Per-library provenance manifests ---"

# Provenance manifests are written inside the assign output subdirectory
# (e.g., .../grna_de/PolyIII/pf_library_provenance.tsv)
GRNA_PROV=$(find "$RUN_A/cr_assign" -path "*/grna_de/*/pf_library_provenance.tsv" 2>/dev/null | head -n1)
LARRY_PROV=$(find "$RUN_A/cr_assign" -path "*/larry_de/*/pf_library_provenance.tsv" 2>/dev/null | head -n1)

if [ -n "$GRNA_PROV" ] && [ -f "$GRNA_PROV" ]; then
    pass "grna_de provenance manifest exists"
    for KEY in library_id feature_type feature_ref status return_code; do
        if grep -q "^${KEY}	" "$GRNA_PROV"; then
            pass "grna_de provenance has $KEY"
        else
            fail "grna_de provenance missing $KEY"
        fi
    done
    if grep -q "^status	OK" "$GRNA_PROV"; then
        pass "grna_de provenance status=OK"
    else
        fail "grna_de provenance status != OK"
    fi
else
    fail "grna_de provenance manifest missing"
fi

if [ -n "$LARRY_PROV" ] && [ -f "$LARRY_PROV" ]; then
    pass "larry_de provenance manifest exists"
    if grep -q "^status	OK" "$LARRY_PROV"; then
        pass "larry_de provenance status=OK"
    else
        fail "larry_de provenance status != OK"
    fi
else
    fail "larry_de provenance manifest missing"
fi

# --- A3: MEX stub outputs (features.tsv, barcodes.tsv) ---
echo ""
echo "--- A3: MEX stub outputs ---"

for LIB_DIR in "$RUN_A/cr_assign/CRISPR_Guide_Capture/grna_de" "$RUN_A/cr_assign/Custom/larry_de"; do
    LIB_NAME=$(basename "$LIB_DIR")
    # Look for features.tsv in any subdirectory (library type subdir)
    FEAT_FILES=$(find "$LIB_DIR" -name "features.tsv" -o -name "features.tsv.gz" 2>/dev/null | head -n1)
    if [ -n "$FEAT_FILES" ]; then
        pass "$LIB_NAME has features.tsv"
    else
        fail "$LIB_NAME missing features.tsv"
    fi
    BARC_FILES=$(find "$LIB_DIR" -name "barcodes.tsv" -o -name "barcodes.tsv.gz" 2>/dev/null | head -n1)
    if [ -n "$BARC_FILES" ]; then
        pass "$LIB_NAME has barcodes.tsv"
    else
        fail "$LIB_NAME missing barcodes.tsv"
    fi
done

# --- A4: Chemistry and feature ref logged correctly ---
echo ""
echo "--- A4: Log.out diagnostics ---"

if grep -q 'grna_de.*CRISPR Guide Capture' "$RUN_A/Log.out"; then
    pass "grna_de chemistry logged"
else
    fail "grna_de chemistry not logged"
fi

if grep -q 'larry_de.*Custom' "$RUN_A/Log.out"; then
    pass "larry_de chemistry logged"
else
    fail "larry_de chemistry not logged"
fi

if grep -q 'star_feature_ref=.*ref_feature_geneBC.csv' "$RUN_A/Log.out"; then
    pass "grna_de per-library feature ref logged"
else
    fail "grna_de per-library feature ref not logged"
fi

if grep -q 'star_feature_ref=.*ref_feature_larryBC.csv' "$RUN_A/Log.out"; then
    pass "larry_de per-library feature ref logged"
else
    fail "larry_de per-library feature ref not logged"
fi

# --- A5: Phase 4 scheduler log ---
echo ""
echo "--- A5: Scheduler log ---"

# Baseline run A does NOT set --dynamicThreadInterface. The Phase 4 scheduler
# must NOT activate; if it does, the gating condition is broken.
if ! grep -q 'pf-multi library scheduler (Phase 4)' "$RUN_A/Log.out"; then
    pass "Scheduler log correctly absent without --dynamicThreadInterface"
else
    fail "Scheduler log present without --dynamicThreadInterface (gating broken)"
fi

# --- A6: Provenance summary table in Log.out ---
echo ""
echo "--- A6: Provenance summary table ---"

if grep -q 'pf-multi library summary:' "$RUN_A/Log.out"; then
    pass "Provenance summary table header present"
else
    fail "Provenance summary table header missing"
fi

# =========================================================================
# Part B: Single-library backward compatibility (GEX only, no feature libs)
# =========================================================================
echo ""
echo "--- Part B: Single-library backward compatibility ---"

RUN_B="$OUT_DIR/run_b_single_lib"
mkdir -p "$RUN_B"

cat > "$RUN_B/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU
EOF

echo "  Running STAR (single-lib, no feature libraries)..."

"$STAR" \
    --runMode alignReads \
    --runThreadN 4 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$RUN_B/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$RUN_B/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    > "$RUN_B/star_stdout.log" 2>&1

if [ $? -eq 0 ]; then
    pass "Single-lib (GEX only) run completes successfully"
else
    fail "Single-lib (GEX only) run failed"
fi

# For GEX-only config, cr_assign must be absent or completely empty.
# Any entries indicate feature routing leaked into a GEX-only run.
if [ ! -d "$RUN_B/cr_assign" ]; then
    pass "No cr_assign directory for GEX-only config (backward compat)"
elif [ -z "$(find "$RUN_B/cr_assign" -mindepth 1 -maxdepth 1 2>/dev/null)" ]; then
    pass "cr_assign exists but is empty for GEX-only config (backward compat)"
else
    fail "cr_assign contains entries for GEX-only config (backward compat broken)"
fi

# =========================================================================
# Part C: runThreadN=1 single-library (no clamp to 2)
# =========================================================================
echo ""
echo "--- Part C: runThreadN=1 edge case ---"

RUN_C="$OUT_DIR/run_c_1thread"
mkdir -p "$RUN_C"

cat > "$RUN_C/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$GRNA_REF,grna_de
EOF

echo "  Running STAR (1 thread, 1 feature library)..."

"$STAR" \
    --runMode alignReads \
    --runThreadN 1 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$RUN_C/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$RUN_C/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    > "$RUN_C/star_stdout.log" 2>&1

if [ $? -eq 0 ]; then
    pass "runThreadN=1 completes without deadlock"
else
    fail "runThreadN=1 failed (possible deadlock)"
fi

PROV_C=$(find "$RUN_C/cr_assign" -name "pf_library_provenance.tsv" 2>/dev/null | head -n1)
if [ -n "$PROV_C" ] && grep -q "^status	OK" "$PROV_C"; then
    pass "runThreadN=1 produces provenance manifest with status=OK"
else
    fail "runThreadN=1 missing provenance or status != OK"
fi

# =========================================================================
# Part D: runThreadN=2, two feature libraries (min guarantee)
# =========================================================================
echo ""
echo "--- Part D: runThreadN=2, two feature libraries ---"

RUN_D="$OUT_DIR/run_d_2thread_2lib"
mkdir -p "$RUN_D"

cat > "$RUN_D/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$GRNA_REF,grna_de
$FIXTURE/LARRY,DE_30KO,Custom,Custom,TRU,$LARRY_REF,larry_de
EOF

echo "  Running STAR (2 threads, 2 feature libraries)..."

"$STAR" \
    --runMode alignReads \
    --runThreadN 2 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$RUN_D/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --outFileNamePrefix "$RUN_D/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    > "$RUN_D/star_stdout.log" 2>&1

if [ $? -eq 0 ]; then
    pass "runThreadN=2, 2 feature libs completes"
else
    fail "runThreadN=2, 2 feature libs failed"
fi

# Both provenance files should exist (under sample subdirs)
for LIB in grna_de larry_de; do
    PROV=$(find "$RUN_D/cr_assign" -path "*/$LIB/*/pf_library_provenance.tsv" 2>/dev/null | head -n1)
    if [ -n "$PROV" ] && grep -q "^status	OK" "$PROV"; then
        pass "runThreadN=2: $LIB provenance status=OK"
    else
        fail "runThreadN=2: $LIB provenance missing or not OK"
    fi
done

# =========================================================================
# Part E: Performance (deferred to after Part H for comparison)
# =========================================================================

# =========================================================================
# Part F: Feature type correctness in MEX stubs
# =========================================================================
echo ""
echo "--- Part F: Feature type correctness ---"

# Assert every row in features.tsv has exactly the expected feature type.
check_all_feature_types() {
    local feat_file="$1" expected_type="$2" lib_label="$3"
    if [ -z "$feat_file" ] || [ ! -s "$feat_file" ]; then
        fail "$lib_label features.tsv not found or empty"
        return
    fi
    local total_rows bad_rows
    total_rows=$(wc -l < "$feat_file")
    bad_rows=$(awk -F'\t' -v et="$expected_type" '$3 != et { count++ } END { print count+0 }' "$feat_file")
    if [ "$total_rows" -gt 0 ] && [ "$bad_rows" -eq 0 ]; then
        pass "$lib_label features.tsv: all $total_rows rows have type '$expected_type'"
    else
        fail "$lib_label features.tsv: $bad_rows/$total_rows rows have wrong type (expected '$expected_type')"
    fi
}

GRNA_FEAT=$(find "$RUN_A/cr_assign/CRISPR_Guide_Capture/grna_de" -name "features.tsv" 2>/dev/null | head -n1)
check_all_feature_types "$GRNA_FEAT" "CRISPR Guide Capture" "grna_de"

LARRY_FEAT=$(find "$RUN_A/cr_assign/Custom/larry_de" -name "features.tsv" 2>/dev/null | head -n1)
check_all_feature_types "$LARRY_FEAT" "Custom" "larry_de"

# =========================================================================
# Part G: No stale .tmp provenance files left behind
# =========================================================================
echo ""
echo "--- Part G: No stale .tmp artifacts ---"

TMP_FILES=$(find "$RUN_A/cr_assign" -name "*.tmp" 2>/dev/null)
if [ -z "$TMP_FILES" ]; then
    pass "No stale .tmp provenance files"
else
    fail "Stale .tmp files found: $TMP_FILES"
fi

# =========================================================================
# Part H: Dynamic permits + scheduler (multi-lib with --dynamicThreadInterface)
# =========================================================================
echo ""
echo "--- Part H: Dynamic permits + Phase 4 scheduler ---"

RUN_H="$OUT_DIR/run_h_dynamic"
mkdir -p "$RUN_H"

cat > "$RUN_H/multi_config.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$GRNA_REF,grna_de
$FIXTURE/LARRY,DE_30KO,Custom,Custom,TRU,$LARRY_REF,larry_de
EOF

echo "  Running STAR (8 threads, dynamic permits active, 2 feature libs)..."
START_H=$(date +%s%N)

"$STAR" \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --pfMultiConfig "$RUN_H/multi_config.csv" \
    --defaultCrCompat yes \
    --crChemistry auto \
    --dynamicThreadInterface 1 \
    --outFileNamePrefix "$RUN_H/" \
    --outSAMtype BAM Unsorted \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures Gene GeneFull \
    --soloBarcodeReadLength 0 \
    --readMapNumber 50000 \
    > "$RUN_H/star_stdout.log" 2>&1

END_H=$(date +%s%N)
ELAPSED_H=$(( (END_H - START_H) / 1000000 ))

if [ $? -eq 0 ]; then
    pass "Dynamic permits run completes"
else
    fail "Dynamic permits run failed"
fi

echo "  Dynamic permits run: ${ELAPSED_H}ms"

if grep -q 'pf-multi library scheduler (Phase 4)' "$RUN_H/Log.out" 2>/dev/null; then
    pass "Phase 4 scheduler log present with dynamic permits"
    if grep -q 'thread_budget=' "$RUN_H/Log.out"; then
        pass "Per-library thread_budget logged"
    else
        fail "thread_budget not in scheduler log"
    fi
else
    fail "Phase 4 scheduler log missing with dynamic permits"
fi

# =========================================================================
# Part E: Performance comparison (baseline vs dynamic)
# =========================================================================
echo ""
echo "--- Part E: Performance comparison ---"

echo "$ELAPSED_A" > "$OUT_DIR/perf_baseline_ms.txt"
echo "$ELAPSED_H" > "$OUT_DIR/perf_dynamic_ms.txt"

pass "Baseline wall time: ${ELAPSED_A}ms"
pass "Dynamic wall time: ${ELAPSED_H}ms"

# Fail if dynamic run exceeds baseline by more than PHASE5_PERF_TOL_PCT%.
THRESHOLD=$(( ELAPSED_A * (100 + PHASE5_PERF_TOL_PCT) / 100 ))
if [ "$ELAPSED_H" -le "$THRESHOLD" ]; then
    pass "Dynamic <= baseline + ${PHASE5_PERF_TOL_PCT}% (${ELAPSED_H}ms <= ${THRESHOLD}ms)"
else
    fail "Dynamic exceeds baseline + ${PHASE5_PERF_TOL_PCT}% (${ELAPSED_H}ms > ${THRESHOLD}ms)"
fi

# =========================================================================
# Summary
# =========================================================================
echo ""
echo "=========================================="
echo "Phase 5 Validation Results: $PASS passed, $FAIL failed"
echo "=========================================="
echo ""
echo "Output: $OUT_DIR"
echo "Baseline log: $RUN_A/Log.out"
echo "Dynamic log:  $RUN_H/Log.out"
echo "Performance:  baseline=${ELAPSED_A}ms  dynamic=${ELAPSED_H}ms  tolerance=${PHASE5_PERF_TOL_PCT}%"

[ "$FAIL" -eq 0 ]
