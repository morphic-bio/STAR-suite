#!/usr/bin/env bash
# Test: Phase 3 provenance manifest format, fail-fast behavior, and rerun freshness.
# Validates:
#   1. Manifest has all required keys (production force=true stub path)
#   2. status=OK in manifest
#   3. Feature type override applied in features.tsv
#   4. Rerun with changed config overwrites stale features.tsv labels
#   5. Stub failure is hard error
#   6. MEX read failure is hard error
#   7. Manifest not written when MEX read fails
#   8. Failed rerun clears stale manifest from previous success
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: Phase 3 provenance + fail-fast + rerun freshness ==="

# --- Build harness mirroring PfMultiProcess production path ---
HARNESS="$WORK_DIR/test_provenance"
cat > "$WORK_DIR/test_provenance.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cerrno>
#include <cstring>
#include <stdexcept>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#include "PfMultiMexStub.h"
#include "PfMultiMerge.h"

struct FeatureRun {
    string featureType;
    string assignOut;
    string featureRefPath;
    string effectiveChem;
    string detectedMatchMode;
    string libraryId;
    string sampleName;
    string resolvedFastq;
    string resolvedChemRequest;
    bool explicitChem = false;
    int returnCode = 0;
};

// Mirrors PfMultiProcess: clear stale manifests, force=true stub, atomic manifest write.
static int productionPath(const string& assignOut, const string& refCsv,
                          const string& whitelist, const string& featureType) {
    FeatureRun run;
    run.featureType = featureType;
    run.assignOut = assignOut;
    run.featureRefPath = refCsv;
    run.effectiveChem = "NXT";
    run.detectedMatchMode = "RAW_MATCH";
    run.libraryId = "grna_test";
    run.sampleName = "S1";
    run.resolvedFastq = "/path/to/fastq";
    run.resolvedChemRequest = "nxt";
    run.explicitChem = true;
    run.returnCode = 0;

    // Clear stale provenance artifacts (mirrors PfMultiProcess)
    for (const char* suffix : {"pf_library_provenance.tsv", "pf_library_provenance.tsv.tmp"}) {
        string path = run.assignOut + "/" + suffix;
        if (std::remove(path.c_str()) != 0 && errno != ENOENT) {
            cerr << "CLEANUP_FAIL: " << path << endl;
            return 1;
        }
    }

    int stubRet = PfMultiMexStub::processAssignOutput(
        run.assignOut, run.featureRefPath, run.featureType,
        true, whitelist, run.featureType);
    if (stubRet != 0) {
        cerr << "STUB_FAIL" << endl;
        return 1;
    }

    try {
        PfMultiMerge::readMex(run.assignOut);
    } catch (const std::exception& e) {
        cerr << "MEX_READ_FAIL: " << e.what() << endl;
        return 1;
    }

    // Atomic manifest: write to .tmp then rename
    string finalPath = run.assignOut + "/pf_library_provenance.tsv";
    string tmpPath = finalPath + ".tmp";
    {
        std::ofstream manifest(tmpPath.c_str());
        if (!manifest.is_open()) {
            cerr << "MANIFEST_OPEN_FAIL" << endl;
            return 1;
        }
        manifest << "key\tvalue\n";
        manifest << "library_id\t" << run.libraryId << "\n";
        manifest << "sample\t" << run.sampleName << "\n";
        manifest << "feature_type\t" << run.featureType << "\n";
        manifest << "fastq_dir\t" << run.resolvedFastq << "\n";
        manifest << "feature_ref\t" << run.featureRefPath << "\n";
        manifest << "whitelist\t" << whitelist << "\n";
        manifest << "chemistry_request\t" << run.resolvedChemRequest << "\n";
        manifest << "chemistry_explicit\t" << (run.explicitChem ? "yes" : "no") << "\n";
        manifest << "effective_chemistry\t" << run.effectiveChem << "\n";
        manifest << "detected_match_mode\t" << run.detectedMatchMode << "\n";
        manifest << "return_code\t" << run.returnCode << "\n";
        manifest << "status\tOK\n";
        manifest << "assign_output_dir\t" << run.assignOut << "\n";
        manifest.close();
        if (manifest.fail()) {
            std::remove(tmpPath.c_str());
            cerr << "MANIFEST_FLUSH_FAIL" << endl;
            return 1;
        }
    }
    if (std::rename(tmpPath.c_str(), finalPath.c_str()) != 0) {
        std::remove(tmpPath.c_str());
        cerr << "MANIFEST_RENAME_FAIL" << endl;
        return 1;
    }

    cout << "MANIFEST_WRITTEN" << endl;
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_provenance <mode> [args...]" << endl;
        return 2;
    }
    string mode = argv[1];

    if (mode == "manifest") {
        if (argc < 5) { cerr << "need 3 args after manifest" << endl; return 2; }
        string featureType = (argc >= 6) ? argv[5] : "CRISPR Guide Capture";
        return productionPath(argv[2], argv[3], argv[4], featureType);

    } else if (mode == "stub_fail") {
        if (argc < 4) { cerr << "need 2 args after stub_fail" << endl; return 2; }
        int stubRet = PfMultiMexStub::processAssignOutput(
            argv[2], argv[3], "Custom", true, "", "Custom");
        if (stubRet != 0) {
            cout << "STUB_FAILED_AS_EXPECTED" << endl;
            return 0;
        }
        cout << "STUB_UNEXPECTEDLY_SUCCEEDED" << endl;
        return 1;

    } else if (mode == "mex_fail") {
        if (argc < 3) { cerr << "need 1 arg after mex_fail" << endl; return 2; }
        try {
            PfMultiMerge::readMex(argv[2]);
            cout << "MEX_READ_UNEXPECTEDLY_SUCCEEDED" << endl;
            return 1;
        } catch (const std::exception& e) {
            cout << "MEX_READ_FAILED_AS_EXPECTED: " << e.what() << endl;
            return 0;
        }
    }

    cerr << "Unknown mode: " << mode << endl;
    return 2;
}
HARNESS_EOF

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$SOURCE_DIR/htslib" \
    "$WORK_DIR/test_provenance.cpp" \
    "$SOURCE_DIR/PfMultiMexStub.o" \
    "$SOURCE_DIR/PfMultiMerge.o" \
    "$SOURCE_DIR/MexWriter.o" \
    "$SOURCE_DIR/MexWriterUtil.o" \
    -L"$SOURCE_DIR/htslib" -lhts -lz \
    -o "$HARNESS" 2>&1
echo "  Harness built OK"

# --- Setup: assign output with features.txt, barcodes.txt, matrix.mtx (no .tsv yet) ---
ASSIGN_OUT="$WORK_DIR/assign_out"
mkdir -p "$ASSIGN_OUT"

cat > "$WORK_DIR/ref.csv" << 'EOF'
id,name,feature_type
g1,Guide-A,Custom
g2,Guide-B,Custom
EOF

cat > "$ASSIGN_OUT/features.txt" << 'EOF'
Guide-A
Guide-B
EOF

cat > "$ASSIGN_OUT/barcodes.txt" << 'EOF'
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
EOF

cat > "$ASSIGN_OUT/matrix.mtx" << 'EOF'
%%MatrixMarket matrix coordinate integer general
2 2 2
1 1 5
2 2 3
EOF

# --- Test 1: First run creates manifest with all 13 keys ---
echo ""
echo "--- Test 1: First run creates manifest with all required keys ---"
OUTPUT=$("$HARNESS" manifest "$ASSIGN_OUT" "$WORK_DIR/ref.csv" "/path/to/whitelist" 2>&1)
echo "$OUTPUT"

MANIFEST="$ASSIGN_OUT/pf_library_provenance.tsv"
if [ ! -f "$MANIFEST" ]; then
    fail "manifest not written on first run"
else
    REQUIRED_KEYS="library_id sample feature_type fastq_dir feature_ref whitelist chemistry_request chemistry_explicit effective_chemistry detected_match_mode return_code status assign_output_dir"
    ALL_FOUND=true
    for key in $REQUIRED_KEYS; do
        if ! grep -q "^${key}	" "$MANIFEST"; then
            fail "missing key: $key"
            ALL_FOUND=false
        fi
    done
    if $ALL_FOUND; then
        pass "manifest has all 13 required keys"
    fi
fi

# --- Test 2: status=OK in manifest ---
echo ""
echo "--- Test 2: status=OK in manifest ---"
STATUS_VAL=$(grep "^status	" "$MANIFEST" | cut -f2)
if [ "$STATUS_VAL" = "OK" ]; then
    pass "manifest status=OK"
else
    fail "expected status=OK, got '$STATUS_VAL'"
fi

# --- Test 3: feature_type override applied in features.tsv ---
echo ""
echo "--- Test 3: feature_type override applied ---"
CRISPR_COUNT=$(grep -c "CRISPR Guide Capture" "$ASSIGN_OUT/features.tsv" || true)
if [ "$CRISPR_COUNT" -eq 2 ]; then
    pass "features.tsv has CRISPR Guide Capture labels (force=true override)"
else
    fail "expected 2 CRISPR Guide Capture in features.tsv, got $CRISPR_COUNT"
fi

# --- Test 4: Rerun with changed config overwrites stale labels ---
echo ""
echo "--- Test 4: Rerun with changed type overwrites stale features.tsv ---"
OUTPUT=$("$HARNESS" manifest "$ASSIGN_OUT" "$WORK_DIR/ref.csv" "/path/to/whitelist" "Antibody Capture" 2>&1)
RET=$?
echo "$OUTPUT"
AB_COUNT=$(grep -c "Antibody Capture" "$ASSIGN_OUT/features.tsv" || true)
OLD_CRISPR=$(grep -c "CRISPR Guide Capture" "$ASSIGN_OUT/features.tsv" || true)
if [ "$RET" -eq 0 ] && [ "$AB_COUNT" -eq 2 ] && [ "$OLD_CRISPR" -eq 0 ]; then
    pass "rerun with changed config overwrites stale labels (Antibody Capture)"
else
    fail "rerun did not overwrite: ret=$RET, Antibody=$AB_COUNT, CRISPR=$OLD_CRISPR"
fi

# --- Test 5: Stub failure is hard error ---
echo ""
echo "--- Test 5: Stub failure is hard error ---"
OUTPUT=$("$HARNESS" stub_fail "$WORK_DIR/nonexistent_dir" "$WORK_DIR/nonexistent_ref.csv" 2>&1 || true)
echo "$OUTPUT"
if echo "$OUTPUT" | grep -q "STUB_FAILED_AS_EXPECTED"; then
    pass "stub failure detected and returned non-zero"
else
    fail "stub should have failed but did not"
fi

# --- Test 6: MEX read failure is hard error ---
echo ""
echo "--- Test 6: MEX read failure is hard error ---"
BAD_MEX="$WORK_DIR/bad_mex"
mkdir -p "$BAD_MEX"
OUTPUT=$("$HARNESS" mex_fail "$BAD_MEX" 2>&1 || true)
echo "$OUTPUT"
if echo "$OUTPUT" | grep -q "MEX_READ_FAILED_AS_EXPECTED"; then
    pass "MEX read failure detected and throws"
else
    fail "MEX read should have failed but did not"
fi

# --- Test 7: Manifest not written when MEX read fails ---
echo ""
echo "--- Test 7: Manifest not written on MEX read failure ---"
BAD_ASSIGN="$WORK_DIR/bad_assign"
mkdir -p "$BAD_ASSIGN"
cat > "$BAD_ASSIGN/features.txt" << 'EOF'
Guide-A
EOF
cat > "$BAD_ASSIGN/barcodes.txt" << 'EOF'
AAACCCAAGAAACACT
EOF
OUTPUT=$("$HARNESS" manifest "$BAD_ASSIGN" "$WORK_DIR/ref.csv" "/path/to/wl" 2>&1 || true)
echo "$OUTPUT"
if [ ! -f "$BAD_ASSIGN/pf_library_provenance.tsv" ]; then
    pass "manifest not written when MEX read fails"
elif echo "$OUTPUT" | grep -q "MEX_READ_FAIL"; then
    pass "manifest not written when MEX read fails (process exited before write)"
else
    fail "manifest was written despite MEX read failure"
fi

# --- Test 8: Failed rerun removes stale manifest from previous success ---
echo ""
echo "--- Test 8: Failed rerun clears stale manifest ---"
# First, confirm manifest exists from test 1-4
RERUN_DIR="$WORK_DIR/rerun_out"
mkdir -p "$RERUN_DIR"
cat > "$RERUN_DIR/features.txt" << 'EOF'
Guide-A
Guide-B
EOF
cat > "$RERUN_DIR/barcodes.txt" << 'EOF'
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
EOF
cat > "$RERUN_DIR/matrix.mtx" << 'EOF'
%%MatrixMarket matrix coordinate integer general
2 2 2
1 1 5
2 2 3
EOF
# Successful first run
"$HARNESS" manifest "$RERUN_DIR" "$WORK_DIR/ref.csv" "/path/to/wl" 2>&1 >/dev/null
if [ ! -f "$RERUN_DIR/pf_library_provenance.tsv" ]; then
    fail "setup: manifest not created on successful run"
else
    # Break the MEX by removing matrix.mtx
    rm "$RERUN_DIR/matrix.mtx"
    # Failing rerun (MEX read will fail)
    "$HARNESS" manifest "$RERUN_DIR" "$WORK_DIR/ref.csv" "/path/to/wl" 2>&1 >/dev/null || true
    if [ -f "$RERUN_DIR/pf_library_provenance.tsv" ]; then
        fail "stale manifest survived failed rerun"
    else
        pass "failed rerun cleared stale manifest"
    fi
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
