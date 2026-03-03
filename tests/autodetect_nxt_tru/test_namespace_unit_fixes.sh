#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

require_file() {
    local p="$1"
    local why="$2"
    if [[ ! -e "$p" ]]; then
        echo "ERROR: missing ${why}: ${p}" >&2
        exit 1
    fi
}

echo "=== Namespace Unit Tests (Fixes 2/3/4/5) ==="

require_file "${STAR_BIN}" "STAR binary"
require_file "${REPO_ROOT}/core/legacy/source/PfMultiAssign.o" "PfMultiAssign.o"
require_file "${REPO_ROOT}/core/legacy/source/PfMultiConfig.o" "PfMultiConfig.o"
require_file "${REPO_ROOT}/core/legacy/source/ThreadControl.o" "ThreadControl.o"

cat > "${WORK_DIR}/wl_unknown.txt" << 'EOF'
AAAAAAAAAAAAAAAA
CCCCCCCCCCCCCCCC
EOF

# col2 is col1 with only positions 8 and 9 complemented (NXT/TRU rule)
cat > "${WORK_DIR}/wl_nxt_2col.txt" << 'EOF'
AAAAAAAAAAAAAAAA	AAAAAAATTAAAAAAA
CCCCCCCCCCCCCCCC	CCCCCCCGGCCCCCCC
EOF

# Same 2-column content but filename has no NXT/TRU hint — orientation ambiguous
cat > "${WORK_DIR}/wl_custom_2col.txt" << 'EOF'
AAAAAAAAAAAAAAAA	AAAAAAATTAAAAAAA
CCCCCCCCCCCCCCCC	CCCCCCCGGCCCCCCC
EOF

cat > "${WORK_DIR}/assign_wl.txt" << 'EOF'
AAAAAAAAAAAAAAAA
CCCCCCCCCCCCCCCC
EOF

echo ""
echo "--- Test 1: PfMultiProcess static helpers (detect + explicit filtered normalization) ---"
cat > "${WORK_DIR}/test_pf_multi_helpers.cpp" << 'EOF'
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include "PfMultiProcess.cpp"

static int assert_true(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "ASSERT_FAIL: " << msg << "\n";
        return 1;
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: test_pf_multi_helpers <unknown_wl> <nxt_2col_wl> <assign_wl> <ambiguous_2col_wl>\n";
        return 2;
    }

    int rc = 0;
    std::string reason;
    bool confident = true;

    std::string chem = detectChemistryFromWhitelistPath(argv[1], reason, confident);
    rc |= assert_true(chem == "UNKNOWN", "unknown whitelist should resolve to UNKNOWN");
    rc |= assert_true(!confident, "unknown whitelist should be non-confident");

    reason.clear();
    confident = false;
    chem = detectChemistryFromWhitelistPath(argv[2], reason, confident);
    rc |= assert_true(chem == "NXT", "2-column whitelist with 'nxt' filename should resolve to NXT");
    rc |= assert_true(confident, "2-column whitelist with 'nxt' filename should be confident");

    // Ambiguous 2-column: complement rule matches but no filename hint
    reason.clear();
    confident = true;
    chem = detectChemistryFromWhitelistPath(argv[4], reason, confident);
    rc |= assert_true(chem == "NXT", "ambiguous 2-col should still resolve to NXT");
    rc |= assert_true(!confident, "ambiguous 2-col (no filename hint) should be non-confident");
    rc |= assert_true(reason.find("WARNING") != std::string::npos,
                      "ambiguous 2-col reason should contain WARNING");

    std::unordered_set<std::string> whitelistSet;
    uint64_t rowsRead = 0;
    uint64_t invalidRows = 0;
    rc |= assert_true(loadWhitelistBarcodeSet(argv[3], whitelistSet, rowsRead, invalidRows),
                      "loadWhitelistBarcodeSet failed");
    rc |= assert_true(whitelistSet.size() == 2, "expected two whitelist barcodes");

    std::vector<std::string> src;
    src.push_back("AAAAAAAAAAAAAAAA");                         // in_set
    src.push_back(translateNxtMiddleTwoBases("CCCCCCCCCCCCCCCC")); // translated_to_set
    src.push_back(translateNxtMiddleTwoBases("CCCCCCCCCCCCCCCC")); // translated_to_set duplicate
    src.push_back("GGGGGGGGGGGGGGGG");                         // unmatched

    FilteredBarcodeNormalizationStats stats;
    std::vector<std::string> out =
        normalizeFilteredBarcodesForAssignNamespace(src, whitelistSet, stats);

    rc |= assert_true(stats.inSet == 1, "inSet should be 1");
    rc |= assert_true(stats.translatedToSet == 2, "translatedToSet should be 2");
    rc |= assert_true(stats.unmatched == 1, "unmatched should be 1");
    rc |= assert_true(stats.dedupDropped == 1, "dedupDropped should be 1");
    rc |= assert_true(stats.outputCount == 2, "outputCount should be 2");
    rc |= assert_true(out.size() == 2, "normalized output size should be 2");

    return rc;
}
EOF

INC_FLAGS="-I${REPO_ROOT}/core/legacy/source -I${REPO_ROOT}/core/features/libscrna/include -I${REPO_ROOT}/core/features/vbem/source -I${REPO_ROOT}/core/features/vbem/source/libem -I${REPO_ROOT}/slam/source -I${REPO_ROOT}/core/features/bamsort/source -I${REPO_ROOT}/core/features/process_features/include -I${REPO_ROOT}/flex/source -I${REPO_ROOT}/flex/source/libflex -I${REPO_ROOT}/core/legacy/source/htslib"

g++ -std=c++11 -O2 -fopenmp -ffunction-sections -fdata-sections \
    ${INC_FLAGS} \
    "${WORK_DIR}/test_pf_multi_helpers.cpp" \
    -Wl,--gc-sections \
    -L"${REPO_ROOT}/core/features/process_features" -lprocess_features \
    -L"${REPO_ROOT}/core/features/libscrna" -lscrna \
    -lpthread -lz -lstdc++ -lhts -lglib-2.0 \
    -o "${WORK_DIR}/test_pf_multi_helpers"

if "${WORK_DIR}/test_pf_multi_helpers" "${WORK_DIR}/wl_unknown.txt" "${WORK_DIR}/wl_nxt_2col.txt" "${WORK_DIR}/assign_wl.txt" "${WORK_DIR}/wl_custom_2col.txt"; then
    pass "PfMultiProcess helper behavior"
else
    fail "PfMultiProcess helper behavior"
fi

echo ""
echo "--- Test 2: PfMultiAssign whitelist normalization metadata ---"
cat > "${WORK_DIR}/test_pf_multi_assign_norm.cpp" << 'EOF'
#include "PfMultiAssign.h"
#include "ThreadControl.h"
#include <fstream>
#include <iostream>
#include <string>

ThreadControl g_threadChunks;

static int assert_true(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "ASSERT_FAIL: " << msg << "\n";
        return 1;
    }
    return 0;
}

static int count_nonempty_lines(const std::string& path) {
    std::ifstream in(path.c_str());
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        if (!line.empty()) n++;
    }
    return n;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: test_pf_multi_assign_norm <wl_2col_nxt> <wl_1col_unknown> <outdir> <wl_2col_nohint>\n";
        return 2;
    }
    int rc = 0;
    auto two = PfMultiAssign::normalizeWhitelistForAssign(argv[1], argv[3]);
    rc |= assert_true(two.hasTwoColumnSource, "two-column whitelist must set hasTwoColumnSource");
    rc |= assert_true(two.assignmentNamespace == "NXT", "two-column translation whitelist must infer NXT");
    rc |= assert_true(two.namespaceConfidence,
                      "two-column whitelist with 'nxt' in filename should be confident");
    rc |= assert_true(two.normalizedRowCount == 2, "two-column normalizedRowCount should be 2");
    rc |= assert_true(two.normalizedPath.find("whitelist.normalized.txt") != std::string::npos,
                      "two-column normalization path should be whitelist.normalized.txt");
    rc |= assert_true(count_nonempty_lines(two.normalizedPath) == 2, "normalized whitelist file should have 2 lines");

    auto one = PfMultiAssign::normalizeWhitelistForAssign(argv[2], argv[3]);
    rc |= assert_true(!one.hasTwoColumnSource, "one-column whitelist must not set hasTwoColumnSource");
    rc |= assert_true(one.assignmentNamespace == "UNKNOWN", "one-column unknown whitelist should infer UNKNOWN");
    rc |= assert_true(!one.namespaceConfidence, "one-column unknown whitelist should be non-confident");
    rc |= assert_true(one.normalizedRowCount == 2, "one-column normalizedRowCount should be 2");
    rc |= assert_true(one.normalizedPath == std::string(argv[2]), "one-column normalized path should stay source path");

    // Ambiguous 2-column: complement rule matches but no filename hint
    auto ambig = PfMultiAssign::normalizeWhitelistForAssign(argv[4], argv[3]);
    rc |= assert_true(ambig.hasTwoColumnSource, "ambiguous 2-col must set hasTwoColumnSource");
    rc |= assert_true(ambig.assignmentNamespace == "NXT", "ambiguous 2-col must still infer NXT");
    rc |= assert_true(!ambig.namespaceConfidence,
                      "ambiguous 2-col (no filename hint) must be non-confident");
    rc |= assert_true(ambig.normalizedRowCount == 2, "ambiguous 2-col normalizedRowCount should be 2");

    return rc;
}
EOF

g++ -std=c++11 -O2 -fopenmp \
    ${INC_FLAGS} \
    "${WORK_DIR}/test_pf_multi_assign_norm.cpp" \
    "${REPO_ROOT}/core/legacy/source/PfMultiAssign.o" \
    "${REPO_ROOT}/core/legacy/source/ThreadControl.o" \
    "${REPO_ROOT}/core/legacy/source/PfMultiConfig.o" \
    -L"${REPO_ROOT}/core/features/process_features" -lprocess_features \
    -L"${REPO_ROOT}/core/features/libscrna" -lscrna \
    -lpthread -lz -lstdc++ -lhts -lglib-2.0 \
    -o "${WORK_DIR}/test_pf_multi_assign_norm"

mkdir -p "${WORK_DIR}/assign_out"
if "${WORK_DIR}/test_pf_multi_assign_norm" "${WORK_DIR}/wl_nxt_2col.txt" "${WORK_DIR}/wl_unknown.txt" "${WORK_DIR}/assign_out" "${WORK_DIR}/wl_custom_2col.txt"; then
    pass "PfMultiAssign normalization metadata"
else
    fail "PfMultiAssign normalization metadata"
fi

echo ""
echo "--- Test 3: Deterministic normalization with COL1/COL2 overlap barcodes ---"
cat > "${WORK_DIR}/test_deterministic_norm.cpp" << 'CPPEOF'
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
#include "PfMultiProcess.cpp"

static int assert_true(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "ASSERT_FAIL: " << msg << "\n";
        return 1;
    }
    return 0;
}

int main() {
    // Build a whitelist set (COL1 = NXT) with deliberate overlap:
    //   col1_a = AAAAAAAACAAAAAAA  →  col2_a = AAAAAAATGAAAAAAA
    //   col1_b = AAAAAAATGAAAAAAA  →  col2_b = AAAAAAAACAAAAAAA
    //   col1_c = CCCCCCCGCCCCCCCC  →  col2_c = CCCCCCCCGCCCCCCC
    // Overlap: col1_b == col2_a and col2_b == col1_a

    std::string col1_a = "AAAAAAAACAAAAAAA";
    std::string col2_a = translateNxtMiddleTwoBases(col1_a); // AAAAAAATGAAAAAAA
    std::string col1_b = col2_a;                             // overlaps with col2_a
    std::string col2_b = translateNxtMiddleTwoBases(col1_b); // should == col1_a
    std::string col1_c = "CCCCCCCGCCCCCCCC";
    std::string col2_c = translateNxtMiddleTwoBases(col1_c);

    int rc = 0;
    rc |= assert_true(col1_b == col2_a, "setup: col1_b must overlap col2_a");
    rc |= assert_true(col2_b == col1_a, "setup: col2_b must overlap col1_a");

    // NXT whitelist set (COL1 only)
    std::unordered_set<std::string> nxtSet;
    nxtSet.insert(col1_a);
    nxtSet.insert(col1_b);
    nxtSet.insert(col1_c);

    // Input: TRU barcodes (Solo output for 2-col NXT)
    std::vector<std::string> truBarcodes;
    truBarcodes.push_back(col2_a);  // AAAAAAATGAAAAAAA — overlaps col1_b!
    truBarcodes.push_back(col2_b);  // AAAAAAAACAAAAAAA — overlaps col1_a!
    truBarcodes.push_back(col2_c);  // CCCCCCCCGCCCCCCC — no overlap

    // --- Heuristic path (BUG): no sourceNamespace, membership decides ---
    FilteredBarcodeNormalizationStats heuristicStats;
    std::vector<std::string> heuristicOut =
        normalizeFilteredBarcodesForAssignNamespace(truBarcodes, nxtSet, heuristicStats);

    // With heuristic, col2_a matches col1_b in set → inSet (WRONG: kept as TRU)
    // col2_b matches col1_a in set → inSet (WRONG: kept as TRU)
    // col2_c doesn't match → translated to NXT
    rc |= assert_true(heuristicStats.inSet == 2,
        "heuristic: expected inSet=2 (overlap collisions)");
    rc |= assert_true(heuristicStats.translatedToSet == 1,
        "heuristic: expected translatedToSet=1");

    // --- Deterministic path (FIX): sourceNS=TRU, whitelistNS=NXT ---
    FilteredBarcodeNormalizationStats deterministicStats;
    std::vector<std::string> deterministicOut =
        normalizeFilteredBarcodesForAssignNamespace(
            truBarcodes, nxtSet, deterministicStats, "TRU", "NXT");

    // Deterministic: always translates TRU→NXT, never uses membership guess
    rc |= assert_true(deterministicStats.inSet == 0,
        "deterministic: expected inSet=0 (no membership guessing)");
    rc |= assert_true(deterministicStats.translatedToSet == 3,
        "deterministic: expected translatedToSet=3 (all translated)");
    rc |= assert_true(deterministicStats.unmatched == 0,
        "deterministic: expected unmatched=0");
    rc |= assert_true(deterministicStats.dedupDropped == 0,
        "deterministic: expected dedupDropped=0");
    rc |= assert_true(deterministicStats.outputCount == 3,
        "deterministic: expected outputCount=3");

    // Verify output barcodes are COL1 (NXT), not COL2 (TRU)
    rc |= assert_true(deterministicOut.size() == 3,
        "deterministic: output size should be 3");
    rc |= assert_true(nxtSet.count(deterministicOut[0]) == 1,
        "deterministic: output[0] must be in NXT whitelist set");
    rc |= assert_true(nxtSet.count(deterministicOut[1]) == 1,
        "deterministic: output[1] must be in NXT whitelist set");
    rc |= assert_true(nxtSet.count(deterministicOut[2]) == 1,
        "deterministic: output[2] must be in NXT whitelist set");

    // The translated barcodes should map back to the correct COL1 entries:
    // TRU col2_a → NXT col1_a (NOT col1_b which it overlaps)
    rc |= assert_true(deterministicOut[0] == col1_a,
        "deterministic: col2_a should translate to col1_a, not col1_b");
    rc |= assert_true(deterministicOut[1] == col1_b,
        "deterministic: col2_b should translate to col1_b");
    rc |= assert_true(deterministicOut[2] == col1_c,
        "deterministic: col2_c should translate to col1_c");

    // --- Same-namespace deterministic (sourceNS=NXT, whitelistNS=NXT): no translation ---
    FilteredBarcodeNormalizationStats sameNsStats;
    std::vector<std::string> nxtBarcodes;
    nxtBarcodes.push_back(col1_a);
    nxtBarcodes.push_back(col1_b);
    nxtBarcodes.push_back(col1_c);
    std::vector<std::string> sameNsOut =
        normalizeFilteredBarcodesForAssignNamespace(
            nxtBarcodes, nxtSet, sameNsStats, "NXT", "NXT");

    rc |= assert_true(sameNsStats.inSet == 3,
        "same-NS deterministic: expected inSet=3 (no translation needed)");
    rc |= assert_true(sameNsStats.translatedToSet == 0,
        "same-NS deterministic: expected translatedToSet=0");
    rc |= assert_true(sameNsOut.size() == 3,
        "same-NS deterministic: output size should be 3");

    if (rc == 0) {
        std::cout << "All deterministic normalization assertions passed\n";
    }
    return rc;
}
CPPEOF

g++ -std=c++11 -O2 -fopenmp -ffunction-sections -fdata-sections \
    ${INC_FLAGS} \
    "${WORK_DIR}/test_deterministic_norm.cpp" \
    -Wl,--gc-sections \
    -L"${REPO_ROOT}/core/features/process_features" -lprocess_features \
    -L"${REPO_ROOT}/core/features/libscrna" -lscrna \
    -lpthread -lz -lstdc++ -lhts -lglib-2.0 \
    -o "${WORK_DIR}/test_deterministic_norm"

if "${WORK_DIR}/test_deterministic_norm"; then
    pass "Deterministic normalization with overlap barcodes"
else
    fail "Deterministic normalization with overlap barcodes"
fi

echo ""
echo "=========================================="
echo "Results: ${PASS_COUNT} passed, ${FAIL_COUNT} failed"
echo "=========================================="

if [[ ${FAIL_COUNT} -ne 0 ]]; then
    exit 1
fi

