#!/usr/bin/env bash
# Test: mixed-chemistry preparation and precedence in PfMultiProcess
#
# Exercises the actual buildPfMultiPreparedContext() logic rather than only
# config parsing. Covers:
#   - GEX star_chemistry re-anchoring of inferred/effective chemistry
#   - per-library override isolation in multi-library configs
#   - distinction between empty star_chemistry (inherit global flag explicitly)
#     and explicit "auto" (remain auto-detect eligible)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
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

echo "=== Test: mixed-chemistry prepared context ==="

require_file "${REPO_ROOT}/core/legacy/source/STAR" "STAR binary"

mkdir -p "${WORK_DIR}/gex" "${WORK_DIR}/grna_a" "${WORK_DIR}/grna_b"
mkdir -p "${WORK_DIR}/out_case1" "${WORK_DIR}/out_case2" "${WORK_DIR}/out_case3"

cat > "${WORK_DIR}/feature_ref_grna.csv" << 'EOF'
id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name
Guide_A,Guide_A,R2,TTCCAGCTTAGCTCTTAAAC(BC),GGGTGGTGCCCATCCTGGTC,CRISPR Guide Capture,target_1,Guide_A
EOF

# 2-column NXT-first whitelist: col1=NXT, col2=TRU
cat > "${WORK_DIR}/whitelist_nxt_2col.txt" << 'EOF'
AAAAAAAAAAAAAAAA	AAAAAAATTAAAAAAA
CCCCCCCCCCCCCCCC	CCCCCCCGGCCCCCCC
EOF

# 2-column TRU-first whitelist: col1=TRU, col2=NXT
cat > "${WORK_DIR}/whitelist_tru_2col.txt" << 'EOF'
AAAAAAATTAAAAAAA	AAAAAAAAAAAAAAAA
CCCCCCCGGCCCCCCC	CCCCCCCCCCCCCCCC
EOF

HARNESS="${WORK_DIR}/test_prepared_context"
cat > "${WORK_DIR}/test_prepared_context.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <stdexcept>
#include <string>
#include "PfMultiProcess.cpp"

static void dumpContext(const PfMultiPreparedContext& ctx) {
    std::cout << "CTX"
              << " requested=" << ctx.requestedChem
              << " inferred=" << ctx.inferredChem
              << " confident=" << (ctx.inferredChemConfident ? "yes" : "no")
              << " effective=" << ctx.effectiveChem
              << " hasTwoColumnWL=" << (ctx.hasTwoColumnWhitelist ? "yes" : "no")
              << " soloOutputNamespace=" << ctx.soloOutputNamespace
              << " outputChem=" << ctx.outputChem
              << "\n";
    for (size_t i = 0; i < ctx.featureLibraries.size(); ++i) {
        const auto& lib = ctx.featureLibraries[i];
        std::cout << "LIB:" << i
                  << " id=" << lib.libraryId
                  << " type=" << lib.libraryType
                  << " resolved=" << lib.resolvedChemRequest
                  << " effective=" << lib.effectiveChem
                  << " explicit=" << (lib.explicitChem ? 1 : 0)
                  << "\n";
    }
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: test_prepared_context <config.csv> <whitelist> <crChemistry> <outPrefix>\n";
        return 2;
    }

    PfMultiPreloadInput input;
    input.pfMultiConfig = argv[1];
    input.crWhitelist = argv[2];
    input.crChemistry = argv[3];
    input.crOutputChemistry = "";
    input.outFileNamePrefix = argv[4];

    try {
        PfMultiPreparedContext ctx = buildPfMultiPreparedContext(input);
        dumpContext(ctx);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
HARNESS_EOF

INC_FLAGS="-I${REPO_ROOT}/core/legacy/source -I${REPO_ROOT}/core/features/libscrna/include -I${REPO_ROOT}/core/features/vbem/source -I${REPO_ROOT}/core/features/vbem/source/libem -I${REPO_ROOT}/slam/source -I${REPO_ROOT}/core/features/bamsort/source -I${REPO_ROOT}/core/features/process_features/include -I${REPO_ROOT}/flex/source -I${REPO_ROOT}/flex/source/libflex -I${REPO_ROOT}/core/legacy/source/htslib"

g++ -std=c++11 -O2 -fopenmp -ffunction-sections -fdata-sections \
    ${INC_FLAGS} \
    "${WORK_DIR}/test_prepared_context.cpp" \
    "${REPO_ROOT}/core/legacy/source/PfMultiConfig.o" \
    -Wl,--gc-sections \
    -L"${REPO_ROOT}/core/features/process_features" -lprocess_features \
    -L"${REPO_ROOT}/core/features/libscrna" -lscrna \
    -lpthread -lz -lstdc++ -lhts -lglib-2.0 \
    -o "${HARNESS}"

echo "  Harness built OK"

echo ""
echo "--- Test 1: GEX override re-anchors auto feature libs and isolates explicit overrides ---"
cat > "${WORK_DIR}/config_case1.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
${WORK_DIR}/gex,S1,Gene Expression,Gene Expression,TRU,,gex_main
${WORK_DIR}/grna_a,S1,CRISPR Guide Capture,CRISPR Guide Capture,,$WORK_DIR/feature_ref_grna.csv,grna_auto
${WORK_DIR}/grna_b,S1,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$WORK_DIR/feature_ref_grna.csv,grna_explicit
EOF

OUTPUT=$("${HARNESS}" "${WORK_DIR}/config_case1.csv" "${WORK_DIR}/whitelist_nxt_2col.txt" auto "${WORK_DIR}/out_case1" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'CTX requested=auto inferred=TRU confident=yes effective=TRU .*soloOutputNamespace=NXT'; then
    pass "GEX star_chemistry=TRU re-anchors inferred/effective chemistry and solo output namespace"
else
    fail "GEX override should re-anchor context from NXT whitelist to TRU effective chemistry"
fi

if echo "$OUTPUT" | grep -q 'LIB:0 id=grna_auto .* resolved=auto effective=TRU explicit=0'; then
    pass "auto feature lib follows re-anchored TRU context without becoming explicit"
else
    fail "auto feature lib should inherit re-anchored TRU context"
fi

if echo "$OUTPUT" | grep -q 'LIB:1 id=grna_explicit .* resolved=nxt effective=NXT explicit=1'; then
    pass "explicit feature NXT override remains isolated from GEX TRU anchor"
else
    fail "explicit feature override should remain NXT and explicit"
fi

echo ""
echo "--- Test 2: explicit star_chemistry=auto stays auto-detect eligible under global flag ---"
cat > "${WORK_DIR}/config_case2.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
${WORK_DIR}/gex,S1,Gene Expression,Gene Expression,,,
${WORK_DIR}/grna_a,S1,CRISPR Guide Capture,CRISPR Guide Capture,auto,$WORK_DIR/feature_ref_grna.csv,grna_auto_flag
EOF

OUTPUT=$("${HARNESS}" "${WORK_DIR}/config_case2.csv" "${WORK_DIR}/whitelist_tru_2col.txt" tru "${WORK_DIR}/out_case2" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'CTX requested=tru inferred=TRU confident=yes effective=TRU .*soloOutputNamespace=NXT'; then
    pass "global TRU flag keeps TRU effective chemistry with TRU-first whitelist"
else
    fail "global TRU flag should yield TRU effective chemistry"
fi

if echo "$OUTPUT" | grep -q 'LIB:0 id=grna_auto_flag .* resolved=auto effective=TRU explicit=0'; then
    pass "explicit auto column remains auto-detect eligible under global TRU flag"
else
    fail "star_chemistry=auto should not become explicit under global TRU"
fi

echo ""
echo "--- Test 3: empty star_chemistry inherits global flag explicitly ---"
cat > "${WORK_DIR}/config_case3.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
${WORK_DIR}/gex,S1,Gene Expression,Gene Expression,,,
${WORK_DIR}/grna_a,S1,CRISPR Guide Capture,CRISPR Guide Capture,,$WORK_DIR/feature_ref_grna.csv,grna_inherit_flag
EOF

OUTPUT=$("${HARNESS}" "${WORK_DIR}/config_case3.csv" "${WORK_DIR}/whitelist_tru_2col.txt" tru "${WORK_DIR}/out_case3" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0 id=grna_inherit_flag .* resolved=tru effective=TRU explicit=1'; then
    pass "empty star_chemistry inherits global flag as explicit chemistry"
else
    fail "empty star_chemistry should inherit global TRU explicitly"
fi

echo ""
echo "=========================================="
echo "Results: ${PASS_COUNT} passed, ${FAIL_COUNT} failed"
echo "=========================================="

if [[ ${FAIL_COUNT} -ne 0 ]]; then
    exit 1
fi
