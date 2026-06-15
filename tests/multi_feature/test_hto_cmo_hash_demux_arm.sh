#!/usr/bin/env bash
# Test: pf-multi HTO/CMO hash demux config columns and adt_mex routing
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d /tmp/hto_cmo_hash_demux_arm.XXXXXX)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: HTO/CMO hash demux pf-multi arm ==="

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
HARNESS="$WORK_DIR/test_hto_hash_arm"
cat > "$WORK_DIR/test_hto_hash_arm.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;
#include "PfMultiConfig.h"
#include "PfMultiFeatureSpecs.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_hto_hash_arm <config.csv>" << endl;
        return 2;
    }
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        ostringstream log;
        vector<PfMultiFeatureSpecs::FeatureSpec> specs =
            PfMultiFeatureSpecs::buildFeatureSpecsFromConfig(cfg, log);
        for (const auto& spec : specs) {
            vector<PfMultiConfig::LibraryEntry> libs =
                cfg.getFeatureLibraries(spec.libraryType);
            for (const auto& lib : libs) {
                const bool adtMex = PfMultiFeatureSpecs::shouldEmitAdtMexOutput(
                    spec.featureRefType, spec.libraryType);
                cout << "LIB id=" << lib.starLibraryId
                     << " feature_types=" << spec.libraryType
                     << " featureRefType=" << spec.featureRefType
                     << " adt_mex=" << (adtMex ? "yes" : "no")
                     << " star_hash_demux=" << lib.starHashDemux
                     << " star_hash_feature_selector=" << lib.starHashFeatureSelector
                     << " star_hash_demux_method=" << lib.starHashDemuxMethod
                     << " star_hash_min_total=" << lib.starHashMinTotal
                     << " star_hash_min_top=" << lib.starHashMinTop
                     << " star_hash_min_ratio=" << lib.starHashMinRatio
                     << endl;
            }
        }
        return 0;
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
}
HARNESS_EOF

g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    "$WORK_DIR/test_hto_hash_arm.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS"

touch "$WORK_DIR/hto_ref.csv"
touch "$WORK_DIR/adt_hto_ref.csv"
touch "$WORK_DIR/whitelist.txt"

cat > "$WORK_DIR/config_hto.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_hash_demux,star_hash_min_total,star_hash_min_top,star_hash_min_ratio
$WORK_DIR/hto_fastqs,S1,Multiplexing Capture,Multiplexing Capture,TRU,hto_s1,$WORK_DIR/hto_ref.csv,$WORK_DIR/whitelist.txt,yes,3,3,2.0
$WORK_DIR/adt_fastqs,S1,Antibody Capture,Protein,TRU,adt_s1,$WORK_DIR/adt_hto_ref.csv,$WORK_DIR/whitelist.txt,auto,,, ,id_prefix:hashtag

[feature]
ref,$WORK_DIR/hto_ref.csv
EOF

# Fix trailing selector row (avoid broken CSV from spaces)
cat > "$WORK_DIR/config_hto.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_hash_demux,star_hash_feature_selector,star_hash_demux_method,star_hash_min_total,star_hash_min_top,star_hash_min_ratio
$WORK_DIR/hto_fastqs,S1,Multiplexing Capture,Multiplexing Capture,TRU,hto_s1,$WORK_DIR/hto_ref.csv,$WORK_DIR/whitelist.txt,yes,,ratio,3,3,2.0
$WORK_DIR/adt_fastqs,S1,Antibody Capture,Protein,TRU,adt_s1,$WORK_DIR/adt_hto_ref.csv,$WORK_DIR/whitelist.txt,auto,id_prefix:hashtag,ratio,3,3,2.0

[feature]
ref,$WORK_DIR/hto_ref.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_hto.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'id=hto_s1 .*adt_mex=yes .*star_hash_demux=yes .*star_hash_min_total=3'; then
    pass "hash-only library routes adt_mex with star_hash_demux=yes"
else
    fail "hash-only library routing"
fi

if echo "$OUTPUT" | grep -q 'id=adt_s1 .*adt_mex=yes .*star_hash_feature_selector=id_prefix:hashtag'; then
    pass "mixed ADT library routes adt_mex with hash selector"
else
    fail "mixed ADT library routing"
fi

if echo "$OUTPUT" | grep -q 'featureRefType=Multiplexing Capture'; then
    pass "HTO feature_types maps to Multiplexing Capture"
else
    fail "HTO featureRefType mapping"
fi

echo ""
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
if [ "$FAIL_COUNT" -gt 0 ]; then
    exit 1
fi
