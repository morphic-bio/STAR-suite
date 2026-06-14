#!/usr/bin/env bash
# Tiny multi-feature smoke: gRNA + Custom + ADT/protein arms share the pf-multi
# library routing path. ADT/protein libraries trigger assignBarcodes ADT MEX mode.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PF_DIR="$REPO_ROOT/core/features/process_features"
ASSIGN="$PF_DIR/assignBarcodes"
WORK_DIR="$(mktemp -d /tmp/multi_adt_arm_test.XXXXXX)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: multi-feature ADT/protein library arm ==="

if [[ ! -x "$ASSIGN" ]]; then
    echo "Building assignBarcodes..."
    make -C "$PF_DIR" -j8 assignBarcodes
fi

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
if [[ ! -f "$SOURCE_DIR/PfMultiConfig.o" ]]; then
    make -C "$SOURCE_DIR" -j8 PfMultiConfig.o
fi

HARNESS="$WORK_DIR/test_adt_arm"
cat > "$WORK_DIR/test_adt_arm.cpp" << 'HARNESS_EOF'
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

static string sanitizeDirName(const string& input) {
    string out = input;
    for (char& c : out) {
        if (!((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
              (c >= '0' && c <= '9') || c == '-' || c == '_')) {
            c = '_';
        }
    }
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_adt_arm <config.csv> <out_prefix>" << endl;
        return 2;
    }
    const string outPrefix = argv[2];
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        ostringstream log;
        vector<PfMultiFeatureSpecs::FeatureSpec> specs =
            PfMultiFeatureSpecs::buildFeatureSpecsFromConfig(cfg, log);

        for (size_t i = 0; i < specs.size(); ++i) {
            const auto& spec = specs[i];
            vector<PfMultiConfig::LibraryEntry> libs =
                cfg.getFeatureLibraries(spec.libraryType);
            for (size_t j = 0; j < libs.size(); ++j) {
                const auto& lib = libs[j];
                const string featureDir = sanitizeDirName(spec.libraryType);
                const string assignOut = outPrefix + "/cr_assign/" + featureDir + "/"
                    + sanitizeDirName(lib.starLibraryId);
                const bool adtMex = PfMultiFeatureSpecs::shouldEmitAdtMexOutput(
                    spec.featureRefType, spec.libraryType);
                cout << "LIBRARY"
                     << " id=" << lib.starLibraryId
                     << " feature_types=" << spec.libraryType
                     << " featureRefType=" << spec.featureRefType
                     << " assignOut=" << assignOut
                     << " adt_mex=" << (adtMex ? "yes" : "no")
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
    "$WORK_DIR/test_adt_arm.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS"

touch "$WORK_DIR/ref_grna.csv" "$WORK_DIR/ref_larry.csv" "$WORK_DIR/ref_adt.csv"

cat > "$WORK_DIR/multi_config.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref,star_library_id
$WORK_DIR/fastq_grna,PolyIII,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv,grna_de
$WORK_DIR/fastq_larry,LARRY,Custom,$WORK_DIR/ref_larry.csv,larry_de
$WORK_DIR/fastq_adt,ADT,Protein,$WORK_DIR/ref_adt.csv,adt_de
EOF

ROUTING="$("$HARNESS" "$WORK_DIR/multi_config.csv" "$WORK_DIR/out" 2>&1)"
echo "$ROUTING"

if echo "$ROUTING" | grep -q 'id=grna_de .*adt_mex=no'; then
    pass "gRNA library does not request ADT MEX"
else
    fail "gRNA library should not request ADT MEX"
fi

if echo "$ROUTING" | grep -q 'id=larry_de .*adt_mex=no'; then
    pass "Custom/LARRY library does not request ADT MEX"
else
    fail "Custom library should not request ADT MEX"
fi

if echo "$ROUTING" | grep -q 'id=adt_de .*featureRefType=Antibody Capture .*adt_mex=yes'; then
    pass "Protein feature_types maps to Antibody Capture + ADT MEX"
else
    fail "Protein library should map to Antibody Capture and ADT MEX"
fi

GRNA_OUT="$WORK_DIR/out/cr_assign/CRISPR_Guide_Capture/grna_de"
LARRY_OUT="$WORK_DIR/out/cr_assign/Custom/larry_de"
ADT_OUT="$WORK_DIR/out/cr_assign/Protein/adt_de"

for path in "$GRNA_OUT" "$LARRY_OUT" "$ADT_OUT"; do
    if echo "$ROUTING" | grep -q "assignOut=$path"; then
        pass "unique assignOut path $path"
    else
        fail "missing expected assignOut path $path"
    fi
done

# Shared whitelist
cat > "$WORK_DIR/whitelist.txt" <<'EOF'
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCT
EOF

make_fastqs() {
    local dir="$1"
    local feat_seq="$2"
    mkdir -p "$dir"
    python3 - <<PY "$dir" "$feat_seq"
import os, sys
out_dir, feat = sys.argv[1], sys.argv[2]
barcodes = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
]
reads = [
    (0, feat, "AAAAAAAAAAAA"),
    (1, feat, "CCCCCCCCCCCC"),
]
qual = lambda n: "I" * n
with open(os.path.join(out_dir, "sample_R1_001.fastq"), "w") as r1, \
     open(os.path.join(out_dir, "sample_R2_001.fastq"), "w") as r2:
    for i, (bc_idx, seq, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc)+len(umi))}\n")
        r2.write(f"@read{i}\n{seq}\n+\n{qual(len(seq))}\n")
PY
}

cat > "$WORK_DIR/ref_grna.csv" <<'EOF'
id,name,sequence,feature_type
g1,Guide-A,ATCGATCG,CRISPR Guide Capture
EOF

cat > "$WORK_DIR/ref_larry.csv" <<'EOF'
id,name,sequence,feature_type
lb1,LARRY-1,TTAATTAATTAATTAA,Custom
EOF

cat > "$WORK_DIR/ref_adt.csv" <<'EOF'
id,name,sequence,feature_type
CD29,CD29,ATCGATCGATCGATCG,ADT
EOF

make_fastqs "$WORK_DIR/fastq_grna" "ATCGATCG"
make_fastqs "$WORK_DIR/fastq_larry" "TTAATTAATTAATTAA"
make_fastqs "$WORK_DIR/fastq_adt" "ATCGATCGATCGATCG"

mkdir -p "$GRNA_OUT" "$LARRY_OUT" "$ADT_OUT"

"$ASSIGN" \
    --whitelist "$WORK_DIR/whitelist.txt" \
    --featurelist "$WORK_DIR/ref_grna.csv" \
    --directory "$GRNA_OUT" \
    --skip_empty_drops \
    "$WORK_DIR/fastq_grna" -b 16 -u 12

"$ASSIGN" \
    --whitelist "$WORK_DIR/whitelist.txt" \
    --featurelist "$WORK_DIR/ref_larry.csv" \
    --directory "$LARRY_OUT" \
    --skip_empty_drops \
    "$WORK_DIR/fastq_larry" -b 16 -u 12

"$ASSIGN" \
    --whitelist "$WORK_DIR/whitelist.txt" \
    --featurelist "$WORK_DIR/ref_adt.csv" \
    --directory "$ADT_OUT" \
    --output-mode adt_mex \
    --skip_empty_drops \
    "$WORK_DIR/fastq_adt" -b 16 -u 12

find_sample_dir() {
    local root="$1"
    if [[ -f "$root/barcodes.tsv.gz" || -f "$root/matrix.mtx" ]]; then
        echo "$root"
        return
    fi
    for sub in "$root"/*; do
        if [[ -d "$sub" ]]; then
            find_sample_dir "$sub"
            return
        fi
    done
}

ADT_SAMPLE="$(find_sample_dir "$ADT_OUT")"
if [[ -f "$ADT_SAMPLE/protein_quant_summary.json" && -f "$ADT_SAMPLE/features.tsv.gz" ]]; then
    pass "ADT assign output has protein MEX sidecars"
else
    fail "ADT assign output missing protein MEX sidecars"
fi

if zgrep -q "Antibody Capture" "$ADT_SAMPLE/features.tsv.gz"; then
    pass "ADT features.tsv.gz uses Antibody Capture"
else
    fail "ADT features.tsv.gz missing Antibody Capture"
fi

GRNA_SAMPLE="$(find_sample_dir "$GRNA_OUT")"
LARRY_SAMPLE="$(find_sample_dir "$LARRY_OUT")"

if [[ -f "$GRNA_SAMPLE/matrix.mtx" ]]; then
    pass "gRNA legacy matrix.mtx present"
else
    fail "gRNA matrix.mtx missing"
fi

if [[ -f "$LARRY_SAMPLE/matrix.mtx" ]]; then
    pass "Custom/LARRY legacy matrix.mtx present"
else
    fail "LARRY matrix.mtx missing"
fi

if [[ ! -f "$GRNA_SAMPLE/protein_quant_summary.json" && ! -f "$LARRY_SAMPLE/protein_quant_summary.json" ]]; then
    pass "gRNA and LARRY outputs omit protein MEX sidecars"
else
    fail "non-ADT libraries should not emit protein_quant_summary.json"
fi

echo ""
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
if [[ "$FAIL_COUNT" -ne 0 ]]; then
    exit 1
fi
