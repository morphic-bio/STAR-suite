#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PF="$ROOT/core/features/process_features"
OUT="$(mktemp -d /tmp/hash_pair_test.XXXXXX)"
trap 'rm -rf "$OUT"' EXIT
cc -std=gnu11 -O2 -I"$PF/include" -c "$ROOT/tests/multi_feature/hash_pair_harness.c" -o "$OUT/test.o"
g++ -o "$OUT/test" "$OUT/test.o" "$PF/libprocess_features.a" "$ROOT/core/features/libscrna/libscrna.a" -lm -lpthread -lz -fopenmp -lhts
printf 'sample\thash_b\thash_a\textra\nplate01\tB\tA\tignored\n' > "$OUT/pairs.tsv"
mkdir "$OUT/pair" "$OUT/ratio"
"$OUT/test" "$OUT/pair" pair "$OUT/pairs.tsv"
"$OUT/test" "$OUT/ratio" ratio "$OUT/pairs.tsv"
python3 - "$OUT" <<'PY'
import csv,json,sys
from pathlib import Path
out=Path(sys.argv[1])
rows=list(csv.DictReader((out/'pair/hash_demux_assignments.tsv').open(),delimiter='\t'))
assert [r['hash_classification'] for r in rows]==['singlet','singlet','multiplet','negative','unknown_pair','singlet']
assert rows[0]['hash_top_feature']=='B' and rows[0]['hash_pair']=='A|B'
assert rows[1]['hash_third_count']=='4' and rows[1]['hash_pair_ratio']=='2.500000'
summary=json.loads((out/'pair/hash_demux_summary.json').read_text())
assert summary['per_sample']=={'plate01':3} and summary['n_unknown_pair']==1
header=(out/'ratio/hash_demux_assignments.tsv').read_text().splitlines()[0]
assert header=='barcode\thash_assignment\thash_classification\thash_total_umis\thash_top_feature\thash_top_count\thash_second_feature\thash_second_count\thash_top_ratio'
PY
# Invalid table rows must fail (unordered duplicates, unknown tags, self pairs).
for row in $'second\tA\tB' $'second\tA\tZ' $'second\tA\tA'; do
  printf 'sample\thash_a\thash_b\nfirst\tA\tB\n%s\n' "$row" > "$OUT/bad.tsv"
  if "$OUT/test" "$OUT/pair" pair "$OUT/bad.tsv"; then echo 'Invalid pair accepted' >&2; exit 1; fi
done
# Exercise both supported pf-multi header blocks and relative path resolution.
cat > "$OUT/config.cpp" <<'CPP'
#include "PfMultiConfig.h"
#include <cassert>
int main(int argc, char **argv) {
    auto cfg=PfMultiConfig::parseConfig(argv[1]);
    for(const auto& lib:cfg.libraries) {
        assert(lib.starHashDemuxMethod=="pair");
        assert(lib.starHashMinPairRatio==2.5);
        assert(lib.starHashSampleTable.find("/pairs.tsv")!=std::string::npos);
    }
}
CPP
g++ -std=c++11 -I"$ROOT/core/legacy/source" -I"$PF/include" "$OUT/config.cpp" "$ROOT/core/legacy/source/PfMultiConfig.o" -o "$OUT/config_test"
printf '[libraries]\nfastqs,sample,library_type,star_hash_demux_method,star_hash_sample_table,star_hash_min_pair_ratio\n/tmp/fastqs,S1,Multiplexing Capture,pair,pairs.tsv,2.5\n' > "$OUT/config.csv"
"$OUT/config_test" "$OUT/config.csv"
tail -n +2 "$OUT/config.csv" > "$OUT/config_simple.csv"
"$OUT/config_test" "$OUT/config_simple.csv"
echo 'PASS: native pair calls, malformed tables, legacy header, pf-multi columns'
