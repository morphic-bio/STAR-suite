#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:-$(mktemp -d /tmp/flex-khash-test.XXXXXX)}"
mkdir -p "$out"
g++ -O2 -std=c++11 -Wall -Wextra -I"$repo/core/legacy/source" -I"$repo/flex/source" \
    "$repo/tests/flex_khash_storage_test.cpp" "$repo/core/legacy/source/FlexHashCacheStorage.cpp" \
    -o "$out/storage_test"
python3 - "$out" <<'PY'
from pathlib import Path
import struct,sys
root=Path(sys.argv[1])
# Zero and all-ones keys, both tiers, H2, sample duplicates, region bits,
# maximum gene/sample IDs, and ignored legacy-version fields.
records=[(0,0,1|(2<<30),0,0,1),(0,0,1|(2<<30),0,0,2),
         (1,0,32767|(1<<30),1,0,0),(2,0,0,2,1,0),
         (3,0,2|(3<<30),4,0,0),(4,0,9,3,0,0),
         ((1<<64)-1,(1<<36)-1,5,0,1,65535)]
for version in [1,2,3]:
    with (root/f'v{version}.bin').open('xb') as f:
        f.write(struct.pack('<8sHHIQ',b'FH01SEQ1',version,50,24,len(records)))
        for r in records:f.write(struct.pack('<QQIBBH',*r))
PY
for version in 1 2 3; do
    "$out/storage_test" "$out/v$version.bin" "$out/v$version.khash" "$out/corrupt$version.khash"
done
g++ -O2 -std=c++11 -I"$repo/flex/tools/hash_screen_replay" \
    -I"$repo/core/legacy/source" -I"$repo/flex/source" \
    "$repo/tests/flex_probe_pair_test.cpp" "$repo/core/legacy/source/FlexHashCacheStorage.cpp" "$repo/core/legacy/source/FlexProbePairKhash.cpp" \
    -o "$out/pair_test"
"$out/pair_test" "$out/pair_fixture.bin"
g++ -O2 -std=c++11 -I"$repo/core/legacy/source" -I"$repo/flex/source" \
    "$repo/tests/flex_probe_pair_storage_test.cpp" "$repo/core/legacy/source/FlexHashCacheStorage.cpp" \
    "$repo/core/legacy/source/FlexProbePairKhash.cpp" -o "$out/pair_storage_test"
"$out/pair_storage_test" "$out/pair_storage"
echo "PASS: Flex khash cache tests; artifacts: $out"
