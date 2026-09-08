#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

g++ -O2 -Wall -Wextra -std=c++11 \
    -I"$repo_root/flex/source" \
    "$repo_root/tests/flex_sample_tag_cache_test.cpp" \
    -o "$tmp_dir/flex_sample_tag_cache_test"
"$tmp_dir/flex_sample_tag_cache_test"

defaults="$repo_root/core/legacy/source/parametersDefault"
params_h="$repo_root/core/legacy/source/ParametersSolo.h"
params_cpp="$repo_root/core/legacy/source/Parameters.cpp"
pipeline="$repo_root/core/legacy/source/FlexPipeline.cpp"
detector="$repo_root/flex/source/SampleDetector.cpp"

grep -Eq '^soloSampleProbeOffset[[:space:]]+68$' "$defaults"
grep -Eq '^soloSampleSearchNearby[[:space:]]+no$' "$defaults"
grep -Eq '^soloSampleTagMismatch[[:space:]]+1$' "$defaults"
grep -q 'sampleProbeOffset = 68' "$params_h"
grep -q 'sampleSearchNearbyStr = "no"' "$params_h"
grep -q 'sampleSearchNearby = false' "$params_h"
grep -q 'setStringIfDefault("soloSampleSearchNearby", "no")' "$params_cpp"

if grep -Eq 'deltas.*-1|primary[[:space:]]*\+[[:space:]]*delta|primary[[:space:]]*-[[:space:]]*delta' \
    "$pipeline" "$detector"; then
    echo "FAIL: neighboring-offset sample-tag search remains in a runtime path" >&2
    exit 1
fi

echo "PASS: Flex sample-tag policy is fixed-offset H0 then H1"
