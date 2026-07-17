#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
make_dir="$repo_root/core/legacy/source"

amd64_config=$(make -s -C "$make_dir" print-simd-config \
  CXX_TARGET_MACHINE=x86_64-linux-gnu)
arm64_config=$(make -s -C "$make_dir" print-simd-config \
  CXX_TARGET_MACHINE=aarch64-linux-gnu)

grep -qx 'target_machine=x86_64-linux-gnu' <<<"$amd64_config"
grep -qx 'simd_flags=-mavx2' <<<"$amd64_config"
grep -qx 'target_machine=aarch64-linux-gnu' <<<"$arm64_config"
grep -qx 'simd_flags=' <<<"$arm64_config"

printf 'architecture-aware SIMD flags: PASS\n'
