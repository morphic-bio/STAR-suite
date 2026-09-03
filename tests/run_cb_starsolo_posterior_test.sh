#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=$(mktemp -d)
cleanup() {
    rm -f "$build_dir/test_cb_starsolo_posterior"
    rmdir "$build_dir"
}
trap cleanup EXIT

"${CXX:-g++}" -std=c++11 -O2 -Wall -Wextra -Werror \
    -I"$repo_root/flex/source" \
    "$repo_root/tests/test_cb_starsolo_posterior.cpp" \
    "$repo_root/flex/source/solo/CbBayesianResolver.cpp" \
    -o "$build_dir/test_cb_starsolo_posterior"

"$build_dir/test_cb_starsolo_posterior"
