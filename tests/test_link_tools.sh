#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

mkdir -p "$TMP_DIR/targets" "$TMP_DIR/links"

echo "hello" > "$TMP_DIR/targets/one.txt"

echo -e "links/one.txt\ttargets/one.txt" > "$TMP_DIR/manifest.tsv"

echo "Running sync_links..."
"$ROOT/scripts/sync_links.sh" --manifest "$TMP_DIR/manifest.tsv" --root "$TMP_DIR"

if [[ ! -L "$TMP_DIR/links/one.txt" ]]; then
  echo "FAIL: expected symlink not found"
  exit 1
fi

resolved=$(python3 - <<PY
import os
print(os.path.realpath("$TMP_DIR/links/one.txt"))
PY
)
if [[ "$resolved" != "$TMP_DIR/targets/one.txt" ]]; then
  echo "FAIL: symlink target mismatch"
  exit 1
fi

echo "Running check_links (should pass)..."
"$ROOT/scripts/check_links.sh" --manifest "$TMP_DIR/manifest.tsv" --root "$TMP_DIR"

# Negative test: break target
rm -f "$TMP_DIR/targets/one.txt"

echo "Running check_links (should fail)..."
if "$ROOT/scripts/check_links.sh" --manifest "$TMP_DIR/manifest.tsv" --root "$TMP_DIR"; then
  echo "FAIL: expected check_links to fail on missing target"
  exit 1
fi

echo "PASS"
