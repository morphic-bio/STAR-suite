#!/usr/bin/env bash

set -euo pipefail

TARGET_DIR="${1:-}"
if [[ -z "${TARGET_DIR}" ]]; then
  echo "Usage: $0 <target-dir>" >&2
  exit 1
fi

if [[ ! -d "${TARGET_DIR}" ]]; then
  echo "ERROR: directory not found: ${TARGET_DIR}" >&2
  exit 1
fi

checksum_file="${TARGET_DIR}/SHA256SUMS"
tmp_file="$(mktemp)"
trap 'rm -f "${tmp_file}"' EXIT

(
  cd "${TARGET_DIR}"
  find . -type f ! -name "SHA256SUMS" -print0 \
    | sort -z \
    | xargs -0 sha256sum
) > "${tmp_file}"

mv "${tmp_file}" "${checksum_file}"
echo "Created ${checksum_file}"
