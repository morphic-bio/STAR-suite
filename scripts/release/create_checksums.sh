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
names_file="$(mktemp)"
trap 'rm -f "${tmp_file}"' EXIT
trap 'rm -f "${tmp_file}" "${names_file}"' EXIT

(
  cd "${TARGET_DIR}"
  find . -type f ! -name "SHA256SUMS" -printf '%P\n' \
    | LC_ALL=C sort > "${names_file}"

  awk -F/ '
    {
      base = $NF
      seen[base]++
    }
    END {
      for (base in seen) {
        if (seen[base] > 1) {
          printf("ERROR: duplicate release asset basename: %s\n", base) > "/dev/stderr"
          dup = 1
        }
      }
      exit dup
    }
  ' "${names_file}"

  while IFS= read -r relpath; do
    hash="$(sha256sum "${relpath}" | awk '{print $1}')"
    printf '%s  %s\n' "${hash}" "$(basename "${relpath}")"
  done < "${names_file}"
) > "${tmp_file}"

mv "${tmp_file}" "${checksum_file}"
echo "Created ${checksum_file}"
