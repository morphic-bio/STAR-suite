#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

IN_DIR="${REPO_ROOT}/dist/release/source"
CHANGES_FILE=""
KEY_ID="${DEB_SIGNING_KEY_ID:-}"
PASSPHRASE_ENV="DEB_SIGNING_KEY_PASSPHRASE"
VERIFY=1

usage() {
  cat <<EOF
Usage: $0 [--in-dir <dir>] [--changes <file>] --key-id <gpg-key-id> [--passphrase-env <var>] [--no-verify]

Sign Debian source package artifacts (.changes and .dsc) using debsign.

Notes:
  - If --changes is omitted, the latest *_source.changes under --in-dir is used.
  - If --key-id is omitted, DEB_SIGNING_KEY_ID is used.
  - Passphrase is read from the environment variable named by --passphrase-env
    (default: DEB_SIGNING_KEY_PASSPHRASE). If empty, standard gpg agent flow is used.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-dir)
      IN_DIR="$2"
      shift 2
      ;;
    --changes)
      CHANGES_FILE="$2"
      shift 2
      ;;
    --key-id)
      KEY_ID="$2"
      shift 2
      ;;
    --passphrase-env)
      PASSPHRASE_ENV="$2"
      shift 2
      ;;
    --no-verify)
      VERIFY=0
      shift 1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if ! command -v debsign >/dev/null 2>&1; then
  echo "ERROR: debsign is required (install devscripts)." >&2
  exit 1
fi

if [[ -z "${KEY_ID}" ]]; then
  echo "ERROR: --key-id is required (or set DEB_SIGNING_KEY_ID)." >&2
  exit 1
fi

if [[ -z "${CHANGES_FILE}" ]]; then
  CHANGES_FILE="$(ls -1t "${IN_DIR}"/*_source.changes 2>/dev/null | head -n1 || true)"
fi

if [[ -z "${CHANGES_FILE}" ]]; then
  echo "ERROR: no *_source.changes file found in ${IN_DIR}" >&2
  exit 1
fi

if [[ ! -f "${CHANGES_FILE}" ]]; then
  echo "ERROR: changes file not found: ${CHANGES_FILE}" >&2
  exit 1
fi

CHANGES_FILE="$(realpath "${CHANGES_FILE}")"
SOURCE_DSC="${CHANGES_FILE%_source.changes}.dsc"

PASSPHRASE="${!PASSPHRASE_ENV:-}"

if [[ -n "${PASSPHRASE}" ]]; then
  export DEBSIGN_GPG_PASSPHRASE="${PASSPHRASE}"
  GPG_WRAPPER="$(mktemp)"
  trap 'rm -f "${GPG_WRAPPER}"' EXIT
  cat > "${GPG_WRAPPER}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
exec gpg --batch --yes --pinentry-mode loopback --passphrase "${DEBSIGN_GPG_PASSPHRASE}" "$@"
EOF
  chmod 700 "${GPG_WRAPPER}"
  debsign -k"${KEY_ID}" -p"${GPG_WRAPPER}" "${CHANGES_FILE}"
else
  debsign -k"${KEY_ID}" "${CHANGES_FILE}"
fi

if [[ "${VERIFY}" -eq 1 ]]; then
  gpg --batch --verify "${CHANGES_FILE}" >/dev/null
  if [[ -f "${SOURCE_DSC}" ]]; then
    gpg --batch --verify "${SOURCE_DSC}" >/dev/null
  fi
fi

echo "Signed source artifacts:"
echo "  ${CHANGES_FILE}"
if [[ -f "${SOURCE_DSC}" ]]; then
  echo "  ${SOURCE_DSC}"
fi
