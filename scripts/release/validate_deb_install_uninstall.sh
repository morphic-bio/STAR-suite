#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

DEB_PATH=""
DOCKER_IMAGE="ubuntu:24.04"
EXPECTED_COMMIT=""

usage() {
  cat <<EOF
Usage: $0 [--deb <path>] [--image <docker-image>] [--expected-commit <sha>]

Validate install/uninstall of star-suite Debian package in a clean Ubuntu container.

Defaults:
  --deb   latest dist/release/deb/star-suite_*_amd64.deb
  --image ubuntu:24.04
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --deb)
      DEB_PATH="$2"
      shift 2
      ;;
    --image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --expected-commit)
      EXPECTED_COMMIT="$2"
      shift 2
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

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker is required for install/uninstall validation." >&2
  exit 1
fi

if [[ -z "${DEB_PATH}" ]]; then
  DEB_PATH="$(ls -1t "${REPO_ROOT}"/dist/release/deb/star-suite_*_amd64.deb 2>/dev/null | head -n1 || true)"
fi

if [[ -z "${DEB_PATH}" ]]; then
  echo "ERROR: no Debian package found. Pass --deb or build one first." >&2
  exit 1
fi

if [[ ! -f "${DEB_PATH}" ]]; then
  echo "ERROR: Debian package not found: ${DEB_PATH}" >&2
  exit 1
fi

DEB_PATH="$(realpath "${DEB_PATH}")"
TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

cp -f "${DEB_PATH}" "${TMPDIR}/star-suite.deb"

cat > "${TMPDIR}/validate.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y --no-install-recommends ca-certificates
apt-get install -y --no-install-recommends /work/star-suite.deb

for bin in STAR star_feature_call run_flexfilter_mex flexfilter molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize slam_requant pileup_snp remove_y_reads; do
  command -v "${bin}" >/dev/null
done

STAR --version
source_revision="$(STAR --source-revision)"
if [[ ! "${source_revision}" =~ ^[0-9a-f]{40}$ ]]; then
  echo "ERROR: installed Debian STAR has invalid source revision: ${source_revision}" >&2
  exit 1
fi
if [[ -n "${EXPECTED_COMMIT:-}" && "${source_revision}" != "${EXPECTED_COMMIT,,}" ]]; then
  echo "ERROR: installed Debian STAR source revision ${source_revision} does not match ${EXPECTED_COMMIT,,}" >&2
  exit 1
fi

dpkg -s star-suite >/dev/null
apt-get purge -y star-suite
apt-get autoremove -y --purge

hash -r

if dpkg-query -W -f='${Status}' star-suite 2>/dev/null | grep -q "install ok installed"; then
  echo "ERROR: package star-suite still installed after purge." >&2
  exit 1
fi

if [[ -x /usr/bin/STAR ]]; then
  echo "ERROR: /usr/bin/STAR still present after purge." >&2
  exit 1
fi

echo "Install/uninstall validation passed."
EOF

chmod +x "${TMPDIR}/validate.sh"

echo "Validating package: ${DEB_PATH}"
echo "Container image: ${DOCKER_IMAGE}"

docker run --rm \
  -e EXPECTED_COMMIT="${EXPECTED_COMMIT}" \
  -v "${TMPDIR}:/work:rw" \
  "${DOCKER_IMAGE}" \
  bash /work/validate.sh

echo "OK: Debian install/uninstall validation complete."
