#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/dist/release/deb"
ALLOW_MISSING_DEPS=0
DOCKER_IMAGE=""
DISTRO_SUFFIX=""
NATIVE_BUILD=0
SNAPSHOT_DIR=""

cleanup() {
  if [[ -n "${SNAPSHOT_DIR}" && -d "${SNAPSHOT_DIR}" ]]; then
    rm -rf "${SNAPSHOT_DIR}"
  fi
}
trap cleanup EXIT

usage() {
  cat <<EOF
Usage: $0 [--out-dir <dir>] [--allow-missing-build-deps]
          [--docker-image <image>] [--distro-suffix <suffix>]

Build Debian binary package(s) from a clean HEAD snapshot and copy artifacts to out-dir.

Options:
  --docker-image <image>    build inside an existing release-build image
  --distro-suffix <suffix> append ~<suffix> to the Debian version
                            (for example ubuntu22.04.1)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --allow-missing-build-deps)
      ALLOW_MISSING_DEPS=1
      shift 1
      ;;
    --docker-image)
      DOCKER_IMAGE="$2"
      shift 2
      ;;
    --distro-suffix)
      DISTRO_SUFFIX="$2"
      shift 2
      ;;
    --native-build)
      NATIVE_BUILD=1
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

if [[ -n "${DISTRO_SUFFIX}" && ! "${DISTRO_SUFFIX}" =~ ^[A-Za-z0-9.+]+$ ]]; then
  echo "ERROR: --distro-suffix may contain only letters, digits, '.', and '+'" >&2
  exit 1
fi

resolve_commit_sha() {
  if git -C "${REPO_ROOT}" rev-parse HEAD >/dev/null 2>&1; then
    git -C "${REPO_ROOT}" rev-parse HEAD
  elif [[ -n "${STAR_SUITE_COMMIT_SHA:-}" ]]; then
    printf '%s\n' "${STAR_SUITE_COMMIT_SHA}"
  else
    return 1
  fi
}

require_commit_sha() {
  local commit_sha
  if ! commit_sha="$(resolve_commit_sha)" \
      || [[ ! "${commit_sha}" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "ERROR: Debian release builds require an exact 40-character source commit SHA" >&2
    exit 1
  fi
  printf '%s\n' "${commit_sha,,}"
}

run_container_build() {
  if ! command -v docker >/dev/null 2>&1; then
    echo "ERROR: docker is required for --docker-image builds" >&2
    exit 1
  fi
  if ! git -C "${REPO_ROOT}" rev-parse --verify HEAD >/dev/null 2>&1; then
    echo "ERROR: --docker-image mode requires a git checkout with a committed HEAD" >&2
    exit 1
  fi

  local out_dir_abs commit_sha suffix_arg
  SNAPSHOT_DIR="$(mktemp -d)"
  mkdir -p "${OUT_DIR}"
  out_dir_abs="$(cd "${OUT_DIR}" && pwd)"
  commit_sha="$(require_commit_sha)"
  suffix_arg=""
  if [[ -n "${DISTRO_SUFFIX}" ]]; then
    suffix_arg="--distro-suffix '${DISTRO_SUFFIX}'"
  fi

  git -C "${REPO_ROOT}" archive --format=tar HEAD | tar -xf - -C "${SNAPSHOT_DIR}"

  docker run --rm \
    -v "${SNAPSHOT_DIR}:/src-ro:ro" \
    -v "${out_dir_abs}:/out" \
    "${DOCKER_IMAGE}" \
    bash -lc "
      set -euo pipefail
      export DEBIAN_FRONTEND=noninteractive
      if [[ \"\${STAR_SUITE_DEB_BUILD_CONTAINER_READY:-0}\" != \"1\" ]]; then
        apt-get update
        apt-get install -y --no-install-recommends \
          build-essential xxd debhelper devscripts dpkg-dev fakeroot \
          pkg-config zlib1g-dev libbz2-dev liblzma-dev \
          libcurl4-gnutls-dev libssl-dev libglib2.0-dev libhts-dev \
          ca-certificates git
      fi
      mkdir -p /src
      cp -a /src-ro/. /src/
      cd /src
      export STAR_SUITE_COMMIT_SHA='${commit_sha}'
      scripts/release/build_deb_binary_package.sh \
        --native-build \
        --out-dir /out \
        ${suffix_arg}
    "
}

if [[ -n "${DOCKER_IMAGE}" && "${NATIVE_BUILD}" -eq 0 ]]; then
  run_container_build
  exit 0
fi

OUT_DIR="$(cd "${REPO_ROOT}" && mkdir -p "${OUT_DIR}" && cd "${OUT_DIR}" && pwd)"

cd "${REPO_ROOT}"

if [[ ! -f debian/control ]]; then
  echo "ERROR: debian/control not found. Debian packaging skeleton is required." >&2
  exit 1
fi

if ! command -v dh >/dev/null 2>&1; then
  cat >&2 <<EOF
ERROR: debhelper (dh) is required but not installed.
Install with: sudo apt-get install -y xxd debhelper devscripts dpkg-dev fakeroot
EOF
  exit 1
fi

build_args=(-b -us -uc)
if [[ "${ALLOW_MISSING_DEPS}" -eq 1 ]]; then
  build_args+=(-d)
fi

source_name="$(dpkg-parsechangelog -S Source)"
version="$(dpkg-parsechangelog -S Version)"
arch="$(dpkg --print-architecture)"
package_version="${version}"
if [[ -n "${DISTRO_SUFFIX}" ]]; then
  package_version+="~${DISTRO_SUFFIX}"
fi
base="${source_name}_${package_version}_${arch}"

tmp_root="$(mktemp -d)"
trap 'rm -rf "${tmp_root}"' EXIT
work_src="${tmp_root}/${source_name}-${version}"
mkdir -p "${work_src}"

echo "Preparing clean source snapshot..."
if [[ "${NATIVE_BUILD}" -eq 1 ]]; then
  tar \
    --exclude='./.git' \
    --exclude='./dist' \
    -cf - . | tar -xf - -C "${work_src}"
else
  git archive --format=tar HEAD | tar -xf - -C "${work_src}"
fi
commit_sha="$(require_commit_sha)"
printf '%s\n' "${commit_sha}" > "${work_src}/debian/source-revision"

if [[ -n "${DISTRO_SUFFIX}" ]]; then
  escaped_version="${version//./\\.}"
  sed -i "1s/(${escaped_version})/(${package_version})/" \
    "${work_src}/debian/changelog"
  if [[ "$(dpkg-parsechangelog -l "${work_src}/debian/changelog" -S Version)" != "${package_version}" ]]; then
    echo "ERROR: failed to apply Debian distribution suffix ${DISTRO_SUFFIX}" >&2
    exit 1
  fi
fi

echo "Building Debian binary package..."
cd "${work_src}"
if ! dpkg-buildpackage "${build_args[@]}"; then
  if [[ "${ALLOW_MISSING_DEPS}" -eq 0 ]]; then
    cat >&2 <<EOF
ERROR: dpkg-buildpackage failed.
Hint: install build deps or retry with --allow-missing-build-deps for local smoke runs.
EOF
  fi
  exit 1
fi

mkdir -p "${OUT_DIR}"
for ext in deb changes buildinfo; do
  for f in "${tmp_root}/${base}.${ext}"; do
    [[ -f "${f}" ]] && cp -f "${f}" "${OUT_DIR}/"
  done
done

echo "Binary package artifacts copied to ${OUT_DIR}"
