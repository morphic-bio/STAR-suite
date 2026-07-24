#!/usr/bin/env bash
set -euo pipefail

# bqtools is a test/fixture oracle. STAR Suite production CBQ input is native
# C++; the Rust tool must not become a runtime dependency of STAR.
BQTOOLS_VERSION="${BQTOOLS_VERSION:-0.5.6}"
RUST_IMAGE="${RUST_IMAGE:-rust:1.88.0-bookworm@sha256:af306cfa71d987911a781c37b59d7d67d934f49684058f96cf72079c3626bfe0}"
INSTALL_ROOT="${BQTOOLS_INSTALL_ROOT:-${HOME}/.local}"
BQTOOLS_BIN="${INSTALL_ROOT}/bin/bqtools"
PROVENANCE="${INSTALL_ROOT}/share/star-suite/bqtools-oracle-${BQTOOLS_VERSION}.tsv"

version_ok() {
    [[ -x "${BQTOOLS_BIN}" ]] || return 1
    "${BQTOOLS_BIN}" --version 2>/dev/null | grep -Fq "${BQTOOLS_VERSION}"
}

if version_ok; then
    echo "bqtools ${BQTOOLS_VERSION} is already installed: ${BQTOOLS_BIN}"
    echo "export BQTOOLS=${BQTOOLS_BIN}"
    exit 0
fi

command -v docker >/dev/null 2>&1 || {
    echo "ERROR: Docker is required to provision the pinned bqtools oracle" >&2
    exit 1
}
if [[ -e "${BQTOOLS_BIN}" ]]; then
    echo "ERROR: refusing to replace an unexpected bqtools binary: ${BQTOOLS_BIN}" >&2
    exit 1
fi

mkdir -p "${INSTALL_ROOT}/bin" "$(dirname "${PROVENANCE}")"
docker run --rm \
    --env CARGO_HOME=/tmp/cargo \
    --env BQTOOLS_VERSION="${BQTOOLS_VERSION}" \
    --env INSTALL_UID="$(id -u)" \
    --env INSTALL_GID="$(id -g)" \
    --volume "${INSTALL_ROOT}:/install" \
    "${RUST_IMAGE}" \
    bash -c '
        apt-get update &&
        apt-get install --no-install-recommends -y cmake clang libclang-dev &&
        cargo install --locked --version "${BQTOOLS_VERSION}" --root /install bqtools &&
        chown -R "${INSTALL_UID}:${INSTALL_GID}" /install
    '

version_ok || {
    echo "ERROR: installed bqtools does not report expected version ${BQTOOLS_VERSION}" >&2
    exit 1
}

{
    printf 'tool\tbqtools\n'
    printf 'version\t%s\n' "${BQTOOLS_VERSION}"
    printf 'binary\t%s\n' "${BQTOOLS_BIN}"
    printf 'binary_sha256\t%s\n' "$(sha256sum "${BQTOOLS_BIN}" | awk '{print $1}')"
    printf 'builder_image\t%s\n' "${RUST_IMAGE}"
    printf 'installed_utc\t%s\n' "$(date -u +%FT%TZ)"
} > "${PROVENANCE}"

echo "Installed bqtools ${BQTOOLS_VERSION}: ${BQTOOLS_BIN}"
echo "Provenance: ${PROVENANCE}"
echo "export BQTOOLS=${BQTOOLS_BIN}"
