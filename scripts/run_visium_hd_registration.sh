#!/usr/bin/env bash
# Register the microscope image of a Visium HD slide to its capture grid, the
# second stage of STARgrid. This is a thin wrapper: the registration itself is
# scripts/run_registration.sh of BioDepot/stargrid. Every option except
# --stargrid is passed to it unchanged.
set -euo pipefail

usage() {
    cat <<'EOF2'
usage: run_visium_hd_registration.sh [--stargrid DIR] --he FILE --cytassist FILE
                                     --layout FILE --area NAME --counts DIR
                                     --out-dir DIR [options]

  --stargrid DIR   checkout of https://github.com/BioDepot/stargrid (default:
                   the environment variable STARGRID_REPO)

All other options are those of DIR/scripts/run_registration.sh and are passed to
it unchanged; run this wrapper with --stargrid DIR --help to list them.
EOF2
}

die() {
    echo "error: $*" >&2
    exit 2
}

stargrid=${STARGRID_REPO:-}
help=0
args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stargrid) [[ $# -ge 2 ]] || die "--stargrid needs a folder"; stargrid=$2; shift 2 ;;
        -h|--help) help=1; args+=("$1"); shift ;;
        *) args+=("$1"); shift ;;
    esac
done

if [[ $help -eq 1 ]]; then
    usage
    if [[ -n "$stargrid" && -f "$stargrid/scripts/run_registration.sh" ]]; then
        echo
        "$stargrid/scripts/run_registration.sh" --help
    fi
    exit 0
fi

[[ -n "$stargrid" ]] || die "give --stargrid DIR or set STARGRID_REPO: a checkout of https://github.com/BioDepot/stargrid"
[[ -d "$stargrid" ]] || die "--stargrid $stargrid is not a folder"
script=$stargrid/scripts/run_registration.sh
[[ -f "$script" ]] || die "$script not found; $stargrid is not a checkout of BioDepot/stargrid"

exec bash "$script" "${args[@]}"
