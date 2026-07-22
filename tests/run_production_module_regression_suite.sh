#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${MANIFEST:-${ROOT_DIR}/tests/production_module_regression_manifest.tsv}"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_production_module_regression_${STAMP}}"

MODE="preflight"
FILTER_MODULE=""
FILTER_CASE=""
FILTER_SCALE=""
CONTINUE_ON_FAIL=0
STRICT_PREFLIGHT=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Production-module regression suite driver. Defaults to --preflight.

Modes:
  --list                  Print selected manifest rows.
  --preflight             Check scripts, tools, and fixture paths without running.
  --run                   Execute selected runnable cases.

Filters:
  --module NAME           Select one module from the manifest.
  --case NAME             Select one case from the manifest.
  --scale NAME            Select one scale, e.g. contract or production-100k.

Other:
  --manifest PATH         Override manifest path.
  --out-root DIR          Wrapper log root for --run.
  --continue-on-fail      Keep running remaining cases after a failure.
  --strict-preflight      Make --preflight exit nonzero when selected cases miss requirements.
  -h, --help              Show this help.
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

trim() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "$value"
}

resolve_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        printf '%s' "$path"
    else
        printf '%s/%s' "$ROOT_DIR" "$path"
    fi
}

requirement_ok() {
    local req
    req="$(trim "$1")"
    [[ -z "$req" || "$req" == "-" ]] && return 0

    if [[ "$req" == any:* ]]; then
        local alternatives="${req#any:}"
        local alt
        IFS='|' read -r -a _alts <<< "$alternatives"
        for alt in "${_alts[@]}"; do
            if requirement_ok "$alt"; then
                return 0
            fi
        done
        return 1
    fi

    case "$req" in
        cmd:*)
            command -v "${req#cmd:}" >/dev/null 2>&1
            ;;
        env:*)
            local var="${req#env:}"
            [[ -n "${!var:-}" ]]
            ;;
        file:*)
            [[ -f "$(resolve_path "${req#file:}")" ]]
            ;;
        dir:*)
            [[ -d "$(resolve_path "${req#dir:}")" ]]
            ;;
        path:*)
            [[ -e "$(resolve_path "${req#path:}")" ]]
            ;;
        *)
            [[ -e "$(resolve_path "$req")" ]]
            ;;
    esac
}

script_path_ok() {
    local script="$1"
    [[ -f "$(resolve_path "$script")" ]]
}

selected_case() {
    local module="$1"
    local case_name="$2"
    local scale="$3"
    [[ -z "$FILTER_MODULE" || "$module" == "$FILTER_MODULE" ]] || return 1
    [[ -z "$FILTER_CASE" || "$case_name" == "$FILTER_CASE" ]] || return 1
    [[ -z "$FILTER_SCALE" || "$scale" == "$FILTER_SCALE" ]] || return 1
}

missing_requirements() {
    local script="$1"
    local reqs="$2"
    local missing=()

    if ! script_path_ok "$script"; then
        missing+=("script:${script}")
    fi

    if [[ -n "$reqs" && "$reqs" != "-" ]]; then
        local req
        IFS=';' read -r -a _reqs <<< "$reqs"
        for req in "${_reqs[@]}"; do
            req="$(trim "$req")"
            [[ -z "$req" || "$req" == "-" ]] && continue
            if ! requirement_ok "$req"; then
                missing+=("$req")
            fi
        done
    fi

    if [[ "${#missing[@]}" -gt 0 ]]; then
        printf '%s\n' "${missing[@]}"
        return 1
    fi
    return 0
}

path_value_available() {
    local value="$1"
    if [[ "$value" == /* || "$value" == ./* || "$value" == ../* || "$value" == *"/"* ]]; then
        [[ -e "$(resolve_path "$value")" ]]
    else
        return 0
    fi
}

env_tokens_for_case() {
    local env_spec="$1"
    ENV_TOKENS=()
    [[ -z "$env_spec" || "$env_spec" == "-" ]] && return 0

    local entry var value
    IFS=';' read -r -a _env_entries <<< "$env_spec"
    for entry in "${_env_entries[@]}"; do
        entry="$(trim "$entry")"
        [[ -z "$entry" || "$entry" == "-" ]] && continue
        if [[ "$entry" == *"?="* ]]; then
            var="${entry%%\?=*}"
            value="${entry#*\?=}"
            if [[ -z "${!var:-}" ]] && path_value_available "$value"; then
                ENV_TOKENS+=("${var}=${value}")
            fi
        else
            ENV_TOKENS+=("$entry")
        fi
    done
}

safe_name() {
    printf '%s' "$1" | tr -c '[:alnum:]_.-' '_'
}

print_header() {
    printf '%-18s %-34s %-15s %s\n' "module" "case" "scale" "script"
}

print_case() {
    local module="$1"
    local case_name="$2"
    local scale="$3"
    local script="$4"
    printf '%-18s %-34s %-15s %s\n' "$module" "$case_name" "$scale" "$script"
}

run_selected_cases() {
    local selected=0
    local runnable=0
    local skipped=0
    local failed=0

    mkdir -p "$OUT_ROOT"
    printf 'mode=%s\nmanifest=%s\nout_root=%s\n' "$MODE" "$MANIFEST" "$OUT_ROOT" > "${OUT_ROOT}/suite.env"

    while IFS=$'\t' read -r module case_name scale script reqs env_spec description; do
        [[ -z "${module:-}" || "$module" == \#* || "$module" == "module" ]] && continue
        selected_case "$module" "$case_name" "$scale" || continue
        selected=$((selected + 1))

        mapfile -t missing < <(missing_requirements "$script" "$reqs" || true)
        if [[ "${#missing[@]}" -gt 0 ]]; then
            skipped=$((skipped + 1))
            printf 'SKIP\t%s\t%s\tmissing: %s\n' "$module" "$case_name" "$(IFS=,; printf '%s' "${missing[*]}")" | tee -a "${OUT_ROOT}/summary.tsv"
            continue
        fi

        runnable=$((runnable + 1))
        env_tokens_for_case "$env_spec"
        local case_dir="${OUT_ROOT}/$(safe_name "${module}_${case_name}")"
        mkdir -p "$case_dir"

        {
            printf 'module=%s\ncase=%s\nscale=%s\nscript=%s\n' "$module" "$case_name" "$scale" "$(resolve_path "$script")"
            printf 'description=%s\n' "$description"
            printf 'env_tokens='
            printf '%q ' "${ENV_TOKENS[@]}"
            printf '\n'
        } > "${case_dir}/case.env"

        printf 'RUN\t%s\t%s\n' "$module" "$case_name" | tee -a "${OUT_ROOT}/summary.tsv"
        set +e
        (
            cd "$ROOT_DIR"
            env OUT_ROOT="${case_dir}/work" "${ENV_TOKENS[@]}" bash "$(resolve_path "$script")"
        ) > "${case_dir}/stdout.log" 2> "${case_dir}/stderr.log"
        local status=$?
        set -e

        printf '%s\n' "$status" > "${case_dir}/exit_status"
        if [[ "$status" -eq 0 ]]; then
            printf 'PASS\t%s\t%s\n' "$module" "$case_name" | tee -a "${OUT_ROOT}/summary.tsv"
        else
            failed=$((failed + 1))
            printf 'FAIL\t%s\t%s\texit=%s\n' "$module" "$case_name" "$status" | tee -a "${OUT_ROOT}/summary.tsv"
            if [[ "$CONTINUE_ON_FAIL" != "1" ]]; then
                break
            fi
        fi
    done < "$MANIFEST"

    [[ "$selected" -gt 0 ]] || die "No cases selected"
    printf 'Selected=%s Runnable=%s Skipped=%s Failed=%s\n' "$selected" "$runnable" "$skipped" "$failed" | tee -a "${OUT_ROOT}/summary.tsv"
    [[ "$failed" -eq 0 ]]
}

preflight_selected_cases() {
    local selected=0
    local runnable=0
    local missing_cases=0

    print_header
    while IFS=$'\t' read -r module case_name scale script reqs env_spec description; do
        [[ -z "${module:-}" || "$module" == \#* || "$module" == "module" ]] && continue
        selected_case "$module" "$case_name" "$scale" || continue
        selected=$((selected + 1))
        print_case "$module" "$case_name" "$scale" "$script"
        mapfile -t missing < <(missing_requirements "$script" "$reqs" || true)
        if [[ "${#missing[@]}" -eq 0 ]]; then
            runnable=$((runnable + 1))
            printf '  OK\n'
        else
            missing_cases=$((missing_cases + 1))
            printf '  MISSING: %s\n' "$(IFS=,; printf '%s' "${missing[*]}")"
        fi
    done < "$MANIFEST"

    [[ "$selected" -gt 0 ]] || die "No cases selected"
    printf '\nSelected=%s Runnable=%s Missing=%s\n' "$selected" "$runnable" "$missing_cases"
    if [[ "$STRICT_PREFLIGHT" == "1" && "$missing_cases" -gt 0 ]]; then
        return 1
    fi
    return 0
}

list_selected_cases() {
    local selected=0
    print_header
    while IFS=$'\t' read -r module case_name scale script reqs env_spec description; do
        [[ -z "${module:-}" || "$module" == \#* || "$module" == "module" ]] && continue
        selected_case "$module" "$case_name" "$scale" || continue
        selected=$((selected + 1))
        print_case "$module" "$case_name" "$scale" "$script"
    done < "$MANIFEST"
    [[ "$selected" -gt 0 ]] || die "No cases selected"
}

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --list)
            MODE="list"
            shift
            ;;
        --preflight)
            MODE="preflight"
            shift
            ;;
        --run)
            MODE="run"
            shift
            ;;
        --module)
            [[ "$#" -ge 2 ]] || die "--module requires a value"
            FILTER_MODULE="$2"
            shift 2
            ;;
        --case)
            [[ "$#" -ge 2 ]] || die "--case requires a value"
            FILTER_CASE="$2"
            shift 2
            ;;
        --scale)
            [[ "$#" -ge 2 ]] || die "--scale requires a value"
            FILTER_SCALE="$2"
            shift 2
            ;;
        --manifest)
            [[ "$#" -ge 2 ]] || die "--manifest requires a value"
            MANIFEST="$2"
            shift 2
            ;;
        --out-root)
            [[ "$#" -ge 2 ]] || die "--out-root requires a value"
            OUT_ROOT="$2"
            shift 2
            ;;
        --continue-on-fail)
            CONTINUE_ON_FAIL=1
            shift
            ;;
        --strict-preflight)
            STRICT_PREFLIGHT=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -f "$MANIFEST" ]] || die "Missing manifest: $MANIFEST"

case "$MODE" in
    list)
        list_selected_cases
        ;;
    preflight)
        preflight_selected_cases
        ;;
    run)
        run_selected_cases
        ;;
    *)
        die "Unknown mode: $MODE"
        ;;
esac
