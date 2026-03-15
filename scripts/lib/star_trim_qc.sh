#!/usr/bin/env bash

star_trim_qc_is_enabled() {
  local value="${STAR_TRIM_QC_ENABLE:-0}"
  case "${value,,}" in
    1|true|yes|on)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

star_trim_qc_basename() {
  printf '%s\n' "${STAR_TRIM_QC_BASENAME:-read_qc}"
}

star_trim_qc_max_reads() {
  printf '%s\n' "${STAR_TRIM_QC_MAX_READS:-250000}"
}

star_trim_qc_report_prefix() {
  local out_dir="$1"
  local basename="${2:-$(star_trim_qc_basename)}"
  printf '%s/%s\n' "${out_dir%/}" "${basename}"
}

star_trim_qc_append_args() {
  local array_name="$1"
  local out_dir="$2"
  local basename="${3:-$(star_trim_qc_basename)}"
  local max_reads="${4:-$(star_trim_qc_max_reads)}"
  local prefix
  local -n cmd_ref="${array_name}"

  if ! star_trim_qc_is_enabled; then
    return 0
  fi

  prefix="$(star_trim_qc_report_prefix "${out_dir}" "${basename}")"
  cmd_ref+=(--trimQcReport "${prefix}")
  if [[ -n "${max_reads}" ]]; then
    cmd_ref+=(--trimQcMaxReads "${max_reads}")
  fi
}

star_trim_qc_write_manifest() {
  local manifest="$1"
  local out_dir="$2"
  local basename="${3:-$(star_trim_qc_basename)}"
  local max_reads="${4:-$(star_trim_qc_max_reads)}"

  if star_trim_qc_is_enabled; then
    printf 'trim_qc_enabled=1\n' >> "${manifest}"
    printf 'trim_qc_prefix=%s\n' "$(star_trim_qc_report_prefix "${out_dir}" "${basename}")" >> "${manifest}"
    printf 'trim_qc_max_reads=%s\n' "${max_reads}" >> "${manifest}"
  else
    printf 'trim_qc_enabled=0\n' >> "${manifest}"
  fi
}

star_trim_qc_list_outputs() {
  local out_dir="$1"
  local basename="${2:-$(star_trim_qc_basename)}"
  local prefix
  prefix="$(star_trim_qc_report_prefix "${out_dir}" "${basename}")"
  printf '%s\n' "${prefix}.trim_qc.json"
  printf '%s\n' "${prefix}.trim_qc.html"
}
