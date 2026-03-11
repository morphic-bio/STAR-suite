#!/usr/bin/env python3
import argparse
import gzip
import pathlib
import re
import sys

SUFFIX_RE = re.compile(r"-(\d+)$")
COUNT_RE = re.compile(r"CR-compat MEX filtering: Using GEX barcodes only \((\d+) barcodes\)")


def open_text(path: pathlib.Path):
    with path.open('rb') as fh:
        magic = fh.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return path.open('r')


def normalize_barcode(raw: str) -> str:
    raw = raw.strip()
    if not raw:
        return raw
    return SUFFIX_RE.sub('', raw)


def complement_base(base: str) -> str:
    return {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }.get(base, base)


def translate_nxt_middle_two_bases(barcode: str) -> str:
    chars = list(barcode)
    if len(chars) >= 9:
        chars[7] = complement_base(chars[7].upper())
        chars[8] = complement_base(chars[8].upper())
    return ''.join(chars)


def read_barcode_set(path: pathlib.Path):
    out = []
    seen = set()
    with open_text(path) as handle:
        for line in handle:
            bc = normalize_barcode(line)
            if not bc:
                continue
            if bc not in seen:
                seen.add(bc)
                out.append(bc)
    return out, seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--solo-barcodes', required=True)
    ap.add_argument('--merged-barcodes', required=True)
    ap.add_argument('--subset-barcodes', action='append', default=[])
    ap.add_argument('--log', required=False)
    ap.add_argument('--report', required=False)
    args = ap.parse_args()

    solo_path = pathlib.Path(args.solo_barcodes)
    merged_path = pathlib.Path(args.merged_barcodes)
    solo_list, solo_set = read_barcode_set(solo_path)
    merged_list, merged_set = read_barcode_set(merged_path)

    lines = []
    lines.append(f"solo_count\t{len(solo_set)}")
    lines.append(f"merged_count\t{len(merged_set)}")

    if solo_set != merged_set:
        solo_only = sorted(solo_set - merged_set)
        merged_only = sorted(merged_set - solo_set)
        lines.append(f"solo_only_count\t{len(solo_only)}")
        lines.append(f"merged_only_count\t{len(merged_only)}")
        if solo_only:
            lines.append(f"solo_only_example\t{solo_only[0]}")
        if merged_only:
            lines.append(f"merged_only_example\t{merged_only[0]}")
        emit(lines, args.report)
        print('Filtered GEX barcode mismatch between Solo and merged filtered MEX', file=sys.stderr)
        return 1

    lines.append('solo_vs_merged\tMATCH')

    for subset_arg in args.subset_barcodes:
        subset_path_str, subset_source_chem = parse_subset_arg(subset_arg)
        subset_path = pathlib.Path(subset_path_str)
        _, subset_set = read_barcode_set(subset_path)
        if subset_source_chem == 'NXT':
            subset_set = {translate_nxt_middle_two_bases(x) for x in subset_set}
        extra = sorted(subset_set - merged_set)
        lines.append(f"subset_count[{subset_path}]\t{len(subset_set)}")
        lines.append(f"subset_source_chem[{subset_path}]\t{subset_source_chem}")
        if extra:
            lines.append(f"subset_status[{subset_path}]\tFAIL")
            lines.append(f"subset_extra_count[{subset_path}]\t{len(extra)}")
            lines.append(f"subset_extra_example[{subset_path}]\t{extra[0]}")
            emit(lines, args.report)
            print(f'Subset barcode mismatch for {subset_path}', file=sys.stderr)
            return 1
        lines.append(f"subset_status[{subset_path}]\tOK")

    if args.log:
        log_text = pathlib.Path(args.log).read_text()
        counts = {int(m.group(1)) for m in COUNT_RE.finditer(log_text)}
        lines.append(f"log_filter_counts\t{','.join(str(x) for x in sorted(counts)) if counts else 'NONE'}")
        if not counts:
            emit(lines, args.report)
            print('No CR-compat GEX filtering count found in log', file=sys.stderr)
            return 1
        if counts != {len(merged_set)}:
            emit(lines, args.report)
            print('Logged CR-compat GEX barcode count does not match merged filtered barcodes', file=sys.stderr)
            return 1
        lines.append('log_count_match\tOK')

    emit(lines, args.report)
    return 0


def emit(lines, report_path):
    text = '\n'.join(lines) + '\n'
    if report_path:
        pathlib.Path(report_path).write_text(text)
    else:
        sys.stdout.write(text)


def parse_subset_arg(raw: str):
    if '::' in raw:
        path, chem = raw.rsplit('::', 1)
        chem = chem.strip().upper()
        if chem not in ('TRU', 'NXT'):
            raise ValueError(f'Unsupported subset chemistry {chem!r} in {raw!r}')
        return path, chem
    return raw, 'TRU'


if __name__ == '__main__':
    raise SystemExit(main())
