#!/usr/bin/env python3
"""Strict biological output equality for two completed scheduler benchmark arms."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def content(path):
    return gzip.open(path, 'rb') if path.suffix == '.gz' else path.open('rb')


def digest(path):
    h = hashlib.sha256()
    with content(path) as source:
        for part in iter(lambda: source.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def locate(directory, name):
    path = directory / name
    if path.is_file():
        return path
    path = directory / (name + '.gz')
    if path.is_file():
        return path
    raise FileNotFoundError(directory / name)


def labels(directory, name):
    with content(locate(directory, name)) as source:
        return source.read().decode().splitlines()


def compare_mex(left, right):
    lf, rf = labels(left, 'features.tsv'), labels(right, 'features.tsv')
    lb, rb = labels(left, 'barcodes.tsv'), labels(right, 'barcodes.tsv')
    lm, rm = locate(left, 'matrix.mtx'), locate(right, 'matrix.mtx')
    if lf == rf and lb == rb and digest(lm) == digest(rm):
        return {'equal': True, 'method': 'identical_uncompressed_files', 'barcodes': len(lb), 'features': len(lf)}
    if sorted(lf) != sorted(rf) or sorted(lb) != sorted(rb):
        return {'equal': False, 'reason': 'label_sets_differ',
                'left_barcodes': len(lb), 'right_barcodes': len(rb),
                'left_only_barcodes': sorted(set(lb)-set(rb))[:20],
                'right_only_barcodes': sorted(set(rb)-set(lb))[:20]}
    from scipy.io import mmread
    import numpy as np
    with content(lm) as source:
        a = mmread(source).tocsr()
    with content(rm) as source:
        b = mmread(source).tocsr()
    if lf != rf or lb != rb:
        assert len(set(lf)) == len(lf) and len(set(lb)) == len(lb), 'Ambiguous duplicate labels'
        fi, bi = {v: i for i, v in enumerate(rf)}, {v: i for i, v in enumerate(rb)}
        b = b[np.array([fi[v] for v in lf]), :][:, np.array([bi[v] for v in lb])]
    difference = (a - b).tocoo()
    difference.eliminate_zeros()
    return {'equal': difference.nnz == 0, 'method': 'exact_sparse_difference',
            'differing_entries': int(difference.nnz), 'total_count_delta': float(difference.sum()),
            'examples': [{'feature': lf[r], 'barcode': lb[c], 'left_minus_right': float(v)}
                         for r, c, v in zip(difference.row[:20], difference.col[:20], difference.data[:20])]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--left', type=Path, required=True)
    parser.add_argument('--right', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--allow-different-binaries', action='store_true',
                        help='Compare a code change against its preserved control; still require exact output equality')
    args = parser.parse_args()
    records = [json.loads((p / 'execution.json').read_text()) for p in (args.left, args.right)]
    assert all(r['status'] == 'passed' and r['exit_status'] == 0 for r in records), 'Incomplete arm'
    same_binary = records[0]['binary_sha256'] == records[1]['binary_sha256']
    assert same_binary or args.allow_different_binaries, 'Matched control requires identical binary (use --allow-different-binaries for code-change validation)'
    roots = [p / 'out' for p in (args.left, args.right)]
    dirs = [{m.parent.relative_to(p) for m in p.rglob('matrix.mtx*')} for p in roots]
    result = {'left': str(args.left), 'right': str(args.right),
              'binary_sha256': records[0]['binary_sha256'] if same_binary else None,
              'left_binary_sha256': records[0]['binary_sha256'], 'right_binary_sha256': records[1]['binary_sha256'],
              'comparison_mode': 'same_binary' if same_binary else 'cross_build',
              'matrix_directory_sets_equal': dirs[0] == dirs[1], 'matrices': {}, 'guide_calls': {}}
    for directory in sorted(dirs[0] & dirs[1]):
        result['matrices'][str(directory)] = compare_mex(roots[0] / directory, roots[1] / directory)
    guides = [{p.relative_to(root) for p in (root / 'outs/crispr_analysis').glob('*.csv')}
              for root in roots]
    result['guide_directory_sets_equal'] = guides[0] == guides[1]
    for relative in sorted(guides[0] & guides[1]):
        path = roots[0] / relative
        other = roots[1] / relative
        with path.open() as source:
            left = list(csv.reader(source))
        with other.open() as source:
            right = list(csv.reader(source))
        equal = left[:1] == right[:1] and sorted(left[1:]) == sorted(right[1:])
        result['guide_calls'][str(relative)] = {'equal': equal, 'left_rows': len(left)-1, 'right_rows': len(right)-1}
    assert result['matrices'] and result['guide_calls'], 'Missing comparison surfaces'
    result['equal'] = result['matrix_directory_sets_equal'] and result['guide_directory_sets_equal'] and all(
        x['equal'] for group in (result['matrices'], result['guide_calls']) for x in group.values())
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'equal': result['equal'], 'matrices': len(result['matrices']),
                      'guide_tables': len(result['guide_calls'])}))
    raise SystemExit(0 if result['equal'] else 1)


if __name__ == '__main__':
    main()
