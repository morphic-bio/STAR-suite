#!/usr/bin/env python3
"""Compare every HDF5 dataset and attribute, ignoring only container layout.

Explicit ROOT replacements may normalize provenance paths. Numeric values,
axes, categorical codes, shapes, dtypes, group membership and attributes must
match; no observations or metadata columns are dropped.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def normalized(value, roots):
    value = np.asarray(value)
    if value.dtype.kind not in 'OSU':
        return value

    def scalar(item):
        if isinstance(item, bytes):
            for root in sorted(roots, key=len, reverse=True):
                item = item.replace(root.encode(), b'<ROOT>')
        elif isinstance(item, str):
            for root in sorted(roots, key=len, reverse=True):
                item = item.replace(root, '<ROOT>')
        return item

    return np.array([scalar(item) for item in value.flat], dtype=object).reshape(value.shape)


def compare(a, b, roots_a=(), roots_b=()):
    result = {'a': str(a), 'b': str(b), 'different': [], 'datasets': 0, 'attributes': 0}

    def equal(x, y):
        try:
            np.testing.assert_equal(normalized(x, roots_a), normalized(y, roots_b))
            return True
        except AssertionError:
            return False

    with h5py.File(a) as aa, h5py.File(b) as bb:
        na, nb = {'/': aa}, {'/': bb}
        aa.visititems(lambda name, obj: na.__setitem__(name, obj))
        bb.visititems(lambda name, obj: nb.__setitem__(name, obj))
        result['only_a'] = sorted(na.keys() - nb.keys())
        result['only_b'] = sorted(nb.keys() - na.keys())
        for name in sorted(na.keys() & nb.keys()):
            x, y = na[name], nb[name]
            if type(x) is not type(y):
                result['different'].append(name + ':node_type')
                continue
            if x.attrs.keys() != y.attrs.keys():
                result['different'].append(name + ':attribute_keys')
            for key in x.attrs.keys() & y.attrs.keys():
                result['attributes'] += 1
                if not equal(x.attrs[key], y.attrs[key]):
                    result['different'].append(name + ':attribute:' + key)
            if isinstance(x, h5py.Dataset):
                result['datasets'] += 1
                if x.dtype != y.dtype or x.shape != y.shape or not equal(x[()], y[()]):
                    result['different'].append(name + ':dataset')
    result['pass'] = not (result['different'] or result['only_a'] or result['only_b'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('a', type=Path)
    parser.add_argument('b', type=Path)
    parser.add_argument('--roots-a', nargs='*', default=[])
    parser.add_argument('--roots-b', nargs='*', default=[])
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        parser.error('Refusing to overwrite report')
    result = compare(args.a, args.b, args.roots_a, args.roots_b)
    with args.report.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(result, indent=2))
    return int(not result['pass'])


if __name__ == '__main__':
    raise SystemExit(main())
