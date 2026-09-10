#!/usr/bin/env python3
"""Exact OCM routing/caller regression; accepts saved controls to avoid reruns.

Build STAR's ocm-multi-unit-tests target first. Each execution has a completion
record; outputs must be fresh. No reference index or FASTQ alignment is used.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def payload(path):
    return gzip.open(path, 'rb').read() if path.suffix == '.gz' else path.read_bytes()


def compare(control, candidate):
    def files(root):
        return {p.relative_to(root): p for name in ('outs', 'samples')
                for p in (root / name).rglob('*') if p.is_file()}
    a, b = files(control), files(candidate)
    assert a.keys() == b.keys(), (sorted(a.keys() - b.keys()), sorted(b.keys() - a.keys()))
    for name, path in a.items():
        x, y = payload(path), payload(b[name])
        if name.name == 'ocm_materialization_summary.json':
            x, y = json.loads(x), json.loads(y)
            # The copied fixture's config location is the only allowed change.
            x.pop('config_path'); y.pop('config_path')
        assert x == y, f'Output mismatch: {name}'
    return len(a)


def prepare(source, dest, native):
    dest.mkdir(parents=True, exist_ok=False)
    shutil.copy2(source / 'config.csv', dest / 'config.csv')
    shutil.copytree(source / 'run/Solo.out', dest / 'run/Solo.out')
    if native:
        filtered = dest / 'run/Solo.out/GeneFull/filtered'
        if filtered.exists():
            shutil.rmtree(filtered)
    return dest


def run(binary, root, threads=4, memory=1 << 30, success=True):
    env = dict(os.environ)
    for key in ('OCM_TEST_RUN_DIR', 'OCM_TEST_CONFIG'):
        env.pop(key, None)
    overrides = dict(OCM_TEST_FIXTURE_ROOT=str(root), OCM_TEST_THREADS=str(threads),
                     OCM_TEST_MEMORY=str(memory), OCM_TEST_LOG=str(root / 'caller.log'),
                     STAR_SOLO_MEMORY_PROFILE='1', OMP_NUM_THREADS=str(threads))
    env.update(overrides)
    args = [str(binary), 'materialize']
    start = time.monotonic()
    with (root / 'stdout.log').open('w') as output:
        result = subprocess.run(args, env=env, stdout=output, stderr=subprocess.STDOUT)
    record = dict(argv=args, env=overrides, exit=result.returncode,
                  wall_seconds=time.monotonic() - start,
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest())
    (root / 'COMPLETE.json').write_text(json.dumps(record, indent=2) + '\n')
    assert (result.returncode == 0) == success, (root, record)
    if not success:
        assert not list((root / 'outs/per_sample_outs').rglob('matrix.mtx.gz')), root
    assert not list(root.rglob('ocm-counts-*')), root
    assert not list(root.rglob('*.tmp')), root
    assert not list(root.rglob('*.body')), root
    return root


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--saved-native', type=Path)
    parser.add_argument('--saved-existing', type=Path)
    parser.add_argument('--saved-ordered', type=Path,
                        help='Previous regression directory containing ordered_*_control results')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixture = Path(__file__).resolve().parent / 'fixtures/ocm_multi_tiny'
    checks = {}
    for native, saved in ((True, args.saved_native), (False, args.saved_existing)):
        label = 'native' if native else 'existing'
        control = saved or run(args.baseline, prepare(fixture, root / (label + '_control'), native))
        changed = run(args.candidate, prepare(fixture, root / label, native))
        checks[label] = compare(control, changed)
        if native:
            assert 'caller_mex_reads=0' in (changed / 'caller.log').read_text()
            assert 'spilled_samples=0' in (changed / 'caller.log').read_text()
            spilled = run(args.candidate, prepare(fixture, root / 'native_spilled', True), memory=0)
            checks['native_spilled'] = compare(control, spilled)
            assert 'spilled_samples=5' in (spilled / 'caller.log').read_text()
        else:
            assert 'raw_matrix_passes=1 filtered_matrix_passes=1' in (changed / 'caller.log').read_text()
    # Interleaved columns, duplicate entries, zeros, fractional values, and a
    # zero-UMI known-tag column expose order/remap and native/export policy drift.
    special = prepare(fixture, root / 'special_input', True)
    matrix = special / 'run/Solo.out/GeneFull/raw/matrix.mtx'
    matrix.write_text('%%MatrixMarket matrix coordinate real general\n%\n2 5 9\n'
                      '2 4 0\n1 1 1.6\n1 2 0\n1 3 3\n2 2 2\n'
                      '2 1 2\n1 1 1\n2 3 0\n2 5 0\n')
    for native in (True, False):
        if not native:
            shutil.copytree(fixture / 'run/Solo.out/GeneFull/filtered',
                            special / 'run/Solo.out/GeneFull/filtered')
        label = 'ordered_native' if native else 'ordered_existing'
        control = (args.saved_ordered / (label + '_control') if args.saved_ordered else
                   run(args.baseline, prepare(special, root / (label + '_control'), native)))
        changed = run(args.candidate, prepare(special, root / label, native), memory=0)
        checks[label] = compare(control, changed)
        compressed = prepare(special, root / (label + '_gzip'), native)
        path = compressed / 'run/Solo.out/GeneFull/raw/matrix.mtx'
        with gzip.open(str(path) + '.gz', 'wb') as out:
            out.write(path.read_bytes())
        path.unlink()
        run(args.candidate, compressed, threads=1)
        checks[label + '_gzip'] = compare(control, compressed)
    # Count/shape/truncation failures must stop before publishing sample MEXs.
    bad_entries = {
        'bad_column': '2 5 1\n1 6 1\n',
        'bad_row': '2 5 1\n3 1 1\n',
        'bad_shape': '3 5 1\n1 1 1\n',
        'truncated': '2 5 2\n1 1 1\n',
        'umi_overflow': '2 5 2\n1 1 4294967295\n2 1 1\n',
        'count_overflow': '2 5 1\n1 1 4294967296\n',
        'bad_syntax': '2 5 1\nnot a count\n',
    }
    for label, body in bad_entries.items():
        changed = prepare(fixture, root / label, True)
        (changed / 'run/Solo.out/GeneFull/raw/matrix.mtx').write_text(
            '%%MatrixMarket matrix coordinate integer general\n%\n' + body)
        run(args.candidate, changed, success=False, memory=0)
        checks[label] = 'rejected'
    changed = prepare(fixture, root / 'bad_shape_existing', False)
    (changed / 'run/Solo.out/GeneFull/raw/matrix.mtx').write_text(
        '%%MatrixMarket matrix coordinate integer general\n%\n3 5 1\n1 1 1\n')
    run(args.candidate, changed, success=False)
    checks['bad_shape_existing'] = 'rejected'
    for label in ('bad_crc', 'truncated_gzip'):
        changed = prepare(fixture, root / label, True)
        path = changed / 'run/Solo.out/GeneFull/raw/matrix.mtx'
        data = bytearray(gzip.compress(path.read_bytes()))
        if label == 'bad_crc':
            data[-8] ^= 1
        else:
            data = data[:-6]
        Path(str(path) + '.gz').write_bytes(data)
        path.unlink()
        run(args.candidate, changed, success=False)
        checks[label] = 'rejected'
    (root / 'PARITY.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
