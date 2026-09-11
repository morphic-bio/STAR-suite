#!/usr/bin/env python3
"""Compare decisions and learning state with a preserved pre-change PF library."""
import argparse
import csv
import gzip
import hashlib
import json
import random
import subprocess
import time
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-dir', type=Path, required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--larry-subset', type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    pf = repo / 'core/features/process_features'
    args.out.mkdir(parents=True, exist_ok=False)
    record = {'status': 'building', 'started_unix': time.time(), 'commands': [], 'cases': {},
              'driver_sha256': sha(pf / 'tests/bootstrap_anchor_probe.c'),
              'libraries': {kind: sha(directory / 'libprocess_features.a')
                            for kind, directory in (('before', args.baseline_dir), ('after', args.build_dir))}}
    def save():
        (args.out / 'execution.json').write_text(json.dumps(record, indent=2) + '\n')
    def execute(command, log):
        record['commands'].append([str(p) for p in command]); save()
        with log.open('w') as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
        if result.returncode:
            record.update(status='failed', exit_status=result.returncode, failed_log=str(log)); save()
            raise RuntimeError(str(log))
    for kind, directory in (('before', args.baseline_dir), ('after', args.build_dir)):
        execute(['gcc', '-O2', '-g', '-I' + str(pf / 'include'), '-I' + str(pf.parent / 'libscrna/include'),
                 '-c', pf / 'tests/bootstrap_anchor_probe.c', '-o', args.out / (kind + '.o')], args.out / (kind + '_compile.log'))
        execute(['g++', args.out / (kind + '.o'), directory / 'libprocess_features.a', directory / 'libscrna.a',
                 '-lm', '-lpthread', '-lz', '-fopenmp', '-lhts', '-o', args.out / kind], args.out / (kind + '_link.log'))
    record['status'] = 'running'; save()
    random.seed(34017)
    cases = []
    def dna(n): return ''.join(random.choice('ACGT') for _ in range(n))
    for shape, length, count in [('suffix', 12, 160), ('prefix', 12, 160), ('suffix128', 40, 160),
                                  ('mixed_length', 12, 160), ('dual', 12, 160), ('singleton', 12, 1),
                                  ('short_h2', 8, 8), ('large_h2', 8, 160)]:
        reference = args.out / (shape + '.csv')
        features = []
        while len(features) < count:
            s = dna(length)
            if s not in features: features.append(s)
        features[0] = features[0][:-3] + 'AAA'
        if count > 1:
            features[1] = features[0][:-1] + next(b for b in 'ACGT' if b != features[0][-1])
        if shape == 'mixed_length': features[-1] += 'AC'
        pattern = 'GCTA(BC)' if shape == 'prefix' else 'GCTA(BC)TGAC' if shape == 'dual' else '(BC)TGAC'
        with reference.open('w') as f:
            writer = csv.writer(f); writer.writerow(['name', 'sequence', 'pattern'])
            writer.writerows((f'f{i}', sequence, pattern) for i, sequence in enumerate(features))
        def anchored(s):
            return ('GCTA' if shape in ('prefix', 'dual') else '') + s + ('' if shape == 'prefix' else 'TGAC')
        queries = []
        for i in range(900):
            s = random.choice(features)
            if i % 4 == 1:
                pos = random.randrange(len(s)); s = s[:pos] + random.choice('ACGTN') + s[pos+1:]
            q = dna(i % 5) + anchored(s) + dna(i % 7)
            if i % 5 == 0: q += anchored(random.choice(features))
            if i % 9 == 0: q += anchored(s)
            if i % 11 == 0: q = dna(20) + s + dna(20)  # No required anchor: broad learning fallback.
            queries.append(q)
        queries += ['', 'TGAC', 'GCTA', 'N' * length, anchored(features[0]) + anchored(features[-1]),
                    anchored(features[-1]) + anchored(features[0]),
                    features[0][:-1], features[0][:-2], features[0][:-3]]
        query_path = args.out / (shape + '.sequences')
        query_path.write_text('\n'.join(queries) + '\n')
        for learning, hamming, max_n in [(100000, 1, 0), (50, 1, 1), (100000, 0, 0)]:
            cases.append((f'{shape}_{learning}_{hamming}_{max_n}', reference, query_path, learning, hamming, max_n))
        if shape in ('short_h2', 'large_h2'):
            cases.append((shape + '_distance2', reference, query_path, 100000, 2, 1))
    if args.larry_subset:
        reference = args.larry_subset / 'references/ref_feature_larryBC.csv'
        paths = sorted((args.larry_subset / 'fastqs/LARRY').glob('*_R2_*.fastq.gz'))
        assert len(paths) == 8
        queries = []
        streams = [gzip.open(p, 'rt') for p in paths]
        try:
            for _ in range(200):
                for stream in streams:
                    lines = [stream.readline() for _ in range(4)]
                    assert all(lines); queries.append(lines[1].rstrip('\r\n'))
        finally:
            for stream in streams: stream.close()
        path = args.out / 'larry_real_1600.sequences'
        path.write_text('\n'.join(queries) + '\n')
        cases.append(('larry_real_1600', reference, path, 100000, 1, 0))
    for name, reference, queries, learning, hamming, max_n in cases:
        row = {'reference_sha256': sha(reference), 'queries_sha256': sha(queries),
               'learning': learning, 'hamming': hamming, 'max_n': max_n}
        for kind in ('before', 'after'):
            ledger = args.out / (name + '_' + kind + '.tsv')
            start = time.time()
            execute([args.out / kind, reference, queries, ledger, str(learning), str(hamming), str(max_n)],
                    args.out / (name + '_' + kind + '.log'))
            row[kind] = {'ledger_sha256': sha(ledger), 'seconds': time.time() - start}
        row['equal'] = row['before']['ledger_sha256'] == row['after']['ledger_sha256']
        record['cases'][name] = row; save()
        if not row['equal']:
            record['status'] = 'mismatch'; save(); raise RuntimeError(name)
        print(name, 'exact decisions/histograms/modes PASS', flush=True)
    record.update(status='passed', finished_unix=time.time()); save()


if __name__ == '__main__':
    main()
