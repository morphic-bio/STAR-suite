#!/usr/bin/env python3
"""Record count/merge/dedup decisions from one PF build; reuse a preserved ledger."""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--library-dir', type=Path, required=True)
    parser.add_argument('--include-dir', type=Path, help='Headers matching the tested library')
    parser.add_argument('--expected-ledger', type=Path)
    parser.add_argument('--sanitize', action='store_true')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    pf = Path(__file__).resolve().parents[1] / 'core/features/process_features'
    include = args.include_dir or pf / 'include'
    args.out.mkdir(parents=True, exist_ok=False)
    flags = ['-O2', '-g']
    if args.sanitize:
        flags += ['-fsanitize=address', '-fno-omit-frame-pointer']
    archive = args.library_dir / 'libprocess_features.a'
    driver = pf / 'tests/umi_counts_probe.c'
    record = {'status': 'building', 'started_unix': time.time(), 'commands': [],
              'library_sha256': sha(archive), 'driver_sha256': sha(driver),
              'common_header_sha256': sha(include / 'common.h')}
    manifest = args.out / 'execution.json'
    def save():
        manifest.write_text(json.dumps(record, indent=2) + '\n')
    def run(command, log):
        record['commands'].append(list(map(str, command))); save()
        with log.open('w') as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
        if result.returncode:
            record.update(status='failed', exit_status=result.returncode, failed_log=str(log))
            save(); raise RuntimeError(str(log))
    run(['gcc', *flags, '-I' + str(include), '-I' + str(pf.parent / 'libscrna/include'),
         '-c', driver, '-o', args.out / 'probe.o'], args.out / 'compile.log')
    run(['g++', *flags, args.out / 'probe.o', archive, args.library_dir / 'libscrna.a',
         '-lm', '-lpthread', '-lz', '-fopenmp', '-lhts', '-o', args.out / 'probe'], args.out / 'link.log')
    record['status'] = 'running'; save()
    command = [str(args.out / 'probe')]; record['commands'].append(command); save()
    ledger = args.out / 'ledger.tsv'
    with ledger.open('w') as output, (args.out / 'stderr.log').open('w') as error:
        result = subprocess.run(command, stdout=output, stderr=error)
    record.update(status='passed' if result.returncode == 0 else 'failed',
                  exit_status=result.returncode, finished_unix=time.time(), ledger_sha256=sha(ledger))
    if args.expected_ledger:
        record['expected_ledger_sha256'] = sha(args.expected_ledger)
        record['equal'] = record['ledger_sha256'] == record['expected_ledger_sha256']
        if not record['equal']: record['status'] = 'mismatch'
    save()
    print(json.dumps({key: record[key] for key in ('status', 'exit_status', 'ledger_sha256')}))
    return 0 if record['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
