#!/usr/bin/env python3
"""Check the production spatial half-probe route against small paired inputs.

Requires a real half-probe cache and the matched feature/3350x3350 barcode axes.
Tests the actual STAR binary with an empty genome directory. Every case saves
argv, environment changes, exit status and logs in a fresh output directory.
"""
import argparse
import gzip
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda: f.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('star', 'cache', 'probe-axis', 'contract', 'bc1-oligos', 'bc2-oligos', 'out'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--r1', required=True, help='Comma-separated original R1 lanes')
    p.add_argument('--r2', required=True, help='Comma-separated original R2 lanes')
    p.add_argument('--pairs-per-lane', type=int, default=512)
    p.add_argument('--bgzip', default='bgzip')
    p.add_argument('--resume', action='store_true', help='Reuse verified successful cases; give failed attempts fresh directories')
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=args.resume)
    # Hold one execution owner even when a failed suite is resumed.
    ownership = (args.out / '.execution.lock').open('a')
    try:
        fcntl.flock(ownership.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('Another test harness owns this output directory')
    inputs = args.out / 'inputs'
    inputs.mkdir(exist_ok=args.resume)
    empty = args.out / 'empty_genome'
    empty.mkdir(exist_ok=args.resume)
    paths = {'plain': [[], []], 'gzip': [[], []], 'bgzf': [[], []]}
    count = 0
    source_counts = []
    sources = [args.r2.split(','), args.r1.split(',')]
    assert len(sources[0]) == len(sources[1])
    for mate, lanes in enumerate(sources):
        counts = []
        for lane, source in enumerate(lanes):
            with gzip.open(source, 'rb') as f:
                records = []
                for _ in range(args.pairs_per_lane):
                    lines = [f.readline() for _ in range(4)]
                    if not lines[0]:
                        break
                    assert all(lines), f'Truncated source fixture {source}'
                    records.append(b''.join(lines))
            counts.append(len(records))
            raw = inputs / f'lane{lane}.mate{mate}.fastq'
            raw.write_bytes(b''.join(records))
            gz = raw.with_suffix('.fastq.gz')
            with gzip.GzipFile(filename=str(gz), mode='wb', mtime=0) as f:
                f.write(raw.read_bytes())
            bgzf = raw.with_suffix('.fastq.bgzf')
            with bgzf.open('wb') as f:
                subprocess.run([args.bgzip, '-c', str(raw)], stdout=f, check=True)
            for kind, file in [('plain', raw), ('gzip', gz), ('bgzf', bgzf)]:
                paths[kind][mate].append(str(file))
        source_counts.append(counts)
    assert source_counts[0] == source_counts[1]
    count = sum(source_counts[0])
    assert count > 0
    base = [str(args.star), '--genomeDir', str(empty), '--outSAMtype', 'None',
            '--outSAMattributes', 'None', '--outSJtype', 'None', '--chimSegmentMin', '0',
            '--soloType', 'CB_UMI_Complex', '--soloBarcodeMate', '0',
            '--soloCBwhitelist', str(args.contract / 'bc1_whitelist.txt'), str(args.contract / 'bc2_whitelist.txt'),
            '--soloStrand', 'Unstranded', '--soloFeatures', 'Gene', '--soloUMIdedup', '1MM_CR',
            '--soloCBmatchWLtype', 'EditDist_2', '--soloCBposition', '0_11_0_24', '0_25_0_38',
            '--soloUMIposition', '0_0_0_8', '--soloCellFilter', 'None',
            '--soloUMIfiltering', 'MultiGeneUMI_CR', '--soloMultiMappers', 'Rescue',
            '--soloSkipProcessing', 'yes', '--soloKeysCompat', 'cr',
            '--soloProbeList', str(args.probe_axis), '--soloHashScreenFile', str(args.cache),
            '--soloRemoveDeprecated', 'No', '--soloProbeMismatch', '1', '--flex', 'yes',
            '--soloRunFlexFilter', 'no', '--soloFlexMinimalMemory', 'no',
            '--soloSpatialFlexIntegrated', 'yes', '--soloSpatialBarcodeContract', str(args.contract),
            '--soloSpatialBc1Oligos', str(args.bc1_oligos), '--soloSpatialBc2Oligos', str(args.bc2_oligos),
            '--soloSpatialAssignmentProducts', 'all', '--soloSpatialBinSizes', '2,8,16',
            '--soloSpatialExpectedReads', str(count + 1), '--soloSpatialExpectedCandidates', str(count * 8),
            '--soloSpatialOverflowPolicy', 'Fail', '--readFilesBgzfMode', 'auto']
    results = []
    def run(name, kind='gzip', threads=4, options=(), env_changes=None, mates=None, expected_failure=None):
        effective = list(base)
        appended = []
        # STAR rejects duplicate CLI definitions, so replace scalar defaults.
        for i in range(0, len(options), 2):
            key, value = str(options[i]), str(options[i + 1])
            if key in effective:
                effective[effective.index(key) + 1] = value
            else:
                appended.extend([key, value])
        out = args.out / name
        def command_for(directory):
            return effective + ['--runThreadN', str(threads), '--outFileNamePrefix', str(directory) + '/',
                                '--readFilesIn', *[','.join(m) for m in (mates or paths[kind])], *appended]
        argv = command_for(out)
        previous = None
        if out.exists() and args.resume:
            attempt = out / 'attempt.json'
            previous = json.loads(attempt.read_text()) if attempt.exists() else None
            if not (previous and (previous.get('passed') or previous.get('exit_code') == 0) and
                    previous['argv'] == argv and previous['binary_sha256'] == digest(args.star)):
                previous = None
                retry = 2
                while out.exists():
                    out = args.out / f'{name}_retry{retry}'
                    retry += 1
                argv = command_for(out)
        out.mkdir(exist_ok=previous is not None)
        env = os.environ.copy()
        for key in ('STAR_DISABLE_FLEX_NO_GENOME', 'STAR_FLEX_HASH_H0_ONLY', 'STAR_FLEX_HASH_SCREEN_CACHE'):
            env.pop(key, None)
        env.update(env_changes or {})
        record = previous or {'case': name, 'argv': argv, 'binary_sha256': digest(args.star),
                              'env_overrides': env_changes or {}, 'status': 'running'}
        manifest = out / 'attempt.json'
        if previous is None:
            manifest.write_text(json.dumps(record, indent=2) + '\n')
            start = time.monotonic()
            with (out / 'console.log').open('w') as log:
                proc = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, env=env, timeout=300)
            record.update(exit_code=proc.returncode, wall_seconds=time.monotonic() - start,
                          status='finished')
            manifest.write_text(json.dumps(record, indent=2) + '\n')
        else:
            assert record['env_overrides'] == (env_changes or {})
            print('REUSE completed case', name, flush=True)
        exit_code = record['exit_code']
        logs = '\n'.join(f.read_text(errors='replace') for f in (out / 'console.log', out / 'Log.out') if f.exists())
        assert '..... loading genome' not in logs, name + ' loaded the reference'
        if expected_failure:
            assert exit_code != 0 and expected_failure in logs, (name, exit_code, logs[-2000:])
            assert not (out / 'SpatialFlex.out/run_summary.tsv').exists()
            record['passed'] = True
        else:
            assert exit_code == 0, (name, logs[-3000:])
            assert 'Flex count-only no-genome: active' in logs
            summary = dict(line.split('\t', 1) for line in (out / 'SpatialFlex.out/run_summary.tsv').read_text().splitlines() if '\t' in line)
            assert summary['feature_route'] == 'half_probe_no_alignment'
            assert int(summary['reads_decoded']) == count
            assert sum(int(summary[k]) for k in ('feature_hash_h0', 'feature_hash_h1', 'feature_hash_h1x2', 'feature_hash_deny', 'feature_hash_miss')) == count
            assert int(summary['feature_alignment_resolved']) == int(summary['feature_alignment_unresolved']) == 0
            mex = {str(f.relative_to(out / 'SpatialFlex.out')): digest(f)
                   for f in (out / 'SpatialFlex.out').rglob('*')
                   if f.name in ('matrix.mtx', 'barcodes.tsv', 'features.tsv')}
            assert len(mex) == 36, (name, len(mex))
            assert not list(out.rglob('*.bam'))
            record.update(passed=True, summary=summary, mex=mex)
        manifest.write_text(json.dumps(record, indent=2) + '\n')
        results.append(record)
        (args.out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print('PASS', name, flush=True)
        return record
    first = run('plain_one_worker', 'plain', 1)
    for name, kind, options in [('gzip_four_workers', 'gzip', []),
                               ('gzip_zcat', 'gzip', ['--readFilesCommand', 'zcat']),
                               ('bgzf_four_workers', 'bgzf', []),
                               ('gzip_forced_spill', 'gzip', ['--soloSpatialOverflowPolicy', 'Spill', '--soloSpatialSpillHighWaterCandidates', '16'])]:
        result = run(name, kind, options=options)
        assert result['mex'] == first['mex'], name + ' changed MEX components'
        if kind == 'bgzf':
            assert 'BGZF parallel range readers: active' in (args.out / name / 'Log.out').read_text()
    run('reject_legacy_spatial', options=['--flexLegacy', 'yes'], expected_failure='no longer supported')
    run('reject_no_genome_disabled', env_changes={'STAR_DISABLE_FLEX_NO_GENOME': '1'}, expected_failure='spatial Flex requires the no-genome route')
    for bad in ('short_r1', 'long_r1', 'truncated_r1', 'mismatched_names', 'bad_r2_quality'):
        mates = [list(m) for m in paths['plain']]
        mate = 0 if bad == 'bad_r2_quality' else 1
        lines = Path(mates[mate][0]).read_bytes().splitlines(keepends=True)
        if bad == 'short_r1': del lines[-4:]
        elif bad == 'long_r1': lines.extend(lines[:4])
        elif bad == 'truncated_r1': del lines[-1]
        elif bad == 'mismatched_names': lines[0] = b'@different_name\n'
        else: lines[3] = b'I\n'
        bad_file = inputs / (bad + '.fastq')
        bad_file.write_bytes(b''.join(lines))
        mates[mate][0] = str(bad_file)
        run(bad, mates=mates, expected_failure='spatial FASTQ')
    print(f'{len(results)} cases passed; all five input/thread/spill variants have identical MEX components.')


if __name__ == '__main__':
    main()
