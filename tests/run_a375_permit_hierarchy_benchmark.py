#!/usr/bin/env python3
"""One full A375 arm, with durable execution provenance and process-tree CPU samples."""
import argparse
import fcntl
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for data in iter(lambda: source.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def write_json(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def process_tree(pid):
    result = []
    pending = [pid]
    while pending:
        current = pending.pop()
        base = Path('/proc') / str(current)
        try:
            stat = (base / 'stat').read_text()
            fields = stat[stat.rfind(')') + 2:].split()
            result.append({'pid': current, 'name': stat[stat.find('(')+1:stat.rfind(')')],
                           'start_ticks': int(fields[19]), 'user_ticks': int(fields[11]),
                           'system_ticks': int(fields[12]), 'reaped_child_user_ticks': int(fields[13]),
                           'reaped_child_system_ticks': int(fields[14]), 'rss_pages': int(fields[21])})
            pending.extend(int(x) for x in (base / 'task' / str(current) / 'children').read_text().split())
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--hierarchy', type=int, choices=(0, 1), required=True)
    parser.add_argument('--balance', type=int, choices=(0, 1))
    parser.add_argument('--bgzf-reader-threads', type=int)
    parser.add_argument('--build-provenance', type=Path,
                        help='Earlier arm on this exact binary, preserving its compiled source provenance')
    parser.add_argument('--timeout', type=float, default=7200)
    args = parser.parse_args()
    if args.bgzf_reader_threads is not None and args.bgzf_reader_threads < 1:
        parser.error('--bgzf-reader-threads must be positive')
    repo = Path(__file__).resolve().parents[1]
    lock = Path('/tmp/star-suite-benchmark.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for process in Path('/proc').iterdir():
        if not process.name.isdigit():
            continue
        try:
            exe = (process / 'exe').resolve(strict=True)
            if exe.name == 'STAR':
                raise RuntimeError(f'STAR already running: PID {process.name}, {exe}')
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
    args.out.mkdir(parents=True, exist_ok=False)
    out = args.out / 'out'
    out.mkdir()
    references = json.loads((repo / 'docs/ADAPTIVE_PERMIT_A375_REFERENCE_RUNS_20260911.json').read_text())
    reference = next(x for x in references['runs'] if x['run_id'] == 'a375_inproc_featfloor8_fifo')
    argv = list(reference['argv'])
    argv[0] = str(args.binary.resolve())

    def set_arg(name, value):
        if name in argv:
            argv[argv.index(name) + 1] = str(value)
        else:
            argv.extend([name, str(value)])

    original_config = Path(argv[argv.index('--pfMultiConfig') + 1])
    config = out / 'multi_config.csv'
    shutil.copy2(original_config, config)
    set_arg('--outFileNamePrefix', str(out) + '/')
    set_arg('--pfMultiConfig', config)
    set_arg('--dynamicThreadMapFloor', 24)
    set_arg('--dynamicThreadFeatureFloor', 8)
    set_arg('--dynamicThreadPfControllerMode', 'off')
    set_arg('--dynamicThreadAtacController', 0)
    set_arg('--dynamicThreadBgzfHierarchy', args.hierarchy)
    if args.balance is not None:
        set_arg('--dynamicThreadBalance', args.balance)
    if args.bgzf_reader_threads is not None:
        set_arg('--bgzfReaderThreads', args.bgzf_reader_threads)
    env = dict(os.environ)
    for name in ('STAR_PF_FEATURE_BOOTSTRAP_READS', 'STAR_PF_ASYNC_ASSIGN'):
        env.pop(name, None)
    env['STAR_SOLO_NONFLEX_HASH_BRIDGE'] = '1'
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    patch = subprocess.check_output(['git', 'diff', '--binary'], cwd=repo)
    (args.out / 'source.patch').write_bytes(patch)
    changed = subprocess.check_output(['git', 'diff', '--name-only'], cwd=repo, text=True).splitlines()
    untracked = subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard'], cwd=repo, text=True).splitlines()
    source_files = {}
    for name in sorted(set(changed + untracked)):
        path = repo / name
        if path.is_file() and path.stat().st_size < 2 * 1024 * 1024:
            source_files[name] = sha(path)
            if name in untracked:
                dest = args.out / 'untracked_source' / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
    inputs = {}
    for flag in ('--readFilesIn', '--genomeDir', '--soloCBwhitelist', '--crFeatureRef'):
        index = argv.index(flag) + 1
        values = []
        while index < len(argv) and not argv[index].startswith('--'):
            values.extend(argv[index].split(','))
            index += 1
        for value in values:
            p = Path(value)
            stat = p.stat()
            inputs[value] = {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                             'inode': stat.st_ino, 'device': stat.st_dev}
            if p.is_file() and stat.st_size < 16 * 1024 * 1024:
                inputs[value]['sha256'] = sha(p)
    build_origin = None
    binary_hash = sha(args.binary)
    if args.build_provenance:
        build_origin = json.loads((args.build_provenance / 'execution.json').read_text())
        assert build_origin['binary_sha256'] == binary_hash, 'Build provenance is for a different binary'
        revision = build_origin['source_revision']
        source_files = build_origin['source_files']
        shutil.copy2(args.build_provenance / 'source.patch', args.out / 'source.patch')
        shutil.rmtree(args.out / 'untracked_source', ignore_errors=True)
        if (args.build_provenance / 'untracked_source').is_dir():
            shutil.copytree(args.build_provenance / 'untracked_source', args.out / 'untracked_source')
    record = {'status': 'prepared', 'argv': argv, 'binary_sha256': binary_hash,
              'source_revision': revision, 'source_patch_sha256': sha(args.out / 'source.patch'),
              'source_files': source_files, 'input_identities': inputs,
              'multi_config_sha256': sha(config), 'hostname': os.uname().nodename,
              'cpu_count': os.cpu_count(), 'environment': {k: v for k, v in env.items()
                  if k.startswith(('STAR_', 'OMP_', 'GOMP_', 'PF_TEST_'))},
              'reference_run': reference['path'], 'parent_floors': {'map': 24, 'feature': 8},
              'purpose': 'full A375 matched hierarchy comparison', 'prepared_unix': time.time(),
              'cpu_accounting': 'GNU time records whole-job CPU; /proc samples retain own and reaped-child CPU. Per-process sampled CPU is a lower bound.'}
    if build_origin:
        record['build_provenance_from'] = str(args.build_provenance.resolve())
    record['benchmark_driver_sha256'] = sha(Path(__file__))
    record['bgzf_reader_threads_override'] = args.bgzf_reader_threads
    provenance = args.out / 'execution.json'
    write_json(provenance, record)
    start = time.monotonic()
    totals = {}
    with (args.out / 'stdout.log').open('w') as stdout, (args.out / 'stderr.log').open('w') as stderr, \
            (args.out / 'cpu_samples.jsonl').open('w') as cpu:
        child = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(args.out / 'time.txt')] + argv,
                                 stdout=stdout, stderr=stderr, env=env, start_new_session=True)
        record.update(status='running', pid=child.pid, started_unix=time.time())
        write_json(provenance, record)
        try:
            while child.poll() is None:
                elapsed = time.monotonic() - start
                rows = process_tree(child.pid)
                cpu.write(json.dumps({'seconds': elapsed, 'processes': rows}) + '\n')
                cpu.flush()
                for row in rows:
                    totals[f"{row['pid']}:{row['start_ticks']}"] = row
                if elapsed > args.timeout:
                    raise TimeoutError('Benchmark exceeded its wall-time limit')
                time.sleep(0.25)
        except BaseException as error:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            record['interruption'] = repr(error)
            raise
        finally:
            record.update(status='complete' if child.returncode == 0 else 'failed',
                          exit_status=child.returncode, elapsed_seconds=time.monotonic() - start,
                          finished_unix=time.time(), sampled_process_cpu=totals,
                          clock_ticks_per_second=os.sysconf('SC_CLK_TCK'))
            write_json(provenance, record)
    required = ['Log.final.out', 'Solo.out/GeneFull/filtered/matrix.mtx',
                'outs/filtered_feature_bc_matrix/matrix.mtx.gz',
                'outs/crispr_analysis/protospacer_calls_per_cell.csv']
    missing = [name for name in required if not (out / name).is_file()]
    record.update(missing_outputs=missing,
                  finished_marker='finished successfully' in (args.out / 'stdout.log').read_text())
    success = record['exit_status'] == 0 and not missing and record['finished_marker']
    record['status'] = 'passed' if success else 'failed'
    write_json(provenance, record)
    print(json.dumps({k: record[k] for k in ('status', 'exit_status', 'elapsed_seconds', 'missing_outputs')}))
    raise SystemExit(0 if success else 1)


if __name__ == '__main__':
    main()
