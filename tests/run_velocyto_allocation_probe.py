"""Run one 100K-GEX Velocyto validation using a preserved A375 run's inputs."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('--reference-run', required=True, type=Path)
parser.add_argument('--binary', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--mode', choices=['stream', 'sorted', 'buckets'], default='stream')
parser.add_argument('--gex-counting', choices=['legacy', 'bridge'], default='legacy',
                    help='Legacy supplies Velocyto read information without requiring BAM tags')
args = parser.parse_args()
reference = json.loads((args.reference_run / 'execution.json').read_text())
assert reference['status'] == 'passed'
lock = Path('/tmp/star-suite-benchmark.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
out = args.out.resolve()
out.mkdir(exist_ok=False)
(out / 'out').mkdir()
binary = args.binary.resolve()
command = list(reference['argv'])
command[0] = str(binary)
command[command.index('--outFileNamePrefix') + 1] = str(out / 'out') + '/'
config_index = command.index('--pfMultiConfig') + 1
shutil.copy2(command[config_index], out / 'out/multi_config.csv')
command[config_index] = str(out / 'out/multi_config.csv')
command.insert(command.index('--soloFeatures') + 2, 'Velocyto')
command.extend(['--readMapNumber', '100000'])
env = {k: v for k, v in os.environ.items() if not k.startswith(('STAR_', 'PF_', 'OMP_'))}
env.update(reference['environment'])
if args.gex_counting == 'legacy':
    command[command.index('--soloInlineHashMode') + 1] = 'no'
    env.pop('STAR_SOLO_NONFLEX_HASH_BRIDGE', None) # The switch tests presence, not its value.
if args.mode != 'stream':
    env['STAR_VELOCYTO_DETERMINISTIC_REPLAY'] = '1'
if args.mode == 'buckets':
    env['STAR_VELOCYTO_INTEGRATED_HASH'] = '1'
build = json.loads((binary.parent / 'execution.json').read_text())
manifest = {
    'status': 'running', 'argv': command,
    'purpose': '100K GEX reads plus full A375 guide library; Velocyto allocation validation',
    'mode': args.mode, 'reference_run': str(args.reference_run.resolve()),
    'gex_counting': args.gex_counting,
    'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'input_identities': reference['input_identities'],
    'environment': {k: v for k, v in env.items() if k.startswith(('STAR_', 'PF_', 'OMP_'))},
    'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
    'build_manifest': build,
    'source_patch_sha256': hashlib.sha256((binary.parent / 'source.patch').read_bytes()).hexdigest(),
    'started_unix': time.time(),
}
status = out / 'execution.json'
status.write_text(json.dumps(manifest, indent=2) + '\n')
with (out / 'stdout.log').open('w') as stdout, (out / 'stderr.log').open('w') as stderr:
    result = subprocess.run(['/usr/bin/time', '-v', '-o', str(out / 'time.txt'), *command],
                            env=env, stdout=stdout, stderr=stderr)
manifest.update(exit_status=result.returncode, finished_unix=time.time())
manifest['elapsed_seconds'] = manifest['finished_unix'] - manifest['started_unix']
manifest['status'] = 'passed' if result.returncode == 0 and (out / 'out/Log.final.out').exists() else 'failed'
if manifest['status'] == 'passed':
    entries = 0
    for layer in ('spliced', 'unspliced', 'ambiguous'):
        with (out / 'out/Solo.out/Velocyto/raw' / (layer + '.mtx')).open() as matrix:
            dims = next(line.split() for line in matrix if line.strip() and not line.startswith('%'))
        entries += int(dims[2])
    manifest['velocyto_nonzero_entries'] = entries
    if entries == 0:
        manifest['status'] = 'unusable_empty_velocyto_layers'
status.write_text(json.dumps(manifest, indent=2) + '\n')
print({k: manifest[k] for k in ('status', 'exit_status', 'elapsed_seconds')})
raise SystemExit(0 if manifest['status'] == 'passed' else 1)
