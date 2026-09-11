#!/usr/bin/env python3
"""Real BGZF ordered input, early cancellation and malformed input with one permit."""
import argparse
import hashlib
import json
import subprocess
import time
import zlib
from pathlib import Path
from run_pf_bgzf_tests import bgzf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = [b'', b'']
    crc = 0
    prefix_crc = 0
    for n in range(513):
        for mate in range(2):
            seq = (('ACGT'[(n+mate) % 4] + 'ACGT' * 20).encode())
            qual = b'I' * len(seq)
            raw[mate] += f'@read{n}/{mate+1}\n'.encode() + seq + b'\n+\n' + qual + b'\n'
            crc = zlib.crc32(seq + qual, crc)
        if n == 16:
            prefix_crc = crc
    inputs = {
        'R1.gz': bgzf(raw[0], 211),
        'R2.gz': bgzf(raw[1], 379),
        'short.gz': bgzf(b'\n'.join(raw[1].split(b'\n')[:-5]) + b'\n', 239),
        'names.gz': bgzf(raw[1].replace(b'@read17/2', b'@other17/2'), 239),
    }
    bad = bytearray(inputs['R2.gz'])
    first_size = int.from_bytes(bad[16:18], 'little') + 1
    bad[first_size - 8] ^= 1
    inputs['crc.gz'] = bytes(bad)
    manifest = {}
    for name, data in inputs.items():
        (args.out / name).write_bytes(data)
        manifest[name] = {'size': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    (args.out / 'inputs.json').write_text(json.dumps(manifest, indent=2))
    binary_hash = hashlib.sha256(args.binary.read_bytes()).hexdigest()
    cases = [('ordered', 'R2.gz', 0, 0, 513, crc),
             ('cancel', 'R2.gz', 17, 0, 17, prefix_crc),
             ('short', 'short.gz', 0, 1, None, None),
             ('names', 'names.gz', 0, 1, None, None),
             ('crc', 'crc.gz', 0, 1, None, None)]
    for name, mate, limit, expected_exit, count, expected_crc in cases:
        argv = [str(args.binary), str(args.out / 'R1.gz'), str(args.out / mate),
                str(limit), str(args.out / name)]
        record = {'argv': argv, 'binary_sha256': binary_hash, 'status': 'running',
                  'started_unix': time.time()}
        path = args.out / f'{name}.execution.json'
        path.write_text(json.dumps(record, indent=2))
        with (args.out / f'{name}.stdout').open('w') as out, (args.out / f'{name}.stderr').open('w') as err:
            try:
                result = subprocess.run(argv, stdout=out, stderr=err, timeout=30)
                record.update(status='complete', exit_status=result.returncode)
            except subprocess.TimeoutExpired:
                record.update(status='timeout')
                raise
            finally:
                record['finished_unix'] = time.time()
                path.write_text(json.dumps(record, indent=2))
        stdout = (args.out / f'{name}.stdout').read_text()
        assert result.returncode == expected_exit, (name, result.returncode)
        assert 'balanced=1' in stdout, (name, stdout)
        if count is not None:
            assert f'pairs={count} ' in stdout and f'crc={expected_crc:08x}' in stdout, stdout
        print(name, 'PASS')


if __name__ == '__main__':
    main()
