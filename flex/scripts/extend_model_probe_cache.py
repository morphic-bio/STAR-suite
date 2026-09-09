#!/usr/bin/env python3
"""Extend an existing H0/H1X2 cache to every probe in a public panel CSV.

The input gene-ID order is preserved. Missing active genes and DEPRECATED_
features are appended separately; deprecated features are never merged into
their active parents. Writes model_gene_ids.txt, included_gene_ids.txt,
model_h01x2_cache.bin, and a provenance manifest to a fresh directory.
No genome index, aligner, or external count matrix is used. Requires numpy.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import struct

import numpy as np

HEADER = struct.Struct('<8sHHIQ')
DT = np.dtype([('lo', '<u8'), ('hi', '<u8'), ('gene', '<u4'),
               ('cls', 'u1'), ('neg', 'u1'), ('sample', '<u2')])
KEY = np.dtype({'names': ['hi', 'lo', 'sample'],
                'formats': ['<u8', '<u8', '<u2'],
                'offsets': [8, 0, 22], 'itemsize': 24})
MASK64 = (1 << 64) - 1


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(16 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def packed(seq):
    require(len(seq) == 50 and set(seq) <= set('ACGT'),
            'Expected a 50-base A/C/G/T probe: ' + seq)
    value = 0
    for b in seq:
        value = (value << 2) | 'ACGT'.index(b)
    return value & MASK64, value >> 64


def ordered(data):
    a, b = data[:-1], data[1:]
    return np.all((a['hi'] < b['hi']) | ((a['hi'] == b['hi']) &
                  ((a['lo'] < b['lo']) | ((a['lo'] == b['lo']) &
                   (a['sample'] < b['sample'])))))


def extend(cache, gene_file, panel, out):
    require(not out.exists(), 'Output directory must be fresh: ' + str(out))
    genes = [s.strip() for s in gene_file.read_text().splitlines()
             if s.strip() and not s.startswith('#')]
    require(genes and len(genes) == len(set(genes)), 'Gene IDs must be unique')
    with panel.open() as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith('#')))
    require(bool(rows), 'Empty probe panel')
    for row in rows:
        packed(row['probe_seq'])
        require(row['included'].upper() in ('TRUE', 'FALSE'), 'Invalid included flag')
        require(row['region'] in ('spliced', 'unspliced'), 'Invalid probe region')
        require(bool(row['gene_id']), 'Empty probe gene ID')
        if 'DEPRECATED' in row['probe_id'].upper():
            require(row['gene_id'].startswith('DEPRECATED_'),
                    'Deprecated probes must have separate DEPRECATED_ feature IDs')
    panel_genes = {row['gene_id'] for row in rows}
    require(set(genes) <= panel_genes, 'Input gene IDs do not match this panel')
    added_genes = panel_genes - set(genes)
    genes += sorted(g for g in added_genes if not g.startswith('DEPRECATED_'))
    genes += sorted(g for g in added_genes if g.startswith('DEPRECATED_'))
    require(len(genes) <= 32767, 'Probe feature IDs exceed the 15-bit cache limit')
    gene_index = {g: i + 1 for i, g in enumerate(genes)}
    included = {r['gene_id'] for r in rows if r['included'].upper() == 'TRUE'
                and not r['gene_id'].startswith('DEPRECATED_')}
    require(bool(included), 'No included export features')
    parents = {}
    for row in rows:
        key = packed(row['probe_seq'])
        encoded = gene_index[row['gene_id']] | (
            {'spliced': 1, 'unspliced': 2}[row['region']] << 30)
        require(key not in parents or parents[key] == encoded,
                'Exact probe sequence has conflicting gene/region assignments')
        parents[key] = encoded
    with cache.open('rb') as f:
        magic, version, k, size, count = HEADER.unpack(f.read(24))
    require((magic, version, k, size) == (b'FH01SEQ1', 3, 50, 24),
            'Expected version-3 50-base Flex sequence cache')
    require(cache.stat().st_size == 24 + count * 24, 'Cache length mismatch')
    data = np.memmap(cache, mode='r', dtype=DT, offset=24, shape=(count,))
    old_parents, samples, class_counts = {}, set(), {}
    previous = None
    for start in range(0, count, 2000000):
        block = data[start:start + 2000000]
        require(ordered(block), 'Cache keys must be strictly sorted and unique')
        first = (int(block[0]['hi']), int(block[0]['lo']), int(block[0]['sample']))
        require(previous is None or previous < first, 'Cache block ordering mismatch')
        last = block[-1]
        previous = (int(last['hi']), int(last['lo']), int(last['sample']))
        classes, numbers = np.unique(block['cls'], return_counts=True)
        for cls, n in zip(classes, numbers):
            class_counts[int(cls)] = class_counts.get(int(cls), 0) + int(n)
        for rec in block[block['cls'] == 0]:
            key = int(rec['lo']), int(rec['hi'])
            gene = int(rec['gene']) & 0x7fff
            require(key in parents and parents[key] & 0x7fff == gene,
                    'Existing H0 parents do not match panel and input gene order')
            require(int(rec['sample']) > 0, 'Expected sample-specific H0 parents')
            old_parents[key] = gene
            samples.add(int(rec['sample']))
    require(samples and class_counts.get(4, 0) > 0, 'Input requires H0 and H1X2 records')
    added = {key: gene for key, gene in parents.items() if key not in old_parents}
    require(bool(added), 'All panel parents already exist; use the existing cache')
    variants = np.zeros(len(added) * 5775, dtype=DT)
    i = 0
    for (lo, hi), gene in added.items():
        value = (hi << 64) | lo
        changes = [[(ref ^ alt) << (2 * (49-p)) for alt in range(4) if alt != ref]
                   for p in range(50) for ref in [(value >> (2*(49-p))) & 3]]
        for p in range(50):
            for delta in changes[p]:
                v = value ^ delta
                variants[i] = (v & MASK64, v >> 64, gene, 4, 0, 0)
                i += 1
        for p in range(25):
            for delta in changes[p]:
                for q in range(25, 50):
                    for second in changes[q]:
                        v = value ^ delta ^ second
                        variants[i] = (v & MASK64, v >> 64, gene, 4, 0, 0)
                        i += 1
    require(i == len(variants), 'Variant generation count mismatch')
    variants.sort(order=['hi', 'lo', 'sample'])
    keys = variants.view(KEY)
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    multiplicity = np.diff(np.r_[starts, len(variants)])
    variants = variants[starts].copy()
    duplicate = multiplicity > 1
    variants['gene'][duplicate] = 0
    variants['cls'][duplicate] = 2
    variants['neg'][duplicate] = 1
    exact = np.array([(lo, hi, gene, 0, 0, s)
                      for (lo, hi), gene in added.items() for s in sorted(samples)], dtype=DT)
    additions = np.concatenate([variants, exact])
    additions.sort(order=['hi', 'lo', 'sample'])
    require(ordered(additions), 'Duplicate generated cache keys')
    positions = np.searchsorted(data.view(KEY), additions.view(KEY))
    match = np.zeros(len(additions), dtype=bool)
    valid = positions < count
    match[valid] = data.view(KEY)[positions[valid]] == additions.view(KEY)[valid]
    variant_match = match & (additions['sample'] == 0)
    require(np.all(np.isin(data['cls'][positions[variant_match]], [2, 4])),
            'Unexpected existing global class at a generated variant')
    additions['gene'][variant_match] = 0
    additions['cls'][variant_match] = 2
    additions['neg'][variant_match] = 1
    # A unique exact panel probe supersedes an old non-exact verdict.
    exact_match = match & (additions['sample'] != 0)
    require(np.all(data['cls'][positions[exact_match]] != 0), 'Conflicting exact parent')
    out.mkdir(parents=True)
    destination = out / 'model_h01x2_cache.bin'
    new_count = count + int((~match).sum())
    previous = None
    with destination.open('xb') as f:
        f.write(HEADER.pack(magic, version, k, size, new_count))
        for start in range(0, count, 2000000):
            end = min(start + 2000000, count)
            a = np.searchsorted(positions, start)
            b = np.searchsorted(positions, end, side='right' if end == count else 'left')
            block = np.array(data[start:end])
            part, pos, matches = additions[a:b], positions[a:b]-start, match[a:b]
            block[pos[matches]] = part[matches]
            insert = part[~matches]
            idx = pos[~matches] + np.arange(len(insert))
            merged = np.empty(len(block)+len(insert), dtype=DT)
            keep = np.ones(len(merged), dtype=bool)
            keep[idx] = False
            merged[idx], merged[keep] = insert, block
            require(ordered(merged), 'Output cache ordering failed')
            first = (int(merged[0]['hi']), int(merged[0]['lo']), int(merged[0]['sample']))
            require(previous is None or previous < first, 'Output block ordering failed')
            last = merged[-1]
            previous = (int(last['hi']), int(last['lo']), int(last['sample']))
            merged.tofile(f)
    require(destination.stat().st_size == 24 + new_count*24, 'Output length mismatch')
    (out/'model_gene_ids.txt').write_text('\n'.join(genes)+'\n')
    (out/'included_gene_ids.txt').write_text('\n'.join(g for g in genes if g in included)+'\n')
    report = dict(source=str(cache), source_sha256=sha(cache),
                  source_gene_ids=str(gene_file), source_gene_ids_sha256=sha(gene_file),
                  panel=str(panel), panel_sha256=sha(panel),
                  output=str(destination), output_sha256=sha(destination),
                  source_records=count, output_records=new_count,
                  source_class_counts=class_counts, sample_indices=sorted(samples),
                  existing_parent_probes=len(old_parents), added_parent_probes=len(added),
                  generated_variants=i, ambiguous_new_variants=int(duplicate.sum()),
                  existing_variant_collisions=int(variant_match.sum()),
                  exact_overrides_nonexact=int(exact_match.sum()),
                  model_features=len(genes), export_features=len(included),
                  deprecated_features=sum(g.startswith('DEPRECATED_') for g in genes),
                  policy='All panel probes for calling; included active genes for export; '
                         'separate deprecated features; preserve old feature indices; '
                         'H1X2 multiple-parent variants are DENY; no alignment or external counts')
    (out/'manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True, type=Path)
    parser.add_argument('--gene-ids', required=True, type=Path)
    parser.add_argument('--probe-csv', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path)
    args = parser.parse_args()
    extend(args.cache, args.gene_ids, args.probe_csv, args.out_dir)
