#!/usr/bin/env python3
"""Small cache integration fixture: feature identity, collisions, and export."""
import csv
import importlib.util
from pathlib import Path
import tempfile

import numpy as np

path = Path(__file__).resolve().parents[2]/'flex/scripts/extend_model_probe_cache.py'
spec = importlib.util.spec_from_file_location('model_cache', path)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)

with tempfile.TemporaryDirectory(prefix='star_model_probe_cache_') as tmp:
    r = Path(tmp)
    active = 'ACGT'*12+'AC'
    deprecated = 'TGCATGCA'*6+'TG'
    excluded = 'A'+deprecated[1:]
    variant = deprecated[:25]+'A'+deprecated[26:]
    rows = [('ACTIVE', active, 'ACTIVE|a', 'TRUE', 'spliced'),
            ('EXCLUDED', excluded, 'EXCLUDED|b', 'FALSE', 'unspliced'),
            ('DEPRECATED_ACTIVE', deprecated, 'DEPRECATED_ACTIVE|c', 'FALSE', 'spliced')]
    with (r/'panel.csv').open('w') as f:
        w = csv.writer(f);w.writerow(['gene_id','probe_seq','probe_id','included','region']);w.writerows(rows)
    (r/'genes.txt').write_text('ACTIVE\n')
    def rec(seq, gene, cls, sample):
        lo, hi = builder.packed(seq)
        return (lo, hi, gene | (1<<30), cls, 0, sample)
    data = np.array([rec(active,1,0,s) for s in [1,2]] +
                    [rec(deprecated,1,1,s) for s in [1,2]] +
                    [rec(variant,1,4,0)], dtype=builder.DT)
    data.sort(order=['hi','lo','sample'])
    with (r/'input.bin').open('wb') as f:
        f.write(builder.HEADER.pack(b'FH01SEQ1',3,50,24,len(data)));data.tofile(f)
    old_bytes = (r/'input.bin').read_bytes()
    report = builder.extend(r/'input.bin',r/'genes.txt',r/'panel.csv',r/'out')
    assert (r/'input.bin').read_bytes() == old_bytes
    assert (r/'out/model_gene_ids.txt').read_text().splitlines() == ['ACTIVE','EXCLUDED','DEPRECATED_ACTIVE']
    assert (r/'out/included_gene_ids.txt').read_text() == 'ACTIVE\n'
    output = np.fromfile(r/'out/model_h01x2_cache.bin',dtype=builder.DT,offset=24)
    assert builder.ordered(output)
    assert report['exact_overrides_nonexact'] == 2 and report['existing_variant_collisions'] == 1
    assert report['ambiguous_new_variants'] > 0
    assert report['model_features'] == 3 and report['deprecated_features'] == 1
    lookup = {(int(a['lo']),int(a['hi']),int(a['sample'])):a for a in output}
    for seq, gene in [(active,1),(excluded,2),(deprecated,3)]:
        for sample in [1,2]:
            v = lookup[(*builder.packed(seq),sample)]
            assert v['cls'] == 0 and int(v['gene']) & 0x7fff == gene
    denied = lookup[(*builder.packed(variant),0)]
    assert denied['cls'] == 2 and denied['neg'] == 1 and denied['gene'] == 0
    try:
        builder.extend(r/'input.bin',r/'genes.txt',r/'panel.csv',r/'out')
    except ValueError:
        pass
    else:
        raise AssertionError('Existing outputs must never be overwritten')
print('PASS: model features stay separate, export mask, exact priority, variant ambiguity, source preservation')
