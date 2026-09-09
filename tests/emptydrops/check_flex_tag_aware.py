#!/usr/bin/env python3
"""Compare the saved in-process fixture with one external call per sample."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--fixture', type=Path, required=True)
p.add_argument('--caller', type=Path, required=True)
p.add_argument('--historical-caller', type=Path)
a = p.parse_args()
root = a.fixture

def calls(path):
    return set(path.read_text().splitlines())

def ledger(path):
    with path.open() as f:
        return {r['barcode']: r for r in csv.DictReader(f, delimiter='\t')}

records = []
for label, tags in [('single', ['ACTTTAGG']), ('paired', ['ACAGTCTG', 'AGTGAGTG'])]:
    target = root/'external'/label
    target.mkdir(parents=True, exist_ok=False)
    common = ['--matrix', str(root/'matrix.mtx'), '--barcodes', str(root/'barcodes.tsv'),
        '--features', str(root/'features.tsv'), '--mitochondrial-genes', str(root/'mt.tsv'),
        '--flex-tag-aware', '--sim-n', '128', '--mc-threads', '2', '--fdr', '.01']
    for tag in tags:
        common += ['--barcode-tag', tag]
    cmd = [str(a.caller)] + common + ['--bootstrap-threads', '2', '--out', str(target/'calls.txt'), '--out-dir', str(target/'detail')]
    r = subprocess.run(cmd, env=dict(os.environ, OMP_NUM_THREADS='2'), capture_output=True, text=True)
    (target/'stdout.txt').write_text(r.stdout); (target/'stderr.txt').write_text(r.stderr)
    assert r.returncode == 0, (label, r.stderr)
    internal = root/'internal'/label
    assert calls(target/'calls.txt') == calls(internal/'calls.txt'), label
    before = ledger(internal/'EmptyDrops/emptydrops_results.tsv')
    after = ledger(target/'detail/EmptyDrops/emptydrops_results.tsv')
    assert before == after, (label, 'ledger mismatch')
    cfg = json.loads((internal/'caller_diagnostics.json').read_text())
    assert cfg['ambient_start_parameter'] == 45000*len(tags)
    assert cfg['ambient_base_end_parameter'] == 90000*len(tags)
    assert cfg['max_expected_cells'] == 22500*len(tags)
    assert cfg['ambient_umi_target'] == 0
    assert cfg['ambient_cells'] == cfg['input_cells'] - 45000*len(tags)
    assert cfg['bootstrap_threads'] == 2 and cfg['apply_bh'] == cfg['gate_on_fdr'] == 1
    assert cfg['candidate_umi_floor'] == 500
    assert any(int(row['umi_count']) == 500 for row in after.values())
    assert all(int(row['umi_count']) >= 500 for row in after.values() if row['is_simple_cell'] == '0')
    selected = {bc for bc, row in after.items() if row['is_simple_cell'] == '1' or row['passes_fdr'] == '1'}
    assert selected == calls(target/'calls.txt')
    assert sum(row['is_simple_cell'] == '1' for row in after.values()) == cfg['simple_cells']
    assert any(row['is_simple_cell'] == '0' and row['passes_fdr'] == '0' for row in after.values())
    records.append(dict(sample=label, tags=len(tags), calls=len(selected), candidates=len(after), exact_ledger_parity=True))
    if label == 'paired' and a.historical_caller:
        old = root/'historical'; old.mkdir(exist_ok=False)
        r = subprocess.run([str(a.historical_caller)]+common+['--out',str(old/'calls.txt'),'--out-dir',str(old/'detail')],
            env=dict(os.environ, OMP_NUM_THREADS='2'),capture_output=True,text=True)
        (old/'stdout.txt').write_text(r.stdout); (old/'stderr.txt').write_text(r.stderr)
        assert r.returncode == 0, r.stderr
        assert calls(old/'calls.txt') == selected
        prior = ledger(old/'detail/EmptyDrops/emptydrops_results.tsv')
        assert prior.keys() == after.keys()
        # The old writer falsely marked every candidate as non-simple. Compare
        # numerical evidence independently of that corrected diagnostic column.
        for bc in prior:
            for key in ['umi_count','p_value','p_adjusted','passes_raw_p','passes_fdr','obs_log_prob']:
                assert prior[bc][key] == after[bc][key], (bc, key)
        records[-1]['historical_numeric_and_call_parity'] = True
(root/'parity_results.json').write_text(json.dumps(records, indent=2))
print(json.dumps(records, indent=2))
