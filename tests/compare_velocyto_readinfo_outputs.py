"""Check all GEX/guide counts and Velocyto layers for the readInfo regression."""
import argparse
import csv
import json
from pathlib import Path
from compare_permit_hierarchy_outputs import compare_mex, content, digest, locate

p = argparse.ArgumentParser()
p.add_argument('--left', type=Path, required=True)
p.add_argument('--right', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--expect-restored-layers', action='store_true',
               help='Require empty old layers, positive new layers and unchanged GEX/guide counts')
a = p.parse_args()
for run in (a.left, a.right):
    execution = json.loads((run / 'execution.json').read_text())
    assert execution['exit_status'] == 0 and execution['status'] == 'passed'
roots = [r / 'out' for r in (a.left, a.right)]
dirs = [{m.parent.relative_to(root) for m in root.rglob('matrix.mtx*')
         if not a.expect_restored_layers or 'velocyto' not in str(m.relative_to(root)).lower()}
        for root in roots]
assert dirs[0] == dirs[1] and dirs[0], 'Matrix surfaces differ'
report = {'left': str(a.left), 'right': str(a.right), 'matrices': {}, 'guides': {}, 'layers': {}}
for d in sorted(dirs[0]):
    report['matrices'][str(d)] = compare_mex(roots[0] / d, roots[1] / d)
guides = [{f.relative_to(root) for f in (root / 'outs/crispr_analysis').glob('*.csv')} for root in roots]
assert guides[0] == guides[1] and guides[0]
for f in sorted(guides[0]):
    tables = []
    for root in roots:
        with (root / f).open() as stream:
            rows = list(csv.reader(stream))
            tables.append((rows[:1], sorted(rows[1:])))
    report['guides'][str(f)] = {'equal': tables[0] == tables[1]}
for kind in ('raw', 'filtered'):
    layer_roots = [root / 'Solo.out/Velocyto' / kind for root in roots]
    for layer in ('spliced', 'unspliced', 'ambiguous'):
        files = [locate(root, layer + '.mtx') for root in layer_roots]
        sums = []
        for f in files:
            with content(f) as stream:
                lines = (line for line in stream if line.strip() and not line.startswith(b'%'))
                next(lines)  # dimensions; stored rows may include zero counts
                sums.append(sum(int(line.split()[2]) for line in lines))
        report['layers'][kind + '/' + layer] = {
            'equal': digest(files[0]) == digest(files[1]), 'umi_sums': sums}
    if not a.expect_restored_layers:
        for name in ('features.tsv', 'barcodes.tsv'):
            report['layers'][kind + '/' + name] = {
                'equal': digest(locate(layer_roots[0], name)) == digest(locate(layer_roots[1], name))}
report['equal_counts'] = all(v['equal'] for group in ('matrices', 'guides') for v in report[group].values())
if a.expect_restored_layers:
    report['passed'] = report['equal_counts'] and all(
        report['layers']['raw/' + layer]['umi_sums'][0] == 0 and
        report['layers']['raw/' + layer]['umi_sums'][1] > 0
        for layer in ('spliced', 'unspliced', 'ambiguous'))
else:
    report['passed'] = report['equal_counts'] and all(v['equal'] for v in report['layers'].values())
a.out.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'passed': report['passed'], 'matrices': len(report['matrices']),
                  'guide_tables': len(report['guides']), 'layers': report['layers']}))
raise SystemExit(0 if report['passed'] else 1)
