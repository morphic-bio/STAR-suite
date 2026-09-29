#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PF="$ROOT/core/features/process_features"
OUT="$(mktemp -d /tmp/celltag_ambient_test.XXXXXX)"
trap 'rm -rf "$OUT"' EXIT
cat > "$OUT/call.c" <<'C'
#include "call_features.h"
int main(int argc, char **argv) {
    if(argc!=4) return 2;
    cf_ambient_fdr_config *cfg=cf_ambient_fdr_config_create();
    if(!cfg) return 2;
    cfg->feature_type="CellTag";
    int ret=cf_process_mex_dir_ambient_fdr(argv[1],argv[2],argv[3],cfg);
    cf_ambient_fdr_config_destroy(cfg);
    return ret!=0;
}
C
cc -I"$PF/include" -c "$OUT/call.c" -o "$OUT/call.o"
g++ -o "$OUT/call" "$OUT/call.o" "$PF/libprocess_features.a" "$ROOT/core/features/libscrna/libscrna.a" -lm -lpthread -lz -fopenmp -lhts
python3 - "$OUT" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])
features='guideA\tguideA\tCRISPR Guide Capture\ntagA\ttagA\tCellTag\ntagB\ttagB\tCellTag\n'
for name in ('raw','cells'):
 path=root/name;path.mkdir()
 (path/'features.tsv').write_text(features)
 cells=['cell1','cell2']+(['empty1','empty2'] if name=='raw' else [])
 (path/'barcodes.tsv').write_text('\n'.join(cells)+'\n')
 values=[(1,1,100),(2,1,25),(3,2,1)]
 if name=='raw': values += [(1,3,100),(2,3,1),(3,3,99),(2,4,1),(3,4,99)]
 (path/'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n'+f'3 {len(cells)} {len(values)}\n'+''.join(f'{r} {c} {v}\n' for r,c,v in values))
PY
"$OUT/call" "$OUT/raw" "$OUT/cells" "$OUT/calls"
python3 - "$OUT" <<'PY'
from pathlib import Path
import csv,json,sys
path=Path(sys.argv[1])/'calls'
features=(path/'guide_qvalues_features.tsv').read_text()
assert 'CRISPR' not in features and features.count('CellTag')==2
calls=list(csv.DictReader((path/'guide_fdr_calls_per_cell.csv').open()))
assert calls[0]['feature_call']=='tagA' and calls[0]['num_features']=='1'
assert calls[1]['num_features']=='0'
summary=json.loads((path/'guide_fdr_summary.json').read_text())
assert summary['ambient_total_umis']==200
PY
echo 'PASS: explicit CellTag selection excludes real guides and preserves feature types'
