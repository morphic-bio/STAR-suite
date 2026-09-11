from pathlib import Path
import argparse,subprocess,json,hashlib,time
parser=argparse.ArgumentParser(description="Compare one Flex MEX loader build against a preserved ledger")
parser.add_argument('--out', required=True, type=Path)
parser.add_argument('--source-dir', type=Path)
parser.add_argument('--expected-ledger', type=Path)
parser.add_argument('--sanitize', action='store_true')
parser.add_argument('--fixtures', required=True, type=Path)
args=parser.parse_args()
repo=Path(__file__).resolve().parents[1]
src=args.source_dir.resolve() if args.source_dir else repo/'flex/source/libflex'
out=args.out.resolve();out.mkdir(exist_ok=False)
inc=[src,repo/'core/legacy/source',repo/'flex/source',repo/'flex/source/libflex',repo/'core/features/libscrna/include',repo/'core/features/process_features/include',repo/'core/features/vbem/source',repo/'core/features/vbem/source/libem',repo/'core/features/bamsort/source',repo/'slam/source',repo/'core/legacy/source/htslib']
flags=['-std=c++11','-O2','-g','-fopenmp','-ffunction-sections','-fdata-sections']+(['-fsanitize=address,undefined','-fno-sanitize-recover=all','-fno-omit-frame-pointer'] if args.sanitize else [])
cmd=['g++',*flags,*['-I'+str(p) for p in inc],'-DLOADER_SOURCE="'+str(src/'FlexFilter.cpp')+'"',str(repo/'tests/test_flex_mex_storage.cpp'),'-Wl,--gc-sections','-o',str(out/'probe')]
run=[str(out/'probe'),*[str(p) for p in sorted(args.fixtures.resolve().iterdir()) if p.is_dir()]]
rec={'status':'building','commands':[cmd],'source_sha256':hashlib.sha256((src/'FlexFilter.cpp').read_bytes()).hexdigest(),'started_unix':time.time()};p=out/'execution.json'
def save():p.write_text(json.dumps(rec,indent=2)+'\n')
save()
with (out/'build.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
if r.returncode:rec.update(status='build_failed',exit_status=r.returncode);save();raise SystemExit(r.returncode)
rec['commands'].append(run);rec['status']='running';save()
with (out/'ledger.tsv').open('w') as f,(out/'stderr.log').open('w') as e:r=subprocess.run(run,stdout=f,stderr=e)
rec.update(status='passed' if r.returncode==0 else 'failed',exit_status=r.returncode,finished_unix=time.time(),ledger_sha256=hashlib.sha256((out/'ledger.tsv').read_bytes()).hexdigest())
if args.expected_ledger:
 rec['equal']=rec['ledger_sha256']==hashlib.sha256(args.expected_ledger.read_bytes()).hexdigest()
 if not rec['equal']:rec['status']='mismatch'
save();print({k:rec[k] for k in ('status','exit_status','ledger_sha256')});raise SystemExit(0 if rec['status']=='passed' else 1)
