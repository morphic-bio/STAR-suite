#!/usr/bin/env python3
"""Compare current barcode storage with a preserved pre-change CbCorrector source."""
import argparse, hashlib, json, re, shutil, subprocess, time
from pathlib import Path

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--baseline-dir',type=Path,required=True,help='Directory containing solo/CbCorrector.{h,cpp}')
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--sanitize',action='store_true')
    args=ap.parse_args(); repo=Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True,exist_ok=False)
    legacy=args.out/'legacy/solo';legacy.mkdir(parents=True)
    hashes={}
    for ext in ('h','cpp'):
        src=args.baseline_dir/'solo'/('CbCorrector.'+ext)
        hashes[str(src)]=sha(src)
        text=src.read_text().replace('CODE_CbCorrector','CODE_LegacyCbCorrector')
        text=re.sub(r'\bCbCorrector\b','LegacyCbCorrector',text)
        text=re.sub(r'\bCbMatch\b','LegacyCbMatch',text)
        (legacy/('LegacyCbCorrector.'+ext)).write_text(text)
    sources=[repo/'tests/test_cb_corrector_storage.cpp',repo/'flex/source/solo/CbCorrector.cpp',legacy/'LegacyCbCorrector.cpp']
    flags=['-std=c++11','-O2','-g','-pthread','-I'+str(repo/'flex/source'),'-I'+str(repo/'core/legacy/source'),'-I'+str(repo/'core/legacy/source/htslib'),'-I'+str(args.out/'legacy')]
    if args.sanitize:flags+=['-fsanitize=address,undefined','-fno-omit-frame-pointer','-fno-sanitize-recover=all']
    binary=args.out/'test'
    command=['g++']+flags+[str(p) for p in sources]+['-o',str(binary)]
    record={'status':'building','started_unix':time.time(),'command':command,'baseline_hashes':hashes,'current_hashes':{str(p):sha(p) for p in sources[:2]+[repo/'flex/source/solo/CbCorrector.h']},'sanitizers':args.sanitize}
    manifest=args.out/'execution.json'
    def save():manifest.write_text(json.dumps(record,indent=2)+'\n')
    save()
    with (args.out/'build.log').open('w') as f:r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT)
    if r.returncode:record.update(status='build_failed',exit_status=r.returncode);save();return r.returncode
    record.update(status='running',binary_sha256=sha(binary));save()
    with (args.out/'stdout.log').open('w') as out,(args.out/'stderr.log').open('w') as err:
        r=subprocess.run([str(binary)],stdout=out,stderr=err)
    record.update(status='passed' if r.returncode==0 else 'failed',exit_status=r.returncode,finished_unix=time.time());save()
    print(json.dumps({k:record[k] for k in ('status','exit_status')}))
    return r.returncode
if __name__=='__main__':raise SystemExit(main())
