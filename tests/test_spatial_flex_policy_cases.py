#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,hashlib
import argparse
p=argparse.ArgumentParser(description="Hand-constructed half-probe decisions through the spatial FASTQ adapter")
p.add_argument('--template', type=Path, required=True, help='Accepted small spatial attempt.json with argv')
p.add_argument('--cache', type=Path, required=True, help='pair_fixture.bin.half.khash from test_flex_khash_storage.sh')
p.add_argument('--r1-source', type=Path, required=True, help='Small plain R1 fixture')
p.add_argument('--out', type=Path, required=True)
aargs=p.parse_args()
out=aargs.out;out.mkdir(parents=True,exist_ok=False)
a,c,g,t=('A'*25,'C'*25,'G'*25,'T'*25)
def mutate(s,changes):
 s=list(s)
 for i,b in changes:s[i]=b
 return ''.join(s)
fourth=mutate(a+c,[(0,'C'),(25,'G')])
queries=[a+c,mutate(a+g,[(10,'C')]),mutate(a+g,[(30,'T')]),mutate(a+g,[(10,'C'),(30,'T')]),mutate(a+g,[(i,'T') for i in range(9)]),mutate(a+g,[(10,'N')]),mutate(a+c,[(0,'G'),(25,'A')]),t+g,mutate(a+g,[(i,'T') for i in range(11)]),'A'*10,'N'*50,fourth]
raw=aargs.r1_source.read_text().splitlines()
r1s=[seq for seq in raw[1::4] if len(seq)>=39 and set(seq[:9])<=set('ACGT')][:len(queries)]
assert len(r1s)==len(queries);r1s[3]=mutate(r1s[3],[(12,'N')]);r1s[4]=mutate(r1s[4],[(0,'N')])
for mate,seqs in [('r1',r1s),('r2',queries)]:
 # Identical names intentionally do not identify these independent reads.
 (out/(mate+'.fastq')).write_text(''.join('@repeated_name\n'+s+'\n+\n'+'I'*len(s)+'\n' for s in seqs))
(out/'genes.txt').write_text('ENSG_TEST_1\nENSG_TEST_2\nENSG_TEST_3\n')
old=json.loads(aargs.template.read_text())
cmd=list(old['argv'])
for key,value in {'--runThreadN':'4','--outFileNamePrefix':str(out/'star')+'/','--soloHashScreenFile':str(aargs.cache),'--soloProbeList':str(out/'genes.txt'),'--soloSpatialExpectedReads':str(len(queries)),'--soloSpatialExpectedCandidates':'1000'}.items():cmd[cmd.index(key)+1]=value
j=cmd.index('--readFilesIn');cmd[j+1:j+3]=[str(out/'r2.fastq'),str(out/'r1.fastq')]
p=out/'attempt.json';r={'argv':cmd,'purpose':'hand-constructed shared half-probe policy through actual spatial FASTQ adapter, including repeated names and barcode/UMI N','expected':{'feature_hash_h0':2,'feature_hash_h1':1,'feature_hash_h1x2':4,'feature_hash_deny':5,'feature_hash_miss':0,'feature_assigned_reads':7,'reads_decoded':12},'short_read_policy':'shared seed classifier reports terminal HalfNoAnchor denial; no alignment handoff', 'test_only_cache':'four synthetic H0 parents; production half tables built by the storage unit test, not a biological reference','status':'running'};p.write_text(json.dumps(r,indent=2)+'\n')
with (out/'console.log').open('w') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
r['exit_code']=rc
if rc==0:
 s=dict(line.split('\t',1) for line in (out/'star/SpatialFlex.out/run_summary.tsv').read_text().splitlines())
 r['observed']={key:int(s[key]) for key in r['expected']}
 if r['observed']!=r['expected']:rc=1
r['status']='complete' if rc==0 else 'failed';p.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r.get('observed',{})),r['status']);raise SystemExit(rc)
