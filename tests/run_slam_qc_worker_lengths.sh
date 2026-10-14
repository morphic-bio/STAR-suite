#!/usr/bin/env bash
# Constructed inputs: idle workers must not supply SLAM QC mate lengths.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
star=${STAR_BIN:-$root/core/legacy/source/STAR}
threads=${THREADS:-1}
work=${OUT_BASE:-$(mktemp -d /tmp/slam-qc-worker-lengths.XXXXXX)}
mkdir -p "$work"
python3 - "$star" "$threads" "$work" <<'PY'
from pathlib import Path
import json, os, random, subprocess, sys
binary,threads,work=sys.argv[1:]; root=Path(work)
rng=random.Random(614);seq=''.join(rng.choice('ACGT') for _ in range(5000))
(root/'genome.fa').write_text('>chrSynthetic\n'+seq+'\n')
(root/'genes.gtf').write_text('chrSynthetic\tconstructed\texon\t101\t4900\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "G1";\n')
index=root/'index';index.mkdir(exist_ok=False)
cmd=[binary,'--runMode','genomeGenerate','--runThreadN','1','--genomeDir',str(index),'--genomeFastaFiles',str(root/'genome.fa'),'--sjdbGTFfile',str(root/'genes.gtf'),'--genomeSAindexNbases','3','--genomeChrBinNbits','8','--outFileNamePrefix',str(root/'index_')]
with (root/'index.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
for case,reads,expected in [('empty',('', ''),(100,0)),('pair',(seq[200:280],seq[260:320].translate(str.maketrans('ACGT','TGCA'))[::-1]),(80,60)),('single',(seq[200:270],),(70,0))]:
 paths=[]
 for m,read in enumerate(reads):
  p=root/f'{case}-R{m+1}.fastq';p.write_text(f'@synthetic/{m+1}\n{read}\n+\n'+ 'I'*len(read)+'\n' if read else '');paths.append(str(p))
 for fill in ['85','170']:
  d=root/f'{case}-{fill}';d.mkdir(exist_ok=False)
  cmd=[binary,'--runThreadN',threads,'--genomeDir',str(index),'--readFilesIn',*paths,'--outFileNamePrefix',str(d/'star_'),'--outSAMtype','None','--slamQuantMode','1','--slamQcReport',str(d/'qc')]
  (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
  with (d/'driver.log').open('w') as log:subprocess.run(cmd,env=dict(os.environ,MALLOC_PERTURB_=fill),check=True,stdout=log,stderr=subprocess.STDOUT)
  qc=json.loads((d/'qc.slam_qc.json').read_text());got=(qc['mate_len_mate1'],qc['mate_len_mate2'])
  assert got==expected,(case,fill,got,expected)
  print(case,fill,got)
(root/'PASS').write_text('All empty/paired/single input mate lengths match at both allocator fills.\n')
PY
