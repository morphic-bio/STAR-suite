import gzip, itertools, csv, numpy as np
from collections import defaultdict, Counter
FQ="/home/lhhung/jax_stage_20260903/fastq/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_"
PS="/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv"
GENES={"ZIC2","H2BC15","N4BP2","APOE","SOX2","ACTB"}
probes=[]; lhs=defaultdict(list); rhs=defaultdict(list)
for row in csv.reader(l for l in open(PS) if not l.startswith("#")):
    if row[0]=="gene_id" or row[3]!="TRUE" or row[5] not in GENES: continue
    i=len(probes); probes.append((row[2].split("|")[-1],row[1],row[5])); lhs[row[1][:25]].append(i); rhs[row[1][25:50]].append(i)
hist=defaultdict(Counter); mmpos=defaultdict(Counter); nreads=Counter()
with gzip.open(FQ+"R2_001.fastq.gz","rt") as f2:
    for n,l2 in enumerate(itertools.islice(f2,1,None,4)):
        if n>=20_000_000: break
        w=l2[:50]
        c=set(lhs.get(w[:25],()))|set(rhs.get(w[25:],()))
        if not c: continue
        best=min(c,key=lambda i:sum(x!=y for x,y in zip(w,probes[i][1]))); s=probes[best][1]
        d=[k for k,(x,y) in enumerate(zip(w,s)) if x!=y]
        key=(probes[best][2],probes[best][0]); hist[key][min(len(d),9)]+=1; nreads[key]+=1
        if 1<=len(d)<=5:
            for k in d: mmpos[key][k]+=1
for key in sorted(hist):
    g,p=key; n=nreads[key]; h=hist[key]
    top=mmpos[key].most_common(4)
    print(f"{g:7s} {p}: reads {n:>7,}  HD0 {h[0]/n:.1%} HD1 {h[1]/n:.1%} HD2 {h[2]/n:.1%} HD3 {h[3]/n:.1%} HD4-5 {(h[4]+h[5])/n:.1%} HD>=6 {sum(v for k,v in h.items() if k>=6)/n:.1%}  | top mismatch positions (pos:count): "+", ".join(f"{k}:{v}" for k,v in top))
