#!/usr/bin/env python3
"""Simulate our H0/H1 hash verdict on every R2 read of a FASTQ sample using the real cache, aggregate per molecule."""
import sys, gzip, itertools, csv, numpy as np
from collections import defaultdict, Counter
CACHE="/home/lhhung/jax_stage_20260903/ref/h01_cache.bin"; PS="/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv"
FQ="/home/lhhung/jax_stage_20260903/fastq/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_"; NPAIRS=int(sys.argv[1]) if len(sys.argv)>1 else 20_000_000
hdr=np.dtype([("magic","S8"),("version","<u2"),("k","<u2"),("recsize","<u4"),("n","<u8")]); rec=np.dtype([("lo","<u8"),("hi","<u8"),("gene","<u4"),("cls","u1"),("neg","u1"),("res","<u2")])
R=np.fromfile(CACHE,dtype=rec,offset=hdr.itemsize)
keys=(R["hi"].astype(object)<<64)|R["lo"].astype(object)   # python ints
verdict={}   # key -> (cls, neg, gene)
for k,c,ng,g in zip(keys.tolist(),R["cls"].tolist(),R["neg"].tolist(),(R["gene"]&0x7FFF).tolist()):
    v=verdict.get(k)
    if v is None or (v[1]==0 and ng): verdict[k]=(c,ng,g)   # any DENY record for the key denies (H0 records are per-sample copies)
print(f"cache dict {len(verdict):,} keys", flush=True)
tbl=str.maketrans("ACGT","0123")
probes=[]; lhs=defaultdict(list); rhs=defaultdict(list); gene_of=[]
for row in csv.reader(l for l in open(PS) if not l.startswith("#")):
    if row[0]=="gene_id" or row[3]!="TRUE": continue
    i=len(probes); probes.append((row[2],row[1],row[5])); lhs[row[1][:25]].append(i); rhs[row[1][25:50]].append(i)
seqs=[p[1] for p in probes]
def hd(a,b): return sum(x!=y for x,y in zip(a,b))
def tolerant(w):
    c=set(lhs.get(w[:25],()))|set(rhs.get(w[25:],()))
    if not c: return (99,-1)
    best=min(c,key=lambda i:hd(w,seqs[i])); return (hd(w,seqs[best]),best)
# per read: K (keep H0/H1), D (deny), M2 (miss, HD2-5 to a probe -> aligner would keep), N, J (junk: no probe within 5)
mols={}; rc=Counter()
with gzip.open(FQ+"R1_001.fastq.gz","rt") as f1, gzip.open(FQ+"R2_001.fastq.gz","rt") as f2:
    for n,(l1,l2) in enumerate(zip(itertools.islice(f1,1,None,4),itertools.islice(f2,1,None,4))):
        if n>=NPAIRS: break
        w=l2[:50]; k=(l1[:16],l1[16:28],l2[68:76]); best=-1
        if "N" in w or len(w)<50: cls="N"
        else:
            v=verdict.get(int(w.translate(tbl),4))
            if v is None:
                d,best=tolerant(w); cls="M2" if 2<=d<=5 else "J"
            else: cls="D" if v[1] else "K"; 
            if cls=="K": best=v[2]
        rc[cls]+=1
        m=mols.get(k)
        if m is None: mols[k]=[cls, best, 1]
        else:
            m[2]+=1
            rank={"K":0,"M2":1,"D":2,"N":3,"J":4}
            if rank[cls]<rank[m[0]]: m[0]=cls; m[1]=best
tot=sum(rc.values()); print(f"{tot:,} reads: "+"  ".join(f"{k}={v/tot:.3%}" for k,v in rc.items()))
for lab,sel in (("multi-read",lambda m:m[2]>=2),("singleton",lambda m:m[2]==1)):
    ms=[m for m in mols.values() if sel(m)]; c=Counter(m[0] for m in ms); n=len(ms)
    print(f"{lab} molecules {n:,}: best verdict  KEEP {c['K']/n:.3%}  MISS-but-aligner-keeps(HD2-5) {c['M2']/n:.3%}  DENY-only {c['D']/n:.3%}  N-only {c['N']/n:.3%}  junk {c['J']/n:.3%}")
    if lab=="multi-read":
        for cl in ("M2","D"):
            genes=Counter(probes[m[1]][2] for m in ms if m[0]==cl and m[1]>=0); base=Counter()
            for m in ms:
                if m[1]>=0: base[probes[m[1]][2]]+=1
            print(f"  genes with most {cl} multi-read molecules (lost/all): "+", ".join(f"{g} {c}/{base[g]} ({c/base[g]:.0%})" for g,c in genes.most_common(12)))
