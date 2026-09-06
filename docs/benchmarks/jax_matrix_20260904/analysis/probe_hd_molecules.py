#!/usr/bin/env python3
"""Per-read probe Hamming distance with a tolerant matcher (25-mer halves, then 1-mismatch halves),
aggregated per molecule (CB,UMI,tag). Question: how many molecules have NO read within HD1 of a probe
(our hash misses them) but a read that Cell Ranger's aligner would still place (HD<=5, both halves)?"""
import sys, gzip, itertools, csv, numpy as np
from collections import defaultdict, Counter
FQ="/home/lhhung/jax_stage_20260903/fastq/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_"
PS="/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv"
NPAIRS=int(sys.argv[1]) if len(sys.argv)>1 else 20_000_000
probes=[]; lhs=defaultdict(list); rhs=defaultdict(list)
for row in csv.reader(l for l in open(PS) if not l.startswith("#")):
    if row[0]=="gene_id" or row[3]!="TRUE": continue
    i=len(probes); probes.append((row[2],row[1],row[5])); lhs[row[1][:25]].append(i); rhs[row[1][25:50]].append(i)
seqs=[p[1] for p in probes]; gc=np.array([(s.count("G")+s.count("C"))/50 for s in seqs])
def hd(a,b): return sum(x!=y for x,y in zip(a,b))
def variants(s):
    for i in range(len(s)):
        for b in "ACGT":
            if b!=s[i]: yield s[:i]+b+s[i+1:]
def classify(r):
    w=r[:50]
    if len(w)<50: return ("short",-1)
    if "N" in w: return ("N",-1)
    c=set(lhs.get(w[:25],()))|set(rhs.get(w[25:],()))
    if not c:
        for v in variants(w[:25]): c.update(lhs.get(v,()))
        if not c:
            for v in variants(w[25:]): c.update(rhs.get(v,()))
    if not c: return ("none",-1)
    best=min(c,key=lambda i:hd(w,seqs[i])); d=hd(w,seqs[best])
    return (d,best)
mols={}  # (cb,umi,tag) -> [nreads, minHD, bestprobe, nHDle1]
readcls=Counter()
with gzip.open(FQ+"R1_001.fastq.gz","rt") as f1, gzip.open(FQ+"R2_001.fastq.gz","rt") as f2:
    for n,(l1,l2) in enumerate(zip(itertools.islice(f1,1,None,4),itertools.islice(f2,1,None,4))):
        if n>=NPAIRS: break
        cb=l1[:16]; umi=l1[16:28]; tag=l2[68:76]
        cls,best=classify(l2.rstrip())
        readcls[cls if isinstance(cls,str) else min(cls,6)]+=1
        k=(cb,umi,tag); m=mols.get(k)
        d=cls if not isinstance(cls,str) else 99
        if m is None: mols[k]=[1,d,best,int(d<=1)]
        else:
            m[0]+=1; m[3]+=int(d<=1)
            if d<m[1]: m[1]=d; m[2]=best
tot=sum(readcls.values()); print(f"{tot:,} reads; per-read probe class:", {k:f"{v/tot:.3%}" for k,v in sorted(readcls.items(),key=lambda x:str(x[0]))})
multi=[m for m in mols.values() if m[0]>=2]; single=[m for m in mols.values() if m[0]==1]
print(f"molecules {len(mols):,}; multi-read {len(multi):,}; singletons {len(single):,}")
def summarize(label,ms):
    n=len(ms); lost=[m for m in ms if m[3]==0]; recov=[m for m in lost if 2<=m[1]<=5]
    print(f"  {label}: no read within HD1 (we miss the molecule): {len(lost):,} = {len(lost)/n:.3%}; of which best read HD2-5 (CR aligner keeps): {len(recov):,} = {len(recov)/n:.3%}; HD2 exactly {sum(1 for m in lost if m[1]==2)/n:.3%}; N/none/HD>5 {sum(1 for m in lost if m[1]>5)/n:.3%}")
    return recov
r_multi=summarize("multi-read molecules",multi); r_single=summarize("singletons",single)
# GC dependence and gene attribution of recoverable multi-read losses
allgc=np.array([gc[m[2]] for m in multi if m[2]>=0]); lostgc=np.array([gc[m[2]] for m in r_multi if m[2]>=0])
print(f"  probe GC: all multi-read molecules mean {allgc.mean():.3f}; lost-but-recoverable mean {lostgc.mean():.3f}")
for lo,hi in [(0,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,1.01)]:
    sel=[m for m in multi if m[2]>=0 and lo<=gc[m[2]]<hi]; l=sum(1 for m in sel if m[3]==0 and 2<=m[1]<=5)
    print(f"    GC {lo:.1f}-{hi:.1f}: molecules {len(sel):>9,}  lost-recoverable {l/max(1,len(sel)):.3%}")
genes=Counter(probes[m[2]][2] for m in r_multi if m[2]>=0); base=Counter(probes[m[2]][2] for m in multi if m[2]>=0)
print("  genes with most lost-recoverable multi-read molecules (lost / all, loss%):")
for g,c in genes.most_common(15): print(f"    {g:10s} {c:>6,} / {base[g]:>8,}  {c/base[g]:.1%}")
