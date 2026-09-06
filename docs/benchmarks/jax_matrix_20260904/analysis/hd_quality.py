"""Quality profile of reads by probe Hamming distance: are HD2-5 reads low-quality reads with correlated UMI/CB errors?"""
import gzip, itertools, csv, numpy as np
from collections import defaultdict, Counter
FQ="/home/lhhung/jax_stage_20260903/fastq/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_"; PS="/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv"
seqs=[]; lhs=defaultdict(list); rhs=defaultdict(list); quarters=[defaultdict(list) for _ in range(4)]; cuts=[(0,12),(12,25),(25,37),(37,50)]
for row in csv.reader(l for l in open(PS) if not l.startswith("#")):
    if row[0]=="gene_id" or row[3]!="TRUE": continue
    i=len(seqs); seqs.append(row[1])
    for q,(a,b) in enumerate(cuts): quarters[q][row[1][a:b]].append(i)
def hd(a,b): return sum(x!=y for x,y in zip(a,b))
def best(w):
    c=set()
    for q,(a,b) in enumerate(cuts): c.update(quarters[q].get(w[a:b],()))
    if not c: return 99
    return min(hd(w,seqs[i]) for i in c)
acc=defaultdict(lambda: {"n":0,"qprobe":0.0,"qumi":0.0,"qcb":0.0,"umi_lowq":0,"cb_lowq":0,"umi_min":0.0})
WL=set(l.strip() for l in open("/home/lhhung/jax_stage_20260903/ref/737K-fixed-rna-profiling.txt"))
cbwl=Counter()
with gzip.open(FQ+"R1_001.fastq.gz","rt") as f1, gzip.open(FQ+"R2_001.fastq.gz","rt") as f2:
    it1=itertools.islice(f1,1,None,1); it2=itertools.islice(f2,1,None,1)
    n=0
    while n<5_000_000:
        try:
            s1=next(it1); next(it1); q1=next(it1).rstrip(); next(it1)
            s2=next(it2); next(it2); q2=next(it2).rstrip(); next(it2)
        except StopIteration: break
        n+=1; w=s2[:50]
        if "N" in w: cls="N"
        else:
            d=best(w); cls=("HD%d"%d) if d<=5 else ("HD6+" if d<99 else "none")
        a=acc[cls]; a["n"]+=1
        qp=np.frombuffer(q2[:50].encode(),dtype=np.uint8)-33; qu=np.frombuffer(q1[16:28].encode(),dtype=np.uint8)-33; qc=np.frombuffer(q1[:16].encode(),dtype=np.uint8)-33
        a["qprobe"]+=qp.mean(); a["qumi"]+=qu.mean(); a["qcb"]+=qc.mean(); a["umi_lowq"]+=int((qu<20).any()); a["cb_lowq"]+=int((qc<20).any()); a["umi_min"]+=qu.min()
        cbwl[(cls, s1[:16] in WL)]+=1
print(f"{n:,} reads")
print(f"{'class':6s} {'reads':>9s} {'frac':>7s} {'meanQ probe':>11s} {'meanQ UMI':>9s} {'minQ UMI':>8s} {'UMI has Q<20':>12s} {'CB has Q<20':>11s} {'CB in whitelist':>15s}")
for cls in ("HD0","HD1","HD2","HD3","HD4","HD5","HD6+","N","none"):
    a=acc.get(cls); 
    if not a or a["n"]==0: continue
    k=a["n"]; inwl=cbwl[(cls,True)]/(cbwl[(cls,True)]+cbwl[(cls,False)])
    print(f"{cls:6s} {k:>9,} {k/n:>7.3%} {a['qprobe']/k:>11.1f} {a['qumi']/k:>9.1f} {a['umi_min']/k:>8.1f} {a['umi_lowq']/k:>12.3f} {a['cb_lowq']/k:>11.3f} {inwl:>15.3f}")
