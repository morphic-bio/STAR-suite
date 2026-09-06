import numpy as np, csv, h5py, gzip
from pathlib import Path
CACHE="/home/lhhung/jax_stage_20260903/ref/h01_cache.bin"; PS="/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv"
hdr=np.dtype([("magic","S8"),("version","<u2"),("k","<u2"),("recsize","<u4"),("n","<u8")]); h=np.fromfile(CACHE,dtype=hdr,count=1)[0]
rec=np.dtype([("lo","<u8"),("hi","<u8"),("gene","<u4"),("cls","u1"),("neg","u1"),("res","<u2")]); R=np.fromfile(CACHE,dtype=rec,offset=hdr.itemsize)
print(f"cache {h['magic']} v{h['version']} k={h['k']} records {len(R):,} (header says {h['n']:,})")
for c in (0,1,2): m=R["cls"]==c; print(f"  class {c}: {m.sum():>10,} records, DENY(neg!=0) {int((R['neg'][m]!=0).sum()):>9,}")
key=(R["hi"].astype(np.uint64)<<np.uint64(64)) if False else None
# lookup structure: sort by (hi,lo)
order=np.lexsort((R["lo"],R["hi"])); hi=R["hi"][order]; lo=R["lo"][order]; neg=R["neg"][order]; gene=R["gene"][order]&0x7FFF; cls=R["cls"][order]
def enc(s):
    lo=0; hi=0
    for ch in s:
        c={"A":0,"C":1,"G":2,"T":3}.get(ch)
        if c is None: return None
        hi=((hi<<2)|(lo>>62))&((1<<64)-1); lo=((lo<<2)|c)&((1<<64)-1)
    return lo,hi
def find(s):
    e=enc(s)
    if e is None: return None
    l,hh=e; i=np.searchsorted(hi,hh); 
    while i<len(hi) and hi[i]==hh:
        if lo[i]==l: return ("DENY" if neg[i] else "KEEP", int(gene[i]), int(cls[i]))
        i+=1
    return None
comp={"A":"T","C":"G","G":"C","T":"A"}
rc=lambda s:"".join(comp.get(c,"N") for c in reversed(s))
status={}; n=0
for row in csv.reader(l for l in open(PS) if not l.startswith("#")):
    if row[0]=="gene_id": continue
    gid,seq,pid,inc=row[0],row[1],row[2],row[3]
    if inc!="TRUE": continue
    r=find(seq); orient="fwd"
    if r is None: r=find(rc(seq)); orient="rc"
    status[pid]=(r,orient); n+=1
vals=[v[0][0] if v[0] else "ABSENT" for v in status.values()]
from collections import Counter; print("included probes:",n, Counter(vals), "orientation:",Counter(v[1] for v in status.values() if v[0]))
deny=[p for p,v in status.items() if v[0] and v[0][0]=="DENY"]; absent=[p for p,v in status.items() if v[0] is None]
print("DENY probes (H0 window denied):",len(deny)); print("ABSENT probes (no H0 record):",len(absent), absent[:5])
# price them with CR molecule_info (WT-Day-8, common cells)
CR=Path("/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs/WT-Day-8/count")
f=h5py.File(CR/"sample_molecule_info.h5"); pidx=f["probe_idx"][:]; bidx=f["barcode_idx"][:]
pid=np.array([x.decode() for x in f["probes"]["probe_id"][:]]); bcs=np.array([b.decode()[:16] for b in f["barcodes"][:]])
crcells=set(l.strip().split("-")[0][:16] for l in gzip.open(CR/"sample_filtered_feature_bc_matrix"/"barcodes.tsv.gz","rt"))
ours=set(l.strip()[:16] for l in open("/storage/jax_tagfix/A_fix_on/per_sample/BC007/Gene/filtered/barcodes.tsv"))
common=crcells&ours; incell=np.isin(bcs[bidx],np.array(sorted(common)))
tot=np.bincount(pidx[incell],minlength=len(pid)); alltot=tot.sum()
dset=set(deny); aset=set(absent)
dm=sum(tot[i] for i,p in enumerate(pid) if p in dset); am=sum(tot[i] for i,p in enumerate(pid) if p in aset)
print(f"WT-Day-8 common cells: CR molecules {alltot:,}; on DENY probes {dm:,} ({dm/alltot:.3%}); on ABSENT probes {am:,} ({am/alltot:.3%}); observed total deficit 486,788 (0.646%)")
top=sorted(((tot[i],p) for i,p in enumerate(pid) if p in dset),reverse=True)[:25]
print("DENY probes by CR molecules (common cells):"); [print(f"   {p:34s} {t:>8,}") for t,p in top]
open("/mnt/pikachu/nvme_jax_v184/deny_probes_WT-Day-8.txt","w").write("\n".join(f"{p}\t{tot[list(pid).index(p)] if p in set(pid) else 0}" for p in deny))
