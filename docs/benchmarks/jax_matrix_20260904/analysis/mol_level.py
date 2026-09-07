import sys, gzip, numpy as np, scipy.io, scipy.sparse as sp, h5py
from pathlib import Path
CR=Path("/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs"); ST=Path(sys.argv[1]); crname=sys.argv[2]; tag=sys.argv[3]
def opn(d,n):
    p=Path(d)/n; return open(p) if p.exists() else gzip.open(str(p)+".gz","rt")
def read_mex(d):
    d=Path(d); m=d/"matrix.mtx"; m=m if m.exists() else d/"matrix.mtx.gz"
    X=sp.csr_matrix(scipy.io.mmread(str(m))).T.tocsr(); bc=[l.strip().split("-")[0][:16] for l in opn(d,"barcodes.tsv")]
    rows=[l.rstrip("\n").split("\t") for l in opn(d,"features.tsv")]; return X,bc,[r[0] for r in rows],{r[0]:(r[1] if len(r)>1 else r[0]) for r in rows}
Xc,bc_c,g_c,sym=read_mex(CR/crname/"count"/"sample_filtered_feature_bc_matrix"); Xs,bc_s,g_s,_=read_mex(ST/tag/"Gene"/"filtered")
common=sorted(set(bc_c)&set(bc_s)); ic={b:i for i,b in enumerate(bc_c)}; is_={b:i for i,b in enumerate(bc_s)}
gset=set(g_s); genes=[g for g in g_c if g in gset]; gi_c={g:i for i,g in enumerate(g_c)}; gi_s={g:i for i,g in enumerate(g_s)}
C=Xc[[ic[b] for b in common]][:,[gi_c[g] for g in genes]].astype(np.int64); S=Xs[[is_[b] for b in common]][:,[gi_s[g] for g in genes]].astype(np.int64)
f=h5py.File(CR/crname/"count"/"sample_molecule_info.h5"); bidx=f["barcode_idx"][:]; fidx=f["feature_idx"][:]; cnt=f["count"][:]; pidx=f["probe_idx"][:]
bcs=np.array([b.decode()[:16] for b in f["barcodes"][:]]); fids=np.array([x.decode() for x in f["features"]["id"][:]])
ci={b:i for i,b in enumerate(common)}; brow=np.array([ci.get(b,-1) for b in bcs]); crow=brow[bidx]; keep=crow>=0
gi={g:i for i,g in enumerate(genes)}; gmap=np.array([gi.get(x,-1) for x in fids]); gcol=gmap[fidx]; keep&=gcol>=0
crow=crow[keep]; gcol=gcol[keep]; cnt=cnt[keep]; pidx=pidx[keep]
Sarr=np.asarray(S[crow,gcol]).ravel(); Carr=np.asarray(C[crow,gcol]).ravel(); dp=Carr-Sarr
print(f"{crname}: molecules in common cells x matched genes: {len(cnt):,} (CR matrix total {C.sum():,}); singleton fraction {np.mean(cnt==1):.3f}; median reads/mol {np.median(cnt):.0f}")
for lab,m in [("pairs STAR==CR",dp==0),("pairs STAR<CR",dp>0),("pairs STAR>CR",dp<0)]:
    print(f"  {lab:16s} molecules {m.sum():>10,} ({m.mean():.3%})  singleton frac {np.mean(cnt[m]==1):.3f}  reads/mol median {np.median(cnt[m]):.0f} mean {cnt[m].mean():.1f}")
print("  read-depth profile, equal pairs vs deficit pairs (fraction of molecules with reads<=k):")
for k in (1,2,3,5,10,20): print(f"    k={k:<3} equal {np.mean(cnt[dp==0]<=k):.3f}   deficit {np.mean(cnt[dp>0]<=k):.3f}")
# Deficit-pair molecules by probe: which probes carry the missing molecules? (top probes by molecules in deficit pairs relative to their totals)
pid=np.array([x.decode() for x in f["probes"]["probe_id"][:]]) if "probe_id" in f["probes"] else None
print("  probe fields:", list(f["probes"].keys()))
if pid is not None:
    tot=np.bincount(pidx,minlength=len(pid)); dft=np.bincount(pidx[dp>0],minlength=len(pid))
    excess=dft/np.maximum(tot,1); big=np.where(tot>=2000)[0]; order=big[np.argsort(-excess[big])][:15]
    print("  probes with >=2000 molecules and the highest share of molecules sitting in deficit pairs (share; deficit-pair mols / total):")
    for i in order: print(f"    {pid[i]:32s} {excess[i]:.3f}  {dft[i]:>7,}/{tot[i]:<8,}")
    for g in ("ZIC2","N4BP2","APOE","NNAT","GOLM1"):
        gid=[k for k,v in sym.items() if v==g]
        if not gid: continue
        sel=np.array([p.startswith(gid[0]) or p.split("|")[0]==g for p in pid]); 
        idx=np.where(sel)[0]
        print(f"  {g}: probes {len(idx)}: "+", ".join(f"{pid[i].split('|')[-1]} {tot[i]:,} (in deficit pairs {dft[i]:,})" for i in idx))
