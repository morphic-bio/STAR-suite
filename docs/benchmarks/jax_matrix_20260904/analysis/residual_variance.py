#!/usr/bin/env python3
"""Decompose what still differs between STAR-Flex (fixed) and Cell Ranger, per sample.
Matrix level (filtered matrices, common cells, Ensembl-matched genes):
  per-cell total ratio, per-gene deficits, count-bin ratio (dedup signature),
  presence/absence pairs, cells called by one tool only.
Molecule level (CR molecule_info): whether our missing UMIs are CR singletons."""
import sys, gzip, numpy as np, scipy.io, scipy.sparse as sp, h5py
from pathlib import Path
CR=Path("/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs"); ST=Path(sys.argv[1])
GROUPS=[("WT-Day-7","BC004"),("PAX6-PTC-D9-Day7","BC006"),("WT-Day-8","BC007"),("PAX6-PTC-D9-Day8","BC008")]
MOLSAMPLE=sys.argv[2] if len(sys.argv)>2 else "WT-Day-8"
def opn(d,n):
    p=d/n; return gzip.open(p+".gz" if False else str(p)+".gz","rt") if not p.exists() else open(p)
def read_mex(d):
    d=Path(d); m=d/"matrix.mtx"; m=m if m.exists() else d/"matrix.mtx.gz"
    X=sp.csr_matrix(scipy.io.mmread(str(m))).T.tocsr()          # cells x genes
    bc=[l.strip().split("-")[0][:16] for l in opn(d,"barcodes.tsv")]
    rows=[l.rstrip("\n").split("\t") for l in opn(d,"features.tsv")]
    return X,bc,[r[0] for r in rows],{r[0]:(r[1] if len(r)>1 else r[0]) for r in rows}
def align(Xa,ga,Xb,gb):
    common=[g for g in ga if g in set(gb)]; ia={g:i for i,g in enumerate(ga)}; ib={g:i for i,g in enumerate(gb)}
    return Xa[:,[ia[g] for g in common]], Xb[:,[ib[g] for g in common]], common
for crname,tag in GROUPS:
    Xc,bc_c,g_c,sym=read_mex(CR/crname/"count"/"sample_filtered_feature_bc_matrix")
    Xs,bc_s,g_s,_=read_mex(ST/tag/"Gene"/"filtered")
    cs=set(bc_c); ss=set(bc_s); common=sorted(cs&ss)
    ic={b:i for i,b in enumerate(bc_c)}; is_={b:i for i,b in enumerate(bc_s)}
    C=Xc[[ic[b] for b in common]]; S=Xs[[is_[b] for b in common]]
    S,C,genes=align(S,g_s,C,g_c); S=S.astype(np.int64); C=C.astype(np.int64)
    print(f"\n===== {crname} ({tag}): CR cells {len(bc_c)}, STAR cells {len(bc_s)}, common {len(common)}, genes matched {len(genes)} =====")
    tc=np.asarray(C.sum(1)).ravel(); ts=np.asarray(S.sum(1)).ravel()
    print(f"UMIs in common cells: STAR/CR = {ts.sum()/tc.sum():.4f}  (STAR {ts.sum():,}  CR {tc.sum():,}; deficit {tc.sum()-ts.sum():,})")
    r=ts/tc; q=np.percentile(r,[1,5,25,50,75,95,99])
    print("per-cell total ratio STAR/CR percentiles 1/5/25/50/75/95/99:", " ".join(f"{v:.4f}" for v in q), f" min {r.min():.4f} max {r.max():.4f}")
    # per-gene deficits
    gc=np.asarray(C.sum(0)).ravel(); gs=np.asarray(S.sum(0)).ravel(); d=gs-gc
    order=np.argsort(d)
    tot_abs=np.abs(d).sum(); net=d.sum()
    print(f"sum over genes of |STAR-CR| = {tot_abs:,}; net = {net:,}; top-20 deficit genes account for {(-d[order[:20]]).sum()/max(1,-(d[d<0]).sum()):.1%} of all gene-level deficit")
    print("  largest deficits (STAR-CR): "+", ".join(f"{sym[genes[i]]} {d[i]:+,} ({gs[i]}/{gc[i]})" for i in order[:12]))
    print("  largest excesses (STAR-CR): "+", ".join(f"{sym[genes[i]]} {d[i]:+,} ({gs[i]}/{gc[i]})" for i in order[::-1][:8]))
    both=(gc>0)&(gs>0); print(f"  genes CR>0 & STAR==0: {((gc>0)&(gs==0)).sum()} (CR UMIs {gc[(gc>0)&(gs==0)].sum():,});  STAR>0 & CR==0: {((gs>0)&(gc==0)).sum()} (STAR UMIs {gs[(gs>0)&(gc==0)].sum():,})")
    rg=gs[both]/gc[both]; w=gc[both]
    print(f"  per-gene ratio (genes with >=100 CR UMIs): median {np.median(rg[w>=100]):.4f}, 5th pct {np.percentile(rg[w>=100],5):.4f}, 95th {np.percentile(rg[w>=100],95):.4f}")
    # count-bin analysis on (cell,gene) pairs: dedup signature
    Cc=C.tocoo(); Sd=S.tocsr()
    sv=np.asarray(Sd[Cc.row,Cc.col]).ravel()
    bins=[(1,1),(2,2),(3,5),(6,10),(11,20),(21,50),(51,100),(101,10**9)]
    print("  (cell,gene) pairs binned by CR count -> STAR sum / CR sum, and fraction of pairs where STAR==0:")
    for lo,hi in bins:
        m=(Cc.data>=lo)&(Cc.data<=hi); 
        if m.sum()==0: continue
        print(f"    CR count {lo:>3}-{hi if hi<10**9 else 'inf':<4} pairs {m.sum():>9,}  STAR/CR {sv[m].sum()/Cc.data[m].sum():.4f}  STAR==0 {np.mean(sv[m]==0):.4f}")
    Sc=S.tocoo(); cv=np.asarray(C.tocsr()[Sc.row,Sc.col]).ravel()
    print(f"  pairs STAR>0 & CR==0: {np.sum(cv==0):,} pairs, {Sc.data[cv==0].sum():,} UMIs;  pairs CR>0 & STAR==0: {np.sum(sv==0):,} pairs, {Cc.data[sv==0].sum():,} UMIs")
    # cells called by one tool only
    Xcr,bc_r,g_r,_=read_mex(CR/crname/"count"/"sample_raw_feature_bc_matrix"); ir={b:i for i,b in enumerate(bc_r)}
    only_cr=sorted(cs-ss); only_st=sorted(ss-cs)
    tcr=np.asarray(Xc.sum(1)).ravel(); t_only_cr=np.array([tcr[ic[b]] for b in only_cr]) if only_cr else np.array([])
    st_tot=np.asarray(Xs.sum(1)).ravel(); t_only_st=np.array([st_tot[is_[b]] for b in only_st]) if only_st else np.array([])
    t_only_st_inCRraw=np.array([np.asarray(Xcr[ir[b]].sum()).item() if b in ir else -1 for b in only_st]) if only_st else np.array([])
    thr_cr=tcr.min(); thr_st=st_tot.min()
    print(f"  cells only in CR: {len(only_cr)} (their CR totals: median {np.median(t_only_cr) if len(only_cr) else 0:.0f}, min {t_only_cr.min() if len(only_cr) else 0:.0f}; CR's smallest called cell {thr_cr:.0f})")
    print(f"  cells only in STAR: {len(only_st)} (their STAR totals: median {np.median(t_only_st) if len(only_st) else 0:.0f}; in CR raw: median {np.median(t_only_st_inCRraw[t_only_st_inCRraw>=0]) if len(only_st) else 0:.0f}, absent from CR raw {np.sum(t_only_st_inCRraw<0)}; STAR's smallest called cell {thr_st:.0f})")
    if crname==MOLSAMPLE:
        f=h5py.File(CR/crname/"count"/"sample_molecule_info.h5"); bidx=f["barcode_idx"][:]; fidx=f["feature_idx"][:]; cnt=f["count"][:]
        bcs=np.array([b.decode()[:16] for b in f["barcodes"][:]]); fids=np.array([x.decode() for x in f["features"]["id"][:]])
        keep=np.isin(bcs[bidx],np.array(common)); bidx=bidx[keep]; fidx=fidx[keep]; cnt=cnt[keep]
        gi={g:i for i,g in enumerate(genes)}; gmap=np.array([gi.get(x,-1) for x in fids]); gcol=gmap[fidx]; ok=gcol>=0
        ci={b:i for i,b in enumerate(common)}; crow=np.array([ci[b] for b in bcs])[bidx]
        crow=crow[ok]; gcol=gcol[ok]; cnt=cnt[ok]
        key=crow.astype(np.int64)*len(genes)+gcol
        Sarr=np.asarray(S[crow,gcol]).ravel(); Carr=np.asarray(C[crow,gcol]).ravel()
        deficit_pair=Carr-Sarr   # per molecule: the (cell,gene) pair's deficit
        print(f"  molecule_info ({crname}, common cells, matched genes): {len(cnt):,} molecules; singletons (1 read) {np.mean(cnt==1):.3f}")
        for lab,m in [("pairs with STAR==CR",deficit_pair==0),("pairs with STAR<CR",deficit_pair>0),("pairs with STAR>CR",deficit_pair<0)]:
            print(f"    {lab:22s}: molecules {m.sum():>10,}  singleton frac {np.mean(cnt[m]==1):.3f}  median reads/mol {np.median(cnt[m]):.0f}")
        # in deficit pairs, are the molecules low-read? compare read-count distribution of deficit-pair molecules vs equal pairs
        for k in (1,2,3,5,10):
            print(f"    fraction of molecules with reads<={k}: equal-pairs {np.mean(cnt[deficit_pair==0]<=k):.3f}  deficit-pairs {np.mean(cnt[deficit_pair>0]<=k):.3f}")
