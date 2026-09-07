#!/usr/bin/env python3
"""Small-sample test set for sample-tag tolerance.
Takes real JAX read pairs whose R2 carries an exact BC004 tag at offset 68 (so the
probe region is genuine and passes the hash screen), and rewrites the tag region into
five known classes, 2,000 pairs each:
  exact        the canonical BC004 tag
  variant      a 10x-listed BC004 variant (first table must accept)
  hd1_unique   one mismatch from the canonical BC004 tag, unique to BC004
  hd1_ambig    an 8-mer one mismatch from listed tags of two different samples
  shifted      the canonical tag one base later (offset 69) -> nearby search
Expected with fix on : all classes except hd1_ambig assigned to BC004.
Expected with fix off: only exact and variant assigned; the rest rejected."""
import gzip, sys, itertools, random
rows=[l.split() for l in open("/home/lhhung/jax_stage_20260903/ref/probe-barcodes-fixed-rna-profiling-rna.txt") if l.strip()]
acc={r[0]:r[2] for r in rows}; canon={r[1]:r[2] for r in rows}
bc4=[r[1] for r in rows if r[2]=="BC004"][0]; variants=[r[0] for r in rows if r[2]=="BC004" and r[0]!=bc4]
def hd(a,b): return sum(x!=y for x,y in zip(a,b))
# HD1-unique neighbour of the BC004 canonical
hd1u=None
for i in range(8):
    for b in "ACGT":
        v=bc4[:i]+b+bc4[i+1:]
        if v==bc4 or v in acc: continue
        if {acc[s] for s in acc if hd(v,s)==1}=={"BC004"}: hd1u=v; break
    if hd1u: break
# HD1-ambiguous: an 8-mer one mismatch from listed sequences of two samples (BC004 and another)
amb=None
for s in acc:
    for i in range(8):
        for b in "ACGT":
            v=s[:i]+b+s[i+1:]
            if v in acc: continue
            near={acc[t] for t in acc if hd(v,t)==1}
            if len(near)>=2 and "BC004" in near: amb=v; break
        if amb: break
    if amb: break
print("canonical BC004", bc4, "| variant", variants[0], "| hd1_unique", hd1u, "| hd1_ambig", amb, "->", {acc[t] for t in acc if hd(amb,t)==1})
S="/home/lhhung/jax_stage_20260903/fastq/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L001_"
out=sys.argv[1]; per=2000; classes=[("exact",bc4,68),("variant",variants[0],68),("hd1_unique",hd1u,68),("hd1_ambig",amb,68),("shifted",bc4,69)]
o1=gzip.open(f"{out}/test_R1.fastq.gz","wt"); o2=gzip.open(f"{out}/test_R2.fastq.gz","wt"); n=0; ci=0; k=0
with gzip.open(S+"R1_001.fastq.gz","rt") as f1, gzip.open(S+"R2_001.fastq.gz","rt") as f2:
    for (h1,s1,p1,q1),(h2,s2,p2,q2) in zip(zip(f1,f1,f1,f1),zip(f2,f2,f2,f2)):
        if s2[68:76]!=bc4: continue
        name,tag,off=classes[ci]; s=s2.rstrip("\n")
        if off==68: s=s[:68]+tag+s[76:]
        else: s=s[:68]+"A"+tag+s[77:]   # insert one base, tag now at 69, read stays 90
        o1.write(f"@{name}:{k} 1\n{s1}{p1}{q1}"); o2.write(f"@{name}:{k} 2\n{s}\n{p2}{q2}")
        k+=1
        if k%per==0:
            ci+=1
            if ci==len(classes): break
o1.close(); o2.close(); print(f"wrote {k} pairs to {out}")
