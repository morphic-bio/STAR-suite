#!/usr/bin/env python3
"""Subset a (optionally gzipped) FASTQ to reads whose ID is in a set.

usage: subset_fastq_by_id.py ids.txt in.fastq[.gz] out.fastq.gz
IDs are bare read names (no leading @, no /1 /2 mate suffix), as emitted by
`samtools view ... | cut -f1`. Mate suffixes in the FASTQ header are stripped
before matching, so the same ID set subsets both R1 and R2.
"""
import sys, gzip

ids = {l.strip() for l in open(sys.argv[1]) if l.strip()}
inp, out = sys.argv[2], sys.argv[3]
opener = (lambda p: gzip.open(p, "rt")) if inp.endswith(".gz") else (lambda p: open(p))

kept = tot = 0
with opener(inp) as fh, gzip.open(out, "wt") as o:
    while True:
        h = fh.readline()
        if not h:
            break
        s = fh.readline(); plus = fh.readline(); q = fh.readline()
        tot += 1
        name = h[1:].split(None, 1)[0]
        if name[-2:] in ("/1", "/2"):
            name = name[:-2]
        if name in ids:
            o.write(h); o.write(s); o.write(plus); o.write(q); kept += 1
sys.stderr.write(f"{out}: kept {kept}/{tot} reads\n")
