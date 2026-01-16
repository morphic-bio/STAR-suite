#!/usr/bin/env python3
import argparse, pandas as pd, gzip, re

def read_gtf(path):
    fh = gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')
    for line in fh:
        if line.startswith('#'):
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
        if feature != 'exon':
            continue
        gd = {}
        for m in re.finditer(r'(\S+) \"([^\"]+)\";', attrs):
            gd[m.group(1)] = m.group(2)
        gene_id = gd.get('gene_id')
        gene_name = gd.get('gene_name', gene_id)
        yield (chrom, int(start)-1, int(end), gene_id, gene_name, strand)
    fh.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--probeset', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    print(f"Building probe_genes_exons.bed from {args.probeset} and {args.gtf}")
    probes = pd.read_csv(args.probeset, comment='#')
    gene_col = None
    for c in probes.columns:
        if c.lower() in ('gene_id','gene','ensgene','ensgene_id'):
            gene_col = c; break
    if gene_col is None:
        for c in probes.columns:
            if c.lower() in ('gene_name','gene_symbol','symbol'):
                gene_col = c; break
    if gene_col is None:
        raise SystemExit("Could not find gene_id or gene_name in probeset.csv")

    genes = set(probes[gene_col].dropna().astype(str).unique())

    with open(args.out, 'w') as out:
        for chrom, start, end, gene_id, gene_name, strand in read_gtf(args.gtf):
            key = gene_id if any(k.startswith('ENSG') for k in genes) else gene_name
            if key in genes:
                out.write(f\"{chrom}\\t{start}\\t{end}\\t{gene_id}\\t0\\t{strand}\\n\")
