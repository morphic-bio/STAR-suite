#!/usr/bin/env python3
import argparse, gzip, os, glob

PROBE_LEN = 8

def load_probe_kmers(probe_file):
    kmers = set()
    with open(probe_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            k = parts[0].upper()
            if len(k) == PROBE_LEN:
                kmers.add(k)
    return kmers

def first_fastq(test_dir, read_tag):
    pattern = f"*_{read_tag}_*fastq.gz"
    candidates = sorted(glob.glob(os.path.join(test_dir, pattern)))
    return candidates[0] if candidates else None

def iter_fastq_seqs(fq_path, nreads):
    read = 0
    with gzip.open(fq_path, 'rt') if fq_path.endswith('.gz') else open(fq_path, 'r') as fh:
        while read < nreads:
            h = fh.readline()
            if not h: break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual: break
            yield seq.strip().upper()
            read += 1

def sweep_offsets(seq_iter, kmers, center, window):
    results = {}
    offs = range(center - window, center + window + 1)
    seqs = list(seq_iter)
    total = len(seqs)
    for off in offs:
        if off < 0:
            results[off] = (0, total)
            continue
        matched = 0
        for s in seqs:
            if len(s) >= off + PROBE_LEN:
                k = s[off:off+PROBE_LEN]
                if k in kmers:
                    matched += 1
        results[off] = (matched, total)
    return results

def main():
    ap = argparse.ArgumentParser(description="Offset sweep for exact 8-mer probe matches")
    ap.add_argument("--test_dir", default=os.path.join(os.path.dirname(__file__), "..", "test_files", "downsampled"))
    ap.add_argument("--probe_file", default=os.path.join(os.path.dirname(__file__), "..", "tables", "probe-barcodes-fixed-rna-profiling-rna.txt"))
    ap.add_argument("--fastq", help="Explicit FASTQ path; overrides --read/--test_dir")
    ap.add_argument("--read", choices=["R1","R2","R3"], default="R2", help="Which read to scan (default R2)")
    ap.add_argument("--offset", type=int, default=68, help="Center offset (0-based)")
    ap.add_argument("--window", type=int, default=5, help="+/- window around center")
    ap.add_argument("--n", type=int, default=1000, help="Number of reads to sample")
    args = ap.parse_args()

    probe_file = os.path.abspath(args.probe_file)
    test_dir = os.path.abspath(args.test_dir)
    fastq = os.path.abspath(args.fastq) if args.fastq else first_fastq(test_dir, args.read)

    if not os.path.isfile(probe_file):
        raise SystemExit(f"Probe table not found: {probe_file}")
    if not fastq or not os.path.isfile(fastq):
        raise SystemExit(f"No FASTQ found for {args.read} under {test_dir}. Provide --fastq.")

    kmers = load_probe_kmers(probe_file)
    seq_iter = iter_fastq_seqs(fastq, args.n)
    results = sweep_offsets(seq_iter, kmers, args.offset, args.window)

    print(f"# FASTQ: {fastq}")
    print(f"# Reads sampled: {args.n}")
    print(f"# Center offset: {args.offset}, window: +/-{args.window}")
    print("offset\tmatched\ttotal\tpct")
    for off in sorted(results.keys()):
        matched, total = results[off]
        pct = (100.0 * matched / total) if total else 0.0
        print(f"{off}\t{matched}\t{total}\t{pct:.3f}")

if __name__ == "__main__":
    main()
