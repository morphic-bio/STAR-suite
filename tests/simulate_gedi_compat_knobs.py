#!/usr/bin/env python3
"""
Cheap simulation to test whether "GEDI SNP compat" knobs are worth turning, without rerunning STAR/GEDI.

Idea:
  - Take divergent loci (STAR-only and GEDI-only BEDs).
  - For a sample of loci, compute pileup counts from the SAME BAM (samtools mpileup):
      n = depth
      k_any = number of non-ref base observations (mismatch-agnostic)
      k_conv = number of "conversion" observations:
              if ref==T: count C/c
              if ref==A: count G/g
              else: 0
  - Recompute binomial tail p-value under different models:
      STAR-current-ish: use k_conv with p_err_tc (user-provided)
      GEDI-compat-ish:  use k_any with p_err_any (either user-provided or estimated from random loci)
  - Report how many loci would flip call status between the two.

This is not a perfect emulation of GEDI (transcript-conditioned counting etc.), but it is a fast
sanity check to see if switching from conversion-only -> mismatch-agnostic could plausibly close the gap.
"""

from __future__ import annotations

import argparse
import math
import os
import random
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class Locus:
    chrom: str
    start0: int  # 0-based
    end: int     # 0-based half-open (start0+1)

    @property
    def pos1(self) -> int:
        return self.start0 + 1


@dataclass
class PileupCounts:
    chrom: str
    pos1: int
    ref: str
    depth: int
    k_any: int
    k_conv: int


def logsumexp(vals: List[float]) -> float:
    m = max(vals)
    if math.isinf(m):
        return m
    s = sum(math.exp(v - m) for v in vals)
    return m + math.log(s)


def log_binom_pmf(n: int, k: int, p: float) -> float:
    # log( nCk p^k (1-p)^(n-k) )
    if k < 0 or k > n:
        return float("-inf")
    if p <= 0.0:
        return 0.0 if k == 0 else float("-inf")
    if p >= 1.0:
        return 0.0 if k == n else float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * math.log(p) + (n - k) * math.log(1.0 - p)


def binom_sf(n: int, k: int, p: float) -> float:
    """Upper tail P[X >= k] for Bin(n,p), computed in log-space by summing pmf terms."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    # sum from k..n; for cheapness, do direct sum in log-space (n will be modest for sampled loci)
    logs = [log_binom_pmf(n, i, p) for i in range(k, n + 1)]
    return math.exp(logsumexp(logs))


def parse_bed3(path: str, max_rows: Optional[int] = None) -> List[Locus]:
    loci: List[Locus] = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start0 = int(parts[1])
            end = int(parts[2])
            if end != start0 + 1:
                # force to 1bp loci for mpileup -l
                end = start0 + 1
            loci.append(Locus(chrom=chrom, start0=start0, end=end))
            if max_rows is not None and len(loci) >= max_rows:
                break
    return loci


def write_bed_for_mpileup(loci: List[Locus], out_path: str) -> None:
    with open(out_path, "w") as out:
        for l in loci:
            out.write(f"{l.chrom}\t{l.start0}\t{l.end}\n")


def count_bases_from_mpileup_readbases(readbases: str, ref_base: str) -> Tuple[int, int]:
    """
    Returns (k_any, k_conv) counts from mpileup readbases string.
    mpileup readbases semantics:
      - '.' ',' are matches to ref on forward/reverse
      - 'ACGTNacgtn' are mismatching observed bases
      - '^' + next char: start of read + mapping quality
      - '$': end of read
      - '+'/'-' followed by number and inserted/deleted sequence
      - '*' is a placeholder for a deletion in earlier read
    """
    ref = ref_base.upper()
    k_any = 0
    k_conv = 0
    i = 0
    while i < len(readbases):
        c = readbases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            j = i
            while j < len(readbases) and readbases[j].isdigit():
                j += 1
            if j == i:
                # malformed; bail
                break
            indel_len = int(readbases[i:j])
            i = j + indel_len
            continue
        if c == "*":
            # deletion placeholder; not a base observation
            i += 1
            continue
        if c in ".,":  # match to ref
            i += 1
            continue
        if c.isalpha():
            base = c.upper()
            # mismatch-agnostic count
            if base != ref and base in ("A", "C", "G", "T", "N"):
                k_any += 1
            # conversion count for ref T/A only
            if ref == "T" and base == "C":
                k_conv += 1
            elif ref == "A" and base == "G":
                k_conv += 1
            i += 1
            continue
        # unknown char
        i += 1
    return k_any, k_conv


def run_mpileup(bam: str, ref_fa: str, bed_path: str, max_depth: int = 1000000) -> List[PileupCounts]:
    cmd = [
        "samtools",
        "mpileup",
        "-f",
        ref_fa,
        "-l",
        bed_path,
        "-q",
        "0",
        "-Q",
        "0",
        "-d",
        str(max_depth),
        bam,
    ]
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    out: List[PileupCounts] = []
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom = parts[0]
        pos1 = int(parts[1])
        ref = parts[2]
        depth = int(parts[3])
        readbases = parts[4]
        k_any, k_conv = count_bases_from_mpileup_readbases(readbases, ref)
        out.append(PileupCounts(chrom=chrom, pos1=pos1, ref=ref, depth=depth, k_any=k_any, k_conv=k_conv))
    return out


def estimate_p_any_mismatch(bam: str, ref_fa: str, genome_fai: str, sample_n: int, seed: int) -> float:
    """
    Very cheap baseline estimate for mismatch rate p_any by sampling random positions genome-wide.
    This is not perfect (includes repeats, introns, etc.) but it's enough to see if the knob matters.
    """
    # load contigs and lengths
    contigs: List[Tuple[str, int]] = []
    with open(genome_fai, "r") as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contigs.append((parts[0], int(parts[1])))
    if not contigs:
        raise RuntimeError(f"Failed to load contigs from {genome_fai}")

    rng = random.Random(seed)
    loci: List[Locus] = []
    for _ in range(sample_n):
        chrom, ln = contigs[rng.randrange(len(contigs))]
        if ln <= 2:
            continue
        pos0 = rng.randrange(0, ln - 1)
        loci.append(Locus(chrom=chrom, start0=pos0, end=pos0 + 1))
    tmp_bed = f"/tmp/p_any_sample_{os.getpid()}.bed"
    try:
        write_bed_for_mpileup(loci, tmp_bed)
        counts = run_mpileup(bam, ref_fa, tmp_bed)
        n_tot = 0
        k_tot = 0
        for c in counts:
            n_tot += c.depth
            k_tot += c.k_any
        return 0.0 if n_tot == 0 else (k_tot / n_tot)
    finally:
        try:
            os.remove(tmp_bed)
        except OSError:
            pass


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True, help="Reference FASTA for samtools mpileup (-f)")
    ap.add_argument("--ref-fai", required=True, help="Reference FASTA .fai for random p_any estimation")
    ap.add_argument("--star-only-bed", required=True)
    ap.add_argument("--gedi-only-bed", required=True)
    ap.add_argument("--sample-per-set", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--p-tc", type=float, default=0.001, help="Baseline p for conversion-only model")
    ap.add_argument("--p-any", type=float, default=-1.0, help="Baseline p for mismatch-agnostic model; -1 to estimate")
    ap.add_argument("--p-any-sample-n", type=int, default=2000, help="How many random positions to use for p_any estimation")
    ap.add_argument("--pval", type=float, default=0.001)
    ap.add_argument("--min-ratio", type=float, default=0.3)
    ap.add_argument("--min-cov", type=int, default=6)
    ap.add_argument("--min-alt", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    if args.p_any < 0.0:
        p_any = estimate_p_any_mismatch(args.bam, args.ref, args.ref_fai, args.p_any_sample_n, args.seed)
    else:
        p_any = args.p_any

    star_only = parse_bed3(args.star_only_bed)
    gedi_only = parse_bed3(args.gedi_only_bed)
    rng.shuffle(star_only)
    rng.shuffle(gedi_only)
    star_s = star_only[: args.sample_per_set]
    gedi_s = gedi_only[: args.sample_per_set]

    loci_all = star_s + gedi_s
    bed_path = args.out_prefix + ".sampled_loci.bed"
    write_bed_for_mpileup(loci_all, bed_path)
    counts = run_mpileup(args.bam, args.ref, bed_path)

    # index by (chrom,pos1)
    by_key: Dict[Tuple[str, int], PileupCounts] = {(c.chrom, c.pos1): c for c in counts}

    def called(n: int, k: int, p: float) -> Tuple[bool, float]:
        if n < args.min_cov or k < args.min_alt:
            return False, 1.0
        ratio = k / n if n else 0.0
        if ratio < args.min_ratio:
            return False, 1.0
        pv = binom_sf(n, k, p)
        return pv < args.pval, pv

    rows_out = []
    for l in loci_all:
        key = (l.chrom, l.pos1)
        c = by_key.get(key)
        if c is None:
            continue
        star_call, star_pv = called(c.depth, c.k_conv, args.p_tc)
        gedi_call, gedi_pv = called(c.depth, c.k_any, p_any)
        group = "star_only" if l in star_s else "gedi_only"
        rows_out.append((group, c.chrom, c.pos1, c.ref, c.depth, c.k_any, c.k_conv, star_call, star_pv, gedi_call, gedi_pv))

    tsv_path = args.out_prefix + ".tsv"
    with open(tsv_path, "w") as out:
        out.write("group\tchrom\tpos1\tref\tdepth\tk_any\tk_conv\tcall_conv\tp_conv\tcall_any\tp_any\n")
        for r in rows_out:
            out.write(
                f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t"
                f"{1 if r[7] else 0}\t{r[8]:.3e}\t{1 if r[9] else 0}\t{r[10]:.3e}\n"
            )

    def summarize(group: str):
        gr = [r for r in rows_out if r[0] == group]
        if not gr:
            return
        flips_to_any = sum((not r[7]) and r[9] for r in gr)   # not called under conv, called under any
        flips_to_conv = sum(r[7] and (not r[9]) for r in gr)   # called under conv, not under any
        called_conv = sum(r[7] for r in gr)
        called_any = sum(r[9] for r in gr)
        total = len(gr)
        print(f"{group}: N={total}  called(conv)={called_conv}  called(any)={called_any}  flip conv->any={flips_to_any}  flip any->conv={flips_to_conv}")

    print("=== Simulation summary ===")
    print(f"bam={args.bam}")
    print(f"p_tc={args.p_tc}  p_any={p_any:.6f}  (p_any_estimated={args.p_any < 0.0})")
    print(f"filters: minCov={args.min_cov} minAlt={args.min_alt} minRatio={args.min_ratio} pval<{args.pval}")
    summarize("star_only")
    summarize("gedi_only")
    print(f"wrote: {tsv_path}")
    print(f"wrote: {bed_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

