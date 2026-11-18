#!/usr/bin/env python3
"""Tier allowlists for --soloPipseqMode, derived from the barcode reads of a
library.

This tool reads the four tier lists off the data: a read is used when exactly
one phase (0-3 leading bases) puts the three constant linkers at their
positions with no mismatch and all four tier words are A/C/G/T; the words at
each tier position are then ranked by read count and the most frequent ones
(96 by default) are written. No reference list selects reads or sequences.
Lists given with --compare are only compared with the result afterwards. (A
set made this way from one public library is bundled in
share/star-suite/barcode_lists/pipseq_readderived/.)

  make_pipseq_tier_lists.py --barcode-fastq R1.fastq.gz --out-dir DIR
      [--top 96] [--compare LABEL=tier1,tier2,tier3,tier4 ...]

The barcode FASTQ may be gzipped or plain. When no read can be used, the tool
stops with a message and writes nothing.

Output: readderived_tier{1..4}.tsv (pass these to --soloPipseqWhitelist in
this order), tier{1..4}.counts.tsv and summary.json. Check in summary.json
that the read count of the last listed word is well above the next one
(min_over_next_ratio); a few hundred thousand reads are enough. A library
may use fewer than 96 sequences per tier; then lower --top or supply the
lists of another library of the same kit.
"""
import argparse, collections, gzip, hashlib, json, os, sys

TIERS = [(0, 8), (11, 17), (20, 26), (31, 39)]
LINKERS = [(8, "ATG"), (17, "GAG"), (26, "TCGAG")]
ACGT = set("ACGT")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def open_reads(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path, "rt")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--barcode-fastq", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top", type=int, default=96)
    ap.add_argument("--compare", action="append", default=[])
    a = ap.parse_args()
    if a.top < 1:
        sys.exit("make_pipseq_tier_lists.py: --top must be at least 1")
    compare = {}
    for item in a.compare:
        label, _, paths = item.partition("=")
        files = paths.split(",")
        if not label or len(files) != len(TIERS):
            sys.exit("make_pipseq_tier_lists.py: --compare takes LABEL=tier1,tier2,tier3,tier4, got " + item)
        try:
            compare[label] = [open(p).read().split() for p in files]
        except OSError as e:
            sys.exit("make_pipseq_tier_lists.py: cannot read a --compare list: %s" % e)
    counts = [collections.Counter() for _ in TIERS]
    status = collections.Counter()
    phases = collections.Counter()
    with open_reads(a.barcode_fastq) as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            seq = line.strip()
            status["reads"] += 1
            ok = [p for p in range(4)
                  if all(seq[p + s:p + s + len(l)] == l for s, l in LINKERS) and len(seq) >= p + 39]
            if len(ok) != 1:
                status["no_linker_phase" if not ok else "several_linker_phases"] += 1
                continue
            p = ok[0]
            words = [seq[p + s:p + e] for s, e in TIERS]
            if any(set(w) - ACGT for w in words):
                status["non_acgt_tier"] += 1
                continue
            status["used"] += 1
            phases[p] += 1
            for c, w in zip(counts, words):
                c[w] += 1
    if not status["used"]:
        sys.exit("make_pipseq_tier_lists.py: no read has the three linkers at exactly one phase with A/C/G/T tier words "
                 "(reads %d, no linker phase %d, several phases %d, non-ACGT tier %d); nothing written"
                 % (status["reads"], status["no_linker_phase"], status["several_linker_phases"], status["non_acgt_tier"]))
    os.makedirs(a.out_dir, exist_ok=True)
    summary = {"schema": "pipseq_tier_lists_from_reads/1", "barcode_fastq": a.barcode_fastq,
               "barcode_fastq_sha256": sha256(a.barcode_fastq), "script_sha256": sha256(__file__),
               "selection": "exactly one phase with all three linkers exact; tier words A/C/G/T only; no reference list",
               "status": dict(status), "phase_counts": {str(k): v for k, v in sorted(phases.items())},
               "top": a.top, "tiers": []}
    for t, c in enumerate(counts, 1):
        ranked = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))
        with open(os.path.join(a.out_dir, f"tier{t}.counts.tsv"), "w") as f:
            f.write("rank\tsequence\treads\n")
            for r, (w, n) in enumerate(ranked, 1):
                f.write(f"{r}\t{w}\t{n}\n")
        top = [w for w, _ in ranked[:a.top]]
        topset = set(top)
        with open(os.path.join(a.out_dir, f"readderived_tier{t}.tsv"), "w") as f:
            f.write("\n".join(sorted(top)) + "\n")
        n_top = [n for _, n in ranked[:a.top]]
        tail = []
        for r, (w, n) in enumerate(ranked[a.top:a.top + 12], a.top + 1):
            near = sorted(((hamming(w, x), -c[x], x) for x in top))[0]
            tail.append({"rank": r, "sequence": w, "reads": n, "nearest_top_distance": near[0],
                         "nearest_top_sequence": near[2], "nearest_top_reads": -near[1]})
        row = {"tier": t, "distinct_words": len(c), "reads_in_top": sum(n_top),
               "fraction_of_used_reads_in_top": sum(n_top) / max(1, status["used"]),
               "top_reads_max": n_top[0], "top_reads_median": sorted(n_top)[len(n_top) // 2],
               "top_reads_min": n_top[-1], "next_reads": ranked[a.top][1] if len(ranked) > a.top else 0,
               "min_over_next_ratio": n_top[-1] / max(1, ranked[a.top][1]) if len(ranked) > a.top else None,
               "words_with_at_least_half_of_top_min": sum(1 for _, n in ranked if n >= n_top[-1] / 2),
               "pairs_one_substitution_apart_in_top": sum(1 for i, x in enumerate(top) for y in top[i + 1:] if hamming(x, y) == 1),
               "following_words": tail,
               "sha256_sorted_list": hashlib.sha256(("\n".join(sorted(top)) + "\n").encode()).hexdigest(),
               "compare": {label: {"entries": len(lists[t - 1]), "same_set": set(lists[t - 1]) == topset,
                                   "shared": len(set(lists[t - 1]) & topset)} for label, lists in compare.items()}}
        summary["tiers"].append(row)
    with open(os.path.join(a.out_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    for row in summary["tiers"]:
        print(row["tier"], "top min/next", row["top_reads_min"], row["next_reads"], "ratio %.1f" % (row["min_over_next_ratio"] or 0),
              "frac in top %.4f" % row["fraction_of_used_reads_in_top"], row["compare"])
    print(summary["status"], summary["phase_counts"])

if __name__ == "__main__":
    main()
