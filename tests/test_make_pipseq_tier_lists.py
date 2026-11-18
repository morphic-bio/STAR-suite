#!/usr/bin/env python3
"""Checks scripts/make_pipseq_tier_lists.py on generated barcode reads.

  python3 tests/test_make_pipseq_tier_lists.py --out tests/pipseq_tier_lists_output

The reads carry known tier words in the layout of the newer PIP chemistry
(tier 1, ATG, tier 2, GAG, tier 3, TCGAG, tier 4, after 0 to 3 leading
bases), with rarer words one substitution away, reads without linkers and
reads with N in a tier. The tool must return exactly the known words.
No external data is read. --tool PATH runs the checks on another copy of the
tool.
"""
import argparse
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
TOOL = ROOT / "scripts" / "make_pipseq_tier_lists.py"
LENGTHS = (8, 6, 6, 8)
LINKERS = ("ATG", "GAG", "TCGAG")


def run(argv, ok=True):
    proc = subprocess.run([sys.executable, str(run.tool)] + [str(x) for x in argv], capture_output=True, text=True)
    assert (proc.returncode == 0) == ok, proc.stdout + proc.stderr
    return proc


def words(rng, length, n, taken):
    """n words of this length, every pair at least two substitutions apart."""
    out = []
    while len(out) < n:
        w = "".join(rng.choice("ACGT") for _ in range(length))
        if all(sum(a != b for a, b in zip(w, x)) >= 2 for x in out + taken):
            out.append(w)
    return out


def read(tiers, phase, rng):
    body = tiers[0] + LINKERS[0] + tiers[1] + LINKERS[1] + tiers[2] + LINKERS[2] + tiers[3]
    return "".join(rng.choice("CT") for _ in range(phase)) + body + "ACG" + "TTTT"


def write_fastq(path, sequences, gz=True):
    opener = gzip.open if gz else open
    with opener(path, "wt") as handle:
        for i, sequence in enumerate(sequences):
            handle.write("@r%d\n%s\n+\n%s\n" % (i, sequence, "I" * len(sequence)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--tool", type=Path, default=TOOL)
    args = ap.parse_args()
    run.tool = args.tool.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rng = random.Random(4127)
    top = 5
    known = [words(rng, length, top, []) for length in LENGTHS]

    reads, phases = [], [0, 0, 0, 0]
    for i in range(4000):
        phase = i % 4
        phases[phase] += 1
        reads.append(read([rng.choice(tier) for tier in known], phase, rng))
    # rarer words one substitution from a known word, as sequencing errors give
    errors = 0
    for tier in range(4):
        for word in known[tier]:
            changed = ("A" if word[0] != "A" else "C") + word[1:]
            for _ in range(6):
                tiers = [rng.choice(t) for t in known]
                tiers[tier] = changed
                reads.append(read(tiers, 0, rng))
                errors += 1
    phases[0] += errors
    no_linker = ["".join(rng.choice("ACGT") for _ in range(46)) for _ in range(300)]
    no_linker = [s for s in no_linker if not any(all(s[p + o:p + o + len(l)] == l for o, l in zip((8, 17, 26), LINKERS)) for p in range(4))]
    short = [read([t[0] for t in known], 0, rng)[:30]] * 5
    with_n = [read(["N" + known[0][0][1:], known[1][0], known[2][0], known[3][0]], 1, rng)] * 7
    everything = reads + no_linker + short + with_n
    rng.shuffle(everything)
    fastq = out / "barcode.fastq.gz"
    write_fastq(fastq, everything)

    # the known words come back, whatever the order of the reads
    d = out / "lists"
    run(["--barcode-fastq", fastq, "--out-dir", d, "--top", top])
    for tier in range(4):
        assert (d / ("readderived_tier%d.tsv" % (tier + 1))).read_text().split() == sorted(known[tier]), tier
    summary = json.loads((d / "summary.json").read_text())
    assert summary["status"] == {"reads": len(everything), "used": len(reads), "no_linker_phase": len(no_linker) + len(short),
                                 "non_acgt_tier": len(with_n)}, summary["status"]
    assert summary["phase_counts"] == {str(p): n for p, n in enumerate(phases)}, summary["phase_counts"]
    for row in summary["tiers"]:
        assert row["next_reads"] == 6 and row["min_over_next_ratio"] > 50, row
        assert row["pairs_one_substitution_apart_in_top"] == 0
        assert all(f["nearest_top_distance"] == 1 for f in row["following_words"][:top]), row["following_words"]
        assert abs(row["fraction_of_used_reads_in_top"] - (len(reads) - errors / 4) / len(reads)) < 1e-9, row
    header = (d / "tier1.counts.tsv").read_text().splitlines()[0]
    assert header == "rank\tsequence\treads"
    rng.shuffle(everything)
    write_fastq(out / "shuffled.fastq.gz", everything)
    run(["--barcode-fastq", out / "shuffled.fastq.gz", "--out-dir", out / "shuffled", "--top", top])
    for tier in range(1, 5):
        name = "readderived_tier%d.tsv" % tier
        assert (out / "shuffled" / name).read_bytes() == (d / name).read_bytes(), name

    # a plain FASTQ gives the same lists
    write_fastq(out / "plain.fastq", everything, gz=False)
    run(["--barcode-fastq", out / "plain.fastq", "--out-dir", out / "plain", "--top", top])
    for tier in range(1, 5):
        name = "readderived_tier%d.tsv" % tier
        assert (out / "plain" / name).read_bytes() == (d / name).read_bytes(), name

    # --compare reports the same set for the lists just written and a different set otherwise
    other = out / "other_tier1.tsv"
    other.write_text("\n".join(sorted(known[0][:-1] + ["AAAAAAAA"])) + "\n")
    same = ",".join(str(d / ("readderived_tier%d.tsv" % t)) for t in range(1, 5))
    differs = ",".join([str(other)] + [str(d / ("readderived_tier%d.tsv" % t)) for t in range(2, 5)])
    run(["--barcode-fastq", fastq, "--out-dir", out / "compare", "--top", top, "--compare", "same=" + same, "--compare", "differs=" + differs])
    compared = json.loads((out / "compare" / "summary.json").read_text())["tiers"]
    assert all(row["compare"]["same"]["same_set"] for row in compared)
    assert compared[0]["compare"]["differs"] == {"entries": top, "same_set": False, "shared": top - 1}
    assert all(row["compare"]["differs"]["same_set"] for row in compared[1:])

    # fewer distinct words than --top: all of them are listed
    run(["--barcode-fastq", fastq, "--out-dir", out / "top_large", "--top", 5000])
    large = json.loads((out / "top_large" / "summary.json").read_text())
    assert all(row["next_reads"] == 0 and row["min_over_next_ratio"] is None for row in large["tiers"])
    assert set(known[0]) < set((out / "top_large" / "readderived_tier1.tsv").read_text().split())

    # errors: no usable read, a missing file, a bad --top, a bad --compare; nothing is written
    write_fastq(out / "no_linker.fastq.gz", no_linker)
    proc = run(["--barcode-fastq", out / "no_linker.fastq.gz", "--out-dir", out / "e1"], ok=False)
    assert "nothing written" in proc.stderr and "Traceback" not in proc.stderr, proc.stderr
    proc = run(["--barcode-fastq", fastq, "--out-dir", out / "e2", "--top", 0], ok=False)
    assert "--top" in proc.stderr and "Traceback" not in proc.stderr, proc.stderr
    proc = run(["--barcode-fastq", fastq, "--out-dir", out / "e3", "--compare", "x=" + str(other)], ok=False)
    assert "--compare" in proc.stderr and "Traceback" not in proc.stderr, proc.stderr
    proc = run(["--barcode-fastq", fastq, "--out-dir", out / "e4", "--compare", "x=" + ",".join([str(out / "missing.tsv")] * 4)], ok=False)
    assert "cannot read" in proc.stderr and "Traceback" not in proc.stderr, proc.stderr
    for name in ("e1", "e2", "e3", "e4"):
        assert not (out / name).exists(), name
    print("PASS make_pipseq_tier_lists: known words recovered, read order, plain FASTQ, counts, compare, errors")


if __name__ == "__main__":
    main()
